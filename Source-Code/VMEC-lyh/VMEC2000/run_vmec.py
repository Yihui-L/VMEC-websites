#!/usr/bin/env python3
"""Allocate one node's CPU quota between MPI and OpenMP, then run VMEC."""
import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys

def cpu_budget():
    """Respect affinity and all visible parent cgroup CPU quotas."""
    limits = [len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 1]
    cgroups = Path("/proc/self/cgroup")
    if cgroups.exists():
        for entry in cgroups.read_text().splitlines():
            _, controllers, relative = entry.split(":", 2)
            if not controllers:
                base = Path("/sys/fs/cgroup")
                paths = [base, base / relative.lstrip("/")]
                for path in list(paths):
                    paths.extend(parent for parent in path.parents if parent == base or base in parent.parents)
                for path in set(paths):
                    quota_file = path / "cpu.max"
                    if quota_file.exists():
                        quota, period = quota_file.read_text().split()
                        if quota != "max":
                            limits.append(max(1, int(quota) // int(period)))
            elif "cpu" in controllers.split(","):
                for base in (Path("/sys/fs/cgroup/cpu"), Path("/sys/fs/cgroup/cpu,cpuacct")):
                    path = base / relative.lstrip("/")
                    while path == base or base in path.parents:
                        quota_file = path / "cpu.cfs_quota_us"
                        period_file = path / "cpu.cfs_period_us"
                        if quota_file.exists() and period_file.exists():
                            quota = int(quota_file.read_text())
                            if quota > 0:
                                limits.append(max(1, quota // int(period_file.read_text())))
                        path = path.parent
    if os.environ.get("SLURM_CPUS_PER_TASK"):
        limits.append(int(os.environ["SLURM_CPUS_PER_TASK"]))
    return max(1, min(limits))

def allocation(budget, enabled, threads, ranks, first_ns, thread_limit):
    if min(threads, ranks) < 0:
        raise ValueError("CPU_THREADS and MPI_RANKS must be nonnegative")
    if not enabled:
        threads = 1
        ranks = ranks or 1
    threads = min(threads, thread_limit)
    if ranks == 0:
        # Most VMEC stages distribute work by radial surface across MPI ranks.
        target_threads = threads or 1
        ranks = max(1, min(budget // target_threads, max(1, first_ns // 2)))
    if threads == 0:
        threads = min(thread_limit, max(1, budget // ranks))
    if ranks * threads > budget:
        raise ValueError(f"Requested {ranks} ranks x {threads} threads exceeds CPU budget {budget}")
    return ranks, threads

def run_rank():
    """Bind each MPI process before the OpenMP runtime is initialized."""
    cpus = list(map(int, os.environ["VMEC_CPUSET"].split(",")))
    threads = int(os.environ["OMP_NUM_THREADS"])
    rank = int(os.environ.get("OMPI_COMM_WORLD_LOCAL_RANK", "0"))
    selected = cpus[rank * threads:(rank + 1) * threads]
    if len(selected) != threads:
        raise RuntimeError("MPI allocation exceeds the selected CPU set")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, selected)
    os.execve(sys.argv[2], sys.argv[2:], os.environ)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="input.<case> file")
    parser.add_argument("--executable", type=Path, default=Path(__file__).absolute().parent / "xvmec2000")
    parser.add_argument("--dry-run", action="store_true", help="parse input and print allocation only")
    args = parser.parse_args()
    input_file = args.input.absolute()
    executable = args.executable.absolute()
    helper = executable.parent / "vmec_cpu_input"
    if not input_file.name.startswith("input."):
        parser.error("The input filename must start with input.")
    for key in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_RANK"):
        if key in os.environ:
            parser.error("Invoke this launcher once outside mpirun; use xvmec2000 directly inside a scheduled MPI job.")
    result = subprocess.run([str(helper), str(input_file)], text=True, capture_output=True)
    if result.returncode:
        parser.error((result.stdout + result.stderr).strip())
    records = [line.split()[1:] for line in result.stdout.splitlines() if line.startswith("VMEC_CPU_INPUT ")]
    if len(records) != 1:
        raise RuntimeError("Could not read VMEC resource controls")
    enabled, threads, ranks, first_ns = map(int, records[0])
    # Match MyEnvVariables: GETENV stores PARVMEC in CHARACTER*100.
    parvmec = os.environ.get("PARVMEC", "")[:100].rstrip(" ")
    if parvmec in ("F", "f", "FALSE", "false"):
        if ranks > 1:
            parser.error("PARVMEC=F requires a single MPI rank")
        ranks = 1
    budget = cpu_budget()
    try:
        thread_limit = int(os.environ.get("OMP_THREAD_LIMIT", sys.maxsize))
        if thread_limit <= 0:
            raise ValueError("OMP_THREAD_LIMIT must be positive")
        ranks, threads = allocation(budget, enabled, threads, ranks, first_ns, thread_limit)
    except ValueError as exc:
        parser.error(str(exc))
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(threads), VMEC_CPU_BUDGET=str(budget),
                       OMP_DYNAMIC="FALSE", OMP_PLACES="threads", OMP_PROC_BIND="spread",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")
    environment.setdefault("OMP_STACKSIZE", "64M")
    # CPU-only entry point: avoid graphics-device discovery (which may
    # block on a remote display/driver), retaining CPU and PCI topology.
    environment.setdefault("HWLOC_COMPONENTS", "-gl,-opencl")
    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else list(range(os.cpu_count() or 1))
    environment["VMEC_CPUSET"] = ",".join(map(str, cpus[:ranks * threads]))
    if os.geteuid() == 0:
        environment.update(OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    command = [sys.executable, str(Path(__file__).absolute()), "--rank-worker",
               str(executable), input_file.name]
    if ranks > 1:
        mpi = shutil.which("mpirun")
        if not mpi:
            raise RuntimeError("mpirun is required for the selected MPI allocation")
        command = [mpi, "--host", f"localhost:{ranks}", "-np", str(ranks), "--bind-to", "none"] + command
    print(json.dumps({"cpu_budget": budget, "mpi_ranks": ranks, "threads_per_rank": threads,
                      "cwd": str(input_file.parent), "command": command}, indent=2), flush=True)
    if args.dry_run:
        return
    _, hard = resource.getrlimit(resource.RLIMIT_STACK)
    resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))
    os.chdir(input_file.parent)
    os.execvpe(command[0], command, environment)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--rank-worker":
        run_rank()
    else:
        main()
