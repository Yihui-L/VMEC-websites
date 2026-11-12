# VMEC-lyh 改动说明

VMEC2000 · 三维磁流体平衡

原版 VMEC2000 对应 STELLOPT 提交 `780f5db`。本仓库的 `Source-Code/VMEC-lyh/` 保存维护版源码。本页按功能列出当前 VMEC-lyh 相对原版的全部改动。[原版源码](https://github.com/PrincetonUniversity/STELLOPT/tree/780f5db557eacc7ca1160f4ca5a10a80922e4f26/VMEC2000) · [VMEC-lyh 源码](https://github.com/Yihui-L/VMEC-websites/tree/main/Source-Code/VMEC-lyh)

维护此页时，先判断新的代码改动是否属于已有条目。同一功能或影响范围的改动并入对应条目，独立主题另列一项。每个条目随源码更新，说明当前行为、适用范围及相关验证。

[源码版本说明](SOURCE_VERSION.md) · [NCSX 算例网站](https://yihui-l.github.io/VMEC-websites/VMEC-lyh-example/)

## 章节目录

- [01 · MGRID 缓存与 MPI 读入](#mgrid)
- [02 · 自由边界真空场更新间隔](#vacuum-update)
- [03 · FTOL 与平衡输出](#convergence)
- [04 · 磁通和旋转变换的径向端点](#flux-endpoints)
- [05 · `currvmns` 径向差分索引](#current-spectrum)
- [06 · `&ADDITIONAL_LYH` 参数参考](#input)
- [07 · Fourier 变换的 OpenMP 并行](#fourier-openmp)
- [08 · NESTOR 真空积分与线性求解](#nestor)
- [09 · MPI 标量归约](#mpi-reduction)
- [10 · CPU 资源分配与启动](#cpu-launch)
- [11 · 运行计时](#timing)
- [12 · 可选 GMRES 求解路径](#solver-options)
- [13 · WOUT 种子重启](#wout-restart)

<a name="mgrid"></a>

## 01 · MGRID 缓存与 MPI 读入

- **缓存文件身份。** 缓存同时核对 MGRID 实际读取路径、`EXTCUR`、`NFP`、环向采样点数 `NV` 及诊断配置。Linux 构建比较文件的设备、inode、大小及纳秒级修改时间和状态变更时间，识别同路径覆盖或替换。读取前后身份均可确认且一致时才登记缓存。身份无法确认时重新读取。命中时恢复线圈数量、网格模式与诊断计数。 [MGRID 读入](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/LIBSTELL/Sources/Modules/mgrid_mod.f)
- **MPI 缓存和共享外场。** 求解通信器内所有秩均命中时才复用缓存，MPI 缓存也核对参与组与秩顺序。同组的新通信器可复用缓存。二进制和 NetCDF 读入按阶段同步错误，共享内存区由组内首秩初始化。广播前为线圈元数据分配空间。
- **私有数组和共享窗口。** 无通信器的读取使用私有数组，MPI 读取使用共享窗口，两者可以在同一进程切换。共享窗口按其分配组集体释放。如果新 MPI 组未覆盖原窗口的成员，读取返回状态 9，可先由原组集体释放后再切换。 [存储管理](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/LIBSTELL/Sources/Modules/mgrid_mod.f)
- **输入校验。** 读入时检查网格维度、场周期、环向采样相容性和线圈场数组。NetCDF 数值读入检查非有限值和缺失值标记。只有完整读入并通过检查的数据才能复用。
- **输入错误返回。** 缺失、损坏或与输入不相容的 MGRID 会使自由边界求解以状态 9 结束。本次求解不会生成 WOUT。MPI 各秩共同返回错误，连续算例之间也应用同一规则。 [输入处理](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Input_Output/readin.f) · [求解入口](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/TimeStep/runvmec.f) · [命令行入口](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/TimeStep/vmec.f)
- **运行时网格覆盖。** NESTOR 检查实际插值点的 LCFS 位置。越出 MGRID 的 R/Z 范围时以状态 20 结束，避免用网格边缘场替代越界点的场并输出平衡。NCSX 裁切网格在单秩和双 MPI 秩均复现了该错误返回。 [外场插值](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/NESTOR_vacuum/becoil.f)

<a name="vacuum-update"></a>

## 02 · 自由边界真空场更新间隔

- **原版行为。** 自由边界求解中的 `NVACSKIP` 控制两次完整真空场更新之间的迭代间隔。原版串行和并行 `funct3d` 各有一个只在首次使用时初始化的 `nvskip0`。同一进程读入新算例后，自适应更新间隔仍可能被上个算例的下限限制。
- **修法。** 串行和并行路径共用同一个 `nvskip0`，每次读入新算例时根据当前输入重置。`NVACSKIP` 非正时仍先按原有规则取 `NFP`，下限至少为 1。同一算例续算或切换径向网格时沿用该下限，自适应更新逻辑继续生效。[新输入初始化](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Input_Output/readin.f) · [下限状态](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/vmec_main.f) · [并行更新](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/funct3d.f) · [串行更新](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/funct3d.f)

<a name="convergence"></a>

## 03 · FTOL 与平衡输出

- **容差优先级。** 全部为零的 `FTOL_ARRAY` 表示按 `FTOL` 生成容差。`FTOL` 须为有限正数，单阶段直接使用它，多阶段从 `1e-8` 按几何序列过渡到它。显式 `FTOL_ARRAY` 的每个有效阶段须有有限正容差。`NS_ARRAY` 的有效阶段须连续、非递减且 `NS>=3`，允许同一分辨率重复求解。省略 `NS_ARRAY` 时沿用 `NSIN` 入口。WOUT 的 `ftolv` 记录实际采用的阶段容差。 [输入检查](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Input_Output/readin.f) · [WOUT 容差](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Input_Output/wrout.f)
- **平衡输出资格。** 标准 WOUT 及相应的平衡输出要求原生求解成功，当前阶段的 `fsqr`、`fsqz`、`fsql` 有限、非负且均满足 `ftolv`。自由边界在当前候选几何上强制刷新完整真空场和力归一化，再按同一原生判据复验。 [停机复验](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/TimeStep/evolve.f) · [输出检查](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Input_Output/fileout.f)
- **原子发布和旧结果归档。** 正式输出先写入同目录临时文件，成功关闭后通过重命名发布。先发布相应诊断，最后发布标准 NetCDF WOUT。发布中途失败时清理本轮已发布文件。开始新求解前将已有同名结果移到 `wout_<算例>.previous.<随机串>.nc` 等归档路径，避免将旧结果误读为本次结果。 [文件发布和归档](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Input_Output/vmec_io.f)
- **常用失败状态。** 命令行返回非零状态。`2` 表示未达收敛要求，`5` 表示输入或运行配置不合法，`9` 表示 MGRID 读入失败，`20` 表示真空场覆盖或求解失败，`21` 表示输出归档、写入或发布失败。失败时不签发本轮平衡。 [状态定义](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/vmec_params.f) · [命令行返回](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/TimeStep/vmec.f)
- **串行 m=1 力缩放。** 原版 `scale_m1` 将含额外存储点的 `NS+1` 长度系数整列赋给长度为 `NS` 的缩放数组，在开启数组检查时会触发维度错误。现在仅使用有效径向点 `1:NS`，与 MPI 路径的作用范围及原缩放公式一致。 [力缩放实现](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/residue.f90)

<a name="flux-endpoints"></a>

## 04 · 磁通和旋转变换的径向端点

- **原版行为。** 并行例程先把径向索引裁剪到内部网格 `2:ns-1`，随后用裁剪后的索引判断是否应更新旋转变换 `iotaf` 的磁轴和边界值，以及极向磁通径向导数 `chipf` 的边界值。这些分支无法进入。`chipf` 的轴值另有分支，但首进程的径向分区过短时会退化为直接使用单个半网格值。
- **修法。** 改用 `tlglob == 1` 与 `trglob == ns` 判断磁轴和最外闭合磁面（LCFS）的所属进程。各端点按串行定义，以最近和次近两个半网格值分别乘以 `3/2` 与 `-1/2` 作线性外推。`LRFP` 分支对 `1/iota` 外推后再取倒数。[并行端点](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/add_fluxes.f90#L63-L88) · [串行公式](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/add_fluxes.f90#L151-L169)
- **对 `wout` 的影响。** 原版标准 `wout` 在写出前已通过串行后处理重建 `iotaf` 和 `chipf`，其端点已按串行公式计算。这项修复使并行阶段直接读取全网格端点时也获得同样定义的值。[原版输出前重算](https://github.com/PrincetonUniversity/STELLOPT/blob/780f5db557eacc7ca1160f4ca5a10a80922e4f26/VMEC2000/Sources/Input_Output/fileout.f) · [数组写出](https://github.com/PrincetonUniversity/STELLOPT/blob/780f5db557eacc7ca1160f4ca5a10a80922e4f26/VMEC2000/Sources/Input_Output/wrout.f)

<a name="current-spectrum"></a>

## 05 · `currvmns` 径向差分索引

- **触发条件。** `LASYM` 启用时，后处理为 Nyquist 谱中的奇数极向模 `m` 计算 `currvmns`，即雅可比加权环向电流的正弦谱。构造全网格上的径向导数时，左右两个半网格的协变极向磁场系数 `B_u` 应分别使用本点的径向因子 `shalf` 归一化。
- **修法。** 原版两项都除以右侧的 `shalf(js+1)`。左侧项 `bu0` 改为 `bsubumns_(:,js)/shalf(js)`，右侧项 `bu1` 仍使用 `shalf(js+1)`。
- **影响范围。** 修正奇数 `m` 的内部电流谱及由这些值外推的相关端点。该差分在输出电流的后处理中执行。[电流谱计算](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/LIBSTELL/Sources/Modules/read_wout_mod.f90) · [输出调用](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Input_Output/wrout.f)

<a name="input"></a>

## 06 · `&ADDITIONAL_LYH` 参数参考

- **读入与写回。** CPU 资源和新增 GMRES1D 的控制项使用同一个可选输入组。省略该组时使用默认值，每份新输入重新初始化。重复组和非法参数均返回输入错误。MPI 广播和写回包含该组。[输入实现](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/LIBSTELL/Sources/Modules/vmec_input.f)

### 完整参数表

通常将可选的 `&ADDITIONAL_LYH` 放在 `&INDATA` 的结束符 `/` 之后。每个输入文件可包含一个该组。省略某一项时使用下表的默认值。 新增加速和 GMRES1D 开关均默认关闭。

Fourier 几何与平衡变量使用 `&INDATA` 的 `MPOL` 和 `NTOR`。谱凝聚、角向网格和 Nyquist 谱沿用 VMEC 的既有规则。

<a name="additional-lyh-parameters"></a>

| 变量 | 类型与默认值 | 取值与实际作用 |
| --- | --- | --- |
| `L_CPU_ACCEL` | 逻辑 `.FALSE.` | 设为真时启用 Fourier 变换和 NESTOR Green 积分的 OpenMP 路径，并允许启动器自动分配 CPU。关闭时每秩一个线程，显式 MPI 仍可用。 |
| `CPU_THREADS` | 整数 `0` | `L_CPU_ACCEL=T` 时，`0` 自动选择，正整数请求指定的每秩线程数，实际值受线程上限和 CPU 资源限制。`L_CPU_ACCEL=F` 时请求超过一个线程会报输入错误。 |
| `MPI_RANKS` | 整数 `0` | `run_vmec.py` 的 MPI 秩数。`0` 在 CPU 加速关闭时选 1，在启用时自动分配。正整数显式指定秩数。直接运行 `xvmec2000` 时，由外部启动命令决定秩数。 |
| `L_GMRES1D` | 逻辑 `.FALSE.` | 设为真时启用可选的 1D 预条件 GMRES，要求 `PRECON_TYPE=NONE` 或 `DEFAULT`，并设置有限正阈值。 |
| `GMRES1D_THRESHOLD` | 实数 `0` | GMRES1D 的残差和切换阈值，仅在 `L_GMRES1D=T` 时使用。 |

### 完整默认配置

```fortran
&ADDITIONAL_LYH
  L_CPU_ACCEL      = .FALSE.
  CPU_THREADS      = 0
  MPI_RANKS        = 0
  L_GMRES1D        = .FALSE.
  GMRES1D_THRESHOLD = 0
/
```

<a name="fourier-openmp"></a>

## 07 · Fourier 变换的 OpenMP 并行

- **径向面并行。** 四个并行布局的正、逆 Fourier 变换按径向面静态分配线程。各线程处理独立径向切片，并使用私有循环变量和临时标量。 [逆变换](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/totzsp_mod.f) · [正变换](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/tomnsp_mod.f)
- **单线程执行。** 每秩只有一个工作线程、局部只有一个径向面或关闭 CPU 加速时，转入同一源码单独编译的单线程内核。该内核使用普通 Fortran 循环，编译时关闭 OpenMP。两条路径沿用相同的变换公式和调用次序，每个径向面上的运算顺序一致。

<a name="nestor"></a>

## 08 · NESTOR 真空积分与线性求解

- **并行积分与工作区复用。** Green 核按最多 64 个源点分块，各线程独立计算输出列。每线程复用临时数组，按原源点顺序累加。单线程使用单列缓冲直接累加。 [真空势](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/NESTOR_vacuum/scalpot.f) · [Green 内核](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/NESTOR_vacuum/greenf.f)
- **核函数验证。** 合成输入在 1 和 4 线程下与原内核逐项一致。固定合成数据在 `-O2` 和 4 核亲和性下，3 组交错测量的中位数从 0.534 s 降至 0.167 s，约 3.20 倍。
- **OpenBLAS 兼容。** GNU/OpenBLAS 构建识别到 OpenBLAS 0.3.20 的 Cooperlake 内核时，NESTOR 和原生二维预条件的 LU 分解使用 `DGETF2`，求解按单个右端调用 `DGETRS`。其他后端沿用常规 LAPACK 调用。 [兼容实现](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Hessian/hessian_lapack.f90) · [真空分解](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/NESTOR_vacuum/scalpot.f)
- **真空线性求解错误。** 检查 LU 分解与回代状态，以及矩阵和求解结果是否有限。分解失败的秩跳过回代，真空参与秩同步状态，再向所有运行秩传播，失败返回 20。双秩作业中向一个秩注入分解错误，两个秩均以状态 20 结束。 [分解和缓存](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/NESTOR_vacuum/scalpot.f) · [回代和错误传播](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/NESTOR_vacuum/vacuum.f)
- **直接读取 LU 缓存。** 完整真空更新时分配矩阵构造工作区。随后 `DGETRS` 直接读取已分解的 LU 因子，省去每次求解的整矩阵复制和增量更新时的矩阵分配。Green 积分、Fourier 变换和回代的浮点运算次序保持一致。单核 OpenBLAS 合成测试中，矩阵阶数为 169、338、676 时，每次复制加回代的配对中位时间分别减少约 5.8、32、250 µs。这个测试只衡量局部线性求解。双秩 NCSX 对照中，去除复制前后的 WOUT 和 JXBOUT 逐变量一致，真空压力与边界磁场诊断表相同。 [真空回代](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/NESTOR_vacuum/vacuum.f)

<a name="mpi-reduction"></a>

## 09 · MPI 标量归约

- **MPI 通信。** 能量、归一化量和 R/Z/lambda 力残差以局部求和及打包 `Allreduce` 计算，替代仅用于求和的全径向收集。保留所用径向范围、边界残差及完整 `vp` 剖面。MPI 求和顺序变化可能产生浮点舍入差异。 [磁场与能量](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/bcovar.f) · [力残差](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/residue.f90)
- **回归验证。** 能量、归一化量和力残差通过 1、2、4 秩的合成归约测试。
- **MPI 通信器生命周期。** 在重复运行及 VSEQ 中，替换前释放旧的运行和真空通信器。参与秩数变化时调整真空分区和计时数组容量。求解完成或清理时回收本次运行资源。 [通信器与分区](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/TimeStep/parallel_vmec_module.f90) · [运行控制](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/TimeStep/runvmec.f)

<a name="cpu-launch"></a>

## 10 · CPU 资源分配与启动

CPU 控制项为 `L_CPU_ACCEL`、`CPU_THREADS` 和 `MPI_RANKS`。三项均放在 [`&ADDITIONAL_LYH`](#additional-lyh-parameters)，线程数与秩数都要求为非负整数。

### 两种启动入口的线程设置

- **自动启动器。** `run_vmec.py` 根据输入参数、可用 CPU 与 `OMP_THREAD_LIMIT` 分配资源，再设置子进程的 `OMP_NUM_THREADS`。调用者原先导出的 `OMP_NUM_THREADS` 由启动器计算值覆盖。 `L_CPU_ACCEL=F` 且 `MPI_RANKS=0` 时默认启动一秩一线程。明确设置正整数 `MPI_RANKS` 或直接使用外部 MPI 启动仍可运行原生 MPI 布局，此时每秩使用一个线程。
- **直接 MPI 启动。** 设定正整数 `CPU_THREADS` 时优先使用该值。取 `0` 时读取环境变量 `OMP_NUM_THREADS` 的首个整数。两项均未指定时，Linux 原生入口识别本节点各秩互不重叠的 CPU 亲和掩码，直接使用各秩的可见 CPU 数。亲和掩码重叠或无法可靠识别时，按本节点秩数保守分配。结果受本秩亲和性、可得的节点 CPU 配额、`OMP_THREAD_LIMIT` 和 MPI 线程级别限制。
- **旧输入项。** `&INDATA` 中的同名 `OMP_NUM_THREADS` 仍能读入。当前加速路径的线程数由 `CPU_THREADS`、环境变量及运行资源决定。

### 环境变量与资源限制

| 变量 | 作用 |
| --- | --- |
| `HWLOC_COMPONENTS` | CPU 启动器在未设置时默认 `-gl,-opencl`，跳过图形/OpenCL 拓扑发现，保留 CPU/PCI 拓扑；用户显式设置原样保留。直接运行原生入口时由外部环境决定。 |
| `OMP_THREAD_LIMIT` | OpenMP 线程硬上限，启动器和原生入口都遵循。 |
| `OMP_NUM_THREADS` | 直接启动时为 `CPU_THREADS=0` 提供线程数。自动启动器为子进程重设此变量。 |
| `SLURM_CPUS_PER_TASK` | 自动启动器将其作为可用 CPU 预算上限之一。 |
| `PARVMEC` | 默认使用并行布局。取 `F` 时使用单秩串行布局，启动器拒绝同时请求多个 MPI 秩，尾随空格按原生入口的规则处理。 |
| `OMP_DYNAMIC`、`OMP_PLACES`、`OMP_PROC_BIND` | 自动启动器分别设为 `FALSE`、`threads` 和 `spread`，配合每秩独立的 CPU 集合。 |
| `OPENBLAS_NUM_THREADS`、`MKL_NUM_THREADS`、`BLIS_NUM_THREADS` | 自动启动器设为 `1`，避免线性代数库与 OpenMP 嵌套占用线程。直接启动时由外部环境配置。 |
| `OMP_STACKSIZE` | 线程栈大小。自动启动器在用户未设置时取 `64M`。 |
| `VMEC_CPU_BUDGET`、`VMEC_CPUSET` | 启动器在进程间传递预算与 CPU 集合的内部变量。 |

### 两秩四线程示例

以下配置在可用 CPU 预算至少为 8、且线程上限至少为 4 时分配两秩、每秩四线程。使用 `run_vmec.py --dry-run` 检查实际秩数和每秩线程数。

```fortran
&ADDITIONAL_LYH
  L_CPU_ACCEL = .TRUE.
  CPU_THREADS = 4
  MPI_RANKS   = 2
/
```

- **自动分配策略。** `run_vmec.py` 通过编译后的原生输入模块读取参数，结合 CPU 亲和性、可见的 cgroup v1/v2 配额及 Slurm 每任务 CPU 限额计算可用预算。线程数与秩数都自动选择时优先分配 MPI 秩。秩数由可用预算及求解器首级径向网格的 `max(1, floor(NS/2))` 上限约束，剩余预算按每秩等量分配给 OpenMP。显式设置线程数时，自动选择的秩数还受每秩线程请求约束。超出预算的显式组合会报错。
- **单节点自动启动。** 启动器为各 MPI 进程分配互不重叠的 CPU 集合，并将 OpenBLAS、MKL 和 BLIS 内部线程设为 1。程序以 `MPI_THREAD_FUNNELED` 初始化，MPI 通信与计时由主线程执行。外部库调用只提供 `MPI_THREAD_SINGLE` 时，OpenMP 内核使用单线程。
- **直接 MPI 调用。** 多节点作业由调度器或 `mpirun` 配置进程。原生入口结合本秩的亲和掩码与共享节点配额裁剪线程数。`PARVMEC=F` 仅允许单秩，原生入口在创建输出前拒绝多 MPI 秩。VSEQ 日志与正式平衡文件由全局秩 0 写入。自动启动器面向单节点，调用应在 MPI 作业外进行。 [原生启动与输出](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/TimeStep/vmec.f)
- **构建和启动。** 提供 GNU Fortran、Open MPI、OpenBLAS、ScaLAPACK、NetCDF Fortran 的构建配置。构建前解析真实路径，拒绝通过符号链接把生成文件写回源码目录。`build_cpu.sh` 在指定的外部目录构建，并生成 `xvmec2000`、原生参数解析辅助程序和 `run_vmec.py`。Release 使用 `-O3 -fopenmp`，Debug 加入运行时检查。原生内核使用双精度。

```bash
./Source-Code/VMEC-lyh/build_cpu.sh /path/to/vmec-build release
/path/to/vmec-build/bin/run_vmec.py --dry-run \
  VMEC-lyh-example/NCSX/input.ncsx_c09r00_free_additional_lyh
/path/to/vmec-build/bin/run_vmec.py \
  VMEC-lyh-example/NCSX/input.ncsx_c09r00_free_additional_lyh
```

- **运行目录。** 启动器在输入文件所在目录运行，`MGRID_FILE` 的相对路径按该目录解释。NCSX 的 `mgrid_c09r00.nc` 需要由同目录的 MAKEGRID 输入提前生成。`--dry-run` 仅解析输入和输出资源计划。[NCSX 示例](https://github.com/Yihui-L/VMEC-websites/blob/main/VMEC-lyh-example/NCSX/input.ncsx_c09r00_free_additional_lyh)
- **构建与验证范围。** GNU Fortran、Open MPI、OpenMP 的 Debug 与 Release 构建通过。NCSX 自由边界的端到端结果见下表。

### NCSX 自由边界回归

使用同一 NCSX 自由边界输入和 MGRID，最终径向网格为 99 个磁面，最终 FTOL 为 1e-12。原版与对照版在同一核心集合上交错重复测试，当前版本用相同布局复跑。混合布局明确启用 ADDITIONAL_LYH 的 CPU 加速。各组固定 CPU 亲和性并使用单线程 OpenBLAS Haswell。共享节点的计时以中位数和范围报告。

| 版本 | CPU 配置 | 次数 | 端到端中位数（s） | 端到端范围（s） | 求解中位数（s） |
| --- | --- | --- | --- | --- | --- |
| 原版 VMEC | 串行，1 核 | 3 | 58.87 | 55.88–66.59 | 57.26 |
| 对照版 VMEC-lyh（133d504） | 串行，1 核 | 3 | 56.33 | 55.22–57.32 | 54.70 |
| 当前 VMEC-lyh | 串行，1 核 | 2 | 62.13 | 59.07–65.19 | 59.83 |
| 原版 VMEC | 4 秩 × 1 线程 | 3 | 18.20 | 17.41–19.16 | 15.83 |
| 对照版 VMEC-lyh（133d504） | 4 秩 × 1 线程 | 3 | 20.06 | 18.76–20.91 | 18.19 |
| 当前 VMEC-lyh | 4 秩 × 1 线程 | 3 | 19.81 | 18.75–21.76 | 17.12 |
| 当前 VMEC-lyh | 2 秩 × 2 线程，显式启用加速 | 2 | 26.84 | 25.81–27.87 | 24.77 |

当前版本的三种 CPU 布局均达到原 FTOL。相对同工况原版，LCFS 和磁轴最大位置差为 1.05e-07 m 和 2.37e-09 m，iotaf 最大差 3.68e-08，b0 最大差 2.12e-09 T，环向电流剖面最大差 2.97e-09 A，equif 最大差 2.65e-06。全网格压强剖面逐点相同。

在保持候选谱系数固定的独立真空复算中，默认、二维 GMRES 和 GMRES1D 三种路径的原生三项残差均达到原 FTOL。验收使用原收敛容差和相同的离散平衡方程。

<a name="timing"></a>

## 11 · 运行计时

- **MPI 计时文件。** 命令行运行由 MPI 全局秩 0 写入 `timings.txt`，记录该秩的阶段用时及通信等待。`PARVMEC=F` 时也记录实际力评估和正、逆 Fourier 变换的用时，正向对应 `tomnsp`，逆向对应 `totzsp`。 [计时输出](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Input_Output/timer_sub.f) · [力评估与 Fourier 计时](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/funct3d.f)
- **真空场积分。** 按结束时间减去起始时间累计 `scalpot` 用时。 [真空积分计时](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/NESTOR_vacuum/scalpot.f)
- **串行二维预条件器。** 修正计时日志的输出单元与实数格式。 [二维预条件计时](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Hessian/precon2d.f)

<a name="solver-options"></a>

## 12 · 可选 GMRES 求解路径

- **输入参数。** `&INDATA` 的原有 `PRECON_TYPE` 默认为 `NONE`，使用原一维预条件松弛。原有 `GMRES` 使用二维块预条件，沿用 `PREC2D_THRESHOLD`。新增 `GMRES1D` 由 `&ADDITIONAL_LYH` 的 `L_GMRES1D` 和 `GMRES1D_THRESHOLD` 独立启用，与原有二维预条件互斥。两种路径均在最终径向网格且三项原生残差之和小于对应阈值时尝试切换。自由边界还要求真空已启用。`GMRESR` 和 `TFQMR` 仅支持 `PARVMEC=F` 的单秩真串行布局。旧的 `PRECON_TYPE=GMRES1D` 会提示迁移并报输入错误。 [输入分派](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Input_Output/readin.f) · [算法切换](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/TimeStep/evolve.f)
- **独立试步评价。** GMRES 用独立的试步评价量判断残差下降，自由边界的评价量包含 LCFS 的 R/Z 力贡献一次。自由边界重新计算基点力，基点与试探点采用相同的力评价规则。拒绝试步时恢复基点状态并转回原预条件松弛，收敛仍按原三项残差和 `FTOL` 判定。 [试步接受与恢复](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/TimeStep/gmres_mod.f) · [试步残差](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/General/residue.f90)

<a name="wout-restart"></a>

## 13 · WOUT 种子重启

- **几何初值与新输入。** 使用 `reset=` 读取已有 WOUT 的 R/Z/λ 作为初值。NetCDF 只读取种子几何、λ 和相容性校验所需字段。旧文本由独立种子读取路径兼容。自由边界开关、MGRID、线圈和剖面来自新输入。种子须与目标的 `NS`、`MPOL`、`NTOR`、`NFP` 相符，并包含完整、无重复的有限谱系数。由 MPI 秩 0 读取和重构，再向各秩分发。 [种子校验和读取](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Initialization_Cleanup/load_xc_from_wout.f)
- **λ 的径向重构。** 从磁轴向外按原网格顺序反解 WOUT 的半网格 λ 插值。偶极向模使用普通径向递推，奇极向模带上近轴径向权重。先恢复全网格 λ，再按 `phipf/lamscale` 转回内部变量，使重构独立于 MPI 径向分区。 [径向递推与内部尺度](https://github.com/Yihui-L/VMEC-websites/blob/main/Source-Code/VMEC-lyh/VMEC2000/Sources/Initialization_Cleanup/load_xc_from_wout.f)
