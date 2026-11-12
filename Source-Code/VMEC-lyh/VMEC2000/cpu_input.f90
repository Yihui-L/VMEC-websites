program vmec_cpu_input
  use vmec_input, only: read_indata_namelist, cpu_threads, mpi_ranks, l_cpu_accel, ns_array, nsin
  use vparams, only: nsd
  implicit none
  integer :: unit, status, first_ns
  character(1024) :: filename
  call get_command_argument(1, filename)
  open(newunit=unit, file=trim(filename), status='old', action='read', iostat=status)
  if (status /= 0) stop 1
  call read_indata_namelist(unit, status)
  close(unit)
  if (status /= 0) stop 2
  ! Match readin's legacy NSIN branch when NS_ARRAY(1) is zero.
  if (ns_array(1) == 0) then
    first_ns = min(nsd, max(3, nsin))
  else if (ns_array(1) >= 3) then
    first_ns = min(nsd, ns_array(1))
  else
    write(*,'(a)') 'NS_ARRAY(1) must be zero for NSIN or at least 3'
    stop 2
  endif
  write(*,'(a,4(1x,i0))') 'VMEC_CPU_INPUT', merge(1,0,l_cpu_accel), cpu_threads, mpi_ranks, first_ns
end program vmec_cpu_input
