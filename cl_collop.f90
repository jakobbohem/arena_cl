!! Copyright (C) 2011 X. Andrade
!!
!! FortranCL is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! FortranCL is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

!! CHANGE THE NAME OF THIS?? THIS IS THE TESTBED OF THE OPENCL COLLISION OPERATOR
module cl_collop
contains
subroutine runClCompare(nloop, outp_times)
  use cl
!  use comek
  use collision_operator ! for in line loop compare
!  use system_clock

  implicit none
  ! input params:
  integer, intent(in) :: nloop
  real, dimension(4), intent(out) :: outp_times

  type(cl_platform_id)   :: platform
  type(cl_device_id)     :: device
  type(cl_context)       :: context
  type(cl_command_queue) :: command_queue
  type(cl_program)       :: prog
  type(cl_kernel)        :: kernel

  integer    :: num, ierr, irec, size, ivar
  integer(8) :: size_in_bytes, globalsize, localsize
  character(len = 100)  :: info
  character(len=100) :: inputfile
  integer, parameter :: iunit = 10
  integer, parameter :: source_length = 10000
  character(len = source_length) :: source
  real, allocatable  :: vec1(:), vec2(:), fsum(:)
  type(cl_mem)       :: cl_vec1, cl_vec2
  
!  real, dimension(500) :: du1, du2, du3, du4
  ! Collision operator parameters
  
  real, allocatable :: pcdot_output(:), lcdot_output(:), ac_lam_output(:), ac_p_output(:)! for old compare:
  
  ! Collision operator CL parameters
  type(cl_mem)  :: cl_peps, cl_pp, cl_plam, cl_ppcdot, cl_pacp, cl_plcdot, cl_paclam, cl_I1, cl_I2
  real, allocatable  :: peps(:), pp(:), plam(:), ppcdot(:), pacp(:), plcdot(:), paclam(:), sigma(:)
  real, allocatable :: I1(:), I2(:)
  
  real :: old_loop_count, test1, test2
  real :: temperature, Z_eff, lnLam
  
  real :: sum1, sum2, sumsum, r, itemp1, itemp2
  integer :: ri,k
  
  double precision :: fsecond ! a C-function to measure time!
  double precision :: start_time
  
  ! Variables for timing;
  integer :: t0,t1,clock_rate, clock_max, tick_time, tock_time, clock_rateT, clock_rate0, tick, tock
  real :: time0, time1, time_old0, time_old1
  
  !! TEST SUBROUTINES
  test1 = 1
  test2 = 2
  call addtwo(test1,test2,itemp1);
  if (itemp1 .ne. 3.0) then
    print *, "sum of ",test1," and ",test2," is: ",itemp1
    print *, "ERROR in function call!!"
  end if
  
  
  !=====================
  ! SET UP THE INPUT READER AND LOOP THE INTEGRALS
  !=====================
  
  call system_clock(t0,clock_rateT,clock_max) ! Measure total run, real time
  call initek
  
  size = 10000
  size_in_bytes = int(size, 8)*4_8
  allocate(vec1(1:size))
  allocate(vec2(1:size))
  allocate(fsum(1:size))
  allocate(sigma(1:size))
  
  allocate(peps(1:size))
  allocate(pp(1:size))
  allocate(plam(1:size))
  allocate(I1(1:size))
  allocate(I2(1:size))
  
  allocate(ppcdot(1:size))
  allocate(pacp(1:size))
  allocate(plcdot(1:size))
  allocate(paclam(1:size))
  
  inputfile = "arena_initial/pl4.res"
  call readState("arena_initial/pl4.res",size, pp,plam,peps,sigma)
    ! OLD SCHOOL LOOP FOR GETTING THE IFUNC:S (i1 and i2 elliptic integral)
  ! this should also be parallelised

!  A BIT STRANGE WITH INDEXING HERE...
!  print *, "peps: ",peps(0)," plam: ", plam(0)," I1: ", I1(0),"I2: ", I2(0)
!  print *, "peps: ",peps(2)," plam: ", plam(2)," I1: ", I1(2),"I2: ", I2(2)
  do ivar = 1,size
!    print *, "vec1(i): ",vec1(ivar)," vec2(i): ",vec2(ivar)
!    print *, "peps: ",peps(ivar)," plam: ", plam(ivar)," I1: ", I1(ivar),"I2: ", I2(ivar)
    call ifunc(peps(ivar), plam(ivar) ,itemp1,itemp2);
    I1(ivar) = itemp1
    I2(ivar) = itemp2
!    print *, 'called ifunc'
  end do
  
!  print *, "A, peps: ",peps(2)," plam: ", plam(2)," I1: ", I1(2),"I2: ", I2(2)
  
    !=====================
  ! INITIALIZATION
  !=====================
  ! get the platform ID
  call clGetPlatformIDs(platform, num, ierr)
  if(ierr /= CL_SUCCESS) stop "Cannot get CL platform."

  ! get the device ID
  call clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, device, num, ierr)
  if(ierr /= CL_SUCCESS) stop "Cannot get CL device."

  ! get the device name and print it
  call clGetDeviceInfo(device, CL_DEVICE_NAME, info, ierr)
  print*, "CL device: ", info

  
  ! create the context and the command queue
  context = clCreateContext(platform, device, ierr)
  command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, ierr)
  
  !=====================
  ! BUILD THE KERNEL
  !=====================

  ! read the source file
  open(unit = iunit, file = 'coll_op.cl', access='direct', status = 'old', action = 'read', iostat = ierr, recl = 1)
  if (ierr /= 0) stop 'Cannot open file coll_op.cl'

  source = ''
  irec = 1
  do
    read(unit = iunit, rec = irec, iostat = ierr) source(irec:irec)
    if (ierr /= 0) exit
    if(irec == source_length) stop 'Error: CL source file is too big'
    irec = irec + 1
  end do
  close(unit = iunit)

  ! create the program
  prog = clCreateProgramWithSource(context, source, ierr)
  if(ierr /= CL_SUCCESS) stop 'Error: cannot create program from source.'

  ! build
  call clBuildProgram(prog, '-cl-mad-enable', ierr)

  ! get the compilation log
  call clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, source, irec)
  if(len(trim(source)) > 0) print*, trim(source)
!  print *, CL_PROGRAM_BUILD_LOG

  if(ierr /= CL_SUCCESS) stop 'Error: program build failed.'

  ! finally get the kernel and release the program
  kernel = clCreateKernel(prog, 'coll_op', ierr)
  call clReleaseProgram(prog, ierr)

  !=====================
  ! RUN THE KERNEL
  !=====================
  
  ! read in the ARENA distribution AND SIZE from file (for now, just empties)
  ! funky fortran vector initialization
!  peps = 1.0
!  pp = 1.0
!  plam = 0.3
  
  vec1 = 1.0
  vec2 = 2.0
  
  temperature = 1.0
  Z_eff = 1.0
  lnLam = 18
    
  ! loop here?
  call cpu_time(time0)
  call system_clock(tick, clock_rate0, clock_max)
  do k = 1,nloop
  
  ! allocate device memory
  cl_vec1 = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, ierr)
  cl_vec2 = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
  cl_peps = clCreateBuffer(context,CL_MEM_READ_ONLY, size_in_bytes, ierr)
  cl_pp = clCreateBuffer(context,CL_MEM_READ_ONLY, size_in_bytes, ierr)
  cl_plam = clCreateBuffer(context,CL_MEM_READ_ONLY, size_in_bytes, ierr)
  
  cl_I1 = clCreateBuffer(context,CL_MEM_READ_ONLY, size_in_bytes, ierr)
  cl_I2 = clCreateBuffer(context,CL_MEM_READ_ONLY, size_in_bytes, ierr)
  
  cl_ppcdot = clCreateBuffer(context,CL_MEM_READ_WRITE, size_in_bytes, ierr)
  cl_pacp = clCreateBuffer(context,CL_MEM_READ_WRITE, size_in_bytes, ierr)
  cl_plcdot = clCreateBuffer(context,CL_MEM_READ_WRITE, size_in_bytes, ierr)
  cl_paclam = clCreateBuffer(context,CL_MEM_READ_WRITE, size_in_bytes, ierr)
  

  ! copy data to device memory
  call clEnqueueWriteBuffer(command_queue, cl_vec1, cl_bool(.true.), 0_8, size_in_bytes, vec1(1), ierr)
  call clEnqueueWriteBuffer(command_queue, cl_vec2, cl_bool(.true.), 0_8, size_in_bytes, vec2(1), ierr)
  
  call clEnqueueWriteBuffer(command_queue, cl_I1, cl_bool(.true.), 0_8, size_in_bytes, I1(1), ierr)
  call clEnqueueWriteBuffer(command_queue, cl_I2, cl_bool(.true.), 0_8, size_in_bytes, I2(1), ierr)
  
  call clEnqueueWriteBuffer(command_queue, cl_peps, cl_bool(.true.), 0_8, size_in_bytes, peps(1), ierr)
  call clEnqueueWriteBuffer(command_queue, cl_pp, cl_bool(.true.), 0_8, size_in_bytes, pp(1), ierr)
  call clEnqueueWriteBuffer(command_queue, cl_plam, cl_bool(.true.), 0_8, size_in_bytes, plam(1), ierr)
  

  ! set the kernel arguments - REQUIRED FOR ALL PARAMETERS, in order?
  call clSetKernelArg(kernel, 0, size, ierr)
!  call clSetKernelArg(kernel, 1, cl_vec1, ierr)
!  call clSetKernelArg(kernel, 2, cl_vec2, ierr)
  call clSetKernelArg(kernel, 1, cl_peps, ierr)
  call clSetKernelArg(kernel, 2, cl_pp, ierr)
  call clSetKernelArg(kernel, 3, cl_plam, ierr)
  call clSetKernelArg(kernel, 4, cl_I1, ierr)
  call clSetKernelArg(kernel, 5, cl_I2, ierr)
  call clSetKernelArg(kernel, 6, cl_ppcdot, ierr)
  call clSetKernelArg(kernel, 7, cl_pacp, ierr)
  call clSetKernelArg(kernel, 8, cl_plcdot, ierr)
  call clSetKernelArg(kernel, 9, cl_paclam, ierr)
  call clSetKernelArg(kernel, 10, temperature, ierr)
  call clSetKernelArg(kernel, 11, Z_eff, ierr)
  call clSetKernelArg(kernel, 12, lnLam, ierr) ! 13 args, oK!

  ! get the localsize for the kernel (note that the sizes are integer(8) variable)
  call clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, localsize, ierr)
  globalsize = int(size, 8)
  if(mod(globalsize, localsize) /= 0) globalsize = globalsize + localsize - mod(globalsize, localsize) ! rounding function...

  
    ! execute the kernel
    call clEnqueueNDRangeKernel(command_queue, kernel, (/globalsize/), (/localsize/), ierr)
    call clFinish(command_queue, ierr)

    ! read the resulting vector from device memory
  !  call clEnqueueReadBuffer(command_queue, cl_vec2, cl_bool(.true.), 0_8, size_in_bytes, fsum(1), ierr)
    call clEnqueueReadBuffer(command_queue, cl_ppcdot, cl_bool(.true.), 0_8, size_in_bytes, ppcdot(1), ierr)
    call clEnqueueReadBuffer(command_queue, cl_pacp, cl_bool(.true.), 0_8, size_in_bytes, pacp(1), ierr)
    call clEnqueueReadBuffer(command_queue, cl_plcdot, cl_bool(.true.), 0_8, size_in_bytes, plcdot(1), ierr)
    call clEnqueueReadBuffer(command_queue, cl_paclam, cl_bool(.true.), 0_8, size_in_bytes, paclam(1), ierr)

  end do ! the time-step loop

  !=====================
  ! RELEASE EVERYTHING
  !=====================

  call clReleaseKernel(kernel, ierr)
  call clReleaseCommandQueue(command_queue, ierr)
  call clReleaseContext(context, ierr)

  call cpu_time(time1) ! tock for OpenCL loop.
  call system_clock(tock,clock_rate0,clock_max)
  ! Slow old school fortran loop: (time this and compare)
  
  
  sum1 = 0
  sum2 = 0
  sumsum = 0
  

  call system_clock(tick_time, clock_rate, clock_max) ! not currently using clock_max, the truncation value.
  call cpu_time(time_old0)
  
  old_loop_count = 0
  allocate(pcdot_output(1:size))
  allocate(lcdot_output(1:size))
  allocate(ac_p_output(1:size))
  allocate(ac_lam_output(1:size))
  
  do k = 1,nloop ! loop over time-steps (old)
    
    do ivar = 0,size
  !    fsum(ivar) = 0
  !    if (ivar .ne. 0) then
  !      fsum(ivar) = fsum(ivar-1) + vec2(ivar)
  !    end if
  !    sum1 = sum1+vec1(ivar)
  !    sum2 = sum2+vec2(ivar)
  !    sumsum = sumsum + fsum(ivar)

      call mccoll_v1(peps(ivar), pp(ivar), plam(ivar),I1(ivar),I2(ivar), &
      pcdot_output(ivar), lcdot_output(ivar), ac_p_output(ivar), ac_lam_output(ivar), temperature, & 
      Z_eff, lnLam)
  !    pcdot_output
  !    lcdot_output
  !    ac_p_output
  !    ac_lam_output
      old_loop_count = old_loop_count + 1
    end do
  
  end do
  
  print *, "called old_loop ",old_loop_count," times."
  
  call system_clock(tock_time, clock_rate, clock_max)
  call cpu_time(time_old1)  
  
  call random_number(r)
  ri = nint(r*size)
  
  
  ! Print the output:
  print *, 'row 5 in coll.op;  ppcdot(',ppcdot(5),'), pacp (',pacp(5),'), plcdot (',plcdot(5),'), paclam (',paclam(5),')'
  print *, 'random i=',ri,' in coll.op; ppcdot(',ppcdot(ri),'), pacp (',pacp(ri),'), plcdot (',plcdot(ri),'), paclam (',&
  paclam(ri),')'
  print *, 'compare random i in old loop; ppcdot(',ppcdot(ri),'), pacp (',pacp(ri),'), plcdot (',plcdot(ri), &
    '), paclam (',paclam(ri),')'
!  print *, 'DBG OLD CHECK: (at ',ri,') ' ,vec1(ri),' plus ',vec2(ri),' is ', fsum(ri)  
    
  print *, "OpenCL real time: ", real(tock-tick) / real(clock_rate0)
  print *, "OpenCL routine CPU time: ",time1-time0
  print *, "old loop real time: ",real(tock_time-tick_time) / real(clock_rate), " (called ",old_loop_count," times)"
  print *, "old loop CPU time: ", time_old1-time_old0
  
  call system_clock(t1, clock_rateT,clock_max) ! end of run, timing (maybe use a separate variable here)
  print *, '' !newline.
  print *, 'total elapsed time: ', real(t1-t0) / real(clock_rateT)
  
  ! Assign outputs for writing to file in containing program ! CL_real, CL_CPU, old_real, old_CPU
  outp_times(1) = real(tock-tick) / real(clock_rate0)
  outp_times(2) = time1-time0
  outp_times(3) = real(tock_time-tick_time) / real(clock_rate)
  outp_times(4) = time_old1-time_old0
  
  !=====================
  ! END OF PROGRAM 
  !===================== 
end subroutine runClCompare



!> the latter three should output the appropriate vectors!
subroutine readState(full_file, file_length, p, lam, eps, sigma)
  implicit none
  character(len=100), intent(in) :: full_file
  integer, intent(in) :: file_length
  real :: pt, lamt, epst, sigmat
  real, dimension(file_length) , intent(out) :: p, lam, eps, sigma ! change to read LENGTH OF FILE!!
  integer :: i, le, rand_row, io, c 
  
  ! DEBUG PARAMETERS:
  real :: psum, lamsum, epssum, sigmasum
  double precision :: r ! random no
  
  open(unit=11,file='arena_initial/pl4.res', status='old', iostat = io)
!  open (11,file='Godea.txt',status='old')
!  open (12,file='Zbirno.txt',status='unknown')

  ! switch end=100, to iostat = io and if io < 0 ; break;
  c = 0
  do
    read (11,*,iostat = io) p(c),lam(c),eps(c), sigma(c)
    c = c+1
    if (io < 0) then
      exit
    end if
  end do
  
  
  print *,"read ",c," rows from file " !,trim(full_file)
  close(11)
  
  !DEBUG IT:
!  psum = 0;lamsum=0;epssum = 0; sigmasum = 0
!  do i = 1,file_length
!    psum = psum + p(i)
!    lamsum = lamsum + lam(i)
!    epssum = epssum + eps(i)
!    sigmasum = sigmasum + sigma(i)
!  end do
!  
!  call random_number(r)
!  rand_row = nint(r*file_length)
!  print *, "debug: random row: ", p(rand_row),lam(rand_row), eps(rand_row), sigma(rand_row)
!  print *, "debug: row 5: ",p(5),lam(5), eps(5), sigma(5)
!  print *, "debig: sums: ", psum, lamsum, epssum, sigmasum

  
end subroutine readState

subroutine addtwo(a,b, res)
  implicit none
  real, intent(in):: a,b
  real, intent(out) :: res
  
  res = a+b
end subroutine addtwo
subroutine initek
  !> 1.02 calculate elliptic integrals and store them in an array
    use comek
    implicit none
    integer :: j, elliptic_resolution
    real :: zarg
    real :: ze, zk, zt
!    real :: ef(:), ek(:)

    zt = 1.0 / float (max_num_of_elliptic_integral_points-2)
    do j = 1, max_num_of_elliptic_integral_points - 2
      zarg = float (j-1) * zt
      call ellipf(zarg, ze, zk)
      ef (j) = ze
      kf (j) = zk
    end do

  !> fixing the last two items in kf and ef
  !> \todo why do we need this??
    ef (max_num_of_elliptic_integral_points-1) = 1.0
    kf (max_num_of_elliptic_integral_points-1) = 1.0e20
    ef (max_num_of_elliptic_integral_points) = 1.0
    kf (max_num_of_elliptic_integral_points) = 1.0e20
    return
end subroutine initek
! - - ARE ALL THESE FUNCTIONS NECESSARY or just REDUNDANT CRAP?
!> 2.04 calculate the bounce integrals i1 and i2
subroutine ifunc (peps, plam, pi1, pi2)
!    use numeric_constants ! using only pi!!!
!    use globals
    implicit none
    real ::  pi = 3.1416; ! constant
    real, intent(in) :: peps, plam
    real, intent(out) :: pi1, pi2
    real :: zarg, ze, zk, zsqr
    
    ! trapped particles
    if (plam .ge. 1.0) then
      zarg = (1.0+plam*(2.0*peps-1.0)) / (2.0*peps*plam)
      if (zarg .lt. 0.0) zarg = 0.0

      call ek(zarg, ze, zk)

      zsqr = sqrt (2.0*peps*plam)
      pi1 = 2.0 * zk / (pi*zsqr)
      pi2 = 2.0 * zsqr * (ze-(1.0-zarg)*zk) / pi

    ! passing particles
    else
      zarg = 2.0 * peps * plam / (1.0+plam*(2.0*peps-1.0))
      if (zarg .lt. 0.0) zarg = 0.0

      call ek (zarg, ze, zk)
      zsqr = sqrt (1.0-plam+2.0*peps*plam)

      pi1 = 2.0 / (pi*zsqr) * zk
      pi2 = 2.0 * zsqr / pi * ze
    end if
    ! note: E and K are the full elliptic integrals that enter into i1, i2 expressions.
!    print *, "DBG: E: ",ze,", K: ",zk
    return
end subroutine ifunc
subroutine ek (parg, pe, pk)
    use comek
    implicit none
    integer :: i1, i2
    real, intent(in) :: parg
    real, intent(out) :: pe, pk
    real :: zarg1, zt
    
    real :: maxpt !max_num_of_elliptic_integral_points number of points at which to evaluate the integral(s)
    ! This SHOULD be read from comek??
    maxpt = max_num_of_elliptic_integral_points
!    print *, "DBG: maxpt =",maxpt
    
    zt = 1.0 / float (max_num_of_elliptic_integral_points-1)
    i1 = int ((max_num_of_elliptic_integral_points-1)*parg+0.000001) + 1
    i1 = min  (max_num_of_elliptic_integral_points-1, i1)
    i2 = i1 + 1

    zarg1 = float (i1-1) * zt

    zt = (parg-zarg1) * zt

   !> \todo why do we have item 501 and 502?
    pe = ef (i1) + (ef(i2)-ef(i1)) * zt
    pk = kf (i1) + (kf(i2)-kf(i1)) * zt

    return
end subroutine ek

end module cl_collop