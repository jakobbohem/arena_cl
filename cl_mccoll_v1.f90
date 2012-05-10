!     
! File:   cl_mccoll_v1.f90
! Author: jakob
!
! Created on 02 May 2012, 14:41
!

MODULE cl_collision_operator
subroutine collisionOperator
  ! each call will run a collision operator calculation, this could be restructured again and again..
  use cl
  use comek
  use collision_operator ! for in line loop compare
!  use system_clock

  implicit none

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
  integer :: tick_time, tock_time, ri

  !! TEST SUBROUTINES
  test1 = 1
  test2 = 2
  call addtwo(test1,test2,itemp1);
  if (itemp1 .ne. 3.0) then
    print *, "sum of ",test1," and ",test2," is: ",itemp1
    print *, "ERROR in function call!!"
  end if
  
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
  call readState("arena_initial/pl4.res",pp,plam,peps,sigma)
  
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
  
  print *, "peps: ",peps(0)," plam: ", plam(0)," I1: ", I1(0),"I2: ", I2(0)
  print *, "peps: ",peps(2)," plam: ", plam(2)," I1: ", I1(2),"I2: ", I2(2)
  
  ! OLD SCHOOL LOOP FOR GETTING THE IFUNC:S (i1 and i2 elliptic integral)
  ! this should also be parallelised
  do ivar = 1,size
!    print *, "vec1(i): ",vec1(ivar)," vec2(i): ",vec2(ivar)
!    print *, "peps: ",peps(ivar)," plam: ", plam(ivar)," I1: ", I1(ivar),"I2: ", I2(ivar)
    call ifunc(peps(ivar), plam(ivar) ,itemp1,itemp2);
!    print *, 'called ifunc'
  end do
  
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
  

  !=====================
  ! RELEASE EVERYTHING
  !=====================

  call clReleaseKernel(kernel, ierr)
  call clReleaseCommandQueue(command_queue, ierr)
  call clReleaseContext(context, ierr)

  
  ! Slow old school fortran loop: (time this and compare)
  
  
  sum1 = 0
  sum2 = 0
  sumsum = 0
  

  call system_clock(tick_time)
  old_loop_count = 0
  allocate(pcdot_output(1:size))
  allocate(lcdot_output(1:size))
  allocate(ac_p_output(1:size))
  allocate(ac_lam_output(1:size))
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
  print *, "called old_loop ",old_loop_count," times."
  
  call system_clock(tock_time)
  
  call random_number(r)
  ri = nint(r*size)
  
  ! Print the output:
  print *, 'first in coll.op;  ppcdot(',ppcdot(1),'), pacp (',pacp(1),'), plcdot (',plcdot(1),'), paclam (',paclam(1),')'
  print *, 'random i in coll.op; ppcdot(',ppcdot(ri),'), pacp (',pacp(ri),'), plcdot (',plcdot(ri),'), paclam (',paclam(ri),')'
!  print *, 'DBG OLD CHECK: (at ',ri,') ' ,vec1(ri),' plus ',vec2(ri),' is ', fsum(ri)  
    
  print *, "old loop time: ",tock_time-tick_time, " (called ",old_loop_count," times)"
  print *, 'same rand i in old loop; ppcdot(',ppcdot(ri),'), pacp (',pacp(ri),'), plcdot (',plcdot(ri),'), paclam (',paclam(ri),')'

  
  !=====================
  ! END OF PROGRAM 
  !===================== 
end subroutine cl_mccoll_v1

END MODULE cl_collision_operator

