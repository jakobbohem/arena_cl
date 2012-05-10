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

program sum
  use cl
!  use system_clock

  implicit none

  type(cl_platform_id)   :: platform
  type(cl_device_id)     :: device
  type(cl_context)       :: context
  type(cl_command_queue) :: command_queue
  type(cl_program)       :: prog
  type(cl_kernel)        :: kernel

  integer    :: num, ierr, irec, size
  integer(8) :: size_in_bytes, globalsize, localsize
  character(len = 100)  :: info, filename
  integer, parameter :: iunit = 10
  integer, parameter :: source_length = 10000 ! why is this exreplaceplicitly set?
  character(len = source_length) :: source
  
  ! remember to add the CL_variables here:
  real, allocatable  :: vec1(:), vec2(:), vec2_old_output(:), fsum(:), vec3(:), fsum_old_output(:)
  type(cl_mem)       :: cl_vec1, cl_vec2, cl_multip
  
  real :: sum1, sum2, sumsum, r, lightspeed
  integer :: tick_time, tock_time, ri
  integer :: k, ivar! loop variables

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
!  filename = 'square_add_one.cl'
  open(unit = iunit, file = 'sum.cl', access='direct', status = 'old', action = 'read', iostat = ierr, recl = 1)
  if (ierr /= 0) stop 'Cannot open the .cl program file!!'

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

  if(ierr /= CL_SUCCESS) stop 'Error: program build failed.'

  ! finally get the kernel and release the program
  kernel = clCreateKernel(prog, 'sum', ierr) ! note, the name here needs to match the call to the C-function
  call clReleaseProgram(prog, ierr)

  !=====================
  ! RUN THE KERNEL
  !=====================
  
  size = 5000
  lightspeed = 3e8
  size_in_bytes = int(size, 8)*4_8
  allocate(vec1(1:size))
  allocate(vec2(1:size))
  allocate(vec2_old_output(1:size))
  allocate(vec3(1:size))
  
  allocate(fsum(1:size))
  
  vec1 = 1.5 ! this is funky fortran initialisation of vectors!!
  vec2 = 2.0
  vec2_old_output = vec2; ! saved for the old loop comparison since vec2 is overwritten by the CLprogam
!  vec3 = 2.0; ! this is funky fortran initialisation of vectors!!
  
  ! can I pass in an integer or real as a parameter?
  
  ! allocate device memory
  cl_vec1 = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, ierr)
  cl_vec2 = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
  cl_multip = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)

  ! copy data to device memory
  call clEnqueueWriteBuffer(command_queue, cl_vec1, cl_bool(.true.), 0_8, size_in_bytes, vec1(1), ierr)
  call clEnqueueWriteBuffer(command_queue, cl_vec2, cl_bool(.true.), 0_8, size_in_bytes, vec2(1), ierr)
  call clEnqueueWriteBuffer(command_queue, cl_multip, cl_bool(.true.), 0_8, size_in_bytes, vec3(1), ierr);
  
  ! set the kernel arguments - REQUIRED FOR ALL PARAMETERS
  call clSetKernelArg(kernel, 0, size, ierr)
  call clSetKernelArg(kernel, 1, cl_vec1, ierr)
  call clSetKernelArg(kernel, 2, cl_vec2, ierr)
  call clSetKernelArg(kernel, 3, cl_multip, ierr)
  call clSetKernelArg(kernel, 4, lightspeed, ierr)

  ! get the localsize for the kernel (note that the sizes are integer(8) variable)
  call clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, localsize, ierr)
  globalsize = int(size, 8)
  
  print *, "DBG: modulus of global size and local size: " ,globalsize, localsize, mod(globalsize, localsize)
  if(mod(globalsize, localsize) /= 0) globalsize = globalsize + localsize - mod(globalsize, localsize) 

  ! execute the kernel
  call clEnqueueNDRangeKernel(command_queue, kernel, (/globalsize/), (/localsize/), ierr)
  call clFinish(command_queue, ierr)

  ! read the resulting vector from device memory
  call clEnqueueReadBuffer(command_queue, cl_vec2, cl_bool(.true.), 0_8, size_in_bytes, fsum(1), ierr)

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
  
 ! - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! OLD SCHOOL LOOP FOR COMPARISON:
 ! - - - - - - - - - - - - - - - - - - - - - - - - - -
 allocate(fsum_old_output(1:size))
  call system_clock(tick_time)
  do ivar = 0,size
      fsum_old_output(ivar) = sqrt(vec1(ivar)+vec2_old_output(ivar))*lightspeed
  end do
  call system_clock(tock_time) ! doesn't this measure milliseconds?
  
 ! - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! SUM THE OUTPUT FOR EASY DISPLAY:
 ! - - - - - - - - - - - - - - - - - - - - - - - - - -
  do ivar = 0,size
!    fsum(ivar) = 0
!    if (ivar .ne. 0) then
!      fsum(ivar) = fsum(ivar-1) + vec2(ivar)
!    end if
    sum1 = sum1+vec1(ivar)
    sum2 = sum2+vec2(ivar)
    sumsum = sumsum + fsum(ivar)
    
  end do
  
  
  
  call init_random_seed()
  call random_number(r)
  ri = nint(r*size)
  
  ! Print the output:
  print *, 'sqrt(sum( vec1 (',sum1,') and vec2 (',sum2,') ) is :', sumsum
  print *, 'output of random sum (at ',ri,') ' ,vec1(ri),' plus ',vec2(ri),' is ', fsum(ri)
  print *, 'old minus new: ',fsum_old_output(ri) - fsum(ri) 
    
  print *, "old loop time: ",tock_time-tick_time
  
  
end program sum

SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE