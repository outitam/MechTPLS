program mainprogram

  use tpls_constants
  use tpls_mpi
  use tpls_configuration
  use tpls_pressure_solver
  use tpls_userchk
  use tpls_io
  use tpls_solver

  implicit none
  include 'mpif.h'

  !----- To be determined -----!

  double precision, allocatable, dimension(:,:,:) :: vx, vy
  double precision, allocatable, dimension(:,:,:) :: vz
  double precision, allocatable, dimension(:,:,:) :: pres
  double precision, allocatable, dimension(:,:,:) :: phi

  
  !----- Miscellaneous -----!
  
  integer :: i, j, k

  double precision :: pi
  parameter ( pi = 4.0D+00*atan(1.0D+00) )

  double precision :: cfl
  double precision, dimension(1000000) :: dummy

  character*80 filename

  ! ****************************************************************************************
  ! Parameters

  call read_parameters()

  call initialize_tpls_mpi
  
  if ( my_id == 0 ) then
     call SRJ_weights(dummy)
  endif

  call mpi_bcast(max_srj, 1, mpi_integer, 0, comm2d_quasiperiodic, ierr)
  allocate(relax_srj(1:max_srj))
  relax_srj = 0

  call mpi_bcast(dummy, 100000, mpi_double_precision, 0, comm2d_quasiperiodic, ierr)

  relax_srj = dummy(1:max_srj)
  
  total_time    = mpi_wtime()
  time_levelset = 0.0D+00
  time_fluid    = 0.0D+00
  time_pressure = 0.0D+00
  t_temp        = 0.0D+00
    
  ! ****************************** Allocate all variables.  ******************************
      
  allocate(pres(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p))
  allocate(phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p))
  allocate(vx(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
  allocate(vy(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
  allocate(vz(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w))
  
  call logfile_header
  call dataload_channel_binary(vx,vy,vz,pres,phi)

  if ( my_id == master_id ) then
     write(*,*)
     write(*,*) 'Calling user_postprocess.'
     write(*,*)
  endif

  call user_postprocess(phi,pres,vx,vy,vz)

  call exchange2d(phi,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
  call exchange2d(vx,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
  call exchange2d(vy,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
  call exchange2d(vz, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,comm2d_quasiperiodic)
  call exchange2d(pres,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

  istep = 0
  call outpost_exact(vx,vy,vz,pres,phi,'   ')

  if ( my_id == master_id ) then
     write(*,*) '!----- Beginning of the time loop -----!'
     write(*,*)
  endif

  if (.not. if_linear_stability ) then

     call tpls_dns(vx,vy,vz,pres,phi)

  else

     call linear_stability(vx,vy,vz,pres,phi)

  endif
  
  total_time = mpi_wtime() - total_time  
  
  if ( my_id==master_id ) then    
     write(*,*) 'Runtime (excluding initialisation) :', total_time
     write(*,*) 'Total time in the levelset solver  :', time_levelset
     write(*,*) 'Total time in the fluid solver     :', time_fluid
     write(*,*) 'Total time in the pressure solver  :', time_pressure
     write(*,*) 'Finalising and deallocating.'
  end if
  
  deallocate(vx,vy,vz,pres,phi)
  
  call mpi_comm_free (comm2d_quasiperiodic, ierr)
  call mpi_finalize(ierr)
  
end program mainprogram
