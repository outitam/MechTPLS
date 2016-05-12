module tpls_io

  implicit none

contains

  ! ****************************************************************************************
  subroutine logfile_header

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    integer :: date_time(8)
    character*10 :: time_aux(3)

    if ( my_id == master_id ) then
       
       !----- RUN DETAILS -----!
       
       write(*,*) "-----     JOB DETAILS     -----"
       write(*,*)
       write(*,*) "Authors : Jean-Christophe Loiseau (loiseau@mech.kth.se)"
       write(*,*) "        & Outi Tammisola          (outi@mech.kth.se)"
       write(*,*) "Institution : KTH Mechanics, Stockholm, Sweden"
       write(*,*)
       
       call date_and_time(time_aux(1),time_aux(2),time_aux(3),date_time)
       
       write(*,'(4x,a,i2.2,a,i2.2,a,i4.4)') "Date : ", date_time(3), "/", date_time(2), "/", date_time(1)
       write(*,'(4x,a,i2.2,a,i2.2,a,i2.2)') "Time : ", date_time(5), ":", date_time(6), ":", date_time(7)
       
       write(*,*)
       write(*,*) 'Total number of processors         : ', num_procs
       write(*,*) 'Number of procs in the x-direction : ', num_procs_x
       write(*,*) 'Number of procs in the y-direction : ', num_procs_y
       write(*,*) 'Number of procs in the z-direction : ', num_procs_z
       write(*,*)
       
       if (num_procs_x*num_procs_y*num_procs_z /= num_procs) then
          write(*,*) 'Error 1: domain decomposition inconsistent with number of processors, exiting'
          stop
          
       end if
       
       if ( num_procs_z /= 1 ) then
          write(*,*) 'Error 2: domain decomposition inconsistent with number of processors, exiting'
          stop
       end if
       
       if ( (mod(maxl-1,num_procs_x) /= 0).or.(mod(maxm-1,num_procs_y) /= 0) ) then
          write(*,*) 'Error 3: number of processors must evenly divide (grid dimension-1), exiting'
          stop
       end if
       
       
       write(*,*)
       write(*,*) '*******************************************'
       write(*,*) '** For channel BCs: ex_max=', ex_max
       write(*,*) '*******************************************'
       
       write(*,*)
       write(*,*) '----- Physical parameters -----'
       write(*,*)
       write(*,*) 'Reynolds number              : ', Re
       write(*,*) 'Inverse Capillary number     : ', Scap
       write(*,*) 'Density of the upper layer   : ', rho_plus
       write(*,*) 'Viscosity of the upper layer : ', mu_plus
       write(*,*) 'Density of the lower layer   : ', rho_minus
       write(*,*) 'Viscosity of the lower layer : ', mu_minus
       write(*,*) 'Height of the liquid film    : ', height
       write(*,*) 'Gravitational number         : ', Grav
       write(*,*)
       
       write(*,*) '----- Numerical parameters -----'
       write(*,*)
       write(*,*) 'Time step                           : ', dt
       write(*,*) 'Number of time steps                : ', nsteps
       write(*,*) 'Number of points in x (streamwise)  : ', maxl-1
       write(*,*) 'Number of points in y (spanwise)    : ', maxm-1
       write(*,*) 'Number of points in z (wall-normal) : ', maxn-1
       write(*,*) 'Total number of points              : ', (maxl-1)*(maxm-1)*(maxn-1)
       write(*,*) 'Grid sizes (dx,dy,dz)               : ', dx, dy, dz
       write(*,*) 'Smooth width for levelset           : ', smooth_width/dx
       write(*,*) 'Tolerance of the levelset solver    : ', tolerance_levelset
       write(*,*) 'Tolerance of the Helmholtz solver   : ', tolerance_helmholtz
       write(*,*) 'Tolerance of the Poisson solver     : ', tolerance_poisson
       write(*,*)
       
       write(*,*) '----- Logicals -----'
       write(*,*)
       write(*,*) 'Density match                     : ', density_matched
       write(*,*) 'Viscosity match                   : ', viscosity_matched
       write(*,*) 'Implicit time integration         : ', if_implicit
       write(*,*) 'Restart from an existing file     : ', if_restart
       write(*,*) '                   Exact restart? : ', if_exact_restart
       write(*,*) 'Save files for an exact restart?  : ', if_exact_restart_save
       write(*,*) 'Selective Frequency Damping       : ', if_sfd
       write(*,*) '         Restart SFD from a file? : ', if_exact_restart_sfd
       write(*,*) 'Linear stability mode             : ', if_linear_stability
       write(*,*)

       if ( if_restart ) then

          write(*,*) '----- Restarting procedure -----'
          write(*,*)
          write(*,*) 'Restart basic fields from file : ', restart_handle
          write(*,*)
          if (if_exact_restart) then
               write(*,*) 'Exact restart without initial Euler step' 
               write(*,*) 'Additional restart file: fieldexactrestart_conv.bin' 
               write(*,*)
          endif
          if (if_exact_restart_sfd) then
               write(*,*) 'Exact restart of SFD variables' 
               write(*,*) 'Additional restart file: fieldexactrestart_sfd.bin'
               write(*,*)
          endif 
       endif

       if ( if_linear_stability ) then

          write(*,*) '----- Linear stability setup -----'
          write(*,*)
          write(*,*) 'Dimension of the Krylov subspace : ', krylov_dim
          write(*,*) 'Baseflow restarting file         : ', restart_handle
          write(*,*)
          
       endif
       
       write(*,*) '----- Initialization -----'
       write(*,*)
       write(*,*) 'Calling user_initial_conditions.'

    endif
       
    return
  end subroutine logfile_header

  ! ****************************************************************************************
  subroutine outpost_exact(vx,vy,vz,pr,phi,prefix)
!   This creates a file in exactly the same format as outpost.
!   However, no BC:s are imposed. All variables are stored including 
!   their boundary values at the time of call of outpost_exact.
    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w), intent(in) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(in) :: pr, phi
    character*3 :: prefix

    !----- Miscellaneous -----!

    double precision, dimension(0:maxl,0:maxm,0:maxn-2) :: u_global, v_global
    double precision, dimension(0:maxl,0:maxm,0:maxn-1) :: w_global
    double precision, dimension(0:maxl,0:maxm,0:maxn)   :: pr_global, phi_global

    integer :: i, j, k, msg_length
    character*80 filename

    !----- Gather all the global variables -----!

    u_global   = 0.0D+00
    v_global   = 0.0D+00
    w_global   = 0.0D+00
    pr_global  = 0.0D+00
    phi_global = 0.0D+00

    if (num_procs==1) then
       u_global = vx
       v_global = vy
       w_global = vz
       pr_global = pr
       phi_global = phi
    else
       call gather_tables_exact(u_global,v_global,w_global,pr_global,phi_global,vx,vy,vz,pr,phi)
    endif

    !----- Outpost the resulting arrays -----!
    
    if ( my_id == master_id ) then
       
       write(*,*) 'Iteration : ', istep, 'Writing to file on master node'
       
       if(istep.lt.10) then
          write(filename,'(A,I1,A)')'channel_00000000',istep,'.vtk'
       elseif((istep.GE.10).and.(istep.lt.100)) then
          write(filename,'(A,I2,A)')'channel_0000000',istep,'.vtk'
       elseif((istep.GE.100).and.(istep.lt.1000)) then
          write(filename,'(A,I3,A)')'channel_000000',istep,'.vtk'
       elseif((istep.GE.1000).and.(istep.lt.10000)) then
          write(filename,'(A,I4,A)')'channel_00000',istep,'.vtk'
       elseif((istep.GE.10000).and.(istep.lt.100000)) then
          write(filename,'(A,I5,A)')'channel_0000',istep,'.vtk'
       elseif((istep.GE.100000).and.(istep.lt.1000000)) then
          write(filename,'(A,I6,A)')'channel_000',istep,'.vtk'
       endif

       if ( prefix /= '   ' ) filename = prefix // filename

       call vtk_output(u_global, v_global, w_global, pr_global, phi_global, filename)
       
       if ( istep.GT.0 ) call backup_channel_binary(u_global, v_global, w_global, pr_global, phi_global)
       
 
    end if
    
    return
  end subroutine outpost_exact
  
  ! ****************************************************************************************
  subroutine outpost(vx,vy,vz,pr,phi,prefix)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w), intent(in) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(in) :: pr, phi
    character*3 :: prefix

    !----- Miscellaneous -----!

    double precision, dimension(0:maxl,0:maxm,0:maxn-2) :: u_global, v_global
    double precision, dimension(0:maxl,0:maxm,0:maxn-1) :: w_global
    double precision, dimension(0:maxl,0:maxm,0:maxn)   :: pr_global, phi_global

    integer :: i, j, k, msg_length
    character*80 filename

    !----- Gather all the global variables -----!

    u_global   = 0.0D+00
    v_global   = 0.0D+00
    w_global   = 0.0D+00
    pr_global  = 0.0D+00
    phi_global = 0.0D+00

    call gather_tables(u_global,v_global,w_global,pr_global,phi_global,vx,vy,vz,pr,phi)

    !----- Outpost the resulting arrays -----!
    
    if ( my_id == master_id ) then
       
       write(*,*) 'Iteration : ', istep, 'Writing to file on master node'
       
       if(istep.lt.10) then
          write(filename,'(A,I1,A)')'channel_00000000',istep,'.vtk'
       elseif((istep.GE.10).and.(istep.lt.100)) then
          write(filename,'(A,I2,A)')'channel_0000000',istep,'.vtk'
       elseif((istep.GE.100).and.(istep.lt.1000)) then
          write(filename,'(A,I3,A)')'channel_000000',istep,'.vtk'
       elseif((istep.GE.1000).and.(istep.lt.10000)) then
          write(filename,'(A,I4,A)')'channel_00000',istep,'.vtk'
       elseif((istep.GE.10000).and.(istep.lt.100000)) then
          write(filename,'(A,I5,A)')'channel_0000',istep,'.vtk'
       elseif((istep.GE.100000).and.(istep.lt.1000000)) then
          write(filename,'(A,I6,A)')'channel_000',istep,'.vtk'
       endif

       if ( prefix /= '   ' ) filename = prefix // filename

       call vtk_output(u_global, v_global, w_global, pr_global, phi_global, filename)
       
       if ( istep.GT.0 ) call backup_channel_binary(u_global, v_global, w_global, pr_global, phi_global)
       
 
    end if
    
    return
  end subroutine outpost
  
  ! ****************************************************************************************
  subroutine vtk_output(vx,vy,vz,pr,phi,filename)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2), intent(in) :: vx, vy
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1), intent(in) :: vz
    double precision, dimension(0:maxl, 0:maxm, 0:maxn),   intent(in) :: phi, pr

    character(len=80), intent(in) :: filename

    !----- Miscellaneous -----!

    integer :: i, j, k, ifich
    integer :: nx, ny, nz
    integer :: iostatus
    character :: q
    character(len=1), parameter :: newline = char(10)
    character(len=100) :: s_buffer

    integer :: output_unit, input_unit
    parameter ( input_unit = 8, output_unit = 9 )

    !----- Writing to file -----!

    nx = maxl-1
    ny = maxm-1
    nz = maxn-1

    ifich = 10
    q = char(34)

    if ( maxm.GT.2 ) then

       open( unit = ifich , file = filename , form = 'unformatted' , access = 'stream' , &
            action = 'write' , convert = 'BIG_ENDIAN' , iostat = iostatus )

       write(unit = ifich, iostat = iostatus) '# vtk DataFile Version 3.0' // newline
       write(unit = ifich, iostat = iostatus) 'test file' // newline
       write(unit = ifich, iostat = iostatus) 'BINARY' // newline
       write(unit = ifich, iostat = iostatus) newline
       write(unit = ifich, iostat = iostatus) 'DATASET ' // 'STRUCTURED_POINTS' // newline

       write(s_buffer, FMT = '(A,3I12)', iostat = iostatus) 'DIMENSIONS', nx, ny, nz
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'ORIGIN', dx/2.0D+00, dy/2.0D+00, dz/2.0D+00
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'SPACING', dx, dy, dz
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

       write(unit = ifich, iostat = iostatus) newline
       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'POINT_DATA', nx*ny*nz
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS U double', 1
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

       do k = 0, nz-1
          do j = 0, ny-1
             do i = 1,nx

                write(unit = ifich, iostat = iostatus) (vx(i,j,k) + vx(i+1,j,k))/2.0D+00

             enddo
          enddo
       enddo

       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS V double', 1
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

       do k = 0, nz-1
          do j = 0, ny-1
             do i = 1,nx

                write(unit = ifich, iostat = iostatus) (vy(i,j,k) + vy(i,j+1,k))/2.0D+00

             enddo
          enddo
       enddo


       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS W double', 1
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

       do k = 0, nz-1
          do j = 0, ny-1
             do i = 1,nx

                write(unit = ifich, iostat = iostatus) (vz(i,j,k) + vz(i,j,k+1))/2.0D+00

             enddo
          enddo
       enddo

       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS Pressure double', 1
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

       do k = 0, nz-1
          do j = 0, ny-1
             do i = 1,nx

                write(unit = ifich, iostat = iostatus) pr(i,j,k+1)

             enddo
          enddo
       enddo

       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS Phi double', 1
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

       do k = 0, nz-1
          do j = 0, ny-1
             do i = 1,nx

                write(unit = ifich, iostat = iostatus) phi(i,j,k+1)

             enddo
          enddo
       enddo

       write(unit = ifich, iostat = iostatus) newline
       close(ifich)

    elseif ( maxm .EQ. 2) then

       ny = 1

       open( unit = ifich , file = filename , form = 'unformatted' , access = 'stream' , &
            action = 'write' , convert = 'BIG_ENDIAN' , iostat = iostatus )

       write(unit = ifich, iostat = iostatus) '# vtk DataFile Version 3.0' // newline
       write(unit = ifich, iostat = iostatus) 'test file' // newline
       write(unit = ifich, iostat = iostatus) 'BINARY' // newline
       write(unit = ifich, iostat = iostatus) newline
       write(unit = ifich, iostat = iostatus) 'DATASET ' // 'STRUCTURED_POINTS' // newline

       write(s_buffer, FMT = '(A,3I12)', iostat = iostatus) 'DIMENSIONS', nx, nz, ny
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'ORIGIN', dx/2.0D+00, dz/2.0D+00, 0.0D+00
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'SPACING', dx, dy, dz
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

       write(unit = ifich, iostat = iostatus) newline
       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'POINT_DATA', nx*ny*nz
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS U double', 1
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

       do k = 0,nz-1
          do j = 0,ny-1
             do i = 1,nx

                write(unit = ifich, iostat = iostatus) (vx(i,j,k) + vx(i+1,j,k))/2.0D+00
                
             enddo
          enddo
       enddo

       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS W double', 1
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

       do k = 0, nz-1
          do j = 0, ny-1
             do i = 1,nx

                write(unit = ifich, iostat = iostatus) (vz(i,j,k) + vz(i,j,k+1))/2.0D+00
                
             enddo
          enddo
       enddo

       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS Pressure double', 1
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

       do k = 1, nz
          do j = 0, ny-1
             do i = 1,nx

                write(unit = ifich, iostat = iostatus) pr(i,j,k)

             enddo
          enddo
       enddo

       write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS Phi double', 1
       write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
       write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

       do k = 1, nz
          do j = 0, ny-1
             do i = 1,nx

                write(unit = ifich, iostat = iostatus) phi(i,j,k)

             enddo
          enddo
       enddo

       write(unit = ifich, iostat = iostatus) newline
       close(ifich)


    endif
    
    return
  end subroutine vtk_output
  
  ! ****************************************************************************************

  subroutine backup_channel_binary(vx,vy,vz,pr,phi)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2), intent(in) :: vx, vy
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1), intent(in) :: vz
    double precision, dimension(0:maxl, 0:maxm, 0:maxn),   intent(in) :: phi, pr

    !----- Miscellaneous -----!

    character(len=100) :: filename
    integer :: i, j, k, ifich

    !----- Output data for possible restart -----!

    ifich = 10

    filename = 'fieldrestart_backup.bin'
    open( unit = ifich, file = filename, form = 'unformatted')

    write(ifich) vx
    write(ifich) vy
    write(ifich) vz
    write(ifich) pr
    write(ifich) phi

    close(ifich)

    return
  end subroutine backup_channel_binary
  
  ! ****************************************************************************************
  subroutine  backup_exactrest_sfd(vx_bar,vy_bar,vz_bar,fx_sfd,fy_sfd,fz_sfd)
    
    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx_bar, vy_bar, fx_sfd, fy_sfd
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(in) :: vz_bar, fz_sfd

    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2) :: vx_bar_gl, vy_bar_gl, fx_sfd_gl, fy_sfd_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1) :: vz_bar_gl, fz_sfd_gl

    !----- Miscellaneous -----!

    character(len=100) :: filename
    integer :: i, j, k, ifich

    call gather_tables_exactsfd(vx_bar_gl,vy_bar_gl,vz_bar_gl,fx_sfd_gl,fy_sfd_gl,fz_sfd_gl, &
                                vx_bar,vy_bar,vz_bar,fx_sfd,fy_sfd,fz_sfd)
    
    !----- Output data for possible restart -----!


    if (my_id==0) then
       
        write(*,*) 'Iteration : ', istep, 'Writing sfd variables to file on master node'
    ifich = 101

    filename = 'fieldexactrestart_sfd.bin'

    open( unit = ifich, file = filename, form = 'unformatted')

    write(ifich) vx_bar_gl
    write(ifich) vy_bar_gl
    write(ifich) vz_bar_gl
    write(ifich) fx_sfd_gl
    write(ifich) fy_sfd_gl
    write(ifich) fz_sfd_gl

    close(ifich)
    endif

    return
  end subroutine backup_exactrest_sfd
  
  ! ****************************************************************************************
  subroutine backup_exactrest_conv(conv1_u,conv1_v,conv1_w,csf_u1,csf_v1,csf_w1,   &
                              diffusion1_u,diffusion1_v,diffusion1_w)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: conv1_u, conv1_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: csf_u1, csf_v1
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: diffusion1_u, diffusion1_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(in) :: conv1_w
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(in) :: csf_w1
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(in) :: diffusion1_w

    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2) :: conv1_u_gl, conv1_v_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2) :: csf_u1_gl, csf_v1_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2) :: diffusion1_u_gl, diffusion1_v_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1) :: conv1_w_gl, csf_w1_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1) :: diffusion1_w_gl

    !----- Miscellaneous -----!

    character(len=100) :: filename
    integer :: i, j, k, ifich

    call gather_tables_exactconv(conv1_u,conv1_v,conv1_w,csf_u1,csf_v1,csf_w1,  &
            diffusion1_u,diffusion1_v,diffusion1_w,  &
             conv1_u_gl,conv1_v_gl,conv1_w_gl,csf_u1_gl,csf_v1_gl,csf_w1_gl, &
            diffusion1_u_gl,diffusion1_v_gl,diffusion1_w_gl)
    !----- Output data for possible restart -----!


    if (my_id==0) then
    ifich = 10

    filename = 'fieldexactrestart_conv.bin'

    open( unit = ifich, file = filename, form = 'unformatted')

    write(ifich) conv1_u_gl
    write(ifich) csf_u1_gl
    write(ifich) diffusion1_u_gl
    write(ifich) conv1_v_gl
    write(ifich) csf_v1_gl
    write(ifich) diffusion1_v_gl
    write(ifich) conv1_w_gl
    write(ifich) csf_w1_gl
    write(ifich) diffusion1_w_gl

    close(ifich)
    endif

    return
  end subroutine backup_exactrest_conv
  
  ! ****************************************************************************************

  subroutine gather_tables(u_global,v_global,w_global,pres_global,phi_global,vx,vy,vz,pr,phi)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(in) :: pr, phi

    !----- Outputs -----!

    double precision, dimension(0:maxl,0:maxm,0:maxn-2), intent(out) :: u_global, v_global
    double precision, dimension(0:maxl,0:maxm,0:maxn-1), intent(out) :: w_global
    double precision, dimension(0:maxl,0:maxm,0:maxn)  , intent(out) :: pres_global, phi_global

    !----- Miscellaneous -----!

    integer :: i, j, k
    double precision, dimension(1:n_local_x,1:n_local_y,1:n_local_z_uv) :: tempr_u, tempr_v
    double precision, dimension(1:n_local_x,1:n_local_y,1:n_local_z_w)  :: tempr_w
    double precision, dimension(1:n_local_x,1:n_local_y,1:n_local_z_p)  :: tempr_p, tempr_phi

    !----- MPI related variables -----!
    
    integer :: status(MPI_STATUS_SIZE)

    u_global    = 0.0D+00
    v_global    = 0.0D+00
    w_global    = 0.0D+00
    pres_global = 0.0D+00
    phi_global  = 0.0D+00

!!!!! First: impose the same BC as after velocity correction

    if ( my_id == master_id ) then          
       
       pz = 0

       do py = 0,dims(2)-1

          do px = 0,dims(1)-1
             
             if( (pz==0).and.(py==0).and.(px==0) )then
                u_global(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv) = vx(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)    
                v_global(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv) = vy(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)    
             else
                
                coords(1) = px
                coords(2) = py
                coords(3) = pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                Call mpi_recv(sx_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_u,n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_v,n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)

                
                do k = 1,n_local_z_uv
                   do j = 1,n_local_y
                      do i = 1,n_local_x
                         u_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_u(i,j,k)
                         v_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_v(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do

    else
       
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(vx(sx:ex,sy:ey,sz_uv:ez_uv),n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(vy(sx:ex,sy:ey,sz_uv:ez_uv),n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)


    endif

    ! enforcing  boundary conditions on ghost cells
!!  Changing to the same as after pressure correction!! 
    
    if (my_id==master_id) then

!!!    inflow u
       u_global(0,:,:) = 2.0D+00*u_global(1,:,:)-u_global(2,:,:);

!!!    no-slip  u
       u_global(:,:,sz_uv) = (1.0D+00/3.0D+00)*u_global(:,:,sz_uv+1);
       u_global(:,:,ez_uv) = (1.0D+00/3.0D+00)*u_global(:,:,ez_uv-1);
!!!    outflow  u
       u_global(maxl,:,:) = u_global(maxl-1,:,:);

       u_global(maxl-1,:,:) = u_global(maxl-2,:,:)
       u_global(maxl,:,:)   = u_global(maxl-1,:,:)
       
       u_global(:,0,:)    = u_global(:,maxm-1,:)
       u_global(:,maxm,:) = u_global(:,1,:)

       v_global(0,:,:)    = vy(0,:,:)
       v_global(maxl,:,:) = v_global(maxl-1,:,:)
       v_global(:,0,:)    = v_global(:,maxm-1,:)
       v_global(:,maxm,:) = v_global(:,1,:)
       
    end if


    If (my_id==master_id) Then

       pz=0

       do py=0,dims(2)-1
          do px=0,dims(1)-1
             
             if((pz==0).and.(py==0).and.(px==0))then
                w_global(sx_o:ex_o,sy_o:ey_o,sz_w:ez_w) = vz(sx_o:ex_o,sy_o:ey_o,sz_w:ez_w)    
             else
                
                coords(1)=px
                coords(2)=py
                coords(3)=pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_w,n_local_x*n_local_y*n_local_z_w,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                
                do k=1,n_local_z_w
                   do j=1,n_local_y
                      do i=1,n_local_x
                         w_global(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_w(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do
    Else
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(vz(sx:ex,sy:ey,sz_w:ez_w),n_local_x*n_local_y*n_local_z_w,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
    End If
    
    ! enforcing conditions on ghost cells
    
    if (my_id==master_id) then
       
!!$       w_global(:,:,0)      = vz(:,:,0)
!!$       w_global(:,:,maxn-1) = vz(:,:,maxn-1)
       
       w_global(0,:,:)    = vz(0,:,:)
       w_global(maxl,:,:) = w_global(maxl-1,:,:)
       
       w_global(:,0,:)    = w_global(:,maxm-1,:)
       w_global(:,maxm,:) = w_global(:,1,:)
       
    end if




    If (my_id==master_id) Then
       pz=0
       do py=0,dims(2)-1
          do px=0,dims(1)-1
             
             if((pz==0).and.(py==0).and.(px==0))then
                pres_global(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p) = pr(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)    
             else
                
                coords(1)=px
                coords(2)=py
                coords(3)=pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_p,n_local_x*n_local_y*n_local_z_p,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                
                do k=1,n_local_z_p
                   do j=1,n_local_y
                      do i=1,n_local_x
                         pres_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_p(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do
    Else
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(pr(sx:ex,sy:ey,sz_p:ez_p),n_local_x*n_local_y*n_local_z_p,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       
    End If
    
    ! enforcing boundary conditions on ghost cells
    
    if (my_id==master_id) then            
       pres_global(0,:,:)    = pr(0,:,:)
       pres_global(maxl,:,:) = pres_global(maxl-1,:,:)
       
       pres_global(:,0,:)=pres_global(:,maxm-1,:)
       pres_global(:,maxm,:)=pres_global(:,1,:)
       pres_global(:,:,maxn)=pres_global(:,:,maxn-1)
       pres_global(:,:,0)=pres_global(:,:,1)
    end if





    If (my_id==master_id) Then
       phi_global = 0.d0
       pz=0
       do py=0,dims(2)-1
          do px=0,dims(1)-1
             
             if((pz==0).and.(py==0).and.(px==0))then
                phi_global(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)=phi(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)    
             else
                
                coords(1)=px
                coords(2)=py
                coords(3)=pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_phi,n_local_x*n_local_y*n_local_z_p,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                
                do k=1,n_local_z_p
                   do j=1,n_local_y
                      do i=1,n_local_x
                         phi_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_phi(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do
    Else
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(phi(sx:ex,sy:ey,sz_p:ez_p),n_local_x*n_local_y*n_local_z_p,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       
    End If
    
    ! enforcing boundary conditions on ghost cells
    
    if (my_id==master_id) then
       
       phi_global(0,:,:) = phi(0,:,:)
       
       phi_global(maxl,:,:) = phi_global(maxl-1,:,:)
       
       phi_global(:,0,:) = phi_global(:,maxm-1,:)
       phi_global(:,maxm,:) = phi_global(:,1,:)
       phi_global(:,:,maxn) = phi_global(:,:,maxn-1)
       phi_global(:,:,0) = phi_global(:,:,1)

    end if
    
    
    return
  end subroutine gather_tables
  
  ! ****************************************************************************************

  subroutine gather_tables_exact(u_global,v_global,w_global,pres_global,phi_global,vx,vy,vz,pr,phi)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(in) :: pr, phi

    !----- Outputs -----!

    double precision, dimension(0:maxl,0:maxm,0:maxn-2), intent(out) :: u_global, v_global
    double precision, dimension(0:maxl,0:maxm,0:maxn-1), intent(out) :: w_global
    double precision, dimension(0:maxl,0:maxm,0:maxn)  , intent(out) :: pres_global, phi_global

    !----- Miscellaneous -----!

    integer :: i, j, k, istart, iend, jstart, jend, kstart, kend
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_uv) :: tempr_u, tempr_v
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_w)  :: tempr_w
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_p)  :: tempr_p, tempr_phi

    !----- MPI related variables -----!
    
    integer :: status(MPI_STATUS_SIZE)

    u_global    = 0.0D+00
    v_global    = 0.0D+00
    w_global    = 0.0D+00
    pres_global = 0.0D+00
    phi_global  = 0.0D+00


    if ( my_id == master_id ) then          
       
       pz = 0

       do py = 0,dims(2)-1

          do px = 0,dims(1)-1
             
             if( (pz==0).and.(py==0).and.(px==0) )then
!         Since this is the first processor in x and y, save from sx-1 & sy-1 
                u_global(sx-1:ex,sy-1:ey,sz_uv:ez_uv) = vx(sx-1:ex,sy-1:ey,sz_uv:ez_uv)    
                v_global(sx-1:ex,sy-1:ey,sz_uv:ez_uv) = vy(sx-1:ex,sy-1:ey,sz_uv:ez_uv)    
             else
                
                coords(1) = px
                coords(2) = py
                coords(3) = pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                Call mpi_recv(sx_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_u,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_v,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)


                jstart=1
                jend=n_local_y
                istart=1
                iend=n_local_x
                
!         If this is the first processor in x (or y), save from sx-1 (or sy-1) 
                if (px==0) istart=istart-1
                if (py==0) jstart=jstart-1
!         If this is the last processor in x (or y), save until ex+1 (or ey+1) 
                if (px==(dims(1)-1)) iend=iend+1
                if (py==(dims(2)-1)) jend=jend+1
                

                do k = 1,n_local_z_uv
                   do j = jstart,jend
                      do i = istart,iend
                         u_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_u(i,j,k)
                         v_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_v(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do

    else
       
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(vx(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(vy(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)


    endif



    If (my_id==master_id) Then

       pz=0

       do py=0,dims(2)-1
          do px=0,dims(1)-1
             
             if((pz==0).and.(py==0).and.(px==0))then
!         Since this is the first processor in x and y, save from sx-1 & sy-1 
                w_global(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w) = vz(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w)    
             else
                
                coords(1)=px
                coords(2)=py
                coords(3)=pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_w,(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                
                jstart=1
                jend=n_local_y
                istart=1
                iend=n_local_x
!         If this is the first processor in x (or y), save from sx-1 (or sy-1) 
                if (px==0) istart=istart-1
                if (py==0) jstart=jstart-1
!         If this is the last processor in x (or y), save until ex+1 (or ey+1) 
                if (px==(dims(1)-1)) iend=iend+1
                if (py==(dims(2)-1)) jend=jend+1

                do k=1,n_local_z_w
                   do j=jstart,jend
                      do i=istart,iend
                         w_global(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_w(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do
    Else
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(vz(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
    End If
    

    If (my_id==master_id) Then
       pz=0
       do py=0,dims(2)-1
          do px=0,dims(1)-1
             
             if((pz==0).and.(py==0).and.(px==0))then
!         Since this is the first processor in x and y, save from sx-1 & sy-1 
                pres_global(sx_o-1:ex_o,sy_o-1:ey_o,sz_p:ez_p) = pr(sx_o-1:ex_o,sy_o-1:ey_o,sz_p:ez_p)    
             else
                
                coords(1)=px
                coords(2)=py
                coords(3)=pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_p,(n_local_x+2)*(n_local_y+2)*n_local_z_p,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                
                jstart=1
                jend=n_local_y
                istart=1
                iend=n_local_x
!         If this is the first processor in x (or y), save from sx-1 (or sy-1) 
                if (px==0) istart=istart-1
                if (py==0) jstart=jstart-1
!         If this is the last processor in x (or y), save until ex+1 (or ey+1) 
                if (px==(dims(1)-1)) iend=iend+1
                if (py==(dims(2)-1)) jend=jend+1
                
                do k=1,n_local_z_p
                   do j=jstart,jend
                      do i=istart,iend
                         pres_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_p(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do
    Else
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(pr(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),(n_local_x+2)*(n_local_y+2)*n_local_z_p,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       
    End If
    

    If (my_id==master_id) Then
       phi_global = 0.d0
       pz=0
       do py=0,dims(2)-1
          do px=0,dims(1)-1
             
             if((pz==0).and.(py==0).and.(px==0))then
!         Since this is the first processor in x and y, save from sx-1 & sy-1 
                phi_global(sx_o-1:ex_o,sy_o-1:ey_o,sz_p:ez_p)=phi(sx_o-1:ex_o,sy_o-1:ey_o,sz_p:ez_p)    
             else
                
                coords(1)=px
                coords(2)=py
                coords(3)=pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_phi,(n_local_x+2)*(n_local_y+2)*n_local_z_p,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                
                jstart=1
                jend=n_local_y
                istart=1
                iend=n_local_x
!         If this is the first processor in x (or y), save from sx-1 (or sy-1) 
                if (px==0) istart=istart-1
                if (py==0) jstart=jstart-1
!         If this is the last processor in x (or y), save until ex+1 (or ey+1) 
                if (px==(dims(1)-1)) iend=iend+1
                if (py==(dims(2)-1)) jend=jend+1
                
                do k=1,n_local_z_p
                   do j=jstart,jend
                      do i=istart,iend
                         phi_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_phi(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do
    Else
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),(n_local_x+2)*(n_local_y+2)*n_local_z_p,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       
    End If
    
    
    
    return
  end subroutine gather_tables_exact
  
  ! ****************************************************************************************

  subroutine gather_tables_exactconv(conv1_u,conv1_v,conv1_w,csf_u1,csf_v1,csf_w1, &
             diffusion1_u,diffusion1_v,diffusion1_w,  &
              conv1_ugl,conv1_vgl,conv1_wgl,csf1_u1gl,csf1_v1gl,csf1_w1gl, &
             diffusion1_ugl,diffusion1_vgl,diffusion1_wgl) 

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: conv1_u, conv1_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: csf_u1, csf_v1
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: diffusion1_u, diffusion1_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(in) :: conv1_w
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(in) :: csf_w1
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(in) :: diffusion1_w

    !----- Outputs -----!

    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2), intent(out) :: conv1_ugl, conv1_vgl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2), intent(out) :: csf1_u1gl, csf1_v1gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2), intent(out) :: diffusion1_ugl, diffusion1_vgl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1), intent(out) :: conv1_wgl, csf1_w1gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1), intent(out) :: diffusion1_wgl

    !----- Miscellaneous -----!

    integer :: i, j, k, istart, iend, jstart, jend
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_uv) :: tempr_convu, tempr_convv
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_uv) :: tempr_csfu, tempr_csfv
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_uv) :: tempr_diffu, tempr_diffv
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_w)  :: tempr_convw
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_w)  :: tempr_csfw
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_w)  :: tempr_diffw

    !----- MPI related variables -----!
    
    integer :: status(MPI_STATUS_SIZE)

    conv1_ugl    = 0.0D+00
    conv1_vgl    = 0.0D+00
    conv1_wgl    = 0.0D+00
    csf1_u1gl    = 0.0D+00
    csf1_v1gl    = 0.0D+00
    csf1_w1gl    = 0.0D+00
    diffusion1_ugl    = 0.0D+00
    diffusion1_vgl    = 0.0D+00
    diffusion1_wgl    = 0.0D+00

    if ( my_id == master_id ) then          
       
       pz = 0

       do py = 0,dims(2)-1

          do px = 0,dims(1)-1
            
!         Since this is the first processor in x and y, save from sx-1 & sy-1 
             if( (pz==0).and.(py==0).and.(px==0) )then
                conv1_ugl(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv) = conv1_u(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv)    
                conv1_vgl(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv) = conv1_v(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv)    
                csf1_u1gl(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv) = csf_u1(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv)    
                csf1_v1gl(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv) = csf_v1(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv)    
                diffusion1_ugl(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv) = &
                 diffusion1_u(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv)    
                diffusion1_vgl(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv) = &
                 diffusion1_v(sx_o-1:ex_o,sy_o-1:ey_o,sz_uv:ez_uv)    
             else
                
                coords(1) = px
                coords(2) = py
                coords(3) = pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                Call mpi_recv(sx_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_convu,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_convv,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_csfu,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_csfv,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_diffu,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_diffv,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)

               
                jstart=1
                jend=n_local_y
                istart=1
                iend=n_local_x
                
!         If this is the first processor in x (or y), save from sx-1 (or sy-1) 
                if (px==0) istart=istart-1
                if (py==0) jstart=jstart-1
!         If this is the last processor in x (or y), save until ex+1 (or ey+1) 
                if (px==(dims(1)-1)) iend=iend+1
                if (py==(dims(2)-1)) jend=jend+1
 
                do k = 1,n_local_z_uv
                   do j = jstart,jend
                      do i = istart,iend
                         conv1_ugl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_convu(i,j,k)
                         conv1_vgl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_convv(i,j,k)
                         csf1_u1gl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_csfu(i,j,k)
                         csf1_v1gl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_csfv(i,j,k)
                         diffusion1_ugl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_diffu(i,j,k)
                         diffusion1_vgl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_diffv(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do

    else
       
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(conv1_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), &
            (n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(conv1_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), &
            (n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(csf_u1(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), &
            (n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(csf_v1(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), &
            (n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(diffusion1_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), &
            (n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(diffusion1_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), &
            (n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       


    endif

    If (my_id==master_id) Then

       pz=0

       do py=0,dims(2)-1
          do px=0,dims(1)-1
             
             if((pz==0).and.(py==0).and.(px==0))then
!         Since this is the first processor in x and y, save from sx-1 & sy-1 
                conv1_wgl(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w) = conv1_w(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w)    
                csf1_w1gl(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w) = csf_w1(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w)    
                diffusion1_wgl(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w) = &
                 diffusion1_w(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w)    
             else
                
                coords(1)=px
                coords(2)=py
                coords(3)=pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_convw,(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_csfw,(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_diffw,(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                

                  
                jstart=1
                jend=n_local_y
                istart=1
                iend=n_local_x
                
!         If this is the first processor in x (or y), save from sx-1 (or sy-1) 
                if (px==0) istart=istart-1
                if (py==0) jstart=jstart-1
!         If this is the last processor in x (or y), save until ex+1 (or ey+1) 
                if (px==dims(1)-1) iend=iend+1
                if (py==dims(2)-1) jend=jend+1

                do k=1,n_local_z_w
                   do j=jstart,jend
                      do i=istart,iend
                         conv1_wgl(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_convw(i,j,k)
                         csf1_w1gl(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_csfw(i,j,k)
                         diffusion1_wgl(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_diffw(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do
    Else
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
      
       Call mpi_ssend(conv1_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(csf_w1(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(diffusion1_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       

    End If
    

    
    return
  end subroutine gather_tables_exactconv
  
  ! ****************************************************************************************

  subroutine gather_tables_exactsfd(vx_bar_gl,vy_bar_gl,vz_bar_gl,fx_sfd_gl,fy_sfd_gl,fz_sfd_gl, &
                                vx_bar,vy_bar,vz_bar,fx_sfd,fy_sfd,fz_sfd)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx_bar, vy_bar, fx_sfd, fy_sfd
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: vz_bar, fz_sfd

    !----- Outputs -----!

    double precision, dimension(0:maxl,0:maxm,0:maxn-2), intent(out) :: vx_bar_gl, vy_bar_gl, fx_sfd_gl, fy_sfd_gl
    double precision, dimension(0:maxl,0:maxm,0:maxn-1), intent(out) :: vz_bar_gl, fz_sfd_gl

    !----- Miscellaneous -----!

    integer :: i, j, k, istart, iend, jstart, jend
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_uv) :: tempr_vx_bar, tempr_vy_bar
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_w)  :: tempr_vz_bar
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_uv) :: tempr_fx_sfd, tempr_fy_sfd
    double precision, dimension(0:n_local_x+1,0:n_local_y+1,1:n_local_z_w)  :: tempr_fz_sfd

    !----- MPI related variables -----!
    
    integer :: status(MPI_STATUS_SIZE)

    vx_bar_gl    = 0.0D+00
    vy_bar_gl    = 0.0D+00
    vz_bar_gl    = 0.0D+00
    fx_sfd_gl    = 0.0D+00
    fy_sfd_gl    = 0.0D+00
    fz_sfd_gl    = 0.0D+00


    if ( my_id == master_id ) then          
       
       pz = 0

       do py = 0,dims(2)-1

          do px = 0,dims(1)-1
             
             if( (pz==0).and.(py==0).and.(px==0) )then
!         Since this is the first processor in x and y, save from sx-1 & sy-1 
                vx_bar_gl(sx-1:ex,sy-1:ey,sz_uv:ez_uv) = vx_bar(sx-1:ex,sy-1:ey,sz_uv:ez_uv)    
                vy_bar_gl(sx-1:ex,sy-1:ey,sz_uv:ez_uv) = vy_bar(sx-1:ex,sy-1:ey,sz_uv:ez_uv)    
                fx_sfd_gl(sx-1:ex,sy-1:ey,sz_uv:ez_uv) = fx_sfd(sx-1:ex,sy-1:ey,sz_uv:ez_uv)    
                fy_sfd_gl(sx-1:ex,sy-1:ey,sz_uv:ez_uv) = fy_sfd(sx-1:ex,sy-1:ey,sz_uv:ez_uv)    
             else
                
                coords(1) = px
                coords(2) = py
                coords(3) = pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                Call mpi_recv(sx_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_vx_bar,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_vy_bar,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_fx_sfd,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_fy_sfd,(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)


                jstart=1
                jend=n_local_y
                istart=1
                iend=n_local_x
                
!         If this is the first processor in x (or y), save from sx-1 (or sy-1) 
                if (px==0) istart=istart-1
                if (py==0) jstart=jstart-1
!         If this is the last processor in x (or y), save until ex+1 (or ey+1) 
                if (px==(dims(1)-1)) iend=iend+1
                if (py==(dims(2)-1)) jend=jend+1
                

                do k = 1,n_local_z_uv
                   do j = jstart,jend
                      do i = istart,iend
                         vx_bar_gl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_vx_bar(i,j,k)
                         vy_bar_gl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_vy_bar(i,j,k)
                         fx_sfd_gl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_fx_sfd(i,j,k)
                         fy_sfd_gl(sx_f-1+i,sy_f-1+j,sz_uv-1+k) = tempr_fy_sfd(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do

    else
       
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(vx_bar(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(vy_bar(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(fx_sfd(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)

       Call mpi_ssend(fy_sfd(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),(n_local_x+2)*(n_local_y+2)*n_local_z_uv,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)


    endif



    If (my_id==master_id) Then

       pz=0

       do py=0,dims(2)-1
          do px=0,dims(1)-1
             
             if((pz==0).and.(py==0).and.(px==0))then
!         Since this is the first processor in x and y, save from sx-1 & sy-1 
                vz_bar_gl(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w) = vz_bar(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w)    
                fz_sfd_gl(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w) = fz_sfd(sx_o-1:ex_o,sy_o-1:ey_o,sz_w:ez_w)    
             else
                
                coords(1)=px
                coords(2)=py
                coords(3)=pz
                
                !! Find the sender's rank from the senders virtual Cartesian coords.
                Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                
                Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                
                Call mpi_recv(tempr_vz_bar,(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                Call mpi_recv(tempr_fz_sfd,(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,sender_id,0, &
                     comm2d_quasiperiodic,status,ierr)
                
                jstart=1
                jend=n_local_y
                istart=1
                iend=n_local_x
!         If this is the first processor in x (or y), save from sx-1 (or sy-1) 
                if (px==0) istart=istart-1
                if (py==0) jstart=jstart-1
!         If this is the last processor in x (or y), save until ex+1 (or ey+1) 
                if (px==(dims(1)-1)) iend=iend+1
                if (py==(dims(2)-1)) jend=jend+1

                do k=1,n_local_z_w
                   do j=jstart,jend
                      do i=istart,iend
                         vz_bar_gl(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_vz_bar(i,j,k)
                         fz_sfd_gl(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_fz_sfd(i,j,k)
                      end do
                   end do
                end do
                
             end if
             
          end do
       end do
    Else
       Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
       
       Call mpi_ssend(vz_bar(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
       Call mpi_ssend(fz_sfd(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),(n_local_x+2)*(n_local_y+2)*n_local_z_w,mpi_double_precision,   &
            master_id,0,comm2d_quasiperiodic,ierr)
    End If

    

    
    return
  end subroutine gather_tables_exactsfd
  
  ! ****************************************************************************************
end module tpls_io
