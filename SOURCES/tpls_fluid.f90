module tpls_fluid

  implicit none

  double precision, dimension(:,:,:), allocatable :: csf_u1, csf_u0
  double precision, dimension(:,:,:), allocatable :: csf_v1, csf_v0
  double precision, dimension(:,:,:), allocatable :: csf_w1, csf_w0

  double precision, dimension(:,:,:), allocatable :: conv1_u, conv0_u
  double precision, dimension(:,:,:), allocatable :: conv1_v, conv0_v
  double precision, dimension(:,:,:), allocatable :: conv1_w, conv0_w

  double precision, dimension(:,:,:), allocatable :: diffusion1_u
  double precision, dimension(:,:,:), allocatable :: diffusion1_v
  double precision, dimension(:,:,:), allocatable :: diffusion1_w

  private
  public :: advance_fluid, pressure_correction

contains

  subroutine advance_fluid(phi,vx,vy,vz,pressure)

    use tpls_constants
    use tpls_mpi
    use tpls_userchk
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(in)    :: phi
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(inout) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(inout) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(inout) :: pressure

    !----- Miscellaneous -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: RHS_u, RHS_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: RHS_w
    
    integer :: i, j, k
    double precision :: dummy, max_dummy

    t_temp = mpi_wtime()

    !----- Compute the right hand side for the velocity equation -----

    call make_rhs_velocity(RHS_u,RHS_v,RHS_w,phi,vx,vy,vz,pressure)

    !----- Update the velocity ------!
    
    vx = RHS_u
    vy = RHS_v
    vz = RHS_w

    call velocity_bc(vx, vy, vz)

    dummy = maxval(vx)
    call mpi_allreduce(dummy,max_dummy,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d_quasiperiodic,ierr)

    if ( my_id.eq.0 ) then
       write(*,*)'Iteration : ', istep, 'max-u vel    is ', max_dummy
    endif
    dummy = maxval(vy)
    call mpi_allreduce(dummy,max_dummy,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d_quasiperiodic,ierr)

    if ( my_id.eq.0 ) then
       write(*,*)'Iteration : ', istep, 'max-v vel    is ', max_dummy
    endif
    dummy = maxval(vz)
    call mpi_allreduce(dummy,max_dummy,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d_quasiperiodic,ierr)

    if ( my_id.eq.0 ) then
       write(*,*)'Iteration : ', istep, 'max-w vel    is ', max_dummy
    endif

    time_fluid = time_fluid + (mpi_wtime() - t_temp)
    
    return
  end subroutine advance_fluid
  
  subroutine make_rhs_velocity(RHS_u,RHS_v,RHS_w,phi,vx,vy,vz,pressure)

    use tpls_constants
    use tpls_mpi
    use tpls_levelset
    use tpls_configuration
    use tpls_selective_frequency_damping, only : selective_frequency_damping
    use mpi
    use tpls_io
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(in) :: phi
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w), intent(in) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(in) :: pressure

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(out) :: RHS_u, RHS_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w), intent(out) :: RHS_w

    !----- Miscellaneous -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: csf_u2, csf_v2
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: csf_w2

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: density
    double precision, dimension(sx-2:ex+2,sy-2:ey+2,sz_p:ez_p)   :: viscosity

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: diffusion2_u, diffusion2_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: diffusion2_w

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: conv2_u, conv2_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: conv2_w

    integer :: i, j, k

    !----- Initialization ( istep == 1 ) -----!
    
    if (istep.eq.1) then
       allocate(conv0_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(conv1_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       
       allocate(csf_u0(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(csf_u1(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       
       allocate(diffusion1_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       
       allocate(conv0_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(conv1_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       
       allocate(csf_v0(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(csf_v1(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       
       allocate(diffusion1_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       
       allocate(conv0_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w))
       allocate(conv1_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w))
       
       allocate(csf_w0(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w))
       allocate(csf_w1(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w))
       
       allocate(diffusion1_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w))

       conv0_u = 0.d0
       conv1_u = 0.d0
       
       diffusion1_u = 0.d0
       
       csf_u0 = 0.d0
       csf_u1 = 0.d0
       
       conv0_v = 0.d0
       conv1_v = 0.d0
       
       diffusion1_v = 0.d0
       
       csf_v0 = 0.d0
       csf_v1 = 0.d0
       
       conv0_w = 0.d0
       conv1_w = 0.d0
       
       diffusion1_w = 0.d0
       
       csf_w0 = 0.d0
       csf_w1 = 0.d0
       
    endif

    !----- Compute the density array -----!

    call get_density(density,phi)

    !----- Compute the viscosity array -----!
    
    call get_viscosity(viscosity,phi)
    
    !----- Compute the surface tension force -----!

    if ( scap .ne. 0.0D+00 ) then
       call compute_csf(csf_u2,csf_v2,csf_w2,phi,density)
    else
       csf_u2 = 0
       csf_v2 = 0
       csf_w2 = 0
    endif

    !----- Compute the diffusive terms -----!
    
    call get_diffusion_all(diffusion2_u,diffusion2_v,diffusion2_w,viscosity,density,vx,vy,vz)
    
    !----- Compute the convective terms -----!
    
    call get_conv_all(conv2_u,conv2_v,conv2_w,viscosity,density,vx,vy,vz)
    
    !----- Compute the right-hand side for the velocity -----!
    
    if((istep.EQ.1) .and. (.not. if_exact_restart)) then
       if (my_id==master_id) write(*,*) 'Performing forward Euler start for conv/diff terms' 
       RHS_u = vx - dt*conv2_u + dt*csf_u2 + dt*diffusion2_u
       RHS_v = vy - dt*conv2_v + dt*csf_v2 + dt*diffusion2_v
       RHS_w = vz - dt*conv2_w + dt*csf_w2 + dt*diffusion2_w
    
    elseif((istep.EQ.1) .and. (if_exact_restart)) then
       if (my_id==master_id) write(*,*) 'Performing exact restart for conv/diff terms' 
       call dataload_exact_conv(conv1_u,conv1_v,conv1_w,csf_u1,csf_v1,csf_w1, &
                                diffusion1_u,diffusion1_v,diffusion1_w)

    endif
   
    if((istep.GE.2) .or. (if_exact_restart)) then


       RHS_u = vx                                                 &
            - dt*(3.d0/2.d0*conv2_u      - 1.d0/2.d0*conv1_u)     &
            + dt*(3.d0/2.d0*csf_u2       - 1.d0/2.d0*csf_u1 )     &
            + dt*(3.d0/2.d0*diffusion2_u - 1.d0/2.d0*diffusion1_u)
       
       RHS_v = vy                                                 &
            - dt*(3.d0/2.d0*conv2_v      - 1.d0/2.d0*conv1_v)     &
            + dt*(3.d0/2.d0*csf_v2       - 1.d0/2.d0*csf_v1 )     &
            + dt*(3.d0/2.d0*diffusion2_v - 1.d0/2.d0*diffusion1_v)
       
       RHS_w = vz                                                 &
            - dt*(3.d0/2.d0*conv2_w      - 1.d0/2.d0*conv1_w)     &
            + dt*(3.d0/2.d0*csf_w2       - 1.d0/2.d0*csf_w1 )     &
            + dt*(3.d0/2.d0*diffusion2_w - 1.d0/2.d0*diffusion1_w)
       
    endif

    if ( if_sfd ) call selective_frequency_damping(RHS_u,RHS_v,RHS_w &
         , vx, vy, vz)
    
    call exchange2d(RHS_u,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(RHS_v,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(RHS_w,stride_w_xz, stride_w_yz, neighbours,ex,ey, ez_w,sx,sy, sz_w,comm2d_quasiperiodic)
    
    diffusion1_u = diffusion2_u
    diffusion1_v = diffusion2_v
    diffusion1_w = diffusion2_w
    
    conv0_u = conv1_u
    conv1_u = conv2_u
    
    conv0_v = conv1_v
    conv1_v = conv2_v
    
    conv0_w = conv1_w
    conv1_w = conv2_w
    
    csf_u0 = csf_u1
    csf_u1 = csf_u2
    
    csf_v0 = csf_v1
    csf_v1 = csf_v2
    
    csf_w0 = csf_w1
    csf_w1 = csf_w2

    if (if_exact_restart_save) then 
       if( mod(istep,iostep).eq.0 .and. (istep .gt. 0)) then

       call  backup_exactrest_conv(conv1_u,conv1_v,conv1_w,csf_u1,csf_v1,csf_w1, &
                                diffusion1_u,diffusion1_v,diffusion1_w)
       endif
    endif

    if (istep.eq.nsteps) then
       deallocate(conv0_u,conv1_u)
       deallocate(csf_u0,csf_u1)
       deallocate(diffusion1_u)
       deallocate(conv0_v,conv1_v)
       deallocate(csf_v0,csf_v1)
       deallocate(diffusion1_v)
       deallocate(conv0_w,conv1_w)
       deallocate(csf_w0,csf_w1)
       deallocate(diffusion1_w)
    endif
    
    return
  end subroutine make_rhs_velocity
  

  subroutine compute_csf(csf_u,csf_v,csf_w,phi,density)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(in) :: phi, density

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(out) :: csf_u, csf_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w), intent(out) :: csf_w

    !----- Miscellaneous -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p) :: fx_csf, fy_csf, fz_csf
    double precision :: rho_ugrid, rho_vgrid, rho_wgrid

    integer :: i, j, k

    call get_csf(fx_csf,fy_csf,fz_csf,phi)

    !Ajout JC

    forall (i = sx:ex, j = sy:ey , k = 0:maxn-2)
       csf_u(i,j,k) = (fx_csf(i,j,k+1) + fx_csf(i-1,j,k+1))/2.0D+00
       csf_v(i,j,k) = (fy_csf(i,j,k+1) + fy_csf(i,j-1,k+1))/2.0D+00
    end forall
    
    forall (i = sx:ex, j = sy:ey, k = 0:maxn-1)
       csf_w(i,j,k) = (fz_csf(i,j,k+1) + fz_csf(i,j,k))/2.0D+00
    end forall
    
    if(.not.density_matched) then
       
       do k = 0,maxn-2
          do j = sy,ey
             do i = sx,ex
                
                rho_ugrid     = (density(i+1,j+1,k+1)+density(i,j+1,k+1))/2.d0
                csf_u(i,j,k) = csf_u(i,j,k)/rho_ugrid
                
                rho_vgrid     = (density(i,j+1,k+1)+density(i,j,  k+1))/2.d0
                csf_v(i,j,k) = csf_v(i,j,k)/rho_vgrid
                
             end do
          end do
       end do
       
       do k = 0,maxn-1
          do j = sy,ey
             do i = sx,ex
                
                rho_wgrid     = (density(i,j+1,k+1)+density(i,j+1,k))/2.d0
                csf_w(i,j,k) = csf_w(i,j,k)/rho_wgrid
                
             end do
          end do
       end do
       
    endif

    if ( sx == 1 ) then
       csf_u(0,:,:) = 2.0D+00*csf_u(1,:,:) - csf_u(2,:,:)     
       csf_v(0,:,:) = 2.0D+00*csf_v(1,:,:) - csf_v(2,:,:)
       csf_w(0,:,:) = 2.0D+00*csf_w(1,:,:) - csf_w(2,:,:)
    endif

    if ( ex == ex_max ) then
       csf_u(ex_max+1,:,:) = 2*csf_u(ex_max,:,:) - csf_u(ex_max-1,:,:)
       csf_v(ex_max+1,:,:) = 2*csf_v(ex_max,:,:) - csf_v(ex_max-1,:,:)
       csf_w(ex_max+1,:,:) = 2*csf_w(ex_max,:,:) - csf_w(ex_max-1,:,:)
    endif
    
    return
  end subroutine compute_csf


  subroutine get_csf(fx_csf,fy_csf,fz_csf,phi)
    
    !---
    !
    !   Computation of the continuous surface tension force on the phi-grid
    !
    !---
    
    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none

    !----- Input -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)  ::  phi

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(out) :: fx_csf, fy_csf, fz_csf

    !----- Miscellaneous -----!

    double precision :: curvature


    integer :: i, j, k, n
    integer :: dist_x, dist_y, dist_z, counter

    double precision :: pi, max_kappa, norm_val, term1, term2, term3, dirac_phi

    double precision :: dphidx, dphidy, dphidz
    double precision :: d2phidx, d2phidy, d2phidz
    double precision :: d2phidxdy, d2phidxdz, d2phidydz

    double precision, dimension(27,10) :: A
    double precision, dimension(10)    :: x
    double precision, dimension(27)    :: b

    curvature = 0.

    n = maxn
    pi = 4.0D+00*atan(1.0D+00)
    max_kappa = 1.0D+00/(dmin1(dx,dy,dz))

    fx_csf = 0
    fy_csf = 0
    fz_csf = 0

    do k = 1,n-1
       do j = sy,ey
          do i = sx,ex

             if ( abs(phi(i,j,k)) .lt. 5*smooth_width ) then
                
                !----- Commpute curvature -----!
                
                A(:,1)  = 1.0D+00
                
                counter = 1
                do dist_z = -1,1
                   do dist_y = -1,1
                      do dist_x = -1,1
                         
                         b(counter)   = phi( i + dist_x , j + dist_y , k + dist_z )
                         
                         A(counter,2) = dx*dist_x
                         A(counter,3) = dy*dist_y
                         A(counter,4) = dz*dist_z
                         
                         counter = counter + 1
                         
                      enddo
                   enddo
                enddo
                
                A(:,5) = 0.5D+00 * A(:,2)**2.00D+00
                A(:,6) = 0.5D+00 * A(:,3)**2.00D+00
                A(:,7) = 0.5D+00 * A(:,4)**2.00D+00
                
                A(:,8)  = A(:,2)*A(:,3)
                A(:,9)  = A(:,2)*A(:,4)
                A(:,10) = A(:,3)*A(:,4)
                
                x = 0.0D+00
                call Least_Squares(A,x,b,27,10)
                
                dphidx = x(2)
                dphidy = x(3)
                dphidz = x(4)
                
                d2phidx = x(5)
                d2phidy = x(6)
                d2phidz = x(7)
                
                d2phidxdy = x(8)
                d2phidxdz = x(9)
                d2phidydz = x(10)
                
                norm_val = dphidx**2.0D+00 &
                     + dphidy**2.0D+00 &
                     + dphidz**2.0D+00
                
                if (norm_val.lt.1.0D-10) then
                   
                   curvature = 0.0D+00
                   
                else
                   
                   term1 = (d2phidx + d2phidy + d2phidz)/(norm_val**0.5D+00)
                   
                   term2 = dphidx*dphidx*d2phidx  + dphidy*dphidy*d2phidy  + dphidz*dphidz*d2phidz
                   term2 = term2/(norm_val**1.5D+00)
                   
                   term3 = dphidx*dphidy*d2phidxdy + dphidx*dphidz*d2phidxdz + dphidy*dphidz*d2phidydz
                   term3 = 2.0D+00*term3/(norm_val**1.5D+00)
                   
                   curvature = -(term1-term2-term3)
                   
                end if

                if ( curvature.lt.-max_kappa ) then
                   curvature = -max_kappa
                end if
                
                dirac_phi = Dirac_function(phi(i,j,k))
                fx_csf(i,j,k) = dirac_phi*scap*curvature*dphidx / dsqrt( norm_val )
                fy_csf(i,j,k) = dirac_phi*scap*curvature*dphidy / dsqrt( norm_val )
                fz_csf(i,j,k) = dirac_phi*scap*curvature*dphidz / dsqrt( norm_val )

             endif
             
          enddo
       enddo
    enddo
    
    fx_csf(:,:,0) = 3.0D+00 * fx_csf(:,:,1) &
         - 3.0D+00/2.D0+00 * fx_csf(:,:,2)  &
         + fx_csf(:,:,3)

    fx_csf(:,:,n) = 3.0D+00 * fx_csf(:,:,n-1) &
         - 3.0D+00/2.0D+00 * fx_csf(:,:,n-2)  &
         + fx_csf(:,:,n-3)

    fy_csf(:,:,0) = 3.0D+00 * fy_csf(:,:,1) &
         - 3.0D+00/2.D0+00 * fy_csf(:,:,2)  &
         + fy_csf(:,:,3)

    fy_csf(:,:,n) = 3.0D+00 * fy_csf(:,:,n-1) &
         - 3.0D+00/2.0D+00 * fy_csf(:,:,n-2)  &
         + fy_csf(:,:,n-3)

    fz_csf(:,:,0) = 3.0D+00 * fz_csf(:,:,1) &
         - 3.0D+00/2.D0+00 * fz_csf(:,:,2)  &
         + fz_csf(:,:,3)

    fz_csf(:,:,n) = 3.0D+00 * fz_csf(:,:,n-1) &
         - 3.0D+00/2.0D+00 * fz_csf(:,:,n-2)  &
         + fz_csf(:,:,n-3)
    
    
    if ( sx == 1 ) then
       fx_csf(0,:,:) = 2.0D+00*fx_csf(1,:,:) - fx_csf(2,:,:)
       fy_csf(0,:,:) = 2.0D+00*fy_csf(1,:,:) - fy_csf(2,:,:)
       fz_csf(0,:,:) = 2.0D+00*fz_csf(1,:,:) - fz_csf(2,:,:)
    end if
    
    if ( ex == ex_max ) then
       fx_csf(ex_max+1,:,:) = 2.0D+00*fx_csf(ex_max,:,:) - fx_csf(ex_max-1,:,:)
       fy_csf(ex_max+1,:,:) = 2.0D+00*fy_csf(ex_max,:,:) - fy_csf(ex_max-1,:,:)
       fz_csf(ex_max+1,:,:) = 2.0D+00*fz_csf(ex_max,:,:) - fz_csf(ex_max-1,:,:)
    end if

    call exchange2d(fx_csf,stride_p_xz,stride_p_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    call exchange2d(fy_csf,stride_p_xz,stride_p_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    call exchange2d(fz_csf,stride_p_xz,stride_p_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    
  end subroutine get_csf

  
  subroutine get_diffusion_all(diff_u,diff_v,diff_w,viscosity,density,u2,v2,w2)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(in) :: u2, v2
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-1), intent(in) :: w2
    double precision, dimension(sx-2:ex+2,sy-2:ey+2,0:maxn)  , intent(in) :: viscosity
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn)  , intent(in) :: density

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(out) :: diff_u, diff_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-1), intent(out) :: diff_w

    !----- Miscellaneous -----!
    
    integer :: i, j, k
    integer :: ip1, im1, jp1, jm1
    integer :: ex_loc
    
    double precision :: u_minusz, u_plusz
    double precision :: mu_duxE, mu_duxW, mu_duyN, mu_duyS, mu_duzI, mu_duzO
    double precision :: v_minusz, v_plusz
    double precision :: mu_dvxE, mu_dvxW, mu_dvyN, mu_dvyS, mu_dvzI, mu_dvzO
    double precision :: mu_dwxE, mu_dwxW, mu_dwyN, mu_dwyS, mu_dwzI, mu_dwzO
    double precision :: v_plusx, v_minusx
    double precision :: w_plusx, w_minusx
    double precision :: mu_plushalf_x_val, mu_minushalf_x_val
    double precision :: mu_plushalf_y_val, mu_minushalf_y_val
    double precision :: mu_plushalf_z_val,mu_minushalf_z_val
    
    double precision :: rho_ugrid, rho_vgrid, rho_wgrid

    double precision :: dudxp, dudxm, dudyp, dudym, dudzp, dudzm
    double precision :: dvdxp, dvdxm, dvdyp, dvdym, dvdzp, dvdzm
    double precision :: dwdxp, dwdxm, dwdyp, dwdym, dwdzp, dwdzm

    !----- Computation of the diffusion for the x-velocity component -----

    diff_u = 0.0D+00
    
    do k = 0,maxn-2
       do j = sy,ey
          do i = sx,ex
             
             ip1 = i+1
             im1 = i-1
             jp1 = j+1
             jm1 = j-1
             
             if(k.eq.0)then
                u_minusz = -u2(i,j,0)
             else
                u_minusz = u2(i,j,k-1)
             end if
             
             if(k.eq.(maxn-2))then
                u_plusz = -u2(i,j,maxn-2)
             else
                u_plusz = u2(i,j,k+1)
             end if

             !----- Compute gradients of u -----!

             dudxp = ( u2(i+1,j,k) - u2(i,j,k) ) / dx
             dudxm = ( u2(i,j,k) - u2(i-1,j,k) ) / dx

             dudyp = ( u2(i,j+1,k) - u2(i,j,k) ) / dy
             dudym = ( u2(i,j,k) - u2(i,j-1,k) ) / dy

             dudzp = ( u_plusz - u2(i,j,k) ) / dz
             dudzm = ( u2(i,j,k) - u_minusz ) /dz

             !----- Compute average viscosity on u_grid -----

             mu_plushalf_x_val  = viscosity(i,j,k+1)
             mu_minushalf_x_val = viscosity(i-1,j,k+1)
             
             mu_plushalf_y_val  = ( viscosity(i,j,k+1) + viscosity(i,j+1,k+1) &
                  + viscosity(i-1,j,k+1) + viscosity(i-1,j+1,k+1) )/4.0D+00
             mu_minushalf_y_val = ( viscosity(i,j,k+1) + viscosity(i,j-1,k+1) &
                  + viscosity(i-1,j,k+1) + viscosity(i-1,j-1,k+1) )/4.0D+00

             mu_plushalf_z_val  = ( viscosity(i,j,k+1) + viscosity(i,j,k+2) &
                  + viscosity(i-1,j,k+1) + viscosity(i-1,j,k+2) )/4.0D+00
             mu_minushalf_z_val = ( viscosity(i,j,k+1) + viscosity(i,j,k) &
                  + viscosity(i-1,j,k) + viscosity(i-1,j,k+1) )/4.0D+00

             !----- Compute div( mu*grad(u) ) -----!

             mu_duxE = mu_plushalf_x_val  * dudxp
             mu_duxW = mu_minushalf_x_val * dudxm

             mu_duyN = mu_plushalf_y_val  * dudyp
             mu_duyS = mu_minushalf_z_val * dudym

             mu_duzO = mu_plushalf_z_val  * dudzp
             mu_duzI = mu_minushalf_z_val * dudzm

             rho_ugrid = (density(i,j,k+1)+density(i-1,j,k+1))/2.0D+00
             
             diff_u(i,j,k)= ( ((mu_duxE-mu_duxW)/dx) + ((mu_duyN-mu_duyS)/dy) + ((mu_duzO-mu_duzI)/dz) ) /rho_ugrid
             
          end do
       end do
    end do

    !----- Computation of the diffusion for the y-velocity component -----

    diff_v = 0.0D+00

    do k = 0,maxn-2
       do j = sy,ey
          do i = sx,ex
             
             ip1 = i+1
             im1 = i-1
             jp1 = j+1
             jm1 = j-1
             
             if(k.eq.0)then
                v_minusz = -v2(i,j,0)
             else
                v_minusz = v2(i,j,k-1)
             end if
             
             if(k.eq.(maxn-2))then
                v_plusz = -v2(i,j,maxn-2)
             else
                v_plusz = v2(i,j,k+1)
             end if

             !----- Compute gradients of v -----!

             dvdxp = ( v2(i+1,j,k) - v2(i,j,k) ) / dx
             dvdxm = ( v2(i,j,k) - v2(i-1,j,k) ) / dx

             dvdyp = ( v2(i,j+1,k) - v2(i,j,k) ) / dy
             dvdym = ( v2(i,j,k) - v2(i,j-1,k) ) / dy

             dvdzp = ( v_plusz - v2(i,j,k) ) / dz
             dvdzm = ( v2(i,j,k) - v_minusz ) / dz

             !----- Compute average viscosity on v-grid -----!

             mu_plushalf_x_val  = ( viscosity(i,j,k+1) + viscosity(i+1,j,k+1) &
                  + viscosity(i,j-1,k+1) + viscosity(i+1,j-1,k+1) )/4.0D+00
             mu_minushalf_x_val = ( viscosity(i,j,k+1) + viscosity(i-1,j,k+1) &
                  + viscosity(i,j-1,k+1) + viscosity(i-1,j-1,k+1) )/4.0D+00

             mu_plushalf_y_val  = viscosity(i,j,k+1)
             mu_minushalf_y_val = viscosity(i,j-1,k+1)

             mu_plushalf_z_val  = ( viscosity(i,j,k+1) + viscosity(i,j,k+2) &
                  + viscosity(i,j-1,k+1) + viscosity(i,j-1,k+2))/4.0D+00
             mu_minushalf_z_val = ( viscosity(i,j,k+1) + viscosity(i,j,k) &
                  + viscosity(i,j-1,k+1) + viscosity(i,j-1,k) )/4.0D+00

             !----- Compute div(mu*grad(v)) -----!

             mu_dvxE = mu_plushalf_x_val  * dvdxp
             mu_dvxW = mu_minushalf_x_val * dvdxm
             
             mu_dvyN = mu_plushalf_y_val  * dvdyp
             mu_dvyS = mu_minushalf_y_val * dvdym
             
             mu_dvzO= mu_plushalf_z_val * dvdzp
             mu_dvzI=mu_minushalf_z_val * dvdzm
             
             rho_vgrid=(density(i,j+1,k+1)+density(i,j,  k+1))/2.0D+00
             
             diff_v(i,j,k)=(((mu_dvxE-mu_dvxW)/dx)+((mu_dvyN-mu_dvyS)/dy)+((mu_dvzO-mu_dvzI)/dz))/rho_vgrid
             
             
          end do
       end do
    end do

    !----- Computation of the diffusion for the z-velocity component -----

    diff_w = 0.0D+00

    do k=1,maxn-2
       do j=sy,ey
          do i=sx,ex
             
             ip1=i+1
             im1=i-1
             jp1=j+1
             jm1=j-1

             !----- Compute gradients of w -----!

             dwdxp = ( w2(i+1,j,k) - w2(i,j,k) ) / dx
             dwdxm = ( w2(i,j,k) - w2(i-1,j,k) ) / dx

             dwdyp = ( w2(i,j+1,k) - w2(i,j,k) ) / dy
             dwdym = ( w2(i,j,k) - w2(i,j-1,k) ) / dy

             dwdzp = ( w2(i,j,k+1) - w2(i,j,k) ) / dz
             dwdzm = ( w2(i,j,k) - w2(i,j,k-1) ) / dz

             !----- Compute average viscosity on w-grid -----!

             mu_plushalf_x_val  = ( viscosity(i,j,k) + viscosity(i+1,j,k) &
                  + viscosity(i,j,k+1) + viscosity(i+1,j,k+1) )/4.0D+00
             mu_minushalf_x_val = ( viscosity(i,j,k) + viscosity(i,j,k+1) &
                  + viscosity(i-1,j,k) + viscosity(i-1,j,k+1) )/4.0D+00

             mu_plushalf_y_val  = ( viscosity(i,j,k) + viscosity(i,j+1,k) &
                  + viscosity(i,j,k+1) + viscosity(i,j+1,k+1) )/4.0D+00
             mu_minushalf_y_val = ( viscosity(i,j,k) + viscosity(i,j-1,k) &
                  + viscosity(i,j,k+1) + viscosity(i,j-1,k+1) )/4.0D+00

             mu_plushalf_z_val  = viscosity(i,j,k+1)
             mu_minushalf_z_val = viscosity(i,j,k)

             !----- Compute div(mu*grad(w)) -----!

             mu_dwxE = mu_plushalf_x_val  * dwdxp
             mu_dwxW = mu_minushalf_x_val * dwdxm
             
             mu_dwyN = mu_plushalf_y_val  * dwdyp
             mu_dwyS = mu_minushalf_y_val * dwdym
             
             mu_dwzO = mu_plushalf_z_val  * dwdzp
             mu_dwzI = mu_minushalf_z_val * dwdzm
             
             rho_wgrid = (density(i,j+1,k+1)+density(i,j+1,k))/2.0D+00
             
             diff_w(i,j,k)=(((mu_dwxE-mu_dwxW)/dx)+((mu_dwyN-mu_dwyS)/dy)+((mu_dwzO-mu_dwzI)/dz))/rho_wgrid
             
          end do
       end do
    end do
    
    return
  end subroutine get_diffusion_all
  
  !**********************************************************************************************
  !
  ! Subroutine to compute the convective term at level maxn.  
  ! There are two parts to the convective term: the ordinary convective term, and a part involving
  ! gradients of the viscosity (the transposed part of the stress tensor).
  ! The differencing is fully flux-conservative.
  ! To compute the cell-averaged viscosities, the augmented viscosity array is needed (MPI stuff).
  !
  !**********************************************************************************************
  
  subroutine get_conv_all(conv_u2,conv_v2,conv_w2,viscosity,density,u2,v2,w2)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(in) :: u2, v2
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-1), intent(in) :: w2
    double precision, dimension(sx-2:ex+2,sy-2:ey+2,0:maxn)  , intent(in) :: viscosity
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn)  , intent(in) :: density

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(out) :: conv_u2, conv_v2
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-1), intent(out) :: conv_w2

    !----- Miscellaneous -----!
    
    integer :: i, j, k
    integer :: im1, ip1
    integer :: jm1, jp1
    integer :: ex_loc
    
    double precision :: ub, vb
    double precision :: tempu, tempv, tempw
    double precision :: temp_conv1, temp_conv2
    double precision :: dudx, dudy, dudz
    double precision :: dvdx, dvdy, dvdz
    double precision :: dwdx, dwdy, dwdz
    double precision :: Xconv, Yconv, Zconv
    double precision :: mu_dudxE, mu_dudxW, mu_dvdxN, mu_dvdxS, mu_dwdxO, mu_dwdxI, &
         mu_dudyE, mu_dudyW, mu_dvdyN, mu_dvdyS, mu_dwdyO, mu_dwdyI, &
         mu_dudzE, mu_dudzW, mu_dvdzN, mu_dvdzS, mu_dwdzO, mu_dwdzI
    double precision :: mu_plushalf_x_val, mu_minushalf_x_val, mu_plushalf_y_val, &
         mu_minushalf_y_val, mu_plushalf_z_val, mu_minushalf_z_val
    double precision :: v_plusx, v_minusx, w_plusx, w_minusx
    double precision :: rho_ugrid, rho_vgrid, rho_wgrid

    double precision :: dudxp, dudxm, dudyp, dudym, dudzp, dudzm
    double precision :: dvdxp, dvdxm, dvdyp, dvdym, dvdzp, dvdzm
    double precision :: dwdxp, dwdxm, dwdyp, dwdym, dwdzp, dwdzm
    
    ! ******************************  X direction

    conv_u2=0.d0

    do k = 0,maxn-2
       do j = sy, ey
          do i = sx, ex
             
             im1 = i-1
             ip1 = i+1

             tempv = (v2(i,j,k) + v2(im1,j,k) + v2(i,j+1,k) + v2(im1,j+1,k))/4.d0
             tempw = (w2(i,j,k) + w2(im1,j,k) + w2(i,j,k+1) + w2(im1,j,k+1))/4.d0
             
             dudx  = (u2(i+1,j,k)-u2(i-1,j,k))/(2.0D+00*dx)
             Xconv = dudx*u2(i,j,k)
             
             ub   = u2(i,j+1,k)
             vb   = u2(i,j-1,k)
             dudy = 0.5d0*(ub-vb)/dy
             
             Yconv = dudy*tempv
             
             if( k == 0 ) then
                dudz = (u2(i,j,k)+u2(i,j,k+1))/(2.0D+00*dz)
             elseif( k == maxn-2 ) then
                dudz = -(u2(i,j,k)+u2(i,j,k-1))/(2.0D+00*dz)
             else
                ub = u2(i,j,k+1)
                vb = u2(i,j,k-1)
                dudz = 0.5d0*(ub-vb)/dz
             endif
             
             Zconv = dudz*tempw
             temp_conv1 = Xconv+Yconv+Zconv

             !----- Compute gradients along x -----!

             dudxp = ( u2(i+1,j,k) - u2(i,j,k) ) / dx
             dudxm = ( u2(i,j,k) - u2(i-1,j,k) ) / dx

             dvdxp = ( v2(i,j+1,k) - v2(i-1,j+1,k) ) / dx
             dvdxm = ( v2(i,j,k) - v2(i-1,j,k) ) / dx

             dwdxp = ( w2(i,j,k+1) - w2(i-1,j,k+1) ) / dx
             dwdxm = ( w2(i,j,k) - w2(i-1,j,k) ) / dx

             !----- Compute average viscosity on u-grid -----!

             mu_plushalf_x_val  = viscosity(i,j,k+1)
             mu_minushalf_x_val = viscosity(i-1,j,k+1)
             
             mu_plushalf_y_val  = ( viscosity(i,j,k+1) + viscosity(i,j+1,k+1) &
                  + viscosity(i-1,j,k+1) + viscosity(i-1,j+1,k+1) )/4.0D+00
             mu_minushalf_y_val = ( viscosity(i,j,k+1) + viscosity(i,j-1,k+1) &
                  + viscosity(i-1,j,k+1) + viscosity(i-1,j-1,k+1) )/4.0D+00

             mu_plushalf_z_val  = ( viscosity(i,j,k+1) + viscosity(i,j,k+2) &
                  + viscosity(i-1,j,k+1) + viscosity(i-1,j,k+2) )/4.0D+00
             mu_minushalf_z_val = ( viscosity(i,j,k+1) + viscosity(i,j,k) &
                  + viscosity(i-1,j,k) + viscosity(i-1,j,k+1) )/4.0D+00

             !----- Compute div(mu*grad(u).T) -----!
             
             mu_dudxE = mu_plushalf_x_val  * dudxp
             mu_dudxW = mu_minushalf_x_val * dudxm
             
             mu_dvdxN = mu_plushalf_y_val  * dvdxp
             mu_dvdxS = mu_minushalf_y_val * dvdxm
             
             mu_dwdxO = mu_plushalf_z_val  * dwdxp
             mu_dwdxI = mu_minushalf_z_val * dwdxm
             
             rho_ugrid=(density(i+1,j,k+1)+density(i,j,k+1))/2.d0
             
             temp_conv2=(((mu_dudxE-mu_dudxW)/dx)+((mu_dvdxN-mu_dvdxS)/dy)+((mu_dwdxO-mu_dwdxI)/dz))/rho_ugrid
             
             conv_u2(i,j,k) = temp_conv1-temp_conv2
          enddo
       enddo
    enddo

    conv_v2 = 0.d0

    do k=0, maxn-2
       do j=sy, ey
          do i=sx, ex
             
             jm1=j-1

!!$             tempu=(u2(i-1,jm1,k)+u2(i-1,j,k)+u2(i,jm1,k)+u2(i,j,k))/4.d0
!!$             tempw=(w2(i,jm1,k)+w2(i,j,k)+w2(i,jm1,k+1)+w2(i,j,k+1))/4.d0

             ! Ajout JC

             tempu=(u2(i,jm1,k)+u2(i,j,k)+u2(i+1,jm1,k)+u2(i+1,j,k))/4.d0
             tempw=(w2(i,jm1,k)+w2(i,j,k)+w2(i,jm1,k+1)+w2(i,j,k+1))/4.d0

             
             dvdy  = 0.5d0*(v2(i,j+1,k)-v2(i,j-1,k))/dy
             Yconv = dvdy*v2(i,j,k)
             
             v_minusx = v2(i-1,j,k)
             v_plusx  = v2(i+1,j,k)
             
             dvdx = 0.5d0*(v_plusx-v_minusx)/dx
             Xconv=dvdx*tempu
             
             if(k==0)then
                dvdz=(v2(i,j,k)+v2(i,j,k+1))/(2.d0*dz)
             elseif(k==maxn-2)then
                dvdz=-(v2(i,j,k)+v2(i,j,k-1))/(2.d0*dz)
             else
                dvdz=0.5d0*(v2(i,j,k+1)-v2(i,j,k-1))/dz
             endif
             
             Zconv = dvdz*tempw
             temp_conv1 = Xconv + Yconv + Zconv
             
             !----- Compute gradients along y -----!

             dudyp = ( u2(i+1,j,k) - u2(i+1,j-1,k) ) / dy
             dudym = ( u2(i,j,k)   - u2(i,j-1,k) ) / dy

             dvdyp = ( v2(i,j+1,k) - v2(i,j,k) ) / dy
             dvdym = ( v2(i,j,k) - v2(i,j-1,k) ) / dy

             dwdyp = ( w2(i,j,k+1) - w2(i,j-1,k+1) ) / dy
             dwdym = ( w2(i,j,k) - w2(i,j-1,k) ) / dy

             !----- Compute average viscosity on v-grid -----!

             mu_plushalf_x_val  = ( viscosity(i,j,k+1) + viscosity(i+1,j,k+1) &
                  + viscosity(i,j-1,k+1) + viscosity(i+1,j-1,k+1) )/4.0D+00
             mu_minushalf_x_val = ( viscosity(i,j,k+1) + viscosity(i-1,j,k+1) &
                  + viscosity(i,j-1,k+1) + viscosity(i-1,j-1,k+1) )/4.0D+00
             
             mu_plushalf_y_val  = viscosity(i,j,k+1)
             mu_minushalf_y_val = viscosity(i,j-1,k+1)

             mu_plushalf_z_val  = ( viscosity(i,j,k+1) + viscosity(i,j,k+2) &
                  + viscosity(i,j-1,k+1) + viscosity(i,j-1,k+2))/4.0D+00
             mu_minushalf_z_val = ( viscosity(i,j,k+1) + viscosity(i,j,k) &
                  + viscosity(i,j-1,k+1) + viscosity(i,j-1,k) )/4.0D+00

             !----- Compute div(mu*grad(u).T) -----!
             
             mu_dudyE=mu_plushalf_x_val  * dudyp
             mu_dudyW=mu_minushalf_x_val * dudym
             
             mu_dvdyN=mu_plushalf_y_val  * dvdyp
             mu_dvdyS=mu_minushalf_y_val * dvdym
             
             mu_dwdyO=mu_plushalf_z_val  * dwdyp
             mu_dwdyI=mu_minushalf_z_val * dwdym
             
             rho_vgrid=(density(i,j+1,k+1)+density(i,j,k+1))/2.0D+00
             
             temp_conv2 = (((mu_dudyE-mu_dudyW)/dx)+((mu_dvdyN-mu_dvdyS)/dy)+((mu_dwdyO-mu_dwdyI)/dz))/rho_vgrid
             
             
             conv_v2(i,j,k) = temp_conv1 - temp_conv2
          enddo
       enddo
    enddo

    ! ******************************  Z direction  
    
    conv_w2 = 0.0D+00

    do k = 1,maxn-2
       do j = sy,ey
          do i = sx,ex

             tempu = (u2(i,j,k-1)+u2(i,j,k)+u2(i+1,j,k-1)+u2(i+1,j,k))/4.d0
             tempv = (v2(i,j,k-1)+v2(i,j,k)+v2(i,j+1,k-1)+v2(i,j+1,k))/4.d0
             
             ub = w2(i,j,k+1)
             vb = w2(i,j,k-1)

             dwdz = 0.5*(ub-vb)/dz
             Zconv = dwdz*w2(i,j,k)
             
             w_minusx = w2(i-1,j,k)
             w_plusx  = w2(i+1,j,k)
             
             dwdx  = 0.5d0*(w_plusx-w_minusx)/dx
             Xconv = dwdx*tempu
             
             ub = w2(i,j+1,k)
             vb = w2(i,j-1,k)

             dwdy  = 0.5d0*(ub-vb)/dy
             Yconv = dwdy*tempv
             
             temp_conv1 = Xconv + Yconv + Zconv

             !----- Compute gradients along z-direction -----!

             dudzp = ( u2(i+1,j,k) - u2(i+1,j,k-1) ) / dz
             dudzm = ( u2(i,j,k) - u2(i,j,k-1) ) / dz

             dvdzp = ( v2(i,j+1,k) - v2(i,j+1,k-1) ) / dz
             dvdzm = ( v2(i,j,k) - v2(i,j,k-1) ) / dz

             dwdzp = ( w2(i,j,k+1) - w2(i,j,k) ) / dz
             dwdzm = ( w2(i,j,k) - w2(i,j,k-1) ) / dz

             !----- Compute average viscosity on w-grid -----!

             mu_plushalf_x_val  = ( viscosity(i,j,k) + viscosity(i+1,j,k) &
                  + viscosity(i,j,k+1) + viscosity(i+1,j,k+1) )/4.0D+00
             mu_minushalf_x_val = ( viscosity(i,j,k) + viscosity(i,j,k+1) &
                  + viscosity(i-1,j,k) + viscosity(i-1,j,k+1) )/4.0D+00

             mu_plushalf_y_val  = ( viscosity(i,j,k) + viscosity(i,j+1,k) &
                  + viscosity(i,j,k+1) + viscosity(i,j+1,k+1) )/4.0D+00
             mu_minushalf_y_val = ( viscosity(i,j,k) + viscosity(i,j-1,k) &
                  + viscosity(i,j,k+1) + viscosity(i,j-1,k+1) )/4.0D+00

             mu_plushalf_z_val=  viscosity(i,j,k+1)
             mu_minushalf_z_val= viscosity(i,j,k)

             !----- Compute div(mu*grad(u).T) -----!
           
             mu_dudzE = mu_plushalf_x_val  * dudzp
             mu_dudzW = mu_minushalf_x_val * dudzm
           
             mu_dvdzN = mu_plushalf_y_val  * dvdzp
             mu_dvdzS = mu_minushalf_y_val * dvdzm
             
             mu_dwdzO = mu_plushalf_z_val  * dwdzp
             mu_dwdzI = mu_minushalf_z_val * dwdzm
             
             rho_wgrid = (density(i,j,k+1)+density(i,j,k))/2.0D+00
             
             temp_conv2 = ( (mu_dudzE-mu_dudzW)/dx + (mu_dvdzN-mu_dvdzS)/dy + (mu_dwdzO-mu_dwdzI)/dz )/rho_wgrid
             
             conv_w2(i,j,k) = temp_conv1-temp_conv2

          enddo
       enddo
    enddo
    
    return
  end subroutine get_conv_all

  subroutine get_curvature(curvature,phi)
    
    !---
    !
    !   Computation of the continuous surface tension force on the phi-grid
    !
    !---
    
    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none

    !----- Input -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)  ::  phi

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn)  , intent(out) :: curvature

    !----- Miscellaneous -----!

    integer :: i, j, k, n
    integer :: dist_x, dist_y, dist_z, counter

    double precision :: pi, max_kappa, norm_val, term1, term2, term3, dirac_phi

    double precision :: dphidx, dphidy, dphidz
    double precision :: d2phidx, d2phidy, d2phidz
    double precision :: d2phidxdy, d2phidxdz, d2phidydz

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn,3) :: normal_xp, normal_xm
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn,3) :: normal_yp, normal_ym
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn,3) :: normal_zp, normal_zm

    double precision, dimension(27,10) :: A
    double precision, dimension(10)    :: x
    double precision, dimension(27)    :: b
    double precision :: local_error, global_error, x_val, z_val, exact_curvature

    curvature = 0.
    local_error = 0.0D+00
    global_error = 0.0D+00

    n = maxn
    pi = 4.0D+00*atan(1.0D+00)
    max_kappa = 1.0D+00/(dmin1(dx,dy,dz))

    do k = 1,n-1
       do j = sy,ey
          do i = sx,ex

             x_val = i*dx - dx/2.0D+00 - Lx/2.0D+00
             z_val = k*dz - dz/2.0D+00 - Lz/2.0D+00

             !----- Exact curvature ----!

             exact_curvature = dsqrt( x_val**2.0D+00 + z_val**2.0D+00 )
             exact_curvature = -1.0D+00 / exact_curvature

             !----- Commpute curvature -----!

             A(:,1)  = 1.0D+00

             counter = 1
             do dist_z = -1,1
                do dist_y = -1,1
                   do dist_x = -1,1
                      
                      b(counter)   = phi( i + dist_x , j + dist_y , k + dist_z )
                      
                      A(counter,2) = dx*dist_x
                      A(counter,3) = dy*dist_y
                      A(counter,4) = dz*dist_z
                      
                      counter = counter + 1
                      
                   enddo
                enddo
             enddo
             
             A(:,5) = 0.5D+00 * A(:,2)**2.00D+00
             A(:,6) = 0.5D+00 * A(:,3)**2.00D+00
             A(:,7) = 0.5D+00 * A(:,4)**2.00D+00
             
             A(:,8)  = A(:,2)*A(:,3)
             A(:,9)  = A(:,2)*A(:,4)
             A(:,10) = A(:,3)*A(:,4)
             
             x = 0.0D+00
             call Least_Squares(A,x,b,27,10)
             
             dphidx = x(2)
             dphidy = x(3)
             dphidz = x(4)

             d2phidx = x(5)
             d2phidy = x(6)
             d2phidz = x(7)

             d2phidxdy = x(8)
             d2phidxdz = x(9)
             d2phidydz = x(10)

             norm_val = dphidx**2.0D+00 &
                      + dphidy**2.0D+00 &
                      + dphidz**2.0D+00
             
             if (norm_val.lt.1.0D-10) then

                curvature(i,j,k) = 0.0D+00

             else

                term1 = (d2phidx + d2phidy + d2phidz)/(norm_val**0.5D+00)
                
                term2 = dphidx*dphidx*d2phidx  + dphidy*dphidy*d2phidy  + dphidz*dphidz*d2phidz
                term2 = term2/(norm_val**1.5D+00)
                
                term3 = dphidx*dphidy*d2phidxdy + dphidx*dphidz*d2phidxdz + dphidy*dphidz*d2phidydz
                term3 = 2.0D+00*term3/(norm_val**1.5D+00)
                
                curvature(i,j,k) = -(term1-term2-term3)
                
             end if

             
             
             dirac_phi = Dirac_function(phi(i,j,k))
             local_error = local_error + dirac_phi*(exact_curvature - curvature(i,j,k))**2.0D+00

             curvature(i,j,k) = dirac_phi*curvature(i,j,k)

          enddo
       enddo
    enddo

    if ( sx == 1 ) then
       curvature(0,:,:) = 3.0D+00 * curvature(1,:,:) &
            - 3.0D+00/2.0D+00 * curvature(2,:,:) &
            + curvature(3,:,:)
    endif

    if ( ex == ex_max ) then
       curvature(ex_max+1,:,:) = 3.0D+00 * curvature(ex_max,:,:) &
            - 3.0D+00/2.0D+00 * curvature(ex_max-1,:,:) &
            + curvature(ex_max-2,:,:)
    endif

    curvature(:,:,0) = 3.0D+00*curvature(:,:,1) &
         - 3.0D+00/2.0D+00*curvature(:,:,2) &
         + curvature(:,:,3)

    curvature(:,:,maxn) = 3.0D+00*curvature(:,:,maxn-1) &
         - 3.0D+00/2.0D+00*curvature(:,:,maxn-2) &
         + curvature(:,:,maxn-3)

    call exchange2d(curvature,stride_p_xz,stride_p_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    
  end subroutine get_curvature

  subroutine get_curvature_original(curvature, phi)

    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none

    !----- Input -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)  ::  phi

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn)  , intent(out) :: curvature

    !----- Miscellaneous -----!

    integer :: n, i, j, k
    integer :: im1, ip1
    integer :: jm1, jp1
    
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: fi
    double precision :: term1, term2, term3
    double precision :: nrm_val
    double precision :: dfi_x, dfi_y, dfi_z
    double precision :: phi_val
    double precision :: dirac_phi, pi, max_kappa
    
    double precision :: d2fi_x, d2fi_y, d2fi_z
    double precision :: dfi_x_up, dfi_x_down, d2fi_xz
    double precision :: dfi_y_up, dfi_y_down, d2fi_yz
    double precision :: dfi_x_in, dfi_x_out, d2fi_xy

    double precision :: local_error, global_error, x_val, z_val, exact_curvature

    local_error = 0.0D+00
    global_error = 0.0D+00

    n = maxn
    pi = 4.0D+00*atan(1.0D+00)
    max_kappa = 1.0D+00/(dmin1(dx,dy,dz))

    fi = phi
    
    ! First, compute curvature
    do k = 1,n-1  
       do j = sy,ey
          do i = sx,ex

             x_val = i*dx - dx/2.0D+00 - Lx/2.0D+00
             z_val = k*dz - dz/2.0D+00 - Lz/2.0D+00

!!$             !----- Exact curvature ----!
!!$
!!$             exact_curvature = dsqrt( x_val**2.0D+00 + z_val**2.0D+00 )
!!$             exact_curvature = -1.0D+00 / exact_curvature

             
             ip1 = i+1
             im1 = i-1
             jp1 = j+1
             jm1 = j-1
                  
             d2fi_x = (fi(ip1,j,k)-2.d0*fi(i,j,k)+fi(im1,j,k))/(dx**2.d0)
             d2fi_y = (fi(i,jp1,k)-2.d0*fi(i,j,k)+fi(i,jm1,k))/(dy**2.d0)
             d2fi_z = (fi(i,j,k+1)-2.d0*fi(i,j,k)+fi(i,j,k-1))/(dz**2.d0)
             
             dfi_x = (fi(ip1,j,k)-fi(im1,j,k))/(2.d0*dx)
             dfi_y = (fi(i,jp1,k)-fi(i,jm1,k))/(2.d0*dy)
             dfi_z = (fi(i,j,k+1)-fi(i,j,k-1))/(2.d0*dz)
             
             dfi_x_up   = (fi(ip1,j,k+1)-fi(im1,j,k+1))/(2.d0*dx)
             dfi_x_down = (fi(ip1,j,k-1)-fi(im1,j,k-1))/(2.d0*dx)
             d2fi_xz    = (dfi_x_up-dfi_x_down)/(2.d0*dz)
             
             dfi_y_up   = (fi(i,jp1,k+1)-fi(i,jm1,k+1))/(2.d0*dy)
             dfi_y_down = (fi(i,jp1,k-1)-fi(i,jm1,k-1))/(2.d0*dy)
             d2fi_yz    = (dfi_y_up-dfi_y_down)/(2.d0*dz)
             
             dfi_x_out = (fi(ip1,jp1,k)-fi(im1,jp1,k))/(2.d0*dx)
             dfi_x_in  = (fi(ip1,jm1,k)-fi(im1,jm1,k))/(2.d0*dx)
             d2fi_xy   = (dfi_x_out-dfi_x_in)/(2.d0*dy)
             
             nrm_val = dfi_x**2.0D+00 &
                     + dfi_y**2.0D+00 &
                     + dfi_z**2.0D+00
             
             if (nrm_val.lt.1.0D-10) then
                curvature(i,j,k) = 0.0D+00
             else

                term1 = (d2fi_x + d2fi_y + d2fi_z)/(nrm_val**0.5D+00)
                
                term2 = dfi_x*dfi_x*d2fi_x  + dfi_y*dfi_y*d2fi_y  + dfi_z*dfi_z*d2fi_z
                term2 = term2/(nrm_val**1.5D+00)
                
                term3 = dfi_x*dfi_y*d2fi_xy + dfi_x*dfi_z*d2fi_xz + dfi_y*dfi_z*d2fi_yz
                term3 = 2.0D+00*term3/(nrm_val**1.5D+00)
                
                curvature(i,j,k) = -(term1-term2-term3)

!!$                curvature(i,j,k) = (d2fi_y + d2fi_z)*dfi_x**2.D0+00 &
!!$                     + (d2fi_x + d2fi_z)*dfi_y**2.0D+00             &
!!$                     + (d2fi_x + d2fi_y)*dfi_z**2.0D+00             &
!!$                     - 2.0D+00*dfi_x*dfi_y*d2fi_xy                  &
!!$                     - 2.0D+00*dfi_x*dfi_z*d2fi_xz                  &
!!$                     - 2.0D+00*dfi_y*dfi_z*d2fi_yz
!!$
!!$                curvature(i,j,k) = -curvature(i,j,k)/(nrm_val**(3.0D+00/2.0D+00))
                
             end if
!!$             
!!$             if (curvature(i,j,k).gt.max_kappa) then
!!$                curvature(i,j,k) = max_kappa
!!$             end if

             dirac_phi = Dirac_function(phi(i,j,k))
             local_error = local_error + dirac_phi*(exact_curvature - curvature(i,j,k))**2.0D+00
             curvature(i,j,k) = dirac_phi*exact_curvature

             
          end do
       end do
    end do

    call mpi_allreduce(local_error,global_error,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d_quasiperiodic,ierr)
!!$    global_error = global_error * dx*dy*dz
!!$    global_error = dsqrt( global_error )
!!$    if ( my_id == 0 ) then
!!$       write(*,*) ' GLOBAL ERROR (Original) : ', global_error
!!$    endif


    return
  end subroutine get_curvature_original

  subroutine velocity_bc(vx, vy, vz)

    use tpls_constants
    use tpls_userchk
    use tpls_mpi
    use mpi

    !----- Inputs/Outputs -----!

    double precision, dimension(sx-1:ex+1, sy-1:ey+1, sz_uv:ez_uv), intent(inout) :: vx, vy
    double precision, dimension(sx-1:ex+1, sy-1:ey+1, sz_w:ez_w), intent(inout) :: vz

    !----- Miscellaneous -----!

    integer :: i, j, k
    double precision, dimension(0:maxn-2) :: u_inlet
    double precision, dimension(0:maxn) :: phi_inlet

    !----- Inflow boundary condition : u = u_inlet -----!
    
    call user_inflow(u_inlet,phi_inlet)
    
    if ( sx == 1 ) then
       
       do j = sy,ey
          vx(0,j,:) = 2.0D+00*u_inlet(:) - vx(2,j,:)
          vx(1,j,:) = u_inlet(:)
       enddo
       
    end if
    
    !----- Outflow boundary condition : du/dx = 0 -----!
    
    if ( ex == ex_max ) then       
       !vx(ex_max+1,:,:) = vx(ex_max,:,:)
       vx(ex_max+1,:,:) = 2*vx(ex_max,:,:) - vx(ex_max-1, :, :)
    end if
    
    !----- No-slip boundary conditions on the walls : u = 0 -----!

    vx(:,:,sz_uv) = (1.0D+00/3.0D+00)*vx(:,:,sz_uv+1)
    vx(:,:,ez_uv) = (1.0D+00/3.0D+00)*vx(:,:,ez_uv-1)

        !----- Inflow boundary condition : v = 0 -----!
    
    if ( sx == 1 ) then
       vy(0,:,:) = -vy(1,:,:)
    end if
    
    !----- Outflow boundary condition : dv/dx = 0 -----!
    
    if ( ex == ex_max ) then
       !vy(ex_max+1,:,:) = vy(ex_max,:,:)
       vy(ex_max+1,:,:) = 2*vy(ex_max,:,:) - vy(ex_max-1, :, :)
    end if
    
    !----- No-slip boundary condition on the walls : v = 0 -----!
    
    vy(:,:,sz_uv) = (1.0D+00/3.0D+00)*vy(:,:,sz_uv+1)
    vy(:,:,ez_uv) = (1.0D+00/3.0D+00)*vy(:,:,ez_uv-1)
    
    !----- No-slip boundary condition on the walls : w = 0 -----!
    
    vz(:,:,0)      = 0.0D+00
    vz(:,:,maxn-1) = 0.0D+00
    
    !----- Inflow boundary condition : w = 0 -----!
    
    if ( sx == 1 ) then
       vz(0,:,:) = -vz(1,:,:)       
    end if
    
    !----- Outflow boundary condition : dw/dx = 0 -----!
    
    if ( ex == ex_max ) then
       !vz(ex_max+1,:,:) = vz(ex_max,:,:)
       vz(ex_max+1,:,:) = 2*vz(ex_max,:,:) - vz(ex_max-1, :, :)
    end if
    
    
    call exchange2d(vx,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(vy,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(vz,stride_w_xz, stride_w_yz, neighbours,ex,ey, ez_w,sx,sy, sz_w,comm2d_quasiperiodic)
    
    return
  end subroutine velocity_bc

  subroutine pressure_correction(vx, vy, vz, p, phi)
    
    use tpls_constants
    use tpls_levelset
    use tpls_mpi
    use tpls_userchk
    use mpi
    implicit none

    !----- Inputs / Ouputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(inout) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-1), intent(inout) :: vz
    
    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)      :: phi
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)      :: p

    !----- Miscellaneous -----!
    
    integer :: i,j,k
    integer :: sx_loc,ex_loc
    double precision :: rho_ugrid,rho_vgrid,rho_wgrid
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: density
   
    if ( .not. density_matched) then

       call get_density(density, phi)
    
       do k = 0,maxn-2
          do j = sy,ey
             do i = sx,ex
                rho_ugrid = ( density(i+1,j+1,k+1) + density(i,j+1,k+1) ) / 2.0D+00
                vx(i,j,k) = vx(i,j,k) - (dt/dx)*( 1.d0/rho_ugrid ) * ( p(i,j,k+1) - p(i-1,j,k+1) )
             end do
          end do
       end do
       
       do k = 0,maxn-2
          do j = sy,ey
             do i = sx,ex
                rho_vgrid = ( density(i,j+1,k+1) + density(i,j,k+1) ) / 2.d0
                vy(i,j,k) = vy(i,j,k) - (dt/dy)*( 1.d0/rho_vgrid ) * ( p(i,j,k+1) - p(i,j-1,k+1) )             
             end do
          end do
       end do
       
       do k = 0,maxn-1
          do j = sy,ey
             do i = sx,ex
                rho_wgrid = ( density(i,j+1,k+1) + density(i,j+1,k) ) / 2.0D+00
                vz(i,j,k) = vz(i,j,k) - (dt/dz)*( 1.d0/rho_wgrid ) * ( p(i,j,k+1) - p(i,j,k) )
             end do
          end do
       end do

    else

       vx(sx:ex, sy:ey, 0:maxn-2) = vx(sx:ex, sy:ey, 0:maxn-2) &
            - (dt/dx) * ( p(sx:ex, sy:ey, 1:maxn-1) - p(sx-1:ex-1, sy:ey, 1:maxn-1) )
       vy(sx:ex, sy:ey, 0:maxn-2) = vy(sx:ex, sy:ey, 0:maxn-2) &
            - (dt/dy) * ( p(sx:ex, sy:ey, 1:maxn-1) - p(sx:ex, sy-1:ey-1, 1:maxn-1) )
       vz(sx:ex, sy:ey, 0:maxn-1) = vz(sx:ex, sy:ey, 0:maxn-1) &
            - (dt/dz) * ( p(sx:ex, sy:ey, 1:maxn)   - p(sx:ex, sy:ey, 0:maxn-1) )
       
    endif
       
    call velocity_bc(vx, vy, vz)
    
    return
  end subroutine pressure_correction


    
end module tpls_fluid
