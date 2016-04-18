module tpls_selective_frequency_damping

  implicit none

  double precision, allocatable, dimension(:,:,:) :: vx_bar
  double precision, allocatable, dimension(:,:,:) :: vy_bar
  double precision, allocatable, dimension(:,:,:) :: vz_bar

  double precision, allocatable, dimension(:,:,:) :: fx
  double precision, allocatable, dimension(:,:,:) :: fy
  double precision, allocatable, dimension(:,:,:) :: fz

  double precision :: omega, chi
  parameter ( omega = 0.1D+00 )
  parameter ( chi = 2.5D-01 )

  private
  public :: selective_frequency_damping

contains

  subroutine selective_frequency_damping(RHS_u,RHS_v,RHS_w,vx,vy,vz)

    use tpls_constants
    use tpls_io
    use tpls_configuration
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: vz

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(inout) :: RHS_u, RHS_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(inout) :: RHS_w

    !----- Initialization (istep==1) -----!

    if (istep==1) then
       allocate(vx_bar(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(vy_bar(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(vz_bar(sx-1:ex+1,sy-1:ey+1,sz_w :ez_w ))
       
       allocate(fx(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(fy(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(fz(sx-1:ex+1,sy-1:ey+1,sz_w :ez_w ))
       
       if (if_exact_restart_sfd) then
          
          call dataload_exact_sfd(vx_bar,vy_bar,vz_bar,fx,fy,fz)   
          
       else 
          vx_bar = vx
       	  vy_bar = vy
          vz_bar = vz
          
          fx = 0.0D+00
          fy = 0.0D+00
          fz = 0.0D+00
       endif
    endif

    vx_bar = vx_bar + dt*omega*( vx - vx_bar )
    vy_bar = vy_bar + dt*omega*( vy - vy_bar )
    vz_bar = vz_bar + dt*omega*( vz - vz_bar )

    call velocity_bc(vx_bar, vy_bar, vz_bar)

    fx = -chi*( vx - vx_bar )
    fy = -chi*( vy - vy_bar )
    fz = -chi*( vz - vz_bar )

    RHS_u = RHS_u + dt*fx
    RHS_v = RHS_v + dt*fy
    RHS_w = RHS_w + dt*fz

    if( mod(istep,iostep).eq.0 ) then
       call backup_exactrest_sfd(vx_bar,vy_bar,vz_bar,fx,fy,fz)
    endif

    if (istep==nsteps) then
       deallocate(vx_bar,vy_bar,vz_bar)
       deallocate(fx,fy,fz)
    endif
    
    return
  end subroutine selective_frequency_damping

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

end module tpls_selective_frequency_damping
