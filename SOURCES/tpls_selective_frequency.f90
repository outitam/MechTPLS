module tpls_selective_frequency_damping

  implicit none

  double precision, allocatable, dimension(:,:,:) :: vx_bar
  double precision, allocatable, dimension(:,:,:) :: vy_bar
  double precision, allocatable, dimension(:,:,:) :: vz_bar

  double precision, allocatable, dimension(:,:,:) :: fx_sfd
  double precision, allocatable, dimension(:,:,:) :: fy_sfd
  double precision, allocatable, dimension(:,:,:) :: fz_sfd

  double precision :: omega, chi
  parameter ( omega = 0.1D-01 )
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
       
       allocate(fx_sfd(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(fy_sfd(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(fz_sfd(sx-1:ex+1,sy-1:ey+1,sz_w :ez_w ))
       
       if (if_exact_restart_sfd) then
          
          call dataload_exact_sfd(vx_bar,vy_bar,vz_bar,fx_sfd,fy_sfd,fz_sfd)   
          
       else 
          vx_bar = vx
       	  vy_bar = vy
          vz_bar = vz
          
          fx_sfd = 0.0D+00
          fy_sfd = 0.0D+00
          fz_sfd = 0.0D+00
       endif
    endif

    vx_bar = vx_bar + dt*omega*( vx - vx_bar )
    vy_bar = vy_bar + dt*omega*( vy - vy_bar )
    vz_bar = vz_bar + dt*omega*( vz - vz_bar )

    fx_sfd = -chi*( vx - vx_bar )
    fy_sfd = -chi*( vy - vy_bar )
    fz_sfd = -chi*( vz - vz_bar )

    RHS_u = RHS_u + dt*fx_sfd
    RHS_v = RHS_v + dt*fy_sfd
    RHS_w = RHS_w + dt*fz_sfd

    if( mod(istep,iostep).eq.0 ) then
       call backup_exactrest_sfd(vx_bar,vy_bar,vz_bar,fx_sfd,fy_sfd,fz_sfd)
    endif

    if (istep==nsteps) then
       deallocate(vx_bar,vy_bar,vz_bar)
       deallocate(fx_sfd,fy_sfd,fz_sfd)
    endif
    
    return
  end subroutine selective_frequency_damping

end module tpls_selective_frequency_damping
