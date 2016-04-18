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

contains

  subroutine selective_frequency_damping(RHS_u,RHS_v,RHS_w,vx,vy,vz)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: vz

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(inout) :: RHS_u, RHS_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(inout) :: RHS_w

    vx_bar = vx_bar + dt*omega*( vx - vx_bar )
    vy_bar = vy_bar + dt*omega*( vy - vy_bar )
    vz_bar = vz_bar + dt*omega*( vz - vz_bar )

    fx_sfd = -chi*( vx - vx_bar )
    fy_sfd = -chi*( vy - vy_bar )
    fz_sfd = -chi*( vz - vz_bar )

    RHS_u = RHS_u + dt*fx_sfd
    RHS_v = RHS_v + dt*fy_sfd
    RHS_w = RHS_w + dt*fz_sfd

    return
  end subroutine selective_frequency_damping

end module tpls_selective_frequency_damping
