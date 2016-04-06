module tpls_userchk

  use tpls_constants
  use tpls_mpi
  implicit none

  double precision, dimension(:,:,:), allocatable :: vx_old, vy_old, vz_old
  double precision, dimension(:,:,:), allocatable :: vx_bas, vy_bas, vz_bas
  double precision, dimension(:,:,:), allocatable :: vx_pertnorm, vy_pertnorm, vz_pertnorm
  double precision :: pertnorm, basnorm, pertprobe, basprobe, probe_xpos=10.d0, probe_zpos=3.d0
  integer :: probe_xind, probe_zind, probe_procid

contains
  




!-------------------------------------------------------------------------------------------------------------------------------





  subroutine user_postprocess(phi,pr,vx,vy,vz)

    !-----

!!$    This subroutine is called at the end of each time step. The velocity, pressure and levelset functions
!!$    are in-out arguments allowing to do some run-time post-processing if needed.

    !-----

    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use tpls_io
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(inout) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(inout) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(inout) :: pr, phi

    !----- Miscellaneous -----!

    integer :: i, j, k, n, ifich = 111, ifich2=123, ifich3=134
    double precision :: residu, z_val, x_val, y_val, pi

    pi = 4.0*atan(1.0)
    
    !----- Begining of the subroutine -----!

!!$    if ( istep.EQ.0 ) then 
!!$
!!$       do k = 0,maxn
!!$          do j = sy,ey
!!$             do i = sx,ex
!!$
!!$                z_val = k*dz - dz/2.0D+00 - Lz/2.0D+00
!!$
!!$                if ( z_val .GE. 0.0D+00 ) then
!!$                   phi(i,j,k) =  -z_val + height
!!$                else
!!$                   phi(i,j,k) =  z_val + height
!!$                endif
!!$
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$    endif

    !----- Computing the residual for steady state monitoring -----!

    if ( istep.EQ. 0 ) then

       allocate(vx_old(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(vy_old(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
       allocate(vz_old(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w))

       vx_old = 0.0D+00
       vy_old = 0.0D+00
       vz_old = 0.0D+00
      
       if ( my_id == master_id ) open( unit = ifich , file = 'residu.dat' )
       
       if (if_record_perturbation_norm) then 
          
          allocate(vx_pertnorm(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
          allocate(vy_pertnorm(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
          allocate(vz_pertnorm(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w))
          allocate(vx_bas(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
          allocate(vy_bas(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv))
          allocate(vz_bas(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w))
      
          vx_bas = vx 
          vy_bas = vy 
          vz_bas = vz 

          call inner_product(basnorm &
            , phi, phi &
            , vx_bas, vy_bas, vz_bas &
            , vx_bas, vy_bas, vz_bas)

       
          if ( my_id == master_id ) open( unit = ifich2 , file = 'pertnorm.dat' )

          if ( my_id == master_id ) then
            write(ifich2,*) istep, basnorm
            call flush(ifich2)
          endif
        
!         Setting a w-vel. probe approximately at x=10, z=1
          probe_zind=floor(probe_zpos/dz)
!!        Need to find out which processor has this x-point
          probe_procid=ceiling(probe_xpos/Lx*num_procs_x)-1
          probe_xind=floor((probe_xpos-probe_procid*(Lx/num_procs_x))/dz)
          if (my_id==probe_procid) then
                  basprobe=vz(sx+probe_xind,1,probe_zind) 
                  open( unit = ifich3 , file = 'pertprobe.dat' )
                  write(ifich3,*) istep, basprobe
                  call flush(ifich3)
          endif

       endif

    endif

    vx_old = vx - vx_old
    vy_old = vy - vy_old
    vz_old = vz - vz_old
    


    call inner_product(residu &
         , phi, phi &
         , vx_old, vy_old, vz_old &
         , vx_old, vy_old, vz_old)

 
    residu = residu/(Lx*Ly*Lz)
    residu = dsqrt(residu)
    residu = residu/dt

    vx_old = vx
    vy_old = vy
    vz_old = vz


    if ( my_id == master_id ) then
       write(ifich,*) istep, residu
       call flush(ifich)
    endif
    
    if ( istep == nsteps ) then
       if ( my_id == master_id ) close(ifich)
    endif

    if (if_record_perturbation_norm) then   
       vx_pertnorm = vx - vx_bas
       vy_pertnorm = vy - vy_bas
       vz_pertnorm = vz - vz_bas
    
       call inner_product(pertnorm &
            , phi, phi &
            , vx_pertnorm, vy_pertnorm, vz_pertnorm &
            , vx_pertnorm, vy_pertnorm, vz_pertnorm)

       pertnorm = dsqrt(pertnorm)

       if ( my_id == master_id ) then
          write(ifich2,*) istep, pertnorm
          call flush(ifich2)
       endif
       

       if (my_id==probe_procid) then
              pertprobe=vz(sx+probe_xind,1,probe_zind)-basprobe
              write(ifich3,*) istep, pertprobe
              call flush(ifich3)
       endif

       if ( istep == nsteps ) then
           if ( my_id == master_id ) close(ifich2)
           if ( my_id == probe_procid ) close(ifich3)
       endif

    endif

    return
  end subroutine user_postprocess





!-------------------------------------------------------------------------------------------------------------------------------





  subroutine user_inflow(u_inlet,phi_inlet)

    !-----

!!$    This subroutine allows the user to easily defined the inlet profile for the velocity and the
!!$    levelset function.

    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none

    !----- Outputs -----!

    double precision, dimension(0:maxn-2) :: u_inlet
    double precision, dimension(0:maxn)   :: phi_inlet

    !----- Intermediate arrays -----!

    double precision, dimension(0:maxn)   :: eta, z_val
    
    !-----Miscellaneous -----!

    integer          :: i, j, k, n
    double precision :: pi
    double precision :: dummy
    double precision :: L_inv

    pi = 4.d0*atan(1.d0)
    n = maxn
    dummy = 150.d0

    do k = 0,n

       z_val(k) = dz*k - (dz/2.d0) - Lz/2.0D+00

       if ( z_val(k) .GE. 0.0D+00 ) then
          phi_inlet(k) = -z_val(k) + height
       else
          phi_inlet(k) =  z_val(k) + height
       endif

    enddo

    !----- Three Poiseuille streams joining at the inflow -----!

    do k = 0,maxn-2

       if ( phi_inlet(k+1).GE.0.0D+00 ) then
          u_inlet(k) = 3.0D+00/2.0D+00*(1.0D+00 - z_val(k+1)**2.0D+00)
       else
          if ( z_val(k+1).GE.0.0D+00 ) then
             eta(k) = 2.0D+00*(z_val(k+1)-1.5D+00)
             u_inlet(k) = 1.0D+00/4.0D+00*(1.0D+00-eta(k)**2.0D+00)
          else
             eta(k) = 2.0D+00*(z_val(k+1)+1.5D+00)
             u_inlet(k) = 1.0D+00/4.0D+00*(1.0D+00-eta(k)**2.0D+00)
          endif
       endif
!!$       u_inlet(k) = 1.0 - (z_val(k+1)/2.0D+00)**2.0D+00
       
    enddo
!!$
!!$    u_inlet = 0.0D+00
    
    return
  end subroutine user_inflow






!-------------------------------------------------------------------------------------------------------------------------------





  subroutine user_initial_conditions(vx,vy,vz,pressure,phi)

    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none

    !----- Outputs -----!

    double precision, dimension(0:maxl,0:maxm,0:maxn-2) :: vx, vy
    double precision, dimension(0:maxl,0:maxm,0:maxn-1) :: vz
    double precision, dimension(0:maxl,0:maxm,0:maxn)   :: pressure, phi

    !----- Intermediate arrays -----!

    double precision, dimension(0:maxn-2) :: u_inlet
    double precision, dimension(0:maxn)   :: phi_inlet

    !----- Miscellaneous -----!

    integer :: i, j, k, n
    double precision :: heaviside

    n = maxn

    call user_inflow(u_inlet,phi_inlet)
    
    do j = 0,maxm
       do i = 0,maxl
          
          vx(i,j,:)  = u_inlet
          phi(i,j,:) = phi_inlet

       enddo
    enddo

    vy  = 0.0D+00
    vz  = 0.0D+00

!!$    pressure = 0.0D+00

    return
  end subroutine user_initial_conditions





!-------------------------------------------------------------------------------------------------------------------------------





end module tpls_userchk
