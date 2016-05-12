module tpls_levelset

  implicit none

  private
  public :: advance_levelset, get_density, get_viscosity

contains

  
  subroutine advance_levelset(phi,vx,vy,vz)

    !-----
    !
    !     This subroutine advances the levelset equation in time based on a SSP-RK 3 scheme.
    !     It also calls for the redistancing subroutine.
    !
    !     INPUTS
    !     ------
    !
    !     phi        : levelset field at time n
    !     vx, vy, vz : velocity fields at time n
    !
    !     OUTPUT
    !     ------
    !
    !     phi        : levelset field at time n+1
    !
    !-----

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(inout) :: phi
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: vz

    !----- Arrays for the RHS -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: conv2_phi
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: phi_temp

    t_temp = mpi_wtime()

!!$    !----- Advection of the levelset function using SSP-RK3 -----!
!!$
!!$    call convective_terms_phi(conv2_phi,phi,vx,vy,vz)
!!$    phi_temp = phi - dt * conv2_phi
!!$    call boundary_conditions_levelset(phi_temp)
!!$
!!$    call convective_terms_phi(conv2_phi,phi_temp,vx,vy,vz)
!!$    phi_temp = (3./4.)*phi + (1./4.)*(phi_temp - dt*conv2_phi)
!!$    call boundary_conditions_levelset(phi_temp)
!!$    
!!$    call convective_terms_phi(conv2_phi,phi_temp,vx,vy,vz)
!!$    phi = (1./3.)*phi + (2./3.)*(phi_temp - dt*conv2_phi)
!!$    call boundary_conditions_levelset(phi)

    !----- Advection of the levelset function using RK2 -----!

    call convective_terms_phi(conv2_phi,phi,vx,vy,vz)
    phi_temp = phi - dt * conv2_phi
    call boundary_conditions_levelset(phi_temp)
    
    call convective_terms_phi(conv2_phi,phi_temp,vx,vy,vz)
    phi_temp = phi_temp - dt * conv2_phi
    call boundary_conditions_levelset(phi_temp)

    phi = .5*(phi + phi_temp)
    call boundary_conditions_levelset(phi)

    !----- Re-initialization to the signed distance function -----
    
    call redistancing_levelset(phi)

    !----- Miscellaneous -----!

    t_temp2 = mpi_wtime()
    time_levelset = time_levelset + (t_temp2 - t_temp)

    return
  end subroutine advance_levelset





  subroutine convective_terms_phi(conv2_phi,phi,vx,vy,vz)

    !-----
    !
    !     This subroutine computes the advection term in the levelset equation  in a
    !     non-conservative way using either WENO5 or HOUC5.
    !
    !     INPUTS
    !     ------
    !
    !     phi        : levelset field
    !     vx, vy, vz : velocity fields
    !
    !     OUTPUT
    !     ------
    !
    !     conv2_phi  : advection term of the levelset : vx*dphidx + vy*dphidy + vz*dphidxz
    !
    !-----

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(in) :: phi
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: vz

    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(out) :: conv2_phi

    !----- Miscellaneous -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: dfx, dfy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: dfz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p) :: dfx_, dfy_, dfz_

    integer :: i, j, k

    conv2_phi  = 0.0D+00

    !call WENO5_FD(phi,vx,vy,vz,dfx_,dfy_,dfz_)
    call HOUC5(dfx_,dfy_,dfz_,phi,vx,vy,vz)

    conv2_phi(sx:ex,sy:ey,1:maxn-1) = &
           (vx(sx+1:ex+1,sy:ey,0:maxn-2)+vx(sx:ex,sy:ey,0:maxn-2))/2.0D+00 * dfx_(sx:ex,sy:ey,1:maxn-1) &
         + (vy(sx:ex,sy+1:ey+1,0:maxn-2)+vy(sx:ex,sy:ey,0:maxn-2))/2.0D+00 * dfy_(sx:ex,sy:ey,1:maxn-1) &
         + (vz(sx:ex,sy:ey,1:maxn-1)+vz(sx:ex,sy:ey,0:maxn-2))/2.0D+00 * dfz_(sx:ex,sy:ey,1:maxn-1)

    return
  end subroutine convective_terms_phi




  subroutine WENO5(phi_,vx,vy,vz,dfx,dfy,dfz)

    !-----
    !
    !     This subroutine computed the upwinding derivatives of the levelset function
    !     using a standard WENO5 procedure.
    !
    !-----

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   , intent(in) :: phi_
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) , intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   , intent(in) :: vz
    
    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(out) :: dfx, dfy, dfz

    !----- WENO-5 parameters -----!
    
    double precision :: eps_weno
    parameter ( eps_weno = 1.0D-06 )

    double precision :: sigma_1, sigma_2, sigma_3
    parameter ( sigma_1 = 0.1D+00 )
    parameter ( sigma_2 = 0.6D+00 )
    parameter ( sigma_3 = 0.3D+00 )

    double precision :: v1, v2, v3, v4, v5
    double precision :: S1, S2, S3
    double precision :: a1, a2, a3
    double precision :: w1, w2, w3

    !----- Miscellaneous -----!

    integer :: i, j, k, l, n
    double precision, dimension(sx-3:ex+3,sy-3:ey+3,sz_p:ez_p) :: phi

    n = maxn
    l = maxl

    dfx = 0.0D+00
    dfy = 0.0D+00
    dfz = 0.0D+00

    !+++++ Creation of the augmented array +++++!

    phi = 0.0D+00
    phi(sx:ex,sy:ey,:) = phi_(sx:ex,sy:ey,:)

    if ( sx == 1 ) then
       phi(0,sy:ey,:) = phi_(0,sy:ey,:)
    endif
    if ( ex == ex_max ) then
       phi(ex_max+1,sy:ey,:) = phi_(ex_max+1,sy:ey,:)
    endif
    
    ! Augmented phi: Exchange first-, second- and third-order halos
    
    call exchange2d_augaug1(phi,stride_p_augaug1_xz,stride_p_augaug1_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)  
    
    call exchange2d_augaug2(phi,stride_p_augaug2_xz,stride_p_augaug2_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    
    call exchange2d_augaug3(phi,stride_p_augaug3_xz,stride_p_augaug3_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)

    !+++++ Computation of the WENO 5 fluxes in all three directions +++++!

    do k = 1,maxn-1
       do j = sy,ey
          do i = sx,ex

             !----- Computation of the fluxes in the x-direction -----!


             if ( ( i .GE. 4 ) .AND. ( i .LE. l-4 ) ) then

                if ( (vx(i+1,j,k-1)+vx(i,j,k-1)) .GE. 0.0D+00 ) then

                   v1 = ( phi(i-2,j,k) - phi(i-3,j,k) ) / dx
                   v2 = ( phi(i-1,j,k) - phi(i-2,j,k) ) / dx
                   v3 = ( phi(i,j,k)   - phi(i-1,j,k) ) / dx
                   v4 = ( phi(i+1,j,k) - phi(i,j,k)   ) / dx
                   v5 = ( phi(i+2,j,k) - phi(i+1,j,k) ) / dx

                else

                   v1 = ( phi(i+3,j,k) - phi(i+2,j,k) ) / dx
                   v2 = ( phi(i+2,j,k) - phi(i+1,j,k) ) / dx
                   v3 = ( phi(i+1,j,k) - phi(i,j,k)   ) / dx
                   v4 = ( phi(i,j,k)   - phi(i-1,j,k) ) / dx
                   v5 = ( phi(i-1,j,k) - phi(i-2,j,k) ) / dx

                endif

                S1 = (13.0D+00/12.0D+00) * ( v1 - 2.0D+00*v2 + v3 )**2.0D+00 &
                     + (1.0D+00/4.0D+00) * ( v1 - 4.0D+00*v2 + 3.0D+00*v3)**2.0D+00

                S2 = (13.0D+00/12.0D+00) * (v2 - 2.0D+00*v3 + v4)**2.0D+00 &
                     + (1.0D+00/4.0D+00) * (v2 - v4)**2.0D+00

                S3 = (13.0D+00/12.0D+00) * (v3 - 2.0D+00*v4 + v5)**2.0D+00 &
                     + (1.0D+00/4.0D+00) * (3.0D+00*v3 - 4.0D+00*v4 + v5)**2.0D+00

                a1 = sigma_1 / (eps_weno + S1)**2.0D+00
                a2 = sigma_2 / (eps_weno + S2)**2.0D+00
                a3 = sigma_3 / (eps_weno + S3)**2.0D+00

                w1 = a1 / (a1 + a2 + a3)
                w2 = a2 / (a1 + a2 + a3)
                w3 = a3 / (a1 + a2 + a3)

                dfx(i,j,k) = w1 * ( (1.0D+00/3.0D+00)  * v1   &
                                  - (7.0D+00/6.0D+00)  * v2   &
                                  + (11.0D+00/6.0D+00) * v3 ) &
                           + w2 * ( (-1.0D+00/6.0D+00) * v2   &
                                  + (5.0D+00/6.0D+00)  * v3   &
                                  + (1.0D+00/3.0D+00)  * v4 ) &
                           + w3 * ( (1.0D+00/3.0D+00)  * v3   &
                                  + (5.0D+00/6.0D+00)  * v4   &
                                  - (1.0D+00/6.0D+00)  * v5 )

             else

                if ( (vx(i,j,k-1)+vx(i+1,j,k-1)).LT.0.0D+00) then
                   dfx(i,j,k) = ( phi(i+1,j,k) - phi(i,j,k) ) / dx
                else
                   dfx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / dx
                endif

             endif
             
             !----- Computation of the fluxes in the y-direction -----!
             
             if ( (vy(i,j,k-1)+vy(i,j+1,k-1)) .GE. 0.0D+00 ) then
                
                v1 = ( phi(i,j-2,k) - phi(i,j-3,k) ) / dx
                v2 = ( phi(i,j-1,k) - phi(i,j-2,k) ) / dx
                v3 = ( phi(i,j,k)   - phi(i,j-1,k) ) / dx
                v4 = ( phi(i,j+1,k) - phi(i,j,k)   ) / dx
                v5 = ( phi(i,j+2,k) - phi(i,j+1,k) ) / dx
                
             else
                
                v1 = ( phi(i,j+3,k) - phi(i,j+2,k) ) / dx
                v2 = ( phi(i,j+2,k) - phi(i,j+1,k) ) / dx
                v3 = ( phi(i,j+1,k) - phi(i,j,k)   ) / dx
                v4 = ( phi(i,j,k)   - phi(i,j-1,k) ) / dx
                v5 = ( phi(i,j-1,k) - phi(i,j-2,k) ) / dx
                
             endif
             
             S1 = (13.0D+00/12.0D+00) * ( v1 - 2.0D+00*v2 + v3 )**2.0D+00 &
                  + (1.0D+00/4.0D+00) * ( v1 - 4.0D+00*v2 + 3.0D+00*v3)**2.0D+00
             
             S2 = (13.0D+00/12.0D+00) * (v2 - 2.0D+00*v3 + v4)**2.0D+00 &
                  + (1.0D+00/4.0D+00) * (v2 - v4)**2.0D+00
             
             S3 = (13.0D+00/12.0D+00) * (v3 - 2.0D+00*v4 + v5)**2.0D+00 &
                  + (1.0D+00/4.0D+00) * (3.0D+00*v3 - 4.0D+00*v4 + v5)**2.0D+00
             
             a1 = sigma_1 / (eps_weno + S1)**2.0D+00
             a2 = sigma_2 / (eps_weno + S2)**2.0D+00
             a3 = sigma_3 / (eps_weno + S3)**2.0D+00
             
             w1 = a1 / (a1 + a2 + a3)
             w2 = a2 / (a1 + a2 + a3)
             w3 = a3 / (a1 + a2 + a3)
             
             dfy(i,j,k) = w1 * ( (1.0D+00/3.0D+00)  * v1   &
                               - (7.0D+00/6.0D+00)  * v2   &
                               + (11.0D+00/6.0D+00) * v3 ) &
                        + w2 * ( (-1.0D+00/6.0D+00) * v2   &
                               + (5.0D+00/6.0D+00)  * v3   &
                               + (1.0D+00/3.0D+00)  * v4 ) &
                        + w3 * ( (1.0D+00/3.0D+00)  * v3   &
                               + (5.0D+00/6.0D+00)  * v4   &
                              - (1.0D+00/6.0D+00)  * v5 )

             !----- Computation of the fluxes in the z-direction -----!
             
             if ( ( k .GE. 4 ) .AND. ( k .LE. maxn-4 ) ) then
                
                if ( ( vz(i,j,k)+vz(i,j,k-1) ) .GE. 0.0D+00 ) then
                   
                   v1 = ( phi(i,j,k-2) - phi(i,j,k-3) ) / dz
                   v2 = ( phi(i,j,k-1) - phi(i,j,k-2) ) / dz
                   v3 = ( phi(i,j,k)   - phi(i,j,k-1) ) / dz
                   v4 = ( phi(i,j,k+1) - phi(i,j,k)   ) / dz
                   v5 = ( phi(i,j,k+2) - phi(i,j,k+1) ) / dz

                else

                   v1 = ( phi(i,j,k+3) - phi(i,j,k+2) ) / dz
                   v2 = ( phi(i,j,k+2) - phi(i,j,k+1) ) / dz
                   v3 = ( phi(i,j,k+1) - phi(i,j,k)   ) / dz
                   v4 = ( phi(i,j,k)   - phi(i,j,k-1) ) / dz
                   v5 = ( phi(i,j,k-1) - phi(i,j,k-2) ) / dz

                endif

                S1 = (13.0D+00/12.0D+00) * ( v1 - 2.0D+00*v2 + v3 )**2.0D+00 &
                     + (1.0D+00/4.0D+00) * ( v1 - 4.0D+00*v2 + 3.0D+00*v3)**2.0D+00

                S2 = (13.0D+00/12.0D+00) * (v2 - 2.0D+00*v3 + v4)**2.0D+00 &
                     + (1.0D+00/4.0D+00) * (v2 - v4)**2.0D+00

                S3 = (13.0D+00/12.0D+00) * (v3 - 2.0D+00*v4 + v5)**2.0D+00 &
                     + (1.0D+00/4.0D+00) * (3.0D+00*v3 - 4.0D+00*v4 + v5)**2.0D+00

                a1 = sigma_1 / (eps_weno + S1)**2.0D+00
                a2 = sigma_2 / (eps_weno + S2)**2.0D+00
                a3 = sigma_3 / (eps_weno + S3)**2.0D+00

                w1 = a1 / (a1 + a2 + a3)
                w2 = a2 / (a1 + a2 + a3)
                w3 = a3 / (a1 + a2 + a3)

                dfz(i,j,k) = w1 * ( (1.0D+00/3.0D+00)  * v1   &
                                  - (7.0D+00/6.0D+00)  * v2   &
                                  + (11.0D+00/6.0D+00) * v3 ) &
                           + w2 * ( (-1.0D+00/6.0D+00) * v2   &
                                  + (5.0D+00/6.0D+00)  * v3   &
                                  + (1.0D+00/3.0D+00)  * v4 ) &
                           + w3 * ( (1.0D+00/3.0D+00)  * v3   &
                                  + (5.0D+00/6.0D+00)  * v4   &
                                  - (1.0D+00/6.0D+00)  * v5 )

             else

                if ( (vz(i,j,k)+vz(i,j,k-1)).LT.0.0D+00 ) then
                   dfz(i,j,k) = ( phi(i,j,k+1) - phi(i,j,k) ) / dz
                else
                   dfz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dz 
                endif

             endif



          enddo
       enddo
    enddo

!!$    call exchange2d(dfx,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
!!$    call exchange2d(dfy,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
!!$    call exchange2d(dfz,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    return
  end subroutine WENO5

  



  
  subroutine get_levelset_status(phi0,sign_mx,h)

    !-----
    !
    !     Computation of the mollified sign function based on the Heaviside function.
    !
    !     INPUTS
    !     ------
    !
    !     phi : levelset field
    !     h   : deprecated
    !
    !     OUTPUT
    !     ------
    !
    !     sign_mx : sign of phi over the whole domain.
    !
    !-----

    use tpls_constants
    use tpls_maths
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs/Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: phi0
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: sign_mx
    double precision :: h

    !----- Miscellaneous -----!

    integer :: n, i, j, k
    double precision :: a1, a2, a3, a4, a5, a6, a7
    
    n = maxn

    sign_mx = 0.0D+00
    do k = 0,maxn
       do j = sy,ey
          do i = sx,ex
             sign_mx(i,j,k) = 2.*(heaviside_function(phi0(i,j,k))-.5)
          enddo
       enddo
    enddo
    
    return
  end subroutine get_levelset_status    





  subroutine boundary_conditions_levelset(phi)

    !-----
    !
    !     Apply the boundary conditions for the levelset field :
    !          + Neuman boundaries at the wall.
    !          + Classical inflow/outflow boundaries in case needed.
    !
    !     INPUT/OUTPUT
    !     ------------
    !
    !     phi : levelset field
    !
    !-----

    use tpls_constants
    use tpls_userchk
    use tpls_mpi
    use mpi
    implicit none

    !----- Input/Output -----!
    
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(inout) :: phi

    !----- Miscellaneous -----!
    
    double precision, dimension(0:maxn-2) :: u_inlet
    double precision, dimension(0:maxn) :: phi_inlet
    double precision :: xc, zc, pi, t, z_val, x_val

    integer :: i, j, k
    
    !----- Imposing the boundary conditions on phi -----!
    
    call user_inflow(u_inlet,phi_inlet)
    
    !+++ Inlet: prescribed profile +++!
    
    if ( sx == 1 ) then
       do j = sy,ey
          phi(0,j,:) = 2.0D+00*phi_inlet(:) - phi(1,j,:)
       enddo
    end if
    
    !+++ Oulet: Neumann boundary condition +++!
    
    if ( ex == ex_max ) then
       phi(ex_max+1, :, :) = phi(ex_max, :, :)
    end if
    
    phi(:,:,0) = phi(:,:,1)
    phi(:,:,maxn) = phi(:,:,maxn-1)
    
    call exchange2d(phi,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    return
  end subroutine boundary_conditions_levelset
  
  
  
  

  subroutine get_density(density,phi)
    
    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none

    !----- Input -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in) :: phi
    
    !----- Ouput -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(out) :: density

    integer :: i,j,k
  
    double precision :: heaviside_phi
    double precision :: phi_val, temp

    if (density_matched) then

       density = 1.0D+00

    else
       
       do k = 0,maxn
          do j = sy,ey
             do i = sx,ex
                
                phi_val = phi(i,j,k)
                
                heaviside_phi = Heaviside_function(phi_val)
                
                temp = rho_plus*heaviside_phi + rho_minus*(1.d0-heaviside_phi)
                density(i,j,k) = temp
                
             end do
          end do
       end do
       
       ! Initialise augmented density at boundary points.
       
       if ( sx == 1 ) then
          
          do k = 0,maxn
             do j = sy,ey
                
                phi_val = phi(0,j,k)
                
                heaviside_phi  = Heaviside_function(phi_val)
                
                temp           = rho_plus*heaviside_phi+rho_minus*(1.d0-heaviside_phi)
                density(0,j,k) = temp
                
             end do
          end do
       end if
       
       if ( ex == ex_max ) then
          
          do k = 0,maxn
             do j = sy,ey
                
                phi_val = phi(ex_max+1,j,k)
                
                heaviside_phi = Heaviside_function(phi_val)
                
                temp                  = rho_plus*heaviside_phi+rho_minus*(1.d0-heaviside_phi)
                density(ex_max+1,j,k) = temp
                
             end do
          end do
       end if
       
       call exchange2d(density,stride_p_xz,stride_p_yz, &
            neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    endif
    
    return
  end subroutine get_density
  




  subroutine get_viscosity(viscosity,phi)
    
    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none
    
    !----- Input -----!
    
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in) :: phi

    !----- Ouput -----!

    double precision, dimension(sx-2:ex+2,sy-2:ey+2,0:maxn), intent(out) :: viscosity
    
    integer :: i,j,k
    
    double precision :: heaviside_phi
    double precision :: phi_val,temp

    do k=0,maxn
       do j=sy,ey
          do i=sx,ex
             
             phi_val = phi(i,j,k)
             
             heaviside_phi = Heaviside_function(phi_val)
             
             temp = mu_plus*heaviside_phi + mu_minus*(1.d0-heaviside_phi)
             viscosity(i,j,k) = temp/Re
             
          end do
       end do
    end do
    
    ! Initialise augmented viscosity at boundary points.
    
    if ( sx == 1 ) then
       
       do k = 0,maxn
          do j = sy,ey
             
             phi_val = phi(0,j,k)
             
             heaviside_phi = Heaviside_function(phi_val)
             
             temp             = mu_plus*heaviside_phi+mu_minus*(1.d0-heaviside_phi)
             viscosity(0,j,k) = temp/Re
             
          end do
       end do
    end if
    
    if ( ex == ex_max ) then
       
       do k = 0,maxn
          do j = sy,ey
             
             phi_val = phi(ex_max+1,j,k)
             
             heaviside_phi = Heaviside_function(phi_val)
             
             temp                    = mu_plus*heaviside_phi+mu_minus*(1.d0-heaviside_phi)
             viscosity(ex_max+1,j,k) = temp/Re
             
          end do
       end do
    end if

    call exchange2d_aug1(viscosity,stride_p_aug1_xz,stride_p_aug1_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    
    call exchange2d_aug2(viscosity,stride_p_aug2_xz,stride_p_aug2_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    
    return
  end subroutine get_viscosity












  subroutine redistancing_levelset(phi)

    !-----
    !
    !     Wrapper subroutine for the redistancing procedure of the levelset field.
    !
    !     INPUT/OUTPUT
    !     ------------
    !
    !     phi : levelset field.
    !
    !-----

    use tpls_constants
    use tpls_mpi
    use mpi

    !----- Input/Output -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(inout) :: phi

    !----- Miscellaneous -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: sign_mx
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: phi_reinit

    double precision :: erreur
    integer :: iteration

    call get_levelset_status(phi,sign_mx,dx)

    phi_reinit = phi
    erreur = 1.0D+00

    do iteration = 1,max_iteration_levelset
       
       call iteration_levelset(sign_mx,phi,phi_reinit)
       
    end do
    
    phi = phi_reinit
    
    return
  end subroutine redistancing_levelset


  subroutine iteration_levelset(sign_phi0,phi0,phi)

    !-----
    !
    !     Solve the PDE-based equation for the redistancing of the levelset.
    !     Time-discretization is performed using the SSP-RK3 scheme while
    !     spatial discretization relies on a WENO5 procedure.
    !
    !     INPUTS
    !     ------
    !
    !     sign_phi : sign of phi as given by get_levelset_status.
    !     phi0     : deprecated.
    !     phi      : levelset field
    !
    !     OUTPUT
    !     ------
    !
    !     phi      : redistanced levelset field.
    !
    !-----

    use tpls_constants
    use tpls_maths
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs/Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)  :: phi0, sign_phi0
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(inout) :: phi

    !----- Miscellaneous -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: phi_n1, phi_n2, phi_n3
    
!!$    !----- RK3 -----!
!!$
!!$    call advance_redistancing(phi,phi_n1,sign_phi0)
!!$    call boundary_conditions_levelset(phi_n1)
!!$
!!$    call advance_redistancing(phi_n1,phi_n2,sign_phi0)
!!$    phi_n2 = (3./4.)*phi + (1./4.)*phi_n2
!!$    call boundary_conditions_levelset(phi_n2)
!!$
!!$    call advance_redistancing(phi_n2,phi_n3,sign_phi0)
!!$    phi = (1./3.)*phi + (2./3.)*phi_n3
!!$    call boundary_conditions_levelset(phi)
!!$    !where( abs(phi) .GE. 10*smooth_width ) phi = sign_phi0 * 10 * smooth_width

    !----- RK2 -----!

    call advance_redistancing(phi,phi_n1,sign_phi0)
    call boundary_conditions_levelset(phi_n1)

    call advance_redistancing(phi_n1,phi_n2,sign_phi0)
    call boundary_conditions_levelset(phi_n2)
    phi = .5*(phi + phi_n2)
    call boundary_conditions_levelset(phi)

    call exchange2d(phi,stride_p_xz,stride_p_yz,neighbours, &
         ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)


    return
  end subroutine iteration_levelset






  subroutine advance_redistancing(phi_in,phi_out,sign_phi0)

    !-----
    !
    !     Basic temporal step for the redistacing equation of the levelset.
    !     The equation is solved only on a band of 10 smoothing width around the interface.
    !
    !     INPUTS
    !     ------
    !
    !     phi_in   : levelset field at the pseudo time step n.
    !     sign_phi : sign of the levelset as given by get_levelset_status
    !
    !     OUTPUT
    !     ------
    !
    !     phi_out  : levelset field at the pseudo time step n+1
    !
    !-----

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(in) :: phi_in, sign_phi0
    
    !----- Output -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(out) :: phi_out

    integer :: i, j, k
    double precision :: dtau, dphidx, dphidy, dphidz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: rhs
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: A, B, C, D, E, F

    dtau = .3*min(dx,dy,dz)

    call WENO5_derivatives(phi_in,A,B,C,D,E,F)

    rhs = 0

    do k = 1,maxn-1
       do j = sy,ey
          do i = sx,ex

             !if ( abs(phi_in(i,j,k)) .LT. 10*smooth_width ) then

                !+++++ In the x-direction +++++!
                
                if ( ( B(i,j,k)*sign_phi0(i,j,k) .LT. 0.0D+00 ) &
                     .AND. ( (A(i,j,k) + B(i,j,k))*sign_phi0(i,j,k) .LT. 0.0D+00) ) then
                   
                   dphidx = B(i,j,k)
                   
                elseif ( ( A(i,j,k)*sign_phi0(i,j,k) .GT. 0.0D+00 ) &
                     .AND. ( (A(i,j,k) + B(i,j,k))*sign_phi0(i,j,k) .GT. 0.0D+00) ) then
                   
                   dphidx = A(i,j,k)
                   
                else
                   
                   dphidx = .5*(A(i,j,k) + B(i,j,k))
                   
                endif
                
                !+++++ In the y-direction +++++!
                
                if ( ( D(i,j,k)*sign_phi0(i,j,k) .LT. 0.0D+00 ) &
                     .AND. ( (C(i,j,k) + D(i,j,k))*sign_phi0(i,j,k) .LT. 0.0D+00) ) then
                   
                   dphidy = D(i,j,k)
                   
                elseif ( ( C(i,j,k)*sign_phi0(i,j,k) .GT. 0.0D+00 ) &
                     .AND. ( (C(i,j,k) + D(i,j,k))*sign_phi0(i,j,k) .GT. 0.0D+00) ) then
                   
                   dphidy = C(i,j,k)
                   
                else
                   
                   dphidy = .5*(C(i,j,k) + D(i,j,k))
                   
                endif
                
                !+++++ In the z-direction +++++!
                
                if ( ( F(i,j,k)*sign_phi0(i,j,k) .LT. 0.0D+00 ) &
                     .AND. ( (E(i,j,k) + F(i,j,k))*sign_phi0(i,j,k) .LT. 0.0D+00) ) then
                   
                   dphidz = F(i,j,k)
                   
                elseif ( ( E(i,j,k)*sign_phi0(i,j,k) .GT. 0.0D+00 ) &
                     .AND. ( (E(i,j,k) + F(i,j,k))*sign_phi0(i,j,k) .GT. 0.0D+00) ) then
                   
                   dphidz = E(i,j,k)
                   
                else
                   
                   dphidz = .5*(E(i,j,k)+D(i,j,k))
                   
                endif
                
                !+++++ Compute right hand side +++++!
                
                RHS(i,j,k) = sign_phi0(i,j,k) * ( 1.0D+00 - dsqrt( dphidx**2.0D+00 + dphidy**2.0D+00 + dphidz**2.0D+00 ) )

             !endif
             
          enddo
       enddo
    enddo
    
    phi_out = phi_in + dtau * RHS
    call boundary_conditions_levelset(phi_out)
    call exchange2d(phi_out,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

  end subroutine advance_redistancing


  subroutine WENO5_derivatives(phi_,A,B,C,D,E,F)

    !---
    !
    !   This subroutine computed the upwinding derivatives of the levelset function
    !   required for the PDE-based redistancing equation using a standard WENO5 procedure.
    !
    !---

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) , intent(in) :: phi_
    
    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(out) :: A, B, C, D, E, F

    !----- WENO-5 parameters -----!
    
    double precision :: eps_weno
    parameter ( eps_weno = 1.0D-06 )

    double precision :: sigma_1, sigma_2, sigma_3
    parameter ( sigma_1 = 0.1D+00 )
    parameter ( sigma_2 = 0.6D+00 )
    parameter ( sigma_3 = 0.3D+00 )

    double precision :: v1, v2, v3, v4, v5
    double precision :: S1, S2, S3
    double precision :: a1, a2, a3
    double precision :: w1, w2, w3

    !----- Miscellaneous -----!

    integer :: i, j, k, l, n
    double precision, dimension(sx-3:ex+3,sy-3:ey+3,0:maxn) :: phi

    n = maxn
    l = maxl

    !+++++ Creation of the augmented array +++++!

    phi = 0.0D+00
    phi(sx:ex,sy:ey,:) = phi_(sx:ex,sy:ey,:)

    if ( sx == 1 ) then
       phi(0,sy:ey,:) = phi_(0,sy:ey,:)
    endif
    if ( ex == ex_max ) then
       phi(ex_max+1,sy:ey,:) = phi_(ex_max+1,sy:ey,:)
    endif
    
    ! Augmented phi: Exchange first-, second- and third-order halos
    
    call exchange2d_augaug1(phi,stride_p_augaug1_xz,stride_p_augaug1_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)  
    
    call exchange2d_augaug2(phi,stride_p_augaug2_xz,stride_p_augaug2_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    
    call exchange2d_augaug3(phi,stride_p_augaug3_xz,stride_p_augaug3_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)

    !+++++ Computation of the WENO 5 fluxes in all three directions +++++!

    A = 0
    B = 0
    C = 0
    D = 0
    E = 0
    F = 0

    do k = 1,maxn-1
       do j = sy,ey
          do i = sx,ex

             !if ( abs(phi(i,j,k)) .LT. 10*smooth_width ) then

                !----- Computation of the fluxes in the x-direction -----!
                
                if ( ( i .GE. 4 ) .AND. ( i .LE. l-4 ) ) then
                   
                   v1 = ( phi(i-2,j,k) - phi(i-3,j,k) ) / dx
                   v2 = ( phi(i-1,j,k) - phi(i-2,j,k) ) / dx
                   v3 = ( phi(i,j,k)   - phi(i-1,j,k) ) / dx
                   v4 = ( phi(i+1,j,k) - phi(i,j,k)   ) / dx
                   v5 = ( phi(i+2,j,k) - phi(i+1,j,k) ) / dx
                   
                   S1 = (13.0D+00/12.0D+00) * ( v1 - 2.0D+00*v2 + v3 )**2.0D+00 &
                        + (1.0D+00/4.0D+00) * ( v1 - 4.0D+00*v2 + 3.0D+00*v3)**2.0D+00
                   
                   S2 = (13.0D+00/12.0D+00) * (v2 - 2.0D+00*v3 + v4)**2.0D+00 &
                        + (1.0D+00/4.0D+00) * (v2 - v4)**2.0D+00
                   
                   S3 = (13.0D+00/12.0D+00) * (v3 - 2.0D+00*v4 + v5)**2.0D+00 &
                        + (1.0D+00/4.0D+00) * (3.0D+00*v3 - 4.0D+00*v4 + v5)**2.0D+00
                   
                   a1 = sigma_1 / (eps_weno + S1)**2.0D+00
                   a2 = sigma_2 / (eps_weno + S2)**2.0D+00
                   a3 = sigma_3 / (eps_weno + S3)**2.0D+00
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   
                   A(i,j,k) = w1 * ( (1.0D+00/3.0D+00)  * v1   &
                        - (7.0D+00/6.0D+00)  * v2   &
                        + (11.0D+00/6.0D+00) * v3 ) &
                        + w2 * ( (-1.0D+00/6.0D+00) * v2   &
                        + (5.0D+00/6.0D+00)  * v3   &
                        + (1.0D+00/3.0D+00)  * v4 ) &
                        + w3 * ( (1.0D+00/3.0D+00)  * v3   &
                        + (5.0D+00/6.0D+00)  * v4   &
                        - (1.0D+00/6.0D+00)  * v5 )
                   
                   
                   v1 = ( phi(i+3,j,k) - phi(i+2,j,k) ) / dx
                   v2 = ( phi(i+2,j,k) - phi(i+1,j,k) ) / dx
                   v3 = ( phi(i+1,j,k) - phi(i,j,k)   ) / dx
                   v4 = ( phi(i,j,k)   - phi(i-1,j,k) ) / dx
                   v5 = ( phi(i-1,j,k) - phi(i-2,j,k) ) / dx
                   
                   
                   S1 = (13.0D+00/12.0D+00) * ( v1 - 2.0D+00*v2 + v3 )**2.0D+00 &
                        + (1.0D+00/4.0D+00) * ( v1 - 4.0D+00*v2 + 3.0D+00*v3)**2.0D+00
                   
                   S2 = (13.0D+00/12.0D+00) * (v2 - 2.0D+00*v3 + v4)**2.0D+00 &
                        + (1.0D+00/4.0D+00) * (v2 - v4)**2.0D+00
                   
                   S3 = (13.0D+00/12.0D+00) * (v3 - 2.0D+00*v4 + v5)**2.0D+00 &
                        + (1.0D+00/4.0D+00) * (3.0D+00*v3 - 4.0D+00*v4 + v5)**2.0D+00
                   
                   a1 = sigma_1 / (eps_weno + S1)**2.0D+00
                   a2 = sigma_2 / (eps_weno + S2)**2.0D+00
                   a3 = sigma_3 / (eps_weno + S3)**2.0D+00
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   
                   B(i,j,k) = w1 * ( (1.0D+00/3.0D+00)  * v1   &
                        - (7.0D+00/6.0D+00)  * v2   &
                        + (11.0D+00/6.0D+00) * v3 ) &
                        + w2 * ( (-1.0D+00/6.0D+00) * v2   &
                        + (5.0D+00/6.0D+00)  * v3   &
                        + (1.0D+00/3.0D+00)  * v4 ) &
                        + w3 * ( (1.0D+00/3.0D+00)  * v3   &
                        + (5.0D+00/6.0D+00)  * v4   &
                        - (1.0D+00/6.0D+00)  * v5 )
                   
                else
                   
                   A(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / dx
                   B(i,j,k) = ( phi(i+1,j,k) - phi(i,j,k) ) / dx
                   
                endif
                
                !----- Computation of the fluxes in the y-direction -----!
                
                v1 = ( phi(i,j-2,k) - phi(i,j-3,k) ) / dx
                v2 = ( phi(i,j-1,k) - phi(i,j-2,k) ) / dx
                v3 = ( phi(i,j,k)   - phi(i,j-1,k) ) / dx
                v4 = ( phi(i,j+1,k) - phi(i,j,k)   ) / dx
                v5 = ( phi(i,j+2,k) - phi(i,j+1,k) ) / dx
                
                S1 = (13.0D+00/12.0D+00) * ( v1 - 2.0D+00*v2 + v3 )**2.0D+00 &
                     + (1.0D+00/4.0D+00) * ( v1 - 4.0D+00*v2 + 3.0D+00*v3)**2.0D+00
                
                S2 = (13.0D+00/12.0D+00) * (v2 - 2.0D+00*v3 + v4)**2.0D+00 &
                     + (1.0D+00/4.0D+00) * (v2 - v4)**2.0D+00
                
                S3 = (13.0D+00/12.0D+00) * (v3 - 2.0D+00*v4 + v5)**2.0D+00 &
                     + (1.0D+00/4.0D+00) * (3.0D+00*v3 - 4.0D+00*v4 + v5)**2.0D+00
                
                a1 = sigma_1 / (eps_weno + S1)**2.0D+00
                a2 = sigma_2 / (eps_weno + S2)**2.0D+00
                a3 = sigma_3 / (eps_weno + S3)**2.0D+00
                
                w1 = a1 / (a1 + a2 + a3)
                w2 = a2 / (a1 + a2 + a3)
                w3 = a3 / (a1 + a2 + a3)
                
                C(i,j,k) = w1 * ( (1.0D+00/3.0D+00)  * v1   &
                     - (7.0D+00/6.0D+00)  * v2   &
                     + (11.0D+00/6.0D+00) * v3 ) &
                     + w2 * ( (-1.0D+00/6.0D+00) * v2   &
                     + (5.0D+00/6.0D+00)  * v3   &
                     + (1.0D+00/3.0D+00)  * v4 ) &
                     + w3 * ( (1.0D+00/3.0D+00)  * v3   &
                     + (5.0D+00/6.0D+00)  * v4   &
                     - (1.0D+00/6.0D+00)  * v5 )
                
                
                v1 = ( phi(i,j+3,k) - phi(i,j+2,k) ) / dx
                v2 = ( phi(i,j+2,k) - phi(i,j+1,k) ) / dx
                v3 = ( phi(i,j+1,k) - phi(i,j,k)   ) / dx
                v4 = ( phi(i,j,k)   - phi(i,j-1,k) ) / dx
                v5 = ( phi(i,j-1,k) - phi(i,j-2,k) ) / dx
                
                S1 = (13.0D+00/12.0D+00) * ( v1 - 2.0D+00*v2 + v3 )**2.0D+00 &
                     + (1.0D+00/4.0D+00) * ( v1 - 4.0D+00*v2 + 3.0D+00*v3)**2.0D+00
                
                S2 = (13.0D+00/12.0D+00) * (v2 - 2.0D+00*v3 + v4)**2.0D+00 &
                     + (1.0D+00/4.0D+00) * (v2 - v4)**2.0D+00
                
                S3 = (13.0D+00/12.0D+00) * (v3 - 2.0D+00*v4 + v5)**2.0D+00 &
                     + (1.0D+00/4.0D+00) * (3.0D+00*v3 - 4.0D+00*v4 + v5)**2.0D+00
                
                a1 = sigma_1 / (eps_weno + S1)**2.0D+00
                a2 = sigma_2 / (eps_weno + S2)**2.0D+00
                a3 = sigma_3 / (eps_weno + S3)**2.0D+00
                
                w1 = a1 / (a1 + a2 + a3)
                w2 = a2 / (a1 + a2 + a3)
                w3 = a3 / (a1 + a2 + a3)
                
                D(i,j,k) = w1 * ( (1.0D+00/3.0D+00)  * v1   &
                     - (7.0D+00/6.0D+00)  * v2   &
                     + (11.0D+00/6.0D+00) * v3 ) &
                     + w2 * ( (-1.0D+00/6.0D+00) * v2   &
                     + (5.0D+00/6.0D+00)  * v3   &
                     + (1.0D+00/3.0D+00)  * v4 ) &
                     + w3 * ( (1.0D+00/3.0D+00)  * v3   &
                     + (5.0D+00/6.0D+00)  * v4   &
                     - (1.0D+00/6.0D+00)  * v5 )
                
                !----- Computation of the fluxes in the z-direction -----!
                
                if ( ( k .GE. 4 ) .AND. ( k .LE. maxn-4 ) ) then
                   
                   v1 = ( phi(i,j,k-2) - phi(i,j,k-3) ) / dz
                   v2 = ( phi(i,j,k-1) - phi(i,j,k-2) ) / dz
                   v3 = ( phi(i,j,k)   - phi(i,j,k-1) ) / dz
                   v4 = ( phi(i,j,k+1) - phi(i,j,k)   ) / dz
                   v5 = ( phi(i,j,k+2) - phi(i,j,k+1) ) / dz
                   
                   S1 = (13.0D+00/12.0D+00) * ( v1 - 2.0D+00*v2 + v3 )**2.0D+00 &
                        + (1.0D+00/4.0D+00) * ( v1 - 4.0D+00*v2 + 3.0D+00*v3)**2.0D+00
                   
                   S2 = (13.0D+00/12.0D+00) * (v2 - 2.0D+00*v3 + v4)**2.0D+00 &
                        + (1.0D+00/4.0D+00) * (v2 - v4)**2.0D+00
                   
                   S3 = (13.0D+00/12.0D+00) * (v3 - 2.0D+00*v4 + v5)**2.0D+00 &
                        + (1.0D+00/4.0D+00) * (3.0D+00*v3 - 4.0D+00*v4 + v5)**2.0D+00
                   
                   a1 = sigma_1 / (eps_weno + S1)**2.0D+00
                   a2 = sigma_2 / (eps_weno + S2)**2.0D+00
                   a3 = sigma_3 / (eps_weno + S3)**2.0D+00
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   
                   E(i,j,k) = w1 * ( (1.0D+00/3.0D+00)  * v1   &
                        - (7.0D+00/6.0D+00)  * v2   &
                        + (11.0D+00/6.0D+00) * v3 ) &
                        + w2 * ( (-1.0D+00/6.0D+00) * v2   &
                        + (5.0D+00/6.0D+00)  * v3   &
                        + (1.0D+00/3.0D+00)  * v4 ) &
                        + w3 * ( (1.0D+00/3.0D+00)  * v3   &
                        + (5.0D+00/6.0D+00)  * v4   &
                        - (1.0D+00/6.0D+00)  * v5 )
                   
                   
                   v1 = ( phi(i,j,k+3) - phi(i,j,k+2) ) / dz
                   v2 = ( phi(i,j,k+2) - phi(i,j,k+1) ) / dz
                   v3 = ( phi(i,j,k+1) - phi(i,j,k)   ) / dz
                   v4 = ( phi(i,j,k)   - phi(i,j,k-1) ) / dz
                   v5 = ( phi(i,j,k-1) - phi(i,j,k-2) ) / dz
                   
                   S1 = (13.0D+00/12.0D+00) * ( v1 - 2.0D+00*v2 + v3 )**2.0D+00 &
                        + (1.0D+00/4.0D+00) * ( v1 - 4.0D+00*v2 + 3.0D+00*v3)**2.0D+00
                   
                   S2 = (13.0D+00/12.0D+00) * (v2 - 2.0D+00*v3 + v4)**2.0D+00 &
                        + (1.0D+00/4.0D+00) * (v2 - v4)**2.0D+00
                   
                   S3 = (13.0D+00/12.0D+00) * (v3 - 2.0D+00*v4 + v5)**2.0D+00 &
                        + (1.0D+00/4.0D+00) * (3.0D+00*v3 - 4.0D+00*v4 + v5)**2.0D+00
                   
                   a1 = sigma_1 / (eps_weno + S1)**2.0D+00
                   a2 = sigma_2 / (eps_weno + S2)**2.0D+00
                   a3 = sigma_3 / (eps_weno + S3)**2.0D+00
                   
                   w1 = a1 / (a1 + a2 + a3)
                   w2 = a2 / (a1 + a2 + a3)
                   w3 = a3 / (a1 + a2 + a3)
                   
                   F(i,j,k) = w1 * ( (1.0D+00/3.0D+00)  * v1   &
                        - (7.0D+00/6.0D+00)  * v2   &
                        + (11.0D+00/6.0D+00) * v3 ) &
                        + w2 * ( (-1.0D+00/6.0D+00) * v2   &
                        + (5.0D+00/6.0D+00)  * v3   &
                        + (1.0D+00/3.0D+00)  * v4 ) &
                        + w3 * ( (1.0D+00/3.0D+00)  * v3   &
                        + (5.0D+00/6.0D+00)  * v4   &
                        - (1.0D+00/6.0D+00)  * v5 )
                   
                else
                   
                   E(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dz
                   F(i,j,k) = ( phi(i,j,k+1) - phi(i,j,k) ) / dz
                   
                endif

             !endif
             
          enddo
       enddo
    enddo
    
    return
  end subroutine WENO5_derivatives





  subroutine HOUC5(dfx,dfy,dfz,phi_,vx,vy,vz)

    !---
    !
    !   This subroutine computed the upwinding derivatives of the levelset function using
    !   the High Order Upstream Central scheme of order 5. For further details please read:
    !   ???
    !
    !---

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   , intent(in) :: phi_
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) , intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   , intent(in) :: vz
    
    !----- Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(out) :: dfx, dfy, dfz

    !----- Miscellaneous -----!

    integer :: i, j, k, l, n
    double precision, dimension(sx-3:ex+3,sy-3:ey+3,sz_p:ez_p) :: phi

    n = maxn
    l = maxl

    dfx = 0.0D+00
    dfy = 0.0D+00
    dfz = 0.0D+00

    !+++++ Creation of the augmented array +++++!

    phi = 0.0D+00
    phi(sx:ex,sy:ey,:) = phi_(sx:ex,sy:ey,:)

    if ( sx == 1 ) then
       phi(0,sy:ey,:) = phi_(0,sy:ey,:)
    endif
    if ( ex == ex_max ) then
       phi(ex_max+1,sy:ey,:) = phi_(ex_max+1,sy:ey,:)
    endif
    
    ! Augmented phi: Exchange first-, second- and third-order halos
    
    call exchange2d_augaug1(phi,stride_p_augaug1_xz,stride_p_augaug1_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)  
    
    call exchange2d_augaug2(phi,stride_p_augaug2_xz,stride_p_augaug2_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    
    call exchange2d_augaug3(phi,stride_p_augaug3_xz,stride_p_augaug3_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)

    !+++++ Computation of the WENO 5 fluxes in all three directions +++++!

    do k = 1,maxn-1
       do j = sy,ey
          do i = sx,ex

             !----- Computation of the fluxes in the x-direction -----!


             if ( ( i .GE. 4 ) .AND. ( i .LE. l-4 ) ) then

                if ( (vx(i+1,j,k-1)+vx(i,j,k-1)) .GE. 0.0D+00 ) then

                   dfx(i,j,k) = -2.*phi(i-3,j,k) + 15.*phi(i-2,j,k) &
                        - 60.*phi(i-1,j,k) + 20.*phi(i,j,k) + 30.*phi(i+1,j,k) -3.*phi(i+2,j,k)

                else

                   dfx(i,j,k) = 2.*phi(i+3,j,k) - 15.*phi(i+2,j,k) &
                        + 60.*phi(i+1,j,k) - 20.*phi(i,j,k) -30.*phi(i-1,j,k) + 3.*phi(i-2,j,k)
                   
                endif

                dfx(i,j,k) = dfx(i,j,k)/(60.*dx)

             else

                if ( (vx(i,j,k-1)+vx(i+1,j,k-1)).LT.0.0D+00) then
                   dfx(i,j,k) = ( phi(i+1,j,k) - phi(i,j,k) ) / dx
                else
                   dfx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / dx
                endif

             endif
             
             !----- Computation of the fluxes in the y-direction -----!
             
             if ( (vy(i,j,k-1)+vy(i,j+1,k-1)) .GE. 0.0D+00 ) then
                
                dfy(i,j,k) = -2.*phi(i,j-3,k) + 15.*phi(i,j-2,k) &
                     - 60.*phi(i,j-1,k) + 20.*phi(i,j,k) + 30.*phi(i,j+1,k) - 3.*phi(i,j+2,k)
                
             else

                dfy(i,j,k) = 2.*phi(i,j+3,k) - 15.*phi(i,j+2,k) &
                     + 60.*phi(i,j+1,k) - 20.*phi(i,j,k) - 30.*phi(i,j-1,k) + 3.*phi(i,j-2,k)
                
             endif

             dfy(i,j,k) = dfy(i,j,k)/(60.*dy)

             !----- Computation of the fluxes in the z-direction -----!
             
             if ( ( k .GE. 4 ) .AND. ( k .LE. maxn-4 ) ) then
                
                if ( ( vz(i,j,k)+vz(i,j,k-1) ) .GE. 0.0D+00 ) then

                dfz(i,j,k) = -2.*phi(i,j,k-3) + 15.*phi(i,j,k-2) &
                     - 60.*phi(i,j,k-1) + 20.*phi(i,j,k) + 30.*phi(i,j,k+1) - 3.*phi(i,j,k+2)                   

                else

                dfz(i,j,k) = 2.*phi(i,j,k+3) - 15.*phi(i,j,k+2) &
                     + 60.*phi(i,j,k+1) - 20.*phi(i,j,k) - 30.*phi(i,j,k-1) + 3.*phi(i,j,k-2)

                endif
                
                dfz(i,j,k) = dfz(i,j,k)/(60.*dz)

             else

                if ( (vz(i,j,k)+vz(i,j,k-1)).LT.0.0D+00 ) then
                   dfz(i,j,k) = ( phi(i,j,k+1) - phi(i,j,k) ) / dz
                else
                   dfz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dz 
                endif

             endif



          enddo
       enddo
    enddo
    
    return
  end subroutine HOUC5



end module tpls_levelset
