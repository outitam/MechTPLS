module tpls_levelset

  implicit none

contains

  
  subroutine advance_levelset(phi,vx,vy,vz)

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
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: phi1, phi2, phi3

    t_temp = mpi_wtime()

    !----- Advection of the levelset function using RK-2 -----!


    call convective_terms_phi(conv2_phi,phi,vx,vy,vz)
    phi1 = phi - dt * conv2_phi
    call boundary_conditions_levelset(phi1)
    call exchange2d(phi1,stride_p_xz,stride_p_yz,neighbours, &
         ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    call convective_terms_phi(conv2_phi,phi1,vx,vy,vz)
    phi2 = phi1 - dt*conv2_phi
    call boundary_conditions_levelset(phi2)
    call exchange2d(phi2,stride_p_xz,stride_p_yz,neighbours, &
         ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    phi = (1.0D+00/2.0D+00)*(phi + phi2)
    call boundary_conditions_levelset(phi)
    call exchange2d(phi,stride_p_xz,stride_p_yz,neighbours, &
         ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    !----- Re-initialization to the signed distance function -----
    
    call redistancing_levelset(phi)

    !----- Miscellaneous -----!

    t_temp2 = mpi_wtime()
    time_levelset = time_levelset + (t_temp2 - t_temp)

    return
  end subroutine advance_levelset





  subroutine convective_terms_phi(conv2_phi,phi,vx,vy,vz)

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

    call WENO5_FD(phi,vx,vy,vz,dfx_,dfy_,dfz_)

    conv2_phi(sx:ex,sy:ey,1:maxn-1) = (vx(sx+1:ex+1,sy:ey,0:maxn-2)+vx(sx:ex,sy:ey,0:maxn-2))/2.0D+00 * dfx_(sx:ex,sy:ey,1:maxn-1) &
         + (vy(sx:ex,sy+1:ey+1,0:maxn-2)+vy(sx:ex,sy:ey,0:maxn-2))/2.0D+00 * dfy_(sx:ex,sy:ey,1:maxn-1) &
         + (vz(sx:ex,sy:ey,1:maxn-1)+vz(sx:ex,sy:ey,0:maxn-2))/2.0D+00 * dfz_(sx:ex,sy:ey,1:maxn-1)

    call exchange2d(conv2_phi,stride_p_xz,stride_p_yz,neighbours, &
         ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    return
  end subroutine convective_terms_phi




  subroutine WENO5_FD(phi_,vx,vy,vz,dfx,dfy,dfz)

    !---
    !
    !   This subroutine computed the upwinding convective term for the levelset function.
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


             if ( ( i .GE. 3 ) .AND. ( i .LE. l-3 ) ) then

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
!!$
!!$                w1 = sigma_1
!!$                w2 = sigma_2
!!$                w3 = sigma_3

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
!!$                w1 = sigma_1
!!$                w2 = sigma_2
!!$                w3 = sigma_3
             
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
             
             if ( ( k .GE. 3 ) .AND. ( k .LE. maxn-3 ) ) then
                
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
!!$
!!$                w1 = sigma_1
!!$                w2 = sigma_2
!!$                w3 = sigma_3

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

    call exchange2d(dfx,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(dfy,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(dfz,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    return
  end subroutine WENO5_FD

  



  
  subroutine get_levelset_status(phi0,D_mx,istatus,sign_mx,h)

    !---
    !
    !   Re-initialization step following Russo & Smereka, JCP 163, 51
    !
    !---

    use tpls_constants
    use tpls_maths
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs/Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: phi0
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: sign_mx
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: D_mx
    integer, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn)          :: istatus
    double precision :: h

    !----- Miscellaneous -----!

    integer :: n, i, j, k
    double precision :: a1, a2, a3, a4, a5, a6, a7
    
    n = maxn

    D_mx    = 0.0D+00
    sign_mx = 0.0D+00
    istatus = 0
     
    do k = 1,n-1
       do j = sy,ey
          do i = sx,ex
             
             ! Compute sign function of phi0

!!$             if (phi0(i,j,k).lt. 0.0D+00) then
!!$                sign_mx(i,j,k) = -1.0D+00
!!$             elseif (phi0(i,j,k).gt.0.0D+00) then
!!$                sign_mx(i,j,k) = 1.0D+00
!!$             else
!!$                sign_mx(i,j,k) = 0.0D+00
!!$             endif
!!$
             if ( phi0(i,j,k) .LT. -dx ) then
                sign_mx(i,j,k) = -1.0D+00
             elseif ( phi0(i,j,k) .GT. dx ) then
                sign_mx(i,j,k) = 1.0D+00
             else
                sign_mx(i,j,k) = phi0(i,j,k)/sqrt(phi0(i,j,k)**2. + dx**2.)
             endif
!!$             sign_mx(i,j,k) = 2.*(Heaviside_function(phi0(i,j,k)) - 1./2.)
             
             ! Compute indicator function \Sigma_{\Delta x} (=istatus)
             
             a1 = phi0(i,j,k)*phi0(i+1,j,k)
             a2 = phi0(i,j,k)*phi0(i-1,j,k)
             a3 = phi0(i,j,k)*phi0(i,j+1,k)
             a4 = phi0(i,j,k)*phi0(i,j-1,k)  
             a5 = phi0(i,j,k)*phi0(i,j,k+1)
             a6 = phi0(i,j,k)*phi0(i,j,k-1)
             
             if (dmin1(a1,a2,a3,a4,a5,a6).le.0.0D+00) then
                istatus(i,j,k) = 1
             else
                istatus(i,j,k) = 0
             endif
             
             ! Compute distance function (D_mx)
             
             if (istatus(i,j,k).eq.1) then
                a1 = dsqrt( (phi0(i+1,j,k)-phi0(i-1,j,k))**2.0D+00 &
                     + (phi0(i,j+1,k)-phi0(i,j-1,k))**2.0D+00 &
                     + (phi0(i,j,k+1)-phi0(i,j,k-1))**2.0D+00 ) /2.0D+00
                a2 = dabs(phi0(i+1,j,k)-phi0(i,j,k))
                a3 = dabs(phi0(i-1,j,k)-phi0(i,j,k))
                a4 = dabs(phi0(i,j+1,k)-phi0(i,j,k))
                a5 = dabs(phi0(i,j-1,k)-phi0(i,j,k))
                a6 = dabs(phi0(i,j,k+1)-phi0(i,j,k))
                a7 = dabs(phi0(i,j,k-1)-phi0(i,j,k))

                D_mx(i,j,k) = h*phi0(i,j,k)/dmax1(a1,a2,a3,a4,a5,a6,a7,h)

             endif
             
          end do
       end do
    end do
    
    return
  end subroutine get_levelset_status    





  subroutine boundary_conditions_levelset(phi)

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
       phi(ex_max+1,:,:) = phi(ex_max,:,:)
    end if

    !+++ Linear extrapolation for the levelset on the walls +++!

    phi(:,:,0) = 3.0D+00*phi(:,:,1) &
         - 3.0D+00*phi(:,:,2) &
         + phi(:,:,3)
    
    phi(:,:,maxn) = 3.0D+00*phi(:,:,maxn-1) &
         - 3.0D+00*phi(:,:,maxn-2) &
         + phi(:,:,maxn-3)
!!$
!!$    
!!$    phi(:,:,0)    = 2.d0*phi(:,:,1)      - phi(:,:,2)
!!$    phi(:,:,maxn) = 2.d0*phi(:,:,maxn-1) - phi(:,:,maxn-2)
!!$
!!$    !----- Boundary conditions for Zalesk disk -----!
!!$
!!$    pi = 4.0D+00*atan(1.0D+00)
!!$    t = dt*istep
!!$
!!$    xc = 50.0D+00 + 25.0D+00*cos(pi*t/314.0D+00 + pi/2.0D+00)
!!$    zc = 50.0D+00 + 25.0D+00*sin(pi*t/314.0D+00 + pi/2.0D+00)
!!$
!!$    if ( sx == 1 ) then
!!$
!!$       do k = 0,maxn
!!$          z_val = k*dz - dz/2.0D+00
!!$          x_val = dx/2.0D+00
!!$          phi(1,:,k) = sqrt( (x_val-xc)**2.0D+00 + (z_val-zc)**2.0D+00 ) - 15.0D+00
!!$
!!$          x_val = -dx/2.0D+00
!!$          phi(0,:,k) = sqrt( (x_val-xc)**2.0D+00 + (z_val-zc)**2.0D+00 ) - 15.0D+00
!!$       enddo
!!$    endif
!!$
!!$    if ( ex == ex_max ) then
!!$
!!$       do k = 0,maxn
!!$          z_val = k*dz - dz/2.0D+00
!!$          x_val = (maxl-1)*dx - dx/2.0D+00
!!$          phi(ex_max,:,k) = sqrt( (x_val-xc)**2.0D+00 + (z_val-zc)**2.0D+00 ) - 15.0D+00
!!$          phi(ex_max+1,:,k) = 2.0D+00*phi(ex_max,:,k) - phi(ex_max-1,:,k)
!!$       enddo
!!$
!!$    endif
!!$
!!$    do i = sx,ex
!!$
!!$       x_val = i*dx - dx/2.0D+00
!!$
!!$       z_val = dz/2.0D+00
!!$       phi(i,:,1) = sqrt( (x_val-xc)**2. + (z_val-zc)**2. ) - 15.0D+00
!!$       z_val = -dz/2.0D+00
!!$       phi(i,:,0) = sqrt( (x_val-xc)**2. + (z_val-zc)**2. ) - 15.0D+00
!!$
!!$       z_val = maxn*dz - dz/2.0D+00
!!$       phi(i,:,maxn-1) = sqrt( (x_val-xc)**2. + (z_val-zc)**2. ) - 15.0D+00
!!$       z_val = maxn*dz + dz/2.0D+00
!!$       phi(i,:,maxn) = sqrt( (x_val-xc)**2. + (z_val-zc)**2. ) - 15.0D+00
!!$
!!$    enddo
!!$       
!!$    !    call mpi_barrier(comm2d_quasiperiodic,ierr)  
    
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
    
    call mpi_barrier(comm2d_quasiperiodic,ierr) 
    
    ! Exchange second-order halos 
    
    call exchange2d_aug2(viscosity,stride_p_aug2_xz,stride_p_aug2_yz, &
         neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
    
    return
  end subroutine get_viscosity





  subroutine iteration_levelset(phi0,phi_reinit,D_mx,istatus,sign_mx,h,error)
    
    use tpls_constants
    use tpls_maths
    use tpls_mpi
    use tpls_userchk
    use mpi
    implicit none

    !----- Inputs/Outputs -----!

    integer, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)             :: istatus

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)    :: phi0, sign_mx, D_mx
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(inout) :: phi_reinit

    double precision, intent(out)                                          :: error
    double precision, intent(in)                                           :: h


    !----- Miscellaneous -----!
        
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: phi_reinit_updated
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: phi_n, phi_n1, phi_n2
    double precision :: dtau, dt_
    double precision :: a1, a2, a3, a4, a5, a6, a7
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: DC
    double precision :: A, B, C, D, E, F
    double precision :: aP,aM
    double precision :: bP,bM
    double precision :: cP,cM
    double precision :: dP,dM
    double precision :: eP,eM
    double precision :: fP,fM
    double precision :: z,y
    double precision :: pi
    integer          :: i, j, k, n
    double precision :: error_val
    double precision :: laplacian, laplacian_m, laplacian_p
    double precision :: dxp, dyp, dzp, dxm, dym, dzm
    double precision :: discriminant, phi0xx, epsilon, dummy
    double precision :: dirac
    
    pi = 4.D0*atan(1.D0)

    n = maxn
    epsilon = dx/3.0
    
    ! Fictitious time
    dt_ = 0.3D+00 * dmin1(dx,dy,dz)

    phi_n = phi_reinit
    
    do k = 1,n-1
       do j = sy,ey
          do i = sx,ex

             if ( istatus(i,j,k) == 1 ) then

                DC(i,j,k) = sign_mx(i,j,k) * dabs( phi_n(i,j,k) ) - D_mx(i,j,k)
                DC(i,j,k) = DC(i,j,k)/dx


             else

                !----- Along x-direction -----!
                
                A = ( phi_n(i,j,k) - phi_n(i-1,j,k) ) / dx
                B = ( phi_n(i+1,j,k) - phi_n(i,j,k) ) / dx

!!$                if ( (i .GT. sx+2 ) .AND. ( i .LT. ex-2 ) ) then
!!$                   
!!$                   laplacian_m = ( phi_n(i-2,j,k) - 2.0D+00*phi_n(i-1,j,k) + phi_n(i,j,k) ) / dx**2.0D+00
!!$                   laplacian   = ( phi_n(i-1,j,k) - 2.0D+00*phi_n(i,j,k) + phi_n(i+1,j,k) ) / dx**2.0D+00
!!$                   laplacian_p = ( phi_n(i,j,k) - 2.0D+00*phi_n(i+1,j,k) + phi_n(i+2,j,k) ) / dx**2.0D+00
!!$
!!$                   A = A + (dx/2.0D+00) * minmod( laplacian , laplacian_m )
!!$                   B = B - (dx/2.0D+00) * minmod( laplacian_p , laplacian )
!!$
!!$                endif
                
                !----- Along y-direction -----!
                
                C = ( phi_n(i,j,k) - phi_n(i,j-1,k) ) / dy
                D = ( phi_n(i,j+1,k) - phi_n(i,j,k) ) / dy

!!$                if ( ( j .GT. sy+2 ) .AND. ( j .LT. ey-2 ) ) then
!!$                   
!!$                   laplacian_m = ( phi_n(i,j-2,k) - 2.0D+00*phi_n(i,j-1,k) + phi_n(i,j,k) ) / dy**2.0D+00
!!$                   laplacian   = ( phi_n(i,j-1,k) - 2.0D+00*phi_n(i,j,k) + phi_n(i,j+1,k) ) / dy**2.0D+00
!!$                   laplacian_p = ( phi_n(i,j,k) - 2.0D+00*phi_n(i,j+1,k) + phi_n(i,j+2,k) ) / dy**2.0D+00
!!$
!!$                   C = C + (dy/2.0D+00) * minmod( laplacian , laplacian_m )
!!$                   D = D - (dy/2.0D+00) * minmod( laplacian_p , laplacian )
!!$
!!$                endif

                                
                !----- Along z-direction -----!
                
                E = ( phi_n(i,j,k) - phi_n(i,j,k-1) ) / dz
                F = ( phi_n(i,j,k+1) - phi_n(i,j,k) ) / dz

!!$                if ( ( k .GT. sz_p+2 ) .AND. ( k .LT. ez_p-2 ) ) then
!!$                   
!!$                   laplacian_m = ( phi_n(i,j,k-2) - 2.0D+00*phi_n(i,j,k-1) + phi_n(i,j,k) ) / dz**2.0D+00
!!$                   laplacian   = ( phi_n(i,j,k-1) - 2.0D+00*phi_n(i,j,k) + phi_n(i,j,k+1) ) / dz**2.0D+00
!!$                   laplacian_p = ( phi_n(i,j,k) - 2.0D+00*phi_n(i,j,k+1) + phi_n(i,j,k+2) ) / dz**2.0D+00
!!$
!!$                   E = E + (dz/2.0D+00) * minmod( laplacian , laplacian_m )
!!$                   F = F - (dz/2.0D+00) * minmod( laplacian_p , laplacian )
!!$
!!$                endif


                if ( phi0(i,j,k) .GE. 0.0D+00 ) then
                   
                   DC(i,j,k) = dsqrt( dmax1( dmax1(A,0.0D+00)**2.0D+00 , dmin1(B,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(C,0.0D+00)**2.0D+00 , dmin1(D,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(E,0.0D+00)**2.0D+00 , dmin1(F,0.0D+00)**2.0D+00 ) ) &
                        - 1.0D+00

                   DC(i,j,k) = DC(i,j,k) * sign_mx(i,j,k)
                   
                elseif ( phi0(i,j,k) .LT. 0.0D+00 ) then
                   
                   DC(i,j,k) = dsqrt( dmax1( dmax1(B,0.0D+00)**2.0D+00 , dmin1(A,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(D,0.0D+00)**2.0D+00 , dmin1(C,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(F,0.0D+00)**2.0D+00 , dmin1(E,0.0D+00)**2.0D+00 ) ) &
                        - 1.0D+00
                   
                   DC(i,j,k) = DC(i,j,k) * sign_mx(i,j,k)
                   
                endif
                
             endif

          enddo
       enddo
    enddo

    phi_n1 = phi_n - dt_*DC

!!$    if ( sx == 1 ) then
!!$       phi_n1(0,:,:) = 3.0D+00*phi_n1(1,:,:) &
!!$            - 3.0D+00/2.0D+00*phi_n1(2,:,:) &
!!$            + phi_n1(3,:,:)
!!$    endif
!!$    
!!$    if ( ex == ex_max ) then
!!$       phi_n1(ex_max+1,:,:) = 3.0D+00*phi_n1(ex_max,:,:) &
!!$            - 3.0D+00/2.0D+00*phi_n1(ex_max-1,:,:) &
!!$            + phi_n1(ex_max-2,:,:)
!!$    endif
!!$
!!$    phi_n1(:,:,0) = 3.0D+00*phi_n1(:,:,1) &
!!$         - 3.0D+00/2.0D+00*phi_n1(:,:,2) &
!!$         + phi_n1(:,:,3)
!!$
!!$    phi_n1(:,:,maxn) = 3.0D+00*phi_n1(:,:,maxn-1) &
!!$         - 3.0D+00/2.0D+00*phi_n1(:,:,maxn-2) &
!!$         + phi_n1(:,:,maxn-3)

    call boundary_conditions_levelset(phi_n1)
    
    call exchange2d(phi_n1,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    phi_reinit_updated = phi_n1
    call exchange2d(phi_reinit_updated,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    error_val = 0.0D+00
    
    do k = 0,n
       do j = sy,ey
          do i = sx,ex
             dirac = Dirac_function(phi_reinit_updated(i,j,k))
             error_val = error_val + dirac*( phi_reinit(i,j,k)-phi_reinit_updated(i,j,k) )**2.0D+00
          enddo
       enddo
    enddo
    
    phi_reinit = phi_reinit_updated
    
    call exchange2d(phi_reinit,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    call mpi_allreduce(error_val, error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm2d_quasiperiodic, ierr)
    error = error*dx*dy*dz/(Lx*Ly*Lz)
    error = dsqrt(error)/dt_

    return
  end subroutine iteration_levelset







  subroutine redistancing_levelset(phi)

    use tpls_constants
    use tpls_mpi
    use mpi

    !----- Input/Output -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), intent(inout) :: phi

    !----- Miscellaneous -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: D_mx, sign_mx,sign_mx_bis
    integer, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)            :: istatus_mx, istatus_bis
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: phi_reinit

    double precision :: erreur
    integer :: iteration

    call get_levelset_status(phi,D_mx,istatus_mx,sign_mx,dx)

    phi_reinit = phi
    erreur = 1.0D+00

    do iteration = 1,max_iteration_levelset
       
!!$       call iteration_levelset(phi,phi_reinit,D_mx,istatus_mx,sign_mx,dx,erreur)
!!$       call iteration_levelset_WENO5(phi,phi_reinit,D_mx,istatus_mx,sign_mx,dx,erreur)
       call iteration_levelset_WENO5_SussmanFatemi(sign_mx,phi,phi_reinit)
!!$       if(erreur.LT.tolerance_levelset) EXIT Iteration_phi
       
    end do
    
    ! Now update phi3 to its redistanced value.
    
    phi = phi_reinit
    
    call exchange2d(phi,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    if ( my_id == master_id ) then
       write(*,*) 'Iteration : ', istep, 'phi-residual is ', erreur, '(', iteration, ')'
    endif
    
    return
  end subroutine redistancing_levelset

  subroutine iteration_levelset_WENO5_SussmanFatemi(sign_phi0,phi0,phi)

    use tpls_constants
    use tpls_maths
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs/Outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)  :: phi0, sign_phi0
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(inout) :: phi

    !----- Miscellaneous -----!

    double precision :: dtau, dphidx, dphidy, dphidz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: phi_n1, phi_n2, rhs
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: A, B, C, D, E, F

    integer :: i, j, k

    dtau = 0.3D+00 * dmin1(dx,dy,dz)

    RHS = 0.0D+00

    !----- RK 2 -----!
    
    call WENO5_derivatives(phi,A,B,C,D,E,F)

    do k = 1,maxn-1
       do j = sy,ey
          do i = sx,ex

             !+++++ In the x-direction +++++!

             if ( ( B(i,j,k)*sign_phi0(i,j,k) .LT. 0.0D+00 ) &
                  .AND. ( (A(i,j,k) + B(i,j,k))*sign_phi0(i,j,k) .LT. 0.0D+00) ) then

                dphidx = B(i,j,k)

             elseif ( ( A(i,j,k)*sign_phi0(i,j,k) .GT. 0.0D+00 ) &
                  .AND. ( (A(i,j,k) + B(i,j,k))*sign_phi0(i,j,k) .GT. 0.0D+00) ) then
                
                dphidx = A(i,j,k)
                
             else

                dphidx = 0.0D+00

             endif

             !+++++ In the y-direction +++++!

             if ( ( D(i,j,k)*sign_phi0(i,j,k) .LT. 0.0D+00 ) &
                  .AND. ( (C(i,j,k) + D(i,j,k))*sign_phi0(i,j,k) .LT. 0.0D+00) ) then

                dphidy = D(i,j,k)

             elseif ( ( C(i,j,k)*sign_phi0(i,j,k) .GT. 0.0D+00 ) &
                  .AND. ( (C(i,j,k) + D(i,j,k))*sign_phi0(i,j,k) .GT. 0.0D+00) ) then
                
                dphidy = C(i,j,k)
                
             else

                dphidy = 0.0D+00

             endif

             !+++++ In the z-direction +++++!

             if ( ( F(i,j,k)*sign_phi0(i,j,k) .LT. 0.0D+00 ) &
                  .AND. ( (E(i,j,k) + F(i,j,k))*sign_phi0(i,j,k) .LT. 0.0D+00) ) then

                dphidz = F(i,j,k)

             elseif ( ( E(i,j,k)*sign_phi0(i,j,k) .GT. 0.0D+00 ) &
                  .AND. ( (E(i,j,k) + F(i,j,k))*sign_phi0(i,j,k) .GT. 0.0D+00) ) then
                
                dphidz = E(i,j,k)
                
             else

                dphidz = 0.0D+00

             endif

             !+++++ Compute right hand side +++++!

             RHS(i,j,k) = sign_phi0(i,j,k) * ( 1.0D+00 - dsqrt( dphidx**2.0D+00 + dphidy**2.0D+00 + dphidz**2.0D+00 ) )
!!$             RHS(i,j,k) = -(dphidx**2./sqrt(dphidx**2. + dphidy**2. + dphidz**2.) &
!!$                  + dphidy**2./sqrt(dphidx**2. + dphidy**2. + dphidz**2.) &
!!$                  + dphidz**2./sqrt(dphidx**2. + dphidy**2. + dphidz**2.) ) + 1.
!!$             RHS(i,j,k) = sign_phi0(i,j,k)*RHS(i,j,k)

          enddo
       enddo
    enddo

    phi_n1 = phi + dt * RHS
    call boundary_conditions_levelset(phi_n1)
    call exchange2d(phi_n1,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    call WENO5_derivatives(phi_n1,A,B,C,D,E,F)

    do k = 1,maxn-1
       do j = sy,ey
          do i = sx,ex

             !+++++ In the x-direction +++++!

             if ( ( B(i,j,k)*sign_phi0(i,j,k) .LT. 0.0D+00 ) &
                  .AND. ( (A(i,j,k) + B(i,j,k))*sign_phi0(i,j,k) .LT. 0.0D+00) ) then

                dphidx = B(i,j,k)

             elseif ( ( A(i,j,k)*sign_phi0(i,j,k) .GT. 0.0D+00 ) &
                  .AND. ( (A(i,j,k) + B(i,j,k))*sign_phi0(i,j,k) .GT. 0.0D+00) ) then
                
                dphidx = A(i,j,k)
                
             else

                dphidx = 0.0D+00

             endif

             !+++++ In the y-direction +++++!

             if ( ( D(i,j,k)*sign_phi0(i,j,k) .LT. 0.0D+00 ) &
                  .AND. ( (C(i,j,k) + D(i,j,k))*sign_phi0(i,j,k) .LT. 0.0D+00) ) then

                dphidy = D(i,j,k)

             elseif ( ( C(i,j,k)*sign_phi0(i,j,k) .GT. 0.0D+00 ) &
                  .AND. ( (C(i,j,k) + D(i,j,k))*sign_phi0(i,j,k) .GT. 0.0D+00) ) then
                
                dphidy = C(i,j,k)
                
             else

                dphidy = 0.0D+00

             endif

             !+++++ In the z-direction +++++!

             if ( ( F(i,j,k)*sign_phi0(i,j,k) .LT. 0.0D+00 ) &
                  .AND. ( (E(i,j,k) + F(i,j,k))*sign_phi0(i,j,k) .LT. 0.0D+00) ) then

                dphidz = F(i,j,k)

             elseif ( ( F(i,j,k)*sign_phi0(i,j,k) .GT. 0.0D+00 ) &
                  .AND. ( (E(i,j,k) + F(i,j,k))*sign_phi0(i,j,k) .GT. 0.0D+00) ) then
                
                dphidz = E(i,j,k)
                
             else

                dphidz = 0.0D+00

             endif

             !+++++ Compute right hand side +++++!

!!$             RHS(i,j,k) = -(dphidx**2./sqrt(dphidx**2. + dphidy**2. + dphidz**2.) &
!!$                  + dphidy**2./sqrt(dphidx**2. + dphidy**2. + dphidz**2.) &
!!$                  + dphidz**2./sqrt(dphidx**2. + dphidy**2. + dphidz**2.) ) + 1.
!!$             RHS(i,j,k) = sign_phi0(i,j,k)*RHS(i,j,k)
             RHS(i,j,k) = sign_phi0(i,j,k) * ( 1.0D+00 - dsqrt( dphidx**2.0D+00 + dphidy**2.0D+00 + dphidz**2.0D+00 ) )

          enddo
       enddo
    enddo

    phi_n2 = phi_n1 + dt * RHS
    call boundary_conditions_levelset(phi_n2)
    call exchange2d(phi_n2,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    phi = (phi + phi_n2)/2.0D+00
    call boundary_conditions_levelset(phi)
    call exchange2d(phi,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)    

    return
  end subroutine iteration_levelset_WENO5_SussmanFatemi




  subroutine iteration_levelset_WENO5(phi0,phi_reinit,D_mx,istatus,sign_mx,h,error)
    
    use tpls_constants
    use tpls_maths
    use tpls_mpi
    use tpls_userchk
    use mpi
    implicit none

    !----- Inputs/Outputs -----!

    integer, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)             :: istatus

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(inout)    :: phi0, sign_mx, D_mx
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(inout) :: phi_reinit

    double precision, intent(out)                                          :: error
    double precision, intent(in)                                           :: h


    !----- Miscellaneous -----!
        
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: phi_reinit_updated
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: phi_n, phi_n1, phi_n2
    double precision :: dtau, dt_
    double precision :: a1, a2, a3, a4, a5, a6, a7
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: DC

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: A_, B_, C_, D_, E_, F_

    double precision :: A, B, C, D, E, F
    double precision :: aP,aM
    double precision :: bP,bM
    double precision :: cP,cM
    double precision :: dP,dM
    double precision :: eP,eM
    double precision :: fP,fM
    double precision :: z,y
    double precision :: pi
    integer          :: i, j, k, n
    double precision :: error_val
    double precision :: laplacian, laplacian_m, laplacian_p
    double precision :: dxp, dyp, dzp, dxm, dym, dzm
    double precision :: discriminant, phi0xx, epsilon, dummy
    double precision :: dirac
    
    pi = 4.D0*atan(1.D0)

    n = maxn
    epsilon = dx/3.0
    
    ! Fictitious time
    dt_ = 0.3D+00 * dmin1(dx,dy,dz)

    phi_n = phi_reinit

    call WENO5_derivatives(phi_n,A_,B_,C_,D_,E_,F_)
    
    do k = 1,n-1
       do j = sy,ey
          do i = sx,ex

             if ( istatus(i,j,k) == 1 ) then

                DC(i,j,k) = sign_mx(i,j,k) * dabs( phi_n(i,j,k) ) - D_mx(i,j,k)
                DC(i,j,k) = DC(i,j,k)/dx


             else

                !----- Along x-direction -----!
                
                A = A_(i,j,k)
                B = B_(i,j,k)
                
                !----- Along y-direction -----!
                
                C = C_(i,j,k)
                D = D_(i,j,k)
                                
                !----- Along z-direction -----!
                
                E = E_(i,j,k)
                F = F_(i,j,k)


                if ( phi0(i,j,k) .GE. 0.0D+00 ) then
                   
                   DC(i,j,k) = dsqrt( dmax1( dmax1(A,0.0D+00)**2.0D+00 , dmin1(B,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(C,0.0D+00)**2.0D+00 , dmin1(D,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(E,0.0D+00)**2.0D+00 , dmin1(F,0.0D+00)**2.0D+00 ) ) &
                        - 1.0D+00

                   DC(i,j,k) = DC(i,j,k) * sign_mx(i,j,k)
                   
                elseif ( phi0(i,j,k) .LT. 0.0D+00 ) then
                   
                   DC(i,j,k) = dsqrt( dmax1( dmax1(B,0.0D+00)**2.0D+00 , dmin1(A,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(D,0.0D+00)**2.0D+00 , dmin1(C,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(F,0.0D+00)**2.0D+00 , dmin1(E,0.0D+00)**2.0D+00 ) ) &
                        - 1.0D+00
                   
                   DC(i,j,k) = DC(i,j,k) * sign_mx(i,j,k)
                   
                endif

             endif

          enddo
       enddo
    enddo

    phi_n1 = phi_n - dt_*DC

    if ( sx == 1 ) then
       phi_n1(0,:,:) = 3.0D+00*phi_n1(1,:,:) &
            - 3.0D+00/2.0D+00*phi_n1(2,:,:) &
            + phi_n1(3,:,:)
    endif
    
    if ( ex == ex_max ) then
       phi_n1(ex_max+1,:,:) = 3.0D+00*phi_n1(ex_max,:,:) &
            - 3.0D+00/2.0D+00*phi_n1(ex_max-1,:,:) &
            + phi_n1(ex_max-2,:,:)
    endif

    phi_n1(:,:,0) = 3.0D+00*phi_n1(:,:,1) &
         - 3.0D+00/2.0D+00*phi_n1(:,:,2) &
         + phi_n1(:,:,3)

    phi_n1(:,:,maxn) = 3.0D+00*phi_n1(:,:,maxn-1) &
         - 3.0D+00/2.0D+00*phi_n1(:,:,maxn-2) &
         + phi_n1(:,:,maxn-3)
!!$
!!$    call boundary_conditions_levelset(phi_n1)

    call exchange2d(phi_n1,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    call WENO5_derivatives(phi_n1,A_,B_,C_,D_,E_,F_)
    
    do k = 1,n-1
       do j = sy,ey
          do i = sx,ex

             if ( istatus(i,j,k) == 1 ) then

                DC(i,j,k) = sign_mx(i,j,k) * dabs( phi_n(i,j,k) ) - D_mx(i,j,k)
                DC(i,j,k) = DC(i,j,k)/dx


             else

                !----- Along x-direction -----!
                
                A = A_(i,j,k)
                B = B_(i,j,k)
                
                !----- Along y-direction -----!
                
                C = C_(i,j,k)
                D = D_(i,j,k)
                                
                !----- Along z-direction -----!
                
                E = E_(i,j,k)
                F = F_(i,j,k)


                if ( phi0(i,j,k) .GE. 0.0D+00 ) then
                   
                   DC(i,j,k) = dsqrt( dmax1( dmax1(A,0.0D+00)**2.0D+00 , dmin1(B,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(C,0.0D+00)**2.0D+00 , dmin1(D,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(E,0.0D+00)**2.0D+00 , dmin1(F,0.0D+00)**2.0D+00 ) ) &
                        - 1.0D+00

                   DC(i,j,k) = DC(i,j,k) * sign_mx(i,j,k)
                   
                elseif ( phi0(i,j,k) .LT. 0.0D+00 ) then
                   
                   DC(i,j,k) = dsqrt( dmax1( dmax1(B,0.0D+00)**2.0D+00 , dmin1(A,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(D,0.0D+00)**2.0D+00 , dmin1(C,0.0D+00)**2.0D+00 ) &
                        + dmax1( dmax1(F,0.0D+00)**2.0D+00 , dmin1(E,0.0D+00)**2.0D+00 ) ) &
                        - 1.0D+00
                   
                   DC(i,j,k) = DC(i,j,k) * sign_mx(i,j,k)
                   
                endif
!!$                
             endif

          enddo
       enddo
    enddo

    phi_n2 = phi_n1 - dt_*DC

    if ( sx == 1 ) then
       phi_n2(0,:,:) = 3.0D+00*phi_n2(1,:,:) &
            - 3.0D+00/2.0D+00*phi_n2(2,:,:) &
            + phi_n2(3,:,:)
    endif
    
    if ( ex == ex_max ) then
       phi_n2(ex_max+1,:,:) = 3.0D+00*phi_n2(ex_max,:,:) &
            - 3.0D+00/2.0D+00*phi_n2(ex_max-1,:,:) &
            + phi_n2(ex_max-2,:,:)
    endif

    phi_n2(:,:,0) = 3.0D+00*phi_n2(:,:,1) &
         - 3.0D+00/2.0D+00*phi_n2(:,:,2) &
         + phi_n2(:,:,3)

    phi_n2(:,:,maxn) = 3.0D+00*phi_n2(:,:,maxn-1) &
         - 3.0D+00/2.0D+00*phi_n2(:,:,maxn-2) &
         + phi_n2(:,:,maxn-3)
!!$
!!$    call boundary_conditions_levelset(phi_n2)

    call exchange2d(phi_n2,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    phi_reinit_updated = ( phi_n + phi_n2)/2.0D+00
    call boundary_conditions_levelset(phi_reinit_updated)
    call exchange2d(phi_reinit_updated,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    error_val = 0.0D+00
    
    do k = 0,n
       do j = sy,ey
          do i = sx,ex
             dirac = Dirac_function(phi_reinit_updated(i,j,k))
             error_val = error_val + dirac * ( phi_reinit(i,j,k)-phi_reinit_updated(i,j,k) )**2.0D+00
          enddo
       enddo
    enddo
    
    phi_reinit = phi_reinit_updated
    
    call exchange2d(phi_reinit,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    call mpi_allreduce(error_val, error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm2d_quasiperiodic, ierr)
    error = error*dx*dy*dz/(Lx*Ly*Lz)
    error = dsqrt(error)/dt_

    return
  end subroutine iteration_levelset_WENO5


  subroutine WENO5_derivatives(phi_,A,B,C,D,E,F)

    !---
    !
    !   This subroutine computed the upwinding convective term for the levelset function.
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

    do k = 1,maxn-1
       do j = sy,ey
          do i = sx,ex

             !----- Computation of the fluxes in the x-direction -----!

             if ( ( i .GE. 3 ) .AND. ( i .LE. l-3 ) ) then

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

             if ( ( k .GE. 3 ) .AND. ( k .LE. maxn-3 ) ) then

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

          enddo
       enddo
    enddo

    A(:,:,0) = 2.*A(:,:,1) - A(:,:,2)
    B(:,:,0) = 2.*B(:,:,1) - B(:,:,2)
    C(:,:,0) = 2.*C(:,:,1) - C(:,:,2)
    D(:,:,0) = 2.*D(:,:,1) - D(:,:,2)
    E(:,:,0) = 2.*E(:,:,1) - E(:,:,2)
    F(:,:,0) = 2.*F(:,:,1) - F(:,:,2)

    A(:,:,maxn) = 2.*A(:,:,maxn-1) - A(:,:,maxn-2)
    B(:,:,maxn) = 2.*B(:,:,maxn-1) - B(:,:,maxn-2)
    C(:,:,maxn) = 2.*C(:,:,maxn-1) - C(:,:,maxn-2)
    D(:,:,maxn) = 2.*D(:,:,maxn-1) - D(:,:,maxn-2)
    E(:,:,maxn) = 2.*E(:,:,maxn-1) - E(:,:,maxn-2)
    F(:,:,maxn) = 2.*F(:,:,maxn-1) - F(:,:,maxn-2)
    

    ! Non-augmented phi fluxes: exchange first-order halos only.  
    
    call exchange2d(A,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(C,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(E,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(B,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(D,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(F,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    return
  end subroutine WENO5_derivatives




end module tpls_levelset