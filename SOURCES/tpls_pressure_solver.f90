module tpls_pressure_solver

  use tpls_levelset
  implicit none

contains





!-----------------------------------------------------------------------------------------------------




  subroutine Poisson_solver(phi,vx,vy,vz,pressure)

    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(in)    :: phi
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(inout) :: vx
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(inout) :: vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(inout) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(inout) :: pressure

    !----- Intermediate arrays -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p) :: density
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p) :: RHS_p
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p) :: pressure_temp, delta_pressure

    !----- Miscellaneous -----!

    double precision :: residual, erreur, norme
    integer          :: iteration, sub_iteration
    integer          :: i, j, k, n

    n = maxn
    t_temp = mpi_wtime()

    !----- Compute the density array -----!

    call get_density(density,phi)

    !----- Compute the RHS of the Poisson equation -----!

    call get_source_pres(RHS_p,vx,vy,vz)

    !----- Begining of the Poisson solver -----!

    residual = 1.0D+00
    delta_pressure = pressure

    Outer_iteration : do iteration = 1, max_iteration_poisson, max_srj
       Inner_iteration : do sub_iteration = 1, max_srj

          pressure_temp = delta_pressure
          
          !----- Scheduled-Relaxation Jacobi -----!

          call Scheduled_relaxation_jacobi(delta_pressure,RHS_p,density,relax_srj(sub_iteration))
          
          !----- Check for residual -----!

          call get_difference(delta_pressure,pressure_temp,ex,ey,ez_p,sx,sy,sz_p,erreur)
          call mpi_allreduce(erreur,residual,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d_quasiperiodic,ierr)
          residual = residual*dx*dy*dz
!!$          residual = residual / (Lx*Ly*Lz)
          residual = dsqrt(residual)
          
          if ( residual.LT.tolerance_poisson ) EXIT Inner_iteration
          if ( isnan(residual) )               STOP 'The Poisson solver blew up!'
          
       enddo Inner_iteration

       if ( residual.LT.tolerance_poisson ) EXIT Outer_iteration
       
    enddo Outer_iteration
    
    pressure = delta_pressure

    call exchange2d(pressure,stride_p_xz,stride_p_yz,neighbours, &
         ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    call get_uvw_pres(vx,vy,vz,pressure,density)
    
    
    if ( my_id == master_id ) then
       write(*,*)'Iteration : ', istep, 'p---residual is ', residual, '(', iteration + sub_iteration, ')'
       write(*,*)
    end if


    time_pressure = time_pressure + (mpi_wtime() - t_temp)
    
  end subroutine Poisson_solver





!-----------------------------------------------------------------------------------------------------





  subroutine get_source_pres(RHS,vx,vy,vz)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none
    
    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(in) :: vx
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(in) :: vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-1), intent(in) :: vz

    !----- Output -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(out) :: RHS

    !----- Intermediate arrays -----!

    integer :: i,j,k
    double precision :: scalar_plusx,scalar_minusx

    RHS = 0.0D+00
    
!!$    do k = 1,maxn-1
!!$       do j = sy,ey
!!$          do i = sx,ex
!!$
!!$             rhs(i,j,k) = ((vx(i,j-1,k-1)-vx(i-1,j-1,k-1))/dx) &
!!$                  + ((vy(i,j,k-1)-vy(i,j-1,k-1))/dy)           &
!!$                  + ((vz(i,j-1,k)-vz(i,j-1,k-1))/dz)
!!$
!!$
!!$          end do
!!$       end do
!!$    end do

!!$    do k = 1,maxn-1
!!$       do j = sy,ey
!!$          do i = sx,ex
!!$
!!$             rhs(i,j,k) = ((vx(i+1,j,k-1) - vx(i,j,k-1))/dx) &
!!$                  + ((vy(i,j+1,k-1) - vy(i,j,k-1))/dy) &
!!$                  + ((vz(i,j,k)     - vz(i,j,k-1))/dz)
!!$                  
!!$          enddo
!!$       enddo
!!$    enddo
!!$             
!!$
    RHS(sx:ex,sy:ey,1:maxn-1) =                                           &
           ( vx(sx+1:ex+1,sy:ey,0:maxn-2) - vx(sx:ex,sy:ey,0:maxn-2) )/dx &
         + ( vy(sx:ex,sy+1:ey+1,0:maxn-2) - vy(sx:ex,sy:ey,0:maxn-2) )/dy &
         + ( vz(sx:ex,sy:ey,1:maxn-1)     - vz(sx:ex,sy:ey,0:maxn-2) )/dz
    
    RHS = RHS/dt

    call exchange2d(RHS,stride_p_xz,stride_p_yz,neighbours, &
         ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    return
  end subroutine get_source_pres





!-----------------------------------------------------------------------------------------------------





  subroutine Scheduled_relaxation_jacobi(pres,RHS,density,relaxation)
    
    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none
    
    !----- Inputs and outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: pres
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: RHS
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: density

    !----- Intermediate arrays -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: residual
    double precision :: diag_val, ax, ay, az
    double precision :: relaxation    
    double precision :: rho_plusx,rho_minusx,rho_plusy,rho_minusy,rho_plusz,rho_minusz

    integer :: i,j,k,n
    
    n = maxn

    ax = 1.d0/(dx*dx)
    ay = 1.d0/(dy*dy)
    az = 1.d0/(dz*dz)

    residual = 0.d0

    if(.not.density_matched) then
       
       do k=1,n-1
          do j=sy,ey
             do i=sx,ex
                
                rho_plusx  = (density(i+1,j,k)+density(i,j,k))/2.d0
                rho_minusx = (density(i-1,j,k)+density(i,j,k))/2.d0
                
                rho_plusy  = (density(i,j+1,k)+density(i,j,k))/2.d0
                rho_minusy = (density(i,j-1,k)+density(i,j,k))/2.d0
                
                rho_plusz  = (density(i,j,k+1)+density(i,j,k))/2.d0
                rho_minusz = (density(i,j,k-1)+density(i,j,k))/2.d0
                
                diag_val = ax*( (1.d0/rho_plusx)+(1.d0/rho_minusx)) + &
                     ay*( (1.d0/rho_plusy)+(1.d0/rho_minusy)) + &
                     az*( (1.d0/rho_plusz)+(1.d0/rho_minusz))
                
                residual(i,j,k) = (1.d0/diag_val) * ( &
                     ax*( (1.d0/rho_minusx)*pres(i-1,j,k) + (1.d0/rho_plusx)*pres(i+1,j,k))  + &
                     ay*( (1.d0/rho_minusy)*pres(i,j-1,k) + (1.d0/rho_plusy)*pres(i,j+1,k))  + &
                     az*( (1.d0/rho_minusz)*pres(i,j,k-1) + (1.d0/rho_plusz)*pres(i,j,k+1))  - RHS(i,j,k) )
                
             end do
          end do
       end do

    elseif(density_matched) then

!!$       do k = 1,n-1
!!$          do j = sy,ey
!!$             do i = sx,ex
!!$
!!$                residual(i,j,k) = (1.d0/6.d0) * (    &
!!$                       pres(i+1,j,k) + pres(i-1,j,k) &
!!$                     + pres(i,j+1,k) + pres(i,j-1,k) &
!!$                     + pres(i,j,k+1) + pres(i,j,k-1) &
!!$                     - dz*dz*rhs(i,j,k)              )
!!$
!!$             enddo
!!$          enddo
!!$       enddo
!!$
       residual(sx:ex,sy:ey,1:n-1) = (1.0D+00/6.0D+00) * ( &
              pres(sx+1:ex+1,sy:ey,1:n-1) + pres(sx-1:ex-1,sy:ey,1:n-1) &
            + pres(sx:ex,sy+1:ey+1,1:n-1) + pres(sx:ex,sy-1:ey-1,1:n-1) &
            + pres(sx:ex,sy:ey,2:n)       + pres(sx:ex,sy:ey,0:n-2)     &
            - dz*dz*rhs(sx:ex,sy:ey,1:n-1)                             )
       
    endif
       
    pres = (1.0D+00-relaxation)*pres + relaxation*residual
    
    !----- Neumann boundary condition dp/dx = 0 at the inflow -----!
    
    if( sx == 1 ) then
       pres(0,:,:) = pres(1,:,:)
    endif

    !----- Neumann boundary condition dp/dz = 0 on the upper and lower walls -----!
    
    pres(:,:,0)    = pres(:,:,1)
    pres(:,:,maxn) = pres(:,:,maxn-1)
    
    !----- Dirichlet boundary condition p = 0 at the outflow -----!
    
    if ( ex == ex_max ) then
       pres(ex_max+1,:,:) = -pres(ex_max,:,:)
    end if
    
    call exchange2d(pres,stride_p_xz,stride_p_yz,neighbours, &
         ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    return
  end subroutine Scheduled_relaxation_jacobi





!-----------------------------------------------------------------------------------------------------




  
  subroutine get_uvw_pres(u3,v3,w3,p,density)
    
    use tpls_constants
    use tpls_mpi
    use tpls_userchk
    use mpi
    implicit none

    !----- Inputs / Ouputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(inout) :: u3
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(inout) :: v3
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-1), intent(inout) :: w3
    
    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)      :: density
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in)      :: p

    !----- Miscellaneous -----!
    
    integer :: i,j,k
    integer :: sx_loc,ex_loc
    double precision :: rho_ugrid,rho_vgrid,rho_wgrid
    double precision, dimension(0:maxn-2) :: u_inlet
    double precision, dimension(0:maxn) :: phi_inlet

    call user_inflow(u_inlet,phi_inlet)
    
    do k = 0,maxn-2
       do j = sy,ey
          do i = sx,ex
             rho_ugrid = ( density(i+1,j+1,k+1) + density(i,j+1,k+1) ) / 2.0D+00
             u3(i,j,k) = u3(i,j,k) - (dt/dx)*( 1.d0/rho_ugrid ) * ( p(i,j,k+1) - p(i-1,j,k+1) )
          end do
       end do
    end do

    if ( sx == 1 ) then
       do j = sy,ey
          u3(1,j,:) = u_inlet(:)
       enddo
       u3(0,:,:) = 2.0D+00*u3(1,:,:) - u3(2,:,:)
    endif

    !----- No-slip boundary conditions on the walls : u = 0 -----!
    
    u3(:,:,sz_uv) = (1.0D+00/3.0D+00)*u3(:,:,sz_uv+1)
    u3(:,:,ez_uv) = (1.0D+00/3.0D+00)*u3(:,:,ez_uv-1)
    
    !----- Outflow boundary condition : du/dx = 0 -----!
    
    if ( ex == ex_max ) then       
       u3(ex_max+1,:,:) = u3(ex_max,:,:)
    end if
    
    do k = 0,maxn-2
       do j = sy,ey
          do i = sx,ex
             rho_vgrid = ( density(i,j+1,k+1) + density(i,j,k+1) ) / 2.d0
             v3(i,j,k) = v3(i,j,k) - (dt/dy)*( 1.d0/rho_vgrid ) * ( p(i,j,k+1) - p(i,j-1,k+1) )             
          end do
       end do
    end do

    if ( sx == 1 ) then
       v3(0,:,:) = -v3(1,:,:)
    endif
    
    !----- No-slip boundary condition on the walls : v = 0 -----!
    
    v3(:,:,sz_uv) = (1.0D+00/3.0D+00)*v3(:,:,sz_uv+1)
    v3(:,:,ez_uv) = (1.0D+00/3.0D+00)*v3(:,:,ez_uv-1)
    
    !----- Outflow boundary condition : dv/dx = 0 -----!
    
    if ( ex == ex_max ) then
       v3(ex_max+1,:,:) = v3(ex_max,:,:)
    end if
    
    do k = 0,maxn-1
       do j = sy,ey
          do i = sx,ex
             rho_wgrid = ( density(i,j+1,k+1) + density(i,j+1,k) ) / 2.0D+00
             w3(i,j,k) = w3(i,j,k) - (dt/dz)*( 1.d0/rho_wgrid ) * ( p(i,j,k+1) - p(i,j,k) )
          end do
       end do
    end do

    if ( sx==1 ) then
       w3(0,:,:) = -w3(1,:,:)
    endif

    !----- No-slip boundary condition on the walls : w = 0 -----!
    
    w3(:,:,0)      = 0.d0
    w3(:,:,maxn-1) = 0.d0
    
    !----- Outflow boundary condition : dw/dx = 0 -----!
    
    if ( ex == ex_max ) then
       w3(ex_max+1,:,:) = w3(ex_max,:,:)
    end if

    call exchange2d(u3,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(v3,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(w3,stride_w_xz ,stride_w_yz ,neighbours,ex,ey,ez_w ,sx,sy,sz_w ,comm2d_quasiperiodic)
    
    return
  end subroutine get_uvw_pres





!-----------------------------------------------------------------------------------------------------





end module tpls_pressure_solver
