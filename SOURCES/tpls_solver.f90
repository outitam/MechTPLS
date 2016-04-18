module tpls_solver

  implicit none

contains

  subroutine tpls_dns(vx,vy,vz,pres,phi)

    use tpls_constants
    use tpls_configuration
    use tpls_mpi
    use tpls_selective_frequency_damping
    use tpls_levelset
    use tpls_fluid
    use tpls_pressure_solver
    use tpls_io
    use tpls_userchk
    use tpls_maths
    use mpi
    implicit none

    double precision, dimension(sx-1:ex+1,sy-1:sy+1,sz_uv:ez_uv), intent(inout) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(inout) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(inout) :: pres, phi
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p) :: curvature, interface_loc

    double precision :: CFL, volume
    
    !----- Allocation -----!
    
    
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
    
    
    !----- Arrays for the Selective Frequency Damping -----!
    
    if ( if_sfd ) then
       
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

    !----- Initialization -----!

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
    
    Time_Loop : do istep = 1,nsteps
       
       !----- Computing the CFL number -----!
       
       cfl = get_cfl(vx,vy,vz)
       volume = get_volume(phi)
       
       if ( my_id == master_id ) then
          write(*,*) 'Iteration : ', istep, 'CFL : ', CFL, 'Time : ', dt*istep
          write(*,*) 'Iteration : ', istep, 'Volume phase 1 : ', volume
       endif

       volume = get_volume(-phi)
       if ( my_id == master_id ) then
          write(*,*) 'Iteration : ', istep, 'Volume phase 2 : ', volume
       endif

       
       !----- Advance the levelset equation -----!
       
       call advance_levelset(phi,vx,vy,vz)
       
       !----- Advance the momentum equation -----!
       
       call advance_fluid(phi,vx,vy,vz,pres)     
       
       !----- Impose the divergence-free constraint -----!
       
       call Poisson_solver(phi,vx,vy,vz,pres)
       
       !----- Performs any required runtime post-processing -----!
       
       call user_postprocess(phi,pres,vx,vy,vz)
       
       call exchange2d(phi,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
       call exchange2d(vx,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(vy,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(vz, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,comm2d_quasiperiodic)
       call exchange2d(pres,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

!!$       call get_curvature_original(curvature,phi)
       call get_curvature(curvature,phi)
       
       !----- Periodic output to files -----!
       
       if( mod(istep,iostep).eq.0 ) then
          
          call outpost_exact(vx,vy,vz,pres,phi,'   ')
          if (if_sfd .and. if_exact_restart_save) then
              call backup_exactrest_sfd(vx_bar,vy_bar,vz_bar,fx_sfd,fy_sfd,fz_sfd)
          endif
       end if
       
    end do Time_Loop

    !----- Deallocation -----!

    deallocate(conv0_u,conv1_u)
    deallocate(csf_u0,csf_u1)
    deallocate(diffusion1_u)
    deallocate(conv0_v,conv1_v)
    deallocate(csf_v0,csf_v1)
    deallocate(diffusion1_v)
    deallocate(conv0_w,conv1_w)
    deallocate(csf_w0,csf_w1)
    deallocate(diffusion1_w)
    
    !----- Arrays for the Selective Frequency Damping -----!
    
    if ( if_sfd ) then
       
       deallocate(vx_bar,vy_bar,vz_bar)
       deallocate(fx_sfd,fy_sfd,fz_sfd)
       
    endif

    
    return
  end subroutine tpls_dns

  subroutine tpls_linear_stability(ub,vb,wb,prb,phi_b)

    use tpls_constants
    use tpls_maths
    use tpls_mpi
    use tpls_io
    use mpi
    implicit none

    !----- Inputs: Base flow arrays -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: ub, vb
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: wb
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(in) :: phi_b, prb

    !----- Perturbation array -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: vxp, vyp
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: vzp
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: phi_p, prp

    !----- Krylov basis -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv,krylov_dim+1) :: vx_, vy_
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w,krylov_dim+1)   :: vz_
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p,krylov_dim+1)   :: phi_, pr_

    !----- Upper Hessenberg matrix -----!

    double precision, dimension(krylov_dim,krylov_dim) :: H_mat
    double precision, dimension(1,krylov_dim)          :: b_vec

    !----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----!

    double complex, dimension(krylov_dim)            :: VP
    double complex, dimension(krylov_dim,krylov_dim) :: FP

    !----- Arrays for storage/output of a given eigenmode of the NS operator -----
    double complex, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: FP_Cx, FP_Cy
    double complex, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: FP_Cz
    double complex, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: FP_Cp, FP_Cphi

    !----- Miscellaneous -----!

    integer :: ifich1 = 111, ifich2 = 222, ifich3 = 333
    integer :: mstart

    double precision :: sampling_period
    double precision, dimension(krylov_dim ) :: residual
    double precision :: alpha, beta
    logical :: is_converged
    integer :: check_convergence
    double precision :: amplitude = 1.0D-05

    !----- Working arrays -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: workx, worky
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: workz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: workp, workphi

    !----- Schur decomposition -----!

    double precision, dimension(krylov_dim,krylov_dim) :: Q1
    double precision, dimension(:,:), allocatable :: Q2, Q
    integer :: have_been_selected

    !----- Miscellaneous -----!

    integer :: i, j, k, index_mode
    double precision :: tolerance
    parameter ( tolerance = 1.0D-06 )
    
    !----- Initialization -----!

    sampling_period = dt*nsteps

    vx_ = 0.0D+00
    vy_ = 0.0D+00
    vz_ = 0.0D+00

    H_mat = 0.0D+00
    b_vec = 0.0D+00

    VP = ( 0.0D+00 , 0.0D+00 )
    FP = ( 0.0D+00 , 0.0D+00 )

    FP_Cx   = ( 0.0D+00 , 0.0D+00 )
    FP_Cy   = ( 0.0D+00 , 0.0D+00 )
    FP_Cz   = ( 0.0D+00 , 0.0D+00 )
    FP_Cp   = ( 0.0D+00 , 0.0D+00 )
    FP_Cphi = ( 0.0D+00 , 0.0D+00 )
    
    residual          = 1.0D0+00
    is_converged      = .false.
    check_convergence = 0

    workx   = 0.0D+00
    worky   = 0.0D+00
    workz   = 0.0D+00
    workp   = 0.0D+00
    workphi = 0.0D+00

    Q1 = 0.0D+00

    have_been_selected = 0
    index_mode         = 0
    
    !-------------------------------------------------
    !-----                                       -----
    !----- Creation of the initial Krylov vector -----
    !-----                                       -----
    !-------------------------------------------------

    if ( my_id == master_id) then
       write(*,*)
       write(*,*) 'Creating the intial Krylov vector'
       write(*,*)
    endif
    call mpi_barrier(comm2d_quasiperiodic,ierr)

    call random_number(vxp)
    call random_number(vyp)
    call random_number(vzp)

    phi_p = 0.0D+00
    prp   = 0.0D+00

    vxp = ub
    vyp = vb
    vzp = wb

    prp = prb
    phi_p = phi_b
    
    !----- Inflow boundary condition : u = 0 -----!

    if ( my_id == master_id ) then
       write(*,*) '--> Imposing the boundary conditions.'
    endif
    call mpi_barrier(comm2d_quasiperiodic,ierr)

    if ( sx == 1 ) then       
       vxp(1,:,:) = 0.0D+00
       vxp(0,:,:) = 2.0D+00*vxp(1,:,:)-vxp(2,:,:)
    end if
    
    !----- Outflow boundary condition : du/dx = 0 -----!
    
    if ( ex == ex_max ) then
       vxp(ex_max,sy:ey,:)   = vxp(ex_max-1,sy:ey,:)
       vxp(ex_max+1,sy:ey,:) = 2.d0*vxp(ex_max,sy:ey,:) - vxp(ex_max-1,sy:ey,:)
    end if
    
    !----- No-slip boundary conditions on the walls : u = 0 -----!
    
    vxp(:,:,sz_uv) = (1.0D+00/3.0D+00)*vxp(:,:,sz_uv+1)
    vxp(:,:,ez_uv) = (1.0D+00/3.0D+00)*vxp(:,:,ez_uv-1)

    !----- Inflow boundary condition : v = 0 -----!
    
    if ( sx == 1 ) then
       vyp(0,:,:) = -vyp(1,:,:)
    end if
    
    !----- Outflow boundary condition : dv/dx = 0 -----!
    
    if ( ex == ex_max ) then
       vyp(ex_max+1,sy:ey,:) = vyp(ex_max,sy:ey,:)
    end if
    
    !----- No-slip boundary condition on the walls : v = 0 -----!
    
    vyp(:,:,sz_uv) = (1.0D+00/3.0D+00)*vyp(:,:,sz_uv+1)
    vyp(:,:,ez_uv) = (1.0D+00/3.0D+00)*vyp(:,:,ez_uv-1)
    
    !----- No-slip boundary condition on the walls : w = 0 -----!
    
    vzp(:,:,0)      = 0.0D+00
    vzp(:,:,maxn-1) = 0.0D+00
    
    !----- Inflow boundary condition : w = 0 -----!
    
    if ( sx == 1 ) then
       vzp(0,:,:) = -vzp(1,:,:)
    end if
    
    !----- Outflow boundary condition : dw/dx = 0 -----!
    
    if ( ex == ex_max ) then
       vzp(ex_max+1,sy:ey,:) = vzp(ex_max,sy:ey,:)
    end if

    !----- Normalized to unit-norm -----!

    call inner_product(alpha &
         , phi_p, phi_p      &
         , vxp, vyp, vzp     &
         , vxp, vyp, vzp)

    alpha = 1.0D+00/dsqrt(alpha)

    vxp   = alpha*vxp
    vyp   = alpha*vyp
    vzp   = alpha*vzp

    prp   = alpha*prp
    phi_p = alpha*phi_p

    call exchange2d(phi_p,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(vxp,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(vyp,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(vzp, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,comm2d_quasiperiodic)
    call exchange2d(prp,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

    workx   = ub + amplitude*vxp
    worky   = vb + amplitude*vyp
    workz   = wb + amplitude*vzp
    workp   = prb
    workphi = phi_b
    
    call exchange2d(workphi,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(workx,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(worky,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(workz, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,comm2d_quasiperiodic)
    call exchange2d(workp,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    if ( my_id == master_id ) then
       write(*,*) 'Calling TPLS_DNS'
       write(*,*)
    endif

    call tpls_dns( workx, worky, workz, workp, workphi )
    
    vxp   = (workx   - ub)    / amplitude
    vyp   = (worky   - vb)    / amplitude
    vzp   = (workz   - wb)    / amplitude
    prp   = (workp   - prb)   / amplitude
    phi_p = (workphi - phi_b) / amplitude

    !----- Inflow boundary condition : u = 0 -----!

    if ( my_id == master_id ) then
       write(*,*) '--> Imposing the boundary conditions.'
    endif
    call mpi_barrier(comm2d_quasiperiodic,ierr)
    
    if ( sx == 1 ) then       
       vxp(1,:,:) = 0.0D+00
       vxp(0,:,:) = 2.0D+00*vxp(1,:,:) - vxp(2,:,:)
    end if
    
    !----- Outflow boundary condition : du/dx = 0 -----!
    
    if ( ex == ex_max ) then
       vxp(ex_max,sy:ey,:)   = vxp(ex_max-1,sy:ey,:)
       vxp(ex_max+1,sy:ey,:) = 2.d0*vxp(ex_max,sy:ey,:) - vxp(ex_max-1,sy:ey,:)
    end if
    
    !----- No-slip boundary conditions on the walls : u = 0 -----!
    
    vxp(:,:,sz_uv) = (1.0D+00/3.0D+00)*vxp(:,:,sz_uv+1)
    vxp(:,:,ez_uv) = (1.0D+00/3.0D+00)*vxp(:,:,ez_uv-1)

    !----- Inflow boundary condition : v = 0 -----!
    
    if ( sx == 1 ) then
       vyp(0,:,:) = -vyp(1,:,:)
    end if
    
    !----- Outflow boundary condition : dv/dx = 0 -----!
    
    if ( ex == ex_max ) then
       vyp(ex_max+1,sy:ey,:) = vyp(ex_max,sy:ey,:)
    end if
    
    !----- No-slip boundary condition on the walls : v = 0 -----!
    
    vyp(:,:,sz_uv) = (1.0D+00/3.0D+00)*vyp(:,:,sz_uv+1)
    vyp(:,:,ez_uv) = (1.0D+00/3.0D+00)*vyp(:,:,ez_uv-1)
    
    !----- No-slip boundary condition on the walls : w = 0 -----!
    
    vzp(:,:,0)      = 0.d0
    vzp(:,:,maxn-1) = 0.d0
    
    !----- Inflow boundary condition : w = 0 -----!
    
    if ( sx == 1 ) then
       vzp(0,:,:) = -vzp(1,:,:)
    end if
    
    !----- Outflow boundary condition : dw/dx = 0 -----!
    
    if ( ex == ex_max ) then
       vzp(ex_max+1,sy:ey,:) = vzp(ex_max,sy:ey,:)
    end if

    call inner_product(alpha,phi_p,phi_p,vxp,vyp,vzp,vxp,vyp,vzp)
    
    alpha = 1.0D+00/dsqrt(alpha)
    
    vxp   = alpha*vxp
    vyp   = alpha*vyp
    vzp   = alpha*vzp

    prp   = alpha*prp
    phi_p = alpha*phi_p

    call exchange2d(phi_p,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    call exchange2d(vxp,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(vyp,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
    call exchange2d(vzp, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,comm2d_quasiperiodic)
    call exchange2d(prp,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    !----- Storing the starting vector in the Krylov basis -----!

    vx_(:,:,:,1)  = vxp
    vy_(:,:,:,1)  = vyp
    vz_(:,:,:,1)  = vzp

    pr_(:,:,:,1)  = prp
    phi_(:,:,:,1) = phi_p

    mstart = 1

    if ( my_id == master_id ) then
       write(*,*)
       write(*,*) '-----------------------------------------------------'
       write(*,*) '-----     Begining of the Arnoldi algorithm     -----'
       write(*,*) '-----------------------------------------------------'
       write(*,*)
    endif

    call mpi_barrier(comm2d_quasiperiodic,ierr)

    Krylov_Schur : do while (.not. is_converged )

       !----- m steps Arnoldi factorization -----!

       call Arnoldi_factorization(vx_, vy_, vz_, pr_, phi_, &
            ub, vb, wb, prb, phi_b, &
            H_mat, beta, mstart)

       !----- Check for convergence -----!

       VP = ( 0.0D+00 , 0.0D+00 )
       FP = ( 0.0D+00 , 0.0D+00 )

       call eigendecomposition(H_mat,krylov_dim,VP,FP)

       b_vec               = 0.0D+00
       b_vec(1,krylov_dim) = beta

       check_convergence = 0
       residual = 1.0D+00
       
       do i = 1, krylov_dim
          
          residual(i) = abs(beta*FP(krylov_dim,i))
          if ( my_id == master_id ) write(*,*) 'Residual of eigenvalue', i, ' :', residual(i)
          
          if ( residual(i) .LT. tolerance ) then
             check_convergence = check_convergence + 1
          endif

          if ( my_id == master_id ) write(*,*) check_convergence, ' eigenvalues have converged'
       enddo

       if ( check_convergence.GE.8 ) then

          is_converged = .true.

       else

          VP = ( 0.0D+00 , 0.0D+00 )
          FP = ( 0.0D+00 , 0.0D+00 )

          call Schur_decomposition(H_mat,Q1,VP,krylov_dim,mstart)

          H_mat(1:mstart,mstart+1:krylov_dim)            = 0.0D+00
          H_mat(mstart+1:krylov_dim,1:mstart)            = 0.0D+00
          H_mat(mstart+1:krylov_dim,mstart+1:krylov_dim) = 0.0D+00

          if ( mstart == 0 ) then

             mstart = 1

             H_mat = 0.0D+00
             
             vx_(:,:,:,mstart)  =  vx_(:,:,:,krylov_dim+1)
             vy_(:,:,:,mstart)  =  vy_(:,:,:,krylov_dim+1)
             vz_(:,:,:,mstart)  =  vz_(:,:,:,krylov_dim+1)
             pr_(:,:,:,mstart)  =  pr_(:,:,:,krylov_dim+1)
             phi_(:,:,:,mstart) = phi_(:,:,:,krylov_dim+1)

          else

             if ( my_id == master_id ) then
                write(*,*)
                write(*,*) mstart, 'Ritz eigenpairs have been selected.'
                write(*,*) check_convergence, 'eigenpairs have converged.'
                write(*,*)
             endif

             allocate(Q2(mstart,mstart))
             allocate(Q(krylov_dim,mstart))
             Q2 = 0.0D+00
             Q  = 0.0D+00

             if ( my_id == master_id ) then
                write(*,*) 'Performing Hessenberg decomposition.'
                write(*,*)
             endif

             call Hessenberg_decomposition(H_mat(1:mstart,1:mstart),Q2,mstart)

             call mpi_barrier(comm2d_quasiperiodic, ierr)

             if ( my_id == master_id ) then
                write(*,*) 'Computation of the transformation matrix Q.'
                write(*,*)
             endif

             Q = Q1(:,1:mstart)

             !----- Re-order the Krylov basis  -----!

             if ( my_id == master_id ) then
                write(*,*) 'Re-ordering the Krylov basis :'
                write(*,*)
             endif

             call mpi_barrier(comm2d_quasiperiodic, ierr)

             call mxm_tpls(vx_(:,:,:,1:krylov_dim), &
                  Q, &
                  krylov_dim, mstart, &
                  sx, ex, sy, ey, sz_uv, ez_uv)

             if ( my_id == master_id) write(*,*) '---> vx : done.'
             
             call mxm_tpls(vy_(:,:,:,1:krylov_dim), &
                  Q, &
                  krylov_dim, mstart, &
                  sx, ex, sy, ey, sz_uv, ez_uv)
             
             if ( my_id == master_id) write(*,*) '---> vy : done.'

             call mxm_tpls(vz_(:,:,:,1:krylov_dim), &
                  Q, &
                  krylov_dim, mstart, &
                  sx, ex, sy, ey, sz_w, ez_w)
             
             if ( my_id == master_id) write(*,*) '---> vz : done.'

             call mxm_tpls(pr_(:,:,:,1:krylov_dim), &
                  Q, &
                  krylov_dim, mstart, &
                  sx, ex, sy, ey, sz_p, ez_p)
             
             if ( my_id == master_id) write(*,*) '---> pr : done.'

             call mxm_tpls(phi_(:,:,:,1:krylov_dim), &
                  Q, &
                  krylov_dim, mstart, &
                  sx, ex, sy, ey, sz_p, ez_p)
             
             if ( my_id == master_id) write(*,*) '---> phi : done.'

             !----- Update the matrix with b'*Q -----

             call matrix_matrix(b_vec,Q1,1,krylov_dim,krylov_dim)

             H_mat(mstart+1,1:mstart) = b_vec(1,1:mstart)
             mstart = mstart + 1

             vx_(:,:,:,mstart)  =  vx_(:,:,:,krylov_dim+1)
             vy_(:,:,:,mstart)  =  vy_(:,:,:,krylov_dim+1)
             vz_(:,:,:,mstart)  =  vz_(:,:,:,krylov_dim+1)
             pr_(:,:,:,mstart)  =  pr_(:,:,:,krylov_dim+1)
             phi_(:,:,:,mstart) = phi_(:,:,:,krylov_dim+1)
             
             deallocate(Q)             

          endif

       endif

    enddo Krylov_Schur

    !----- Sorting the eigenvalues and eigenvectors -----!
    
    call sort_eigendecomp(VP,FP,residual,krylov_dim)

    if ( my_id == master_id ) then
       open(unit = ifich1, file = 'Hessenberg.dat')
       do j = 1,krylov_dim
          do i = 1,krylov_dim
             write(ifich1,*) H_mat(i,j)
          enddo
       enddo
       close(ifich1)
    endif
        
    !----- Output all the spectrums and converged eigenmodes -----!

    if ( my_id == master_id ) then

       open(unit = ifich1,          &
            file = 'Spectre_H.dat', &
            form = 'formatted' )

       open(unit = ifich2,           &
            file = 'Spectre_NS.dat', &
            form = 'formatted' )

       open(unit = ifich3,                &
            file = 'Spectre_NS_conv.dat', &
            form = 'formatted' )

       open(unit = ifich1*ifich2,        &
            file = 'Orthonormality.dat', &
            form = 'formatted' )

    endif

    call mpi_barrier(comm2d_quasiperiodic,ierr)

    do i = 1, krylov_dim

       alpha = 0.0D+00
       
       call inner_product(alpha, &
            phi_(:,:,:,i), phi_(:,:,:,i), &
            vx_(:,:,:,i) , vy_(:,:,:,i) , vz_(:,:,:,i) , &
            vx_(:,:,:,i) , vy_(:,:,:,i) , vz_(:,:,:,i) )

       alpha = dsqrt(alpha)

       if ( my_id == master_id ) then
          write(ifich1*ifich2,*) '----- Norm of mode ', i, ' = ', alpha
          write(*,*)
       endif
       
       do j = i+1, krylov_dim

          beta  = 0.0D+00

          call inner_product(beta, &
               phi_(:,:,:,i), phi_(:,:,:,j), &
               vx_(:,:,:,i) , vy_(:,:,:,i) , vz_(:,:,:,i) , &
               vx_(:,:,:,j) , vy_(:,:,:,j) , vz_(:,:,:,j) )

          if ( my_id == master_id ) then
             write(ifich1*ifich2,*) 'Inner product (', i, ',', j,') = ', beta
          endif

       enddo
    enddo

    if ( my_id == master_id ) close(ifich1*ifich2)


    do i = 1, krylov_dim

       if ( my_id == master_id ) then
          
          write(ifich1,*) dble(VP(i)), &
               dimag(VP(i)), &
               residual(i)
          
          write(ifich2,*) log(abs(VP(i)))/sampling_period, &
               ATAN2(dimag(VP(i)),dble(VP(i)))/sampling_period, &
               residual(i)
          
       endif
       call mpi_barrier(comm2d_quasiperiodic,ierr)

       if( residual(i).LT.tolerance ) then

          index_mode = index_mode + 1
          istep = index_mode
          
          if ( my_id == master_id ) then
             
             write(ifich3,*) log(abs(VP(i)))/sampling_period, &
                  ATAN2(aimag(VP(i)),real(VP(i)))/sampling_period
             
          endif

          !----- Computation of the corresponding eigenmode -----!

          FP_Cx   = (0.0D+00,0.0D+00)
          FP_Cy   = (0.0D+00,0.0D+00)
          FP_Cz   = (0.0D+00,0.0D+00)
          FP_Cp   = (0.0D+00,0.0D+00)
          FP_Cphi = (0.0D+00,0.0D+00)

          call matvec(FP_Cx,vx_,FP(:,i),krylov_dim,sx,ex,sy,ey,sz_uv,ez_uv)
          call matvec(FP_Cy,vy_,FP(:,i),krylov_dim,sx,ex,sy,ey,sz_uv,ez_uv)
          call matvec(FP_Cz,vz_,FP(:,i),krylov_dim,sx,ex,sy,ey,sz_w,ez_w)
          call matvec(FP_Cp,pr_,FP(:,i),krylov_dim,sx,ex,sy,ey,sz_p,ez_p)
          call matvec(FP_Cphi,phi_,FP(:,i),krylov_dim,sx,ex,sy,ey,sz_p,ez_p)

          !----- Normalization to unit-norm -----!

          alpha = 0.0D+00
          beta  = 0.0D+00

          call inner_product(alpha,dble(FP_Cphi),dble(FP_Cphi) &
               , dble(FP_Cx), dble(FP_Cy), dble(FP_Cz) &
               , dble(FP_Cx), dble(FP_Cy), dble(Fp_Cz))

          call inner_product(beta,dimag(FP_Cphi),dimag(FP_Cphi) &
               , dimag(FP_Cx), dimag(FP_Cy), dimag(Fp_Cz) &
               , dimag(FP_Cx), dimag(Fp_Cy), dimag(Fp_Cz))

          if ( my_id == master_id ) then
             write(*,*) 'Norm of mode', i, 'before normalization : ', alpha, beta, dsqrt(alpha+beta)
          endif

          alpha = alpha + beta

          alpha = 1.0D+00/dsqrt(alpha)

          FP_Cx   = alpha*FP_Cx
          FP_Cy   = alpha*FP_Cy
          FP_Cz   = alpha*FP_Cz
          FP_Cp   = alpha*FP_Cp
          FP_Cphi = alpha*FP_Cphi

          alpha = 0.0D+00
          beta  = 0.0D+00

          call inner_product(alpha,dble(FP_Cphi),dble(FP_Cphi) &
               , dble(FP_Cx), dble(FP_Cy), dble(FP_Cz)         &
               , dble(FP_Cx), dble(FP_Cy), dble(Fp_Cz))
          
          call inner_product(beta,dimag(FP_Cphi),dimag(FP_Cphi) &
               , dimag(FP_Cx), dimag(FP_Cy), dimag(Fp_Cz)       &
               , dimag(FP_Cx), dimag(FP_Cy), dimag(Fp_Cz))
          
          if ( my_id == master_id ) then
             write(*,*) 'Norm of mode', i, 'after normalization : ', alpha, beta, dsqrt(alpha+beta)
          endif
          
          !----- Output the real part -----!

          vxp   = dble(FP_Cx)
          vyp   = dble(FP_Cy)
          vzp   = dble(FP_Cz)
          prp   = dble(FP_Cp)
          phi_p = dble(FP_Cphi)

          call outpost_exact(vxp,vyp,vzp,prp,phi_p,'Re_')

          !----- Output the imaginary part -----!

          vxp   = dimag(FP_Cx)
          vyp   = dimag(FP_Cy)
          vzp   = dimag(Fp_Cz)
          prp   = dimag(FP_Cp)
          phi_p = dimag(FP_Cphi)

          call outpost_exact(vxp,vyp,vzp,prp,phi_p,'Im_')

       endif

       call mpi_barrier(comm2d_quasiperiodic,ierr)

    enddo
    
    if ( my_id == master_id ) then

       close(ifich1)
       close(ifich2)
       close(ifich3)

    endif
    call mpi_barrier(comm2d_quasiperiodic,ierr)
          
    return
  end subroutine tpls_linear_stability

  subroutine Arnoldi_factorization(Qx,Qy,Qz,Qp,Qphi,ub,vb,wb,prb,phi_b,H,beta,mstart)

    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use tpls_io
    use mpi
    implicit none

    !----- Inputs: Krylov basis -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv,1:krylov_dim+1)  , intent(inout) :: Qx, Qy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w,1:krylov_dim+1)    , intent(inout) :: Qz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p,1:krylov_dim+1)    , intent(inout) :: Qp, Qphi

    !----- Inputs: Base flow -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: ub, vb
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: wb
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(in) :: prb, phi_b

    !----- Inputs: Arnoldi related -----!

    double precision, dimension(krylov_dim,krylov_dim), intent(inout) :: H
    double precision                                  , intent(out)   :: beta
    integer                                           , intent(in)    :: mstart

    !----- Miscellaneous -----!

    double precision :: alpha, norme
    integer :: mstep, i, j, k
    double precision, dimension(krylov_dim) :: h_vec
    double precision :: amplitude
    parameter ( amplitude = 1.0D-05 )

    !----- Array for the matrix-vector product w = M*v -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: wx, wy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: wz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: wp, wphi

    !----- Orthogonal residual f = w - (V,w)*V -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: fx, fy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: fz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: fp, fphi

    !----- Working arrays -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: workx, worky
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: workz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: workp, workphi

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: vxp, vyp
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: vzp
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   :: prp, phi_p

    !=================================================
    !=====                                       =====
    !=====     Initialization of most arrays     =====
    !=====                                       =====
    !=================================================

    wx   = 0.0D+00
    wy   = 0.0D+00
    wz   = 0.0D+00
    wp   = 0.0D+00
    wphi = 0.0D+00

    fx   = 0.0D+00
    fy   = 0.0D+00
    fz   = 0.0D+00
    fp   = 0.0D+00
    fphi = 0.0D+00

    h_vec = 0.0D+00

    workx   = 0.0D+00
    worky   = 0.0D+00
    workz   = 0.0D+00
    workp   = 0.0D+00
    workphi = 0.0D+00

    vxp   = 0.0D+00
    vyp   = 0.0D+00
    vzp   = 0.0D+00
    prp   = 0.0D+00
    phi_p = 0.0D+00

    Arnoldi : do mstep = mstart, krylov_dim

       if ( my_id == master_id ) then
          write(*,*)
          write(*,*) '-----'
          write(*,*) '-----'
          write(*,*) '-----     Begin Arnoldi iteration: ', mstep
          write(*,*) '-----'
          write(*,*) '-----'
          write(*,*)
       endif

       vxp   = Qx(:,:,:,mstep)
       vyp   = Qy(:,:,:,mstep)
       vzp   = Qz(:,:,:,mstep)
       prp   = Qp(:,:,:,mstep)
       phi_p = Qphi(:,:,:,mstep)

       if ( my_id == master_id ) then
          write(*,*) '--> Imposing the boundary conditions.'
       endif
       call mpi_barrier(comm2d_quasiperiodic,ierr)
       
       if ( sx == 1 ) then       
          vxp(1,:,:) = 0.0D+00
          vxp(0,:,:) = 2.0D+00*vxp(1,:,:) - vxp(2,:,:)
       end if
       
       !----- Outflow boundary condition : du/dx = 0 -----!
       
       if ( ex == ex_max ) then
          vxp(ex_max,sy:ey,:)   = vxp(ex_max-1,sy:ey,:)
          vxp(ex_max+1,sy:ey,:) = 2.d0*vxp(ex_max,sy:ey,:) - vxp(ex_max-1,sy:ey,:)
       end if
       
       !----- No-slip boundary conditions on the walls : u = 0 -----!
       
       vxp(:,:,sz_uv) = (1.0D+00/3.0D+00)*vxp(:,:,sz_uv+1)
       vxp(:,:,ez_uv) = (1.0D+00/3.0D+00)*vxp(:,:,ez_uv-1)
       
       !----- Inflow boundary condition : v = 0 -----!
       
       if ( sx == 1 ) then
          vyp(0,:,:) = -vyp(1,:,:)
       end if
       
       !----- Outflow boundary condition : dv/dx = 0 -----!
       
       if ( ex == ex_max ) then
          vyp(ex_max+1,sy:ey,:) = vyp(ex_max,sy:ey,:)
       end if
       
       !----- No-slip boundary condition on the walls : v = 0 -----!
       
       vyp(:,:,sz_uv) = (1.0D+00/3.0D+00)*vyp(:,:,sz_uv+1)
       vyp(:,:,ez_uv) = (1.0D+00/3.0D+00)*vyp(:,:,ez_uv-1)
       
       !----- No-slip boundary condition on the walls : w = 0 -----!
       
       vzp(:,:,0)      = 0.d0
       vzp(:,:,maxn-1) = 0.d0
       
       !----- Inflow boundary condition : w = 0 -----!
       
       if ( sx == 1 ) then
          vzp(0,:,:) = -vzp(1,:,:)
       end if
       
       !----- Outflow boundary condition : dw/dx = 0 -----!
       
       if ( ex == ex_max ) then
          vzp(ex_max+1,sy:ey,:) = vzp(ex_max,sy:ey,:)
       end if
       
       call exchange2d(phi_p,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
       call exchange2d(vxp,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(vyp,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(vzp, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,comm2d_quasiperiodic)
       call exchange2d(prp,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

       !----- Matrix-vector product w = M*v -----

       workx   = ub + amplitude*vxp
       worky   = vb + amplitude*vyp
       workz   = wb + amplitude*vzp
       workp   = prb
       workphi = phi_b

       call exchange2d(workphi,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
       call exchange2d(workx,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(worky,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(workz, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,comm2d_quasiperiodic)
       call exchange2d(workp,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
       
       istep = 0
       call outpost_exact(vxp,vyp,vzp,prp,phi_b,'prt')

       if ( my_id == master_id ) then
          write(*,*) 'Calling TPLS_DNS'
          write(*,*)
       endif

       call tpls_dns( workx, worky, workz, workp, workphi)

       wx   = (workx   - ub)    / amplitude
       wy   = (worky   - vb)    / amplitude
       wz   = (workz   - wb)    / amplitude
       wp   = (workp   - prb)   / amplitude
       wphi = (workphi - phi_b) / amplitude

       call exchange2d(wphi,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
       call exchange2d(wx,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(wy,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(wz, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,comm2d_quasiperiodic)
       call exchange2d(wp,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

       fx   = wx
       fy   = wy
       fz   = wz
       fp   = wp
       fphi = wphi

       !----- Compute the orthogonal residual f -----!

       h_vec = 0.0D+00

       do i = 1, mstep

          workx   = Qx(:,:,:,i)
          worky   = Qy(:,:,:,i)
          workz   = Qz(:,:,:,i)
          workp   = Qp(:,:,:,i)
          workphi = Qphi(:,:,:,i)

          call inner_product(alpha, &
               fphi, workphi, &
               fx, fy, fz, &
               workx, worky, workz)

          h_vec(i) = alpha

          fx   = fx   - h_vec(i) * workx
          fy   = fy   - h_vec(i) * worky
          fz   = fz   - h_vec(i) * workz
          fp   = fp   - h_vec(i) * workp
          fphi = fphi - h_vec(i) * workphi

       enddo

       !----- Fills in the upper Hessenberg matrix -----!
       
       H(1:mstep,mstep) = h_vec(1:mstep)
       
       !----- Compute the norm of the orthogonal residual f -----!
       
       call inner_product(beta &
            , fphi, fphi       &
            , fx, fy, fz       &
            , fx, fy, fz )

       beta = dsqrt(beta)

       if ( mstep.LT.krylov_dim ) then
          
          !----- Normalizes the orthogonal residual and uses it as a new Krylov vector -----!
          
          H(mstep+1,mstep) = beta
          
       endif
       
       beta = 1.0D+00/beta
       
       fx   = beta*fx
       fy   = beta*fy
       fz   = beta*fz
       fp   = beta*fp
       fphi = beta*fphi

       call exchange2d(fphi,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
       call exchange2d(fx,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(fy,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
       call exchange2d(fz, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,comm2d_quasiperiodic)
       call exchange2d(fp,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

       Qx(:,:,:,mstep+1)   = fx
       Qy(:,:,:,mstep+1)   = fy
       Qz(:,:,:,mstep+1)   = fz
       Qp(:,:,:,mstep+1)   = fp
       Qphi(:,:,:,mstep+1) = fphi
       
    end do Arnoldi
    
    return
  end subroutine Arnoldi_factorization

end module tpls_solver
