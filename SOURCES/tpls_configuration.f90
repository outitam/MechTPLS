module tpls_configuration

  implicit none

contains

  ! ****************************************************************************************
  subroutine SRJ_weights
    
    use tpls_constants
    implicit none
    
    double precision, dimension(5) :: W, dis
    integer, dimension(5) :: Q

    double precision, parameter :: pi = 4.D0*atan(1.D0)
    double precision :: kmin, ww
    double precision, dimension(:), allocatable :: k_wavenumber, G
    integer :: i, j, k, m, n
    integer :: istart, iend, index, counter
    integer, dimension(1) :: max_loc, min_loc

    n = int((maxn - 1)/Lz)
    kmin = sin(pi/2.0D+00/N)**2.0D+00

    k = floor(2.d0/kmin)
    allocate(k_wavenumber(k))
    allocate(G(k))
    G = 1.0D+00

    ww = 0.0D+00

    do i = 1, k
       k_wavenumber(i) = i*kmin
    enddo

    !----- Weights for P=5 SRJ schemes -----!

    !Note : For more detail, see Yang & Mittall, "Acceleration of the Jacobi
    ! iterative method by factors exceeding 100 using scheduled relaxation",
    ! J. Comp. Phys., vol. 274, pp 695--708, 2014

    if ( (N.GE.16).AND.(N.LT.32) ) then

       W(1) = 88.190D+00
       W(2) = 30.122D+00
       W(3) = 6.8843D+00
       W(4) = 1.6008D+00
       W(5) = 0.58003D+00

       Q(1) = 1
       Q(2) = 2
       Q(3) = 5
       Q(4) = 12
       Q(5) = 23

    elseif ( (N.GE.32).AND.(N.LT.64) ) then

       W(1) = 330.57D+00
       W(2) = 82.172D+00
       W(3) = 13.441D+00
       W(4) = 2.2402D+00
       W(5) = 0.6081D+00

       Q(1) = 1
       Q(2) = 2
       Q(3) = 7
       Q(4) = 20
       Q(5) = 46

    elseif ( (N.GE.64).AND.(N.LT.128) ) then

       W(1) = 1228.8D+00
       W(2) = 220.14D+00
       W(3) = 26.168D+00
       W(4) = 3.1668D+00
       W(5) = 0.6389D+00

       Q(1) = 1
       Q(2) = 3
       Q(3) = 10
       Q(4) = 38
       Q(5) = 106

    elseif ( (N.GE.128).AND.(N.LT.256) ) then

       W(1) = 4522.0D+00
       W(2) = 580.86D+00
       W(3) = 50.729D+00
       W(4) = 4.5018D+00
       W(5) = 0.67161D+00

       Q(1) = 1
       Q(2) = 3
       Q(3) = 16
       Q(4) = 73
       Q(5) = 250

    elseif ( (N.GE.256).AND.(N.LT.512) ) then

       W(1) = 16459D+00
       W(2) = 1513.4D+00
       W(3) = 97.832D+00
       W(4) = 6.4111D+00
       W(5) = 0.70531D+00

       Q(1) = 1
       Q(2) = 4
       Q(3) = 26
       Q(4) = 142
       Q(5) = 605

    elseif ( N.GE.512 ) then

       W(1) = 59226D+00
       W(2) = 3900.56D+00
       W(3) = 187.53D+00
       W(4) = 9.1194D+00
       W(5) = 0.73905D+00

       Q(1) = 1
       Q(2) = 6
       Q(3) = 40
       Q(4) = 277
       Q(5) = 1500

    endif

    !----- Compute the appropriate sequence of over- and under-relaxation -----

    m = sum(Q)
    max_srj = m
    allocate(relax_srj(m))
    relax_srj = 0.0D+00

    relax_srj(1) = W(1)
    Q(1) = Q(1) - 1
    G = G*abs(1.D+00 - k_wavenumber*relax_srj(1))

    counter = 2

    do while (sum(Q).NE.0)

       do i = 1,5
          if ( Q(i).EQ.0 ) W(i) = -1000
       enddo

       max_loc = maxloc(G)
       ww = 1.0D+00/k_wavenumber(max_loc(1))
       dis = abs(W - ww)
       min_loc = minloc(dis)
       relax_srj(counter) = W(min_loc(1))
       G = G*abs(1.D+00 - k_wavenumber*relax_srj(counter))
       Q(min_loc(1)) = Q(min_loc(1)) - 1
       counter = counter + 1

    enddo

  end subroutine SRJ_weights


  ! ****************************************************************************************
  subroutine read_parameters()

    use tpls_constants
    use tpls_mpi
    implicit none

    integer :: i, restart_int, sfd_int

    open( unit = 10 , file = 'tpls_config.opt' , status = 'old' )

    do i = 1,4
       read(10,*)
    enddo

    read(10,*) dt
    read(10,*) nsteps
    read(10,*) iostep
    read(10,*) Re
    read(10,*) mu_plus
    read(10,*) rho_plus
    read(10,*) mu_minus
    read(10,*) rho_minus
    read(10,*) height
    read(10,*) Scap
    read(10,*) Grav
    read(10,*) gz

    do i = 1,5
       read(10,*)
    enddo

    read(10,*) maxl
    read(10,*) maxm
    read(10,*) maxn
    read(10,*) num_procs_x
    read(10,*) num_procs_y
    read(10,*) num_procs_z

    do i = 1,5
       read(10,*)
    enddo

    read(10,*) smooth_width
    read(10,*) tolerance_levelset
    read(10,*) max_iteration_levelset
    read(10,*) tolerance_poisson
    read(10,*) max_iteration_poisson
    read(10,*) tolerance_helmholtz
    read(10,*) max_iteration_helmholtz

    do i = 1,5
       read(10,*)
    enddo

    read(10,*) if_implicit
    read(10,*) sfd_int
    read(10,*) restart_int
    read(10,*) if_exact_restart_save
    read(10,*) if_record_perturbation_norm
    read(10,*) if_linear_stability

    do i = 1,5
       read(10,*)
    enddo

    read(10,*) restart_handle

    close(10)

    if_restart=.false.
    if_exact_restart=.false.
    if ((restart_int .lt. 0) .or. (restart_int .gt. 2)) then
        if (my_id==master_id) then
        write(*,*) "Illegal nr. in Restart field of tpls_config. opt! Choose either 0, 1 or 2"
        stop
        endif
    elseif (restart_int>0) then
        if_restart=.true.
        if (restart_int>1) if_exact_restart=.true.
    endif
  
    
    if_sfd=.false.
    if_exact_restart_sfd=.false.
    if ((sfd_int .lt. 0) .or. (sfd_int .gt. 2)) then
        if (my_id==master_id) then
        write(*,*) "Illegal nr. in SFD field of tpls_config. opt! Choose either 0, 1 or 2"
        stop
        endif
    elseif (sfd_int>0) then 
        if_sfd=.true.
        if (sfd_int>1) if_exact_restart_sfd=.true.
    endif

    density_matched = .false.
    if( rho_plus .EQ. rho_minus ) density_matched = .true.
    
    viscosity_matched = .false.
    if( mu_plus .EQ. mu_minus ) viscosity_matched = .true.

    Lz = 4.0D+00

    dz = Lz/dble(maxn-1)
    dx = dz
    dy = dz

    Lx = (maxl-1)*dx
    Ly = (maxm-1)*dy
    
    istep = 0
    
    smooth_width = smooth_width*dx
    
  end subroutine read_parameters

  ! ****************************************************************************************
  subroutine dataload_channel_binary(vx,vy,vz,pr,phi)
    
    use tpls_constants
    use tpls_mpi
    use tpls_userchk
    use mpi
    implicit none
    
    !----- Inputs -----!
    
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(out) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(out) :: vz
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(out) :: phi, pr
    
    !----- Miscellaneous -----!
    
    double precision, dimension(0:maxl,0:maxm,0:maxn-2) :: u, v
    double precision, dimension(0:maxl,0:maxm,0:maxn-1) :: w
    double precision, dimension(0:maxl,0:maxm,0:maxn)   :: pr_, phi_

    integer :: i, j, k, ifich

    vx   = 0.0D+00
    vy   = 0.0D+00
    vz   = 0.0D+00

    phi  = 0.0D+00
    pr   = 0.0D+00

    u    = 0.0D+00
    v    = 0.0D+00
    w    = 0.0D+00

    pr_  = 0.0D+00
    phi_ = 0.0D+00


    if ( if_restart ) then
     
       !----- Read data for possible restart -----!
       
       ifich = 10
       
       open( unit = ifich, file = restart_handle, form = 'unformatted')
       
       read(ifich) u
       read(ifich) v
       read(ifich) w
       read(ifich) pr_
       read(ifich) phi_
       
       close(ifich)
       
    elseif ( .not. if_restart ) then
       
       call user_initial_conditions(u,v,w,pr_,phi_)
       
    endif
    
    pr(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)   = pr_(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  = phi_(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    vx(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    vy(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    vz(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   = w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)
    
    call mpi_barrier(mpi_comm_world,ierr)
    
    return
  end subroutine dataload_channel_binary

  ! ****************************************************************************************
    
  subroutine dataload_exact_conv(conv1_u,conv1_v,conv1_w,csf_u1,csf_v1,csf_w1,   &
                              diffusion1_u,diffusion1_v,diffusion1_w)
    use tpls_constants
    use tpls_mpi
    use tpls_userchk
    use mpi
    implicit none
    
    !----- Inputs -----!
    
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(out) :: csf_u1, csf_v1
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(out) :: diffusion1_u, diffusion1_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(out) :: conv1_u, conv1_v
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(out) :: csf_w1
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(out) :: diffusion1_w
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(out) :: conv1_w
    
    !----- Miscellaneous -----!
    

    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2) :: conv1_u_gl, conv1_v_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2) :: csf_u1_gl, csf_v1_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2) :: diffusion1_u_gl, diffusion1_v_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1) :: conv1_w_gl, csf_w1_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1) :: diffusion1_w_gl

    integer :: i, j, k, ifich
    character(len=100) :: filename

    conv1_u=0.d0
    conv1_v=0.d0
    conv1_w=0.d0
    csf_u1=0.d0
    csf_v1=0.d0
    csf_w1=0.d0
    diffusion1_u=0.d0
    diffusion1_v=0.d0
    diffusion1_w=0.d0

    if ( if_restart ) then
       
       !----- Read data for restart -----!
       
       ifich = 10
       
       filename = 'fieldexactrestart_conv.bin'
       
       open( unit = ifich, file = filename, form = 'unformatted')
       
     read(ifich) conv1_u_gl
     read(ifich) csf_u1_gl
     read(ifich) diffusion1_u_gl
     read(ifich) conv1_v_gl
     read(ifich) csf_v1_gl
     read(ifich) diffusion1_v_gl
     read(ifich) conv1_w_gl
     read(ifich) csf_w1_gl
     read(ifich) diffusion1_w_gl
       
     close(ifich)
       
    endif
    
    conv1_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = conv1_u_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    conv1_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = conv1_v_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    csf_u1(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = csf_u1_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    csf_v1(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = csf_v1_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    diffusion1_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = &
     diffusion1_u_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    diffusion1_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = &
     diffusion1_v_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    conv1_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   = conv1_w_gl(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)
    csf_w1(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   = csf_w1_gl(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)
    diffusion1_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   = &
     diffusion1_w_gl(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)
    
    call mpi_barrier(mpi_comm_world,ierr)
    
    return
  end subroutine dataload_exact_conv

  ! ****************************************************************************************
    
  subroutine dataload_exact_sfd(vx_bar,vy_bar,vz_bar,fx_sfd,fy_sfd,fz_sfd)   
    
    use tpls_constants
    use tpls_mpi
    use tpls_userchk
    use mpi
    implicit none
    
    !----- Inputs -----!
    
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(out) :: vx_bar, vy_bar, fx_sfd, fy_sfd
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   intent(out) :: vz_bar, fz_sfd
    
    !----- Miscellaneous -----!
    

    double precision, dimension(0:maxl, 0:maxm, 0:maxn-2) :: vx_bar_gl, vy_bar_gl, fx_sfd_gl, fy_sfd_gl
    double precision, dimension(0:maxl, 0:maxm, 0:maxn-1) :: vz_bar_gl, fz_sfd_gl

    integer :: i, j, k, ifich
    character(len=100) :: filename

    vx_bar=0.d0
    vy_bar=0.d0
    vz_bar=0.d0
    fx_sfd=0.d0
    fy_sfd=0.d0
    fz_sfd=0.d0

    if ( if_restart ) then
       
       !----- Read data for restart -----!
       
       ifich = 103
       
       filename = 'fieldexactrestart_sfd.bin'
       
       open( unit = ifich, file = filename, form = 'unformatted')
       
     read(ifich) vx_bar_gl
     read(ifich) vy_bar_gl
     read(ifich) vz_bar_gl
     read(ifich) fx_sfd_gl
     read(ifich) fy_sfd_gl
     read(ifich) fz_sfd_gl
       
     close(ifich)
       
    endif
    
    vx_bar(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = vx_bar_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    vy_bar(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = vy_bar_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    fx_sfd(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = fx_sfd_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    fy_sfd(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) = fy_sfd_gl(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    vz_bar(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   = vz_bar_gl(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)
    fz_sfd(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   = fz_sfd_gl(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)
    
    call mpi_barrier(mpi_comm_world,ierr)
    
    return
  end subroutine dataload_exact_sfd

  ! ****************************************************************************************
end module tpls_configuration
