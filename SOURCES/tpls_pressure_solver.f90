module tpls_pressure_solver

  use tpls_levelset
  implicit none
  private
  public :: poisson_solver, srj_weights

contains





!-----------------------------------------------------------------------------------------------------




  subroutine poisson_solver(x, b, maxiter, tol)

    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, optional, intent(in) :: tol
    integer, optional, intent(in) :: maxiter

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(in)    :: b
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(inout) :: x

    !----- Miscellaneous -----!

    double precision :: residual, erreur, norme, weight, tolerance
    integer          :: it, sub_it, iterations
    integer          :: i, j, k, n

    if ( present(tol) ) then
       tolerance = tol
    else
       tolerance = 1e-6
    endif

    if ( present(maxiter) ) then
       iterations = maxiter
    else
       iterations = 100000
    endif
    
    n = maxn
    t_temp = mpi_wtime()

!!$    call conjugate_gradient(residual, x, -b, iterations, tolerance)
!!$    call pressure_bc(x)

    !----- SCHEDULED RELAXATION JACOBI METHOD -----!
    
    if ( maxm == 2) then
       call norm(norme, b(sx:ex,1,1:n-1), size(b(sx:ex,1,1:n-1)))
    else
       call norm(norme, b(sx:ex,sy:ey,1:n-1), size(b(sx:ex,sy:ey,1:n-1)))
    endif
    
    !----- Begining of the Poisson solver -----!
    
    residual = 1.0D+00
    
    Outer_iteration : do it = 1, iterations, max_srj
       Inner_iteration : do sub_it = 1, max_srj
          
          !----- Skip the big over-relaxation step as a starting point -----!
          !      Note: this trick tends to significantly decrease the cpu
          !            time required for convergence.
          
          if ( (it == 1) .and. (sub_it==1) ) then
             weight = relax_srj(sub_it+1)
          else
             weight = relax_srj(sub_it)
          endif
          
          call srj_iteration(residual, x, b, weight)
          residual = min(residual, residual/norme)
          
          !----- Perform various checks -----!
          
          if ( residual.LT.tolerance ) EXIT Outer_iteration
          if ( it+sub_it .gt. iterations ) EXIT Outer_iteration
          if ( isnan(residual) ) STOP 'The Poisson solver blew up!'
          
       enddo Inner_iteration
    enddo Outer_iteration
    
    if ( my_id == master_id ) then
       write(*,*)'Iteration : ', istep, 'p---residual is ', residual, '(', it + sub_it, ')'
       write(*,*)
    end if
    
    time_pressure = time_pressure + (mpi_wtime() - t_temp)
    
  end subroutine Poisson_solver



  
!-----------------------------------------------------------------------------------------------------





  subroutine srj_iteration(residual, pres, rhs, relaxation)
    
    use tpls_constants
    use tpls_mpi
    use tpls_maths
    use mpi
    implicit none
    
    !----- Inputs and outputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(inout) :: pres
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(in) :: rhs
    double precision, intent(in) :: relaxation
    double precision, intent(out) :: residual

    !----- Intermediate arrays -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn) :: residual_array
    integer :: n, i, j, k

    n = maxn
    
    residual = 0.0D+00
    residual_array = 0.0D+00

    if ( maxm==2) then
       
       residual_array(sx:ex,sy:ey,1:n-1) = (1.0D+00/6.0D+00) * (        &
            pres(sx+1:ex+1,sy:ey,1:n-1) + pres(sx-1:ex-1,sy:ey,1:n-1) &
            + pres(sx:ex,sy+1:ey+1,1:n-1) + pres(sx:ex,sy-1:ey-1,1:n-1) &
            + pres(sx:ex,sy:ey,2:n)       + pres(sx:ex,sy:ey,0:n-2)     &
            - (dz**2.)*rhs(sx:ex,sy:ey,1:n-1)                           )
       
       residual_array = residual_array-pres
       call norm(residual, residual_array(sx:ex,sy:ey,1:n-1), &
            size(residual_array(sx:ex,sy:ey,1:n-1)))
       residual = residual * (6/dz**2)
       
    else if ( maxm == 2 ) then
       
       residual_array(sx:ex,1,1:n-1) = .25 * (        &
            pres(sx+1:ex+1,1,1:n-1) + pres(sx-1:ex-1,1,1:n-1) &
            + pres(sx:ex,1,2:n)       + pres(sx:ex,1,0:n-2)     &
            - (dz**2.)*rhs(sx:ex,1,1:n-1)                           )
       
       residual_array(sx:ex,1,1:n-1) = residual_array(sx:ex,1,1:n-1)-pres(sx:ex,1,1:n-1)
       call norm(residual, residual_array(sx:ex,1,1:n-1), &
            size(residual_array(sx:ex,1,1:n-1)))
       residual = residual * (4/dz**2)
       do j = sy-1, ey+1
          residual_array(sx:ex,j,1:n-1) = residual_array(sx:ex,1,1:n-1)
       enddo
       
    endif
    
    pres = pres + relaxation*residual_array
    call pressure_bc(pres)
    
    return
  end subroutine srj_iteration


  !-----------------------------------------------------------------------------------------------------
  
  
  subroutine pressure_bc(pressure)
    
    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none
    
    !----- Input/Output -----!

    double precision, dimension(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p), intent(inout) :: pressure

    !----- Neuman boundary condition dp/dx at the inflow -----!

    if ( sx==1 ) then
       pressure(0, :, :) = pressure(1, :, :)
    endif

    !----- Neumann boundary condition dp/dz at the upper and lower walls -----!

    pressure(:, :, 0) = pressure(:, :, 1)
    pressure(:, :, maxn) = pressure(:, :, maxn-1)

    !----- Dirichlet boundary condition p=0 at the outflow -----!

    if ( ex==ex_max ) then
       pressure(ex_max+1, :, :) = -pressure(ex_max, :, :)
    endif

    call exchange2d(pressure,stride_p_xz,stride_p_yz,neighbours, &
         ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
    
    return
  end subroutine pressure_bc

  !-----------------------------------------------------------------------------------------------------
  
  subroutine SRJ_weights(weights)
    
    use tpls_constants
    implicit none
    
    double precision, dimension(1000000), intent(out) :: weights
    
    double precision, dimension(15) :: W, dis, beta
    integer, dimension(15) :: Q

    double precision, parameter :: pi = 4.D0*atan(1.D0)
    double precision :: kmin, ww
    double precision, dimension(:), allocatable :: k_wavenumber, G
    integer :: i, j, k, m, n
    integer :: istart, iend, index, counter
    integer, dimension(1) :: max_loc, min_loc

    n = max(maxl, maxm, maxn)
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

    if ( (N.GE.64).AND.(N.LT.128) ) then
       
       W(1) = 1604.55D+00
       W(2) = 1236.6D+00
       W(3) = 777.72D+00
       W(4) = 429.57D+00
       W(5) = 220.699D+00
       W(6) = 109.268D+00
       W(7) = 53.1653D+00
       W(8) = 25.7023D+00
       W(9) = 12.4395D+00
       W(10) = 6.06839D+00
       W(11) = 3.77684D+00
       W(12) = 2.26342D+00
       W(13) = 1.17188D+00
       W(14) = 0.697364D+00
       W(15) = 0.519746D+00
       
       beta(1) = 0.00324844D+00
       beta(2) = 0.00375019D+00
       beta(3) = 0.00483085D+00
       beta(4) = 0.00665688D+00
       beta(5) = 0.00950942D+00
       beta(6) = 0.0138266D+00
       beta(7) = 0.0202681D+00
       beta(8) = 0.0298105D+00
       beta(9) = 0.0439172D+00
       beta(10) = 0.0661899D+00
       beta(11) = 0.0257826D+00
       beta(12) = 0.120006D+00
       beta(13) = 0.167699D+00
       beta(14) = 0.222552D+00
       beta(15) = 0.261952D+00
       
       Q = floor(beta/beta(1))
       
    elseif ( (N.GE.128).AND.(N.LT.256) ) then
       
       W(1) = 6371.83D+00
       W(2) = 4666.19D+00
       W(3) = 2709.05D+00
       W(4) = 1370.56D+00
       W(5) = 646.134D+00
       W(6) = 294.622D+00
       W(7) = 132.367D+00
       W(8) = 59.1371D+00
       W(9) = 26.4159D+00
       W(10) = 11.8529D+00
       W(11) = 5.94244D+00
       W(12) = 2.89717D+00
       W(13) = 1.36449D+00
       W(14) = 0.743778D+00
       W(15) = 0.523863D+00

       beta(1) = 0.00167178D+00
       beta(2) = 0.00198802D+00
       beta(3) = 0.00268004D+00
       beta(4) = 0.00387851D+00
       beta(5) = 0.00580983D+00
       beta(6) = 0.00883865D+00
       beta(7) = 0.0135361D+00
       beta(8) = 0.0207864D+00
       beta(9) = 0.031965D+00
       beta(10) = 0.0496738D+00
       beta(11) = 0.0510383D+00
       beta(12) = 0.113836D+00
       beta(13) = 0.169489D+00
       beta(14) = 0.236713D+00
       beta(15) = 0.288096D+00
       
       Q = floor(beta/beta(1))
       
    elseif ( (N.GE.256).AND.(N.LT.512) ) then
       
       W(1) = 25234.4D+00
       W(2) = 17233.0D+00
       W(3) = 9009.26D+00
       W(4) = 4080.06D+00
       W(5) = 1729.54D+00
       W(6) = 712.786D+00
       W(7) = 290.428D+00
       W(8) = 117.862D+00
       W(9) = 47.8274D+00
       W(10) = 19.4772D+00
       W(11) = 8.32423D+00
       W(12) = 3.5472D+00
       W(13) = 1.54899D+00
       W(14) = 0.785632D+00
       W(15) = 0.527435D+00
       
       beta(1) = 0.000861253D+00
       beta(2) = 0.0010671D+00
       beta(3) = 0.00152771D+00
       beta(4) = 0.00235298D+00
       beta(5) = 0.00373991D+00
       beta(6) = 0.00601951D+00
       beta(7) = 0.00973541D+00
       beta(8) = 0.0157724D+00
       beta(9) = 0.0255664D+00
       beta(10) = 0.0415604D+00
       beta(11) = 0.0589968D+00
       beta(12) = 0.108078D+00
       beta(13) = 0.169214D+00
       beta(14) = 0.246607D+00
       beta(15) = 0.308902D+00
       
       Q = floor(beta/beta(1))
       
    elseif ( ( N.GE.512 ) .AND. ( N.LT.1024 ) ) then
       
       W(1) = 99805.2D+00
       W(2) = 63101.3D+00
       W(3) = 29545.0D+00
       W(4) = 11959.4D+00
       W(5) = 4558.78D+00
       W(6) = 1698.18D+00
       W(7) = 627.242D+00
       W(8) = 231.042D+00
       W(9) = 85.1043D+00
       W(10) = 31.433D+00
       W(11) = 11.8839D+00
       W(12) = 4.53525D+00
       W(13) = 1.81056D+00
       W(14) = 0.841402D+00
       W(15) = 0.532005D+00
       
       beta(1) = 0.000435073D+00
       beta(2) = 0.000564418D+00
       beta(3) = 0.000861316D+00
       beta(4) = 0.00141386D+00
       beta(5) = 0.00238609D+00
       beta(6) = 0.00406664D+00
       beta(7) = 0.00695417D+00
       beta(8) = 0.0119047D+00
       beta(9) = 0.0203829D+00
       beta(10) = 0.0349166D+00
       beta(11) = 0.0569284D+00
       beta(12) = 0.101306D+00
       beta(13) = 0.167425D+00
       beta(14) = 0.256851D+00
       beta(15) = 0.333604D+00
       
       Q = floor(beta/beta(1))
       
    elseif ( (N.GE.1024) .and. (n .lt. 2048) ) then
       
       W(1) = 394347.D+00
       W(2) = 229799.D+00
       W(3) = 96276.D+00
       W(4) = 34921.9D+00
       W(5) = 12008.9D+00
       W(6) = 4053.99D+00
       W(7) = 1360.11D+00
       W(8) = 455.47D+00
       W(9) = 152.531D+00
       W(10) = 51.1795D+00
       W(11) = 17.388D+00
       W(12) = 5.98513D+00
       W(13) = 2.16481D+00
       W(14) = 0.91159D+00
       W(15) = 0.537479D+00
       
       beta(1) = 0.000215354D+00
       beta(2) = 0.000293494D+00
       beta(3) = 0.00047792D+00
       beta(4) = 0.00083541D+00
       beta(5) = 0.00149547D+00
       beta(6) = 0.0026972D+00
       beta(7) = 0.00487579D+00
       beta(8) = 0.00881986D+00
       beta(9) = 0.0159551D+00
       beta(10) = 0.0288594D+00
       beta(11) = 0.0511968D+00
       beta(12) = 0.0935117D+00
       beta(13) = 0.1673741D+00
       beta(14) = 0.266199D+00
       beta(15) = 0.360827D+00
       
       Q = floor(beta/beta(1))

    elseif ( (n .ge. 2048) .and. (n .lt. 4096) ) then

       W(1) = 1556575
       W(2) = 832736
       W(3) = 312142
       W(4) = 101721
       W(5) = 31639.4
       W(6) = 9698.49
       W(7) = 2959.69
       W(8) = 902.95
       W(9) = 274.961
       W(10) = 83.9203
       W(11) = 25.799
       W(12) = 8.03399
       W(13) = 2.62374
       W(14) = 0.99542
       W(15) = 0.543653
       
       beta(1) = 0.000104643
       beta(2) = 0.000150277
       beta(3) = 0.000261282
       beta(4) = 0.000485959
       beta(5) = 0.000922098
       beta(6) = 0.0017595
       beta(7) = 0.00336258
       beta(8) = 0.00642874
       beta(9) = 0.0122909
       beta(10) = 0.0234935
       beta(11) = 0.0445452
       beta(12) = 0.0851626
       beta(13) = 0.158288
       beta(14) = 0.273725
       beta(15) = 0.38902
       
       Q = floor(beta/beta(1))

    elseif ( (n .ge. 4096) .and. (n .lt. 8192) ) then

       W(1) = 6045072
       W(2) = 2698661
       W(3) = 811234
       W(4) = 214918
       W(5) = 55022.7
       W(6) = 13963.2
       W(7) = 3535.64
       W(8) = 894.898
       W(9) = 226.17
       W(10) = 57.5285
       W(11) = 14.7527
       W(12) = 3.93026
       W(13) = 1.20411
       W(14) = 0.557585
       
       beta(1) = 0.0000505714
       beta(2) = 0.0000812954
       beta(3) = 0.000161234
       beta(4) = 0.000338826
       beta(5) = 0.000721677
       beta(6) = 0.00154274
       beta(7) = 0.00329580
       beta(8) = 0.00704624
       beta(9) = 0.0150625
       beta(10) = 0.0321750
       beta(11) = 0.0684563
       beta(12) = 0.143616
       beta(13) = 0.282283
       beta(14) = 0.445170
       
       Q = floor(beta/beta(1))
       W(15) = -1000
       Q(15) = 0

    elseif ( (n .ge. 8192) .and. (n .lt. 16384) ) then

       W(1) = 20841177
       W(2) = 4339863
       W(3) = 589668
       W(4) = 75210.5
       W(5) = 9514.64
       W(6) = 1202.61
       W(7) = 152.183
       W(8) = 19.4605
       W(9) = 2.70028
       W(10) = 0.624451
       
       beta(1) = 0.0000219770
       beta(2) = 0.0000581897
       beta(3) = 0.000189695
       beta(4) = 0.000632223
       beta(5) = 0.00211144
       beta(6) = 0.00705278
       beta(7) = 0.0235524
       beta(8) = 0.0784280
       beta(9) = 0.253403
       beta(10) = 0.634551
       
       Q = floor(beta/beta(1))
       W(11:15) = -1000
       Q(11:15) = 0

    elseif ( (n .gt. 16384) ) then

       W(1) = 67375460
       W(2) = 6181549.54
       W(3) = 397871.724
       W(4) = 24979.4457
       W(5) = 1566.01561
       W(6) = 98.3988686
       W(7) = 6.42559672
       W(8) = 0.70552038
       
       beta(1) = 0.00000788702
       beta(2) = 0.0000379175
       beta(3) = 0.000202091
       beta(4) = 0.00110994
       beta(5) = 0.00609737
       beta(6) = 0.0334824
       beta(7) = 0.181878
       beta(8) = 0.777185
       
       Q = floor(beta/beta(1))
       W(9:15) = -1000
       Q(9:15) = 0
       
    endif
    
    !----- Compute the appropriate sequence of over- and under-relaxation -----

    m = sum(Q)
    max_srj = m
    weights = 0.0D+00

    weights(1) = W(1)
    Q(1) = Q(1) - 1
    G = G*abs(1.D+00 - k_wavenumber*weights(1))

    counter = 2

    do while (sum(Q).NE.0)

       do i = 1,15
          if ( Q(i).EQ.0 ) W(i) = -1000
       enddo

       max_loc = maxloc(G)
       ww = 1.0D+00/k_wavenumber(max_loc(1))
       dis = abs(W - ww)
       min_loc = minloc(dis)
       weights(counter) = W(min_loc(1))
       G = G*abs(1.D+00 - k_wavenumber*weights(counter))
       Q(min_loc(1)) = Q(min_loc(1)) - 1
       counter = counter + 1

    enddo

    open( unit=10, file='srj.dat')
    do i = 1, m
       write(10, *) weights(i)
    enddo
    close(10)
    
  end subroutine SRJ_weights

  subroutine laplace_operator(Lu, u)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Input -----!

    double precision, dimension(sx-1:ex+1, sy-1:ey+1, 0:maxn), intent(in) :: u

    !----- Output -----!

    double precision, dimension(sx-1:ex+1, sy-1:ey+1, 0:maxn), intent(out) :: Lu

    !----- Miscellaneous -----!

    integer :: i, j, k, n
    double precision, dimension(sx-1:ex+1, sy-1:ey+1, 0:maxn) :: work1, work2, work3

    n = maxn-1
    
    work1 = 0
    work2 = 0
    work3 = 0
    Lu = 0
    
    !----- d2u/dx2 -----!

    work1(sx:ex, :, :) = u(sx+1:ex+1, :, :) - 2*u(sx:ex, :, :) + u(sx-1:ex-1, :, :)
    work1 = work1/dx**2

    !----- d2u/dy2 -----!

    work2(:, sy:ey, :) = u(:, sy+1:ey-1, :) - 2*u(:, sy:ey, :) + u(:, sy-1:ey-1, :)
    work2 = work2/dy**2

    !----- d2u/dz2 -----!

    work3(:, :, 1:n) = u(:, :, 2:n+1) - 2*u(:, :, 1:n) + u(:, :, 0:n-1)     
    work3 = work3/dz**2

    Lu = -(work1 + work2 + work3)
    
    return
  end subroutine laplace_operator

  subroutine conjugate_gradient(residual, u, b, maxiter, tol)

    use tpls_constants
    use tpls_maths
    use tpls_mpi
    use mpi
    implicit none

    !----- Input/Output -----!

    double precision, dimension(sx-1:ex+1, sy-1:ey+1, 0:maxn), intent(inout) :: u
    double precision, dimension(sx-1:ex+1, sy-1:ey+1, 0:maxn), intent(in)    :: b
    double precision, intent(in) :: tol
    integer, intent(in) :: maxiter
    double precision, intent(out) :: residual

    !----- Miscellaneous -----!

    integer :: i, j, k, n, size_vec
    double precision :: alpha, beta, dummy, normb, rTr, pTLp
    double precision, dimension(sx-1:ex+1, sy-1:ey+1, 0:maxn) :: r, p, Lp, dummy_vec, r_old

    n = maxn-1
    residual = 1.
    r = 0
    r_old = 0
    p = 0
    size_vec = size(r(sx:ex, sy:ey, 1:n))
    call norm(normb, b(sx:ex, sy:ey, 1:n), size_vec)
    call laplace_operator(Lp, u)
    r = b - Lp
    p = r
    
    cg_loop : do k = 0, maxiter
       if (residual .lt. tol) then
          exit cg_loop
       else
          r_old = r
          call pressure_bc(p) 
          call laplace_operator(Lp, p)
          
          rTr = euclidean_scalar_product(r(sx:ex, sy:ey, 1:n), &
               r(sx:ex, sy:ey, 1:n), &
               size_vec )
          pTLp = euclidean_scalar_product(p(sx:ex, sy:ey, 1:n), &
               Lp(sx:ex, sy:ey, 1:n), &
               size_vec )
          
          alpha = rTr/pTLp
          
          u = u + alpha*p
          r = r - alpha*Lp
          call exchange2d(r,stride_p_xz,stride_p_yz,neighbours, &
               ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
 
          residual = euclidean_scalar_product(r(sx:ex, sy:ey, 1:n), &
               r(sx:ex, sy:ey, 1:n), &
               size_vec )
          dummy = euclidean_scalar_product(r(sx:ex, sy:ey, 1:n), &
               r_old(sx:ex, sy:ey, 1:n), &
               size_vec )
          beta = (residual-dummy)/rTr
          
          p = r + beta*p
          
          residual = sqrt(residual)/normb
          
       endif
       
    enddo cg_loop

    if (residual .lt. tol) then
       if ( my_id == master_id ) then
          write(*,*)'Iteration : ', istep, 'p---residual is ', residual, '(', k, ')'
          write(*,*)
       endif
    else
       if ( my_id == master_id ) then
          write(*,*)'Iteration : ', istep, 'p---residual is ', residual, '(', k, ')'
          write(*,*) 'CG FAILED TO CONVERGE'
       endif
    endif
    call pressure_bc(u)
      
    return
  end subroutine conjugate_gradient


end module tpls_pressure_solver
