module tpls_maths

  implicit none

contains

    subroutine divergence(div,vx,vy,vz)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none
    
    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(in) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-1), intent(in) :: vz

    !----- Output -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn), intent(out) :: div

    !----- Intermediate arrays -----!

    integer :: i,j,k
    double precision :: scalar_plusx,scalar_minusx

    div = 0.0D+00
    
    div(sx:ex,sy:ey,1:maxn-1) =                                           &
           ( vx(sx+1:ex+1,sy:ey,0:maxn-2) - vx(sx:ex,sy:ey,0:maxn-2) )/dx &
         + ( vy(sx:ex,sy+1:ey+1,0:maxn-2) - vy(sx:ex,sy:ey,0:maxn-2) )/dy &
         + ( vz(sx:ex,sy:ey,1:maxn-1)     - vz(sx:ex,sy:ey,0:maxn-2) )/dz
    
    return
  end subroutine divergence



  subroutine Least_Squares(A,x,b,nr,nc)

    implicit none

    integer :: nr, nc
    double precision, dimension(nr,nc) :: A
    double precision, dimension(nr)    :: b
    double precision, dimension(nc)    :: x

    !----- Variables for LAPACK DGELS -----!

    character*1 :: TRANS = 'N'
    integer :: m
    integer :: n
    integer :: nrhs = 1
    double precision, dimension(nr,nc) :: A_
    integer :: LDA
    integer :: LDB
    integer :: LWORK
    double precision, dimension(:), allocatable :: WORK
    integer INFO

    !----- Solving the least square problem -----!

    m = nr
    n = nc

    A_ = A
    LDA = max(1,m)
    LDB = max(1,m,n)
    LWORK = max( 1, m*n + max(m*n,nrhs) )

    allocate(WORK(max(1,LWORK)))

    call DGELS( TRANS, M, N, NRHS, A_ , LDA, B, LDB, WORK, LWORK, INFO )

    deallocate(WORK)

    x = 0.0D+00
    x = b(1:nc)

  end subroutine Least_Squares

  subroutine Hessenberg_decomposition(A,Q,n)

    implicit none

    !----- Required variables for dgehrd -----

    integer :: n, ILO, IHI, LDA, LWORK, INFO
    double precision, dimension(n,n) :: A
    double precision, dimension(n-1) :: TAU
    double precision, dimension(n)   :: WORK

    !----- Required variables for dorghr -----

    double precision, dimension(n,n) :: Q

    !----- Hessenberg decomposition -----

    ILO = 1
    IHI = n
    LDA = max(1,n)
    LWORK = max(1,n)

    call dgehrd(n,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)

    !----- Compute the unitary transformation matrix Q -----

    Q = A
    LWORK = IHI-ILO
    call dorghr(n,ILO,IHI,Q,LDA,TAU,WORK,LWORK,INFO)

    return
  end subroutine Hessenberg_decomposition

  subroutine eigendecomposition(A,n,VP,VR)

    implicit none
    
    !     ----- Required variables for zgeev (eigenpairs computations) -----
    
    integer    :: INFO
    integer    :: n, LWORK
    double precision, dimension(n,n) :: A
    double precision, dimension(2*n) :: RWORK
    double complex, dimension(1,n) :: VL
    double complex, dimension(n,n) :: VR
    double complex, dimension(4*n) :: WORK
    double complex, dimension(n)   :: VP
    
    !     ----- Computation of eig(A) -----
    
    LWORK = 4*n
    
    call ZGEEV('N','V',n,A*(1.,0.),n,VP,VL,1,VR,n, &
         WORK,LWORK,RWORK,INFO)
    
    return
  end subroutine eigendecomposition

  subroutine matvec(FP_Cx,Qx,FP,k_dim,sx,ex,sy,ey,sz,ez)

    implicit none

    integer :: k_dim, i, sx, ex, sy, ey, sz, ez
  
    double complex, dimension(sx-1:ex+1,sy-1:ey+1,sz:ez)               :: FP_Cx
    double complex, dimension(k_dim)                                   :: FP
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz:ez,1:k_dim) :: Qx
    
    FP_Cx = (0.0D+00,0.0D+00)
    
    do i = 1,k_dim
       FP_Cx = FP_Cx + Qx(:,:,:,i)*FP(i)
    enddo
    
    return
  end subroutine matvec

  subroutine mxm_tpls(A,B,nr,nc,sx,ex,sy,ey,sz,ez)

    implicit none

    !----- Inputs -----!

    integer :: nr, nc, sx, ex, sy, ey, sz, ez
    
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz:ez,1:nr), intent(inout) :: A
    double precision, dimension(1:nr,1:nc), intent(in)                         :: B

    !----- Miscellaneous -----!
    
    integer :: i, j, k
    integer :: num_points
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz:ez,1:nc) :: A_temp
    double precision, dimension(:,:), allocatable :: C, C_temp

    num_points = size(A,1)*size(A,2)*size(A,3)

    allocate( C( 1:num_points , size(A,4) ) )
    allocate( C_temp( 1:num_points , nc ) )

    C = 0.0D+00

    C = reshape( A, (/ num_points, size(A,4) /) )

    C_temp = matmul(C,B)

    A_temp = reshape( C_temp, (/ size(A,1), size(A,2), size(A,3), nc /) )
    A(:,:,:,1:nc) = A_temp

    deallocate(C)

    return
  end subroutine mxm_tpls

  subroutine matrix_matrix(A,B,nra,nc,ncb)

    implicit none
    
    integer :: nra, nc, ncb
    
    !     ----- Required variables for dgeem -----
    
    character :: TRANSA = 'N', TRANSB = 'N'
    integer   :: m, n, k
    integer   :: lda, ldb, ldc
    double precision, dimension(nra,nc)  :: A
    double precision, dimension(nc,ncb)  :: B
    double precision, dimension(nra,ncb) :: C
    double precision :: alpha, beta
    
    !     ----- Matrix-matrix multiplication -----
    
    m = nra
    n = ncb
    k = nc
    
    lda = max(1,m)
    ldb = max(1,k)
    ldc = max(1,m)
    
    alpha = 1.0D0
    beta  = 0.0D0
    
    call dgemm(TRANSA,TRANSB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
    
    A(:,1:ncb) = C
    
    return
  end subroutine matrix_matrix
  
  subroutine Schur_decomposition(A,VS,VP,n,sdim)

    implicit none
    
    !     ----- Required variables for dgees -----
    
    character*1           :: JOBVS = 'V', SORT = 'S'
    integer               :: n, LDA, SDIM, LDVS, LWORK, INFO
    double precision, dimension(n,n)  :: A, VS
    double precision, dimension(n)    :: WR, WI
    double precision, dimension(3*n)  :: WORK
    logical, dimension(n) :: BWORK
    double complex, dimension(n) :: VP
    
!    external select_eigenvalue
    
    !     ----- Schur decomposition -----
    
    LDA   = max(1,n)
    LDVS  = max(1,n)
    LWORK = max(1,3*n)
    
    call dgees(JOBVS, SORT, SELECT_EIGENVALUE, N, A, LDA &
         , SDIM, WR, WI, VS, LDVS, WORK, LWORK, BWORK, INFO)
    
    VP = WR*(1.0D+00,0.0D+00) + WI*(0.0D+00,1.0D+00)
    
    return
  end subroutine Schur_decomposition
  
  subroutine compute_eigenvec_schur(A,FP,n)

    implicit none
    
    double complex, dimension(n,n) :: FP
    
    !     ----- Required variables for ztrevc -----
    
    character*1 :: JOB = 'R', HOWMY = 'A'
    integer :: n
    double precision, dimension(n,n) :: A
    double precision, dimension(n,n) :: T, VR, VL
    integer :: ldt, ldvl, ldvr, mm, m
    double precision, dimension(3*n) :: WORK
    integer :: INFO, i
    double precision :: norme
    
    !     ----- Compute the eigenvectors of T -----
    
    T    = A
    ldt  = max(1,n)
    ldvl = max(1,n)
    ldvr = max(1,n)
    m    = n
    mm   = m
    
    call dtrevc( JOB, HOWMY, select_eigenvalue, N, T, LDT, VL, &
         LDVL, VR, LDVR, MM, M, WORK, INFO )            
    
    
    i = 1
    do while ( i.LT.n )
       
       if ( abs(T(i+1,i)).EQ.0.0D+00 ) then
          FP(:,i) = VR(:,i)*(1.0D+00,0.0D+00)
          
          norme   = dsqrt( dot_product(VR(:,i),VR(:,i)) )
          FP(:,i) = FP(:,i)/norme
          
          i = i + 1
       else
          FP(:,i)   = VR(:,i)*(1.0D+00,0.0D+00) + VR(:,i+1)*(0.0D+00,1.0D+00)
          FP(:,i+1) = VR(:,i)*(1.0D+00,0.0D+00) - VR(:,i+1)*(0.0D+00,1.0D+00)
          
          norme     = dsqrt( dot_product(VR(:,i),VR(:,i)) + dot_product(VR(:,i+1),VR(:,i+1))  )
          FP(:,i)   = FP(:,i)/norme
          FP(:,i+1) = FP(:,i+1)/norme
          
          i = i + 2
       endif
       
    enddo
    
    return
  end subroutine compute_eigenvec_schur

  function select_eigenvalue(wr,wi)

    implicit none
    
    logical :: select_eigenvalue
    double precision :: wr, wi
    double precision :: delta = 0.05
    
    if ( sqrt(wr**2. + wi**2.) .GT. 1.d0-delta ) then
       select_eigenvalue = .true.
    else
       select_eigenvalue = .false.
    endif
    
  end function select_eigenvalue
  
  subroutine sort_eigendecomp(VP,FP,residu,n)

    implicit none
    
    integer                        :: n
    double complex, dimension(n)   :: VP
    double complex, dimension(n,n) :: FP
    double precision, dimension(n) :: residu, norm
    double precision               :: temp_real, temp_residu
    double complex                 :: temp_complex
    double complex, dimension(n)   :: temp_n
    integer                        :: i, j, k, l
    
    !     ----- Sorting the eigenvalues according to their norm -----
    
    temp_n   = (0.0D+00,0.0D+00)
    norm     = 0.0D+00
    
    norm = sqrt(dble(VP)**2. + dimag(VP)**2.)
    
    do k = 1,n-1
       do l = k+1,n
          
          if (norm(k).LT.norm(l)) then
             
             temp_real    = norm(k)
             temp_complex = VP(k)
             temp_n       = FP(:,k)
             temp_residu  = residu(k)
             
             norm(k) = norm(l)
             norm(l) = temp_real
             
             residu(k) = residu(l)
             residu(l) = temp_residu
             
             VP(k)  = VP(l)
             VP(l)  = temp_complex
             
             FP(:,k) = FP(:,l)
             FP(:,l) = temp_n
             
          endif
          
       enddo
    enddo
    
    return
  end subroutine sort_eigendecomp
  
  
  subroutine inner_product(ip,phi1,phi2,vx1,vy1,vz1,vx2,vy2,vz2)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Inputs -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), intent(in) :: vx1, vx2, vy1, vy2
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)  , intent(in) :: vz1, vz2
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)  , intent(in) :: phi1, phi2
    double precision :: tempin(2),tempout(2)

    !----- Output -----!

    double precision, intent(out) :: IP

    !----- Miscellaneous -----!

    integer :: i, j, k
    double precision :: temp
    double precision, dimension(sx:ex, sy:ey, 0:maxn-2) :: work1, work2

    !----- Kinetic energy -----!

    ip = 0
    
    work1 = .5* ( vx1(sx:ex, sy:ey, :) + vx1(sx+1:ex+1, sy:ey, :) )
    work2 = .5* ( vx2(sx:ex, sy:ey, :) + vx2(sx+1:ex+1, sy:ey, :) )

    call scalar_product(temp, work1, work2, size(work1))
    ip = ip + temp

    if ( maxm .gt. 2 ) then
       work1 = .5* ( vy1(sx:ex, sy:ey, :) + vy1(sx:ex, sy+1:ey+1, :) )
       work2 = .5* ( vy2(sx:ex, sy:ey, :) + vy2(sx:ex, sy+1:ey+1, :) )
       
       call scalar_product(temp, work1, work2, size(work1))
       ip = ip + temp
    endif

    work1 = .5* ( vz1(sx:ex, sy:ey, 0:maxn-2) + vz1(sx:ex, sy:ey, 1:maxn-1) )
    work2 = .5* ( vz2(sx:ex, sy:ey, 0:maxn-2) + vz2(sx:ex, sy:ey, 1:maxn-1) )

    call scalar_product(temp, work1, work2, size(work1))
    ip = ip + temp

    ip = ip*(dx*dy*dz)
    
!!$    temp = 0.0D+00
!!$
!!$    do k = 0,maxn-2
!!$       do j = sy,ey
!!$          do i = sx,ex
!!$             temp = temp &
!!$                  + (vx1(i,j,k) + vx1(i+1,j,k)) * (vx2(i,j,k) + vx2(i+1,j,k))/4.0D+00 &
!!$                  + (vy1(i,j,k) + vy1(i,j+1,k)) * (vy2(i,j,k) + vy2(i,j+1,k))/4.0D+00 &
!!$                  + (vz1(i,j,k) + vz1(i,j,k+1)) * (vz2(i,j,k) + vz2(i,j,k+1))/4.0D+00
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    if ((my_id .eq. 127) .and. (istep .gt. 99)) then
!!$       write(*,*) 'temporary value after stupid integral at istep =', istep
!!$       write(*,*) 'temp', temp
!!$    endif
!!$    !----- Interfacial energy -----!
!!$
!!$    if ( scap.ne.0 ) then 
!!$       temp = temp + 0.0D+00
!!$    endif
!!$    
!!$    temp = temp*dx*dz*dy
!!$
!!$    if (istep .gt. 99) then
!!$    tempin(1)=temp
!!$    tempin(2)=my_id   
!!$    call mpi_reduce(tempin, tempout, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, master_id, &
!!$                                                               comm2d_quasiperiodic, ierr);
!!$ 
!!$    if (my_id .eq. master_id) write(*,*) 'maximum temp before mpi_allreduce is', tempout(1)
!!$    if (my_id .eq. master_id) write(*,*) 'proc. ', tempout(2)
!!$    endif
!!$ 
!!$    call mpi_allreduce(temp,IP,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d_quasiperiodic,ierr)
!!$
!!$    if ((my_id .eq. 127) .and. (istep .gt. 99)) then
!!$       write(*,*) 'temporary value just after mpi_allreduce at istep =', istep
!!$       write(*,*) 'temp', temp
!!$    endif

    return
  end subroutine inner_product

  subroutine get_difference(ua,ub,ex,ey,ez,sx,sy,sz,diff)
    
    implicit none

    double precision :: diff, sum
    integer :: ex,ey,ez,sx,sy,sz,i,j,k
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz:ez) :: ua, ub
    
    sum = 0.0D+00
    
    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             sum = sum + (ua(i,j,k)-ub(i,j,k))**2.0D+00
          end do
       end do
    end do
    
    diff = sum
    
    return
  end subroutine get_difference

  subroutine norm(norme, x, n)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: x
    double precision :: norme, dummy

    dummy = dot_product(x,x)
    call mpi_allreduce(dummy, norme, 1, mpi_double_precision, mpi_sum, comm2d_quasiperiodic, ierr)
    norme = sqrt(norme)

    return
  end subroutine norm

  subroutine scalar_product(sp, x, y, n)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: x, y
    double precision :: sp, dummy

    dummy = dot_product(x,y)
    call mpi_allreduce(dummy, sp, 1, mpi_double_precision, mpi_sum, comm2d_quasiperiodic, ierr)
    
    return
  end subroutine scalar_product

  double precision function sign_of(a)

    implicit none
    double precision :: a

    if ( a .LT. 0.0D+00 ) then
       sign_of = -1.0D+00
    else
       sign_of = 1.0D+00
    endif

    return
  end function sign_of

  double precision function minmod(a,b)

    implicit none
    double precision :: a, b

    if ( a*b .LT. 0.0D+00 ) then
       minmod = 0.0D+00
    else
       minmod = dmin1(dabs(a),dabs(b))
    endif
    
    return
  end function minmod

  double precision function Dirac_function(val)

    use tpls_constants
    implicit none

    double precision :: val
    double precision :: pi

    pi = 4.D0*atan(1.D0)

    if ( abs(val).GT.smooth_width ) then
       Dirac_function = 0.0D+00
    else
       Dirac_function = (1.d0/(2.d0*smooth_width))+(1.d0/(2.d0*smooth_width))*dcos(pi*val/smooth_width)
!!$       Dirac_function = (1.0D+00/smooth_width) * (1.0D+00 - dabs(val/smooth_width))
!!$       Dirac_function = 1.0D+00/(2.0D+00 * smooth_width)
    endif

    return
  end function Dirac_function

  double precision function Heaviside_function(val)

    use tpls_constants
    implicit none

    double precision :: val
    double precision :: pi
    
    pi = 4.d0*atan(1.d0)

    if ( val.LT.-smooth_width ) then
       Heaviside_function = 0.d0
    elseif ( val.GT.smooth_width) then
       Heaviside_function = 1.d0
    else
       Heaviside_function = .5d0 + (val/(2.d0*smooth_width)) &
            + (1.d0/(2.d0*pi)*sin(pi*val/smooth_width))
    endif
            
    return
  end function Heaviside_function

  double precision function get_volume(phi)

    use tpls_constants
    use tpls_mpi
    use mpi
    
    !----- Input -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p) :: phi

    !----- Compute the volume of phase 1 -----!

    double precision :: volume
    integer :: i, j, k

    volume = 0.0D+00

    do k = 1,maxn-1
       do j = sy,ey
          do i = sx,ex

             volume = volume + Heaviside_function(phi(i,j,k))

          enddo
       enddo
    enddo

    call mpi_allreduce(volume,get_volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d_quasiperiodic,ierr)
    get_volume = get_volume*dx*dy*dz

    return
  end function get_volume

  double precision function get_cfl(vx,vy,vz)

    use tpls_constants
    use tpls_mpi
    use mpi

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) :: vx, vy
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)   :: vz
    double precision :: dummy, maxu
    
    dummy = max(maxval(vx),maxval(vy),maxval(vz))
    call mpi_allreduce(dummy,maxu,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d_quasiperiodic,ierr)
    get_cfl = maxu*dt/dx

    return
  end function get_cfl
  
end module tpls_maths
