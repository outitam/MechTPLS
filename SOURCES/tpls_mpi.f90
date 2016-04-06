module tpls_mpi
  
  implicit none

  integer :: Ndim, etag, stag
  parameter ( Ndim = 3 )
  parameter ( etag = 1 , stag = 2 )

  integer :: dims(Ndim), coords(Ndim)
  logical :: periodic(Ndim), reorder

  integer :: master_id
  parameter ( master_id = 0 )

  integer :: my_id, sender_id, num_procs
  integer :: comm2d_quasiperiodic
  integer :: neighbours(6)

  integer :: ierr

  integer :: n_local_x, n_local_y, n_local_z_uv, n_local_z_w, n_local_z_p
  integer :: sx, sy, sz_uv, sz_w, sz_p, sx_max
  integer :: ex, ey, ez_uv, ez_w, ez_p, ex_max
  integer :: sx_o, ex_o
  integer :: sy_o, ey_o
  integer :: sx_f, ex_f
  integer :: sy_f, ey_f

  integer :: px, py, pz

  integer :: stride_uv_xz, stride_uv_yz
  integer :: stride_uv_aug1_xz, stride_uv_aug1_yz
  integer :: stride_uv_aug2_xz, stride_uv_aug2_yz

  integer :: stride_w_xz, stride_w_yz
  integer :: stride_w_aug1_xz, stride_w_aug1_yz
  integer :: stride_w_aug2_xz, stride_w_aug2_yz

  integer :: stride_p_xz, stride_p_yz
  integer :: stride_p_aug1_xz, stride_p_aug1_yz
  integer :: stride_p_aug2_xz, stride_p_aug2_yz
  integer :: stride_p_augaug1_xz, stride_p_augaug1_yz
  integer :: stride_p_augaug2_xz, stride_p_augaug2_yz
  integer :: stride_p_augaug3_xz, stride_p_augaug3_yz

contains

  subroutine initialize_tpls_mpi

    use tpls_constants
    use mpi
    implicit none
    
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,my_id,ierr)
    call mpi_comm_size(mpi_comm_world,num_procs,ierr)

    dims(1) = num_procs_x
    dims(2) = num_procs_y
    dims(3) = num_procs_z
    
    periodic(1) = .false.
    periodic(2) = .true.
    periodic(3) = .false.
    reorder     = .false.
    
    call mpi_cart_create(mpi_comm_world,Ndim,dims,periodic,reorder,comm2d_quasiperiodic,ierr)
    call mpi_comm_rank(  comm2d_quasiperiodic,my_id,ierr)
    call mpi_cart_coords(comm2d_quasiperiodic,my_id,Ndim,coords,ierr)
    call get_mpi_neighbours(neighbours,comm2d_quasiperiodic)
    
    call mpi_barrier(mpi_comm_world,ierr)
    
    call mpi_decomp_2d(sx,ex,sy,ey,n_local_x,n_local_y,maxl,maxm,coords,dims,Ndim)
    
    sz_uv        = 0
    ez_uv        = maxn-2
    n_local_z_uv = ez_uv-sz_uv + 1
    
    sz_w        = 0
    ez_w        = maxn-1
    n_local_z_w = ez_w - sz_w + 1 
    
    sz_p        = 0
    ez_p        = maxn
    n_local_z_p = ez_p - sz_p + 1 
    
    call mpi_barrier(mpi_comm_world,ierr)
    
    ! Pick out max value ex -- signal for a right-hand boundary.
    Call mpi_allreduce(ex,ex_max,1,mpi_integer,mpi_max,mpi_comm_world,ierr)

    ! First set of strides - unchanged
    
    call get_stride_p(stride_uv_yz,stride_uv_xz,sx,ex,sy,ey,sz_uv,ez_uv)
    call get_stride_p(stride_w_yz, stride_w_xz, sx,ex,sy,ey,sz_w, ez_w )
    call get_stride_p(stride_p_yz, stride_p_xz, sx,ex,sy,ey,sz_p, ez_p )
    
    ! Second set of strides: first- and second-order halos for augmented arrays
    
    call get_stride_p_aug1(stride_uv_aug1_yz,stride_uv_aug1_xz,sx,ex,sy,ey,sz_uv,ez_uv)
    call get_stride_p_aug2(stride_uv_aug2_yz,stride_uv_aug2_xz,sx,ex,sy,ey,sz_uv,ez_uv)
    
    call get_stride_p_aug1(stride_w_aug1_yz,stride_w_aug1_xz,sx,ex,sy,ey,sz_w,ez_w)
    call get_stride_p_aug2(stride_w_aug2_yz,stride_w_aug2_xz,sx,ex,sy,ey,sz_w,ez_w)
    
    call get_stride_p_aug1(stride_p_aug1_yz,stride_p_aug1_xz,sx,ex,sy,ey,sz_p,ez_p)
    call get_stride_p_aug2(stride_p_aug2_yz,stride_p_aug2_xz,sx,ex,sy,ey,sz_p,ez_p)
    
    ! Third set of strides: first-, second-, and third-order halos for augmented arrays
    
    call get_stride_p_augaug1(stride_p_augaug1_yz,stride_p_augaug1_xz,sx,ex,sy,ey,sz_p,ez_p)
    call get_stride_p_augaug2(stride_p_augaug2_yz,stride_p_augaug2_xz,sx,ex,sy,ey,sz_p,ez_p)
    call get_stride_p_augaug3(stride_p_augaug3_yz,stride_p_augaug3_xz,sx,ex,sy,ey,sz_p,ez_p)

    ! ****************************************************************************************  
    ! Getting coords of base node in Cartesian topology.
    
    if ( my_id == master_id ) then
       
       pz = 0
       do py = 0,dims(2)-1
          do px = 0,dims(1)-1
             
             if ( (pz == 0) .and. (py == 0) .and. (px ==0 ) ) then
                sx_o = sx
                ex_o = ex
                sy_o = sy
                ey_o = ey 
             end if
             
          end do
       end do
       
    end if
    
    call mpi_barrier(mpi_comm_world,ierr)
    
    return
  end subroutine initialize_tpls_mpi

  subroutine mpi_decomp_2d(sx,ex,sy,ey,n_local_x,n_local_y,   &
       maxl,maxm,coords,dims,Ndim)     

    use mpi
    implicit none
    
    integer :: sx,ex,sy,ey,sz,ez,n_local_x,n_local_y,n_local_z
    integer :: n_procs_x,n_procs_y,n_procs_z,rank
    integer :: maxl,maxm,Ndim,ierr,my_id
    integer :: coords(Ndim),dims(Ndim)
    
    integer :: nt,remainder
    
    n_procs_x=dims(1)        
    rank=coords(1)
    n_local_x = (maxl-1)/n_procs_x
    sx = rank*n_local_x + 1
    ex = sx + n_local_x - 1
    
    n_procs_y=dims(2)        
    rank=coords(2)
    n_local_y = (maxm-1)/n_procs_y
    sy = rank*n_local_y + 1
    ey = sy + n_local_y - 1
    
    call mpi_comm_rank(mpi_comm_world,my_id,ierr)
    Write(*,*) 'Rank:',my_id,'Coordinates:',sx,ex,sy,ey
    
    return
  end subroutine mpi_decomp_2d
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  
  subroutine get_mpi_neighbours(neighbours,comm3d)     

    use mpi
    implicit none
    
    integer :: neighbours(6), N,E,S,W,Fr,Bk,comm3d,ierr
    parameter (N=1,E=2,S=3,W=4,Fr=5,Bk=6)
    
    neighbours(1:6)  =  MPI_PROC_NULL
    CALL  MPI_CART_SHIFT(comm3d,  0, 1, neighbours(W),  neighbours(E),  ierr)
    CALL  MPI_CART_SHIFT(comm3d,  1, 1, neighbours(S),  neighbours(N),  ierr)
    CALL  MPI_CART_SHIFT(comm3d,  2, 1, neighbours(Fr), neighbours(Bk), ierr)
    
    return
  end subroutine get_mpi_neighbours
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! For passing 1st-order halos of non-augmented arrays WITH CORNERS (slightly modified code)
  
  subroutine get_stride_p(stride_p_yz,stride_p_xz,sx,ex,sy,ey,sz,ez)

    use mpi
    implicit none
    
    integer :: stride_p_xz,stride_p_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real
    
    ! For passing xz planes in the y direction (North-South).  This one is changed from before to allow
    ! the passing of vertex points.
    CALL  MPI_TYPE_VECTOR (    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+3,              & ! longueur d'un bloc (change this one for bigger blocks)
         (ex-sx+3)*(ey-sy+3),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_xz, ierr  )
    CALL  MPI_TYPE_COMMIT (  stride_p_xz,  ierr  )
    
    ! For passing yz planes in the x direction (East-West).  This one is ALSO changed from before.
    CALL MPI_TYPE_SIZE ( mpi_double_precision,  size_real, ierr )
    
    CALL MPI_TYPE_HVECTOR (    &
         ey-sy+3,              &
         1,                    &
         (ex-sx+3)*size_real,  &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )
    
    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+3)*(ex-sx+3)*size_real, & 
         type_y,    &           
         stride_p_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_yz,  ierr  )
    
    return
  end subroutine get_stride_p
  
  ! ****************************************************************************************
  ! For passing 1st-order halos of non-augmented arrays WITH CORNERS (slightly modified code)
  
  Subroutine exchange2d(uu,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none
    
    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_xz,stride_p_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz:ez), intent(inout) :: uu
    
    integer :: tag1=100,tag2=200,ierr,status(mpi_status_size)
    integer :: N,E,S,W,Fr,Bk
    parameter (N=1,E=2,S=3,W=4,Fr=5,Bk=6)
    
    
    ! ******************************************************************************
    ! Send neighbours "N" and receive neighbours "S" using the stride type stride_p_xz.  
    ! This part is changed to allow for sending extra strips of data.
    
    CALL  MPI_SENDRECV( &
         uu(sx-1,ey,  sz),  1,  stride_p_xz,  neighbours(N),  tag1, &
         uu(sx-1,sy-1,sz),  1,  stride_p_xz,  neighbours(S),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "S" and receive neighbours "N" using the stride type stride_p_xz
    CALL  MPI_SENDRECV( &
         uu(sx-1,sy,  sz),  1,  stride_p_xz,  neighbours(S),  tag2, &
         uu(sx-1,ey+1,sz),  1,  stride_p_xz,  neighbours(N),  tag2, &
         comm_topology,  status,  ierr )
    
    ! ******************************************************************************
    ! This part is also changed to allow for sending extra strips of data.  
    ! Corners need to be exchanged twice becuase in the first (N-S) swap, 
    ! wrong information is received into corners.  
    ! The second (E-W) swap fixes this problem.
    
    ! send neighbours "W" and receive neighbours "E" using the stride type stride_p_yz
    CALL  MPI_SENDRECV( &
         uu(sx,  sy-1,sz),  1,  stride_p_yz,  neighbours(W),  tag1, &
         uu(ex+1,sy-1,sz),  1,  stride_p_yz,  neighbours(E),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "E" and receive neighbours "W" using the stride type stride_p_yz
    CALL  MPI_SENDRECV( &
         uu(ex,  sy-1,sz),  1,  stride_p_yz,  neighbours(E),  tag2, &
         uu(sx-1,sy-1,sz),  1,  stride_p_yz,  neighbours(W),  tag2, &
         comm_topology,  status,  ierr )
    
    return
  end subroutine exchange2d
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! For passing 1st-order halos of AUGMENTED arrays WITH CORNERS (modified code)
  ! ****************************************************************************************
  ! ****************************************************************************************
  
  subroutine get_stride_p_aug1(stride_p_aug1_yz,stride_p_aug1_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none
    
    integer :: stride_p_aug1_xz,stride_p_aug1_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real
    
    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+3,              & ! longueur d'un bloc
         (ex-sx+5)*(ey-sy+5),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_aug1_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_aug1_xz,  ierr  )        
    
    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+3, &
         1, &
         (ex-sx+5)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )
    
    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+5)*(ex-sx+5)*size_real, & 
         type_y,    &           
         stride_p_aug1_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_aug1_yz,  ierr  )
    
    return
  end subroutine get_stride_p_aug1
  
  ! ****************************************************************************************
  ! For passing 1st-order halos of augmented arrays WITH CORNERS (modified code)
  
  Subroutine exchange2d_aug1(uu,stride_p_aug1_xz,stride_p_aug1_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none
    
    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_aug1_xz,stride_p_aug1_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-2:ex+2,sy-2:ey+2,sz:ez), intent(inout) :: uu
    
    integer :: tag1=100,tag2=200,ierr,status(mpi_status_size)
    integer :: N,E,S,W,Fr,Bk
    parameter (N=1,E=2,S=3,W=4,Fr=5,Bk=6)
    
    ! ******************************************************************************
    
    ! Send neighbours "N" and receive neighbours "S" using the stride type stride_p_xz.  
    ! This one is changed to allow for sending extra strips of data.
    CALL  MPI_SENDRECV( &
         uu(sx-1,ey,  sz),  1,  stride_p_aug1_xz,  neighbours(N),  tag1, &
         uu(sx-1,sy-1,sz),  1,  stride_p_aug1_xz,  neighbours(S),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "S" and receive neighbours "N" using the stride type stride_p_aug1_xz
    CALL  MPI_SENDRECV( &
         uu(sx-1,sy,  sz),  1,  stride_p_aug1_xz,  neighbours(S),  tag2, &
         uu(sx-1,ey+1,sz),  1,  stride_p_aug1_xz,  neighbours(N),  tag2, &
         comm_topology,  status,  ierr )
    
    ! ******************************************************************************
    
    ! send neighbours "W" and receive neighbours "E" using the stride type stride_p_yz
    CALL  MPI_SENDRECV( &
         uu(sx,  sy-1,sz),  1,  stride_p_aug1_yz,  neighbours(W),  tag1, &
         uu(ex+1,sy-1,sz),  1,  stride_p_aug1_yz,  neighbours(E),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "E" and receive neighbours "W" using the stride type stride_p_aug1_yz
    CALL  MPI_SENDRECV( &
         uu(ex,  sy-1,sz),  1,  stride_p_aug1_yz,  neighbours(E),  tag2, &
         uu(sx-1,sy-1,sz),  1,  stride_p_aug1_yz,  neighbours(W),  tag2, &
         comm_topology,  status,  ierr )
    
    return
  end subroutine exchange2d_aug1
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! For passing 2nd-order halos of augmented arrays WITHOUT CORNERS (modified code)
  
  subroutine get_stride_p_aug2(stride_p_aug2_yz,stride_p_aug2_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none
    
    integer :: stride_p_aug2_xz,stride_p_aug2_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real
    
    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+3,              & ! longueur d'un bloc
         (ex-sx+5)*(ey-sy+5),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_aug2_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_aug2_xz,  ierr  )        
    
    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+3, &
         1, &
         (ex-sx+5)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )
    
    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+5)*(ex-sx+5)*size_real, & 
         type_y,    &           
         stride_p_aug2_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_aug2_yz,  ierr  )
    
    return
  end subroutine get_stride_p_aug2
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! For passing 2nd-order halos of augmented arrays WITHOUT CORNERS (modified code)
  
  Subroutine exchange2d_aug2(uu,stride_p_aug2_xz,stride_p_aug2_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none
    
    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_aug2_xz,stride_p_aug2_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-2:ex+2,sy-2:ey+2,sz:ez), intent(inout) :: uu
    
    integer :: tag1=100,tag2=200,ierr,status(mpi_status_size)
    integer :: N,E,S,W,Fr,Bk
    parameter (N=1,E=2,S=3,W=4,Fr=5,Bk=6)
    
    ! ******************************************************************************
    
    ! send neighbours "N" and receive neighbours "S" using the stride type stride_p_aug2_xz
    CALL  MPI_SENDRECV( &
         uu(sx-1, ey-1,  sz),  1,  stride_p_aug2_xz,  neighbours(N),  tag1, &
         uu(sx-1, sy-2,sz),  1,  stride_p_aug2_xz,  neighbours(S),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "S" and receive neighbours "N" using the stride type stride_p_aug2_xz
    CALL  MPI_SENDRECV( &
         uu(sx-1, sy+1,  sz),  1,  stride_p_aug2_xz,  neighbours(S),  tag2, &
         uu(sx-1, ey+2,sz),  1,  stride_p_aug2_xz,  neighbours(N),  tag2, &
         comm_topology,  status,  ierr )
    
    ! ******************************************************************************
    
    ! send neighbours "W" and receive neighbours "E" using the stride type stride_p_aug2_yz
    CALL  MPI_SENDRECV( &
         uu(sx+1, sy-1,sz),  1,  stride_p_aug2_yz,  neighbours(W),  tag1, &
         uu(ex+2, sy-1,sz),  1,  stride_p_aug2_yz,  neighbours(E),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "E" and receive neighbours "W" using the stride type stride_p_aug2_yz
    CALL  MPI_SENDRECV( &
         uu(ex-1, sy-1,sz),  1,  stride_p_aug2_yz,  neighbours(E),  tag2, &
         uu(sx-2, sy-1,sz),  1,  stride_p_aug2_yz,  neighbours(W),  tag2, &
         comm_topology,  status,  ierr )
    
    return
  end subroutine exchange2d_aug2
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! Doubly-Augmented arrays
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! For passing 1st-order halos of DOUBLY-augmented arrays WITH CORNERS
  
  subroutine get_stride_p_augaug1(stride_p_augaug1_yz,stride_p_augaug1_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none
    
    integer :: stride_p_augaug1_xz,stride_p_augaug1_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real
    
    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+3,              & ! longueur d'un bloc
         (ex-sx+7)*(ey-sy+7),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_augaug1_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_augaug1_xz,  ierr  )        
    
    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+3, &
         1, &
         (ex-sx+7)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )
    
    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+7)*(ex-sx+7)*size_real, & 
         type_y,    &           
         stride_p_augaug1_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_augaug1_yz,  ierr  )
    
    return
  end subroutine get_stride_p_augaug1
  
  ! ****************************************************************************************
  ! For passing 1st-order halos of DOUBLY-augmented arrays WITH CORNERS
  
  Subroutine exchange2d_augaug1(uu,stride_p_augaug1_xz,stride_p_augaug1_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none
    
    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_augaug1_xz,stride_p_augaug1_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-3:ex+3,sy-3:ey+3,sz:ez), intent(inout) :: uu
    
    integer :: tag1=100,tag2=200,ierr,status(mpi_status_size)
    integer :: N,E,S,W,Fr,Bk
    parameter (N=1,E=2,S=3,W=4,Fr=5,Bk=6)
    
    ! ******************************************************************************
    
    ! Send neighbours "N" and receive neighbours "S" using the stride type stride_p_xz.  
    ! This one is changed to allow for sending extra strips of data.
    CALL  MPI_SENDRECV( &
         uu(sx-1,ey,  sz),  1,  stride_p_augaug1_xz,  neighbours(N),  tag1, &
         uu(sx-1,sy-1,sz),  1,  stride_p_augaug1_xz,  neighbours(S),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "S" and receive neighbours "N" using the stride type stride_p_augaug1_xz
    CALL  MPI_SENDRECV( &
         uu(sx-1,sy,  sz),  1,  stride_p_augaug1_xz,  neighbours(S),  tag2, &
         uu(sx-1,ey+1,sz),  1,  stride_p_augaug1_xz,  neighbours(N),  tag2, &
         comm_topology,  status,  ierr )
    
    ! ******************************************************************************
    
    ! send neighbours "W" and receive neighbours "E" using the stride type stride_p_yz
    CALL  MPI_SENDRECV( &
         uu(sx,  sy-1,sz),  1,  stride_p_augaug1_yz,  neighbours(W),  tag1, &
         uu(ex+1,sy-1,sz),  1,  stride_p_augaug1_yz,  neighbours(E),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "E" and receive neighbours "W" using the stride type stride_p_augaug1_yz
    CALL  MPI_SENDRECV( &
         uu(ex,  sy-1,sz),  1,  stride_p_augaug1_yz,  neighbours(E),  tag2, &
         uu(sx-1,sy-1,sz),  1,  stride_p_augaug1_yz,  neighbours(W),  tag2, &
         comm_topology,  status,  ierr )
    
    return
  end subroutine exchange2d_augaug1
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! For passing 2nd-order halos of DOUBLY-augmented arrays WITH CORNERS
  
  subroutine get_stride_p_augaug2(stride_p_augaug2_yz,stride_p_augaug2_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none
    
    integer :: stride_p_augaug2_xz,stride_p_augaug2_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real
    
    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+5,              & ! longueur d'un bloc
         (ex-sx+7)*(ey-sy+7),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_augaug2_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_augaug2_xz,  ierr  )        
    
    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+5, &
         1, &
         (ex-sx+7)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )
    
    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+7)*(ex-sx+7)*size_real, & 
         type_y,    &           
         stride_p_augaug2_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_augaug2_yz,  ierr  )
    
    return
  end subroutine get_stride_p_augaug2
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! For passing 2nd-order halos of DOUBLY-augmented arrays WITH CORNERS
  
  Subroutine exchange2d_augaug2(uu,stride_p_augaug2_xz,stride_p_augaug2_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none
    
    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_augaug2_xz,stride_p_augaug2_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-3:ex+3,sy-3:ey+3,sz:ez), intent(inout) :: uu
    
    integer :: tag1=100,tag2=200,ierr,status(mpi_status_size)
    integer :: N,E,S,W,Fr,Bk
    parameter (N=1,E=2,S=3,W=4,Fr=5,Bk=6)
    
    ! ******************************************************************************
    
    ! send neighbours "N" and receive neighbours "S" using the stride type stride_p_augaug2_xz
    CALL  MPI_SENDRECV( &
         uu(sx-2, ey-1,  sz),  1,  stride_p_augaug2_xz,  neighbours(N),  tag1, &
         uu(sx-2, sy-2,sz),  1,  stride_p_augaug2_xz,  neighbours(S),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "S" and receive neighbours "N" using the stride type stride_p_augaug2_xz
    CALL  MPI_SENDRECV( &
         uu(sx-2, sy+1,  sz),  1,  stride_p_augaug2_xz,  neighbours(S),  tag2, &
         uu(sx-2, ey+2,sz),  1,  stride_p_augaug2_xz,  neighbours(N),  tag2, &
         comm_topology,  status,  ierr )
    
    ! ******************************************************************************
    
    ! send neighbours "W" and receive neighbours "E" using the stride type stride_p_augaug2_yz
    CALL  MPI_SENDRECV( &
         uu(sx+1, sy-2,sz),  1,  stride_p_augaug2_yz,  neighbours(W),  tag1, &
         uu(ex+2, sy-2,sz),  1,  stride_p_augaug2_yz,  neighbours(E),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "E" and receive neighbours "W" using the stride type stride_p_augaug2_yz
    CALL  MPI_SENDRECV( &
         uu(ex-1, sy-2,sz),  1,  stride_p_augaug2_yz,  neighbours(E),  tag2, &
         uu(sx-2, sy-2,sz),  1,  stride_p_augaug2_yz,  neighbours(W),  tag2, &
         comm_topology,  status,  ierr )
    
    return
  end subroutine exchange2d_augaug2
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! For passing 3rd-order halos of DOUBLY-augmented arrays WITHOUT CORNERS
  
  subroutine get_stride_p_augaug3(stride_p_augaug3_yz,stride_p_augaug3_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none
    
    integer :: stride_p_augaug3_xz,stride_p_augaug3_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real
    
    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+5,              & ! longueur d'un bloc
         (ex-sx+7)*(ey-sy+7),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_augaug3_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_augaug3_xz,  ierr  )        
    
    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+5, &
         1, &
         (ex-sx+7)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )
    
    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+7)*(ex-sx+7)*size_real, & 
         type_y,    &           
         stride_p_augaug3_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_augaug3_yz,  ierr  )
    
    return
  end subroutine get_stride_p_augaug3
  
  ! ****************************************************************************************
  ! ****************************************************************************************
  ! For passing 2nd-order halos of DOUBLY-augmented arrays WITHOUT CORNERS
  
  Subroutine exchange2d_augaug3(uu,stride_p_augaug3_xz,stride_p_augaug3_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none
    
    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_augaug3_xz,stride_p_augaug3_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-3:ex+3,sy-3:ey+3,sz:ez), intent(inout) :: uu
    
    integer :: tag1=100,tag2=200,ierr,status(mpi_status_size)
    integer :: N,E,S,W,Fr,Bk
    parameter (N=1,E=2,S=3,W=4,Fr=5,Bk=6)
    
    ! ******************************************************************************
    
    ! send neighbours "N" and receive neighbours "S" using the stride type stride_p_augaug2_xz
    CALL  MPI_SENDRECV( &
         uu(sx-2, ey-2,  sz),  1,  stride_p_augaug3_xz,  neighbours(N),  tag1, &
         uu(sx-2, sy-3,sz),  1,  stride_p_augaug3_xz,  neighbours(S),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "S" and receive neighbours "N" using the stride type stride_p_augaug2_xz
    CALL  MPI_SENDRECV( &
         uu(sx-2, sy+2,  sz),  1,  stride_p_augaug3_xz,  neighbours(S),  tag2, &
         uu(sx-2, ey+3,sz),  1,  stride_p_augaug3_xz,  neighbours(N),  tag2, &
         comm_topology,  status,  ierr )
    
    ! ******************************************************************************
  
    ! send neighbours "W" and receive neighbours "E" using the stride type stride_p_augaug2_yz
    CALL  MPI_SENDRECV( &
         uu(sx+2, sy-2,sz),  1,  stride_p_augaug3_yz,  neighbours(W),  tag1, &
         uu(ex+3, sy-2,sz),  1,  stride_p_augaug3_yz,  neighbours(E),  tag1, &
         comm_topology,  status,  ierr )
    
    ! send neighbours "E" and receive neighbours "W" using the stride type stride_p_augaug2_yz
    CALL  MPI_SENDRECV( &
         uu(ex-2, sy-2,sz),  1,  stride_p_augaug3_yz,  neighbours(E),  tag2, &
         uu(sx-3, sy-2,sz),  1,  stride_p_augaug3_yz,  neighbours(W),  tag2, &
         comm_topology,  status,  ierr )
    
    return
  end subroutine exchange2d_augaug3
  
end module tpls_mpi
