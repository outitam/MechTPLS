module tpls_constants

  implicit none

  !----- Definition of the non-dimensional numbers -----!

  double precision :: Re
  double precision :: mu_plus, mu_minus, mu_0
  double precision :: rho_plus, rho_minus, rho_0
  double precision :: Scap
  double precision :: height
  double precision :: dpdl
  parameter ( dpdl = -1.d0 )
  double precision :: gz, Grav

  !----- Mesh related general variables -----!

  integer :: maxl, maxm, maxn
  double precision :: dx, dy, dz, dt
  double precision :: Lx, Ly, Lz
  
  !----- Miscellaneous -----!

  integer :: nsteps
  integer :: istep
  integer :: iostep

  !----- MPI related stuff -----!

  integer :: num_procs_x, num_procs_y, num_procs_z

  !----- Levelset related general variables -----!
  
  double precision :: smooth_width
  double precision :: tolerance_levelset
  integer          :: max_iteration_levelset

  !----- Setup of the Poisson and Helmholtz solvers -----!

  double precision :: tolerance_poisson
  integer          :: max_iteration_poisson
  double precision :: tolerance_helmholtz
  integer          :: max_iteration_helmholtz

  !----- Scheduled Relaxation Jacobi -----!

  integer :: max_srj
  double precision, dimension(:), allocatable :: relax_srj

  !----- Logical parameters -----!

  logical :: density_matched
  logical :: viscosity_matched
  logical :: if_implicit
  logical :: if_sfd
  logical :: if_restart
  logical :: if_exact_restart
  logical :: if_exact_restart_sfd
  logical :: if_exact_restart_save
  logical :: if_record_perturbation_norm
  logical :: if_linear_stability

  !----- Filenames for restart -----!

  character*80 :: restart_handle

  !----- Linear stability setup -----!

  integer :: krylov_dim
  parameter ( krylov_dim = 128 )

  !----- Runtime statistics -----!

  double precision :: total_time, time_levelset, time_fluid, time_pressure, t_temp, t_temp2

end module tpls_constants
