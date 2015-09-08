  module Stiffness
  use Types
  use Linkedlist_Handling

  integer, save :: neq                           ! Total no. equations
  logical, save :: unsymmetric_stiffness         ! Set to .true. if stiffness is unsymmetric
    
  real (prec), save, allocatable :: alow(:)      ! Lower triangle of stiffness
  real (prec), save, allocatable :: aupp(:)      ! Upper triangle of stiffness
  real (prec), save, allocatable :: diag(:)      ! Diagonal
  real (prec), save, allocatable :: rhs(:)       ! RHS of equation system; replaced by solution

  integer, save, allocatable :: ieqs(:)          ! Equation number array
  
  !     Data structures for direct solver
  !
  !     Detailed description of storage:
  !     For array
  !     |  K11  K12   K13   0   0  |
  !     |  K21  K22   K23  K24  0  |
  !     |  K31  K32   K33  K34 K35 |
  !     |  0    K42   K43  K44 K45 |
  !     |  0    0     K53  K54 K55 |

  !     AUPP =  [K12   K13 K23   K24 K34   K35 K45]
  !     ALOW = [K21   K31 K32   K42 K43   K53 K54]
  !     DIAG = [K11 K22 K33 K44 K55]
  !     JPOIN = [0  1  3  5  7]

  integer, save, allocatable :: jpoin(:)        ! Pointer to last element in each row/column of AUPP/ALOW respectively
  integer, save, allocatable :: node_numbers(:)  ! Bandwidth minimizing node number for each node
  integer, save, allocatable :: node_order_index(:) ! Index listing nodes in order of ascending node number

  ! CG solver
  real( prec ), parameter :: cg_tolerance=1.0D-17  ! CG solver tolerance
  real( prec ), parameter :: thresh=0.d0           ! Threshold below which stiffness is assumed to be zero
  integer, parameter :: max_cg_iterations = 1000
  integer, save, allocatable :: rows(:)             ! rows(i) gives row number of ith entry in aupp
  integer, save, allocatable :: cols(:)             ! cols(i) gives column number of ith entry in aupp
  
  integer, save :: last_filled_aupp
  integer, save :: length_aupp
  integer, save :: last_filled_equation_adjacency
  integer, save :: length_equation_adjacency
  
  integer, save, allocatable :: equation_adjacency_index(:)  !  Points to first data block 
  type (Integer_Linked_List), save, allocatable :: equation_adjacency(:)  ! Linked data block list of index in aupp 
  
  real (prec), save, allocatable :: p(:)            ! Workspace arrays for CG solver
  real (prec), save, allocatable :: r(:)
  real (prec), save, allocatable :: z(:)
  real (prec), save, allocatable :: precon(:)
  real (prec), save, allocatable :: sol(:)
  
  

  end module
