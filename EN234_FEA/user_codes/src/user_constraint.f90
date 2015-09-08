subroutine user_constraint(constraint, constraint_flag, n_nodes, node_property_list, dof_list,&              ! Input variables
     n_properties, constraint_properties,nodal_coords,length_coord_array, &
     dof_increment, dof_total, length_dof_array, lagrange_multiplier,  &      ! Input variables
     constraint_stiffness,constraint_residual)      ! Output variables
  use Types
  use ParamIO
  use Mesh, only : node
  implicit none

  integer, intent( in )         :: constraint                                             ! Constraint number
  integer, intent( in )         :: constraint_flag                                        ! Flag identifying constraint type
  integer, intent( in )         :: n_nodes                                                ! # nodes in the constraint  
  integer, intent( in )         :: n_properties                                           ! # properties for the element
  integer, intent( in )         :: length_coord_array                                     ! Total # coords for element
  integer, intent( in )         :: length_dof_array                                       ! Total # Dof for the element
  integer, intent( in )         :: dof_list(n_nodes)                                      ! List of DOFs being constrained

  type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
!  type node
!      sequence
!      integer :: flag                          ! Integer identifier
!      integer :: coord_index                   ! Index of first coordinate in coordinate array
!      integer :: n_coords                      ! Total no. coordinates for the node
!      integer :: dof_index                     ! Index of first DOF in dof array
!      integer :: n_dof                         ! Total no. of DOF for node
!   end type node 
!   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

  real( prec ), intent( in )    :: nodal_coords(length_coord_array)                       ! Coordinates, stored as x1,x2,(x3) for each node in constraint in turn
  real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
  real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment
  real( prec ), intent( in )    :: lagrange_multiplier                                    ! Current value of lagrange multiplier

  real( prec ), intent( in )    :: constraint_properties(n_properties)                    ! Constraint properties, stored in order listed in input file
                  
  real( prec ), intent( out )   :: constraint_stiffness(length_dof_array,length_dof_array)   ! Constraint stiffness (ROW,COLUMN)
  real( prec ), intent( out )   :: constraint_residual(length_dof_array)                     ! Constraint residual force (ROW)
  
  integer :: dof,ndof
  
   ! Subroutine to set up stiffness for multi--point constraint.  The constraint has one parameter, specifying a (small) compliance.
  constraint_residual = 0.d0
  constraint_stiffness = 0.d0
!
  ndof = node_property_list(1)%n_dof
  dof = dof_list(1)
  constraint_residual(1) = -lagrange_multiplier
  constraint_residual(2) = lagrange_multiplier
  constraint_residual(3) = dof_total(ndof+dof)+dof_increment(ndof+dof)-dof_total(dof)-dof_increment(dof) &
                          -lagrange_multiplier*constraint_properties(1)
  constraint_stiffness(1, 3) = 1.D0
  constraint_stiffness(3, 1) = 1.D0
  constraint_stiffness(2, 3) = -1.D0
  constraint_stiffness(3, 2) = -1.D0
  constraint_stiffness(3,3) = constraint_properties(1)
  
  
  
  end subroutine user_constraint
