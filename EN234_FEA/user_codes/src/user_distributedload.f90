subroutine user_distributed_load_static(lmn, element_identifier,face,dload_parameters,n_parameters,  &  ! Input variables
     n_nodes, node_property_list, &                                                                     ! Input variables
     n_properties, element_properties,element_coords, length_coord_array, &                             ! Input variables
     dof_increment, dof_total, length_dof_array,  &                                                     ! Input variables
     element_stiffness,element_residual)                           ! Output variables
  use Types
  use ParamIO
  use Mesh, only : node
  use Element_Utilities, only : N => shape_functions_2D
  use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
  use Element_Utilities, only : dNdx => shape_function_spatial_derivatives_2D
  use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
  use Element_Utilities, only : dxdxi => jacobian_2D
  use Element_Utilities, only : initialize_integration_points
  use Element_Utilities, only : calculate_shapefunctions
  use Element_Utilities, only : invert_small
  implicit none

  integer, intent( in )         :: lmn                                                    ! Element number
  integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
  integer, intent( in )         :: face                                                   ! Loaded face
  integer, intent( in )         :: n_parameters                                           ! User supplied parameters
  integer, intent( in )         :: n_nodes                                                ! # nodes on the element  
  integer, intent( in )         :: n_properties                                           ! # properties for the element                
  integer, intent( in )         :: length_coord_array                                     ! Total # coords
  integer, intent( in )         :: length_dof_array                                       ! Total # DOF

  real (prec), intent( in )     :: dload_parameters(n_parameters)                         ! User supplied parameters.

  type (node), intent( in )     :: node_property_list(n_nodes)                            ! Data structure describing storage for nodal variables - see below
!  type node
!      sequence
!      integer :: flag                          ! Integer identifier
!      integer :: coord_index                   ! Index of first coordinate in coordinate array
!      integer :: n_coords                      ! Total no. coordinates for the node
!      integer :: dof_index                     ! Index of first DOF in dof array
!      integer :: n_dof                         ! Total no. of DOF for node
!   end type node 
!   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

  real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
  real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
  real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

  real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
                
  real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
  real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
  element_stiffness = 0.d0
  element_residual = 0.d0

  end subroutine user_distributed_load_static
  
  subroutine user_distributed_load_dynamic(lmn, element_identifier,face,dload_parameters,n_parameters,  &  ! Input variables
     n_nodes, node_property_list, &                                                                        ! Input variables
     n_properties, element_properties,element_coords,length_coord_array, &                                 ! Input variables
     dof_increment, dof_total, length_dof_array, &                                                         ! Input variables
     element_residual)                                                                                     ! Output variables
  use Types
  use ParamIO
  use Mesh, only : node
  use Element_Utilities, only : N => shape_functions_2D
  use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
  use Element_Utilities, only : dNdx => shape_function_spatial_derivatives_2D
  use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_2D
  use Element_Utilities, only : dxdxi => jacobian_2D
  use Element_Utilities, only : initialize_integration_points
  use Element_Utilities, only : calculate_shapefunctions
  use Element_Utilities, only : invert_small
  implicit none

  integer, intent( in )         :: lmn                                                    ! Element number
  integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
  integer, intent( in )         :: face                                                   ! Loaded face
  integer, intent( in )         :: n_parameters                                           ! User supplied parameters
  integer, intent( in )         :: n_nodes                                                ! # nodes on the element  
  integer, intent( in )         :: n_properties                                           ! # properties for the element                
  integer, intent( in )         :: length_coord_array                                     ! Total # coords
  integer, intent( in )         :: length_dof_array                                       ! Total # DOF

  real (prec), intent( in )     :: dload_parameters(n_parameters)                         ! User supplied parameters.

  type (node), intent( in )     :: node_property_list(n_nodes)                            ! Data structure describing storage for nodal variables - see below
!  type node
!      sequence
!      integer :: flag                          ! Integer identifier
!      integer :: coord_index                   ! Index of first coordinate in coordinate array
!      integer :: n_coords                      ! Total no. coordinates for the node
!      integer :: dof_index                     ! Index of first DOF in dof array
!      integer :: n_dof                         ! Total no. of DOF for node
!   end type node 
!   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

  real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
  real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
  real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

  real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
                
  real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
  element_residual = 0.d0

  end subroutine user_distributed_load_dynamic
