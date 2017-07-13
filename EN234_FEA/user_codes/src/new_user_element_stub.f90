!
!     Subroutines to assemble element level matrices for FEA analysis

!=========================== subroutine user_element_stiffness ===================
subroutine new_user_element_static(lmn, element_identifier, n_nodes, node_property_list, &       ! Input variables
    n_properties, element_properties, n_int_properties, int_element_properties, &                ! Input variables
    element_coords, length_coord_array, &                                                        ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables
    use Types
    use ParamIO
    use Globals, only: TIME, DTIME                  ! Total analysis time and time increment
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: n_int_properties                                       ! # integer valued properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    integer, intent( in )         :: int_element_properties(n_int_properties)               ! Integer valued element properties

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

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout )  :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)


    element_stiffness = 0.d0
    element_residual = 0.d0
    fail = .false.
    
    return
end subroutine new_user_element_static


!=========================== subroutine user_element_dynamic ===================
subroutine new_user_element_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, n_int_properties, int_element_properties, &                     ! Input variables
    element_coords, length_coord_array, &                                                             ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                     ! Input variables
    n_state_variables, initial_state_variables, &                                                     ! Input variables
    updated_state_variables,element_residual,element_deleted)      ! Output variables
    use Types
    use ParamIO
    use Globals, only: TIME, DTIME                  ! Total analysis time and time increment
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: n_int_properties                                       ! # integer valued properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    integer, intent( in )         :: int_element_properties(n_int_properties)               ! Integer valued element properties

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
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
              
    real( prec ), intent( inout )  :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    logical, intent( inout )       :: element_deleted                                       ! Set to .true. to delete an element

    element_residual = 0.d0
    
    return

end subroutine new_user_element_dynamic

subroutine new_user_element_fieldvariables(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, n_int_properties, int_element_properties, &                       ! Input variables
    element_coords, length_coord_array, &                                                               ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                            ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                                    ! Input variables
    n_field_variables,field_variable_names, &                                                                ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: n_int_properties                                       ! # integer valued properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! No. projected field variables

    integer, intent( in )         :: int_element_properties(n_int_properties)               ! Integer valued element properties

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

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Element stiffness (ROW,COLUMN)

    nodal_fieldvariables = 0.d0
    
    return

end subroutine new_user_element_fieldvariables






