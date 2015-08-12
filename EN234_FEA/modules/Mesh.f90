module Mesh
   use Types

!  Data type for nodes   
   type node
      sequence
      integer :: flag                          ! Integer identifier
      integer :: coord_index                   ! Index of first coordinate in coordinate array
      integer :: n_coords                      ! Total no. coordinates for the node
      integer :: dof_index                     ! Index of first DOF in dof array
      integer :: n_dof                         ! Total no. of DOF for node
      integer :: displacement_map_index        ! Index of displacement node map
      integer :: n_displacements               ! No. displacement DOFs
   end type node
!  Data type for elements   
   type element
      sequence
      integer :: flag                          ! Integer identifier for element
      integer :: connect_index                 ! Index of first node on element in connectivity(:)
      integer :: n_nodes                       ! No. nodes on the element
      integer :: state_index                   ! Index of first state variable in element_state_variables(:)
      integer :: n_states                      ! No. state variables
      integer :: element_property_index        ! Index of first element property in element_properties(:)
      integer :: n_element_properties          ! No. element properties
      integer :: density_index                 ! Index of density value
   end type element
   
   type zone
      sequence
      integer :: start_element                 ! First element in a zone
      integer :: end_element                   ! Last element in a zone
   end type zone
   
   integer, save :: n_zones                          ! Number of zones
   integer, save :: n_nodes                          ! Total number of nodes
   integer, save :: n_elements                       ! Total number of elements
   integer, save :: length_coords                    ! Length of coordinate array
   integer, save :: length_dofs                      ! Length of nodal DOF array
   integer, save :: length_connectivity              ! Length of connectivity array
   integer, save :: length_element_properties        ! Length of element property array
   integer, save :: length_densities                 ! Length of density array
   integer, save :: length_state_variables           ! Length of state variable array
   integer, save :: length_displacement_map          ! Length of array mapping nodal DOF to displacements
   integer, save :: n_mesh_parameters                ! Parameeters controlling a user-subrouine generated mesh
   
   integer, save, allocatable :: displacement_map(:) ! Array storing mapping of nodal DOF to displacements
   integer, save, allocatable :: connectivity(:)     ! Array storing element connectivity

   real (prec), save :: nodal_force_norm             ! Norm of nodal generalized force
   real (prec), save :: unbalanced_force_norm        ! Norm of out-of-balance forces
   real (prec), save :: correction_norm              ! Norm of solution correction

   real (prec), save, allocatable :: element_properties(:)            ! List of element properties
   real (prec), save, allocatable :: densities(:)                     ! List of density values for zones
   real (prec), save, allocatable :: initial_state_variables(:)       ! Element state variables at the start of a time increment
   real (prec), save, allocatable :: updated_state_variables(:)       ! Element state variables at the end of a time increment
   real (prec), save, allocatable :: coords(:)                        ! List of nodal coordinates
   real (prec), save, allocatable :: dof_total(:)                     ! List of accumulated DOF
   real (prec), save, allocatable :: dof_increment(:)                 ! List of increment in DOF
   real (prec), save, allocatable :: velocity(:)                      ! Velocity (for explicit dynamics)
   real (prec), save, allocatable :: acceleration(:)                  ! Acceleration (for explicit dynamics)
   real (prec), save, allocatable :: lumped_mass(:)                   ! Lumped mass matrix (for explicit dynamics)
   real (prec), save, allocatable :: rforce(:)
   real (prec), save, allocatable :: mesh_subroutine_parameters(:)    ! List of parameters controlling a user-subroutine generated mesh
   
   type (node), save, allocatable :: node_list(:)
   type (element), save, allocatable :: element_list(:)
   type (zone), save, allocatable :: zone_list(:)
   
   character (len=100), allocatable :: zone_namelist(:)         ! Names of zones in mesh



end module
