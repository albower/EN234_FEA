   subroutine user_mesh(n_parameters,parameterlist)
   
   use Types
   use Mesh
   implicit none
   
   integer, intent(in)  :: n_parameters                      ! No. parameters provided by user in input file
   real (prec), intent(in) :: parameterlist(n_parameters)    ! List of parameters provided by user 

   integer :: n_nodes_per_element,n_dof_per_node, butler_volmer_flag
   integer :: lmn,n,iof
   integer :: status
   
   real (prec) :: H
   

 !  User subroutine to generate a mesh.  Note that the user is responsible for allocating all memory  

!   You must specify values for all the following variables:
!   n_zones                     (no. zones, or groups of elements with identical properties, in the mesh)
!   n_nodes                     (No. nodes in the mesh)
!   n_elements                  (No. elements in the mesh)
!   length_dofs                  Total no. DOFs (not including Lagrange multipliers for constraints)
!   length_coords                sum_nodes (# coords for node)
!   length_connectivity          sum_elements (# nodes on element)
!   length_state_variables       sum_elements (# state variables for element) (must be at least 1)
!   length_element_properties    sum_property lists (# properties)
!
!   You must allocate the following memory
!   allocate(connectivity(length_connectivity), stat=status)
!   allocate(element_properties(length_element_properties), stat=status)
!   allocate(coords(length_coords), stat=status)
!   allocate(dof_total(length_dofs), stat=status)
!   allocate(dof_increment(length_dofs), stat=status)
!   allocate(rforce(length_dofs), stat=status)
!   allocate(node_list(n_nodes), stat=status)
!   allocate(element_list(n_elements), stat=status)
!   allocate(zone_list(n_zones), stat=status)
!   allocate(initial_state_variables(length_state_variables), stat = status)
!   allocate(updated_state_variables(length_state_varialbes), stat = status)
!
!  You must define values for the following variables:
!   dof_total               (1D array; stored with u1,u2,u3 for first node, u1, u2,u3 for second node etc, set to initial value of total DOF)
!   dof_increment           (set to initial value for increment in DOF)
!   coords                  (1D array; stored with x1, x2, x3 for first node, x1,x2,x3 for second node, etc. Enter list of corrds for each node)
!   element_properties      (1D array; list properties 1,2,3,etc for first group of elements; props 1,2,3 for second group, etc)
!
!   For each 'zone' (group of elements with identical parameters) you must define
!
!     zone_list(z)%start_element = number of first element in the zone
!     zone_list(z)%end_element = number of last element in the zone.
!
!   For each node you must initialize the following variables
!
!     node_list(n)%flag = (integer identifier for the node)
!     node_list(n)%coord_index = (index of first coord for the node in array coords)
!     node_list(n)%n_coords= (# coords for this node)
!     node_list(n)%dof_index = index of first DOF for the node in array dof_total and dof_increment
!     node_list(n)%n_dof = # DOF for the node
!     node_list(n)%displacement_map_index = 0
!     node_list(n)%n_displacements = 0 
!
!  For each element you must initialize the following variables
!
!     element_list(lmn)%flag = integer flag specifying the element type
!     element_list(lmn)%connect_index = index of first node on the element in the array connectivity
!     element_list(lmn)%n_nodes = # nodes on this element
!     element_list(lmn)%state_index = index of first state variable on the element in the array initial_state_variables (must be a least 1)
!     element_list(lmn)%n_states = # state vars for the element (may be zero)
!     element_list(lmn)%element_property_index = index of first property for the element in the array element_properties
!     element_list(lmn)%n_element_properties = # properties for the element
!
!  You must initialize the following variables:
!     coords
!     connectivity
!     element_properties
!     dof_total
!     dof_increment
!     initial_state_variables
     
      

   
   end subroutine user_mesh
