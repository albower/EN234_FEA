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
     
      
   n_nodes_per_element = int(parameterlist(1))
   n_elements = int(parameterlist(2))
   n_dof_per_node = int(parameterlist(3))
   butler_volmer_flag = int(parameterlist(4))
   H = parameterlist(5)                                    ! Film thickness
   n_nodes = n_elements*(n_nodes_per_element-1) + 1 + butler_volmer_flag
   
   n_zones = 1+butler_volmer_flag
   length_dofs = n_dof_per_node*n_nodes - butler_volmer_flag
   length_connectivity = n_nodes_per_element*n_elements+2*butler_volmer_flag
   length_element_properties = n_parameters-5
   
   allocate(connectivity(length_connectivity), stat=status)
   allocate(element_properties(length_element_properties), stat=status)
   allocate(coords(n_nodes-butler_volmer_flag), stat=status)
   allocate(dof_total(length_dofs), stat=status)
   allocate(dof_increment(length_dofs), stat=status)
   allocate(rforce(length_dofs), stat=status)
   allocate(node_list(n_nodes), stat=status)
   allocate(element_list(n_elements+butler_volmer_flag), stat=status)
   allocate(zone_list(n_zones), stat=status)
   allocate(initial_state_variables(1), stat = status)
   allocate(updated_state_variables(1), stat = status)
   dof_total = 0.d0
   dof_increment = 0.d0
     
   zone_list(1)%start_element = 1
   zone_list(1)%end_element = n_elements

   if (butler_volmer_flag==1) then
      zone_list(2)%start_element = n_elements+1
      zone_list(2)%end_element = n_elements+1
  ! If we are enforcing the Butler-Volmer boundary condition the first node is the voltage node
      node_list(1)%flag = 0
      node_list(1)%coord_index = 0
      node_list(1)%n_coords=0
      node_list(1)%dof_index = 1
      node_list(1)%n_dof = 1
      node_list(1)%displacement_map_index = 0
      node_list(1)%n_displacements = 0
   endif
   
   
   ! Bulk nodes
   do n = 1+butler_volmer_flag,n_nodes
     node_list(n)%flag = 1
     node_list(n)%coord_index = n-butler_volmer_flag
     node_list(n)%n_coords=1
     node_list(n)%dof_index = (n-1)*n_dof_per_node+1-butler_volmer_flag
     node_list(n)%n_dof = n_dof_per_node
     node_list(n)%displacement_map_index = 0
     node_list(n)%n_displacements = 0
     coords(n-1) = (n-1-butler_volmer_flag)*H/(n_nodes-1)
   end do
!   do n = 1,n_nodes/2
!     dof_total(2*n) = parameterlist(8)      ! High concentration phase
!   end do
   do n = 1+butler_volmer_flag,n_nodes
     dof_total(2*n-butler_volmer_flag) = parameterlist(8)/6      ! Low concentration phase
   end do
!  These are the bulk elements   
   do lmn = 1,n_elements
     element_list(lmn)%flag = 500
     element_list(lmn)%connect_index = (lmn-1)*n_nodes_per_element+1
     element_list(lmn)%n_nodes = n_nodes_per_element
     element_list(lmn)%state_index = 1
     element_list(lmn)%n_states = 0
     element_list(lmn)%element_property_index = 1
     element_list(lmn)%n_element_properties = 7
     iof = (lmn-1)*n_nodes_per_element+1
     connectivity(iof) = (lmn-1)*(n_nodes_per_element-1) + 1 + butler_volmer_flag
     connectivity(iof+1) = lmn*(n_nodes_per_element-1) + 1 + butler_volmer_flag
     if (n_nodes_per_element==3) then
       connectivity(iof+2) = lmn*(n_nodes_per_element-1)
     endif
   end do
!  The last element is the butler-volmer element   
   if (butler_volmer_flag==1) then
     n_elements = n_elements + 1
     lmn = n_elements
     element_list(lmn)%flag = 400
     element_list(lmn)%connect_index = (lmn-1)*n_nodes_per_element+1
     element_list(lmn)%n_nodes = 2
     element_list(lmn)%state_index = 1
     element_list(lmn)%n_states = 0
     element_list(lmn)%element_property_index = 8
     element_list(lmn)%n_element_properties = 5
     iof = (lmn-1)*n_nodes_per_element+1
     connectivity(iof) = 1                       ! The first node is the voltage node
     connectivity(iof+1) = 2                     ! Second node is the bulk node
   endif
   element_properties(1:length_element_properties) = parameterlist(6:n_parameters)
   
   
   
   end subroutine user_mesh