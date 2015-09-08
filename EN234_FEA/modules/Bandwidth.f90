module Bandwidth
    use Types
    use Linkedlist_Handling
   
   
    integer, save :: last_filled_node_adjacency
    integer, save :: length_node_adjacency
    integer, save :: max_node_degree
  
    integer, save, allocatable :: iadj(:)                  ! List of nodes adjacent to a single node
  
    integer, save, allocatable :: node_degree(:)
    integer, save, allocatable :: node_adjacency_index(:)  !  Points to first data block of adjacency for each node
    type (Integer_Linked_List), save, allocatable :: node_adjacency(:)  ! Linked data block list of adjacency
   
end module
