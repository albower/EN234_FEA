!============================SUBROUTINE generate_node_numbers ==========================
subroutine generate_node_numbers(element_start,element_end,include_constraints,node_numbers,node_order_index)
    use Types
    use ParamIO
    use Bandwidth
    use Mesh, only : n_nodes
    use Boundaryconditions, only : n_constraints
    !    use Stiffness, only: node_numbers, node_order_index
    implicit none

    integer, intent(in) :: element_start
    integer, intent(in) :: element_end

    logical, intent(in) :: include_constraints

    integer, intent(out) :: node_numbers(*)
    integer, intent(out) :: node_order_index(*)
  
    ! Local Variables
    integer :: ipair, iwumax, iwvmax, n
    integer :: nu, nv
  
    integer :: status
  
    integer :: n_levels_v
    integer, allocatable :: node_level_v(:)
    integer, allocatable :: level_widths_v(:)
    integer, allocatable :: level_list_index_v(:)
    type (Integer_Linked_List), allocatable :: level_list_v(:)
  
    integer :: n_levels_u
    integer, allocatable :: node_level_u(:)
    integer, allocatable :: level_widths_u(:)
    integer, allocatable :: level_list_index_u(:)
    type (Integer_Linked_List), allocatable :: level_list_u(:)

    integer, allocatable :: node_level(:)
    integer, allocatable :: level_widths(:)
    integer, allocatable :: level_list_index(:)
    type (Integer_Linked_List), allocatable :: level_list(:)
  
  
    !     Bandwidth minimization algorithm of Gibbs et al, SIAM J NUMER ANAL 13 (2) 1976

    ! Local variables deallocated after use

    allocate(node_level_u(n_nodes+n_constraints), stat=status)
    allocate(level_widths_u(n_nodes+n_constraints), stat=status)
    allocate(level_list_index_u(n_nodes+n_constraints), stat=status)
    allocate(level_list_u((n_nodes+n_constraints)), stat=status)
    level_widths_u(1:n_nodes+n_constraints) = 0
    level_list_u(1:(n_nodes+n_constraints))%n_entries = 0
    level_list_u(1:(n_nodes+n_constraints))%next = 0
    level_list_index_u(1:n_nodes+n_constraints) = 0


    allocate(node_level_v(n_nodes+n_constraints), stat=status)
    allocate(level_widths_v(n_nodes+n_constraints), stat=status)
    allocate(level_list_index_v(n_nodes+n_constraints), stat=status)
    allocate(level_list_v((n_nodes+n_constraints)), stat=status)
    level_widths_v(1:n_nodes+n_constraints) = 0
    level_list_v(1:(n_nodes+n_constraints))%n_entries = 0
    level_list_v(1:(n_nodes+n_constraints))%next = 0
    level_list_index_v(1:n_nodes+n_constraints) = 0

    allocate(node_level(n_nodes+n_constraints), stat=status)
    allocate(level_widths(n_nodes+n_constraints), stat=status)
    allocate(level_list_index(n_nodes+n_constraints), stat=status)
    allocate(level_list((n_nodes+n_constraints)), stat=status)
    level_widths(1:n_nodes+n_constraints) = 0
    level_list(1:(n_nodes+n_constraints))%n_entries = 0
    level_list(1:(n_nodes+n_constraints))%next = 0
    level_list_index(1:n_nodes+n_constraints) = 0

    if (status/=0) then
        write(IOW,*) ' Error in subroutine generate_node_numbers '
        write(IOW,*) ' Unable to allocate storage for bandwidth minimizer'
        stop
    endif
  
    !     Generate list of adjacent nodes and degree for each node
    call compute_adjacency(element_start,element_end,include_constraints)
    !     Find endpoints of pseudo-diameter
    call diam(n_nodes+n_constraints, nv, n_levels_v, node_level_v, level_list_index_v,  level_list_v, level_widths_v,  &
        iwvmax, nu, n_levels_u, node_level_u, level_list_index_u, level_list_u, level_widths_u, iwumax)
    !     Generate new level structure of minimum width
    call widmin(n_nodes+n_constraints, node_level_v, iwvmax, n_levels_u, node_level_u,  &
        iwumax, node_level, level_list_index, level_list, level_widths, ipair)
    !     Number nodes based on final level structure
    call numgen(n_nodes+n_constraints, ipair, nu, nv, level_list_index, level_list, level_widths, n_levels_v, node_numbers)
    !     Generate an index table ordering the nodes
    call isort(n_nodes+n_constraints, node_numbers, node_order_index)
    !     Check for errors
    do n = 1, n_nodes+n_constraints
        if ( node_numbers(n)>n_nodes+n_constraints ) then
            write (IOW, 99002) n, node_numbers(n)
            stop
        end if
        if ( node_numbers(n)<1.and.node_degree(n)>0 ) then
            write (IOW, 99003) n, node_numbers(n)
            stop
        end if
    end do

    deallocate(node_level_u)
    deallocate(level_widths_u)
    deallocate(level_list_index_u)
    deallocate(level_list_u)

    deallocate(node_level_v)
    deallocate(level_widths_v)
    deallocate(level_list_index_v)
    deallocate(level_list_v)

    deallocate(node_level)
    deallocate(level_widths)
    deallocate(level_list_index)
    deallocate(level_list)
    deallocate(node_degree)
    deallocate(node_adjacency_index)
    deallocate(node_adjacency)

99002 format ( //, '    ***  ERROR DETECTED IN RENUM ', /, &
        '       Node ', I4, ' has number ', I4)
99003 format ( //, '    ***  ERROR DETECTED IN RENUM ', /, &
        '       Node ', I4, ' has number ', I4)

end subroutine generate_node_numbers

!====================SUBROUTINE compute_adjacency ======================
subroutine compute_adjacency(element_start,element_end,include_constraints)
    use Types
    use ParamIO
    use Mesh, only : element, connectivity, element_list, n_elements, n_nodes
    use Boundaryconditions, only : constraint,nodeset, constraint_list, nodeset_list, node_lists, n_constraints
    use Bandwidth
    implicit none

    integer, intent(in) :: element_start
    integer, intent(in) :: element_end

    logical, intent(in) :: include_constraints

    logical, allocatable :: nzstiffness(:,:)

    !  Computes adjacency structure for a mesh and boundary conditions

    ! Local Variables
    integer :: i, j
    integer :: lmn, nnodes
    integer :: ns
    integer :: node1, node2
    integer :: nc, mpc, iof, start_index
    integer :: status
    integer :: n1, n2, nnode, iofc
  
   

    allocate(nzstiffness(n_nodes,n_nodes), stat = status)
    if (status/=0) then
        write(IOW,'(A)') ' *** Error in subroutine compute_adjacency *** '
        write(IOW,'(A)') ' Unable to allocate memory for nzstiffness '
        stop
    endif
    nzstiffness = .false.
    if (allocated(node_degree)) deallocate(node_degree)
    allocate(node_degree(n_nodes+n_constraints), stat = status)

    node_degree = 0
 
    do lmn = element_start,element_end
        do n1 = 1,element_list(lmn)%n_nodes
            node1 = connectivity(element_list(lmn)%connect_index+n1-1)
            do n2 = 1,element_list(lmn)%n_nodes
                node2 = connectivity(element_list(lmn)%connect_index+n2-1)
                if (node1==node2) cycle
                if (nzstiffness(node1,node2)) cycle
                nzstiffness(node1,node2) = .true.
                nzstiffness(node2,node1) = .true.
                node_degree(node1) = node_degree(node1) + 1
                node_degree(node2) = node_degree(node2) + 1
            end do
        end do
    end do
   
    if (include_constraints) then
        do mpc = 1, n_constraints
            if (constraint_list(mpc)%flag<3) then
                node1 = constraint_list(mpc)%node1
                node2 = constraint_list(mpc)%node2
                if (.not.nzstiffness(node1,node2)) then
                    if (node1/=node2) then
                        node_degree(node1) = node_degree(node1) + 1
                        node_degree(node2) = node_degree(node2) + 1
                        nzstiffness(node1,node2) = .true.
                        nzstiffness(node2,node1) = .true.
                    endif
                endif
                node_degree(node1) = node_degree(node1) + 1
                node_degree(node2) = node_degree(node2) + 1
                node_degree(n_nodes+mpc)=node_degree(n_nodes+mpc)+2  ! Degree for Lagrange multiplier
            else
                ns = constraint_list(mpc)%node1
                nnode = nodeset_list(ns)%n_nodes
                iofc = nodeset_list(ns)%index
                do i = 1, nnode           !     ---    Loop over nodes in set
                    node1 = node_lists(i + iofc - 1)
                    node_degree(node1) = node_degree(node1) + 1
                    do j = 1, nnode
                        node2 = node_lists(j+iofc-1)
                        if (nzstiffness(node1,node2)) cycle
                        if (node1==node2) cycle
                        node_degree(node1) = node_degree(node1) + 1
                        node_degree(node2) = node_degree(node2) + 1
                        nzstiffness(node1,node2) = .true.
                        nzstiffness(node2,node1) = .true.
                    end do
                end do
                node_degree(n_nodes+mpc) = node_degree(n_nodes+mpc)+nnode  ! Degree for Lagrange multiplier
            endif
        end do
    endif

    length_node_adjacency = 0
    do i = 1,n_nodes+n_constraints
        if (node_degree(i)>0)  length_node_adjacency = length_node_adjacency + int(node_degree(i)/databinsize)+1
    end do
    deallocate(nzstiffness)
    max_node_degree = maxval(node_degree)
  
    if (allocated(iadj)) deallocate(iadj)
    if (allocated(node_adjacency)) deallocate(node_adjacency)
    if (allocated(node_adjacency_index)) deallocate(node_adjacency_index)
  
    allocate(iadj(max_node_degree), stat = status)
    allocate(node_adjacency(length_node_adjacency), stat=status)
    allocate(node_adjacency_index(n_nodes+n_constraints), stat=status)
    if (status /=0) then
        write(IOW,'(A)') ' *** Error in subroutine compute_adjacency *** '
        write(IOW,'(A)') ' Unable to allocate memory for node adjacency '
        stop
    endif
    node_adjacency(1:length_node_adjacency)%n_entries = 0
    node_adjacency(1:length_node_adjacency)%next = 0
    node_adjacency_index = 0
    last_filled_node_adjacency = 0
  
    !     Loop over all elements
    do lmn = element_start, element_end

        !     For each node on current element add other nodes
        !     on same element to adjacency list

        do i = 1, element_list(lmn)%n_nodes
            iof = element_list(lmn)%connect_index
            node1 = connectivity(i + iof - 1)
            if (node_adjacency_index(node1)==0) then
                last_filled_node_adjacency = last_filled_node_adjacency + 1
                node_adjacency_index(node1) = last_filled_node_adjacency
            endif
            start_index = node_adjacency_index(node1)
            do j = 1, element_list(lmn)%n_nodes
                node2 = connectivity(j + iof - 1)
                if ( node2==node1 ) cycle
                call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,node2) ! Add node2 to adjacency list for node 1 (if not already present)
            end do
        end do
    end do

    if (include_constraints) then
        do nc = 1,n_constraints    !     Loop over constraint sets
            if (constraint_list(nc)%flag<3) then    ! Simple 2 node constraint
                node1 = constraint_list(nc)%node1
                node2 = constraint_list(nc)%node2
                if (node_adjacency_index(n_nodes+nc)==0) then
                    last_filled_node_adjacency = last_filled_node_adjacency + 1
                    node_adjacency_index(n_nodes+nc) = last_filled_node_adjacency
                endif
                start_index = node_adjacency_index(n_nodes+nc)
                call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,node1)
                call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,node2)
                if (node_adjacency_index(node1)==0) then
                    last_filled_node_adjacency = last_filled_node_adjacency + 1
                    node_adjacency_index(node1) = last_filled_node_adjacency
                endif
                start_index = node_adjacency_index(node1)
                call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,n_nodes+nc)
                if (node_adjacency_index(node2)==0) then
                    last_filled_node_adjacency = last_filled_node_adjacency + 1
                    node_adjacency_index(node2) = last_filled_node_adjacency
                endif
                start_index = node_adjacency_index(node2)
                call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,n_nodes+nc)
                if (node1==node2) cycle
                start_index = node_adjacency_index(node1)
                call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,node2) ! Add node2 to adjacency list for node 1 (if not already present)
                start_index = node_adjacency_index(node2)
                call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,node1) ! Add node1 to adjacency list for node 2 (if not already present)
            else if (constraint_list(nc)%flag==3) then   !  General multi-point constraint
                ns = constraint_list(nc)%node1            !  Node set listing nodes
                iof = nodeset_list(ns)%index
                nnodes = nodeset_list(ns)%n_nodes
                do i = 1,nnodes
                    node1 = node_lists(iof+i-1)
                    if (node_adjacency_index(node1)==0) then
                        last_filled_node_adjacency = last_filled_node_adjacency + 1
                        node_adjacency_index(node1) = last_filled_node_adjacency
                    endif
                    start_index = node_adjacency_index(node1)
                    call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,n_nodes+nc)
                    do j = 1,nnodes
                        node2 = node_lists(iof+j-1)
                        if (node2==node1) cycle
                        call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,node2)
                    end do
                end do
                if (node_adjacency_index(n_nodes+nc)==0) then
                    last_filled_node_adjacency = last_filled_node_adjacency + 1
                    node_adjacency_index(n_nodes+nc) = last_filled_node_adjacency
                endif
                start_index = node_adjacency_index(n_nodes+nc)
                do i = 1,nnodes
                    node1 = node_lists(iof+i-1)
                    call adddistinctlistdata(start_index,node_adjacency,last_filled_node_adjacency,node1)
                end do
            end if
    
        end do
    endif
    return

end subroutine compute_adjacency

!=====================SUBROUTINE DIAM ==========================
subroutine diam(nops, nv, n_levels_v, node_level_v, level_list_index_v,  level_list_v, level_widths_v,  &
    iwvmax, nu, n_levels_u, node_level_u, level_list_index_u, level_list_u, level_widths_u, iwumax)
    use Types
    use ParamIO
    use Linkedlist_Handling
    use Bandwidth, only: ideg => node_degree
    implicit none

    integer, intent( in )    :: nops
    integer, intent( out )   :: nv
    integer, intent( out )   :: iwvmax
    integer, intent( out )   :: nu
    integer, intent( out )   :: iwumax

    integer, intent( out )   :: n_levels_v
    integer, intent( out )   :: node_level_v(nops)
    integer, intent( out )   :: level_widths_v(nops)
    integer, intent( out )   :: level_list_index_v(nops)
    type (Integer_Linked_List), intent(out) :: level_list_v(nops)

    integer, intent( out )   :: n_levels_u
    integer, intent( out )   :: node_level_u(nops)
    integer, intent( out )   :: level_widths_u(nops)
    integer, intent( out )   :: level_list_index_u(nops)
    type (Integer_Linked_List), intent(out) :: level_list_u(nops)

    ! Local Variables
    integer :: i, iwsmax, mindeg, minnod, minrem
    integer :: n, ns
    integer :: status
    integer :: n_nodes
    integer :: start_index
  
    logical :: swapped

    integer              :: n_levels_s
    integer, allocatable :: node_level_s(:)
    integer, allocatable :: level_widths_s(:)
    integer, allocatable :: level_list_index_s(:)
    type (Integer_Linked_List), allocatable :: level_list_s(:)
 
    integer, allocatable :: level_node_list(:)
    integer, allocatable :: sort_index(:)
    integer, allocatable :: id(:)
    integer, allocatable :: workspace(:)

    !     Subroutine to find endpoints NV and NU of pseudo-diameter and associated level structures LU and LV

    allocate(node_level_s(nops), stat=status)
    allocate(level_widths_s(nops), stat=status)
    allocate(level_list_index_s(nops), stat=status)
    allocate(level_list_s(nops), stat=status)
    level_widths_s(1:nops) = 0
    level_list_s(1:nops)%n_entries = 0
    level_list_s(1:nops)%next = 0
    level_list_index_s(1:nops) = 0
  
    allocate(level_node_list(nops), stat=status)
    allocate(sort_index(nops), stat=status)
    allocate(id(nops), stat=status)
    allocate(workspace(nops), stat=status)

    if (status /=0) then
        write(IOW,'(A)') ' *** Error in subroutine diam *** '
        write(IOW,'(A)') ' Unable to allocate memory for bandwidth minimizer '
        stop
    endif

    !     Find an arbitrary node of smallest degree
    mindeg = nops
    minnod = 1
    do i = 1, nops
        if ( ideg(i)<mindeg.and.ideg(i)>0 ) then
            mindeg = ideg(i)
            minnod = i
        end if
    end do

    !     Generate a level structure starting at node of lowest degree
    nv = minnod

    call levgen(nops,ideg, nv, node_level_v, level_list_index_v, level_list_v, level_widths_v, n_levels_v, iwvmax, workspace)

    !     Loop over nodes in highest level in order of increasing degree.  For each node, generate
    !     a new level structure.
  
    swapped = .true.

    !---  Generate index table for nodes in highest level
    do while (swapped)
        start_index = level_list_index_v(n_levels_v)                         ! Extract nodes from previous level
        call extractlinkedlistdata(start_index,level_list_v,level_node_list,n_nodes)
        do i = 1, level_widths_v(n_levels_v)
            id(i) = ideg(level_node_list(i))
        end do

        n = level_widths_v(n_levels_v)
        call isort(n, id, sort_index)

        i = 0
        minrem = nops

        do while (i<level_widths_v(n_levels_v))
            i = i + 1
            ns = level_node_list(sort_index(i))
            call levgen(nops,ideg, ns, node_level_s, level_list_index_s, &
                level_list_s, level_widths_s, n_levels_s, iwsmax, workspace)

            !     If depth of new level structure is greater than that of old one, swap over and start again

            swapped = .false.
            if ( n_levels_s>n_levels_v ) then
                nv = ns
                node_level_v = node_level_s
                level_list_index_v = level_list_index_s
                level_list_v = level_list_s
                level_widths_v = level_widths_s
                n_levels_v = n_levels_s
                swapped = .true.
                exit
            end if

            !     Store level structure with minimum width
            if ( iwsmax<=minrem ) then
                minrem = iwsmax
                nu = ns
 
                node_level_u = node_level_s
                level_list_index_u = level_list_index_s
                level_list_u = level_list_s
                level_widths_u = level_widths_s
                n_levels_u = n_levels_s
                iwumax = iwsmax

            end if

        !     Repeat for next node in highest level of v
        end do

    end do

    deallocate(node_level_s)
    deallocate(level_widths_s)
    deallocate(level_list_index_s)
    deallocate(level_list_s)
    deallocate(level_node_list)
    deallocate(sort_index)
    deallocate(id)
end subroutine diam

!====================SUBROUTINE WIDMIN =======================
subroutine widmin(nops, node_level_v, iwvmax, n_levels_u, node_level_u,  &
    iwumax, node_level, level_list_index, level_list, level_widths, ipair)
    use Types
    use ParamIO
    use Linkedlist_Handling
    use Bandwidth, only: ideg => node_degree
    implicit none

    integer, intent( in )    :: nops
    integer, intent( in )    :: iwvmax
    integer, intent( in )    :: iwumax
    integer, intent( out )   :: ipair
                             
    integer, intent( in )    :: node_level_v(nops)

    integer, intent( in )    :: n_levels_u
    integer, intent( in   )  :: node_level_u(nops)
   
    integer, intent( out )   :: node_level(nops)
    integer, intent( out )   :: level_widths(nops)
    integer, intent( out )   :: level_list_index(nops)
    type (Integer_Linked_List), intent(out) :: level_list(nops)

    ! Local Variables
    integer :: i, icount, ih, ihmax, il, ilev, ilmax,  &
        iwsmax, j, jcount
    integer :: n, ncon, node, nstrt, iof
    integer :: status
    integer :: start_index
    integer :: n_nodes

    integer :: last_filled_connected_component
    integer, allocatable :: ideg2(:)
    integer, allocatable :: connected_component_index(:)
    integer, allocatable :: connected_component_nodelist(:)

    integer :: last_filled_level_list
    integer :: n_levels_s
    integer, allocatable :: node_level_s(:)
    integer, allocatable :: level_widths_s(:)
    integer, allocatable :: level_list_index_s(:)
    type (Integer_Linked_List), allocatable :: level_list_s(:)

    integer, allocatable :: levwh(:)
    integer, allocatable :: levwl(:)
    integer, allocatable :: level_node_list(:)
    integer, allocatable :: n_members(:)
    integer, allocatable :: index(:)
    integer, allocatable :: workspace(:)

    !     Subroutine to assemble nodes into new level structure to minimize level width

    allocate(ideg2(nops), stat=status)
    allocate(connected_component_index(nops), stat=status)
    allocate(connected_component_nodelist(nops), stat=status)
   
    allocate(node_level_s(nops), stat=status)
    allocate(level_widths_s(nops), stat=status)
    allocate(level_list_index_s(nops), stat=status)
    allocate(level_list_s(nops), stat=status)
    allocate(level_node_list(nops), stat=status)
    level_widths_s(1:nops) = 0
    level_list_s(1:nops)%n_entries = 0
    level_list_s(1:nops )%next = 0
    level_list_index_s(1:nops) = 0
  
    allocate(levwh(nops), stat=status)
    allocate(levwl(nops), stat=status)
    allocate(level_node_list(nops), stat=status)
    allocate(n_members(nops), stat=status)
    allocate(index(nops), stat=status)
    allocate(workspace(nops), stat=status)
  
  
    if (status /=0) then
        write(IOW,'(A)') ' *** Error in subroutine widmin *** '
        write(IOW,'(A)') ' Unable to allocate memory for node degree '
        stop
    endif
  

    !     Extract nodes from Lu and Lv with identical level pairs
    level_widths = 0
    last_filled_level_list = 0

    icount = 0
    do i = 1,nops
        if (ideg(i) > 0) icount = icount + 1
    end do
    ideg2(1:nops) = ideg(1:nops)
    do i = 1, nops
        if (ideg(i) > 0) then
            if ( node_level_v(i)==n_levels_u + 1 - node_level_u(i) ) then
                n = node_level_v(i)
                node_level(i) = n
                if (level_widths(n)==0) then
                    last_filled_level_list = last_filled_level_list + 1
                    level_list_index(n) = last_filled_level_list
                endif
                level_widths(n) = level_widths(n) + 1
                start_index = level_list_index(n)
                call addlistdata(start_index,level_list,last_filled_level_list,i)
                icount = icount - 1
                ideg2(i) = 0
            else
                ideg2(i) = ideg(i)
            end if
        endif
    end do

    if ( icount==0 ) return
  
    !     Collect nodes into disjoint connected components
    ncon = 0
    last_filled_connected_component = 0
    do i = 1, nops
        if ( ideg2(i)/=0 ) then
            nstrt = i
            !     ---   Find all nodes connected to current one by generating levels
            call levgen(nops,ideg2, nstrt, node_level_s, level_list_index_s, &
                level_list_s, level_widths_s, n_levels_s, iwsmax, workspace)
            !     ---   Add nodes in level structure to connected component
            ncon = ncon + 1
            connected_component_index(ncon) = last_filled_connected_component+1
            jcount = 0
            do n = 1, n_levels_s
                start_index = level_list_index_s(n)                         ! Extract nodes
                call extractlinkedlistdata(start_index,level_list_s,level_node_list,n_nodes)
                do j = 1, level_widths_s(n)
                    jcount = jcount + 1
                    last_filled_connected_component=last_filled_connected_component+1
                    connected_component_nodelist(last_filled_connected_component) = level_node_list(j)
                    ideg2(level_node_list(j)) = 0
                end do
            end do
            n_members(ncon) = jcount
        end if
    end do

    !     Sort connected components by size
    call isort(ncon, n_members, index)

    !     Loop over connected components
    do i = 1, ncon
        n = index(ncon + 1 - i)

        !     Compute estimates of level width based on u and v structures
        levwh = level_widths
        levwl  = level_widths
        iof = connected_component_index(n)
        do j = 1, n_members(n)
            node = connected_component_nodelist(iof+j-1)
            ih = node_level_v(node)
            il = n_levels_u + 1 - node_level_u(node)
            levwh(ih) = levwh(ih) + 1
            levwl(il) = levwl(il) + 1
        end do

        ihmax = 0
        ilmax = 0
        do j = 1, n_levels_u
            if ( levwh(j)>ihmax .and. levwh(j)>level_widths(j) ) ihmax = levwh(j)
            if ( levwl(j)>ilmax .and. levwl(j)>level_widths(j) ) ilmax = levwl(j)
        end do

        !     Assign nodes to levels based on calculated widths
        if ( ilmax>ihmax ) then
            !---  Put nodes in levels of v
            if ( i==1 ) ipair = 1
            iof = connected_component_index(n)
            do j = 1, n_members(n)
                node = connected_component_nodelist(iof+j-1)
                ilev = node_level_v(node)
                node_level(node) = ilev
                if (level_widths(ilev)==0) then
                    last_filled_level_list = last_filled_level_list + 1
                    level_list_index(ilev) = last_filled_level_list
                endif
                start_index = level_list_index(ilev)
                call addlistdata(start_index,level_list,last_filled_level_list,node)
                level_widths(ilev) = level_widths(ilev) + 1
            end do
        else if ( ilmax<ihmax ) then
            !---  Put nodes in levels of u
            if ( i==1 ) ipair = 2
            iof = connected_component_index(n)
            do j = 1, n_members(n)
                node = connected_component_nodelist(iof+j-1)
                ilev = n_levels_u + 1 - node_level_u(node)
                node_level(node) = ilev
                if (level_widths(ilev)==0) then
                    last_filled_level_list = last_filled_level_list + 1
                    level_list_index(ilev) = last_filled_level_list
                endif
                start_index = level_list_index(ilev)
                call addlistdata(start_index,level_list,last_filled_level_list,node)
                level_widths(ilev) = level_widths(ilev) + 1
            end do
        else if ( iwvmax<=iwumax ) then
            !---  Put nodes in levels of v
            if ( i==1 ) ipair = 1
            iof = connected_component_index(n)
            do j = 1, n_members(n)
                node = connected_component_nodelist(iof+j-1)
                ilev = node_level_v(node)
                node_level(node) = ilev
                if (level_widths(ilev)==0) then
                    last_filled_level_list = last_filled_level_list + 1
                    level_list_index(ilev) = last_filled_level_list
                endif
                start_index = level_list_index(ilev)
                call addlistdata(start_index,level_list,last_filled_level_list,node)
                level_widths(ilev) = level_widths(ilev) + 1
            end do
        else
            !---  Put nodes in levels of u
            if ( i==1 ) ipair = 2
            iof = connected_component_index(n)
            do j = 1, n_members(n)
                node = connected_component_nodelist(iof+j-1)
                ilev = n_levels_u + 1 - node_level_u(node)
                node_level(node) = ilev
                if (level_widths(ilev)==0) then
                    last_filled_level_list = last_filled_level_list + 1
                    level_list_index(ilev) = last_filled_level_list
                endif
                start_index = level_list_index(ilev)
                call addlistdata(start_index,level_list,last_filled_level_list,node)
                level_widths(ilev) = level_widths(ilev) + 1
            end do
        end if

    end do
  
    deallocate(ideg2)
    deallocate(connected_component_index)
    deallocate(connected_component_nodelist)
   
    deallocate(node_level_s)
    deallocate(level_widths_s)
    deallocate(level_list_index_s)
    deallocate(level_list_s)
    deallocate(level_node_list)
  
    deallocate(levwh)
    deallocate(levwl)
    deallocate(n_members)
    deallocate(index)

end subroutine widmin

!======================SUBROUTINE NUMGEN ======================
subroutine numgen(nops, ipair, nu, nv, level_list_index, level_list, level_widths, n_levels, numnod)
    use Types
    use ParamIO
    use Linkedlist_Handling
    use Bandwidth, only: ideg => node_degree
    implicit none

    integer, intent( in )    :: nops
    integer, intent( inout ) :: ipair
    integer, intent( inout ) :: nu
    integer, intent( inout ) :: nv
    integer, intent( in )    :: n_levels
    integer, intent( in )    :: level_widths(nops)
    integer, intent( in )    :: level_list_index(nops)
    type (Integer_Linked_List), intent(in) :: level_list(nops)
                                       
    integer, intent( inout ) :: numnod(*)

    ! Local Variables
    integer :: i, irev, l, l1, lcount,   &
        lend, linc,  lstart, mindeg, n, node
    integer :: nodrem, number, status
    integer :: start_index
    integer :: n_nodes1, n_nodes2

    integer, allocatable :: list(:)
    integer, allocatable :: ldeg(:)
    integer, allocatable :: sort_index(:)
  
    integer, allocatable :: node_list_lev1(:)
    integer, allocatable :: node_list_lev2(:)
  
    integer, allocatable :: ltemp(:)

    !     Subroutine to assign node numbers based on final level structure

    allocate(list(nops), stat = status)
    allocate(ldeg(nops), stat=status)
    allocate(sort_index(nops), stat=status)
    allocate(node_list_lev1(nops), stat=status)
    allocate(node_list_lev2(nops), stat=status)
    allocate(ltemp(nops+1), stat=status)
  
    if (status/=0) then
        write(IOW,'(A)') ' *** Error in subroutine numgen *** '
        write(IOW,'(A)') '     Unable to allocate memory for node numbering '
        stop
    endif

    numnod(1:nops) = 0


    if ( ideg(nv)<=ideg(nu) ) then
        numnod(nv) = 1
        irev = 0
    else
        numnod(nu) = 1
        irev = 1
    end if
    number = 1

    !     Number nodes in first (last) level

    l = 1
    if ( irev==1 ) l = n_levels
    do while ( .true. )
        start_index = level_list_index(l)                         ! Extract nodes
        call extractlinkedlistdata(start_index,level_list,node_list_lev1,n_nodes1)
    
        call lowfnd(nops, node_list_lev1, n_nodes1, node_list_lev1, n_nodes1, ltemp, &
            numnod, list, ldeg, lcount)

        if ( lcount>0 ) then
            call isort(lcount, ldeg, sort_index)

            do i = 1, lcount
                n = sort_index(i)
                node = list(n)
                number = number + 1
                numnod(node) = number
            end do
            cycle
        end if

        !     Check for unnumbered nodes in first level

        nodrem = 0
        mindeg = nops
        do i = 1, level_widths(l)
            node = node_list_lev1(i)
            if ( numnod(node)==0 ) then
                if ( ideg(node)<mindeg.and.ideg(node)>0 ) then
                    mindeg = ideg(node)
                    nodrem = node
                end if
            end if
        end do

        if ( nodrem>0 ) then
            number = number + 1
            numnod(nodrem) = number
            cycle
        end if

        if ( irev==0 ) then
            lstart = 2
            lend = n_levels
            linc = 1
        else
            lstart = n_levels - 1
            lend = 1
            linc = -1
        end if

        !     Number nodes in remaining levels

        do l = lstart, lend, linc

            l1 = l - linc
            do while ( .true. )

                start_index = level_list_index(l1)                         ! Extract nodes
                call extractlinkedlistdata(start_index,level_list,node_list_lev1,n_nodes1)

                start_index = level_list_index(l)                         ! Extract nodes
                call extractlinkedlistdata(start_index,level_list,node_list_lev2,n_nodes2)

                call lowfnd(nops, node_list_lev1, n_nodes1, node_list_lev2, n_nodes2, ltemp, &
                    numnod, list, ldeg, lcount)

                if ( lcount>0 ) then
                    call isort(lcount, ldeg, sort_index)

                    do i = 1, lcount
                        n = sort_index(i)
                        node = list(n)
                        number = number + 1
                        numnod(node) = number
                    end do
                    cycle
                end if

                !     Check for unnumbered nodes in current level

                nodrem = 0
                mindeg = nops
                do i = 1, level_widths(l)
                    node = node_list_lev2(i)
                    if ( numnod(node)==0 ) then
                        if ( ideg(node)<=mindeg.and.ideg(node)>0 ) then
                            mindeg = ideg(node)
                            nodrem = node
                        end if
                    end if
                end do

                if ( nodrem>0 ) then
                    number = number + 1
                    numnod(nodrem) = number
                    cycle
                end if
                exit
            end do

        end do

        !     Reverse numbering if appropriate

        nodrem = 0
        do i = 1,nops
            if (ideg(i)>0) nodrem = nodrem + 1
        end do

        if ( ((irev==1.) .and. (ipair==2)) .or.  &
            ((irev==0.) .and. (ipair==1)) ) then
            do i = 1, nops
                if (ideg(i) > 0) numnod(i) = nodrem + 1 - numnod(i)
            end do
        end if

        return
    end do

    deallocate(list)
    deallocate(ldeg)
    deallocate(sort_index)
    deallocate(node_list_lev1)
    deallocate(node_list_lev2)
    deallocate(ltemp)


end subroutine numgen


!=======================Subroutine LOWFND ======================
subroutine lowfnd(nops, lev1, levw1, lev2, levw2, ltemp, numnod, list, ldeg, lcount)
    use Types
    use Bandwidth
    use Linkedlist_Handling
    implicit none

    integer, intent( in )    :: nops
    integer, intent( in )    :: levw1
    integer, intent( in )    :: levw2
    integer, intent( out )   :: lcount
    integer, intent( in )    :: lev1(nops)
    integer, intent( in )    :: lev2(nops)
    integer, intent( inout ) :: numnod(nops)
    integer, intent( out )   :: list(nops)
    integer, intent( out )   :: ldeg(nops)
    integer, intent( inout ) :: ltemp(nops)

    ! Local Variables
    integer :: j, k, lctemp,  low, n, nbr, node1, node2
    integer :: start_index
    integer :: n_datavalues


    !     Function to find lowest numbered node
    !     of level LEV1 which has unnumbered nodes
    !     in level LEV2 adjacent to it, and return list
    !     of adjacent unnumbered nodes

 
    low = nops + 1
    lcount = 0
    !     ---Loop over nodes in level LEV1
    do n = 1, levw1
        node1 = lev1(n)
        !     -- If current node has lower number, then...
        if ( (numnod(node1)<low) .and. (numnod(node1)>0) ) then
            lctemp = 0
            !     --  Check all adjoining nodes
            start_index = node_adjacency_index(node1)
            call extractlinkedlistdata(start_index,node_adjacency,iadj,n_datavalues)
            do j = 1, node_degree(node1)
                nbr = iadj(j)
                if ( numnod(nbr)==0 ) then
                    !     --    If adjoining node is in level LEV2, then...
                    do k = 1, levw2
                        node2 = lev2(k)
                        !     --   add adjoining unnumbered node to list
                        if ( node2==nbr ) then
                            lctemp = lctemp + 1
                            ltemp(lctemp) = node2
                        end if
                    end do
                end if
            end do
            if ( lctemp>0 ) then
                low = numnod(node1)
                lcount = lctemp
                do j = 1, lcount
                    list(j) = ltemp(j)
                    ldeg(j) = node_degree(list(j))
                end do
            end if
        end if
    end do

end subroutine lowfnd

!======================Subroutine LEVGEN ================
subroutine levgen(nops,ideg, start_node, node_levels, level_list_index, level_list, level_widths, n_levels, iwidth, level_node_list)
    use Types
    use ParamIO
    use Linkedlist_Handling
    use Bandwidth, only: node_adjacency_index,node_adjacency,iadj
    implicit none

    integer, intent( in )  :: nops
    integer, intent( in )  :: start_node
    integer, intent( out ) :: n_levels
    integer, intent( out ) :: iwidth
    integer, intent( in )  :: ideg(nops)
    integer, intent( out ) :: node_levels(nops)
    integer, intent( out ) :: level_widths(nops)
    integer, intent( out ) :: level_list_index(nops)
    type (Integer_Linked_List), intent(out) :: level_list(nops)
    integer, intent(inout) :: level_node_list(nops)
                            

    ! Local Variables
    integer :: i, j, nbr, nl1, node
    integer :: start_index, n_datavalues
    integer :: last_filled_level_list
    integer :: n_nodes
  

    !     Subroutine to generate a level structure starting at node start_node


    node_levels = 0
    level_widths = 0
    level_list(1:nops)%n_entries = 0
    level_list(1:nops )%next = 0

    !     Fist level contains only start node

    node_levels(start_node) = 1
    level_widths(1) = 1
    n_levels = 1
    last_filled_level_list = 1
    level_list_index(1) = 1
    start_index = 1
    call addlistdata(start_index,level_list,last_filled_level_list,start_node)
  
    do while ( .true. )

        nl1 = n_levels
        n_levels = n_levels + 1
        if (n_levels>nops) exit
        last_filled_level_list = last_filled_level_list + 1
        level_list_index(n_levels) = last_filled_level_list             ! Create a new level
    
        !     Loop over nodes in preceding level, adding adjoining nodes to next level

        start_index = level_list_index(nl1)                         ! Extract nodes from previous level
        call extractlinkedlistdata(start_index,level_list,level_node_list,n_nodes)

        do i = 1, level_widths(nl1)
            node = level_node_list(i)
            !     --- Loop over adjoining nodes
            start_index = node_adjacency_index(node)
            call extractlinkedlistdata(start_index,node_adjacency,iadj,n_datavalues)
            do j = 1, ideg(node)
                nbr = iadj(j)
                !     ---   If adjoining node has not been assigned a level, add to next level
                if ( node_levels(nbr)==0 .and. ideg(nbr)>0 ) then
                    level_widths(n_levels) = level_widths(n_levels) + 1
                    start_index = level_list_index(n_levels)
                    call addlistdata(start_index,level_list,last_filled_level_list,nbr)
                    node_levels(nbr) = n_levels
                end if
            end do
        end do

        if ( level_widths(n_levels)<=0 ) exit
    end do

    n_levels = n_levels - 1
    iwidth = maxval(level_widths(1:n_levels))

end subroutine levgen

!=======================SUBROUTINE ISORT ========================
subroutine isort(n, in, index)
    use Types
    implicit none

    integer, intent( in )  :: n
    integer, intent( in )  :: in(*)
    integer, intent( out ) :: index(*)

    ! Local Variables
    integer :: i, indxt, ir, j, k, l

    !     Subroutine to generate an index listing
    !     elements of IN in increasing magnitude
    do j = 1, n
        index(j) = j
    end do

    if ( n==1 ) return
    l = n/2 + 1
    ir = n
100 if ( l>1 ) then
        l = l - 1
        indxt = index(l)
        k = in(indxt)
    else
        indxt = index(ir)
        k = in(indxt)
        index(ir) = index(1)
        ir = ir - 1
        if ( ir==1 ) then
            index(1) = indxt
            return
        end if
    end if
    i = l
    j = l + l
    do while ( .true. )
        if ( j<=ir ) then
            if ( j<ir ) then
                if ( in(index(j))<in(index(j + 1)) ) j = j + 1
            end if
            if ( k<in(index(j)) ) then
                index(i) = index(j)
                i = j
                j = j + j
            else
                j = ir + 1
            end if
            cycle
        end if
        index(i) = indxt
        goto 100
    end do

end subroutine isort
