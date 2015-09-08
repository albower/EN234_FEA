subroutine check_stiffness(element_flag)
    use Types
    use ParamIO
    use User_Subroutine_Storage
    use Mesh
    implicit none

    integer, intent(in)   :: element_flag

    ! Local Variables
    integer         :: lmn,n,i,j,k,iof
    integer         :: ix,iu,ns
    integer         :: icount,icount2
    integer         :: status

    real ( prec ) :: err
    logical :: fail

    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
                                                            
    real( prec ), allocatable   :: element_stiffness(:,:)
    real( prec ), allocatable   :: stif1(:,:)
    real( prec ), allocatable   :: numerical_stiffness(:,:)
    real( prec ), allocatable   :: element_residual(:)
    real( prec ), allocatable   :: resid0(:)
    real( prec ), allocatable   :: resid1(:)

    type (node), allocatable ::  local_nodes(:)

    !     Subroutine to assemble global stiffness matrix

    allocate(element_coords(length_coord_array), stat = status)
    allocate(element_dof_increment(length_dof_array), stat = status)
    allocate(element_dof_total(length_dof_array), stat = status)
    allocate(local_nodes(length_node_array), stat = status)
    allocate(element_stiffness(length_dof_array,length_dof_array), stat = status)
    allocate(numerical_stiffness(length_dof_array,length_dof_array), stat = status)
    allocate(stif1(length_dof_array,length_dof_array), stat = status)
    allocate(element_residual(length_dof_array), stat = status)
    allocate(resid0(length_dof_array), stat = status)
    allocate(resid1(length_dof_array), stat = status)

    if (status/=0) then
        write(IOW,*) ' Error in subroutine check_stiffness '
        write(IOW,*) ' Unable to allocate memory for user subroutines '
        stop
    endif


  
    !
    ! Find an element of type element_flag
    do lmn = 1, n_elements
        if (element_list(lmn)%flag==element_flag) exit
    end do
    if (lmn>n_elements) then
        write(IOW,*) ' *** Error in subroutine check_stiffness *** '
        write(IOW,*) ' No element with identifier ',element_flag,' was found in mesh '
        stop
    endif

    !     Extract local coords, DOF and nodal state/props for the element
    ix = 0
    iu = 0
    do j = 1, element_list(lmn)%n_nodes
        n = connectivity(element_list(lmn)%connect_index + j - 1)
        local_nodes(j)%n_coords = node_list(n)%n_coords
        if (ix+node_list(n)%n_coords>length_coord_array) then
            write(IOW,*) ' Error in subroutine assemble_direct_stiffness '
            write(IOW,*) ' Insufficient storage for element dof '
            write(IOW,*) ' Parameter length_dof_array in module User_Element_Storage must be increased'
            stop
        endif
        local_nodes(j)%coord_index = ix+1
        do k = 1, node_list(n)%n_coords
            ix = ix + 1
            element_coords(ix) = coords(node_list(n)%coord_index + k - 1)
        end do
        if (iu+node_list(n)%n_dof>length_dof_array) then
            write(IOW,*) ' Error in subroutine assemble_direct_stiffness '
            write(IOW,*) ' Insufficient storage for element coordinates '
            write(IOW,*) ' Parameter length_coord_array in module User_Element_Storage must be increased'
            stop
        endif
        local_nodes(j)%dof_index = iu+1
        do k = 1, node_list(n)%n_dof
            iu = iu + 1
            element_dof_increment(iu) = dof_increment(node_list(n)%dof_index + k - 1)
            element_dof_total(iu) = dof_total(node_list(n)%dof_index + k - 1)
        end do
    end do

    !     Form element stiffness
    iof = element_list(lmn)%state_index
    if (iof==0) iof = 1
    ns = element_list(lmn)%n_states
    call user_element_static(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes, &
        local_nodes(1:element_list(lmn)%n_nodes), &       ! Input variables
        element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &     ! Input variables
        element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu), iu,      &                              ! Input variables
        ns, initial_state_variables(iof:iof+ns), &                           ! Input variables
        updated_state_variables(iof:iof+ns),element_stiffness(1:iu,1:iu),resid0(1:iu), fail)               ! Output variables

    !     Compute numerical derivative of stiffness

    icount = 0
    do n = 1,element_list(lmn)%n_nodes
        do i = 1,node_list(n)%n_dof
            icount = icount + 1
            element_dof_increment(icount) = element_dof_increment(icount) + 1.D-07

            call user_element_static(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes, &
                 local_nodes(1:element_list(lmn)%n_nodes), &       ! Input variables
                 element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &     ! Input variables
                 element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu), iu,      &                              ! Input variables
                 ns, initial_state_variables(iof:iof+ns), &                           ! Input variables
                 updated_state_variables(iof:iof+ns),stif1(1:iu,1:iu),resid1(1:iu), fail)               ! Output variables
          
            numerical_stiffness(1:iu,icount) = -(resid1(1:iu)-resid0(1:iu))/1.D-07

            write(IOW,*)
            write(IOW,*) ' Column ',icount, ' node ',n,' DOF ',i
         
            icount2 = 0
            do j = 1,element_list(lmn)%n_nodes
                do k = 1,node_list(j)%n_dof
                    icount2 = icount2 + 1
                    write(IOW,1000) icount2,j,k,element_stiffness(icount2,icount),numerical_stiffness(icount2,icount)
1000                format( ' Row ',i4,' node ',i4,' DOF ',i4,' Stiffness ',d15.5,' Numerical deriv ',d15.5 )
                end do
            end do
            element_dof_increment(icount) = element_dof_increment(icount) - 1.D-07
        end do
    end do

    err = sum( (numerical_stiffness-element_stiffness)*(numerical_stiffness-element_stiffness) )
     
    if (err>1.d-06*sum(element_stiffness*element_stiffness)) then
        write(IOW,*)
        write(IOW,*) ' Error detected in stiffness matrix '
    else
        write(IOW,'(//A//)') ' Stiffness matrix is consistent with residual '
    endif
	 
    err = 0.D0
    do i = 1,icount
        do j = 1,icount
            err = element_stiffness(i,j)-element_stiffness(j,i)
            if (err*err>1.D-06*(element_stiffness(i,j)+element_stiffness(j,i))**2) then
                write(IOW,*) ' Stiffness is unsymmetric at col, row ',i,j
            endif
            err = numerical_stiffness(i,j)-numerical_stiffness(j,i)
            if (err*err>1.D-06*(numerical_stiffness(i,j)+numerical_stiffness(j,i))**2) then
                write(IOW,*) ' Numerical stiffness is unsymmetric at col, row ',i,j
            endif

        end do
    end do

    deallocate(element_coords)
    deallocate(element_dof_increment)
    deallocate(element_dof_total)
    deallocate(local_nodes)
    deallocate(element_stiffness)
    deallocate(stif1)
    deallocate(numerical_stiffness)
    deallocate(element_residual)
    deallocate(resid0)
    deallocate(resid1)


    return

end subroutine check_stiffness
