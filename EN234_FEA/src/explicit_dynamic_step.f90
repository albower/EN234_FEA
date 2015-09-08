subroutine explicit_dynamic_step
    use Types
    use ParamIO
    use Globals, only : TIME, DTIME
    use Mesh
    use Dynamicstepparameters
    use Printparameters
    implicit none
    
    integer :: status,neq,istep
    
    neq = node_list(n_nodes)%dof_index + node_list(n_nodes)%n_dof - 1

    if (allocated(velocity)) deallocate(velocity)
    if (allocated(acceleration)) deallocate(acceleration)
    if (allocated(lumped_mass)) deallocate(lumped_mass)
    if (allocated(rforce)) deallocate(rforce)
     
    allocate(velocity(neq), stat=status)
    allocate(acceleration(neq), stat=status)
    allocate(lumped_mass(neq), stat=status)
    allocate(rforce(neq), stat=status)
     
    TIME = 0.d0
    DTIME = dynamic_timestep

    write(IOW,'(//A//)') '  Started explicit dynamic analysis '

    !       Mass matrix assumed to be constant - move this into time loop if not

    call assemble_lumped_mass

    velocity(1:neq) = dof_increment(1:neq)
    acceleration(1:neq) = 0.d0
    do istep = 1,no_dynamic_steps


        dof_increment(1:neq) = DTIME*(velocity(1:neq) + 0.5d0*DTIME*acceleration(1:neq))

        rforce = 0.d0
        call apply_dynamic_boundaryconditions
 
        if (stateprint) then
             write(IOW,*) ' '
             write(IOW,*) '  Step ',istep,' completed'
             write(IOW,*) '  Elapsed time is ',TIME
             write(IOW,*) '  Printing state... '
            if (mod(istep,state_print_steps)==0) call print_state
        endif

        if (userprint) then
             write(IOW,*) '  Step ',istep,' completed'
             write(IOW,*) '  User print activated '
            if (mod(istep,user_print_steps)==0) call user_print(istep)
        endif

        call assemble_dynamic_force

        acceleration(1:neq) = rforce(1:neq)/lumped_mass(1:neq)
        velocity(1:neq) = velocity(1:neq) + DTIME*acceleration(1:neq)
        dof_total(1:neq) = dof_total(1:neq) + dof_increment(1:neq)
        initial_state_variables = updated_state_variables               ! Update element state
        TIME = TIME + DTIME

    end do

end subroutine explicit_dynamic_step


subroutine assemble_lumped_mass
    use Types
    use ParamIO
    use User_Subroutine_Storage
    use Mesh
    implicit none

    ! Local Variables
    integer      :: i1, irow, iu, ix, j,k, status
    integer      :: lmn, n
    integer      :: nn, node1
    integer      :: iof, ns
  

    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
                                                            
    real( prec ), allocatable   :: element_mass(:)
  
    type (node), allocatable ::  local_nodes(:)

    !     Subroutine to assemble global stiffness matrix

    allocate(element_coords(length_coord_array), stat = status)
    allocate(element_dof_increment(length_dof_array), stat = status)
    allocate(element_dof_total(length_dof_array), stat = status)
    allocate(local_nodes(length_node_array), stat = status)
    allocate(element_mass(length_dof_array), stat = status)

  
    if (status/=0) then
        write(IOW,*) ' Error in subroutine assemble_lumped_mass_matrix '
        write(IOW,*) ' Unable to allocate memory for user subroutines '
        stop
    endif
 
    lumped_mass = 0.d0
  
  
    do lmn = 1, n_elements
        !     Extract local coords, DOF for the element
        ix = 0
        iu = 0
        do j = 1, element_list(lmn)%n_nodes
            n = connectivity(element_list(lmn)%connect_index + j - 1)
            local_nodes(j)%n_coords = node_list(n)%n_coords
            local_nodes(j)%coord_index = ix+1
            do k = 1, node_list(n)%n_coords
                ix = ix + 1
                element_coords(ix) = coords(node_list(n)%coord_index + k - 1)
            end do
            local_nodes(j)%dof_index = iu+1
            local_nodes(j)%n_dof = node_list(n)%n_dof
            do k = 1, node_list(n)%n_dof
                iu = iu + 1
                element_dof_increment(iu) = dof_increment(node_list(n)%dof_index + k - 1)
                element_dof_total(iu) = dof_total(node_list(n)%dof_index + k - 1)
            end do
        end do

        !     Form element lumped mass
        iof = element_list(lmn)%state_index
        if (iof==0) iof = 1
        ns = element_list(lmn)%n_states
        if (ns==0) ns=1

        call user_element_lumped_mass(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes,&                           ! Input variables
            local_nodes(1:element_list(lmn)%n_nodes), &                                                                  ! Input variables
            densities(element_list(lmn)%density_index),   &                                                              ! Input variables
            element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &     ! Input variables
            element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu),iu,      &                     ! Input variables
            ns, initial_state_variables(iof:iof+ns-1), &                                                                ! Input variables
            updated_state_variables(iof:iof+ns),element_mass(1:iu))                 ! Output variables

        !     --   Add element lumped mass matrix to global array
        irow = 1
        do i1 = 1, element_list(lmn)%n_nodes
            node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
            iof = node_list(node1)%dof_index
            nn= node_list(node1)%n_dof
            lumped_mass(iof:iof+nn-1) = lumped_mass(iof:iof+nn-1) + element_mass(irow:irow+nn-1)
            irow = irow + nn
        end do
    end do

    return
  
end subroutine assemble_lumped_mass
  
subroutine assemble_dynamic_force
    use Types
    use ParamIO
    use User_Subroutine_Storage
    use Mesh
    use Stiffness
    implicit none

 
    integer      :: status
    integer      :: lmn
    integer      :: ix,iu,j,n,k,iof,ns
    integer      :: irow,i1,node1
  
  
    integer      :: nn

    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
                                                            
    real( prec ), allocatable   :: element_stiffness(:,:)
    real( prec ), allocatable   :: element_residual(:)
  
    type (node), allocatable ::  local_nodes(:)

    !     Subroutine to assemble global stiffness matrix

    allocate(element_coords(length_coord_array), stat = status)
    allocate(element_dof_increment(length_dof_array), stat = status)
    allocate(element_dof_total(length_dof_array), stat = status)
    allocate(local_nodes(length_node_array), stat = status)
    allocate(element_stiffness(length_dof_array,length_dof_array), stat = status)
    allocate(element_residual(length_dof_array), stat = status)
  
    if (status/=0) then
        write(IOW,*) ' Error in subroutine assemble_dynamic_force '
        write(IOW,*) ' Unable to allocate memory for user subroutines '
        stop
    endif


    do lmn = 1, n_elements
        if (element_deleted(lmn)) cycle

          !     Extract local coords, DOF for the element
        ix = 0
        iu = 0
        do j = 1, element_list(lmn)%n_nodes
            n = connectivity(element_list(lmn)%connect_index + j - 1)
            local_nodes(j)%n_coords = node_list(n)%n_coords
            local_nodes(j)%coord_index = ix+1
            do k = 1, node_list(n)%n_coords
                ix = ix + 1
                element_coords(ix) = coords(node_list(n)%coord_index + k - 1)
            end do
            local_nodes(j)%n_dof = node_list(n)%n_dof
            local_nodes(j)%dof_index = iu+1
            do k = 1, node_list(n)%n_dof
                iu = iu + 1
                element_dof_increment(iu) = dof_increment(node_list(n)%dof_index + k - 1)
                element_dof_total(iu) = dof_total(node_list(n)%dof_index + k - 1)
            end do
        end do
        !     Form element contribution to nodal force vector
        iof = element_list(lmn)%state_index
        if (iof==0) iof = 1
        ns = element_list(lmn)%n_states
        if (ns==0) ns=1

        call user_element_dynamic(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes, &                            ! Input variables
            local_nodes(1:element_list(lmn)%n_nodes), &                                                                ! Input variables
            element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &   ! Input variables
            element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu),iu,      &                   ! Input variables
            ns, initial_state_variables(iof:iof+ns-1), &                                                               ! Input variables
            updated_state_variables(iof:iof+ns),element_residual(1:iu),element_deleted(lmn))                           ! Output variables
      

        !     --   Add element force vector to global array
        irow = 1
        do i1 = 1, element_list(lmn)%n_nodes
            node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
            iof = node_list(node1)%dof_index
            nn = node_list(node1)%n_dof
            rforce(iof:iof+nn-1) = rforce(iof:iof+nn-1) + element_residual(irow:irow+nn-1)
            irow = irow + nn
        end do
    end do

  
    deallocate(element_coords)
    deallocate(element_dof_increment)
    deallocate(element_dof_total)
    deallocate(local_nodes)
    deallocate(element_stiffness)
    deallocate(element_residual)
  

end subroutine assemble_dynamic_force
  
  
  
subroutine apply_dynamic_boundaryconditions
    use Types
    use ParamIO
    use User_Subroutine_Storage
    use Mesh
    use Boundaryconditions
    use Element_Utilities, only : facenodes
    use Globals, only : TIME, DTIME
    implicit none


    ! Local Variables
    logical :: ignoredof
    real( prec ) :: ucur, dofvalue, dloadvalue, force_value
    integer :: idof, ix,iu, i,j,k, lmn, n, iof, iofs,nhist, nparam, nnodes, status
    integer :: load,flag,elset,ifac,nel,param_index,ntract,ndims,ndof,nfacenodes
    integer :: i1,node1,irow,nn
    integer :: list(8)
       
       
    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
                                                            
    real( prec ), allocatable   :: element_stiffness(:,:)
    real( prec ), allocatable   :: element_residual(:)
  
    real( prec ), allocatable   :: traction(:)
  


    type (node), allocatable ::  local_nodes(:)

    !     Subroutine to apply boundary conditions

    if (n_distributedloads>0) then
        allocate(element_coords(length_coord_array), stat = status)
        allocate(element_dof_increment(length_dof_array), stat = status)
        allocate(element_dof_total(length_dof_array), stat = status)
        allocate(local_nodes(length_node_array), stat = status)
        allocate(element_stiffness(length_dof_array,length_dof_array), stat = status)
        allocate(element_residual(length_dof_array), stat = status)
        allocate(traction(length_dof_array), stat = status)
    
        if (status /=0) then
            write(IOW,*) ' Error in subroutine apply_direct_boundaryconditions '
            write(IOW,*) ' Unable to allocate memory for distributed forces '
            stop
        endif
    endif


     !  -- Distributed forces on element faces
   
    do load = 1, n_distributedloads
        flag = distributedload_list(load)%flag
        elset = distributedload_list(load)%element_set
        ifac = distributedload_list(load)%face
        nel = elementset_list(elset)%n_elements
        iof = elementset_list(elset)%index
    
        if (flag<4) then
            do k  = 1, nel
                lmn = element_lists(iof+k-1)
                ndims = node_list(connectivity(element_list(lmn)%connect_index))%n_coords
                ndof = node_list(connectivity(element_list(lmn)%connect_index))%n_dof
                traction = 0.d0
                if (flag<3) then
                    ntract = distributedload_list(load)%n_dload_values
                    do i = 1,ntract
                        traction(i) = dload_values(distributedload_list(load)%index_dload_values+i-1)
                    end do
                else
                    ntract = 1
                    traction(1) = 1.d0
                endif
                if (flag==2) traction = traction/dsqrt(dot_product(traction,traction))
                if (flag>1) then
                    iof = history_list(distributedload_list(load)%history_number)%index
                    nhist = history_list(distributedload_list(load)%history_number)%n_timevals
                    call interpolate_history_table(history_data(1,iof),nhist,TIME+DTIME,dloadvalue)
                    traction = traction*dloadvalue
                endif

                call facenodes(ndims,element_list(lmn)%n_nodes,ifac,list,nfacenodes)
                ix = 0
                iu = 0
                do j = 1,nfacenodes
                    n = connectivity(element_list(lmn)%connect_index + list(j) - 1)
                    do i = 1,ndims
                        ix = ix + 1
                        element_coords(ix) = coords(node_list(n)%coord_index+i-1)
                    end do
                    do i = 1, node_list(n)%n_dof
                        iu = iu + 1
                        element_dof_increment(iu) = dof_increment(node_list(n)%dof_index + i - 1)
                        element_dof_total(iu) = dof_total(node_list(n)%dof_index + i - 1)
                    end do
                end do
  
                call traction_boundarycondition_dynamic(flag,ndims,ndof,nfacenodes,element_coords(1:ix),ix,&
                    element_dof_increment(1:iu),element_dof_total(1:iu),traction(1:ntract),ntract,&
                    element_residual(1:iu))               ! Output variables


                irow = 1
                do i1 = 1,nfacenodes
                    node1 = connectivity(element_list(lmn)%connect_index + list(i1) - 1)
                    iof = node_list(node1)%dof_index
                    nn = node_list(node1)%n_dof
                    rforce(iof:iof+nn-1) = rforce(iof:iof+nn-1) + element_residual(irow:irow+nn-1)
                    irow = irow + nn
                end do
            end do
    
        else if (flag==4) then              ! User subroutine controlled distributed load

            param_index = subroutineparameter_list(distributedload_list(load)%subroutine_parameter_number)%index
            nparam = subroutineparameter_list(distributedload_list(load)%subroutine_parameter_number)%index
            if (param_index==0) param_index=1
            if (nparam==0) nparam = 1
            do k  = 1, nel
                lmn = element_lists(iof+k-1)
                !     Extract local coords, DOF for the element
                ix = 0
                iu = 0
                do j = 1, element_list(lmn)%n_nodes
                    n = connectivity(element_list(lmn)%connect_index + j - 1)
                    local_nodes(j)%n_coords = node_list(n)%n_coords
                    local_nodes(j)%coord_index = ix+1
                    do i = 1, node_list(n)%n_coords
                        ix = ix + 1
                        element_coords(ix) = coords(node_list(n)%coord_index + i - 1)
                    end do
                    local_nodes(j)%dof_index = iu+1
                    do i = 1, node_list(n)%n_dof
                        iu = iu + 1
                        element_dof_increment(iu) = dof_increment(node_list(n)%dof_index + i - 1)
                        element_dof_total(iu) = dof_total(node_list(n)%dof_index + i - 1)
                    end do
                end do
     
        
                call user_distributed_load_dynamic(lmn, element_list(lmn)%flag, ifac, &      ! Input variables
                    subroutine_parameters(param_index:param_index+nparam-1),nparam, &       ! Input variables
                    element_list(lmn)%n_nodes, local_nodes(1:element_list(lmn)%n_nodes), &       ! Input variables
                    element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &     ! Input variables
                    element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu),iu,      &                              ! Input variables
                    element_residual(1:iu))               ! Output variables

                !     --   Add element stiffness and residual to global array

                irow = 1
                do i1 = 1, element_list(lmn)%n_nodes
                    node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
                    iof = node_list(node1)%dof_index
                    nn = node_list(node1)%n_dof
                    rforce(iof:iof+nn-1) = rforce(iof:iof+nn-1) + element_residual(irow:irow+nn-1)
                    irow = irow + nn
                end do
            end do
        endif
    end do

    !     -- Prescribed nodal forces
    do k = 1, n_prescribedforces
        if (prescribedforce_list(k)%flag==1) then                                          ! Prescribe dof directly
            force_value = dof_values(prescribedforce_list(k)%index_dof_values)
        else if (prescribedforce_list(k)%flag==2) then                                     ! Interpolate a history table
            iof = history_list(prescribedforce_list(k)%history_number)%index
            nhist = history_list(prescribedforce_list(k)%history_number)%n_timevals
            call interpolate_history_table(history_data(1,iof),nhist,TIME+DTIME,force_value)
        endif
       
        if (prescribedforce_list(k)%node_set==0) then    ! Constrain a single node
            n = prescribedforce_list(k)%node_number
            idof = prescribedforce_list(k)%dof
 
            if (prescribedforce_list(k)%flag==3) then
                iof = subroutineparameter_list(prescribedforce_list(k)%subroutine_parameter_number)%index
                nparam = subroutineparameter_list(prescribedforce_list(k)%subroutine_parameter_number)%n_parameters
                call user_prescribedforce(n,idof,subroutine_parameters(iof),nparam,force_value)
            endif
            irow = idof + node_list(n)%dof_index - 1
            rforce(irow) = rforce(irow) + force_value
        else  ! Constrain a node set
            nnodes = nodeset_list(prescribedforce_list(k)%node_set)%n_nodes
            iof = nodeset_list(prescribedforce_list(k)%node_set)%index
            do i = 1,nnodes
                n = node_lists(iof+i-1)
                idof = prescribedforce_list(k)%dof
                if (prescribedforce_list(k)%flag==3) then
                    iof = subroutineparameter_list(prescribedforce_list(k)%subroutine_parameter_number)%index
                    nparam = subroutineparameter_list(prescribedforce_list(k)%subroutine_parameter_number)%n_parameters
                    call user_prescribedforce(n,idof,subroutine_parameters(iof),nparam,force_value)
                endif
                irow = idof + node_list(n)%dof_index - 1
                rforce(irow) = rforce(irow) + force_value
            end do
        endif
   
    end do


    !     -- Prescribed DOFs
    do k = 1, n_prescribeddof
        ignoredof = .false.
        if (prescribeddof_list(k)%flag==1) then                                          ! Prescribe dof directly
            dofvalue = dof_values(prescribeddof_list(k)%index_dof_values)
        else if (prescribeddof_list(k)%flag==2) then                                     ! Interpolate a history table
            iof = history_list(prescribeddof_list(k)%history_number)%index
            nhist = history_list(prescribeddof_list(k)%history_number)%n_timevals
            call interpolate_history_table(history_data(1,iof),nhist,TIME+DTIME,dofvalue)
        endif
       
        if (prescribeddof_list(k)%node_set==0) then    ! Constrain a single node
            n = prescribeddof_list(k)%node_number
            idof = prescribeddof_list(k)%dof
            if (prescribeddof_list(k)%rate_flag==0) then
                ucur =  dof_total(idof + node_list(n)%dof_index - 1)
            else
                ucur = 0.d0
            endif
            if (prescribeddof_list(k)%flag==3) then
                iofs = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%index
                nparam = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%n_parameters
                call user_prescribeddof(n,idof,subroutine_parameters(iofs),nparam,dofvalue,ignoredof)
            endif
            if (ignoredof) cycle
            dof_increment(idof + node_list(n)%dof_index - 1) = dofvalue - ucur
        else  ! Constrain a node set
            nnodes = nodeset_list(prescribeddof_list(k)%node_set)%n_nodes
            iof = nodeset_list(prescribeddof_list(k)%node_set)%index
            do i = 1,nnodes
                ignoredof = .false.
                n = node_lists(iof+i-1)
                idof = prescribeddof_list(k)%dof
                if (prescribeddof_list(k)%rate_flag==0) then
                    ucur =  dof_total(idof + node_list(n)%dof_index - 1)
                else
                    ucur = 0.d0
                endif
                if (prescribeddof_list(k)%flag==3) then
                    iofs = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%index
                    nparam = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%n_parameters
                    call user_prescribeddof(n,idof,subroutine_parameters(iofs),nparam,dofvalue,ignoredof)
                endif
                if (ignoredof) cycle
                dof_increment(idof + node_list(n)%dof_index - 1) = dofvalue - ucur
            end do
        endif
   
    end do

    if (n_distributedloads>0) then
       deallocate(element_coords)
       deallocate(element_dof_increment)
       deallocate(element_dof_total)
       deallocate(local_nodes)
       deallocate(element_stiffness)
       deallocate(element_residual)
       deallocate(traction)
    endif

end subroutine apply_dynamic_boundaryconditions
