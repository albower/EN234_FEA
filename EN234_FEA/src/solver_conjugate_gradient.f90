subroutine assemble_conjugate_gradient_stiffness(fail)
    use Types
    use ParamIO
    use Globals
    use User_Subroutine_Storage
    use Mesh
    use Boundaryconditions, only: nodeset,nodeset_list, node_lists
    use Boundaryconditions, only: n_constraints,constraint,constraint_list
    use Boundaryconditions, only: constraint_parameters, constraintparameter_list,lagrange_multipliers
    use Staticstepparameters, only : nonlinear,current_step_number,max_total_time,abq_PNEWDT
    use Controlparameters, only : abaqusformat
    use Stiffness
    implicit none

    logical, intent(out) :: fail

    ! Local Variables

    real( prec ) :: diagnorm                         ! measure of average of diagonal in stiffness
   
    integer      :: irow,icol,index
    integer      :: status
    integer      :: lmn
    integer      :: ix,iu,j,n,k,iof,ns
    integer      :: nsr,i1,node1,j1,nsc,i2,node2,j2
    integer      :: nc,nset,nnodes
    integer      :: dof1,iof1,dof2,iof2
    integer      :: ipar,npar,idof
    real (prec)  :: lmult
    integer      :: abq_JTYPE
    integer      :: abq_MDLOAD
    integer      :: abq_NPREDF
    integer      :: abq_LFLAGS(5)

    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
    real( prec ), allocatable    :: element_state_variables(:)

    real( prec ), allocatable    :: abq_PREDEF(:,:,:)
    real( prec ), allocatable    :: abq_ADLMAG(:)
    real( prec ), allocatable    :: abq_DDLMAG(:)
    real( prec ), allocatable    :: abq_V(:)
    real( prec ), allocatable    :: abq_A(:)

    real( prec ), allocatable   :: element_stiffness(:,:)
    real( prec ), allocatable   :: element_residual(:)

    integer, allocatable :: abq_JDLTYP(:)
  
    real( prec ) :: abq_time(2)      ! Time variable for ABAQUS UEL interface
    real( prec ) :: abq_PARAMS(3)    ! Parameter variable for ABAQUS UEL interface
    real( prec ) :: element_PNEWDT   ! ABAQUS time increment control parameter

    type (node), allocatable ::  local_nodes(:)


    !     Subroutine to assemble global stiffness matrix

    if (abaqusformat) then
        call generate_abaqus_dloads           ! Generate_abaqus_dloads is in solver_direct.f90
    endif

    j = length_coord_array
    if (abaqusformat) j = max(length_coord_array,length_dof_array)

    allocate(element_coords(j), stat = status)
    allocate(element_dof_increment(length_dof_array), stat = status)
    allocate(element_dof_total(length_dof_array), stat = status)
    allocate(local_nodes(length_node_array), stat = status)
    allocate(element_stiffness(length_dof_array,length_dof_array), stat = status)
    allocate(element_residual(length_dof_array), stat = status)
    allocate(element_state_variables(length_state_variable_array), stat=status)
  
    if (status/=0) then
        write(IOW,*) ' Error in subroutine assemble_direct_stiffness '
        write(IOW,*) ' Unable to allocate memory for user subroutines '
        stop
    endif

    if (abaqusformat) then
        allocate(abq_PREDEF(2,1,length_node_array), stat = status)
        allocate(abq_V(length_dof_array), stat=status)
        allocate(abq_A(length_dof_array), stat=status)
    endif

    if (length_abq_dlmag_array>0) then
        allocate(abq_ADLMAG(length_abq_dlmag_array), stat = status)
        allocate(abq_DDLMAG(length_abq_dlmag_array), stat = status)
        allocate(abq_JDLTYP(length_abq_dlmag_array), stat = status)
    endif
    if (status/=0) then
        write(IOW,*) ' Error in subroutine assemble_direct_stiffness '
        write(IOW,*) ' Unable to allocate memory for user subroutines '
        stop
    endif

    element_stiffness = 0.d0
    element_residual = 0.d0

    aupp = 0.d0             ! Upper triangle of stiffness
    if (unsymmetric_stiffness) alow = 0.d0    ! Lower triangle of stiffness
    rhs = 0.d0              ! Right hand side
    diag = 0.d0             ! Stiffness diagonal

    last_filled_equation_adjacency = 0               ! Index of last filled databin in equation adjacency linked list
    last_filled_aupp = 0                             ! Index of last filled entry in aupp

    equation_adjacency_index = 0
    equation_adjacency(1:length_equation_adjacency)%n_entries = 0
    equation_adjacency(1:length_equation_adjacency)%next = 0
 
    abq_PNEWDT = 1.d0               ! ABAQUS timestep cutback deactivated

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
        !     Form element stiffness
        iof = element_list(lmn)%state_index
        if (iof==0) iof = 1
        ns = element_list(lmn)%n_states
        if (ns==0) ns=1

        if (element_list(lmn)%flag>99999) then                ! ABAQUS UEL format user subroutine
            !
            abq_V = 0.d0
            abq_A = 0.d0
            abq_time(1:2) = TIME
            abq_JTYPE = element_list(lmn)%flag-99999

            abq_LFLAGS = 0
            if (nonlinear) abq_LFLAGS(2) = 1
            abq_LFLAGS(3) = 1

            !           Change storage of element coords to ABAQUS UEL format
            if (abq_MCRD(lmn) == 0) then
                do j = 1,element_list(lmn)%n_nodes
                    if (local_nodes(j)%n_coords >abq_MCRD(lmn)) then
                        abq_MCRD(lmn) = local_nodes(j)%n_coords
                    endif
                    if (local_nodes(j)%n_dof>abq_MCRD(lmn)) then
                        abq_MCRD(lmn) = node_list(n)%n_dof
                        if (abq_MCRD(lmn)>3) abq_MCRD(lmn) = 3
                    endif
                end do
            endif
            ix = 0
            do j = 1, element_list(lmn)%n_nodes
                n = connectivity(element_list(lmn)%connect_index + j - 1)
                do k = 1, node_list(n)%n_coords
                    ix = ix + 1
                    element_coords(ix) = coords(node_list(n)%coord_index + k - 1)
                end do
                if (node_list(n)%n_coords<abq_MCRD(lmn)) ix = ix + abq_MCRD(lmn)-node_list(n)%n_coords
            end do

            abq_PARAMS = 0

            abq_MDLOAD = abq_uel_bc_list(lmn)%mdload
            if (abq_MDLOAD>0) then
                abq_ADLMAG(1:abq_MDLOAD) = &
                    abq_uel_bc_mag(abq_uel_bc_list(lmn)%mag_index:abq_uel_bc_list(lmn)%mag_index+abq_MDLOAD-1)
                abq_JDLTYP(1:abq_MDLOAD) = &
                    abq_uel_bc_typ(abq_uel_bc_list(lmn)%mag_index:abq_uel_bc_list(lmn)%mag_index+abq_MDLOAD-1)
                abq_DDLMAG(1:abq_MDLOAD) = &
                    abq_uel_bc_dmag(abq_uel_bc_list(lmn)%mag_index:abq_uel_bc_list(lmn)%mag_index+abq_MDLOAD-1)
            endif

            abq_PREDEF(1,1,1:element_list(lmn)%n_nodes) = BTEMP+BTINC
            abq_PREDEF(2,1,1:element_list(lmn)%n_nodes) = BTINC
            abq_NPREDF = 1

            element_PNEWDT = 100.d0
            fail = .false.

            if (element_list(lmn)%n_states==0) then 
               ns=1
               element_state_variables(1) = 0.d0
            else
               element_state_variables(1:ns) = initial_state_variables(iof:iof+ns-1)
            endif    

            call UEL(element_residual(1:iu),element_stiffness(1:iu,1:iu),element_state_variables(1:ns), &
                energy(8*lmn-7:8*lmn),iu,1,ns, &
                element_properties(element_list(lmn)%element_property_index),element_list(lmn)%n_element_properties, &
                element_coords(1:ix),abq_MCRD(lmn),element_list(lmn)%n_nodes,element_dof_increment(1:iu)+element_dof_total(1:iu), &
                element_dof_increment(1:iu),abq_V(1:iu),abq_A(1:iu),abq_JTYPE,abq_time,DTIME, &
                1,current_step_number,lmn,abq_PARAMS,abq_MDLOAD,abq_JDLTYP,abq_ADLMAG, &
                abq_PREDEF(1,1,1:element_list(lmn)%n_nodes),abq_NPREDF, &
                abq_LFLAGS,iu,abq_DDLMAG,abq_MDLOAD,element_PNEWDT, &
                int_element_properties(element_list(lmn)%int_element_property_index),element_list(lmn)%n_int_element_properties, &
                max_total_time)

            if (element_PNEWDT < abq_PNEWDT) then
                abq_PNEWDT = element_PNEWDT
            endif

            if (element_list(lmn)%n_states>0) updated_state_variables(iof:iof+ns-1) = element_state_variables(1:ns)

        else

            call user_element_static(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes, &                               ! Input variables
                local_nodes(1:element_list(lmn)%n_nodes), &                                                                 ! Input variables
                element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &     ! Input variables
                element_list(lmn)%n_int_element_properties,int_element_properties(element_list(lmn)%int_element_property_index), & ! Input variables
                element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu),iu,      &                     ! Input variables
                ns, initial_state_variables(iof:iof+ns-1), &                                                                 ! Input variables
                updated_state_variables(iof:iof+ns),element_stiffness(1:iu,1:iu),element_residual(1:iu), fail)      ! Output variables
      

            if (fail) return
        endif
        !     --   Add element stiffness and residual to global array

        nsr = 0
        do i1 = 1, element_list(lmn)%n_nodes
            node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
            do j1 = 1, node_list(node1)%n_dof
                irow = ieqs(j1 + node_list(node1)%dof_index - 1)
                nsr = nsr + 1
                rhs(irow) = rhs(irow) + element_residual(nsr)
                nsc = 0
                do i2 = 1, element_list(lmn)%n_nodes
                    node2 = connectivity(i2 + element_list(lmn)%connect_index - 1)
                    do j2 = 1, node_list(node2)%n_dof
                        icol = ieqs(j2 + node_list(node2)%dof_index - 1)
                        nsc = nsc + 1
                        if ( icol==irow ) then
                            diag(irow) = diag(irow) + element_stiffness(nsr,nsc)
                        else if ( icol>irow ) then
                            if (dabs(element_stiffness(nsr,nsc))+dabs(element_stiffness(nsc,nsr))>thresh) then
                                call getindex(irow,icol,index)
                                aupp(index) = aupp(index) + element_stiffness(nsr,nsc)
                                if ( unsymmetric_stiffness ) alow(index) = alow(index) + element_stiffness(nsc,nsr)
                            endif
                        end if
                    end do
                end do
            end do
        end do
    end do

    diagnorm = sqrt(dot_product(diag(1:neq),diag(1:neq)))/neq
  
    do nc = 1, n_constraints
  
        if (constraint_list(nc)%flag < 3) then          ! Simple two-node tie constraint
            icol = ieqs(length_dofs+nc)
            node1 = constraint_list(nc)%node1
            dof1 = constraint_list(nc)%dof1
            iof1 = dof1 + node_list(node1)%dof_index - 1
            node2 = constraint_list(nc)%node2
            dof2 = constraint_list(nc)%dof2
            iof2 = dof2 + node_list(node2)%dof_index - 1
        
            irow = ieqs(iof1)
            rhs(irow) = rhs(irow) - lagrange_multipliers(nc)
            if (icol>irow) then
                call getindex(irow,icol,index)
            else
                call getindex(icol,irow,index)
            endif
            aupp(index) = aupp(index) + 1.d0
            if (unsymmetric_stiffness) alow(index) = alow(index) + 1.d0

            irow = ieqs(iof2)
            if (icol>irow) then
                call getindex(irow,icol,index)
            else
                call getindex(icol,irow,index)
            endif
            rhs(irow) = rhs(irow) + lagrange_multipliers(nc)
            aupp(index) = aupp(index) - 1.d0
            if (unsymmetric_stiffness) alow(index) = alow(index) - 1.d0

            rhs(icol) = dof_total(iof2)+dof_increment(iof2)-dof_total(iof1)-&
                dof_increment(iof1)-lagrange_multipliers(nc)*1.d-12*diagnorm
            diag(icol) = 1.d-12*diagnorm
        else                                                                    ! General multi-node constraint
            nset = constraint_list(nc)%node1
            nnodes = nodeset_list(nset)%n_nodes
            ix = 0
            iu = 0
            do j = 1, nnodes
                n = node_lists(nodeset_list(nset)%index + j - 1)
                local_nodes(j)%n_coords = node_list(n)%n_coords
                local_nodes(j)%coord_index = ix+1
                do k = 1, node_list(n)%n_coords
                    ix = ix + 1
                    element_coords(ix) = coords(node_list(n)%coord_index + k - 1)
                end do
                local_nodes(j)%dof_index = iu+1
                do k = 1, node_list(n)%n_dof
                    iu = iu + 1
                    element_dof_increment(iu) = dof_increment(node_list(n)%dof_index + k - 1)
                    element_dof_total(iu) = dof_total(node_list(n)%dof_index + k - 1)
                end do
            end do

            ipar = constraintparameter_list(constraint_list(nc)%index_parameters)%index
            npar = constraintparameter_list(constraint_list(nc)%index_parameters)%n_parameters
            idof = constraint_list(nc)%node2                                 ! Index of node set containing DOF list
            lmult = lagrange_multipliers(nc)
            !
            call user_constraint(nc, constraint_list(nc)%flag, nodeset_list(nset)%n_nodes,&    ! Input variables
                local_nodes(1:nodeset_list(nset)%n_nodes), node_lists(idof:idof+nnodes-1),&    ! Input variables
                npar, constraint_parameters(ipar:ipar+npar-1), &                               ! Input variables
                element_coords(1:ix),ix, &                                                           ! Input variables
                element_dof_increment(1:iu), element_dof_total(1:iu),iu, &                                 ! Input variables
                lmult,  &                                                                      ! Input variables
                element_stiffness(1:iu,1:iu),element_residual(1:iu))      ! Output variables

            nsr = 0
            do i1 = 1, nnodes
                node1 = node_lists(nodeset_list(nset)%index + i1 - 1)
                j1 = node_lists(idof+i1-1)
                irow = ieqs(j1 + node_list(node1)%dof_index - 1)
                nsr = nsr + 1
                rhs(irow) = rhs(irow) + element_residual(nsr)
                nsc = 0
                do i2 = 1, nnodes
                    node2 = node_lists(nodeset_list(nset)%index + i2 - 1)
                    j2 = node_lists(idof+i2-1)
                    icol = ieqs(j2 + node_list(node2)%dof_index - 1)
                    nsc = nsc + 1
                    if ( icol==irow ) then
                        diag(irow) = diag(irow) + element_stiffness(nsr,nsc)
                    else if ( icol>irow ) then
                        if (dabs(element_stiffness(nsr,nsc))+dabs(element_stiffness(nsc,nsr))>thresh) then
                            call getindex(irow,icol,index)
                            aupp(index) = aupp(index) + element_stiffness(nsr,nsc)
                            if ( unsymmetric_stiffness ) alow(index) = alow(index) + element_stiffness(nsc,nsr)
                        endif
                    end if
                end do
            end do
            irow = ieqs(length_dofs+nc)
            nsr = nsr + 1
            rhs(irow) = rhs(irow) + element_residual(nsr)
            do i2 = 1, nnodes
                node2 = node_lists(nodeset_list(nset)%index + i2 - 1)
                j2 = node_lists(idof+i2-1)
                icol = ieqs(j2 + node_list(node2)%dof_index - 1)
                nsc = nsc + 1
                if ( icol==irow ) then
                    diag(irow) = diag(irow) + element_stiffness(nsr,nsc)
                else if ( icol>irow ) then
                    if (dabs(element_stiffness(nsr,nsc))+dabs(element_stiffness(nsc,nsr))>thresh) then
                        call getindex(irow,icol,index)
                        aupp(index) = aupp(index) + element_stiffness(nsr,nsc)
                        if ( unsymmetric_stiffness ) alow(index) = alow(index) + element_stiffness(nsc,nsr)
                    endif
                end if
            end do
            diag(irow) = diag(irow) + element_stiffness(nsr,nsr)
        endif
    end do

    !  Force a cutback if NaN is found in RHS.
    do n = 1, neq
        if (rhs(n)==rhs(n)+1.d0) then
            write(IOW,*) ' NaN detected in residual '
            flush(IOW)
            fail = .true.
            return
        endif
        if (isnan( rhs(n) ))  then
            write(IOW,*) ' NaN detected in residual '
            flush(IOW)
            fail = .true.
            return
        endif
    end do
    ! Store generalized nodal forces
    do n = 1, n_nodes
        do j = 1, node_list(n)%n_dof
            rforce(j + node_list(n)%dof_index - 1) = -rhs(ieqs(j + node_list(n)%dof_index - 1))
        end do
    end do
  
    nodal_force_norm = dsqrt(dot_product(rhs,rhs))
  
    deallocate(element_coords)
    deallocate(element_dof_increment)
    deallocate(element_dof_total)
    deallocate(element_stiffness)
    deallocate(element_residual)

    if (allocated(abq_PREDEF)) deallocate(abq_PREDEF)
    if (allocated(abq_ADLMAG)) deallocate(abq_ADLMAG)
    if (allocated(abq_JDLTYP)) deallocate(abq_JDLTYP)
    if (allocated(abq_uel_bc_typ)) deallocate(abq_uel_bc_typ)
    if (allocated(abq_uel_bc_mag)) deallocate(abq_uel_bc_mag)
    if (allocated(abq_uel_bc_dmag)) deallocate(abq_uel_bc_dmag)

end subroutine assemble_conjugate_gradient_stiffness

!====================Subroutine BCONS_CG =======================
subroutine apply_cojugategradient_boundaryconditions

    use Types
    use ParamIO
    use User_Subroutine_Storage
    use Globals, only: TIME,DTIME
    use Mesh
    use Boundaryconditions
    use Stiffness
    use Element_Utilities, only : facenodes
    use Linkedlist_handling
    implicit none

    ! Local Variables
    logical :: ignoredof
    real( prec ) :: ucur, dofvalue, dofvalue_correction, dloadvalue
    
    integer :: idof, ix,iu, i,j,k, lmn, n, iof, iofs, nhist, nparam, nnodes, status
    integer :: load,flag,elset,ifac,nel,param_index,ntract,ndims,ndof,nfacenodes
    integer :: nsr,i1,node1,j1,irow,nsc,i2,node2,j2,icol, index
    integer :: list(8)
       
       
    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
                                                            
    real( prec ), allocatable   :: element_stiffness(:,:)
    real( prec ), allocatable   :: element_residual(:)
  
    real( prec ), allocatable   :: traction(:)

    type (node), allocatable ::  local_nodes(:)


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
                if (element_list(lmn)%flag>99999) cycle ! Skip boundary conditions on  ABAQUS UEL
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
  
                call traction_boundarycondition_static(flag,ndims,ndof,nfacenodes,element_coords(1:ix),ix,&
                    element_dof_increment(1:iu),element_dof_total(1:iu),iu,traction(1:ntract),ntract,&
                    element_stiffness(1:iu,1:iu),element_residual(1:iu))               ! Output variables

                nsr = 0
                do i1 = 1,nfacenodes
                    node1 = connectivity(element_list(lmn)%connect_index + list(i1) - 1)
                    do j1 = 1, node_list(node1)%n_dof
                        irow = ieqs(j1 + node_list(node1)%dof_index - 1)
                        nsr = nsr + 1
                        rhs(irow) = rhs(irow) + element_residual(nsr)
                        nsc = 0
                        do i2 = 1,nfacenodes
                            node2 = connectivity(element_list(lmn)%connect_index+ list(i2) - 1)
                            do j2 = 1, node_list(node2)%n_dof
                                icol = ieqs(j2 + node_list(node2)%dof_index - 1)
                                nsc = nsc + 1
                                if ( icol==irow ) then
                                    diag(irow) = diag(irow) + element_stiffness(nsr,nsc)
                                else if ( icol>irow ) then
                                    if (dabs(element_stiffness(nsc,nsr))+dabs(element_stiffness(nsr,nsc))>thresh) then
                                        call getindex(irow,icol,index)
                                        aupp(index) = aupp(index) + element_stiffness(nsr,nsc)
                                        if ( unsymmetric_stiffness ) alow(index) = alow(index) + element_stiffness(nsc,nsr)
                                    endif
                                end if
                            end do
                        end do
                    end do
                end do
            end do
    
        else if (flag==4) then              ! User subroutine controlled distributed load


            do k  = 1, nel
                lmn = element_lists(iof+k-1)
                if (element_list(lmn)%flag>99999) cycle ! Skip boundary conditions on  ABAQUS UEL
                param_index = subroutineparameter_list(distributedload_list(load)%subroutine_parameter_number)%index
                nparam = subroutineparameter_list(distributedload_list(load)%subroutine_parameter_number)%index
                if (param_index==0) param_index=1
                if (nparam==0) nparam = 1
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
        
                call user_distributed_load_static(lmn, element_list(lmn)%flag, ifac, &      ! Input variables
                    subroutine_parameters(param_index:param_index+nparam-1),nparam, &       ! Input variables
                    element_list(lmn)%n_nodes, local_nodes(1: element_list(lmn)%n_nodes), &       ! Input variables
                    element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &     ! Input variables
                    element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu), iu,  &   ! Input variables
                    element_stiffness(1:iu,1:iu),element_residual(1:iu) )              ! Output variables

                !     --   Add element stiffness and residual to global array

                nsr = 0
                do i1 = 1, element_list(lmn)%n_nodes
                    node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
                    do j1 = 1, node_list(node1)%n_dof
                        irow = ieqs(j1 + node_list(node1)%dof_index - 1)
                        nsr = nsr + 1
                        rhs(irow) = rhs(irow) + element_residual(nsr)
                        nsc = 0
                        do i2 = 1, element_list(lmn)%n_nodes
                            node2 = connectivity(i2 + element_list(lmn)%connect_index - 1)
                            do j2 = 1, node_list(node2)%n_dof
                                icol = ieqs(j2 + node_list(node2)%dof_index - 1)
                                nsc = nsc + 1
                                if ( icol==irow ) then
                                    diag(irow) = diag(irow) + element_stiffness(nsr,nsc)
                                else if ( icol>irow ) then
                                    if (dabs(element_stiffness(nsr,nsc))+dabs(element_stiffness(nsc,nsr))>thresh) then
                                        call getindex(irow,icol,index)
                                        aupp(index) = aupp(index) + element_stiffness(nsr,nsc)
                                        if ( unsymmetric_stiffness ) alow(index) = alow(index) + element_stiffness(nsc,nsr)
                                    endif
                                end if
                            end do
                        end do
                    end do
                end do
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
                ucur = dof_increment(idof + node_list(n)%dof_index - 1) + dof_total(idof + node_list(n)%dof_index - 1)
            else
                ucur = dof_increment(idof + node_list(n)%dof_index - 1)
            endif
            if (prescribeddof_list(k)%flag==3) then
                iofs = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%index
                nparam = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%n_parameters
                call user_prescribeddof(n,idof,subroutine_parameters(iofs),nparam,dofvalue,ignoredof)
            endif
            if (ignoredof) cycle
            dofvalue_correction = dofvalue - ucur
            call fixdof_conjugategradient(idof, n, dofvalue_correction)
        else  ! Constrain a node set
            nnodes = nodeset_list(prescribeddof_list(k)%node_set)%n_nodes
            iof = nodeset_list(prescribeddof_list(k)%node_set)%index
            do i = 1,nnodes
                ignoredof = .false.
                n = node_lists(iof+i-1)
                idof = prescribeddof_list(k)%dof
                if (prescribeddof_list(k)%rate_flag==0) then
                    ucur = dof_increment(idof + node_list(n)%dof_index - 1) + dof_total(idof + node_list(n)%dof_index - 1)
                else
                    ucur = dof_increment(idof + node_list(n)%dof_index - 1)
                endif
                if (prescribeddof_list(k)%flag==3) then
                    iofs = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%index
                    nparam = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%n_parameters
                    call user_prescribeddof(n,idof,subroutine_parameters(iofs),nparam,dofvalue,ignoredof)
                endif
                if (ignoredof) cycle
                dofvalue_correction = dofvalue - ucur
                call fixdof_conjugategradient(idof, n, dofvalue_correction)
            end do
        endif
   
    end do
  
    unbalanced_force_norm = dsqrt(dot_product(rhs,rhs))


end subroutine apply_cojugategradient_boundaryconditions

subroutine fixdof_conjugategradient(idof,node_number,dofvalue)
    use Types
    use Mesh, only : node, node_list
    use Stiffness, only : unsymmetric_stiffness,ieqs,rows,cols,aupp,alow,diag,rhs,last_filled_aupp
    implicit none

    integer, intent( in )      :: idof            ! Degree of freedom to constrain
    integer, intent( in )      :: node_number     ! Node number to constrain
     
    real( prec ), intent( in )    :: dofvalue     ! Value to apply

    integer :: npresc,i

    npresc = ieqs(idof + node_list(node_number)%dof_index - 1)
    do i = 1,last_filled_aupp
        if (rows(i)==npresc) then
            if (unsymmetric_stiffness) then
                !         Unsymmetric - take lower half of stiffness over to rhs
                rhs(cols(i)) = rhs(cols(i)) - dofvalue*alow(i)
                alow(i) = 0.D0
            else
                !         Symmetric - use fact that aupp=alow to take column irow
                !         from lower half of stiffness over to rhs
                rhs(cols(i)) = rhs(cols(i)) - dofvalue*aupp(i)
            endif
            !       Zero the row for DOF npresc (also column in lower half of stiffness)
            aupp(i) = 0.D0
        endif
        if (cols(i)==npresc)  then
            !       Take column npresc from upper half of stiffness over to rhs
            rhs(rows(i)) = rhs(rows(i)) - dofvalue*aupp(i)
            aupp(i) = 0.D0
            !       zero row npresc in lower half of stiffness
            if (unsymmetric_stiffness) alow(i) = 0.D0
        endif
    end do
    diag(npresc) = 1.D0
    rhs(npresc) = dofvalue


end subroutine fixdof_conjugategradient


!====================Subroutine INSTIF_CG =======================
subroutine solve_conjugategradient(fail)
    use Types
    use Mesh, only : n_nodes, node, node_list, dof_increment, length_dofs,correction_norm
    use Stiffness, only : unsymmetric_stiffness,neq,ieqs,rows,cols,last_filled_aupp
    use Stiffness, only : diag,aupp,alow,rhs,sol
    use Boundaryconditions, only: n_constraints,lagrange_multipliers
    implicit none
    interface

        subroutine cg_iteration_loop(alow, aupp, diag, rhs,icol,irow,size,neq,sol,fail)
            use Types
            use ParamIO
            use Stiffness, only : p,r,z,precon
            implicit none

            integer, intent( in )         :: neq
            integer, intent( in )         :: size
            real( prec ), intent( inout ) :: alow(:)
            real( prec ), intent( inout ) :: aupp(:)
            real( prec ), intent( in )    :: diag(:)
            real( prec ), intent( inout ) :: rhs(:)
            real( prec ), intent( inout ) :: sol(:)
            integer, intent( in )         :: icol(:)
            integer, intent( in )         :: irow(:)
 
            logical, intent( inout )      :: fail

        end subroutine cg_iteration_loop
    end interface

    logical, intent( inout ) :: fail

    ! Local Variables
    integer :: n,iof,k,nc
  
    if ( unsymmetric_stiffness ) then
        call cg_iteration_loop(alow, aupp, diag, rhs,rows,cols,last_filled_aupp,neq,sol,fail)
    else
        call cg_iteration_loop(aupp, aupp, diag, rhs,rows,cols,last_filled_aupp,neq,sol,fail)
    end if

    if (fail) return

    do n = 1,n_nodes
        iof = node_list(n)%dof_index
        do k = 1,node_list(n)%n_dof
            dof_increment(iof+k-1) = dof_increment(iof+k-1) + sol(ieqs(iof+k-1))
        end do
    end do
 
    do nc = 1,n_constraints
        lagrange_multipliers(nc) = lagrange_multipliers(nc) + sol(ieqs(length_dofs+nc))
    end do
 
    correction_norm = dsqrt(dot_product(sol,sol))

end subroutine solve_conjugategradient

!=====================Subroutine SOLV =========================
subroutine cg_iteration_loop(alow, aupp, diag, rhs,irow,icol,size,neq,sol,fail)
    use Types
    use ParamIO
    use Stiffness, only : p,r,z,precon,cg_tolerance,max_cg_iterations
    implicit none

    integer, intent( in )         :: neq
    integer, intent( in )         :: size
    logical, intent( inout )      :: fail
    real( prec ), intent( inout ) :: alow(:)
    real( prec ), intent( inout ) :: aupp(:)
    real( prec ), intent( in )    :: diag(:)
    real( prec ), intent( inout ) :: rhs(:)
    real( prec ), intent( inout ) :: sol(:)
    integer, intent( in )         :: icol(:)
    integer, intent( in )         :: irow(:)

    real( prec ) :: ak, akden, bk, bkden, bknum, bnrm

    integer :: is,iter,i
    real (prec) :: err
    !
    ! Conjugate gradient solver with indexed storage scheme
    !
    fail = .false.
    err = 1.D0

    sol(1:neq) = 0.D0

    is = 0
    do
        is = is + 1
        if (is>neq) then
            write (IOW, 99001)
            sol(1:neq) = 0.D0
            return
        endif
        if ( rhs(is) /= 0.D0 ) exit
    end do

    !  Jacobi preconditioner
    precon(1:neq) = diag(1:neq)
    do i = 1,neq
        precon(i) = 1.d0/precon(i)
    end do
 
    !
    !  r = rhs - [A]*sol
    !
    r(1:neq) = 0.D0
    do i = 1,size
        r(irow(i)) = r(irow(i)) + aupp(i)*sol(icol(i))
        r(icol(i)) = r(icol(i)) + alow(i)*sol(irow(i))
    end do
    r(1:neq) = r(1:neq) + diag(1:neq)*sol(1:neq)

    r(1:neq) = rhs(1:neq) - r(1:neq)
    bnrm = dot_product(rhs(1:neq),rhs(1:neq))
    err = dot_product(r(1:neq),r(1:neq))/bnrm
    if (err<cg_tolerance) then
        write(IOW,9901) iter,neq,size
        return
    endif
    !
    !  ~
    ! [A] = [I]diag([A]) = preconditioning matrix
    !
    !       ~ -1     ~ -1
    !  z = [A]  r = [A]  (rhs - [A]*sol)
    !
    z(1:neq) = r(1:neq)*precon(1:neq)
    !
    !  iterate
    !
    iter = 0
    do while (err>cg_tolerance)

        iter=iter+1
        if (iter>max_cg_iterations) then
            write(IOW,9900) iter
9900        format(' Conjugate gradient solver failed to converge in ',i7,&
                ' iterations ')
            fail = .true.
            return
        endif
        !
        !  update search direction p
        !
        !  bk = z   dot r    / z dot r
        !        k+1     k+1    k     k
        !
        bknum=dot_product(z(1:neq),r(1:neq))
        if (iter == 1) then
            p(1:neq)=z(1:neq)
        else
            bk=bknum/bkden
            p(1:neq)=bk*p(1:neq)+z(1:neq)
        endif
        bkden=bknum
        !
        !  ak = r dot z / p dot ([A]p)
        !
        z(1:neq) = 0.D0
        do i = 1,size
            z(irow(i)) = z(irow(i)) + aupp(i)*p(icol(i))
            z(icol(i)) = z(icol(i)) + alow(i)*p(irow(i))
        end do
        z(1:neq) = z(1:neq) + diag(1:neq)*p(1:neq)
        akden=dot_product(z(1:neq),p(1:neq))
        ak=bknum/akden
        !
        !  update x and r
        !
        sol(1:neq)=sol(1:neq)+ak*p(1:neq)
        r(1:neq)=r(1:neq)-ak*z(1:neq)
        !
        !         ~ -1
        !  z   = [A]  r
        !   k+1

        z(1:neq) = r(1:neq)*precon(1:neq)
        !
        !  check convergence tolerance
        !
        err=dot_product(r(1:neq),r(1:neq))/bnrm
    enddo

    write(IOW,9901) iter,neq,size
9901 format(//'  Conjugate gradient solver converged in ',i7,' iterations'/ &
        '  No. equations:                         ',i7/ &
        '  Size of assembled stiffness matrix:    ',i10)


99001 format ( // ' ***** SOLVER WARNING ***** '/, &
        ' Zero RHS found in SOLV ')

end subroutine cg_iteration_loop
!
!========================== SUBROUTINE GETINDEX ===================
!
subroutine getindex(irow,icol,index)

    use Types
    use ParamIO
    use Linkedlist_Handling, only : databinsize, Integer_Linked_List
    use Stiffness, only : topea => last_filled_equation_adjacency,topaup => last_filled_aupp
    use Stiffness, only : length_equation_adjacency,length_aupp
    use Stiffness, only : cols,rows
    use Stiffness, only : equation_adjacency,equation_adjacency_index
    implicit none

    integer, intent( in )            :: irow
    integer, intent( in )            :: icol
    integer, intent( out )           :: index

    integer :: i,j
    !
    ! Routine to locate position AUPP(index) of entry K(irow,icol) of stiffness matrix
    ! New storage is assigned if the this entry has not yet been filled
    !
  
    !
    ! Search through link list to find entry icol
    !

    i = equation_adjacency_index(irow)

    !
    ! Check whether there are any existing entries for row irow
    !
    if (i==0) then

        !   If not, set the pointer to the top of the eqtoa link list
        topea = topea + 1
        if (topea>length_equation_adjacency) then
            write(IOW,*) ' Error in subroutine getindex '
            write(IOW,*) ' Insufficient storage to assemble row pointer'
            write(IOW,*) ' Increase length_equation_adjacency '
            stop
        endif
        i = topea
        equation_adjacency_index(irow) = topea
        topaup = topaup + 1
        if (topaup>length_aupp) THEN
            write(IOW,9900) length_aupp
9900        format ( // ' *** Error detected in subroutine GETINDEX ***', /,  &
                '   Insufficient storage to assemble indexed matrix ',  &
                /'   Parameter MAXSUP must be increased',  &
                /'   Its current value is ', I10)
            stop
        endif

        equation_adjacency(i)%n_entries = 1
        equation_adjacency(i)%data_list(1) = topaup
        rows(topaup) = irow
        cols(topaup) = icol
        index = topaup
        return

    endif

    ! Keep looping until we either find entry icol, or failing that create a
    ! new entry and quit
    do

        if (equation_adjacency(i)%next==0) then !     we're in the top data block

            do j = 1,equation_adjacency(i)%n_entries
                if (cols(equation_adjacency(i)%data_list(j))==icol) then
                    index = equation_adjacency(i)%data_list(j)
                    return
                endif
            end do
            !     Entry for column icol wasn't in the top data block.  Add it,
            !     creating a new data block if necessary
            if (equation_adjacency(i)%n_entries<databinsize) then   ! Add to an existing data block
                equation_adjacency(i)%n_entries = equation_adjacency(i)%n_entries + 1
                topaup = topaup + 1
                if (topaup>length_aupp) THEN
                    write(IOW,9900) length_aupp
                    stop
                endif
                equation_adjacency(i)%data_list(equation_adjacency(i)%n_entries) = topaup
                rows(topaup) = irow
                cols(topaup) = icol
            else            ! Create a new data bin for the current row
                topaup = topaup + 1
                if (topaup>length_aupp) THEN
                    write(IOW,*) ' Error in subroutine getindex '
                    write(IOW,*) ' Insufficient memory for stiffness '
                    write(IOW,*) ' Increase size of aupp '
                    stop
                endif
                topea = topea + 1
                if (topea>length_equation_adjacency) then
                    write(IOW,*) ' Error in subroutine getindex '
                    write(IOW,*) ' Insufficient storage to assemble equation adjacency'
                    write(IOW,*) ' Increase parameter length_equation_adjacency '
                    stop
                endif
                equation_adjacency(i)%next = topea
                equation_adjacency(topea)%n_entries = 1
                equation_adjacency(topea)%data_list(1) = topaup
                rows(topaup) = irow
                cols(topaup) = icol
            endif
            !     All done
            index = topaup
            return
        else                                        ! Look for current column in filled data block
            do j = 1,equation_adjacency(i)%n_entries
                if (cols(equation_adjacency(i)%data_list(j))==icol) then
                    index = equation_adjacency(i)%data_list(j)
                    return
                endif
            end do
            !     We didnt find column icol in the current data block - go to the next one
            i = equation_adjacency(i)%next
        endif
    end do

end subroutine getindex

subroutine initialize_cg_solver

    use Types
    use ParamIO
    use Linkedlist_Handling
    use Mesh, only: n_elements,element,element_list,connectivity
    use Mesh, only: n_nodes,node,node_list,length_dofs
    use Boundaryconditions, only : n_constraints,constraint,constraint_list
    use Boundaryconditions, only : nodeset,nodeset_list,node_lists
    use Stiffness
    implicit none
  
    logical, allocatable :: nzstiffness(:,:)
    integer, allocatable :: degree(:)

    integer :: n,i,j,ipn,mpc
    integer :: lmn,n1,node1,n2,node2,j1,j2,row,col,status
    integer :: dof1, dof2, dofset
    integer :: ns, nnode, iofc, iofd, iofn1, iofn2

    neq = node_list(n_nodes)%dof_index + node_list(n_nodes)%n_dof - 1 + n_constraints

    if (.not.allocated(ieqs)) then
        allocate(ieqs(neq), stat=status)
        if (status /=0) then
            write(IOW,*) ' Error in subroutine profile_equations '
            write(IOW,*) ' Unable to allocate storage for equation numebers '
            stop
        endif
    else if (size(ieqs)<neq) then
        deallocate(ieqs)
        allocate(ieqs(neq), stat=status)
        if (status /=0) then
            write(IOW,*) ' Error in subroutine profile_equations '
            write(IOW,*) ' Unable to allocate storage for equation numebers '
            stop
        endif
    endif

    ! Compute equation numbers.  Note that nodes are not re-numbered.
    neq = 0
    do n = 1, n_nodes+n_constraints
        if (n<n_nodes+1) then
            do j = 1, node_list(n)%n_dof
                neq = neq + 1
                ipn = node_list(n)%dof_index + j - 1
                ieqs(ipn) = neq
            end do
        else
            neq = neq + 1
            ieqs(length_dofs+n-n_nodes) = neq
        endif
    end do


    allocate(nzstiffness(n_nodes,n_nodes), stat=status)
    allocate(degree(neq), stat=status)
 
    if (status /=0) then
        write(IOW,*) ' Error in subroutine initialize_cg_solver '
        write(IOW,*) ' Unable to allocate memory for CG solver '
        stop
    endif
    ! Initialize storage for conjugate gradient solver

    nzstiffness = .false.
    degree = 0
   
    do lmn = 1,n_elements
        do n1 = 1,element_list(lmn)%n_nodes
            node1 = connectivity(element_list(lmn)%connect_index+n1-1)
            do n2 = 1,element_list(lmn)%n_nodes
                node2 = connectivity(element_list(lmn)%connect_index+n2-1)
                if (nzstiffness(node1,node2)) cycle
                nzstiffness(node1,node2) = .true.
                nzstiffness(node2,node1) = .true.
                do j1 = 1,node_list(node1)%n_dof
                    row = ieqs(node_list(node1)%dof_index+j1-1)
                    do j2 = 1,node_list(node2)%n_dof
                        col = ieqs(node_list(node2)%dof_index+j2-1)
                        if (col>row) then
                            degree(row) = degree(row) + 1
                        else if (col<row) then
                            degree(col) = degree(col) + 1
                        endif
                    end do
                end do
            end do
        end do
    end do
   
    do mpc = 1, n_constraints
        if (constraint_list(mpc)%flag<3) then
            node1 = constraint_list(mpc)%node1
            dof1 = constraint_list(mpc)%dof1
            node2 = constraint_list(mpc)%node2
            dof2 = constraint_list(mpc)%dof2
            row = ieqs(node_list(node1)%dof_index+dof1-1)
            col = ieqs(node_list(node2)%dof_index+dof2-1)
            if (.not.nzstiffness(node1,node2)) then
                if (col>row) then
                    degree(row) = degree(row) + 1
                else if (row>col) then
                    degree(col) = degree(col) + 1
                endif
                nzstiffness(node1,node2) = .true.
                nzstiffness(node2,node1) = .true.
            endif
            degree(row) = degree(row) + 1
            degree(col) = degree(col) + 1
            col = ieqs(length_dofs+mpc)
            degree(col) = degree(col) + 2
        else
            ns = constraint_list(mpc)%node1
            dofset = constraint_list(mpc)%node2
            nnode = nodeset_list(ns)%n_nodes
            iofc = nodeset_list(ns)%index
            iofd = nodeset_list(dofset)%index
            do i = 1, nnode           !     ---    Loop over nodes in set
                node1 = node_lists(i + iofc - 1)
                dof1 = node_lists(i + iofd - 1)
                iofn1 = node_list(node1)%dof_index    !     ---      Find smallest eqn. no. on current element
                row = ieqs(iofn1+dof1-1)
                degree(row) = degree(row) + 1
                do j = 1, nnode
                    node2 = node_lists(j+iofc-1)
                    if (nzstiffness(node1,node2)) cycle
                    dof2 = node_lists(j+iofd-1)
                    iofn2 = node_list(node2)%dof_index
                    col = ieqs(iofn2+dof2-1)
                    if (col>row) then
                        degree(row) = degree(row) + 1
                    else if (col<row) then
                        degree(col) = degree(col) + 1
                    endif
                    nzstiffness(node1,node2) = .true.
                    nzstiffness(node2,node1) = .true.
                end do
 
            end do
            degree(ieqs(length_dofs+mpc)) = degree(ieqs(length_dofs+mpc))+nnode
        endif

    end do
     
   
    length_aupp = sum(degree)
    length_equation_adjacency = 0
    do i = 1,neq
        if (degree(i)>0)  length_equation_adjacency = length_equation_adjacency + int(degree(i)/databinsize)+1
    end do
    deallocate(nzstiffness)
    deallocate(degree)

    if (allocated(diag)) deallocate(diag)
    if (allocated(aupp)) deallocate(aupp)
    if (allocated(alow)) deallocate(alow)
    if (allocated(rhs))  deallocate(rhs)
    if (allocated(sol))  deallocate(sol)
    if (allocated(equation_adjacency)) deallocate(equation_adjacency)
    if (allocated(equation_adjacency_index)) deallocate(equation_adjacency_index)
    if (allocated(p)) deallocate(p)
    if (allocated(r)) deallocate(r)
    if (allocated(z)) deallocate(z)
    if (allocated(precon)) deallocate(precon)
    if (allocated(rows)) deallocate(rows)
    if (allocated(cols)) deallocate(cols)

    allocate(diag(neq), stat = status)
    allocate(rhs(neq), stat = status)
    allocate(sol(neq), stat = status)
    allocate(rows(length_aupp), stat = status)
    allocate(cols(length_aupp), stat = status)
    allocate(aupp(length_aupp), stat = status)
    if (unsymmetric_stiffness) allocate(alow(length_aupp), stat = status)
 
    if (status /=0) then
        write(IOW,*) ' Error in subroutine allocate_direct_stiffness'
        write(IOW,*) ' Unable to allocate storage for equation numebers '
        stop
    endif
 
    aupp = 0.d0
    rhs = 0.d0
    diag = 0.d0

    allocate(equation_adjacency(length_equation_adjacency), stat = status)
    allocate(equation_adjacency_index(neq), stat = status)
    allocate(p(neq), stat = status)
    allocate(r(neq), stat = status)
    allocate(z(neq), stat = status)
    allocate(precon(neq), stat = status)

   
    if (status/=0) then
        write(IOW,*) ' Error in subroutine initialize_cg_solver '
        write(IOW,*) ' Unable to allocate memory for CG solver '
        stop
    endif
   
    equation_adjacency(1:length_equation_adjacency)%n_entries = 0
    equation_adjacency(1:length_equation_adjacency)%next = 0
    equation_adjacency_index = 0
 

end subroutine initialize_cg_solver
