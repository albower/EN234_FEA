subroutine explicit_dynamic_step
    use Types
    use ParamIO
    use Globals, only : TIME, DTIME
    use Mesh
    use Dynamicstepparameters
    use Printparameters
    implicit none
    
    integer :: status,neq,istep,i
    
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
    abq_minStableTimeStep = 1.d99
    abq_prevTimeStep = DTIME

    write(IOW,'(//A//)') '  Started explicit dynamic analysis '

    velocity(1:neq) = dof_increment(1:neq)
    acceleration(1:neq) = 0.d0
    do istep = 1,no_dynamic_steps

        current_dynamic_step = istep

        dof_increment(1:neq) = DTIME*(velocity(1:neq) + 0.5d0*DTIME*acceleration(1:neq))

        rforce = 0.d0
        call assemble_lumped_mass

        call apply_dynamic_boundaryconditions

        if (stateprint) then
            write(IOW,*) ' '
            write(IOW,*) '  Step ',istep,' completed'
            write(IOW,*) '  Elapsed time is ',TIME

            if (mod(istep,state_print_steps)==0) then
               call print_state
               write(IOW,*) '  Printing state... '
            endif
        endif

        if (userprint) then
            write(IOW,*) '  Step ',istep,' completed'
            write(IOW,*) '  User print activated '
            if (istep==1) then
               call user_print(istep)
            else if (mod(istep,user_print_steps)==0) then
               call user_print(istep)
            endif
        endif

        call assemble_dynamic_force

        acceleration(1:neq) = rforce(1:neq)/lumped_mass(1:neq)
        velocity(1:neq) = velocity(1:neq) + DTIME*acceleration(1:neq)
        dof_total(1:neq) = dof_total(1:neq) + dof_increment(1:neq)
        initial_state_variables = updated_state_variables               ! Update element state
        TIME = TIME + DTIME

        if (abq_minStableTimeStep<DTIME) then
            write(IOW,*) ' '
            write(IOW,*) '  *** Abaqus VUEL has cut time step to ',abq_minStableTimeStep
            write(IOW,*) '   '
            abq_prevTimeStep = DTIME
            DTIME = abq_minStableTimeStep
        endif

    end do

end subroutine explicit_dynamic_step


subroutine assemble_lumped_mass
    use Types
    use ParamIO
    use Globals
    use User_Subroutine_Storage
    use Mesh
    use Controlparameters, only : abaqusformat
    use Dynamicstepparameters
    implicit none

    ! Local Variables
    integer      :: i1, irow, iu, ix, j,k, status
    integer      :: lmn, n
    integer      :: nn, node1
    integer      :: iof, ns
  
    integer, parameter  :: abq_length_blocks=1

    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
                                                            
    real( prec ), allocatable   :: element_mass(:)

    real( prec ), allocatable    :: abq_rhs(:)
    real( prec ), allocatable    :: abq_amass(:)
    real( prec ), allocatable    :: abq_dtimeStable(:)
    real( prec ), allocatable    :: abq_dMassScaleFactor(:)
    real( prec ), allocatable    :: abq_svars(:)
    real( prec ), allocatable    :: abq_energy(:)
    real( prec ), allocatable    :: abq_PREDEF(:,:,:,:)
    real( prec ), allocatable    :: abq_ADLMAG(:)
    real( prec ), allocatable    :: abq_V(:)
    real( prec ), allocatable    :: abq_A(:)
    real( prec ), allocatable    :: abq_element_coords(:)
    real( prec ), allocatable    :: abq_element_dof_increment(:)
    real( prec ), allocatable    :: abq_element_dof_total(:)

    integer :: abq_JTYPE
    integer :: abq_JDLTYP
    integer :: abq_NPREDF
    integer :: abq_LFLAGS(3)
    integer :: el_number(1)

    real( prec ) :: abq_time(2)

    type (node), allocatable ::  local_nodes(:)

    !     Subroutine to assemble global stiffness matrix

    allocate(element_coords(length_coord_array), stat = status)
    allocate(element_dof_increment(length_dof_array), stat = status)
    allocate(element_dof_total(length_dof_array), stat = status)
    allocate(local_nodes(length_node_array), stat = status)
    allocate(element_mass(length_dof_array), stat = status)



    if (abaqusformat) then
        allocate(abq_rhs(abq_length_blocks*length_dof_array), stat=status)
        allocate(abq_amass(abq_length_blocks*length_dof_array*length_dof_array), stat=status)
        allocate(abq_dtimeStable(abq_length_blocks), stat=status)
        allocate(abq_dMassScaleFactor(abq_length_blocks), stat=status)
        allocate(abq_svars(abq_length_blocks*max(length_state_variable_array,1)), stat=status)
        allocate(abq_energy(12*abq_length_blocks), stat=status)
        allocate(abq_PREDEF(abq_length_blocks,length_node_array,1,2), stat=status)
        allocate(abq_ADLMAG(abq_length_blocks), stat=status)
        allocate(abq_V(abq_length_blocks*length_dof_array), stat=status)
        allocate(abq_A(abq_length_blocks*length_dof_array), stat=status)
        allocate(abq_element_coords(abq_length_blocks*length_node_array*3), stat=status)
        allocate(abq_element_dof_increment(abq_length_blocks*length_dof_array), stat=status)
        allocate(abq_element_dof_total(abq_length_blocks*length_dof_array), stat=status)
    endif
  
    if (status/=0) then
        write(IOW,*) ' Error in subroutine assemble_lumped_mass_matrix '
        write(IOW,*) ' Unable to allocate memory for user subroutines '
        stop
    endif
 
    lumped_mass = 0.d0

    do lmn = 1, n_elements

        if (element_deleted(lmn)) cycle

        if (element_list(lmn)%flag == 10002.or.element_list(lmn)%flag == 10003 ) then             ! Internal continuum element

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

            if (element_list(lmn)%density_index<1.or.element_list(lmn)%density_index>length_densities) then
                write(IOW,'(A)') ' *** Error in explicit dynamic analysis *** '
                write(IOW,'(A)') '     A density has not been defined for element number ',lmn
                stop
            endif

            call continuum_element_lumped_mass(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes,&                           ! Input variables
                local_nodes(1:element_list(lmn)%n_nodes), &                                                                  ! Input variables
                densities(element_list(lmn)%density_index),   &                                                              ! Input variables
                element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &     ! Input variables
                element_list(lmn)%n_int_element_properties,int_element_properties(element_list(lmn)%int_element_property_index), & ! Input variables
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

        else if (element_list(lmn)%flag>99999) then                ! ABAQUS VUEL format user subroutine

            abq_LFLAGS = 0
            abq_LFLAGS(2) = 1
            abq_LFLAGS(3) = 1                                 ! Requests mass calculation
            abq_time(1:2) = TIME
            abq_JTYPE = element_list(lmn)%flag-99999

            !           Change storage of element coords to ABAQUS UEL format
            if (abq_MCRD(lmn) == 0) then
                do j = 1,element_list(lmn)%n_nodes
                    n = connectivity(element_list(lmn)%connect_index + j - 1)
                    if (node_list(n)%n_coords >abq_MCRD(lmn)) then
                        abq_MCRD(lmn) = node_list(n)%n_coords
                    endif
                    if (node_list(n)%n_dof>abq_MCRD(lmn)) then
                        abq_MCRD(lmn) = node_list(n)%n_dof
                        if (abq_MCRD(lmn)>3) abq_MCRD(lmn) = 3
                    endif
                end do
            endif

            ix = 0
            do k = 1, abq_MCRD(lmn)
                do j = 1, element_list(lmn)%n_nodes
                    n = connectivity(element_list(lmn)%connect_index + j - 1)
                    ix = ix + 1
                    if (k<=node_list(n)%n_coords) then
                        abq_element_coords(ix) = coords(node_list(n)%coord_index + k - 1)
                    else
                        abq_element_coords(ix) = 0.d0
                    endif
                end do
            end do
            iu = 0
            do j = 1, element_list(lmn)%n_nodes
                n = connectivity(element_list(lmn)%connect_index + j - 1)
                do k = 1, node_list(n)%n_dof
                    iu = iu + 1
                    abq_element_dof_total(iu) = dof_increment(node_list(n)%dof_index + k - 1) + &
                        dof_total(node_list(n)%dof_index + k - 1)
                    abq_element_dof_increment(iu) = dof_increment(node_list(n)%dof_index + k - 1)
                    abq_V(iu) = velocity(node_list(n)%dof_index + k - 1)
                    abq_A(iu) = acceleration(node_list(n)%dof_index + k - 1)
                end do
            end do


            abq_PREDEF(1,1:element_list(lmn)%n_nodes,1,2) = BTEMP
            abq_PREDEF(1,1:element_list(lmn)%n_nodes,1,1) = BTEMP+BTINC
            abq_NPREDF = 1

            abq_JDLTYP = 0

            !     Form element contribution to nodal force vector
            iof = element_list(lmn)%state_index
            if (iof==0) iof = 1
            ns = element_list(lmn)%n_states
            if (ns==0) ns=1

            abq_svars(1:ns) = initial_state_variables(iof:iof+ns-1)

            abq_rhs(1:iu) = 0.d0
            abq_amass(1:iu*iu) = 0.d0
            el_number(1) = lmn
            call VUEL(1,abq_rhs(1:iu),abq_amass(1:iu*iu),abq_dtimeStable,&
                abq_svars(1:ns),ns, &
                abq_energy,element_list(lmn)%n_nodes,iu, &
                element_properties(element_list(lmn)%element_property_index),element_list(lmn)%n_element_properties, &
                int_element_properties(element_list(lmn)%int_element_property_index),element_list(lmn)%n_int_element_properties, &
                abq_element_coords,abq_MCRD(lmn),abq_element_dof_total,abq_element_dof_increment,abq_V,abq_A, &
                abq_JTYPE,el_number, &
                abq_time,total_dynamic_time,DTIME,abq_prevTimeStep,1,current_dynamic_step, &
                abq_LFLAGS, &
                abq_dMassScaleFactor, &
                abq_PREDEF,abq_NPREDF, &
                abq_JDLTYP, abq_ADLMAG)

            !     --   Add element lumped mass matrix to global array
            n = element_list(lmn)%n_nodes
            irow = 1
            do i1 = 1, n
                node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
                iof = node_list(node1)%dof_index
                nn= node_list(node1)%n_dof
                do j = 1,nn
                    lumped_mass(iof+j-1) = lumped_mass(iof+j-1) + abq_amass(irow)
                    irow = irow + n*nn+1
                end do
            end do

        else

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

            if (element_list(lmn)%density_index<1.or.element_list(lmn)%density_index>length_densities) then
                write(IOW,'(A)') ' *** Error in explicit dynamic analysis *** '
                write(IOW,'(A)') '     A density has not been defined for element number ',lmn
                stop
            endif

            call user_element_lumped_mass(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes,&                           ! Input variables
                local_nodes(1:element_list(lmn)%n_nodes), &                                                                  ! Input variables
                densities(element_list(lmn)%density_index),   &                                                              ! Input variables
                element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &     ! Input variables
                element_list(lmn)%n_int_element_properties,int_element_properties(element_list(lmn)%int_element_property_index), & ! Input variables
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

        endif

    end do

    if (abaqusformat) then
        deallocate(abq_rhs)
        deallocate(abq_amass)
        deallocate(abq_dtimeStable)
        deallocate(abq_dMassScaleFactor)
        if (allocated(abq_svars)) deallocate(abq_svars)
        deallocate(abq_energy)
        deallocate(abq_PREDEF)
        deallocate(abq_ADLMAG)
        deallocate(abq_V)
        deallocate(abq_A)
        deallocate(abq_element_coords)
        deallocate(abq_element_dof_increment)
        deallocate(abq_element_dof_total)
    endif


    return
  
end subroutine assemble_lumped_mass
  
subroutine assemble_dynamic_force
    use Types
    use ParamIO
    use Globals
    use User_Subroutine_Storage
    use Mesh
    use Stiffness
    use Controlparameters, only : abaqusformat
    use Dynamicstepparameters
    implicit none

 
    integer      :: status
    integer      :: lmn
    integer      :: ix,iu,j,n,k,iof,ns
    integer      :: irow,i1,node1,iofr
    integer      :: nn
    integer      :: mat_prop_index,n_mat_props

    integer, parameter  :: abq_length_blocks=1

    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
                                                            
    real( prec ), allocatable    :: element_stiffness(:,:)
    real( prec ), allocatable    :: element_residual(:)

    real( prec ), allocatable    :: abq_rhs(:)
    real( prec ), allocatable    :: abq_amass(:)
    real( prec ), allocatable    :: abq_dtimeStable(:)
    real( prec ), allocatable    :: abq_dMassScaleFactor(:)
    real( prec ), allocatable    :: abq_svars(:)
    real( prec ), allocatable    :: abq_energy(:)
    real( prec ), allocatable    :: abq_PREDEF(:,:,:,:)
    real( prec ), allocatable    :: abq_ADLMAG(:)
    real( prec ), allocatable    :: abq_V(:)
    real( prec ), allocatable    :: abq_A(:)
    real( prec ), allocatable    :: abq_element_coords(:)
    real( prec ), allocatable    :: abq_element_dof_increment(:)
    real( prec ), allocatable    :: abq_element_dof_total(:)

    integer :: abq_JTYPE
    integer :: abq_JDLTYP
    integer :: abq_NPREDF
    integer :: abq_LFLAGS(3)
    integer :: vuel_elementno(1)

    real( prec ) :: abq_time(2)

    character (len=80) material_name

    type (node), allocatable ::  local_nodes(:)

    !     Subroutine to assemble global stiffness matrix

    allocate(element_coords(length_coord_array), stat = status)
    allocate(element_dof_increment(length_dof_array), stat = status)
    allocate(element_dof_total(length_dof_array), stat = status)
    allocate(local_nodes(length_node_array), stat = status)
    allocate(element_stiffness(length_dof_array,length_dof_array), stat = status)
    allocate(element_residual(length_dof_array), stat = status)

    if (abaqusformat) then
        allocate(abq_rhs(abq_length_blocks*length_dof_array), stat=status)
        allocate(abq_amass(abq_length_blocks*length_dof_array*length_dof_array), stat=status)
        allocate(abq_dtimeStable(abq_length_blocks), stat=status)
        allocate(abq_dMassScaleFactor(abq_length_blocks), stat=status)
        allocate(abq_svars(abq_length_blocks*max(length_state_variable_array,1)), stat=status)
        allocate(abq_energy(12*abq_length_blocks), stat=status)
        allocate(abq_PREDEF(abq_length_blocks,length_node_array,1,2), stat=status)
        allocate(abq_ADLMAG(abq_length_blocks), stat=status)
        allocate(abq_V(abq_length_blocks*length_dof_array), stat=status)
        allocate(abq_A(abq_length_blocks*length_dof_array), stat=status)
        allocate(abq_element_coords(abq_length_blocks*length_node_array*3), stat=status)
        allocate(abq_element_dof_increment(abq_length_blocks*length_dof_array), stat=status)
        allocate(abq_element_dof_total(abq_length_blocks*length_dof_array), stat=status)
    endif


  
    if (status/=0) then
        write(IOW,*) ' Error in subroutine assemble_dynamic_force '
        write(IOW,*) ' Unable to allocate memory for user subroutines '
        stop
    endif

    abq_minStableTimeStep = 1.d99

    do lmn = 1, n_elements
        if (element_deleted(lmn)) cycle


        !     Form element contribution to nodal force vector

        if (element_list(lmn)%flag == 10002) then                  ! 2D continuum element


        else if (element_list(lmn)%flag == 10003) then             ! 3D continuum element

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


            iof = element_list(lmn)%state_index
            if (iof==0) iof = 1
            ns = element_list(lmn)%n_states
            if (ns==0) ns=1

            mat_prop_index = material_list(element_list(lmn)%material_index)%prop_index
            n_mat_props = material_list(element_list(lmn)%material_index)%n_properties
            material_name = material_namelist(element_list(lmn)%material_index)(1:80)

                   !     Form element contribution to nodal force vector
            call continuum_element_dynamic_3D(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes, &                            ! Input variables
                local_nodes(1:element_list(lmn)%n_nodes), &      ! Input variables
                densities(1),n_mat_props, material_properties(mat_prop_index),material_name,  &              ! Input variables
                element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu),iu,      &                   ! Input variables
                ns, initial_state_variables(iof:iof+ns-1), &                                                               ! Input variables
                updated_state_variables(iof:iof+ns-1),element_residual(1:iu),element_deleted(lmn))                           ! Output variables


            !     --   Add element force vector to global array
            irow = 1
            do i1 = 1, element_list(lmn)%n_nodes
                node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
                iofr = node_list(node1)%dof_index
                nn = node_list(node1)%n_dof
                rforce(iofr:iofr+nn-1) = rforce(iofr:iofr+nn-1) + element_residual(irow:irow+nn-1)
                irow = irow + nn
            end do




        else if (element_list(lmn)%flag>99999) then                ! ABAQUS VUEL format user subroutine

            abq_LFLAGS = 0
            abq_LFLAGS(2) = 1
            abq_LFLAGS(3) = 2                                 ! Requests internal force and new stable timestep
            abq_time(1:2) = TIME
            abq_JTYPE = element_list(lmn)%flag-99999

            !           Change storage of element coords to ABAQUS UEL format
            if (abq_MCRD(lmn) == 0) then
                do j = 1,element_list(lmn)%n_nodes
                    n = connectivity(element_list(lmn)%connect_index + j - 1)
                    if (node_list(n)%n_coords >abq_MCRD(lmn)) then
                        abq_MCRD(lmn) = node_list(n)%n_coords
                    endif
                    if (node_list(n)%n_dof>abq_MCRD(lmn)) then
                        abq_MCRD(lmn) = node_list(n)%n_dof
                        if (abq_MCRD(lmn)>3) abq_MCRD(lmn) = 3
                    endif
                end do
            endif
            ix = 0
            do k = 1, abq_MCRD(lmn)
                do j = 1, element_list(lmn)%n_nodes
                    n = connectivity(element_list(lmn)%connect_index + j - 1)
                    ix = ix + 1
                    if (k<=node_list(n)%n_coords) then
                        abq_element_coords(ix) = coords(node_list(n)%coord_index + k - 1)
                    else
                        abq_element_coords(ix) = 0.d0
                    endif
                end do
            end do
            iu = 0
            do j = 1, element_list(lmn)%n_nodes
                n = connectivity(element_list(lmn)%connect_index + j - 1)
                do k = 1, node_list(n)%n_dof
                    iu = iu + 1
                    abq_element_dof_total(iu) = dof_increment(node_list(n)%dof_index + k - 1) + &
                        dof_total(node_list(n)%dof_index + k - 1)
                    abq_element_dof_increment(iu) = dof_increment(node_list(n)%dof_index + k - 1)
                    abq_V(iu) = velocity(node_list(n)%dof_index + k - 1)
                    abq_A(iu) = acceleration(node_list(n)%dof_index + k - 1)
                end do
            end do

            abq_PREDEF(1,1:element_list(lmn)%n_nodes,1,1) = BTEMP
            abq_PREDEF(1,1:element_list(lmn)%n_nodes,1,2) = BTEMP+BTINC
            abq_NPREDF = 1

            abq_JDLTYP = 0

            !     Form element contribution to nodal force vector
            iof = element_list(lmn)%state_index
            if (iof==0) iof = 1
            ns = element_list(lmn)%n_states
            if (ns==0) ns=1

            abq_svars(1:ns) = initial_state_variables(iof:iof+ns-1)

            abq_rhs(1:iu) = 0.d0
            abq_amass(1:iu*iu) = 0.d0

            vuel_elementno(1) = lmn
            call VUEL(1,abq_rhs(1:iu),abq_amass(1:iu*iu),abq_dtimeStable,&
                abq_svars(1:ns),ns, &
                abq_energy,element_list(lmn)%n_nodes,iu, &
                element_properties(element_list(lmn)%element_property_index),element_list(lmn)%n_element_properties, &
                int_element_properties(element_list(lmn)%int_element_property_index),element_list(lmn)%n_int_element_properties, &
                abq_element_coords,abq_MCRD(lmn),abq_element_dof_total,abq_element_dof_increment,abq_V,abq_A, &
                abq_JTYPE,vuel_elementno, &
                abq_time,total_dynamic_time,DTIME,abq_prevTimeStep,1,current_dynamic_step, &
                abq_LFLAGS, &
                abq_dMassScaleFactor, &
                abq_PREDEF,abq_NPREDF, &
                abq_JDLTYP, abq_ADLMAG)

            if (abq_dtimeStable(1)<abq_minStableTimeStep) abq_minStableTimeStep = abq_dtimeStable(1)

            energy(12*lmn-11:12*lmn) = abq_energy(1:12)
            updated_state_variables(iof:iof+ns-1) = abq_svars(1:ns)


            !     --   Add element force vector to global array
            irow = 1
            do i1 = 1, element_list(lmn)%n_nodes
                node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
                iofr = node_list(node1)%dof_index
                nn = node_list(node1)%n_dof
                rforce(iofr:iofr+nn-1) = rforce(iofr:iofr+nn-1) + abq_rhs(irow:irow+nn-1)
                irow = irow + nn
            end do


        else

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


            iof = element_list(lmn)%state_index
            if (iof==0) iof = 1
            ns = element_list(lmn)%n_states
            if (ns==0) ns=1
                   !     Form element contribution to nodal force vector
            call user_element_dynamic(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes, &                            ! Input variables
                local_nodes(1:element_list(lmn)%n_nodes), &                                                                ! Input variables
                element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &   ! Input variables
                element_list(lmn)%n_int_element_properties,int_element_properties(element_list(lmn)%int_element_property_index), & ! Input variables
                element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu),iu,      &                   ! Input variables
                ns, initial_state_variables(iof:iof+ns-1), &                                                               ! Input variables
                updated_state_variables(iof:iof+ns),element_residual(1:iu),element_deleted(lmn))                           ! Output variables
      

            !     --   Add element force vector to global array
            irow = 1
            do i1 = 1, element_list(lmn)%n_nodes
                node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
                iofr = node_list(node1)%dof_index
                nn = node_list(node1)%n_dof
                rforce(iofr:iofr+nn-1) = rforce(iofr:iofr+nn-1) + element_residual(irow:irow+nn-1)
                irow = irow + nn
            end do

        endif

    end do

  
    deallocate(element_coords)
    deallocate(element_dof_increment)
    deallocate(element_dof_total)
    deallocate(local_nodes)
    deallocate(element_stiffness)
    deallocate(element_residual)
  

    if (abaqusformat) then
        deallocate(abq_rhs)
        deallocate(abq_amass)
        deallocate(abq_dtimeStable)
        deallocate(abq_dMassScaleFactor)
        if (allocated(abq_svars)) deallocate(abq_svars)
        deallocate(abq_energy)
        deallocate(abq_PREDEF)
        deallocate(abq_ADLMAG)
        deallocate(abq_V)
        deallocate(abq_A)
        deallocate(abq_element_coords)
        deallocate(abq_element_dof_increment)
        deallocate(abq_element_dof_total)
    endif

end subroutine assemble_dynamic_force
  
  
  
subroutine apply_dynamic_boundaryconditions
    use Types
    use ParamIO
    use Globals
    use User_Subroutine_Storage
    use Mesh
    use Boundaryconditions
    use Element_Utilities, only : facenodes
    use Controlparameters, only : abaqusformat
    use Dynamicstepparameters
    implicit none

    ! Local Variables
    logical :: ignoredof
    real( prec ) :: ucur, dofvalue, dloadvalue, force_value
    integer :: idof, ix,iu, i,j,k,kel, lmn, n,ns, iof, iofs,nhist, nparam, nnodes, status
    integer :: load,flag,elset,ifac,nel,param_index,ntract,ndims,ndof,nfacenodes
    integer :: i1,node1,irow,nn,iofr
    integer :: list(8)

    integer, parameter  :: abq_length_blocks=1
       
    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
                                                            
    real( prec ), allocatable   :: element_stiffness(:,:)
    real( prec ), allocatable   :: element_residual(:)
  
    real( prec ), allocatable   :: traction(:)
  
    real( prec ), allocatable    :: abq_rhs(:)
    real( prec ), allocatable    :: abq_amass(:)
    real( prec ), allocatable    :: abq_dtimeStable(:)
    real( prec ), allocatable    :: abq_dMassScaleFactor(:)
    real( prec ), allocatable    :: abq_svars(:)
    real( prec ), allocatable    :: abq_energy(:)
    real( prec ), allocatable    :: abq_PREDEF(:,:,:,:)
    real( prec ), allocatable    :: abq_ADLMAG(:)
    real( prec ), allocatable    :: abq_V(:)
    real( prec ), allocatable    :: abq_A(:)
    real( prec ), allocatable    :: abq_element_coords(:)
    real( prec ), allocatable    :: abq_element_dof_increment(:)
    real( prec ), allocatable    :: abq_element_dof_total(:)

    integer :: abq_JTYPE
    integer :: abq_JDLTYP
    integer :: abq_NPREDF
    integer :: abq_LFLAGS(3)
    integer :: vuel_elementno(1)

    real( prec ) :: abq_time(2)


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
    


        if (abaqusformat) then
            allocate(abq_rhs(abq_length_blocks*length_dof_array), stat=status)
            allocate(abq_amass(abq_length_blocks*length_dof_array*length_dof_array), stat=status)
            allocate(abq_dtimeStable(abq_length_blocks), stat=status)
            allocate(abq_dMassScaleFactor(abq_length_blocks), stat=status)
            allocate(abq_svars(abq_length_blocks*max(length_state_variable_array,1)), stat=status)
            allocate(abq_energy(12*abq_length_blocks), stat=status)
            allocate(abq_PREDEF(abq_length_blocks,length_node_array,1,2), stat=status)
            allocate(abq_ADLMAG(abq_length_blocks), stat=status)
            allocate(abq_V(abq_length_blocks*length_dof_array), stat=status)
            allocate(abq_A(abq_length_blocks*length_dof_array), stat=status)
            allocate(abq_element_coords(abq_length_blocks*length_node_array*3), stat=status)
            allocate(abq_element_dof_increment(abq_length_blocks*length_dof_array), stat=status)
            allocate(abq_element_dof_total(abq_length_blocks*length_dof_array), stat=status)
        endif

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

        do kel  = 1, nel
            lmn = element_lists(iof+kel-1)
            ndims = node_list(connectivity(element_list(lmn)%connect_index))%n_coords
            ndof = node_list(connectivity(element_list(lmn)%connect_index))%n_dof

            if (element_list(lmn)%flag>99999) then                ! ABAQUS VUEL format user subroutine

                abq_LFLAGS = 0
                abq_LFLAGS(2) = 1
                abq_LFLAGS(3) = 3                                 ! Requests distributed loads
                abq_time(1:2) = TIME

                abq_JTYPE = element_list(lmn)%flag-99999

                    !           Change storage of element coords to ABAQUS UEL format
                if (abq_MCRD(lmn) == 0) then
                    do j = 1,element_list(lmn)%n_nodes
                        n = connectivity(element_list(lmn)%connect_index + j - 1)
                        if (node_list(n)%n_coords >abq_MCRD(lmn)) then
                            abq_MCRD(lmn) = node_list(n)%n_coords
                        endif
                        if (node_list(n)%n_dof>abq_MCRD(lmn)) then
                            abq_MCRD(lmn) = node_list(n)%n_dof
                            if (abq_MCRD(lmn)>3) abq_MCRD(lmn) = 3
                        endif
                    end do
                endif
                ix = 0
                do k = 1, abq_MCRD(lmn)
                    do j = 1, element_list(lmn)%n_nodes
                        n = connectivity(element_list(lmn)%connect_index + j - 1)
                        ix = ix + 1
                        if (k<=node_list(n)%n_coords) then
                            abq_element_coords(ix) = coords(node_list(n)%coord_index + k - 1)
                        else
                            abq_element_coords(ix) = 0.d0
                        endif
                    end do
                end do

                iu = 0
                do j = 1, element_list(lmn)%n_nodes
                    n = connectivity(element_list(lmn)%connect_index + j - 1)
                    do k = 1, node_list(n)%n_dof
                        iu = iu + 1
                        abq_element_dof_total(iu) = dof_increment(node_list(n)%dof_index + k - 1) + &
                            dof_total(node_list(n)%dof_index + k - 1)
                        abq_element_dof_increment(iu) = dof_increment(node_list(n)%dof_index + k - 1)

                        abq_V(iu) = velocity(node_list(n)%dof_index + k - 1)
                        abq_A(iu) = acceleration(node_list(n)%dof_index + k - 1)
                    end do
                end do

                abq_PREDEF(1,1:element_list(lmn)%n_nodes,1,1) = BTEMP
                abq_PREDEF(1,1:element_list(lmn)%n_nodes,1,2) = BTEMP+BTINC
                abq_NPREDF = 1

                abq_JDLTYP = flag

                !     Form element contribution to nodal force vector
                iof = element_list(lmn)%state_index
                if (iof==0) iof = 1
                ns = element_list(lmn)%n_states
                if (ns==0) ns=1

                abq_svars(1:ns) = initial_state_variables(iof:iof+ns-1)

                abq_rhs(1:iu) = 0.d0
                abq_amass(1:iu*iu) = 0.d0

                if (distributedload_list(n_distributedloads)%history_number>0) then
                    iof = history_list(distributedload_list(load)%history_number)%index
                    nhist = history_list(distributedload_list(load)%history_number)%n_timevals
                    call interpolate_history_table(history_data(1,iof),nhist,TIME+DTIME,dloadvalue)
                    abq_adlmag(1) = dloadvalue
                else
                    abq_adlmag(1) = dload_values(distributedload_list(load)%index_dload_values)
                endif

                vuel_elementno(1) = lmn
                call VUEL(1,abq_rhs(1:iu),abq_amass(1:iu*iu),abq_dtimeStable,&
                    abq_svars(1:ns),ns, &
                    abq_energy,element_list(lmn)%n_nodes,iu, &
                    element_properties(element_list(lmn)%element_property_index),element_list(lmn)%n_element_properties, &
                    int_element_properties(element_list(lmn)%int_element_property_index), &
                    element_list(lmn)%n_int_element_properties, &
                    abq_element_coords,abq_MCRD(lmn),abq_element_dof_total+abq_element_dof_increment, &
                    abq_element_dof_increment,abq_V,abq_A, &
                    abq_JTYPE,vuel_elementno, &
                    abq_time,total_dynamic_time,DTIME,abq_prevTimeStep,1,current_dynamic_step, &
                    abq_LFLAGS, &
                    abq_dMassScaleFactor, &
                    abq_PREDEF,abq_NPREDF, &
                    abq_JDLTYP, abq_ADLMAG)


                !     --   Add element force vector to global array
                irow = 1
                do i1 = 1, element_list(lmn)%n_nodes
                    node1 = connectivity(i1 + element_list(lmn)%connect_index - 1)
                    iofr = node_list(node1)%dof_index
                    nn = node_list(node1)%n_dof
                    rforce(iofr:iofr+nn-1) = rforce(iofr:iofr+nn-1) + abq_rhs(irow:irow+nn-1)
                    irow = irow + nn
                end do

            else

                if (flag<4) then
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
                    if (flag==2) traction(1:ntract) = traction(1:ntract)/dsqrt(dot_product(traction(1:ntract),traction(1:ntract)))
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
                        element_dof_increment(1:iu),element_dof_total(1:iu),iu,traction(1:ntract),ntract,&
                        element_residual(1:iu))               ! Output variables
                    irow = 1
                    do i1 = 1,nfacenodes
                        node1 = connectivity(element_list(lmn)%connect_index + list(i1) - 1)
                        iofr = node_list(node1)%dof_index
                        nn = node_list(node1)%n_dof
                        rforce(iofr:iofr+nn-1) = rforce(iofr:iofr+nn-1) + element_residual(irow:irow+nn-1)
                        irow = irow + nn
                    end do

                else if (flag==4) then

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

                endif

            endif
        end do

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
!                if (prescribeddof_list(k)%rate_flag==0) then
!                    ucur =  dof_total(idof + node_list(n)%dof_index - 1)
!                else
!                    ucur = 0.d0
!                endif
                if (prescribeddof_list(k)%flag==3) then
                    iofs = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%index
                    nparam = subroutineparameter_list(prescribeddof_list(k)%subroutine_parameter_number)%n_parameters
                    call user_prescribeddof(n,idof,subroutine_parameters(iofs),nparam,dofvalue,ignoredof)
                endif
                if (ignoredof) cycle
                 if (prescribeddof_list(k)%rate_flag==0) then
                    ucur =  dof_total(idof + node_list(n)%dof_index - 1)
                    dof_increment(idof + node_list(n)%dof_index - 1) = dofvalue-ucur
                 else
                    dof_increment(idof + node_list(n)%dof_index - 1) = DTIME*dofvalue
                 endif
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
        if (abaqusformat) then
            deallocate(abq_rhs)
            deallocate(abq_amass)
            deallocate(abq_dtimeStable)
            deallocate(abq_dMassScaleFactor)
            if (allocated(abq_svars)) deallocate(abq_svars)
            deallocate(abq_energy)
            deallocate(abq_PREDEF)
            deallocate(abq_ADLMAG)
            deallocate(abq_V)
            deallocate(abq_A)
            deallocate(abq_element_coords)
            deallocate(abq_element_dof_increment)
            deallocate(abq_element_dof_total)
        endif
    endif

end subroutine apply_dynamic_boundaryconditions
