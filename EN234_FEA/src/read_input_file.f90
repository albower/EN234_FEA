subroutine read_input_file
    use Types
    use ParamIO
    use Globals
    use Mesh
    use Boundaryconditions
    use Controlparameters
    use Staticstepparameters
    use Dynamicstepparameters
    use Stiffness, only : unsymmetric_stiffness
    use Printparameters
    use Fileparameters
    use User_Subroutine_Storage
    implicit none

    ! Local variables

    character (len=100) uel_string    ! String used to read ABAQUS element type or boundary condition
 
    integer :: ncoor          ! # coords for current node
    integer :: ndof           ! # degrees of freedom for current node
    integer :: nodtyp         ! identifier for current node
    integer :: iopt           ! Flag specifying whether increment of DOF or value of total DOF is read from file
    integer :: nnod           ! no. nodes on an element
    integer :: nsvar          ! no. state vars on an element
    integer :: n_initial_state_vars ! No. initial user defined state vars
    integer :: n_int_pts      ! No. integration points on an internal continuum element
    integer :: lmntyp         ! Element identifier for current element
    integer :: nprops         ! No. properties for elements
    integer :: n_int_props    ! No. integer valued properties for element
    integer :: property_index ! Index to element property array
    integer :: int_property_index ! Index to ingteger element property array
    integer :: material_index ! Index to material list
    integer :: k              ! Counter
    integer :: iof            ! index
    integer :: nodenum        ! node number
    integer :: nodeset1       ! number of first node set in a constraint
    integer :: nodeset2       ! number of second node set in a constraint
    integer :: iof1           ! index of start of node set 1
    integer :: iof2           ! index of start of node set 2
    integer :: izone          ! Zone number to be printed
    integer :: iz             ! counter
    integer :: lmn            ! counter
    integer :: n              ! counter
    integer :: kk             ! counter
    integer :: nn             ! Used to store node number
    integer :: displacement_map_currentindex  ! Index of displacement map for a node set
    integer :: iun            ! Unit number for files
    integer :: ix             ! Coord counter
    integer :: iu             ! DOF counter
    integer :: nc             ! Constraint index
    integer :: nset           ! Node set counter
    integer :: j,i
    integer :: start_node     ! Start node for a node set, generate key
    integer :: end_node       ! end node for a node set, generate
    integer :: increment      ! Increment for a node set, generate key
    integer :: itst           ! Counter
    integer :: status         ! Status flag for memory allocation
   
    real (prec) :: coord_scale_factor   !  coordinates read from file are scaled by this factor

    logical :: strcmp        ! True if two strings agree

    call allocate_mesh_storage
    rewind(IOR)

    ! Procedure flags
    printinitialmesh = .false.
    checkstiffness = .false.
    checktangent = .false.
    staticstep = .false.
    explicitdynamicstep = .false.
    userstep = .false.
    abaqusformat = .false.

    ! Mesh data
    n_zones = 0
    n_nodes = 0
    n_elements = 0
    n_materials = 0
    length_coords = 0
    length_dofs = 0
    length_displacement_map = 0
    displacement_map_currentindex = 1
    length_element_properties = 0
    length_int_element_properties = 0
    length_material_properties = 0
    property_index = 1
    int_property_index = 1
    length_densities = 0
    length_state_variables = 0
    n_initial_state_vars = 0
    length_connectivity = 0
    if (allocated(zone_namelist)) zone_namelist = 'Unnamed Zone'
    if (allocated(node_list)) node_list(1:n_nodes)%displacement_map_index = 0
    if (allocated(node_list)) node_list(1:n_nodes)%n_displacements = 0
    n_mesh_parameters = 0
    n_int_props = 0
  
    ! Boundary condition data
    n_histories = 0
    n_subroutine_parameters = 0
    n_constraint_parameters = 0
    n_nodesets = 0
    n_elementsets = 0
      
    n_prescribeddof = 0
    n_prescribedforces = 0
    n_distributedloads = 0
    n_constraints = 0
  
    length_element_lists = 0
    length_node_lists = 0
    length_history_data = 0
    length_subroutine_parameters = 0
    length_dof_values = 0
    length_dload_values = 0
    length_constraint_parameters = 0

    ! Default values for time stepping parameters
    timestep_max = 0.d0
    timestep_min = 0.d0
    max_no_steps = 1
    max_total_time = 0.d0
    n_user_print_parameters=0
    n_user_print_files=0
    userprint = .false.
    stateprint = .false.
    nonlinear = .true.

    ! Solution data
    if (allocated(initial_state_variables) ) initial_state_variables = 0.d0
    if (allocated(updated_state_variables) ) updated_state_variables = 0.d0
    if (allocated(dof_total) ) dof_total = 0.d0
    if (allocated(dof_increment) )dof_increment = 0.d0
    if (allocated(rforce) ) rforce = 0.d0
    if (allocated(element_properties) ) element_properties = 0.d0
    if (allocated(int_element_properties) ) int_element_properties = 0
    if (allocated(lagrange_multipliers)) lagrange_multipliers = 0.d0

    ! Default print control parameters
    print_displacedmesh = .false.
    zone_print_flag = .true.
    use_lumped_projection_matrix = .false.
    zone_dimension = 0
    zone_ndof = 0
    print_dof = .false.
    n_field_variables = 0
    displacementscalefactor = 1.d0
    n_total_files = 0

    do while (.true.)
  
        iblnk=1
        read (IOR, 99001, ERR = 200, end = 500) strin
        if (echo) write(IOW,*) strin
        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
        if (iblnk==1) cycle
  
        if ( strcmp(strpar(1), 'MESH', 4) ) then
 
            do while ( .true. )
                !     Read a line of the input file
                iblnk =1
                read (IOR, 99001, ERR = 200, end = 500) strin
                if (echo) write(IOW,*) strin
                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                if (iblnk==1) cycle
                !     ---------------------- SOME PORTION OF A MESH IS GENERATED IN A USER SUBROUTINE ---------------------
                if ( strcmp(strpar(1), 'USERSUB',7) ) then
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDUS',5) ) then
                 
                            exit
                        endif
                        do k = 1,nstr
                            if (ityp(k)<2) then
                                n_mesh_parameters = n_mesh_parameters + 1
                                read(strpar(k),*) mesh_subroutine_parameters(n_mesh_parameters)
                            else
                                write(IOW,'(A)') ' *** Error in input file ***'
                                write(IOW,'(A)') ' Expecting a list of parameters defining user-subroutine generated mesh'
                                write(IOW,'(A)') ' Found '
                                write(IOW,'(A)') strin
                            endif
                        end do
                    end do

                    call user_mesh(n_mesh_parameters,mesh_subroutine_parameters)

                !     ---------------------- READ NODAL DATA ---------------------
                else if ( strcmp(strpar(1), 'NODE', 4) ) then

                    do while ( .true. )
                        !     Read a line of the input file
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if ( strcmp(strpar(1), 'PARAM',5) ) then
                            read (strpar(2), *) ncoor
                            read (strpar(3), *) ndof
                            if (nstr==4.and.ityp(4)==0) then
                                read (strpar(4), *) nodtyp
                            else
                                nodtyp = 0
                            endif
                        else if ( strcmp(strpar(1), 'DISPLACEMENTDOF', 15) ) then
                            if (ityp(2)==2) then
                                if (strcmp(strpar(2),'NONE',4) ) then
                                    displacement_map_currentindex = length_displacement_map + 1
                                else
                                    write(IOW,*) ' *** Error in input file ***'
                                    write(IOW,*) ' Expecting a list of DOF for displacements or NONE'
                                    write(IOW,*) ' following a DISPLACEMENT DOF key'
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif
                            else
                                displacement_map_currentindex = length_displacement_map + 1
                                do k = 1,nstr-1
                                    length_displacement_map = length_displacement_map + 1
                                    read (strpar(k+1), *) displacement_map(length_displacement_map)
                                end do
                            endif
                        else if ( strcmp(strpar(1), 'COOR', 4) ) then                                        !      Read nodal coordinates
                            coord_scale_factor = 1.d0
                            if (nstr>1) then
                                Read(strpar(2),*) coord_scale_factor
                            endif
                            if (ndof==0) then
                                write(IOW,*) ' *** Error in input file ***'
                                write(IOW,*) ' Attempt to create nodes with no degrees of freedom '
                                write(IOW,*) ' PARAMETERS must be specified for nodes before COORDINATES '
                                stop
                            endif
                            do while ( .true. )
                                read (IOR, 99001, ERR = 200, end = 500) strin
                                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                if ( iblnk==1 ) cycle
                                if ( strcmp(strpar(1), 'ENDCOOR', 7) ) then
                                    exit
                                endif
                                if ( ityp(1)==2 ) then
                                    write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                    write(IOW,*) ' COORDINATES key must be terminated by an END COORDINATES'
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif
                                if ( nstr>ncoor+1.or.nstr<ncoor ) then
                                    write(IOW,*) '  *** ERROR DETECTED IN INPUT FILE ***'
                                    write(IOW,*) ' Expecting ', ncoor,' coords of a node.  Found ',nstr
                                    stop
                                end if
                                n_nodes = n_nodes + 1
                                node_list(n_nodes)%flag = nodtyp
                                node_list(n_nodes)%coord_index = length_coords + 1
                                node_list(n_nodes)%n_coords = ncoor
                                node_list(n_nodes)%dof_index = length_dofs + 1
                                node_list(n_nodes)%n_dof = ndof
                                node_list(n_nodes)%displacement_map_index = displacement_map_currentindex
                                node_list(n_nodes)%n_displacements = length_displacement_map-displacement_map_currentindex+1
                                do k = nstr-ncoor+1, nstr
                                    if ( ityp(k)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' COORDINATES key must be terminated by an END COORDINATES'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    length_coords = length_coords + 1
                                    read (strpar(k), *) coords(length_coords)
                                    coords(length_coords) = coords(length_coords)*coord_scale_factor
                                end do
                                length_dofs = length_dofs + ndof
                            end do
                        else if ( strcmp(strpar(1), 'CREATE', 6 ) ) then                                               ! Create nodes with no coordinates
                            if (nstr/=2 .or. ityp(2) ==2) then
                                write(IOW,*) ' *** Error in input file ***'
                                write(IOW,*) ' Expecting an integer defining number of nodes following a CREATE NODES key '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            read (strpar(2), *) n
                            do k = 1,n
                                n_nodes = n_nodes + 1
                                node_list(n_nodes)%flag = nodtyp
                                node_list(n_nodes)%coord_index = length_coords + 1
                                node_list(n_nodes)%n_coords = ncoor
                                node_list(n_nodes)%dof_index = length_dofs + 1
                                node_list(n_nodes)%n_dof = ndof
                                node_list(n_nodes)%displacement_map_index = displacement_map_currentindex
                                node_list(n_nodes)%n_displacements = length_displacement_map-displacement_map_currentindex+1
                                length_dofs = length_dofs + ndof
                            end do
                        !     --------- read initial nodal dof  ----------
                        else if ( strcmp(strpar(1), 'INIT', 4) ) then

                            if (n_nodes==0) then
                                write(IOW,*) ' *** Error detected in input file ***'
                                write(IOW,*) ' Nodes must be defined before initial DOF can be created'
                                stop
                            endif

                            iopt = 1
                            if (strcmp(strpar(2),'VEL',3)) iopt = 2
        
                            do while ( .true. )
                                read (IOR, 99001, ERR = 200, end = 500) strin
                                if (echo) write(IOW,*) strin
                                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                if ( iblnk==1 ) cycle

                                if ( strcmp(strpar(1), 'ENDINIT', 7) ) then
                                    exit
                                elseif (strcmp(strpar(1), 'USERSUBR', 8) ) then
                                    iof = n_mesh_parameters+1
                                    do while (.true.)
                                        read (IOR, 99001, ERR = 200, end = 500) strin
                                        if (echo) write(IOW,*) strin
                                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                        if ( iblnk==1 ) cycle
                                        if (strcmp(strpar(1),'ENDUS',5) ) then
                                            exit
                                        endif
                                        do k = 1,nstr
                                            if (ityp(k)<2) then
                                                n_mesh_parameters = n_mesh_parameters + 1
                                                read(strpar(k),*) mesh_subroutine_parameters(n_mesh_parameters)
                                            else
                                                write(IOW,'(A)') ' *** Error in input file ***'
                                                write(IOW,'(A)') ' Expecting a list of parameters for user-subroutine generated DOF'
                                                write(IOW,'(A)') ' Found '
                                                write(IOW,'(A)') strin
                                            endif
                                        end do
                                    end do

                                    if (n_mesh_parameters-iof==0) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' No parameters were supplied for user subroutine generating DOF '
                                        stop
                                    endif

                                    do nodenum = 1,n_nodes
                                        iof1 = node_list(nodenum)%dof_index
                                        iof2 = node_list(nodenum)%coord_index
                                        ndof = node_list(nodenum)%n_dof
                                        ncoor = node_list(nodenum)%n_coords
                                        call user_initialdof(mesh_subroutine_parameters(iof:iof+n_mesh_parameters-1), &
                                            n_mesh_parameters-iof+1, &
                                            nodenum,coords(iof2:iof2+ncoor-1),ncoor, &
                                            dof_total(iof1:iof1+ndof-1),dof_increment(iof1:iof1+ndof-1),ndof)
                                    end do


                                elseif (strcmp(strpar(1), 'ALLNODES', 8) ) then
                                    ! Initialize DOF for all nodes with one line

                                    do nodenum = 1,n_nodes
                                        iof = node_list(nodenum)%dof_index
                                        if (node_list(nodenum)%n_dof /= nstr-1) then
                                            write(IOW,*) ' *** Error detected in input file ***'
                                            write(IOW,*) ' Incorrect number of DOF specified in an INITIAL DOF statement'
                                            write(IOW,*) ' node ',nodenum,' has ',node_list(nodenum)%n_dof,' dof'
                                            stop
                                        endif
                                        if (iopt==1) then
                                            do k = 1, nstr-1
                                                read (strpar(k+1), *) dof_total(iof+k-1)
                                            end do
                                        else
                                            do k = 1, nstr-1
                                                read (strpar(k+1), *) dof_increment(iof+k-1)
                                            end do
                                        endif
                                    end do
                                else
                                    if (ityp(1)/=0) then
                                        write(IOW,*) ' Error detected in input file'
                                        write(IOW,*) ' Expecting node number following INITIAL CONDITIONS key or END INITIAL DOF '
                                        write(IOW,*) 'Found'
                                        write(IOW,*) strin
                                        stop
                                    endif
			
                                    read(strpar(1),*) nodenum
			  
                                    if (nodenum>n_nodes) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' Attempted to define degrees of freedom for node ',nodenum
                                        write(IOW,*) ' Only ',n_nodes,' nodes are defined'
                                        stop
                                    endif
                                    if (nstr/=node_list(nodenum)%n_dof+1) then
                                        write(IOW,*) ' Error detected in input file '
                                        write(IOW,*) ' Expecting ',node_list(nodenum)%n_dof,' initial DOFs for node ',nodenum
                                        write(IOW,*) ' Found only ',nstr-1
                                        stop
                                    endif

                                    iof = node_list(nodenum)%dof_index
                                    if (iopt==1) then
                                        do k = 1, nstr-1
                                            read (strpar(k+1), *) dof_total(iof+k-1)
                                        end do
                                    else
                                        do k = 1, nstr-1
                                            read (strpar(k+1), *) dof_increment(iof+k-1)
                                        end do
                                    endif
                                endif
                            end do
                        else if ( strcmp(strpar(1), 'ENDNODE', 7) ) then
                            exit
                        else
                            write(IOW,*) ' Error detected in input file '
                            write(IOW,*) ' Expecting END NODES keyword'
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                        endif
                    end do
                !   Define a new material
                else if ( strcmp(strpar(1), 'MATE', 4) ) then

                    if (nstr<2) then
                       write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE ***'
                       write(IOW,*) ' MATERIAL key must specify a name for the material '
                       stop
                    endif

                    n_materials = n_materials + 1
                    material_namelist(n_materials) = strpar(2)
                    material_list(n_materials)%n_states = 0
                    nprops = 0
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if (strcmp(strpar(1), 'STAT', 4) ) then
                            if (nstr==2.and.ityp(2)==0) then
                               read(strpar(2),*) material_list(n_materials)%n_states
                            else
                               write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                               write(IOW,*) ' STATE VARIABLES key in a MATERIAL definition must specify no. state variables '
                               write(IOW,*) ' Found '
                               write(IOW,*) strin
                               stop
                            endif

                        else if (strcmp(strpar(1), 'PROP', 4) ) then

                            property_index = length_material_properties+1
                            nprops = 0
                            do while (.true.)
                                iblnk =1
                                read (IOR, 99001, ERR = 200, end = 500) strin
                                if (echo) write(IOW,*) strin
                                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                if (iblnk==1) cycle
                                if ( strcmp(strpar(1), 'ENDPROP', 7) ) then
                                    exit
                                endif
                                if ( ityp(1)==2 ) then
                                    write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                    write(IOW,*) ' PROPERTIES key in a MATERIAL block must be terminated by an END PROPERTIES'
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif
                                do k = 1,nstr
                                    length_material_properties = length_material_properties+1
                                    nprops = nprops + 1
                                    if (ityp(k)==2) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' Material properties must be real or integer values '
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    read (strpar(k),*) material_properties(length_material_properties)
                                end do
                            end do

                        else if ( strcmp(strpar(1), 'ENDMATE', 7) ) then
                            material_list(n_materials)%prop_index = property_index
                            material_list(n_materials)%n_properties = nprops
                            exit
                        else
                            write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                            write(IOW,*) ' MATERIAL key must be terminated by an END MATERIAL'
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                    end do



                !   Read element properties
                else if ( strcmp(strpar(1), 'ELEM', 4) ) then
                    n_initial_state_vars = 0
                    if (nstr<2) then
                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                        write(IOW,*) ' The ELEMENT key must be followed by a USER or INTERNAL option'
                        stop
                    endif

                    if (strcmp(strpar(2),'USER',4) ) then

                        do while (.true.)
                            iblnk =1
                            read (IOR, 99001, ERR = 200, end = 500) strin
                            if (echo) write(IOW,*) strin
                            call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                            if (iblnk==1) cycle

                            if ( strcmp(strpar(1), 'PARAM', 5) ) then

                                if ( nstr<4 .or. (ityp(2)/=0) .or. (ityp(3)/=0) .or.nstr>4 ) then
                                    write(IOW,*) '  *** ERROR DETECTED IN INPUT FILE ***'
                                    write(IOW,*)  ' No. nodes per element, no. state vars and element identifier'
                                    write(IOW,*)  ' must be supplied with element keyword'
                                    write(IOW,*)  ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif

                                read (strpar(2), *) nnod
                                read (strpar(3), *) nsvar
                                if (ityp(4)==0) then
                                    read (strpar(4), *) lmntyp
                                else
                                    if (strcmp(strpar(4),'U',1)) then  ! ABAQUS UEL
                                        abaqusformat = .true.
                                        uel_string = strpar(4)
                                        read(uel_string(2:100), *) lmntyp
                                        lmntyp = lmntyp + 99999
                                    else if (strcmp(strpar(4),'VU',2)) then  ! ABAQUS VUEL
                                        abaqusformat = .true.
                                        uel_string = strpar(4)
                                        read(uel_string(3:100), *) lmntyp
                                        lmntyp = lmntyp + 99999
                                    else
                                        write(IOW,'(A)')  ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,'(A)')  ' Element identifier must be an integer, or Un (for ABAQUS UEL),',&
                                            ' or VUn (for ABAQUS VUEL '
                                        write(IOW,'(A)')  ' Found'
                                        write(IOW,*) strin
                                        stop
                                    endif
                                endif
                            elseif (strcmp(strpar(1), 'DENS',4) ) then
                                length_densities = length_densities + 1
                                if (nstr<2.or.ityp(2)>1) then
                                    write(IOW,'(A)')  ' *** ERROR DETECTED IN INPUT FILE *** '
                                    write(IOW,'(A)')  ' Expecting density value following DENSITY key '
                                    write(IOW,'(A)')  ' Found'
                                    write(IOW,*) strin
                                    stop
                                endif
                                read (strpar(2), *) densities(length_densities)
                            else if (strcmp(strpar(1), 'PROP',4) ) then
                                property_index = length_element_properties+1
                                material_index = 0
                                nprops = 0
                                do while (.true.)
                                    iblnk =1
                                    read (IOR, 99001, ERR = 200, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDPROP', 7) ) then
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' PROPERTIES key must be terminated by an END PROPERTIES'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    do k = 1,nstr
                                        length_element_properties = length_element_properties+1
                                        nprops = nprops + 1
                                        if (ityp(k)==2) then
                                            write(IOW,*) ' *** Error in input file *** '
                                            write(IOW,*) ' Element properties must be real or integer values '
                                            write(IOW,*) ' Found '
                                            write(IOW,*) strin
                                            stop
                                        endif
                                        read (strpar(k),*) element_properties(length_element_properties)
                                    end do
                                end do
                            else if (strcmp(strpar(1), 'INTEGERPROP',11) ) then
                                int_property_index = length_int_element_properties+1
                                material_index = 0
                                do while (.true.)
                                    iblnk =1
                                    read (IOR, 99001, ERR = 200, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDINTE', 7) ) then
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' INTEGER PROPERTIES key must be terminated by an END INTEGER PROPERTIES'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    do k = 1,nstr
                                        length_int_element_properties = length_int_element_properties+1
                                        n_int_props = n_int_props + 1
                                        if (ityp(k)==2) then
                                            write(IOW,*) ' *** Error in input file *** '
                                            write(IOW,*) ' Element properties must be real or integer values '
                                            write(IOW,*) ' Found '
                                            write(IOW,*) strin
                                            stop
                                        endif
                                        read (strpar(k),*) int_element_properties(length_int_element_properties)
                                    end do
                                end do
                            else if (strcmp(strpar(1), 'INITIALS',8) ) then
                                n_initial_state_vars = 0
                                do while (.true.)
                                    iblnk =1
                                    read (IOR, 99001, ERR = 200, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDINIT', 7) ) then
                                        if (n_initial_state_vars==0) exit
                                        if (n_initial_state_vars /= nsvar) then
                                           write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE ***'
                                           write(IOW,*) ' Incorrect number of initial state variables was found for an element set '
                                           write(IOW,*) ' Expecting ',nsvar,' state variables but found only ',n_initial_state_vars
                                           stop
                                        endif
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' INITIAL STATE VARIABLES key must be terminated by', &
                                                     ' an END INITIAL STATE VARIABLES'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    do k = 1,nstr
                                        n_initial_state_vars = n_initial_state_vars + 1
                                        if (ityp(k)==2) then
                                            write(IOW,*) ' *** Error in input file *** '
                                            write(IOW,*) ' Element properties must be real or integer values '
                                            write(IOW,*) ' Found '
                                            write(IOW,*) strin
                                            stop
                                        endif
!                                       Updated state vars is used as temporary storage
                                        read (strpar(k),*) updated_state_variables(n_initial_state_vars)
                                    end do
                                end do
                            else if (strcmp(strpar(1), 'CONN',4) ) then
                                n_zones = n_zones + 1
                                zone_list(n_zones)%start_element = n_elements+1
                                if (nstr>1.and.ityp(2)==2) then
                                    read(strpar(2),*) zone_namelist(n_zones)
                                endif
                                do while (.true.)

                                    if (nnod==0) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' Attempt to create elements with no nodes '
                                        write(IOW,*) ' PARAMETERS must be specified for elements before CONNECTIVITY '
                                        stop
                                    endif
                                    iblnk =1
                                    read (IOR, 99001, ERR = 200, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDCONN', 7) ) then
                                        zone_list(n_zones)%end_element = n_elements
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' CONNECTIVITY key must be terminated by an END CONNECTIVITY'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    if (nstr < nnod.or.nstr > nnod+1) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' Expecting element number (optional) and ',nnod,&
                                            ' nodes for an element but found ',nstr
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    n_elements = n_elements + 1
                                    element_list(n_elements)%flag = lmntyp
                                    element_list(n_elements)%connect_index = length_connectivity + 1
                                    element_list(n_elements)%n_nodes = nnod
                                    element_list(n_elements)%state_index = length_state_variables + 1
                                    element_list(n_elements)%n_states = nsvar
                                    if (n_initial_state_vars>0) then
                                       do k=1,nsvar
                                          initial_state_variables(length_state_variables+k) = updated_state_variables(k)
                                       end do
                                    endif
                                    element_list(n_elements)%element_property_index = property_index
                                    element_list(n_elements)%n_element_properties = nprops
                                    element_list(n_elements)%density_index = length_densities
                                    element_list(n_elements)%int_element_property_index = int_property_index
                                    element_list(n_elements)%n_int_element_properties = n_int_props
                                    element_list(n_elements)%material_index = 0
                                    do k = 1,nnod
                                        length_connectivity = length_connectivity +1
                                        read(strpar(k+nstr-nnod),*) connectivity(length_connectivity)
                                        if (connectivity(length_connectivity)<1.or.connectivity(length_connectivity)>n_nodes) then
                                           write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE ***'
                                           write(IOW,*) ' A nonexistent node number was found in a connectivity list '
                                           write(IOW,*) strin
                                           stop
                                        endif
                                    end do
                                    length_state_variables = length_state_variables + nsvar
                                end do
                            else if ( strcmp(strpar(1), 'ENDELEM', 7) ) then
                                if (n_initial_state_vars>0) updated_state_variables = 0.d0
                                exit
                            else
                                write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                write(IOW,*) ' ELEMENTS key must be terminated by an END ELEMENTS'
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        end do

                    else if (strcmp(strpar(2),'INTER',5) ) then
                        do while (.true.)
                            iblnk =1
                            read (IOR, 99001, ERR = 200, end = 500) strin
                            if (echo) write(IOW,*) strin
                            call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                            if (iblnk==1) cycle
                            if (strcmp(strpar(1), 'TYPE', 4) ) then
                                if ( nstr.ne.2 ) then
                                    write(IOW,*) '  *** ERROR DETECTED IN INPUT FILE ***'
                                    write(IOW,*)  ' TYPE key in an ELEMENT definition must specify the element type'
                                    write(IOW,*)  ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif

                                if (strcmp(strpar(2),'CPE3',4) ) then
                                   lmntyp = 10002
                                   nnod = 3
                                   n_int_pts = 1
                                   nsvar = 11
                                else if (strcmp(strpar(2),'CPE4',4) ) then
                                   lmntyp = 10002
                                   nnod = 4
                                   n_int_pts = 4
                                   nsvar = 11
                                else if (strcmp(strpar(2),'CPE6',4) ) then
                                   lmntyp = 10002
                                   nnod = 6
                                   n_int_pts = 4
                                   nsvar = 11
                                else if (strcmp(strpar(2),'CPE8',4) ) then
                                   lmntyp = 10002
                                   nnod = 8
                                   n_int_pts = 9
                                   nsvar = 11
                                else if (strcmp(strpar(2),'C3D4',4) ) then
                                   lmntyp = 10003
                                   nnod = 4
                                   n_int_pts = 1
                                   nsvar = 15
                                else if (strcmp(strpar(2),'C3D10',4) ) then
                                   lmntyp = 10003
                                   nnod = 10
                                   n_int_pts = 4
                                   nsvar = 15
                                else if (strcmp(strpar(2),'C3D8',4) ) then
                                   lmntyp = 10003
                                   nnod = 8
                                   n_int_pts = 8
                                   nsvar = 15
                                else if (strcmp(strpar(2),'C3D20',4) ) then
                                   lmntyp = 10003
                                   nnod = 20
                                   n_int_pts = 27
                                   nsvar = 15
                                else
                                   write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                   write(IOW,*) ' Unrecognized element type specified with TYPE key in ELEMENT, INTERNAL block'
                                   write(IOW,*) ' Found '
                                   write(IOW,*) strin
                                   stop
                                endif

                            elseif (strcmp(strpar(1), 'DENS',4) ) then
                                length_densities = length_densities + 1
                                if (nstr<2.or.ityp(2)>1) then
                                    write(IOW,'(A)')  ' *** ERROR DETECTED IN INPUT FILE *** '
                                    write(IOW,'(A)')  ' Expecting density value following DENSITY key '
                                    write(IOW,'(A)')  ' Found'
                                    write(IOW,*) strin
                                    stop
                                endif
                                read (strpar(2), *) densities(length_densities)
                            else if (strcmp(strpar(1), 'PROP',4) ) then
                                property_index = 1
                                int_property_index = 1
                                nprops = 0
                                n_int_props = 0
                                material_index = find_name_index(strpar(2),lenstr(2),material_namelist,n_materials)
                                nsvar = nsvar + material_list(material_index)%n_states
                            else if (strcmp(strpar(1), 'INITIALS',8) ) then
                                write(6,*) ' Entered initial state vars '
                                n_initial_state_vars=0
                                do while (.true.)
                                    iblnk =1
                                    read (IOR, 99001, ERR = 200, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDINIT', 7) ) then
                                        if (n_initial_state_vars==0) exit
                                        if (n_initial_state_vars /= material_list(material_index)%n_states) then
                                           write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                           write(IOW,*) ' An incorrect number of initial state', &
                                                        ' variables was specified for a material '
                                           write(IOW,*) ' Expecting ', material_list(material_index)%n_states, &
                                           ' states but found ',n_initial_state_vars
                                           stop
                                        endif
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' INITIAL STATE VARIABLES key must be', &
                                                     ' terminated by an END INITIAL STATE VARIABLES'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    do k = 1,nstr
                                        n_initial_state_vars = n_initial_state_vars+1
                                        if (ityp(k)==2) then
                                            write(IOW,*) ' *** Error in input file *** '
                                            write(IOW,*) ' Initial state variables must be real or integer values '
                                            write(IOW,*) ' Found '
                                            write(IOW,*) strin
                                            stop
                                        endif
                                        !   Updated state variables used as temporary storage
                                        read (strpar(k),*) updated_state_variables(n_initial_state_vars)
                                    end do
                                end do
                            else if (strcmp(strpar(1), 'CONN',4) ) then
                                n_zones = n_zones + 1
                                zone_list(n_zones)%start_element = n_elements+1
                                if (nstr>1.and.ityp(2)==2) then
                                    read(strpar(2),*) zone_namelist(n_zones)
                                endif
                                do while (.true.)

                                    if (nnod==0) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' Attempt to create elements with no nodes '
                                        write(IOW,*) ' PARAMETERS must be specified for elements before CONNECTIVITY '
                                        stop
                                    endif
                                    iblnk =1
                                    read (IOR, 99001, ERR = 200, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDCONN', 7) ) then
                                        zone_list(n_zones)%end_element = n_elements
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' CONNECTIVITY key must be terminated by an END CONNECTIVITY'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    if (nstr < nnod.or.nstr > nnod+1) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' Expecting element number (optional) and ',nnod,&
                                            ' nodes for an element but found ',nstr
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    n_elements = n_elements + 1
                                    element_list(n_elements)%flag = lmntyp
                                    element_list(n_elements)%connect_index = length_connectivity + 1
                                    element_list(n_elements)%n_nodes = nnod
                                    element_list(n_elements)%state_index = length_state_variables + 1
                                    element_list(n_elements)%n_states = nsvar*n_int_pts
                                    if (n_initial_state_vars>0) then
                                       do k = 1,n_int_pts
                                          do n = 1,n_initial_state_vars
                                             iof = length_state_variables+1 + nsvar*k-n_initial_state_vars+n-1
                                             initial_state_variables(iof) = updated_state_variables(n)
                                          end do
                                       end do
                                    endif
                                    element_list(n_elements)%element_property_index = property_index
                                    element_list(n_elements)%n_element_properties = 0
                                    element_list(n_elements)%density_index = length_densities
                                    element_list(n_elements)%int_element_property_index = int_property_index
                                    element_list(n_elements)%n_int_element_properties = 0
                                    element_list(n_elements)%material_index = material_index
                                    do k = 1,nnod
                                        length_connectivity = length_connectivity +1
                                        read(strpar(k+nstr-nnod),*) connectivity(length_connectivity)
                                        if (connectivity(length_connectivity)<1.or.connectivity(length_connectivity)>n_nodes) then
                                           write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE ***'
                                           write(IOW,*) ' A nonexistent node number was found in a connectivity list '
                                           write(IOW,*) strin
                                           stop
                                        endif
                                    end do
                                    length_state_variables = length_state_variables + nsvar*n_int_pts
                                end do
                            else if ( strcmp(strpar(1), 'ENDELEM', 7) ) then
                                if (n_initial_state_vars>0) updated_state_variables = 0.d0
                                exit
                            else
                                write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                write(IOW,*) ' ELEMENTS key must be terminated by an END ELEMENTS'
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        end do




                    else

                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                        write(IOW,*) ' The ELEMENT key must be followed by USER or INTERNAL option '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif

                else if ( strcmp(strpar(1), 'ENDMESH', 7) ) then

                    if (abaqusformat) then
                       allocate(abq_MCRD(n_elements), stat=status)
                       if (status/=0) then
                            write(IOW,*) ' *** Error in subroutine read_input_file *** '
                            write(IOW,*) ' Unable to allocate storage for abq_MCRD array'
                            stop
                       endif
                       abq_MCRD = 0
                    endif
                    exit
        
                else
                    write(IOW,*) ' Error detected in input file '
                    write(IOW,*) ' Expecting an END MESH statement'
                    write(IOW,*) ' Found '
                    write(IOW,*) strin
                    stop
                endif

            end do


        else if ( strcmp(strpar(1), 'BOUNDARYCON', 11) ) then
            do while (.true.)
                iblnk =1
                read (IOR, 99001, ERR = 200, end = 500) strin
                if (echo) write(IOW,*) strin 
                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                if (iblnk==1) cycle

                if ( strcmp(strpar(1), 'HISTORY', 7) ) then
                    n_histories = n_histories + 1
                    read(strpar(2),*) history_namelist(n_histories)
                    history_list(n_histories)%index = length_history_data+1
                    history_list(n_histories)%n_timevals = 0


                    if (nstr<2) then
                        write(IOW,*) ' *** Error in input file '
                        write(IOW,*) ' HISTORY key was used without specifying a name for the history '
                        stop
                    endif

                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDHIST', 7) ) then
                            exit
                        endif
                        length_history_data = length_history_data+1
                        history_list(n_histories)%n_timevals = history_list(n_histories)%n_timevals+1
                        read(strpar(1), *) history_data(1,length_history_data)
                        read(strpar(2), *) history_data(2,length_history_data)
                    end do

                else if ( strcmp(strpar(1), 'SUBROUTINEP', 11) ) then
                    n_subroutine_parameters = n_subroutine_parameters + 1
                    read(strpar(2),*) subroutineparameter_namelist(n_subroutine_parameters)
                    subroutineparameter_list(n_subroutine_parameters)%index = length_subroutine_parameters+1
                    subroutineparameter_list(n_subroutine_parameters)%n_parameters = 0
            
                    do while (.true.)

                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDSUBR', 7) ) then
                            exit
                        endif
                
                        do k = 1,nstr
                            if (ityp(k) ==2) then
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting a list of parameters in a SUBROUTINE PARAMETERS block'
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            length_subroutine_parameters = length_subroutine_parameters + 1
                            subroutineparameter_list(n_subroutine_parameters)%n_parameters = &
                                subroutineparameter_list(n_subroutine_parameters)%n_parameters+1
                            read(strpar(k),*) subroutine_parameters(length_subroutine_parameters)
                        end do
                    end do

                else if ( strcmp(strpar(1), 'CONSTRAINTP', 11) ) then
                    n_constraint_parameters = n_constraint_parameters + 1
                    read(strpar(2),*) constraintparameter_namelist(n_constraint_parameters)
                    constraintparameter_list(n_constraint_parameters)%index = length_constraint_parameters+1
                    constraintparameter_list(n_constraint_parameters)%n_parameters = 0
            
                    do while (.true.)

                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDCONS', 7) ) then
                            exit
                        endif
                
                        do k = 1,nstr
                            if (ityp(k) ==2) then
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting a list of parameters in a CONSTRAINT PARAMETERS block'
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            length_constraint_parameters = length_constraint_parameters + 1
                            constraintparameter_list(n_constraint_parameters)%n_parameters = &
                                constraintparameter_list(n_constraint_parameters)%n_parameters+1
                            read(strpar(k),*) constraint_parameters(length_constraint_parameters)
                        end do
                    end do
       
                else if (strcmp(strpar(1), 'NODESET', 7) ) then
          
                    n_nodesets = n_nodesets + 1
                    read(strpar(2),*) nodeset_namelist(n_nodesets)
                    nodeset_list(n_nodesets)%index = length_node_lists + 1
                    nodeset_list(n_nodesets)%n_nodes = 0
                    if (nstr==3) then
                        do while (.true.)
                            iblnk =1
                            read (IOR, 99001, ERR = 500, end = 500) strin
                            if (echo) write(IOW,*) strin
                            call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                            if (iblnk==1) cycle

                            if ( strcmp(strpar(1), 'ENDNODE', 7) ) then
                                exit
                            endif
                            read(strpar(1), *) start_node
                            read(strpar(2), *) end_node
                            read(strpar(3), *) increment
                            if (increment==0) increment=1
                            do k = start_node,end_node,increment
                                length_node_lists = length_node_lists + 1
                                node_lists(length_node_lists) = k
                                nodeset_list(n_nodesets)%n_nodes = nodeset_list(n_nodesets)%n_nodes + 1
                            end do
                        end do
                    else
                        do while (.true.)
                            iblnk =1
                            read (IOR, 99001, ERR = 200, end = 500) strin
                            if (echo) write(IOW,*) strin
                            call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                            if (iblnk==1) cycle

                            if ( strcmp(strpar(1), 'ENDNODE', 7) ) then
                                exit
                            endif
                            do k = 1,nstr
                                length_node_lists = length_node_lists + 1
                                read(strpar(k),*) node_lists(length_node_lists)
                                nodeset_list(n_nodesets)%n_nodes = nodeset_list(n_nodesets)%n_nodes + 1
                            end do
                        end do
                    endif
                else if (strcmp(strpar(1), 'ELEMENTSET', 10) ) then
          
                    n_elementsets = n_elementsets + 1
                    read(strpar(2),*) elementset_namelist(n_elementsets)
                    elementset_list(n_elementsets)%index = length_element_lists + 1
                    elementset_list(n_elementsets)%n_elements = 0
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDELEM', 7) ) then
                            exit
                        endif
                        do k = 1,nstr
                            length_element_lists = length_element_lists + 1
                            read(strpar(k),*) element_lists(length_element_lists)
                            elementset_list(n_elementsets)%n_elements = elementset_list(n_elementsets)%n_elements + 1
                        end do
                    end do
                else if (strcmp(strpar(1), 'DEGREES', 7) ) then
           
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDDEGR', 7) ) then
                            exit
                        endif
                        n_prescribeddof = n_prescribeddof + 1
                        if (ityp(1)==0) then
                            read(strpar(1),*) prescribeddof_list(n_prescribeddof)%node_number
                            prescribeddof_list(n_prescribeddof)%node_set = 0
                            read(strpar(2),*) prescribeddof_list(n_prescribeddof)%dof
                            if (prescribeddof_list(n_prescribeddof)%dof> &
                                node_list(prescribeddof_list(n_prescribeddof)%node_number)%n_dof) then
                                write(IOW,*) ' *** Error in input file '
                                write(IOW,*) ' Attempt to prescribe dof number ', &
                                    prescribeddof_list(n_prescribeddof)%dof,' at node ', &
                                    prescribeddof_list(n_prescribeddof)%node_number
                                write(IOW,*) ' This node has only ', &
                                    node_list(prescribeddof_list(n_prescribeddof)%node_number)%n_dof,' DOF '
                                stop
                            endif
                        else if (ityp(1)==2) then
                            prescribeddof_list(n_prescribeddof)%node_set = &
                                find_name_index(strpar(1),lenstr(1),nodeset_namelist,n_nodesets)

                            read(strpar(2),*) prescribeddof_list(n_prescribeddof)%dof

                            do kk = 1,nodeset_list(prescribeddof_list(n_prescribeddof)%node_set)%n_nodes

                                nn = node_lists(nodeset_list(prescribeddof_list(n_prescribeddof)%node_set)%index+kk-1)
                                if (prescribeddof_list(n_prescribeddof)%dof>node_list(nn)%n_dof) then
                                    write(IOW,*) ' *** Error in input file '
                                    write(IOW,*) ' Attempt to prescribe dof number ', &
                                        prescribeddof_list(n_prescribeddof)%dof,' at node ',nn
                                    write(IOW,*) ' This node has only ',node_list(nn)%n_dof,' DOF '
                                    stop
                                endif
                            end do
                        else
                            write(IOW,*) ' *** Error in input file ***'
                            write(IOW,*) ' Expecting a node number or node set name in a degree of freedom'
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                
                        if (strcmp(strpar(3),'VALUE',5)) then
                            prescribeddof_list(n_prescribeddof)%flag = 1
                            length_dof_values = length_dof_values + 1
                            prescribeddof_list(n_prescribeddof)%index_dof_values = length_dof_values
                            prescribeddof_list(n_prescribeddof)%history_number = 0
                            prescribeddof_list(n_prescribeddof)%subroutine_parameter_number = 0
                            read(strpar(4),*) dof_values(length_dof_values)
                        else if (strcmp(strpar(3),'HISTORY',5)) then
                            prescribeddof_list(n_prescribeddof)%flag = 2
                            prescribeddof_list(n_prescribeddof)%index_dof_values = 0
                            prescribeddof_list(n_prescribeddof)%history_number = &
                                find_name_index(strpar(4),lenstr(4),history_namelist,n_histories)
                        else if (strcmp(strpar(3),'SUBROUT',5)) then
                            prescribeddof_list(n_prescribeddof)%flag = 3
                            prescribeddof_list(n_prescribeddof)%index_dof_values = 0
                            prescribeddof_list(n_prescribeddof)%history_number = 0
                            prescribeddof_list(n_prescribeddof)%subroutine_parameter_number = &
                                find_name_index(strpar(4),lenstr(4),subroutineparameter_namelist,n_subroutine_parameters)
                        else
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' Expecting key VALUE, HISTORY, or SUBROUTINE in a dof definition '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                
                        if ( strcmp(strpar(nstr),'RATE',4) ) then
                            prescribeddof_list(n_prescribeddof)%rate_flag = 1
                        else
                            prescribeddof_list(n_prescribeddof)%rate_flag = 0
                        endif
                                                 
                    end do


                else if (strcmp(strpar(1), 'FORCES', 6) ) then
           
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDFORC', 7) ) then
                            exit
                        endif
                        n_prescribedforces = n_prescribedforces + 1
                        if (ityp(1)==0) then
                            read(strpar(1),*) prescribedforce_list(n_prescribedforces)%node_number
                            prescribedforce_list(n_prescribedforces)%node_set = 0
                        else if (ityp(1)==2) then
                            prescribedforce_list(n_prescribedforces)%node_set = &
                                find_name_index(strpar(1),lenstr(1),nodeset_namelist,n_nodesets)
                        else
                            write(IOW,*) ' *** Error in input file ***'
                            write(IOW,*) ' Expecting a node number or node set name in a prescribed force'
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        read(strpar(2),*) prescribedforce_list(n_prescribedforces)%dof
                
                        if (strcmp(strpar(3),'VALUE',5)) then
                            prescribedforce_list(n_prescribedforces)%flag = 1
                            length_dof_values = length_dof_values + 1
                            prescribedforce_list(n_prescribedforces)%index_dof_values = length_dof_values
                            prescribedforce_list(n_prescribedforces)%history_number = 0
                            prescribedforce_list(n_prescribedforces)%subroutine_parameter_number = 0
                            read(strpar(4),*) dof_values(length_dof_values)
                        else if (strcmp(strpar(3),'HISTORY',5)) then
                            prescribedforce_list(n_prescribedforces)%flag = 2
                            prescribedforce_list(n_prescribedforces)%index_dof_values = 0
                            prescribedforce_list(n_prescribedforces)%history_number = &
                                find_name_index(strpar(4),lenstr(4),history_namelist,n_histories)
                        else if (strcmp(strpar(3),'SUBROUT',5)) then
                            prescribedforce_list(n_prescribedforces)%flag = 3
                            prescribedforce_list(n_prescribedforces)%index_dof_values = 0
                            prescribedforce_list(n_prescribedforces)%history_number = 0
                            prescribedforce_list(n_prescribedforces)%subroutine_parameter_number = &
                                find_name_index(strpar(4),lenstr(4),subroutineparameter_namelist,n_subroutine_parameters)
                        else
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' Expecting key VALUE, HISTORY, or SUBROUTINE in a dof definition '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                                                           
                    end do

                else if (strcmp(strpar(1), 'DISTRIB', 7) ) then
           
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDDISTRIB', 10) ) then
                            exit
                        endif
                        n_distributedloads = n_distributedloads + 1
                        distributedload_list(n_distributedloads)%element_set = &
                            find_name_index(strpar(1),lenstr(1),elementset_namelist,n_elementsets)


                            if (ityp(2) >0) then  ! Apply a boundary condition to an ABAQUS uel Format is: element set, Un(NU), value/history name
                                if (.not.strcmp(strpar(2),'U',1)) then
                                    write(IOW,*) ' Error in input file '
                                    write(IOW,*) ' Distributed load boundary conditions must specify a face number or Un/UnNU '
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif
                                !
                                !                       Check element type
                                !
                                lmn = element_lists(elementset_list(distributedload_list(n_distributedloads)%element_set)%index)
                                if (element_list(lmn)%flag < 99999) then
                                    write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                    write(IOW,*) ' An ABAQUS format boundary condition can only be applied to an ABAQUS element '
                                    stop
                                endif
                                uel_string = strpar(2)
                            do i = 2, len(uel_string)
                                !     Check for anything other than a number in the string
                                itst = ichar(uel_string(i:i))
                                if ( itst<48 .or. itst>57 ) exit
                            end do
                            read(uel_string(2:i-1), *) distributedload_list(n_distributedloads)%flag
                            if (uel_string(i:i) == 'N') then
                              distributedload_list(n_distributedloads)%flag = -distributedload_list(n_distributedloads)%flag
                            endif
                            distributedload_list(n_distributedloads)%index_dload_values = length_dload_values + 1
                            distributedload_list(n_distributedloads)%n_dload_values = 1
                            length_dload_values = length_dload_values + 1
                            if (ityp(3)<2) then
                                read(strpar(3), *) dload_values(length_dload_values)
                                distributedload_list(n_distributedloads)%history_number = 0
                            else
                                distributedload_list(n_distributedloads)%history_number = &
                                    find_name_index(strpar(3),lenstr(3),history_namelist,n_histories)
                                dload_values(length_dload_values) = 1.d0
                            endif
                        else
                            if (nstr<4) then
                                write(IOW,*) ' Error in input file'
                                write(IOW,*) ' Expecting element set name, face #,'
                                write(IOW,*) '    VALUE/HISTORY/NORMAL or SUBROUTINE,  history/parameter name, or value, values '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            read(strpar(2),*) distributedload_list(n_distributedloads)%face


                            if (strcmp(strpar(3), 'VALUE', 5)) then
                                distributedload_list(n_distributedloads)%flag = 1
                                distributedload_list(n_distributedloads)%index_dload_values = length_dload_values + 1
                                distributedload_list(n_distributedloads)%n_dload_values = 0
                                do k = 4,nstr
                                    if (ityp(k)==2) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' Expecting a list of traction components in a distributed load definition '
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    length_dload_values = length_dload_values + 1
                                    distributedload_list(n_distributedloads)%n_dload_values = &
                                        distributedload_list(n_distributedloads)%n_dload_values + 1
                                    read(strpar(k), *) dload_values(length_dload_values)
                                end do
                            else if (strcmp(strpar(3), 'HISTORY', 5)) then
                                distributedload_list(n_distributedloads)%flag = 2
                                distributedload_list(n_distributedloads)%index_dload_values = length_dload_values + 1
                                distributedload_list(n_distributedloads)%n_dload_values = 0
                                distributedload_list(n_distributedloads)%history_number = &
                                    find_name_index(strpar(4),lenstr(4),history_namelist,n_histories)
                                do k = 5,nstr
                                    if (ityp(k)==2) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' Expecting a list of direction components in a distributed load definition '
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    length_dload_values = length_dload_values + 1
                                    distributedload_list(n_distributedloads)%n_dload_values = &
                                        distributedload_list(n_distributedloads)%n_dload_values + 1
                                    read(strpar(k), *) dload_values(length_dload_values)

                                end do
                            else if (strcmp(strpar(3), 'NORMAL', 5)) then
                                distributedload_list(n_distributedloads)%flag = 3
                                distributedload_list(n_distributedloads)%index_dload_values = 0
                                distributedload_list(n_distributedloads)%n_dload_values = 0
                                distributedload_list(n_distributedloads)%history_number = &
                                    find_name_index(strpar(4),lenstr(4),history_namelist,n_histories)
                            else if (strcmp(strpar(3), 'SUBROUTINE', 10)) then
                                distributedload_list(n_distributedloads)%flag = 4
                                distributedload_list(n_distributedloads)%index_dload_values = 0
                                distributedload_list(n_distributedloads)%n_dload_values = 0
                                distributedload_list(n_distributedloads)%subroutine_parameter_number = &
                                    find_name_index(strpar(4),lenstr(4),subroutineparameter_namelist,n_subroutine_parameters)
                            endif
                        endif
                    end do

                else if (strcmp(strpar(1), 'CONSTRA', 7) ) then
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if ( strcmp(strpar(1), 'CONSTRAINNODEPA',15) ) then
                            if (ityp(3)==0) then
                                n_constraints = n_constraints + 1
                                constraint_list(n_constraints)%flag = 1
                                read(strpar(3),*) constraint_list(n_constraints)%node1
                                read(strpar(4),*) constraint_list(n_constraints)%dof1
                                read(strpar(5),*) constraint_list(n_constraints)%node2
                                read(strpar(6),*) constraint_list(n_constraints)%dof2
                            else
                                nodeset1 = find_name_index(strpar(3),lenstr(3),nodeset_namelist,n_nodesets)
                                nodeset2 = find_name_index(strpar(5),lenstr(5),nodeset_namelist,n_nodesets)
                                iof1 = nodeset_list(nodeset1)%index
                                iof2 = nodeset_list(nodeset2)%index
                                read(strpar(2), *) k
                                if (nodeset_list(nodeset1)%n_nodes /= k.or. nodeset_list(nodeset2)%n_nodes /=k) then
                                    write(IOW,*) ' *** Error in input file *** '
                                    write(IOW,*) ' The node sets in a CONSTRAIN NODE PAIRS constraint do',&
                                        ' not have the correct number of nodes'
                                    write(IOW,*) ' Master set ',strpar(2)
                                    write(IOW,*) ' Slave set ',strpar(3)
                                    stop
                                endif
                                iof = 0
                                if (nstr==7.and.ityp(7)==2) then
                                    iof =  find_name_index(strpar(7),lenstr(7),constraintparameter_namelist,n_constraint_parameters)
                                endif
                                do k = 1,nodeset_list(nodeset1)%n_nodes
                                    n_constraints = n_constraints + 1
                                    constraint_list(n_constraints)%flag = 1
                                    constraint_list(n_constraints)%node1 = node_lists(iof1+k-1)
                                    read(strpar(4),*) constraint_list(n_constraints)%dof1
                                    constraint_list(n_constraints)%node2 = node_lists(iof2+k-1)
                                    read(strpar(6),*) constraint_list(n_constraints)%dof2
                                    constraint_list(n_constraints)%index_parameters = iof
                                end do
                            endif
 
                        else if (strcmp(strpar(1), 'TIENO',5) ) then
                            read(strpar(2),*) k
                            nodeset1 = find_name_index(strpar(3),lenstr(3),nodeset_namelist,n_nodesets)
                            if (nodeset_list(nodeset1)%n_nodes /=k) then
                                write(IOW,*) ' *** Error in input file ***'
                                write(IOW,*) ' The node set in a TIE NODES constraint does not have the correct number of nodes'
                                write(IOW,*) ' Node set name ',strpar(3)
                                write(IOW,*) ' Expecting ',k,' nodes'
                                stop
                            endif
                            iof1 = nodeset_list(nodeset1)%index
                            iof = 0
                            if (nstr==7.and.ityp(7)==2) then
                                iof =  find_name_index(strpar(7),lenstr(7),constraintparameter_namelist,n_constraint_parameters)
                            endif
                            do k = 1,nodeset_list(nodeset1)%n_nodes
                                n_constraints = n_constraints + 1
                                constraint_list(n_constraints)%flag = 2
                                constraint_list(n_constraints)%node1 = node_lists(iof1+k-1)
                                read(strpar(4),*) constraint_list(n_constraints)%dof1
                                read(strpar(5),*) constraint_list(n_constraints)%node2
                                read(strpar(6),*) constraint_list(n_constraints)%dof2
                                constraint_list(n_constraints)%index_parameters = iof
                            end do

                        else if (strcmp(strpar(1), 'CONSTRAINNODESET',16) ) then
                            n_constraints = n_constraints + 1
                            if (nstr>3.and.ityp(2)==0.and.ityp(3)==2.and.ityp(4)==2) then
                                read(strpar(2),*) k
                                nodeset1 = find_name_index(strpar(3),lenstr(3),nodeset_namelist,n_nodesets)
                                nodeset2 = find_name_index(strpar(4),lenstr(4),nodeset_namelist,n_nodesets)
                                if (nodeset_list(nodeset1)%n_nodes /=k .or. nodeset_list(nodeset2)%n_nodes /=k) then
                                    write(IOW,*) ' *** Error in input file *** '
                                    write(IOW,*) ' The node set or DOF list in a CONSTRAIN NODE SET constraint does not',&
                                        ' have the correct number of nodes or dofs'
                                    write(IOW,*) ' Node set names ',strpar(3),strpar(4)
                                    stop
                                endif
                                iof1 = nodeset_list(nodeset1)%index
                                iof2 = nodeset_list(nodeset2)%index
                                constraint_list(n_constraints)%flag = 3
                                constraint_list(n_constraints)%node1 = nodeset1
                                constraint_list(n_constraints)%node2 = nodeset2
                                constraint_list(n_constraints)%dof1 = 0
                                constraint_list(n_constraints)%dof2 = 0
                            else
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting a node set name, a dof set name (generate using a node set),',&
                                    ' parameter name (optional) for a CONSTRAIN NODE SET '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            if (nstr==5) then
                                if (ityp(5)==2) then
                                    constraint_list(n_constraints)%index_parameters = &
                                        find_name_index(strpar(5),lenstr(5),constraintparameter_namelist,n_constraint_parameters)
                                else
                                    write(IOW,*) ' *** Error in input file ***'
                                    write(IOW,*) ' Expecting name of a parameter list in a CONSTRAIN NODE SET constraint'
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif
                            else
                                constraint_list(n_constraints)%index_parameters =0
                            endif
                        else if ( strcmp(strpar(1), 'ENDCONS', 7) ) then
                            exit
                        else
                            write(IOW,*) ' *** Error in input file ***'
                            write(IOW,*) ' CONSTRAINTS block was not terminated by an END CONSTRAINTS'
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                        endif

                    end do

                else if ( strcmp(strpar(1), 'ENDBOUN', 7) ) then
                    exit
                else
                    write(IOW,*) ' *** Error in input file *** '
                    write(IOW,*) ' BOUNDARY CONDITIONS key must be terminated with an END BOUNDARY CONDITIONS'
                    write(IOW,*) ' Found '
                    write(IOW,*) strin
                    stop
                endif
            end do

        else if (strcmp(strpar(1), 'TIME',4) ) then
            if (ityp(3)<2) then
                if (strcmp(strpar(2),'INCREM',6) ) then
                    read(strpar(3),*) DTIME
                else if (strcmp(strpar(2),'VAL',3) ) then
                    read(strpar(3),*) TIME
                else
                    write(IOW,'(A)') ' *** Error in input file *** '
                    write(IOW,'(A)') ' Expecting VALUE or INCREMENT following TIME key '
                    write(IOW,'(A)') ' Found '
                    write(IOW,*) strin
                    stop
                endif
            else
                write(IOW,'(A)') ' *** Error in input file ***'
                write(IOW,'(A)') ' Expecting value for time increment following TIME key'
                write(IOW,'(A)') ' Found'
                write(IOW,*) strin
                stop
            endif
        else if (strcmp(strpar(1), 'PRINTINI', 8) ) then
     
            printinitialmesh = .true.

            if (nstr>1.and.ityp(2)==2) then
                read(strpar(2),'(A)') initial_mesh_print_filename
            else
                write(IOW,*) ' *** Error in input file *** '
                write(IOW,*) ' Expecting a file name following PRINT INITIAL MESH key'
                write(IOW,*) ' Found '
                write(IOW,*) strin
                stop
            endif
                
            call files(initial_mesh_print_filename,1,iun)
            initial_mesh_print_unit = iun


        else if (strcmp(strpar(1), 'CHECKMATER', 10) ) then
            checktangent = .true.
            if (nstr<2.or.nstr>3) then
                write(IOW, *) ' *** Error in input file *** '
                write(IOW,*) ' Expecting an material name following CHECK MATERIAL TANGENT key '
                write(IOW,*) ' Found '
                write(IOW,*) strin
                stop
            endif

            checktangent_materialno = find_name_index(strpar(2),lenstr(2),material_namelist,n_materials)

        else if (strcmp(strpar(1), 'CHECKSTIFF', 10) ) then
            checkstiffness = .true.
            if (nstr<2.or.nstr>3) then
                write(IOW, *) ' *** Error in input file *** '
                write(IOW,*) ' Expecting an element number following CHECK STIFFNESS key '
                write(IOW,*) ' Found '
                write(IOW,*) strin
                stop
            endif

            if (ityp(2)>0) then
                if (strcmp(strpar(2),'U',1)) then  ! ABAQUS UEL
                    if (abaqusformat) then
                        uel_string = strpar(2)
                        read(uel_string(2:100), *) checkstiffness_elementno
                        checkstiffness_elementno = checkstiffness_elementno + 99999
                    else
                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                        write(IOW,*) ' CHECK STIFFNESS has specified an ABAQUS UEL, but no UEL has been created '
                        stop
                    endif
                else if (strcmp(strpar(2),'CPE',3)) then               ! 3D continuum element with UMAT
                    checkstiffness_elementno = 10002
                else if (strcmp(strpar(2),'C3D',3)) then               ! 3D continuum element with UMAT
                    checkstiffness_elementno = 10003
                else
                    write(IOW,'(A)')  ' *** ERROR DETECTED IN INPUT FILE *** '
                    write(IOW,'(A)')  ' Element identifier must be an integer, or Un (for ABAQUS UEL) '
                    write(IOW,'(A)')  ' Found'
                    write(IOW,*) strin
                    stop
                endif
            else

                read(strpar(2),*) checkstiffness_elementno
            endif

        else if (strcmp(strpar(1), 'STATICSTEP',10) ) then
            staticstep = .true.
         
            do while (.true.)
                iblnk =1
                read (IOR, 99001, ERR = 200, end = 500) strin
                if (echo) write(IOW,*) strin
                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                if (iblnk==1) cycle
                if (strcmp(strpar(1),'ENDSTA',6) ) then
                    exit
                else if (strcmp(strpar(1),'INITIALTIMESTEP',15) ) then
                    if (nstr==2.and.ityp(2)<2) then
                        read(strpar(2),*) timestep_initial
                        if (timestep_max==0.d0) timestep_max = timestep_initial
                        if (timestep_min==0.d0) timestep_min = timestep_initial
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for an INITIAL TIME STEP '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'MAXTIMESTEP',11) ) then
                    if (nstr==2.and.ityp(2)<2) then
                        read(strpar(2),*) timestep_max
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for MAX TIME STEP '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'MINTIMESTEP',11) ) then
                    if (nstr==2.and.ityp(2)<2) then
                        read(strpar(2),*) timestep_min
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for MIN TIME STEP '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'MAXN',4) ) then
                    if (nstr==2.and.ityp(2)==0) then
                        read(strpar(2),*) max_no_steps
                        if (max_total_time==0.d0) max_total_time= timestep_max*max_no_steps
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for MAX NO STEPS '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'STOPT',5) ) then
                    if (nstr==2.and.ityp(2)<2) then
                        read(strpar(2),*) max_total_time
                        if(max_no_steps == 0) max_no_steps = int(max_total_time/timestep_min)
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for MAX TOTAL TIME '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'STATEPRINTS',11) ) then
                    if (nstr==2.and.ityp(2)==0) then
                        read(strpar(2),*) state_print_steps
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for STATE PRINT STEPS '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'USERPRINTST',11) ) then
                    if (nstr==2.and.ityp(2)==0) then
                        read(strpar(2),*) user_print_steps
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for USER PRINT STEPS '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'SOLVER',6)) then
                    unsymmetric_stiffness = .false.
                    solvertype = 1
                    nonlinear = .false.
                    if (strcmp(strpar(2),'DIREC',5) ) then
                        solvertype = 1
                    else if (strcmp(strpar(2),'CONJU',5) ) then
                        solvertype = 2
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Unknown solver type ',strpar(2)
                        stop
                    endif
                    if (strcmp(strpar(3),'LINEAR',6) ) then
                        nonlinear = .false.
                        max_newton_iterations = 1
                        newtonraphson_tolerance = 1.d-08
                    else if (strcmp(strpar(3),'NONLINEAR',9) ) then
                        nonlinear = .true.
                        if (nstr>=5.and.ityp(4)<2.and.ityp(5)==0) then
                            read(strpar(4),*) newtonraphson_tolerance
                            read(strpar(5),*) max_newton_iterations
                        else
                            write(IOW,*) ' *** Error in input file ***'
                            write(IOW,*) ' Expecting convergence tolerance and max # iterations following NONLINEAR solver option'
                            write(IOW,*) ' Found'
                            write(IOW,*) strin
                            stop
                        endif
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Unknown solver option ',strpar(3)
                        stop
                    endif
                    if (strcmp(strpar(nstr), 'UNSYMM', 6) ) then
                        unsymmetric_stiffness = .true.
                    endif
                else if (strcmp(strpar(1),'PRINTST',7) ) then
                    if (nstr==2.and.ityp(2)==2) then
                        stateprint = .true.
                        read(strpar(2),'(A)') state_print_filename
                    else
                        write(IOW,*) ' *** Error in input file *** '
                        write(IOW,*) ' Expecting a file name following PRINT STATE key'
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                
                    call files(state_print_filename,1,iun)
                    state_print_unit = iun
                
                
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                
                        if (strcmp(strpar(1),'ENDPRINTS',9) ) then
                      
                            exit
                        endif
                   
                        if (ityp(1)<2) then
                            write(IOW,*) ' *** Error in input file ***'
                            write(IOW,*) ' Unrecognized option in PRINT STATE in a STATIC STEP '
                            write(IOW,*) ' Options are: Zones, degrees of freedom,'
                            write(IOW,*) ' field variables, displaced mesh, displacement scale factor '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                   
                        if (strcmp(strpar(1),'DEGREE',6) ) print_dof = .true.
                        if (strcmp(strpar(1),'DISPLACEDM',10) ) print_displacedmesh = .true.
                        if (strcmp(strpar(1),'DISPLACEMENTS',13) ) then
                            if (nstr<2.or.ityp(2)==2) then
                                write(IOW,*) ' *** Error in input file ***'
                                write(IOW,*) ' Expecting value for displacement scale factor '
                                write(IOW,*) ' Found'
                                write(IOW,*) strin
                                stop
                            endif
                            read(strpar(2),*) displacementscalefactor
                        endif
                        if (strcmp(strpar(1),'FIELD',5) ) then
                            if (nstr<2) then
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting names of field variables '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        
                            do k = 2,nstr
                                n_field_variables = n_field_variables + 1
                                read(strpar(k),*) field_variable_names(k-1)
                            end do
                        endif
                        if (strcmp(strpar(1), 'ZONES',5) ) then
                            zone_print_flag = .false.
                            combinezones = .false.
                            if (nstr==2 .and. ityp(2)==2) then
                                if (strcmp(strpar(2),'COMBINE',7) ) then
                                    combinezones = .true.
                                else
                                    write(IOW,*) ' *** Error in input file *** '
                                    write(IOW,*) ' Unrecognized option when specifying ZONES in a STATE PRINT block'
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif
                            endif
                            do while (.true.)
                                iblnk =1
                                read (IOR, 99001, ERR = 200, end = 500) strin
                                if (echo) write(IOW,*) strin
                                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                if (iblnk==1) cycle
                         
                                if (strcmp(strpar(1),'ENDZ',4) ) exit

                                do k = 1,nstr
                                    izone = find_name_index(strpar(k),lenstr(k),zone_namelist,n_zones)
                                    zone_print_flag(izone) = .true.
                                end do
                            end do
                        endif
              
                    end do
                else if (strcmp(strpar(1),'USERPRINTPA',11) ) then
                    userprint = .true.
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDU',4) ) then
                            exit
                        endif
                        do k = 1,nstr
                            if (ityp(k)<2) then
                                n_user_print_parameters = n_user_print_parameters + 1
                                read(strpar(k),*) user_print_parameters(n_user_print_parameters)
                            else
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting list of numbers following USER PRINT PARAMETERS '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        end do
                    end do
                else if (strcmp(strpar(1),'USERPRINTFI',11) ) then
                    userprint = .true.
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDU',4) ) then
                            exit
                        endif
                        if (ityp(1)==2) then
                            n_user_print_files = n_user_print_files + 1
                            read(strpar(1),'(A)') user_print_filenames(n_user_print_files)
                            call files(user_print_filenames(n_user_print_files),1,iun)
                            user_print_units(n_user_print_files) = iun
                        else
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' Expecting list of file names following USER PRINT FILES '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                    end do
                endif
            end do
        else if (strcmp(strpar(1), 'EXPLICITDYNAMICSTEP',19) ) then
            explicitdynamicstep = .true.
            use_lumped_projection_matrix = .true.

            do while (.true.)
                iblnk =1
                read (IOR, 99001, ERR = 200, end = 500) strin
                if (echo) write(IOW,*) strin
                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                if (iblnk==1) cycle
                if (strcmp(strpar(1),'ENDEXP',6) ) then
                    if (length_densities==0) then
                        if (.not.abaqusformat) then
                        write(IOW,'(A)') ' *** Error in input file *** '
                        write(IOW,'(A)') '     Densities must be specified for an explicit dynamic analysis '
                        stop
                        endif
                    endif
                    if (dynamic_timestep==0.d0) then
                        write(IOW,'(A)') '  ** Error in input file *** '
                        write(IOW,'(A)') ' A value must be provided for the TIME STEP in an explicit dynamic analysis '
                        stop
                    endif
                    if (no_dynamic_steps==0) then
                        write(IOW,'(A)') '  ** Error in input file *** '
                        write(IOW,'(A)') ' A value must be provided for the NUMBER OF STEPS in an explicit dynamic analysis '
                        stop
                    endif
                    total_dynamic_time = dynamic_timestep*no_dynamic_steps
                    exit
                else if (strcmp(strpar(1),'TIMESTEP',8) ) then
                    if (nstr==2.and.ityp(2)<2) then
                        read(strpar(2),*) dynamic_timestep
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for an TIME STEP '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'NUMB',4) ) then
                    if (nstr==2.and.ityp(2)==0) then
                        read(strpar(2),*) no_dynamic_steps
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for MAX NO STEPS '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'STATEPRINTS',11) ) then
                    if (nstr==2.and.ityp(2)==0) then
                        read(strpar(2),*) state_print_steps
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for STATE PRINT STEPS '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'USERPRINTST',11) ) then
                    if (nstr==2.and.ityp(2)==0) then
                        read(strpar(2),*) user_print_steps
                    else
                        write(IOW,*) ' *** Error in input file ***'
                        write(IOW,*) ' Expecting value for USER PRINT STEPS '
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                else if (strcmp(strpar(1),'PRINTST',7) ) then
                    if (nstr==2.and.ityp(2)==2) then
                        stateprint = .true.

                        read(strpar(2),'(A)') state_print_filename
                    else
                        write(IOW,*) ' *** Error in input file *** '
                        write(IOW,*) ' Expecting a file name following PRINT STATE key'
                        write(IOW,*) ' Found '
                        write(IOW,*) strin
                        stop
                    endif
                
                    call files(state_print_filename,1,iun)
                    state_print_unit = iun
                
                
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                
                        if (strcmp(strpar(1),'ENDPRINTS',9) ) then
                      
                            exit
                        endif
                   
                        if (ityp(1)<2) then
                            write(IOW,*) ' *** Error in input file ***'
                            write(IOW,*) ' Unrecognized option in PRINT STATE in a STATIC STEP '
                            write(IOW,*) ' Options are: Zones, degrees of freedom,'
                            write(IOW,*) ' field variables, displaced mesh, displacement scale factor '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                   
                        if (strcmp(strpar(1),'DEGREE',6) ) print_dof = .true.
                        if (strcmp(strpar(1),'DISPLACEDM',10) ) print_displacedmesh = .true.
                        if (strcmp(strpar(1),'DISPLACEMENTS',13) ) then
                            if (nstr<2.or.ityp(2)==2) then
                                write(IOW,*) ' *** Error in input file ***'
                                write(IOW,*) ' Expecting value for displacement scale factor '
                                write(IOW,*) ' Found'
                                write(IOW,*) strin
                                stop
                            endif
                            read(strpar(2),*) displacementscalefactor
                        endif
                        if (strcmp(strpar(1),'FIELD',5) ) then
                            if (nstr<2) then
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting names of field variables '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        
                            do k = 2,nstr
                                n_field_variables = n_field_variables + 1
                                read(strpar(k),*) field_variable_names(k-1)
                            end do
                        endif
                        if (strcmp(strpar(1), 'ZONES',5) ) then
                            zone_print_flag = .false.
                            combinezones = .false.
                            if (nstr==2 .and. ityp(2)==2) then
                                if (strcmp(strpar(2),'COMBINE',7) ) then
                                    combinezones = .true.
                                else
                                    write(IOW,*) ' *** Error in input file *** '
                                    write(IOW,*) ' Unrecognized option when specifying ZONES in a STATE PRINT block'
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif
                            endif
                            do while (.true.)
                                iblnk =1
                                read (IOR, 99001, ERR = 200, end = 500) strin
                                if (echo) write(IOW,*) strin
                                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                if (iblnk==1) cycle
                         
                                if (strcmp(strpar(1),'ENDZ',4) ) exit

                                do k = 1,nstr
                                    izone = find_name_index(strpar(k),lenstr(k),zone_namelist,n_zones)
                                    zone_print_flag(izone) = .true.
                                end do
                            end do
                        endif
              
                    end do
                else if (strcmp(strpar(1),'USERPRINTPA',11) ) then
                    userprint = .true.
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDU',4) ) then
                            exit
                        endif
                        do k = 1,nstr
                            if (ityp(k)<2) then
                                n_user_print_parameters = n_user_print_parameters + 1
                                read(strpar(k),*) user_print_parameters(n_user_print_parameters)
                            else
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting list of numbers following USER PRINT PARAMETERS '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        end do
                    end do
                else if (strcmp(strpar(1),'USERPRINTFI',11) ) then
                    userprint = .true.
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 200, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDU',4) ) then
                            exit
                        endif
                        if (ityp(1)==2) then
                            n_user_print_files = n_user_print_files + 1
                            read(strpar(1),'(A)') user_print_filenames(n_user_print_files)
                            call files(user_print_filenames(n_user_print_files),1,iun)
                            user_print_units(n_user_print_files) = iun
                        else
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' Expecting list of file names following USER PRINT FILES '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                    end do
                endif
            end do


        else if (strcmp(strpar(1),'STOP',4) ) then
      
            do n = 1,n_nodes              ! Check that # of displacements in a displacement map is consistent with # coords
                if (node_list(n)%n_displacements>0) then
                    if (node_list(n)%n_displacements /= node_list(n)%n_coords) then
                        write(IOW,*) ' *** Error in input file *** '
                        write(IOW,*) ' The number of nodes in a displacement map is inconsistent with the number of coordinates '
                        stop
                    endif
                endif
            end do
  
  
            ! Determine storage necessary for user elements
            length_node_array = 0
            length_coord_array = 0
            length_dof_array = 0
            length_property_array = 0
            length_state_variable_array = 0
            do lmn = 1,n_elements
                if (element_list(lmn)%n_nodes>length_node_array) length_node_array = element_list(lmn)%n_nodes
                if (element_list(lmn)%n_element_properties>length_property_array) then
                    length_property_array = element_list(lmn)%n_element_properties
                endif
                if (element_list(lmn)%n_states>length_state_variable_array) length_state_variable_array=element_list(lmn)%n_states
                ix = 0
                iu = 0
                do j = 1,element_list(lmn)%n_nodes
                    n = connectivity(element_list(lmn)%connect_index+j-1)
                    ix = ix + node_list(n)%n_coords
                    iu = iu + node_list(n)%n_dof
                end do
                if (ix>length_coord_array) length_coord_array = ix
                if (iu>length_dof_array)   length_dof_array = iu
            end do
       
            do nc = 1,n_constraints
                if (constraint_list(nc)%flag<3) then
                    if (node_list( constraint_list(nc)%node1 )%n_dof > length_dof_array) then
                        length_dof_array = node_list( constraint_list(nc)%node1 )%n_dof
                    endif
                    if (node_list( constraint_list(nc)%node2 )%n_dof > length_dof_array) then
                        length_dof_array = node_list( constraint_list(nc)%node1 )%n_dof
                    endif
                    if (node_list( constraint_list(nc)%node1 )%n_coords > length_coord_array) then
                        length_coord_array = node_list( constraint_list(nc)%node1 )%n_coords
                    endif
                    if (node_list( constraint_list(nc)%node2 )%n_coords > length_coord_array) then
                        length_coord_array = node_list( constraint_list(nc)%node1 )%n_coords
                    endif
                else if (constraint_list(nc)%flag==3) then
                    nset = constraint_list(nc)%node1
                    ix = 0
                    iu = 0
                    do j = 1, nodeset_list(nset)%n_nodes
                        n = node_lists(nodeset_list(nset)%index + j - 1)
                        do k = 1, node_list(n)%n_coords
                            ix = ix + 1
                        end do
                        do k = 1, node_list(n)%n_dof
                            iu = iu + 1
                        end do
                        if (iu>length_dof_array)   length_dof_array = iu
                        if (ix>length_coord_array) length_coord_array = ix
                    end do
                endif
            end do
     
            if (stateprint) then             ! Check printing for consistency
                do iz = 1,n_zones
                    if (.not.zone_print_flag(iz)) cycle
                    do lmn = zone_list(iz)%start_element,zone_list(iz)%end_element
                        do k = 1,element_list(lmn)%n_nodes
                            n = connectivity(element_list(lmn)%connect_index+k-1)
                            if (node_list(n)%n_coords==2) then
                                if (zone_dimension(iz)==0) then
                                    zone_dimension(iz) = 2
                                else if (zone_dimension(iz) /=2) then
                                    write(IOW,*) ' *** Warning *** '
                                    write(IOW,*) ' A zone to be printed using a STATE PRINT command contains nodes with'
                                    write(IOW,*) ' an inconsistent number of coordinates '
                                    write(IOW,*) ' This zone will not be printed '
                                    zone_print_flag(iz) = .false.
                                    exit
                                endif

                            else if (node_list(n)%n_coords==3) then
                                if (zone_dimension(iz)==0) then
                                    zone_dimension(iz) = 3
                                else if (zone_dimension(iz) /=3) then
                                    write(IOW,*) ' *** Warning *** '
                                    write(IOW,*) ' A zone to be printed using a STATE PRINT command contains nodes with'
                                    write(IOW,*) ' an inconsistent number of coordinates '
                                    write(IOW,*) ' This zone will not be printed '
                                    zone_print_flag(iz) = .false.
                                    exit
                                endif
                            endif
                            if (zone_ndof(iz)==0) then
                                zone_ndof(iz) = node_list(n)%n_dof
                            else if (zone_dimension(iz) /=node_list(n)%n_dof) then
                                write(IOW,*) ' *** Warning *** '
                                write(IOW,*) ' A zone to be printed using a STATE PRINT command contains nodes with'
                                write(IOW,*) ' an inconsistent number of DOF '
                                write(IOW,*) ' This zone will not be printed '
                                zone_print_flag(iz) = .false.
                                exit
                            endif
                        end do
                        if (.not.zone_print_flag(iz)) exit
                        if (zone_dimension(iz)==2) then
                            if (element_list(lmn)%n_nodes == 3) cycle
                            if (element_list(lmn)%n_nodes == 4) cycle
                            if (element_list(lmn)%n_nodes == 6) cycle
                            if (element_list(lmn)%n_nodes == 8) cycle
                            if (element_list(lmn)%n_nodes == 9) cycle
                        else if (zone_dimension(iz)==3) then
                            if (element_list(lmn)%n_nodes == 4) cycle
                            if (element_list(lmn)%n_nodes == 8) cycle
                        endif
                        write(IOW,*) ' *** Warning *** '
                        write(IOW,*) ' A zone to be printed using a STATE PRINT command contains elements with'
                        write(IOW,*) ' a number of nodes that cannot be printed'
                        write(IOW,*) ' This zone will not be printed '
                        zone_print_flag(iz) = .false.
                    end do
                end do

                n = 0
                do iz = 1,n_zones
                    if (zone_dimension(iz)/=3.and.zone_dimension(iz)/=2) zone_print_flag(iz) = .false.
                    if (zone_print_flag(iz)) n = n + 1
                end do
                if (n==0) then
                    write(IOW,*) ' *** Warning *** '
                    write(IOW,*) ' No printable zones were found '
                    write(IOW,*) ' STATE PRINT has been suppressed '
                    stateprint = .false.
                endif
                if (combinezones) then
                    k = 0
                    n = 0
                    do iz = 2,n_zones
                        if (zone_print_flag(iz)) then
                            if (k==0) then
                                k = zone_dimension(iz)
                            else if (k/=zone_dimension(iz)) then
                                write(IOW,*) ' *** Warning *** '
                                write(IOW,*) ' Zones in STATE PRINT option have different number of coordinats'
                                write(IOW,*) ' They cannot be combined '
                                combinezones = .false.
                            endif
                            if (n==0) then
                                n = zone_ndof(iz)
                            else if (n/=zone_ndof(iz)) then
                                write(IOW,*) ' *** Warning *** '
                                write(IOW,*) ' Zones in STATE PRINT option have different number of DOF'
                                write(IOW,*) ' They cannot be combined '
                                combinezones = .false.
                            endif
 
                        endif
                    end do
                endif
            endif
   
            write(IOW,'(//(A))')   ' The input file has been read successfully '
            write(IOW,'(A31,I6)')  ' Number of elements:           ',n_elements
            write(IOW,'(A31,I6)')  ' Number of nodes:              ',n_nodes
            write(IOW,'(A31,I6)')  ' Number of constraints:        ',n_constraints
            write(IOW,'(A31,I6)')  ' Total no. degrees of freedom: ',length_dofs
   
   
            return
        else
            write(IOW,*) ' Error in input file '
            write(IOW,*) ' Unrecognized keyword '
            write(IOW,*) strin
            stop
        endif
    end do

    do n = 1,n_nodes              ! Check that # of displacements in a displacement map is consistent with # coords
        if (node_list(n)%n_displacements>0) then
            if (node_list(n)%n_displacements /= node_list(n)%n_coords) then
                write(IOW,*) ' *** Error in input file *** '
                write(IOW,*) ' The number of nodes in a displacement map is inconsistent with the number of coordinates '
                stop
            endif
        endif
    end do

    ! Determine storage necessary for user elements
    length_node_array = 0
    length_coord_array = 0
    length_dof_array = 0
    do lmn = 1,n_elements
        if (element_list(lmn)%n_nodes>length_node_array) length_node_array = element_list(lmn)%n_nodes
        ix = 0
        iu = 0
        do j = 1,element_list(lmn)%n_nodes
            n = connectivity(element_list(lmn)%connect_index+j-1)
            ix = ix + node_list(n)%n_coords
            iu = iu + node_list(n)%n_dof
        end do
        if (ix>length_coord_array) length_coord_array = ix
        if (iu>length_dof_array)   length_dof_array = iu
    end do
    
    do nc = 1,n_constraints
        if (constraint_list(nc)%flag<3) then
            if (node_list( constraint_list(nc)%node1 )%n_dof > length_dof_array) then
                length_dof_array = node_list( constraint_list(nc)%node1 )%n_dof
            endif
            if (node_list( constraint_list(nc)%node2 )%n_dof > length_dof_array) then
                length_dof_array = node_list( constraint_list(nc)%node1 )%n_dof
            endif
            if (node_list( constraint_list(nc)%node1 )%n_coords > length_coord_array) then
                length_coord_array = node_list( constraint_list(nc)%node1 )%n_coords
            endif
            if (node_list( constraint_list(nc)%node2 )%n_coords > length_coord_array) then
                length_coord_array = node_list( constraint_list(nc)%node1 )%n_coords
            endif
        else if (constraint_list(nc)%flag==3) then
            nset = constraint_list(nc)%node1
            ix = 0
            iu = 0
            do j = 1, nodeset_list(nset)%n_nodes
                n = node_lists(nodeset_list(nset)%index + j - 1)
                do k = 1, node_list(n)%n_coords
                    ix = ix + 1
                end do
                do k = 1, node_list(n)%n_dof
                    iu = iu + 1
                end do
                if (iu>length_dof_array)   length_dof_array = iu
                if (ix>length_coord_array) length_coord_array = ix
            end do
        endif
    end do

    if (stateprint) then             ! Check printing for consistency
        do iz = 1,n_zones
            if (.not.zone_print_flag(iz)) cycle
            do lmn = zone_list(iz)%start_element,zone_list(iz)%end_element
                do k = 1,element_list(lmn)%n_nodes
                    n = connectivity(element_list(lmn)%connect_index+k-1)
                    if (node_list(n)%n_coords==2) then
                        if (zone_dimension(iz)==0) then
                            zone_dimension(iz) = 2
                        else if (zone_dimension(iz) /=2) then
                            write(IOW,*) ' *** Warning *** '
                            write(IOW,*) ' A zone to be printed using a STATE PRINT command contains nodes with'
                            write(IOW,*) ' an inconsistent number of coordinates '
                            write(IOW,*) ' This zone will not be printed '
                            zone_print_flag(iz) = .false.
                            exit
                        endif
                    else if (node_list(n)%n_coords==3) then
                        if (zone_dimension(iz)==0) then
                            zone_dimension(iz) = 3
                        else if (zone_dimension(iz) /=3) then
                            write(IOW,*) ' *** Warning *** '
                            write(IOW,*) ' A zone to be printed using a STATE PRINT command contains nodes with'
                            write(IOW,*) ' an inconsistent number of coordinates '
                            write(IOW,*) ' This zone will not be printed '
                            zone_print_flag(iz) = .false.
                            exit
                        endif
                    endif
                    if (zone_ndof(iz)==0) then
                        zone_ndof(iz) = node_list(n)%n_dof
                    else if (zone_dimension(iz) /=node_list(n)%n_dof) then
                        write(IOW,*) ' *** Warning *** '
                        write(IOW,*) ' A zone to be printed using a STATE PRINT command contains nodes with'
                        write(IOW,*) ' an inconsistent number of DOF '
                        write(IOW,*) ' This zone will not be printed '
                        zone_print_flag(iz) = .false.
                        exit
                    endif
                end do
                if (.not.zone_print_flag(iz)) exit
                if (zone_dimension(iz)==2) then
                    if (element_list(lmn)%n_nodes == 3) cycle
                    if (element_list(lmn)%n_nodes == 4) cycle
                    if (element_list(lmn)%n_nodes == 6) cycle
                    if (element_list(lmn)%n_nodes == 8) cycle
                    if (element_list(lmn)%n_nodes == 9) cycle
                else if (zone_dimension(iz)==3) then
                    if (element_list(lmn)%n_nodes == 4) cycle
                    if (element_list(lmn)%n_nodes == 8) cycle
                endif
                write(IOW,*) ' *** Warning *** '
                write(IOW,*) ' A zone to be printed using a STATE PRINT command contains elements with'
                write(IOW,*) ' a number of nodes that cannot be printed'
                write(IOW,*) ' This zone will not be printed '
                zone_print_flag(iz) = .false.
            end do
        end do

        n = 0
        do iz = 1,n_zones
            if (zone_dimension(iz)/=3.and.zone_dimension(iz)/=2) zone_print_flag(iz) = .false.
            if (zone_print_flag(iz)) n = n + 1
        end do
        if (n==0) then
            write(IOW,*) ' *** Warning *** '
            write(IOW,*) ' No printable zones were found '
            write(IOW,*) ' STATE PRINT has been suppressed '
            stateprint = .false.
        endif
        if (combinezones) then
            k = 0
            n = 0
            do iz = 2,n_zones
                if (zone_print_flag(iz)) then
                    if (k==0) then
                        k = zone_dimension(iz)
                    else if (k/=zone_dimension(iz)) then
                        write(IOW,*) ' *** Warning *** '
                        write(IOW,*) ' Zones in STATE PRINT option have different number of coordinats'
                        write(IOW,*) ' They cannot be combined '
                        combinezones = .false.
                    endif
                    if (n==0) then
                        n = zone_ndof(iz)
                    else if (n/=zone_ndof(iz)) then
                        write(IOW,*) ' *** Warning *** '
                        write(IOW,*) ' Zones in STATE PRINT option have different number of DOF'
                        write(IOW,*) ' They cannot be combined '
                        combinezones = .false.
                    endif
 
                endif
            end do
        endif
    
       
    endif
        
    write(IOW,'(//(A))')   ' The input file has been read successfully '
    write(IOW,'(A31,I6)')  ' Number of elements:           ',n_elements
    write(IOW,'(A31,I6)')  ' Number of nodes:              ',n_nodes
    write(IOW,'(A31,I6)')  ' Number of constraints:        ',n_constraints
    write(IOW,'(A31,I6)')  ' Total no. degrees of freedom: ',length_dofs
   

500 write(IOW,*) ' *** Warning - end of input file was reached before a STOP keyword '
    return

200 write(IOW,*) ' *** An input file error occurred in subroutine readmesh *** '
    stop
99001 format (A100)

end subroutine read_input_file
!
! ==================================================== allocate_mesh_storage ====================================
!
subroutine allocate_mesh_storage
    use Types
    use ParamIO
    use Globals
    use Mesh
    use Boundaryconditions
    use Staticstepparameters
    use Printparameters
    use Fileparameters
    implicit none
  
    integer :: ncoor         ! # coords for current node
    integer :: ndof          ! # degrees of freedom for current node
    integer :: k             ! Counter
    integer :: status        ! status flag for memory allocation
    integer :: nnod          ! no. nodes on an element
    integer :: nsvar         ! no. state vars on an element
    integer :: n_int_pts     ! No. integration points on an internal continuum element
    integer :: n             ! Number of created nodes
    integer :: start_node    ! First node in a generated node set
    integer :: end_node      ! Last node in a generated node set
    integer :: increment     ! Increment in a generated node set
    integer :: n_mater_states ! No. state variables specified for a material
    integer :: max_mater_states ! Max no. state vars on any material
    logical :: strcmp

    n_zones = 0
    n_nodes = 0
    n_elements = 0
    n_materials = 0
    length_coords = 0
    length_dofs = 0
    length_displacement_map = 0
    length_element_properties = 0
    length_int_element_properties = 0
    length_material_properties = 0
    length_densities = 0
    length_state_variables = 0
    length_connectivity = 0
    max_mater_states = 0
  
    n_histories = 0
    n_subroutine_parameters = 0
    n_constraint_parameters = 0
    n_nodesets = 0
    n_elementsets = 0
      
    n_prescribeddof = 0
    n_prescribedforces = 0
    n_distributedloads = 0
    n_constraints = 0
  
    length_element_lists = 0
    length_node_lists = 0
    length_history_data = 0
    length_subroutine_parameters = 0
    length_dof_values = 0
    length_dload_values = 0
    length_constraint_parameters = 0
    n_mesh_parameters = 0

    n_user_print_parameters = 0
    n_user_print_files=0
    n_field_variables = 0
  
    n_total_files = 0
 
    do while (.true.)
  
        iblnk=1
        read (IOR, 99001, ERR = 500, end = 500) strin
        if (echo) write(IOW,*) strin
        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
        if (iblnk==1) cycle
  
        if ( strcmp(strpar(1), 'MESH', 4) ) then
 
  
            do while ( .true. )
                !     Read a line of the input file
                iblnk =1
                read (IOR, 99001, ERR = 500, end = 500) strin
                if (echo) write(IOW,*) strin
                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                if (iblnk==1) cycle
                !     ---------------------- SOME PORTION OF A MESH IS GENERATED IN A USER SUBROUTINE ---------------------
                if ( strcmp(strpar(1), 'USERSUB',7) ) then
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDUS',5) ) exit
                        do k = 1,nstr
                            if (ityp(k)<2) then
                                n_mesh_parameters = n_mesh_parameters + 1
                            else
                                write(IOW,'(A)') ' *** Error in input file ***'
                                write(IOW,'(A)') ' Expecting a list of parameters defining user-subroutine generated mesh'
                                write(IOW,'(A)') ' Found '
                                write(IOW,'(A)') strin
                            endif
                        end do
                    end do
                !     ---------------------- READ NODAL DATA ---------------------

                else if ( strcmp(strpar(1), 'NODE', 4) ) then
                    do while ( .true. )
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        !        Set properties to be assigned to next group of nodes
                        if ( strcmp(strpar(1), 'PARAM',5) ) then
                            if ( nstr<4 .or. (ityp(2)/=0) .or. (ityp(3)/=0) .or. (ityp(4) /=0) ) then
                                write (IOW,*) '  *** ERROR DETECTED IN INPUT FILE ***'
                                write (IOW,*) strin
                                write (IOW,*) ' No. coords, no. dof, and an optional integer',&
                                    '                      identifier must be supplied with node keyword'
                                stop
                            end if
                            read (strpar(2), *) ncoor
                            read (strpar(3), *) ndof
                        !        Set mapping of nodal DOF to displacements
                        else if ( strcmp(strpar(1), 'DISPLACEMENTDOF', 15) ) then
                            if (ityp(2)<2) then
                                do k = 1,nstr-1
                                    length_displacement_map = length_displacement_map + 1
                                end do
                            endif
                        else if ( strcmp(strpar(1), 'CREATE', 6 ) ) then                                               ! Create nodes with no coordinates
                            if (nstr/=2 .or. ityp(2) ==2) then
                                write(IOW,*) ' *** Error in input file ***'
                                write(IOW,*) ' Expecting an integer defining number of nodes following a CREATE NODES key '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            read (strpar(2), *) n
                            do k = 1,n
                                n_nodes = n_nodes + 1
                                length_dofs = length_dofs + ndof
                            end do
                        !        Read coords of nodes
                        else if ( strcmp(strpar(1), 'COOR', 4) ) then

                            do while ( .true. )
                                read (IOR, 99001, ERR = 500, end = 500) strin
                                if (echo) write(IOW,*) strin
                                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                if ( iblnk==1 ) cycle
                                if ( strcmp(strpar(1), 'ENDCOOR', 7) ) then
                                    exit
                                endif
                                if ( ityp(1)==2 ) then
                                    write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                    write(IOW,*) ' COORDINATES key must be terminated by an END COORDINATES'
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                endif
                                if ( nstr>ncoor+1 ) then
                                    write(IOW,*) '  *** ERROR DETECTED IN INPUT FILE ***'
                                    write(IOW,*) ' Expecting ', ncoor,' coords of a node.  Found only ',nstr
                                    stop
                                end if
                                n_nodes = n_nodes + 1
                                do k = nstr-ncoor+1, nstr
                                    length_coords = length_coords + 1
                                end do
                                length_dofs = length_dofs + ndof
                            end do
                        else if ( strcmp(strpar(1), 'INIT', 4) ) then
                            !        Finished reading nodes
                            do while(.true.)
                                read (IOR, 99001, ERR = 500, end = 500) strin
                                if (echo) write(IOW,*) strin
                                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                if ( iblnk==1 ) cycle
                                if ( strcmp(strpar(1), 'ENDINIT', 7) ) then
                                    exit
                                elseif (strcmp(strpar(1), 'USERSUBR', 8) ) then
                                    do while (.true.)
                                        read (IOR, 99001, ERR = 500, end = 500) strin
                                        if (echo) write(IOW,*) strin
                                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                        if ( iblnk==1 ) cycle
                                        if (strcmp(strpar(1),'ENDUS',5) ) then
                                            exit
                                        endif
                                        do k = 1,nstr
                                            if (ityp(k)<2) then
                                                n_mesh_parameters = n_mesh_parameters + 1
                                            else
                                                write(IOW,'(A)') ' *** Error in input file ***'
                                                write(IOW,'(A)') ' Expecting a list of parameters for user-subroutine generated DOF'
                                                write(IOW,'(A)') ' Found '
                                                write(IOW,'(A)') strin
                                            endif
                                        end do
                                    end do
                                endif
                            end do

                        else if ( strcmp(strpar(1), 'ENDNODE', 7) ) then
                            exit
                        else
                            write(IOW,*) ' Error detected in input file '
                            write(IOW,*) ' Expecting END NODES keyword'
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
      
                    end do
                else if ( strcmp(strpar(1), 'MATE', 4) ) then

                    if (nstr<2) then
                       write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE ***'
                       write(IOW,*) ' MATERIAL key must specify a name for the material '
                       stop
                    endif

                    n_materials = n_materials + 1

                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle


                        if (strcmp(strpar(1), 'STAT', 4) ) then
                            if (nstr==2.and.ityp(2)==0) then
                               read(strpar(2),*) n_mater_states
                               if (n_mater_states>max_mater_states) then
                                   max_mater_states = n_mater_states
                               endif
                            else
                               write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                               write(IOW,*) ' STATE VARIABLES key in a MATERIAL definition must specify no. state variables '
                               write(IOW,*) ' Found '
                               write(IOW,*) strin
                               stop
                            endif

                        else if (strcmp(strpar(1), 'PROP', 4) ) then

                            do while (.true.)
                                iblnk =1
                                read (IOR, 99001, ERR = 500, end = 500) strin
                                if (echo) write(IOW,*) strin
                                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                if (iblnk==1) cycle
                                if ( strcmp(strpar(1), 'ENDPROP', 7) ) then
                                    exit
                                endif
                                if ( ityp(1)==2 ) then
                                    write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                    write(IOW,*) ' PROPERTIES key in a MATERIAL block must be terminated by an END PROPERTIES'
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif
                                do k = 1,nstr
                                    length_material_properties = length_material_properties+1
                                end do
                            end do

                        else if ( strcmp(strpar(1), 'ENDMATE', 7) ) then
                            exit
                        else
                            write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                            write(IOW,*) ' MATERIAL key must be terminated by an END MATERIAL'
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                    end do

                !   Read element properties
                else if ( strcmp(strpar(1), 'ELEM', 4) ) then
                    if (nstr<2) then
                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                        write(IOW,*) ' The ELEMENT key must be followed by a USER or INTERNAL option'
                        stop
                    endif

                    if (strcmp(strpar(2),'USER',4) ) then

                        do while (.true.)
                            iblnk =1
                            read (IOR, 99001, ERR = 500, end = 500) strin
                            if (echo) write(IOW,*) strin
                            call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                            if (iblnk==1) cycle
                            if ( strcmp(strpar(1), 'PARAM', 5) ) then

                                if ( nstr<4 .or. (ityp(2)/=0) .or. (ityp(3)/=0) ) then
                                    write(IOW,*) '  *** ERROR DETECTED IN INPUT FILE ***'
                                    write(IOW,*)  ' No. nodes per element, no. state vars and element identifier'
                                    write(IOW,*)  ' must be supplied with element keyword'
                                    stop
                                end if

                                read (strpar(2), *) nnod
                                read (strpar(3), *) nsvar
                            else if (strcmp(strpar(1), 'DENS',4) ) then
                                length_densities = length_densities + 1
                            else if (strcmp(strpar(1), 'PROP',4) ) then
                                do while (.true.)
                                    iblnk =1
                                    read (IOR, 99001, ERR = 500, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDPROP', 7) ) then
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' PROPERTIES key must be terminated by an END PROPERTIES'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    do k = 1,nstr
                                        length_element_properties = length_element_properties+1
                                    end do
                                end do
                            else if (strcmp(strpar(1), 'INTEGERPROP',11) ) then
                                do while (.true.)
                                    iblnk =1
                                    read (IOR, 99001, ERR = 500, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDINTE', 7) ) then
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' INTEGER PROPERTIES key must be terminated by an END INTEGER PROPERTIES'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    do k = 1,nstr
                                        length_int_element_properties = length_int_element_properties+1
                                    end do
                                end do
                            else if (strcmp(strpar(1), 'INITIALS',8)) then
                                do while (.true.)
                                    iblnk =1
                                    read (IOR, 99001, ERR = 500, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDINIT', 7) ) then
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' INITIAL STATE VARIABLES key must be', &
                                                     ' terminated by an END INITIAL STATE VARIABLES'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                 end do
                            else if (strcmp(strpar(1), 'CONN',4) ) then
                                n_zones = n_zones + 1
                                do while (.true.)
                                    iblnk =1
                                    read (IOR, 99001, ERR = 500, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDCONN', 7) ) then
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' CONNECTIVITY key must be terminated by an END CONNECTIVITY'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    n_elements = n_elements + 1
                                    length_connectivity = length_connectivity + nnod
                                    length_state_variables = length_state_variables + nsvar
                                end do
                            else if (strcmp(strpar(1), 'TYPE',4) ) then
                                write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                write(IOW,*) ' TYPE key can only be used with INTERNAL elements '
                                write(IOW,*) ' Use the PARAMETERS key instead '
                                stop
                            else if ( strcmp(strpar(1), 'ENDELEM', 7) ) then
                                exit
                            else
                                write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                write(IOW,*) ' ELEMENTS key must be terminated by an END ELEMENTS'
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        end do

                    else if (strcmp(strpar(2),'INTER',5) ) then
                        do while (.true.)
                            iblnk =1
                            read (IOR, 99001, ERR = 500, end = 500) strin
                            if (echo) write(IOW,*) strin
                            call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                            if (iblnk==1) cycle
                            if (strcmp(strpar(1), 'TYPE', 4) ) then
                                if ( nstr.ne.2 ) then
                                    write(IOW,*) '  *** ERROR DETECTED IN INPUT FILE ***'
                                    write(IOW,*)  ' TYPE key in an ELEMENT definition must specify the element type'
                                    write(IOW,*)  ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif

                                if (strcmp(strpar(2),'CPE3',4) ) then
                                    nnod = 3
                                   n_int_pts = 1
                                   nsvar = 11
                                else if (strcmp(strpar(2),'CPE4',4) ) then
                                   nnod = 4
                                   n_int_pts = 4
                                   nsvar = 11
                                else if (strcmp(strpar(2),'CPE6',4) ) then
                                   nnod = 6
                                   n_int_pts = 4
                                   nsvar = 11
                                else if (strcmp(strpar(2),'CPE8',4) ) then
                                   nnod = 8
                                   n_int_pts = 9
                                   nsvar = 11
                                else if (strcmp(strpar(2),'C3D4',4) ) then
                                   nnod = 4
                                   n_int_pts = 1
                                   nsvar = 15
                                else if (strcmp(strpar(2),'C3D10',4) ) then
                                   nnod = 10
                                   n_int_pts = 4
                                   nsvar = 15
                                else if (strcmp(strpar(2),'C3D8',4) ) then
                                   nnod = 8
                                   n_int_pts = 8
                                   nsvar = 15
                                else if (strcmp(strpar(2),'C3D20',4) ) then
                                   nnod = 20
                                   n_int_pts = 27
                                   nsvar = 15
                                else
                                   write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                   write(IOW,*) ' Unrecognized element type specified with TYPE key in ELEMENT, INTERNAL block'
                                   write(IOW,*) ' Found '
                                   write(IOW,*) strin
                                   stop
                                endif

                            elseif (strcmp(strpar(1), 'DENS',4) ) then
                                length_densities = length_densities + 1
                            else if (strcmp(strpar(1), 'PROP',4) ) then

                                if (nstr<2) then
                                    write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE ***'
                                    write(IOW,*) ' A material name must be specified after a PROPERTY key in an ELEMENT block '
                                    stop
                                endif
                                nsvar = nsvar + max_mater_states
                            else if (strcmp(strpar(1), 'INITIALS',8)) then
                                do while (.true.)
                                    iblnk =1
                                    read (IOR, 99001, ERR = 500, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDINIT', 7) ) then
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' INITIAL STATE VARIABLES key must be', &
                                                     ' terminated by an END INITIAL STATE VARIABLES'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                 end do
                            else if (strcmp(strpar(1), 'CONN',4) ) then
                                n_zones = n_zones + 1
                                do while (.true.)
                                    if (nnod==0) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' Attempt to create elements with no nodes '
                                        write(IOW,*) ' PARAMETERS must be specified for elements before CONNECTIVITY '
                                        stop
                                    endif
                                    iblnk =1
                                    read (IOR, 99001, ERR = 500, end = 500) strin
                                    if (echo) write(IOW,*) strin
                                    call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                    if (iblnk==1) cycle
                                    if ( strcmp(strpar(1), 'ENDCONN', 7) ) then
                                        exit
                                    endif
                                    if ( ityp(1)==2 ) then
                                        write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                        write(IOW,*) ' CONNECTIVITY key must be terminated by an END CONNECTIVITY'
                                        write(IOW,*) ' Found '
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    if (nstr < nnod.or.nstr > nnod+1) then
                                        write(IOW,*) ' *** Error in input file *** '
                                        write(IOW,*) ' Expecting element number (optional) and ',nnod,&
                                            ' nodes for an element but found ',nstr
                                        write(IOW,*) strin
                                        stop
                                    endif
                                    n_elements = n_elements + 1
                                    do k = 1,nnod
                                        length_connectivity = length_connectivity +1
                                    end do
                                    length_state_variables = length_state_variables + nsvar*n_int_pts
                                end do
                            else if ( strcmp(strpar(1), 'ENDELEM', 7) ) then
                                exit
                            else
                                write(IOW,*) ' *** ERROR DETECTED IN INPUT FILE *** '
                                write(IOW,*) ' ELEMENTS key must be terminated by an END ELEMENTS'
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        end do
                    endif

                else if ( strcmp(strpar(1), 'ENDMESH', 7) ) then
                    status = 0
                    if (n_zones>0) then
                        allocate(zone_list(n_zones), stat = status)
                        allocate(zone_namelist(n_zones), stat = status)
                        allocate(zone_print_flag(n_zones), stat = status)
                        allocate(zone_dimension(n_zones), stat = status)
                        allocate(zone_ndof(n_zones), stat = status)
                    endif
                    if (length_displacement_map>0) allocate(displacement_map(length_displacement_map), stat = status)
                    if (n_nodes>0) allocate(node_list(n_nodes), stat=status)
                    if (length_coords>0) allocate(coords(length_coords), stat = status)
                    if (length_dofs>0) then
                        allocate(dof_increment(length_dofs), stat = status)
                        dof_increment = 0.d0
                        allocate(dof_total(length_dofs), stat = status)
                        dof_total = 0.d0
                        allocate(rforce(length_dofs), stat = status)
                        rforce = 0.d0
                    endif
                    if (n_elements>0) then
                        allocate(element_list(n_elements), stat = status)
                        allocate(element_deleted(n_elements), stat=status)
                        allocate(energy(12*n_elements), stat=status)
                        element_deleted = .false.
                        energy = 0.d0
                    endif
                    if (length_connectivity>0) allocate(connectivity(length_connectivity), stat = status)
                    if (n_mesh_parameters>0) allocate(mesh_subroutine_parameters(n_mesh_parameters), stat=status)
                    if (length_state_variables>0) then
                        allocate(initial_state_variables(length_state_variables), stat = status)
                        allocate(updated_state_variables(length_state_variables), stat = status)
                        initial_state_variables = 0.d0
                    else
                        allocate(initial_state_variables(1), stat = status)                      ! Used as dummy variable in user subroutine
                        allocate(updated_state_variables(1), stat = status)
                    endif

                    if (length_element_properties>0) then
                        allocate(element_properties(length_element_properties), stat = status)
                    else
                        allocate(element_properties(1), stat=status)
                    endif
                    if (length_int_element_properties>0) then
                        allocate(int_element_properties(length_int_element_properties), stat = status)
                    else
                        allocate(int_element_properties(1), stat=status)
                    endif
                    if (length_densities>0) allocate(densities(length_densities), stat = status)

                    if (n_materials>0) then
                        allocate(material_list(n_materials), stat=status)
                        allocate(material_namelist(n_materials), stat=status)
                        allocate(material_properties(length_material_properties), stat=status)
                    endif

                    if (status/=0) then
                        write(IOW,*) ' Unable to allocate memory for mesh '
                        stop
                    endif

        
                    exit
                else
                    write(IOW,*) ' Error detected in input file '
                    write(IOW,*) ' Expecting an END MESH statement'
                    write(IOW,*) ' Found '
                    write(IOW,*) strin
                    stop
                endif
            end do

        else if ( strcmp(strpar(1), 'BOUNDARYCON', 11) ) then
            do while (.true.)
                iblnk =1
                read (IOR, 99001, ERR = 500, end = 500) strin
                if (echo) write(IOW,*) strin 
                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                if (iblnk==1) cycle

                if ( strcmp(strpar(1), 'HISTORY', 7) ) then
                    n_histories = n_histories + 1
                    if (nstr<2) then
                        write(IOW,*) ' *** Error in input file '
                        write(IOW,*) ' HISTORY key was used without specifying a name for the history '
                        stop
                    endif

                    do while (.true.)

                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDHIST', 7) ) then
                            exit
                        endif
                
                        if (nstr /=2) then
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' Expecting a pair of real numbers specifying time, value in a HISTORY '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        if (ityp(1) >1 .or. ityp (2) >1) then
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' HISTORY must be terminated with an END HISTORY '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        length_history_data = length_history_data + 1
                    end do

                else if ( strcmp(strpar(1), 'SUBROUTINEP', 11) ) then
                    n_subroutine_parameters = n_subroutine_parameters + 1
                    if (nstr<2) then
                        write(IOW,*) ' *** Error in input file '
                        write(IOW,*) ' SUBROUTINE PARAMETERS key was used without specifying a name for the parameters '
                        stop
                    endif
             
                    do while (.true.)

                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDSUBR', 7) ) then
                            exit
                        endif
                
                        do k = 1,nstr
                            if (ityp(k) ==2) then
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' SUBROUTINE PARAMETERS must be terminated with an ENDSUBROUTINE PARAMETERS '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            length_subroutine_parameters = length_subroutine_parameters + 1
                        end do
                    end do

                else if ( strcmp(strpar(1), 'CONSTRAINTP', 11) ) then
                    n_constraint_parameters = n_constraint_parameters + 1
                    if (nstr<2) then
                        write(IOW,*) ' *** Error in input file '
                        write(IOW,*) ' CONSTRAINT PARAMETERS key was used without specifying a name for the parameters '
                        stop
                    endif
             
                    do while (.true.)

                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDCONS', 7) ) then
                            exit
                        endif
                
                        do k = 1,nstr
                            if (ityp(k) ==2) then
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' CONSTRAINT PARAMETERS must be terminated with an ENDCONSTRAINT PARAMETERS '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            length_constraint_parameters = length_constraint_parameters + 1
                        end do
                    end do
       
                else if (strcmp(strpar(1), 'NODESET', 7) ) then
          
                    n_nodesets = n_nodesets + 1
                    if (nstr<2) then
                        write(IOW,*) ' *** Error in input file '
                        write(IOW,*) ' NODESET key was used without specifying a name for the node set '
                        stop
                    endif
                    if (nstr==3) then
                        if ( .not.strcmp(strpar(3), 'GENERATE',8) ) then
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' Expecting GENERATE option to generate a node set '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        do while (.true.)
                            iblnk =1
                            read (IOR, 99001, ERR = 500, end = 500) strin
                            if (echo) write(IOW,*) strin
                            call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                            if (iblnk==1) cycle

                            if ( strcmp(strpar(1), 'ENDNODE', 7) ) then
                                exit
                            endif
                            if (nstr/=3.or.ityp(1)/=0.or.ityp(2)/=0.or.ityp(3)/=0) then
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting start node, end node, increment for a GENERATE key '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            read(strpar(1), *) start_node
                            read(strpar(2), *) end_node
                            read(strpar(3), *) increment
                            if (increment==0) increment=1
                            do k = start_node,end_node,increment
                                length_node_lists = length_node_lists + 1
                            end do
                        end do
                    else
                        do while (.true.)
                            iblnk =1
                            read (IOR, 99001, ERR = 500, end = 500) strin
                            if (echo) write(IOW,*) strin
                            call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                            if (iblnk==1) cycle

                            if ( strcmp(strpar(1), 'ENDNODE', 7) ) then
                                exit
                            endif
                            do k = 1,nstr
                                if (ityp(k) >1) then
                                    write(IOW,*) ' *** Error in input file *** '
                                    write(IOW,*) ' Expecting a list of nodes or an END NODE SET '
                                    write(IOW,*) ' Found '
                                    write(IOW,*) strin
                                    stop
                                endif
                                length_node_lists = length_node_lists + 1
                            end do
                        end do
                    endif
                else if (strcmp(strpar(1), 'ELEMENTSET', 10) ) then
          
                    n_elementsets = n_elementsets + 1
                    if (nstr<2) then
                        write(IOW,*) ' *** Error in input file '
                        write(IOW,*) ' ELEMENTSET key was used without specifying a name for the element set '
                        stop
                    endif
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDELEM', 7) ) then
                            exit
                        endif
                        do k = 1,nstr
                            if (ityp(k) ==2) then
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting a list of nodes or an END ELEMENT SET '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            length_element_lists = length_element_lists + 1
                        end do
                    end do
                else if (strcmp(strpar(1), 'DEGREES', 7) ) then
           
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDDEGR', 7) ) then
                            exit
                        endif
                        n_prescribeddof = n_prescribeddof + 1
                                   
                        if (ityp(2) /=0.or.nstr<4.or.nstr>5) then
                            write(IOW,*) ' *** Error in input file ***'
                            write(IOW,'(A)') ' Expecting node set name, dof #, VALUE, HISTORY, or '
                            write(IOW,'(A)') '       SUBROUTINE,  history/parameter name, or value, RATE (optional) '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        if (strcmp(strpar(3), 'VALUE', 5)) then
                            length_dof_values = length_dof_values + 1
                        endif
                
                    end do
                else if (strcmp(strpar(1), 'FORCES', 6) ) then
          
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDFORC', 7) ) then
                            exit
                        endif
                        n_prescribedforces = n_prescribedforces + 1
                                   
                        if (ityp(2) /=0.or.nstr<4.or.nstr>4) then
                            write(IOW,'(A)') ' *** Error in input file ***'
                            write(IOW,'(A)') ' Expecting node set name, dof #, VALUE, HISTORY, or SUBROUTINE,'
                            write(IOW,'(A)') '  history/parameter name, or value '
                            write(IOW,'(A)') ' following a FORCE key'
                            write(IOW,'(A)') ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        if (strcmp(strpar(3), 'VALUE', 5)) then
                            length_dof_values = length_dof_values + 1
                        endif
                
                    end do

                else if (strcmp(strpar(1), 'DISTRIB', 7) ) then
           
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'ENDDISTRIB', 10) ) then
                            exit
                        endif
                        n_distributedloads = n_distributedloads + 1
                        if (ityp(2)>0) then
                            length_dload_values = length_dload_values + 1   ! ABAQUS UEL boundary condition
                        else
                            if (strcmp(strpar(3), 'VALUE', 5)) then
                                do k = 4,nstr
                                    length_dload_values = length_dload_values + 1
                                end do
                            else if (strcmp(strpar(3), 'HISTORY', 7)) then
                                do k = 5,nstr
                                    length_dload_values = length_dload_values + 1
                                end do
                            else if (strcmp(strpar(3), 'NORMAL', 5)) then
                                do k = 5,nstr
                                    length_dload_values = length_dload_values + 1
                                end do
                            endif
                        endif
                    end do

                else if (strcmp(strpar(1), 'CONSTRA', 7) ) then
 
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if ( strcmp(strpar(1), 'CONSTRAINNODEPA',15) ) then
                            if (nstr<6.or.(ityp(2)>0).or.ityp(4)>0.or.ityp(6)>0.or.nstr>7) then
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting # of node pairs, nodeset or node, dof,'
                                write(IOW,*) ' nodeset or node, dof, parameter name (optional) '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                            if (ityp(3)==0.and.ityp(5)==0) then
                                n_constraints = n_constraints + 1
                            else if (ityp(3)==2.and.ityp(5)==2) then
                                read(strpar(2),*) k
                                n_constraints = n_constraints + k
                            else
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting # of node pairs, nodeset or node, dof,'
                                write(IOW,*) ' nodeset or node, dof, parameter name (optional) for CONSTRAIN NODE PAIRS constraint'
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        else if (strcmp(strpar(1), 'TIENO',5) ) then
                            if (nstr>5.and.ityp(2)==0.and.ityp(3)==2.and.ityp(4)==0.and.ityp(5)==0.and.ityp(6)==0) then
                                read(strpar(2),*) k
                                n_constraints = n_constraints + k
                            else
                                write(IOW,*) ' *** Error in input file *** '
                                write(IOW,*) ' Expecting nodeset or node, dof, node, dof,'
                                write(IOW,*) ' parameter name (optional)for TIE NODE PAIRS constraint '
                                write(IOW,*) ' Found '
                                write(IOW,*) strin
                                stop
                            endif
                        else if (strcmp(strpar(1), 'CONSTRAINNODESET',16) ) then
                            n_constraints = n_constraints + 1
                        else if ( strcmp(strpar(1), 'ENDCONS', 7) ) then
                            exit
                        else
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' CONSTRAINTS must be terminated by an END CONSTRAINTS '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
               
                    end do

                else if ( strcmp(strpar(1), 'ENDBOUN', 7) ) then

                    if (n_histories>0) allocate(history_data(2,length_history_data), stat = status)
                    if (length_subroutine_parameters>0) then
                        allocate(subroutine_parameters(length_subroutine_parameters), stat = status)
                    else
                        allocate(subroutine_parameters(1), stat = status)
                    endif
                    if (length_constraint_parameters>0) then
                        allocate(constraint_parameters(length_constraint_parameters), stat = status)
                    else
                        allocate(constraint_parameters(1), stat = status)
                    endif
                    if (length_dof_values>0) allocate(dof_values(length_dof_values), stat = status)
                    if (length_dload_values>0) allocate(dload_values(length_dload_values), stat = status)
 
                    if (n_histories>0) then
                        allocate(history_list(n_histories), stat = status)
                        allocate(history_namelist(n_histories), stat = status)
                    endif
                    if (n_subroutine_parameters>0) then
                        allocate(subroutineparameter_list(n_subroutine_parameters), stat = status)
                        allocate(subroutineparameter_namelist(n_subroutine_parameters), stat = status)
                    endif
                    if (n_constraint_parameters>0) then
                        allocate(constraintparameter_list(n_constraint_parameters), stat = status)
                        allocate(constraintparameter_namelist(n_constraint_parameters), stat = status)
                    endif
                    if (n_nodesets>0) then
                        allocate(nodeset_list(n_nodesets), stat = status)
                        allocate(nodeset_namelist(n_nodesets), stat = status)
                    endif
                    if (length_node_lists>0) allocate(node_lists(length_node_lists), stat = status)
                    if (n_elementsets>0) then
                        allocate(elementset_list(n_elementsets), stat = status)
                        allocate(elementset_namelist(n_elementsets), stat = status)
                    endif
                    if (length_element_lists>0) allocate(element_lists(length_element_lists), stat = status)
                    if (n_prescribeddof>0) allocate(prescribeddof_list(n_prescribeddof), stat = status)
                    if (n_prescribedforces>0) allocate(prescribedforce_list(n_prescribedforces), stat = status)
                    if (n_distributedloads>0)  allocate(distributedload_list(n_distributedloads), stat = status)
                    if (n_constraints>0) then
                        allocate(constraint_list(n_constraints), stat = status)
                        allocate(lagrange_multipliers(n_constraints), stat = status)
                    endif
                    exit
                else
                    write(IOW,*) ' *** Error in input file *** '
                    write(IOW,*) ' BOUNDARY CONDITIONS key must be terminated with an END BOUNDARY CONDITIONS'
                    write(IOW,*) ' Found '
                    write(IOW,*) strin
                    stop
                endif
            end do
        else if (strcmp(strpar(1), 'PRINTINI', 8) ) then
            n_total_files = n_total_files + 1
        else if ( strcmp(strpar(1), 'STATICSTEP', 10) ) then
            do while (.true.)
                iblnk =1
                read (IOR, 99001, ERR = 500, end = 500) strin
                if (echo) write(IOW,*) strin
                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                if (iblnk==1) cycle

                if (strcmp(strpar(1),'USERPRINTF',10) ) then
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDU',4) ) then
                            exit
                        endif
                        if (nstr>1.or.ityp(1)<2) then
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' USER PRINT FILES must be terminated by an END USER PRINT FILES'
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        n_user_print_files = n_user_print_files + 1
                        n_total_files = n_total_files + 1
                    end do
                else if (strcmp(strpar(1),'USERPRINTP',10) ) then
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDU',4) ) then
                            exit
                        endif
                        if (ityp(1)==2) then
                            write(IOW,*) ' *** Error in input file ***'
                            write(IOW,*) ' USER PRINT PARAMETERS must be terminated by '
                            write(IOW,*) ' an END USER PRINT PARAMETERS '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        do k = 1,nstr
                            if (ityp(k)==2) then
                                write(IOW,*) ' Expecting an integer or real valued USER PRINT PARAMETER'
                                write(IOW,*) ' Found'
                                write(IOW,*) strin
                                stop
                            endif
                            n_user_print_parameters = n_user_print_parameters+1
                        end do
                    end do
                else if (strcmp(strpar(1),'PRINTSTATE',10) ) then
                    n_total_files = n_total_files + 1
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if (strcmp(strpar(1),'ENDPR',5) ) then
                            exit
                        else if (strcmp(strpar(1),'FIELD',5) ) then
                            n_field_variables = nstr-1
                        else if (strcmp(strpar(1),'DEGR',4) ) then
                            cycle
                        else if (strcmp(strpar(1),'ZON',3) ) then
                            do while (.true.)
                                iblnk =1
                                read (IOR, 99001, ERR = 500, end = 500) strin
                                if (echo) write(IOW,*) strin
                                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                                if (iblnk==1) cycle
                                if (strcmp(strpar(1),'ENDZ',4) ) exit
                            end do
                        else if (strcmp(strpar(1),'DISPL',5) ) then
                            cycle
                        else
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' Unexpected PRINT STATE option '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            write(IOW,*) ' An END PRINT STATE or END ZONES may be missing'
                            stop
                        endif

                    end do
                else if (strcmp(strpar(1),'ENDSTATIC', 9) ) then
                    if (n_user_print_parameters>0) allocate(user_print_parameters(n_user_print_parameters), stat=status)
                    if (n_user_print_files>0) then
                        allocate(user_print_filenames(n_user_print_files), stat=status)
                        allocate(user_print_units(n_user_print_files), stat=status)
                    endif
                    if (n_field_variables>0) allocate(field_variable_names(n_field_variables), stat=status)
                    if (status/=0) then
                        write(IOW,*) ' Unable to allocate memory for static step '
                        stop
                    endif
                    exit
                else
                    cycle
                end if
        
            end do
        else if ( strcmp(strpar(1), 'EXPLICITDYNAMICSTEP', 19) ) then
            do while (.true.)
                iblnk =1
                read (IOR, 99001, ERR = 500, end = 500) strin
                if (echo) write(IOW,*) strin
                call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                if (iblnk==1) cycle

                if (strcmp(strpar(1),'USERPRINTF',10) ) then
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDU',4) ) then
                            exit
                        endif
                        if (nstr>1.or.ityp(1)<2) then
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' USER PRINT FILES must be terminated by an END USER PRINT FILES'
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        n_user_print_files = n_user_print_files + 1
                        n_total_files = n_total_files + 1
                    end do
                else if (strcmp(strpar(1),'USERPRINTP',10) ) then
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle
                        if (strcmp(strpar(1),'ENDU',4) ) then
                            exit
                        endif
                        if (ityp(1)==2) then
                            write(IOW,*) ' *** Error in input file ***'
                            write(IOW,*) ' USER PRINT PARAMETERS must be terminated by '
                            write(IOW,*) ' an END USER PRINT PARAMETERS '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            stop
                        endif
                        do k = 1,nstr
                            if (ityp(k)==2) then
                                write(IOW,*) ' Expecting an integer or real valued USER PRINT PARAMETER'
                                write(IOW,*) ' Found'
                                write(IOW,*) strin
                                stop
                            endif
                            n_user_print_parameters = n_user_print_parameters+1
                        end do
                    end do
                else if (strcmp(strpar(1),'PRINTSTATE',10) ) then
                    n_total_files = n_total_files + 1
                    do while (.true.)
                        iblnk =1
                        read (IOR, 99001, ERR = 500, end = 500) strin
                        if (echo) write(IOW,*) strin
                        call parse(strin, strpar, 100, nstr, lenstr, ityp, iblnk)
                        if (iblnk==1) cycle

                        if (strcmp(strpar(1),'ENDPR',5) ) then
                            exit
                        else if (strcmp(strpar(1),'FIELD',5) ) then
                            n_field_variables = nstr-1
                        else if (strcmp(strpar(1),'DEGR',4) ) then
                            cycle
                        else if (strcmp(strpar(1),'ZON',3) ) then
                            cycle
                        else if (strcmp(strpar(1),'DISPL',5) ) then
                            cycle
                        else if (strcmp(strpar(1),'ENDZ',4) ) then
                            cycle
                        else
                            write(IOW,*) ' *** Error in input file *** '
                            write(IOW,*) ' Unexpected PRINT STATE option '
                            write(IOW,*) ' Found '
                            write(IOW,*) strin
                            write(IOW,*) ' An END PRINT STATE or END ZONES may be missing'
                            stop
                        endif
                    end do

                else if (strcmp(strpar(1),'ENDEXP', 6) ) then
                    if (n_user_print_parameters>0) allocate(user_print_parameters(n_user_print_parameters), stat=status)
                    if (n_user_print_files>0) then
                        allocate(user_print_filenames(n_user_print_files), stat=status)
                        allocate(user_print_units(n_user_print_files), stat=status)
                    endif
                    if (n_field_variables>0) allocate(field_variable_names(n_field_variables), stat=status)

                    if (status/=0) then
                        write(IOW,*) ' Unable to allocate memory for explicit dynamic step '
                        stop
                    endif


                    exit
                else
                    cycle
                end if
        
            end do
            
        else
            cycle
        endif
    end do
  
500 continue                                          ! End of file
  
    allocate(file_units(n_total_files), stat = status)
    allocate(file_list(n_total_files), stat = status)
  
  
    if (status /=0) then
        write(IOW,*) ' Error in subroutine allocate_mesh_storage: unable to allocate memory to file names '
        stop
    endif
   
    return
99001 format (A100)
  
end subroutine allocate_mesh_storage
