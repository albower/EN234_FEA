module Boundaryconditions
    use Types
    !  Data defining a time history of prescribed DOF, generalized forces, or generalized traction
    type history
        sequence
        integer :: index                         ! Index of start of time, value pairs in history_data(:)
        integer :: n_timevals                    ! No. time/value pairs
    end type history
    !  Parameters passed to user-subroutine controlled boundary conditions
    type subroutineparameters
        sequence
        integer :: index                         ! Index of start of parameter list in subroutine_parameters(:)
        integer :: n_parameters                  ! No. subroutine parameters
    end type subroutineparameters
    !  Parameters for constraints
    type constraintparameters
        sequence
        integer :: index                        ! Index of start of parameter list in constraint_parameters(:)
        integer :: n_parameters                 ! No. constraint parameters
    end type constraintparameters
    !  Parameters for node sets
    type nodeset
        sequence
        integer :: index                        ! Index of start of list in node_lists(:)
        integer :: n_nodes                      ! No. nodes in set
    end type nodeset
    !  Parameters for element sets
    type elementset
        sequence
        integer :: index                        ! Index of start of list in element_lists(:)
        integer :: n_elements                   ! No. elements in set
    end type elementset

    type prescribeddof
        sequence
        integer :: flag                          ! Flag specifying nature of DOF - flag=1 prescribed value; flag=2 history, flag=3 user subroutine
        integer :: dof                           ! DOF to prescribe
        integer :: node_number                   ! Node number if single node is constrained, zero otherwise
        integer :: node_set                      ! Node set number if a set of nodes is constrained, zero otherwise
        integer :: history_number                ! Index of DOF history table
        integer :: subroutine_parameter_number   ! Index of user subroutine parameters
        integer :: index_dof_values              ! Index for value of prescribed DOF
        integer :: rate_flag                     ! Set to 1 if DOF rate is to be prescribed; 0 if value is to be prescribed.
    end type prescribeddof

    type prescribedforce
        sequence
        integer :: flag                          ! Flag specifying nature of force - flag=1 prescribed value; flag=2 history, flag=3 user subroutine
        integer :: dof                           ! DOF to apply force to
        integer :: node_number                   ! Node number if single node is loaded, zero otherwise
        integer :: node_set                      ! Node set number if a set of nodes is loaded, zero otherwise
        integer :: history_number                ! Index of DOF history table
        integer :: subroutine_parameter_number   ! Index of user subroutine parameters
        integer :: index_dof_values              ! Index for value of prescribed force
    end type prescribedforce

    type distributedloads
        sequence
        integer :: flag                          ! Flag specifying distributed load type. flag=1, traction value given, flag=2 history and direction give, flag=3 history and normal to element, flag=4 user subroutine
        integer :: element_set                   ! Element set to be loaded
        integer :: face                          ! Face to be loaded
        integer :: history_number                ! History number specifying variation of traction
        integer :: subroutine_parameter_number   !
        integer :: index_dload_values
        integer :: n_dload_values
    end type distributedloads
   
    type constraint
        sequence
        integer :: flag                       ! Flag identifying type of constraint.  flag=1 tie a node pair; flag=2, tie a node set, flag=3 general MPC
        integer :: node1                      ! First node (or node set) in constraint
        integer :: dof1                       ! Dof for first node (not used if flag=3)
        integer :: node2                      ! Second node (or node set, or DOF list) in constraint
        integer :: dof2                       ! Dof for second node (not used if flag=3)
        integer :: index_parameters           ! Index of parameters associated with constraint in parameter list
    end type constraint


   
    integer, save :: n_histories                   ! No. load histories
    integer, save :: n_subroutine_parameters       ! No. lists of subroutine parameters for distributed loads or DOF
    integer, save :: n_constraint_parameters       ! No. lists of parameters for constraints
    integer, save :: n_nodesets                    ! No. node sets
    integer, save :: n_elementsets                 ! No. element sets
    integer, save :: n_prescribeddof               ! No. prescribed DOF
    integer, save :: n_prescribedforces            ! No. prescribed forces
    integer, save :: n_distributedloads            ! No. distributed loads or fluxes
    integer, save :: n_constraints                 ! No. constraints
   
    integer, save :: length_node_lists
    integer, save :: length_element_lists
    integer, save :: length_history_data
    integer, save :: length_subroutine_parameters
    integer, save :: length_dof_values
    integer, save :: length_dload_values
    integer, save :: length_constraint_parameters
   
    integer, save, allocatable :: node_lists(:)
    integer, save, allocatable :: element_lists(:)
   
    real (prec), save, allocatable :: history_data(:,:)
    real (prec), save, allocatable :: subroutine_parameters(:)
    real (prec), save, allocatable :: dload_values(:)
    real (prec), save, allocatable :: dof_values(:)
    real (prec), save, allocatable :: constraint_parameters(:)
   
    real (prec), save, allocatable :: lagrange_multipliers(:)
   
    character (len=100), allocatable :: elementset_namelist(:)
    character (len=100), allocatable :: nodeset_namelist(:)
    character (len=100), allocatable :: history_namelist(:)
    character (len=100), allocatable :: subroutineparameter_namelist(:)
    character (len=100), allocatable :: constraintparameter_namelist(:)
     
    type (history), save, allocatable :: history_list(:)
    type (subroutineparameters), save, allocatable :: subroutineparameter_list(:)
    type (constraintparameters), save, allocatable :: constraintparameter_list(:)
    type (elementset), save, allocatable :: elementset_list(:)
    type (nodeset), save, allocatable :: nodeset_list(:)
    type (prescribeddof), save, allocatable :: prescribeddof_list(:)
    type (prescribedforce), save, allocatable :: prescribedforce_list(:)
    type (distributedloads), save, allocatable :: distributedload_list(:)
    type (constraint), save, allocatable :: constraint_list(:)
   
contains
    function find_name_index(name,len,name_list,n_names)
        use Types
        use ParamIO
        implicit none
        integer                     :: len
        integer                     :: n_names
        integer                     :: find_name_index
        character ( len = 100 )     :: name
        character ( len = 100 )     :: name_list(n_names)

        ! Local Variables
        integer :: n
        logical :: strcmp

        !     Find the pointers for the zone named zname
        !     by searching through the list of named zones

        do n = 1, n_names
            if (len == len_trim(name_list(n)) ) then
                if ( strcmp(name, name_list(n), len) ) then
                    find_name_index = n
                    return
                end if
            endif
        end do

        write(IOW,*) ' *** Error detected in input file ***'
        write(IOW,*) ' The named history, node set, element set, or parameter list ',name
        write(IOW,*) ' was not defined '
        stop

    end function find_name_index
 
 

    subroutine interpolate_history_table(history,nhist,time,dofvalue)

        use Types
   

        integer,  intent( in )         :: nhist

        real( prec ), intent( in )     :: history(2,nhist)
        real( prec ), intent( in )     :: time
        real( prec ), intent( out )    :: dofvalue


        integer :: klo,khi,k

        !   Subroutine to interpolate a load history table
        !
        !  Find positions in history table within which to interpolate

        if (nhist == 1) then
            dofvalue = history(2,1)
            return
        endif

        if (time <= history(1,1) ) then
            dofvalue = history(2,1)
        else if (time >= history(1,nhist) ) then
            dofvalue = history(2,nhist)
        else

            klo = 1
            khi = nhist

            do while ( .true. )

                if ( khi - klo>1 ) then
                    k = (khi + klo)/2
                    if ( history(1, k)>time ) then
                        khi = k
                    else
                        klo = k
                    end if
                    cycle
                else
                    exit
                end if

            end do

            if (history(1,khi) == history(1,klo)) then
                dofvalue = 0.5D0*(history(2,khi) + history(2,klo))
                return
            endif
            dofvalue =( (history(1,khi)-time)*history(2,klo) +  &
                (time-history(1,klo))*history(2,khi)  )/ &
                (history(1,khi)-history(1,klo))

        endif

        return

    end subroutine interpolate_history_table
   
end module
