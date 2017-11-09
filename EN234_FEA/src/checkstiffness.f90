subroutine check_stiffness(element_flag)
    use Types
    use ParamIO
    use Globals
    use User_Subroutine_Storage
    use Mesh
    use Boundaryconditions
    use Stiffness
    use Staticstepparameters, only : nonlinear,current_step_number,max_total_time
    use Controlparameters, only : abaqusformat
    implicit none

    integer, intent(in)   :: element_flag

    ! Local Variables
    integer         :: lmn,n,nn,m,i,j,k,iof
    integer         :: ix,iu,ns
    integer         :: icount,icount2
    integer         :: status
    integer      :: mat_prop_index,n_mat_props
    integer      :: abq_JTYPE
    integer      :: abq_MDLOAD
    integer      :: abq_NPREDF
    integer      :: abq_LFLAGS(5)

    real ( prec ) :: err
    real ( prec ) :: element_PNEWDT
    logical :: fail

    real( prec ), allocatable    :: element_coords(:)
    real( prec ), allocatable    :: element_dof_increment(:)
    real( prec ), allocatable    :: element_dof_total(:)
    real( prec ), allocatable    :: element_state_variables(:)
                                                            
    real( prec ), allocatable   :: element_stiffness(:,:)
    real( prec ), allocatable   :: stif1(:,:)
    real( prec ), allocatable   :: numerical_stiffness(:,:)
    real( prec ), allocatable   :: element_residual(:)
    real( prec ), allocatable   :: resid0(:)
    real( prec ), allocatable   :: resid1(:)

    real( prec ), allocatable    :: abq_PREDEF(:,:,:)
    real( prec ), allocatable    :: abq_ADLMAG(:)
    real( prec ), allocatable    :: abq_DDLMAG(:)
    real( prec ), allocatable    :: abq_V(:)
    real( prec ), allocatable    :: abq_A(:)

    type (node), allocatable ::  local_nodes(:)

    character (len=80) material_name

    integer, allocatable :: abq_JDLTYP(:)

    real( prec ) :: abq_time(2)      ! Time variable for ABAQUS UEL interface
    real( prec ) :: abq_PARAMS(3)    ! Parameter variable for ABAQUS UEL interface
    real( prec ) :: abq_PNEWDT

    !     Subroutine to compare stiffness matrix calculated in a user subroutine with its numerical derivative
    !     It is activated by the CHECK STIFFNESS keyword

    if (abaqusformat) then
       call generate_abaqus_dloads
    endif

    j = length_coord_array
    if (abaqusformat) j = max(length_coord_array,length_dof_array)

    allocate(element_coords(j), stat = status)
    allocate(element_dof_increment(length_dof_array), stat = status)
    allocate(element_dof_total(length_dof_array), stat = status)
    allocate(local_nodes(length_node_array), stat = status)
    allocate(element_state_variables(length_state_variable_array), stat=status)

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

    if (abaqusformat) then
       allocate(abq_PREDEF(2,1,length_node_array), stat = status)
       allocate(abq_V(length_dof_array), stat=status)
       allocate(abq_A(length_dof_array), stat=status)
    endif

    if (length_abq_dlmag_array>0) then
       allocate(abq_ADLMAG(length_abq_dlmag_array), stat = status)
       allocate(abq_DDLMAG(length_abq_dlmag_array), stat = status)
       allocate(abq_JDLTYP(length_abq_dlmag_array), stat = status)
    else
       allocate(abq_ADLMAG(1), stat = status)
       allocate(abq_DDLMAG(1), stat = status)
       allocate(abq_JDLTYP(1), stat = status)
    endif
    if (status/=0) then
        write(IOW,*) ' Error in subroutine assemble_direct_stiffness '
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

    write(IOW,*) ' '
    write(IOW,*) ' '
    write(IOW,*) ' ===================== STIFFNESS CHECK ============================'
    if (element_flag==10002) then
        write(IOW,*)
        write(IOW,*) ' Testing 2D continuum element with UMAT '
        write(IOW,*) ' WARNING: the ABAQUS format continuum elements use an approximate stiffness matrix '
        write(IOW,*) ' Consider checking the material tangent in the UMAT with a CHECK MATERIAL TANGENT key instead'
    else if (element_flag==10003) then
        write(IOW,*)
        write(IOW,*) ' Testing 3D continuum element with UMAT '
        write(IOW,*) ' WARNING: the ABAQUS format continuum elements use an approximate stiffness matrix '
        write(IOW,*) ' Consider checking the material tangent in the UMAT with a CHECK MATERIAL TANGENT key instead '
    else if (element_flag<99999) then
        write(IOW,*)
        write(IOW,*) ' Testing user defined element with flag ',element_flag
    else
        write(IOW,*)
        write(IOW,'(A31,I6)') ' Testing ABAQUS UEL with flag U',element_flag-99999
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
        if (ns==0) ns=1

    if (element_list(lmn)%flag==10002) then            ! 3D continuum element

       element_PNEWDT = 1.d99

        mat_prop_index = material_list(element_list(lmn)%material_index)%prop_index
        n_mat_props = material_list(element_list(lmn)%material_index)%n_properties
        material_name = material_namelist(element_list(lmn)%material_index)(1:80)

        call continuum_element_static_2D(current_step_number,lmn, element_list(lmn)%flag, &
            element_list(lmn)%n_nodes, local_nodes(1:element_list(lmn)%n_nodes), &               ! Input variables
            n_mat_props, material_properties(mat_prop_index),material_name,  &               ! Input variables
            element_coords(1:ix),ix, &                                                       ! Input variables
            element_dof_increment(1:iu), element_dof_total(1:iu),iu,  &                                              ! Input variables
            ns, initial_state_variables(iof:iof+ns-1),&                                               ! Input variables
            updated_state_variables(iof:iof+ns-1),resid0(1:iu),element_stiffness(1:iu,1:iu),element_PNEWDT)


        icount = 0
        do j = 1, element_list(lmn)%n_nodes
            n = connectivity(element_list(lmn)%connect_index + j - 1)
            do i = 1,node_list(n)%n_dof
                icount = icount + 1
                element_dof_increment(icount) = element_dof_increment(icount) + 1.D-07

                call continuum_element_static_2D(current_step_number,lmn, element_list(lmn)%flag, &
                    element_list(lmn)%n_nodes, local_nodes(1:element_list(lmn)%n_nodes), &               ! Input variables
                    n_mat_props, material_properties(mat_prop_index),material_name,  &               ! Input variables
                    element_coords(1:ix),ix, &                                                       ! Input variables
                    element_dof_increment(1:iu), element_dof_total(1:iu),iu,  &                                              ! Input variables
                    ns, initial_state_variables(iof:iof+ns-1),&                                               ! Input variables
                    updated_state_variables(iof:iof+ns-1),resid1(1:iu),stif1(1:iu,1:iu),element_PNEWDT)

                numerical_stiffness(1:iu,icount) = -(resid1(1:iu)-resid0(1:iu))/1.D-07

                write(IOW,*)
                write(IOW,*) ' Column ',icount, ' node ',n,' DOF ',i

                icount2 = 0
                do m = 1,element_list(lmn)%n_nodes
                    nn = connectivity(element_list(lmn)%connect_index + m - 1)
                    do k = 1,node_list(nn)%n_dof
                        icount2 = icount2 + 1
                        write(IOW,4000) icount2,j,k,element_stiffness(icount2,icount),numerical_stiffness(icount2,icount)
4000                    format( ' Row ',i4,' node ',i4,' DOF ',i4,' Stiffness ',d15.5,' Numerical deriv ',d15.5 )
                    end do
                end do
                element_dof_increment(icount) = element_dof_increment(icount) - 1.D-07
            end do
        end do

        err = sum( (numerical_stiffness-element_stiffness)*(numerical_stiffness-element_stiffness) ) / &
              sum(element_stiffness*element_stiffness)

        if (err>1.d-06) then
            write(IOW,*)
            write(IOW,*) ' The normalized difference between the numerical and exact stiffness matrix is ',err
            write(IOW,*)
            write(IOW,*) ' Note that ABAQUS format continuum elements use an approximate stiffness '
            write(IOW,*) ' For ABAQUS elements with a UMAT, consider using the '
            write(IOW,*) ' CHECK MATERIAL TANGENT key instead of CHECK STIFFNESS '
            write(IOW,*) ''
            write(IOW,*) ' If you are testing a user element, check the printed exact and numerical stiffness '
            write(IOW,*) ' for errors.   Small differences may be caused by rounding error.   Large ones '
            write(IOW,*) ' (more than a few %) are likely to be caused by an error in your coded stiffness '
            write(IOW,*) ' '
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



    else if (element_list(lmn)%flag==10003) then            ! 3D continuum element


        element_PNEWDT = 1.d99

        mat_prop_index = material_list(element_list(lmn)%material_index)%prop_index
        n_mat_props = material_list(element_list(lmn)%material_index)%n_properties
        material_name = material_namelist(element_list(lmn)%material_index)(1:80)

        call continuum_element_static_3D(current_step_number,lmn, element_list(lmn)%flag, &
            element_list(lmn)%n_nodes, local_nodes(1:element_list(lmn)%n_nodes), &               ! Input variables
            n_mat_props, material_properties(mat_prop_index),material_name,  &               ! Input variables
            element_coords(1:ix),ix, &                                                       ! Input variables
            element_dof_increment(1:iu), element_dof_total(1:iu),iu,  &                                              ! Input variables
            ns, initial_state_variables(iof:iof+ns-1),&                                               ! Input variables
            updated_state_variables(iof:iof+ns-1),resid0(1:iu),element_stiffness(1:iu,1:iu),element_PNEWDT)

        icount = 0
        do j = 1, element_list(lmn)%n_nodes
            n = connectivity(element_list(lmn)%connect_index + j - 1)
            do i = 1,node_list(n)%n_dof
                icount = icount + 1
                element_dof_increment(icount) = element_dof_increment(icount) + 1.D-07

                call continuum_element_static_3D(current_step_number,lmn, element_list(lmn)%flag, &
                    element_list(lmn)%n_nodes, local_nodes(1:element_list(lmn)%n_nodes), &               ! Input variables
                    n_mat_props, material_properties(mat_prop_index),material_name,  &               ! Input variables
                    element_coords(1:ix),ix, &                                                       ! Input variables
                    element_dof_increment(1:iu), element_dof_total(1:iu),iu,  &                                              ! Input variables
                    ns, initial_state_variables(iof:iof+ns-1),&                                               ! Input variables
                    updated_state_variables(iof:iof+ns-1),resid1(1:iu),stif1(1:iu,1:iu),element_PNEWDT)

                numerical_stiffness(1:iu,icount) = -(resid1(1:iu)-resid0(1:iu))/1.D-07

                write(IOW,*)
                write(IOW,*) ' Column ',icount, ' node ',n,' DOF ',i

                icount2 = 0
                do m = 1,element_list(lmn)%n_nodes
                    nn = connectivity(element_list(lmn)%connect_index + m - 1)
                    do k = 1,node_list(nn)%n_dof
                        icount2 = icount2 + 1
                        write(IOW,3000) icount2,j,k,element_stiffness(icount2,icount),numerical_stiffness(icount2,icount)
3000                    format( ' Row ',i4,' node ',i4,' DOF ',i4,' Stiffness ',d15.5,' Numerical deriv ',d15.5 )
                    end do
                end do
                element_dof_increment(icount) = element_dof_increment(icount) - 1.D-07
            end do
        end do

        err = sum( (numerical_stiffness-element_stiffness)*(numerical_stiffness-element_stiffness) ) / &
              sum(element_stiffness*element_stiffness)

        if (err>1.d-06) then
            write(IOW,*)
            write(IOW,*) ' The normalized difference between the numerical and exact stiffness matrix is ',err
            write(IOW,*)
            write(IOW,*) ' Note that ABAQUS format continuum elements use an approximate stiffness '
            write(IOW,*) ' For ABAQUS elements with a UMAT, consider using the '
            write(IOW,*) ' CHECK MATERIAL TANGENT key instead of CHECK STIFFNESS '
            write(IOW,*) ''
            write(IOW,*) ' If you are testing a user element, check the printed exact and numerical stiffness '
            write(IOW,*) ' for errors.   Small differences may be caused by rounding error.   Large ones '
            write(IOW,*) ' (more than a few %) are likely to be caused by an error in your coded stiffness '
            write(IOW,*) ' '
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



    else if (element_list(lmn)%flag>99999) then                ! ABAQUS UEL format user subroutine
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

        !     Form element stiffness
        iof = element_list(lmn)%state_index
        ns = element_list(lmn)%n_states
        if (ns==0) then 
           ns=1
           element_state_variables(1) = 0.d0
        else
          element_state_variables(1:ns) = initial_state_variables(iof:iof+ns-1)
        endif
        
        call UEL(resid0(1:iu),element_stiffness(1:iu,1:iu),element_state_variables, &
            energy(8*lmn-7:8*lmn),iu,1,ns, &
            element_properties(element_list(lmn)%element_property_index),element_list(lmn)%n_element_properties, &
            element_coords(1:ix),abq_MCRD(lmn),element_list(lmn)%n_nodes,element_dof_increment(1:iu)+element_dof_total(1:iu), &
            element_dof_increment(1:iu),abq_V(1:iu),abq_A(1:iu),abq_JTYPE,abq_time,DTIME, &
            1,current_step_number,lmn,abq_PARAMS,abq_MDLOAD,abq_JDLTYP,abq_ADLMAG, &
            abq_PREDEF(1,1,1:element_list(lmn)%n_nodes),abq_NPREDF, &
            abq_LFLAGS,iu,abq_DDLMAG,abq_MDLOAD,abq_PNEWDT, &
            int_element_properties(element_list(lmn)%int_element_property_index), &
            element_list(lmn)%n_int_element_properties, &
            max_total_time)


        icount = 0
        do j = 1, element_list(lmn)%n_nodes
            n = connectivity(element_list(lmn)%connect_index + j - 1)
            do i = 1,node_list(n)%n_dof
                icount = icount + 1
                element_dof_increment(icount) = element_dof_increment(icount) + 1.D-07

                call UEL(resid1(1:iu),stif1(1:iu,1:iu),element_state_variables, &
                    energy(8*lmn-7:8*lmn),iu,1,ns, &
                    element_properties(element_list(lmn)%element_property_index),element_list(lmn)%n_element_properties, &
                    element_coords(1:ix),abq_MCRD(lmn),element_list(lmn)%n_nodes, &
                    element_dof_increment(1:iu)+element_dof_total(1:iu), &
                    element_dof_increment(1:iu),abq_V(1:iu),abq_A(1:iu),abq_JTYPE,abq_time,DTIME, &
                    1,current_step_number,lmn,abq_PARAMS,abq_MDLOAD,abq_JDLTYP,abq_ADLMAG, &
                    abq_PREDEF(1,1,1:element_list(lmn)%n_nodes),abq_NPREDF, &
                    abq_LFLAGS,iu,abq_DDLMAG,abq_MDLOAD,abq_PNEWDT, &
                    int_element_properties(element_list(lmn)%int_element_property_index), &
                    element_list(lmn)%n_int_element_properties, &
                    max_total_time)


                numerical_stiffness(1:iu,icount) = -(resid1(1:iu)-resid0(1:iu))/1.D-07

                write(IOW,*)
                write(IOW,*) ' Column ',icount, ' node ',n,' DOF ',i

                icount2 = 0
                do m = 1,element_list(lmn)%n_nodes
                    nn = connectivity(element_list(lmn)%connect_index + m - 1)
                    do k = 1,node_list(nn)%n_dof
                        icount2 = icount2 + 1
                        write(IOW,1000) icount2,j,k,element_stiffness(icount2,icount),numerical_stiffness(icount2,icount)
1000                    format( ' Row ',i4,' node ',i4,' DOF ',i4,' Stiffness ',d15.5,' Numerical deriv ',d15.5 )
                    end do
                end do
                element_dof_increment(icount) = element_dof_increment(icount) - 1.D-07
            end do
        end do
    !
    else

        iof = element_list(lmn)%state_index
        if (iof==0) iof = 1
        ns = element_list(lmn)%n_states
        call user_element_static(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes, &
            local_nodes(1:element_list(lmn)%n_nodes), &       ! Input variables
            element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index:),  &     ! Input variables
            element_list(lmn)%n_int_element_properties, int_element_properties(element_list(lmn)%int_element_property_index), &
            element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu), iu,      &                              ! Input variables
            ns, initial_state_variables(iof:iof+ns), &                           ! Input variables
            updated_state_variables(iof:iof+ns),element_stiffness(1:iu,1:iu),resid0(1:iu), fail)               ! Output variables

        !     Compute numerical derivative of stiffness

        icount = 0
        do j = 1, element_list(lmn)%n_nodes
            n = connectivity(element_list(lmn)%connect_index + j - 1)
            do i = 1,node_list(n)%n_dof
                icount = icount + 1
                element_dof_increment(icount) = element_dof_increment(icount) + 1.D-07

                call user_element_static(lmn, element_list(lmn)%flag, element_list(lmn)%n_nodes, &
                    local_nodes(1:element_list(lmn)%n_nodes), &       ! Input variables
                    element_list(lmn)%n_element_properties, element_properties(element_list(lmn)%element_property_index),  &     ! Input variables
                    element_list(lmn)%n_int_element_properties, &
                    int_element_properties(element_list(lmn)%int_element_property_index), &
                    element_coords(1:ix),ix, element_dof_increment(1:iu), element_dof_total(1:iu), iu,      &                              ! Input variables
                    ns, initial_state_variables(iof:iof+ns), &                           ! Input variables
                    updated_state_variables(iof:iof+ns),stif1(1:iu,1:iu),resid1(1:iu), fail)               ! Output variables
          
                numerical_stiffness(1:iu,icount) = -(resid1(1:iu)-resid0(1:iu))/1.D-07

                write(IOW,*)
                write(IOW,*) ' Column ',icount, ' node ',n,' DOF ',i
         
                icount2 = 0
                do m = 1,element_list(lmn)%n_nodes
                    nn = connectivity(element_list(lmn)%connect_index + m - 1)
                    do k = 1,node_list(nn)%n_dof
                        icount2 = icount2 + 1
                        write(IOW,2000) icount2,j,k,element_stiffness(icount2,icount),numerical_stiffness(icount2,icount)
2000                    format( ' Row ',i4,' node ',i4,' DOF ',i4,' Stiffness ',d15.5,' Numerical deriv ',d15.5 )
                    end do
                end do
                element_dof_increment(icount) = element_dof_increment(icount) - 1.D-07
            end do
        end do

        err = sum( (numerical_stiffness-element_stiffness)*(numerical_stiffness-element_stiffness) ) / &
              sum(element_stiffness*element_stiffness)
     
        if (err>1.d-06) then
            write(IOW,*)
            write(IOW,*) ' The normalized difference between the numerical and exact stiffness matrix is ',err
            write(IOW,*)
            write(IOW,*) ' Note that ABAQUS format continuum elements use an approximate stiffness '
            write(IOW,*) ' For ABAQUS elements with a UMAT, consider using the '
            write(IOW,*) ' CHECK MATERIAL TANGENT key instead of CHECK STIFFNESS '
            write(IOW,*) ''
            write(IOW,*) ' If you are testing a user element, check the printed exact and numerical stiffness '
            write(IOW,*) ' for errors.   Small differences may be caused by rounding error.   Large ones '
            write(IOW,*) ' (more than a few %) are likely to be caused by an error in your coded stiffness '
            write(IOW,*) ' '
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

    endif

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

    if (allocated(abq_PREDEF)) deallocate(abq_PREDEF)
    if (allocated(abq_ADLMAG)) deallocate(abq_ADLMAG)
    if (allocated(abq_JDLTYP)) deallocate(abq_JDLTYP)
    if (allocated(abq_uel_bc_typ)) deallocate(abq_uel_bc_typ)
    if (allocated(abq_uel_bc_mag)) deallocate(abq_uel_bc_mag)
    if (allocated(abq_uel_bc_dmag)) deallocate(abq_uel_bc_dmag)

    return

end subroutine check_stiffness
