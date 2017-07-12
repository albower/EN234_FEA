subroutine check_tangent(material_flag)
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

    integer, intent(in)   :: material_flag

    ! Local Variables
    integer         :: lmn,n,j,k,iof
    integer         :: ix,iu,ns
    
    integer         :: status
    integer      :: mat_prop_index,n_mat_props
    
    real ( prec ) :: element_PNEWDT
    
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

    type (node), allocatable ::  local_nodes(:)

    character (len=80) material_name

    !     Subroutine to compare stiffness matrix calculated in a user subroutine with its numerical derivative
    !     It is activated by the CHECK STIFFNESS keyword


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

    !
    ! Find an element of type element_flag
    do lmn = 1, n_elements
        if (element_list(lmn)%material_index==material_flag) exit
    end do
    if (lmn>n_elements) then
        write(IOW,*) ' *** Error in subroutine check_stiffness *** '
        write(IOW,*) ' No element with material identifier ',material_flag,' was found in mesh '
        stop
    endif

    write(IOW,*) ' '
    write(IOW,*) ' '
    write(IOW,*) ' ===================== MATERIAL TANGENT STIFFNESS CHECK ============================'
    if (element_list(lmn)%flag==10002) then
        write(IOW,*)
        write(IOW,*) ' Testing 2D continuum element with UMAT '
        write(IOW,*) 'Material: ',material_namelist(element_list(lmn)%material_index)
    else if (element_list(lmn)%flag==10003) then
        write(IOW,*)
        write(IOW,*) ' Testing 3D continuum element with UMAT '
        write(IOW,*) 'Material: ',material_namelist(element_list(lmn)%material_index)
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

    if (element_list(lmn)%flag==10002) then            ! 2D continuum element

       element_PNEWDT = 1.d99

        mat_prop_index = material_list(element_list(lmn)%material_index)%prop_index
        n_mat_props = material_list(element_list(lmn)%material_index)%n_properties
        material_name = material_namelist(element_list(lmn)%material_index)(1:80)

        call check_material_2D(current_step_number,lmn, element_list(lmn)%flag, &
            element_list(lmn)%n_nodes, local_nodes(1:element_list(lmn)%n_nodes), &               ! Input variables
            n_mat_props, material_properties(mat_prop_index),material_name,  &               ! Input variables
            element_coords(1:ix),ix, &                                                       ! Input variables
            element_dof_increment(1:iu), element_dof_total(1:iu),iu,  &                                              ! Input variables
            ns, initial_state_variables(iof:iof+ns-1),&                                               ! Input variables
            updated_state_variables(iof:iof+ns-1),resid0(1:iu),element_stiffness(1:iu,1:iu),element_PNEWDT)


    else if (element_list(lmn)%flag==10003) then            ! 3D continuum element


        element_PNEWDT = 1.d99

        mat_prop_index = material_list(element_list(lmn)%material_index)%prop_index
        n_mat_props = material_list(element_list(lmn)%material_index)%n_properties
        material_name = material_namelist(element_list(lmn)%material_index)(1:80)

        call check_material_3D(current_step_number,lmn, element_list(lmn)%flag, &
            element_list(lmn)%n_nodes, local_nodes(1:element_list(lmn)%n_nodes), &               ! Input variables
            n_mat_props, material_properties(mat_prop_index),material_name,  &               ! Input variables
            element_coords(1:ix),ix, &                                                       ! Input variables
            element_dof_increment(1:iu), element_dof_total(1:iu),iu,  &                                              ! Input variables
            ns, initial_state_variables(iof:iof+ns-1),&                                               ! Input variables
            updated_state_variables(iof:iof+ns-1),resid0(1:iu),element_stiffness(1:iu,1:iu),element_PNEWDT)

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



    return

end subroutine check_tangent

!==========================SUBROUTINE continuum_element_static_3D ==============================
subroutine check_material_2D(current_step_number,lmn, element_identifier, &
    n_nodes, node_property_list, &               ! Input variables
    n_properties, material_properties, material_name, &               ! Input variables
    element_coords, length_coord_array, &                                                       ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                              ! Input variables
    n_state_variables, initial_state_variables, &                                               ! Input variables
    updated_state_variables,element_residual,element_stiffness,PNEWDT)                                   ! Output variables
    use Types
    use Globals
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : det33
    implicit none

    integer, intent( in )         :: current_step_number                                    ! Current step
    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                            ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element


    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: material_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine

    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness
    real( prec ), intent( inout ) :: PNEWDT                                                 ! User specified time increment


    character (len=80), intent( in ) :: material_name                                       ! Name of the material

    ! Local Variables
    integer      :: n_points,kint
    integer      :: n_state_vars_per_intpt, iof

    real (prec)  ::  stress0(3,3)                             ! Initial stress
    real (prec)  ::  xintpt(3)                                ! Coords of integration point
    real (prec)  ::  dxidx(2,2), determinant                  ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  F0(3,3)                                  ! Def grad at start of increment, global basis
    real (prec)  ::  F1(3,3)                                  ! Def grad at end of increment, global basis
    real (prec)  ::  dF(3,3)                                  ! Def gradient increment, global basis
    real (prec)  ::  Fmidinv(3,3)                             ! Inverse of mid point def gradient
    real (prec)  ::  F1inv(3,3)                               ! Inverse of updated def gradient
    real (prec)  ::  J0                                       ! det(F0)
    real (prec)  ::  J1                                       ! det(F1)
    real (prec)  ::  J0bar                                    ! Vol averaged det(F0)
    real (prec)  ::  J1bar                                    ! Vol averaged det(F1)
    real (prec)  ::  Jmid                                     ! Vol averaged det(Fmid)
    real (prec)  ::  deps(3,3)                                ! Strain increment in global basis
    real (prec)  ::  deps_perturbed(3,3)                      ! Perturbed strain increment for numerical deriv
    real (prec)  ::  dW(3,3)                                  ! Spin increment components in global basis
    real (prec)  ::  dR(3,3)                                  ! Rotation increment components in global basis
    real (prec)  ::  strain0(3,3)                             ! Strain at start of increment; updated for rigid rotation
    real (prec)  ::  el_vol
    real (prec)  ::  temp33a(3,3)                             ! Workspace array
    real (prec)  ::  temp33b(3,3)                             ! Workspace array
    real (prec)  ::  temp33c(3,3)                             ! Workspace array

    real (prec)  ::  dummy                                    ! Dummy variable
    real (prec)  ::  predef(1)                                ! Dummy variable
    real (prec)  ::  dpred(1)                                 ! Dummy variable
    real (prec)  ::  Dcorot(4,4)                              ! Corotational material tangent
    real (prec)  ::  Dcorot_numerical(4,4)                    ! Numerical corotataional material tangent
    real (prec)  ::  Dtherm(4)                                ! Tangent wrt temp change
    !
    real (prec)  ::  stress_initial(4)                        ! Stress
    real (prec)  ::  stress(4)                                ! Stress (in global basis) stored as vector
    real (prec)  ::  strain(4)                                ! Strain (global basis) stored as vector
    real (prec)  ::  strainInc(4)                             ! Strain increment (in global basis) stored as vector
    real (prec)  ::  SSE,  SPD, SCD                           ! Energy dissipation variables
    real (prec)  ::  charLength                               ! Characteristic element length
    real (prec)  ::  rpl                                      ! Rate of plastic heat generation

    real (prec)  ::  drplde(4)                                ! Derivative of plastic work wrt strain (unused)
    real (prec)  ::  drpldt
    real (prec)  ::  abq_TIME(2)

    integer :: ie
    integer :: i,j
    !
    !   This implements a 2D B-Bar continuum element with incremental stress update.
    !   The constitutive response of the solid is specified by an ABAQUS UMAT user subroutine.
    !   Currently coded for plane strain only
    !
    !   The state variable storage for 2D continuum elements is
    !   s11,s22,s33,s12
    !   e11,e22,e33,e12
    !   specific internal energy, specific plastic dissipation, specific creep dissipation,
    !   user defined variables
    !

    x = reshape(element_coords,(/2,length_coord_array/2/))

    if (n_nodes == 1) n_points = 1
    if (n_nodes == 4) n_points = 4
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 8) n_points = 9

    n_state_vars_per_intpt = n_state_variables/n_points

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0

    !   Volume averaged element variables
    J0bar = 0.d0
    J1bar = 0.d0
    el_vol = 0.d0
    F0 = 0.d0
    F1 = 0.d0
    dF = 0.d0
    F0(3,3) = 1.d0
    F1(3,3) = 1.d0
    do kint = 1,n_points

        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        xintpt(1:2) = matmul(x(1:2,1:n_nodes),N(1:n_nodes))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        do i = 1,2
            ie = 2*(n_nodes-1)+i
            F0(i,1:2) = matmul(dof_total(i:ie:2) ,dNdx(1:n_nodes,1:2))
            F0(i,i) = F0(i,i) + 1.d0
            dF(i,1:2) = matmul(dof_increment(i:ie:2),dNdx(1:n_nodes,1:2))
        end do
        F1 = F0 + dF
        call invert_small(F1,F1inv,J1)
        dNdx(1:n_nodes,1:2) = matmul(dNdx(1:n_nodes,1:2),F1inv(1:2,1:2))
        el_vol = el_vol + w(kint)*determinant
        J0bar = J0bar + det33(F0)*w(kint)*determinant
        J1bar = J1bar + J1*w(kint)*determinant
    end do
    J0bar = J0bar/el_vol
    J1bar = J1bar/el_vol
    charLength = el_vol**(1.d0/2.d0)

    kint = 1              ! Check material at the first integration point
    call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
    dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
    xintpt(1:2) = matmul(x(1:2,1:n_nodes),N(1:n_nodes))
    xintpt(3) = 0.d0
    call invert_small(dxdxi,dxidx,determinant)
    dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
    !
    !       Def gradients
    !
    do i = 1,2
        ie = 2*(n_nodes-1)+i
        F0(i,1:2) = matmul(dof_total(i:ie:2),dNdx(1:n_nodes,1:2))
        F0(i,i) = F0(i,i) + 1.d0
        dF(i,1:2) = matmul(dof_increment(i:ie:2),dNdx(1:n_nodes,1:2))
    end do
    F1 = F0 + dF
    call invert_small(F1,F1inv,J1)
    dNdx(1:n_nodes,1:2) = matmul(dNdx(1:n_nodes,1:2),F1inv(1:2,1:2))  !       Spatial shape function derivatives
    J0 = det33(F0)
    F0(1:2,1:2) = F0(1:2,1:2)*(J0bar/J0)**(1.d0/2.d0)
    F1(1:2,1:2) = F1(1:2,1:2)*(J1bar/J1)**(1.d0/2.d0)
    dF = F1-F0
    call invert_small(0.5d0*(F1+F0),Fmidinv,Jmid)

    !       Strain and spin increment
    deps = matmul(dF,Fmidinv)
    dW = 0.5d0*(deps - transpose(deps))
    deps = 0.5d0*(deps + transpose(deps))

    !     Hughes-Winget rotation increment
    temp33a = eye3_D + 0.5d0*dW
    temp33b = eye3_D - 0.5d0*dW
    call invert_small(temp33b,temp33c,dummy)
    dR = matmul(temp33c,temp33a)

    iof = n_state_vars_per_intpt*(kint-1) + 1
    stress0 = reshape([initial_state_variables(iof),initial_state_variables(iof+3),0.d0, &
        initial_state_variables(iof+3),initial_state_variables(iof+1),0.d0, &
        0.d0,0.d0,initial_state_variables(iof+2)], &
        shape(stress0))
    stress0 = matmul(dR,matmul(stress0,transpose(dR)))    !

    stress(1:4) = [stress0(1,1),stress0(2,2),stress0(3,3),stress0(1,2)]            ! UMAT format stress vector

    strain0 = reshape( &
        [initial_state_variables(iof+4),0.5d0*initial_state_variables(iof+7),0.d0, &
        0.5d0*initial_state_variables(iof+7),initial_state_variables(iof+5),0.d0, &
        0.d0,0.d0,initial_state_variables(iof+6)], &
        shape(strain0))
    strain0 =  matmul(dR,matmul(strain0,transpose(dR)))
    strain(1:4) = [strain0(1,1),strain0(2,2),strain0(3,3),strain0(1,2)]            ! UMAT format strain vector

    strainInc(1:4) = [deps(1,1),deps(2,2),deps(3,3),2.d0*deps(1,2)]                    ! UMAT format strain increment vector

    if (n_state_vars_per_intpt>10)   updated_state_variables(iof+11:iof+n_state_vars_per_intpt) = &
        initial_state_variables(iof+11:iof+n_state_vars_per_intpt)
    SSE = initial_state_variables(iof+8)
    SPD = initial_state_variables(iof+9)
    SCD = initial_state_variables(iof+10)

    abq_time(1:2) = TIME

    stress_initial = stress
    Call UMAT(stress_initial,updated_state_variables(iof+11:iof+n_state_vars_per_intpt),Dcorot,SSE,SPD,SCD, &
        RPL,Dtherm,DRPLDE,DRPLDT, &
        strain,strainInc,abq_time,DTIME,BTEMP,BTINC,predef,dpred,material_name, &
        3,1,4,n_state_vars_per_intpt-11,material_properties,n_properties,xintpt,dR,PNEWDT, &
        charLength,F0,F1,lmn,kint,0,0,1,current_step_number)


    do i = 1,4
        strainInc(i) = strainInc(i) + 1.d-07
        deps_perturbed(1,1:3) =  [strainInc(1),0.5d0*strainInc(4),0.d0]
        deps_perturbed(2,1:3) =  [0.5d0*strainInc(4),strainInc(2),0.d0]
        deps_perturbed(3,1:3) =  [ 0.d0, 0.d0, strainInc(3)]

        temp33a = eye3_D - 0.5d0*(deps_perturbed+dW)
        call invert_small(temp33a,temp33b,dummy)
        temp33c = eye3_D + 0.5d0*(deps_perturbed+dW)
        F1 = matmul(temp33b,matmul(temp33c,F0))
        if (n_state_vars_per_intpt>10)   updated_state_variables(iof+11:iof+n_state_vars_per_intpt) = &
            initial_state_variables(iof+11:iof+n_state_vars_per_intpt)
        stress(1:4) = [stress0(1,1),stress0(2,2),stress0(3,3),stress0(1,2)]
        Call UMAT(stress,updated_state_variables(iof+11:iof+n_state_vars_per_intpt),Dcorot,SSE,SPD,SCD, &
            RPL,Dtherm,DRPLDE,DRPLDT, &
            strain,strainInc,abq_time,DTIME,BTEMP,BTINC,predef,dpred,material_name, &
            3,1,4,n_state_vars_per_intpt-11,material_properties,n_properties,xintpt,dR,PNEWDT, &
            charLength,F0,F1,lmn,kint,0,0,1,current_step_number)
        Dcorot_numerical(1:4,i) = (stress(1:4) - stress_initial(1:4))/1.d-07

        write(IOW,*)
        write(IOW,*) ' Column ',i

        do j = 1,4
            write(IOW,2000) j,Dcorot(j,i),Dcorot_numerical(j,i)
2000        format( ' Row ',i4,' Stiffness ',d15.5,' Numerical deriv ',d15.5 )
        end do

        strainInc(i) = strainInc(i) - 1.d-07

    end do



    return
end subroutine check_material_2D


!
!==========================SUBROUTINE continuum_element_static_3D ==============================
subroutine check_material_3D(current_step_number,lmn, element_identifier, &
    n_nodes, node_property_list, &               ! Input variables
    n_properties, material_properties, material_name, &               ! Input variables
    element_coords, length_coord_array, &                                                       ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                              ! Input variables
    n_state_variables, initial_state_variables, &                                               ! Input variables
    updated_state_variables,element_residual,element_stiffness,PNEWDT)                                   ! Output variables
    use Types
    use Globals
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : det33
    use Element_Utilities, only : polardecomp
    implicit none

    integer, intent( in )         :: current_step_number                                    ! Current step
    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                            ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element


    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: material_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine

    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness
    real( prec ), intent( inout ) :: PNEWDT                                                 ! User specified time increment


    character (len=80), intent( in ) :: material_name                                       ! Name of the material

    ! Local Variables
    integer      :: n_points,kint
    integer      :: n_state_vars_per_intpt, iof


    real (prec)  ::  stress0(3,3)                             ! Initial stress
    real (prec)  ::  xintpt(3)                                ! Coords of integration point
    real (prec)  ::  dxidx(3,3), determinant                  ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  F0(3,3)                                  ! Def grad at start of increment, global basis
    real (prec)  ::  F1(3,3)                                  ! Def grad at end of increment, global basis
    real (prec)  ::  dF(3,3)                                  ! Def gradient increment, global basis
    real (prec)  ::  Fmidinv(3,3)                             ! Inverse of mid point def gradient
    real (prec)  ::  F1inv(3,3)                               ! Inverse of updated def gradient
    real (prec)  ::  J0                                       ! det(F0)
    real (prec)  ::  J1                                       ! det(F1)
    real (prec)  ::  J0bar                                    ! Vol averaged det(F0)
    real (prec)  ::  J1bar                                    ! Vol averaged det(F1)
    real (prec)  ::  Jmid                                     ! Vol averaged det(Fmid)
    real (prec)  ::  deps(3,3)                                ! Strain increment in global basis
    real (prec)  ::  deps_perturbed(3,3)                      ! Perturbed strain for numerical derivative
    real (prec)  ::  dW(3,3)                                  ! Spin increment components in global basis
    real (prec)  ::  dR(3,3)                                  ! Rotation increment components in global basis
    real (prec)  ::  strain0(3,3)                             ! Strain at start of increment; updated for rigid rotation
    real (prec)  ::  el_vol
    real (prec)  ::  temp33a(3,3)                             ! Workspace array
    real (prec)  ::  temp33b(3,3)                             ! Workspace array
    real (prec)  ::  temp33c(3,3)                             ! Workspace array

    real (prec)  ::  dummy                                    ! Dummy variable
    real (prec)  ::  Dcorot(6,6)                              ! Corotational material tangent
    real (prec)  ::  Dcorot_numerical(6,6)                    ! Spin tangent
    real (prec)  ::  Dtherm(6)                                ! Thermal tangent
!

    real (prec)  ::  stress_initial(6)                        ! Stress
    real (prec)  ::  stress(6)                                ! Stress (in global basis) stored as vector
    real (prec)  ::  strain(6)                                ! Strain (global basis) stored as vector
    real (prec)  ::  strainInc(6)                             ! Strain increment (in global basis) stored as vector
    real (prec)  ::  SSE,  SPD, SCD                           ! Energy dissipation variables
    real (prec)  ::  charLength                               ! Characteristic element length
    real (prec)  ::  rpl                                      ! Rate of plastic heat generation

    real (prec)  ::  drplde(6)                                ! Derivative of plastic work wrt strain (unused)
    real (prec)  ::  drpldt
    real (prec)  ::  abq_TIME(2)
    real (prec)  ::  predef(1)                                ! Dummy argument for predefined field
    real (prec)  ::  dpred(1)                                 ! Dummy argument for predefined field increment


    integer :: ie
    integer :: i,j
!
!   This implements a 3D B-Bar continuum element with incremental stress update.
!   The constitutive response of the solid is specified by an ABAQUS UMAT user subroutine.
!
!   The state variable storage for 3D continuum elements is
!   s11,s22,s33,s12,s13,s23,
!   e11,e22,e33,e12,e13,e23,
!   specific internal energy, specific plastic dissipation, specific creep dissipation,
!   user defined variables
!

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    n_state_vars_per_intpt = n_state_variables/n_points

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0

    !   Volume averaged element variables
    J0bar = 0.d0
    J1bar = 0.d0
    el_vol = 0.d0
    do kint = 1,n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        xintpt = matmul(x(1:3,1:n_nodes),N(1:n_nodes))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        do i = 1,3
            ie = 3*(n_nodes-1)+i
            F0(i,1:3) = matmul(dof_total(i:ie:3) ,dNdx(1:n_nodes,1:3))
            F0(i,i) = F0(i,i) + 1.d0
            dF(i,1:3) = matmul(dof_increment(i:ie:3),dNdx(1:n_nodes,1:3))
        end do
!        dF = 0.d0
        F1 = F0 + dF
        call invert_small(F1,F1inv,J1)
        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),F1inv)
        el_vol = el_vol + w(kint)*determinant
        J0bar = J0bar + det33(F0)*w(kint)*determinant
        J1bar = J1bar + J1*w(kint)*determinant
    end do
    J0bar = J0bar/el_vol
    J1bar = J1bar/el_vol
    charLength = el_vol**(1.d0/3.d0)


    kint = 1            ! Check material tangent at first integration point
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!
!       Def gradients
!
        do i = 1,3
            ie = 3*(n_nodes-1)+i
            F0(i,1:3) = matmul(dof_total(i:ie:3),dNdx(1:n_nodes,1:3))
            F0(i,i) = F0(i,i) + 1.d0
            dF(i,1:3) = matmul(dof_increment(i:ie:3),dNdx(1:n_nodes,1:3))
        end do
!        dF = 0.d0
        F1 = F0 + dF
        call invert_small(F1,F1inv,J1)
        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),F1inv)  !       Spatial shape function derivatives
        J0 = det33(F0)
        F0 = F0*(J0bar/J0)**(1.d0/3.d0)
        F1 = F1*(J1bar/J1)**(1.d0/3.d0)
        dF = F1-F0
        call invert_small(0.5d0*(F1+F0),Fmidinv,Jmid)


!       Strain and spin increment
        deps = matmul(dF,Fmidinv)
        dW = 0.5d0*(deps - transpose(deps))
        deps = 0.5d0*(deps + transpose(deps))

!     Hughes-Winget rotation increment
        temp33a = eye3_D + 0.5d0*dW
        temp33b = eye3_D - 0.5d0*dW
        call invert_small(temp33b,temp33c,dummy)
        dR = matmul(temp33c,temp33a)

        iof = n_state_vars_per_intpt*(kint-1) + 1
        stress0 = reshape([initial_state_variables(iof),initial_state_variables(iof+3),initial_state_variables(iof+4), &
                           initial_state_variables(iof+3),initial_state_variables(iof+1),initial_state_variables(iof+5), &
                           initial_state_variables(iof+4),initial_state_variables(iof+5),initial_state_variables(iof+2)], &
                           shape(stress0))
        stress0 = matmul(dR,matmul(stress0,transpose(dR)))    !

        stress(1:6) = [stress0(1,1),stress0(2,2),stress0(3,3), &
                       stress0(1,2),stress0(1,3),stress0(2,3)]            ! UMAT format stress vector


        strain0 = reshape( &
               [initial_state_variables(iof+6),0.5d0*initial_state_variables(iof+9),0.5d0*initial_state_variables(iof+10), &
                0.5d0*initial_state_variables(iof+9),initial_state_variables(iof+7),0.5d0*initial_state_variables(iof+11), &
                0.5d0*initial_state_variables(iof+10),0.5d0*initial_state_variables(iof+11),initial_state_variables(iof+8)], &
                           shape(strain0))
        strain0 =  matmul(dR,matmul(strain0,transpose(dR)))
        strain(1:6) = [strain0(1,1),strain0(2,2),strain0(3,3), &
                       strain0(1,2),strain0(1,3),strain0(2,3)]            ! UMAT format strain vector

        strainInc(1:6) = [deps(1,1),deps(2,2),deps(3,3), &
                        2.d0*deps(1,2),2.d0*deps(1,3),2.d0*deps(2,3)]                    ! UMAT format strain increment vector


        if (n_state_vars_per_intpt>15)   updated_state_variables(iof+15:iof+n_state_vars_per_intpt) = &
                                         initial_state_variables(iof+15:iof+n_state_vars_per_intpt)
        SSE = initial_state_variables(iof+12)
        SPD = initial_state_variables(iof+13)
        SCD = initial_state_variables(iof+14)

        abq_time(1:2) = TIME

        temp33a = eye3_D - 0.5d0*(deps+dW)
        call invert_small(temp33a,temp33b,dummy)
        temp33c = eye3_D + 0.5d0*(deps+dW)
        F1 = matmul(temp33b,matmul(temp33c,F0))

        stress_initial = stress
        Call UMAT(stress_initial,updated_state_variables(iof+15:iof+n_state_vars_per_intpt),Dcorot,SSE,SPD,SCD, &
                  RPL,Dtherm,DRPLDE,DRPLDT, &
                  strain,strainInc,abq_time,DTIME,BTEMP,BTINC,predef,dpred,material_name, &
                  3,3,6,n_state_vars_per_intpt-15,material_properties,n_properties,xintpt,dR,PNEWDT, &
                  charLength,F0,F1,lmn,kint,0,0,1,current_step_number)

    do i = 1,6
        strainInc(i) = strainInc(i) + 1.d-07
        deps_perturbed(1,1:3) =  [strainInc(1),0.5d0*strainInc(4),0.5d0*strainInc(5)]
        deps_perturbed(2,1:3) =  [0.5d0*strainInc(4),strainInc(2),0.5d0*strainInc(6)]
        deps_perturbed(3,1:3) =  [ 0.5d0*strainInc(5), 0.5d0*strainInc(6) , strainInc(3)]

        temp33a = eye3_D - 0.5d0*(deps_perturbed+dW)
        call invert_small(temp33a,temp33b,dummy)
        temp33c = eye3_D + 0.5d0*(deps_perturbed+dW)
        F1 = matmul(temp33b,matmul(temp33c,F0))
        if (n_state_vars_per_intpt>14)   updated_state_variables(iof+15:iof+n_state_vars_per_intpt) = &
                                         initial_state_variables(iof+15:iof+n_state_vars_per_intpt)
        stress(1:6) = [stress0(1,1),stress0(2,2),stress0(3,3), &
                       stress0(1,2),stress0(1,3),stress0(2,3)]
        Call UMAT(stress,updated_state_variables(iof+15:iof+n_state_vars_per_intpt),Dcorot,SSE,SPD,SCD, &
                  RPL,Dcorot,DRPLDE,DRPLDT, &
                  strain,strainInc,abq_time,DTIME,BTEMP,BTINC,predef,dpred,material_name, &
                  3,3,6,n_state_vars_per_intpt-15,material_properties,n_properties,xintpt,dR,PNEWDT, &
                  charLength,F0,F1,lmn,kint,0,0,1,current_step_number)
        Dcorot_numerical(1:6,i) = (stress(1:6) - stress_initial(1:6))/1.d-07

        write(IOW,*)
        write(IOW,*) ' Column ',i

        do j = 1,6
            write(IOW,2000) j,Dcorot(j,i),Dcorot_numerical(j,i)
2000        format( ' Row ',i4,' Stiffness ',d15.5,' Numerical deriv ',d15.5 )
        end do

        strainInc(i) = strainInc(i) - 1.d-07

    end do


    return
end subroutine check_material_3D

