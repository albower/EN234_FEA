!
!     Subroutines to assemble element level matrices for FEA analysis

!=========================== subroutine user_element_stiffness ===================
subroutine user_element_static(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, n_int_properties, int_element_properties, &                ! Input variables
    element_coords, length_coord_array, &                                                        ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables
    use Types
    use ParamIO
    use Globals, only: TIME, DTIME                  ! Total analysis time and time increment
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: n_int_properties                                       ! # integer valued properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    integer, intent( in )         :: int_element_properties(n_int_properties)               ! Integer valued element properties

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout )  :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)


    !     Element stiffness routine

    
    element_stiffness = 0.d0
    element_residual = 0.d0

    fail = .false.

    updated_state_variables = initial_state_variables

    if ( element_identifier == 1001 ) then              ! Basic fully integrated 3D linear elastic element

        call el_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, n_int_properties, int_element_properties,    &               ! Input variables
    element_coords, length_coord_array, &                                                          ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables


    else if ( element_identifier ==0) then           ! Stub for a new element
  
        call new_user_element_static(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, n_int_properties, int_element_properties,    &               ! Input variables
    element_coords, length_coord_array, &                                                        ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables
  
  
  
    else
        write (IOW, 99001) element_identifier
        stop

    end if

99001 format ( // ' **** ERROR DETECTED IN SUBROUTINE user_element_static ****'/  &
        '   Invalid element type was specified '/, &
        '   Current element types are: '/  &
        '     IEP=1001     Basic fully integrated 3D linear elastic element       '/&
        '    Subroutine called with IEP = ', I10)

end subroutine user_element_static


!=========================== subroutine user_element_dynamic ===================
subroutine user_element_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, n_int_properties, int_element_properties, &                        ! Input variables
    element_coords, length_coord_array, &                                                                ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                  ! Input variables
    n_state_variables, initial_state_variables, &                                                  ! Input variables
    updated_state_variables,element_residual,element_deleted)      ! Output variables
    use Types
    use ParamIO
    use Globals, only: TIME, DTIME                  ! Total analysis time and time increment
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: n_int_properties                                       ! # integer valued properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    integer, intent( in )         :: int_element_properties(n_int_properties)               ! Integer valued element properties

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
              
    real( prec ), intent( inout )  :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete an element
    !     Element force routine for explicit dynamic analysis

    
    element_residual = 0.d0
    element_deleted = .false.

    updated_state_variables = initial_state_variables

    if ( element_identifier == 1001 ) then              ! Basic fully integrated 3D linear elastic element

        call el_linelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
            n_properties, element_properties, n_int_properties, int_element_properties, &                  ! Input variables
            element_coords, length_coord_array, &                                                          ! Input variables
            dof_increment, dof_total, length_dof_array,  &                                                 ! Input variables
            n_state_variables, initial_state_variables, &                                                  ! Input variables
            updated_state_variables,element_residual,element_deleted)                                      ! Output variables

    else if ( element_identifier ==0) then               ! Stub for a new element
  
        call new_user_element_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
            n_properties, element_properties, n_int_properties, int_element_properties, &               ! Input variables
            element_coords, length_coord_array, &                                                       ! Input variables
            dof_increment, dof_total, length_dof_array,  &                                                 ! Input variables
            n_state_variables, initial_state_variables, &                                                  ! Input variables
            updated_state_variables,element_residual,element_deleted)                                      ! Output variables
  
    else
        write (IOW, 99001) element_identifier
        stop

    end if

99001 format ( // ' **** ERROR DETECTED IN SUBROUTINE user_element_dynamic ****'/  &
        '   Invalid element type was specified '/, &
        '   Current element types are: '/  &
        '     IEP=1001     Basic fully integrated 3D linear elastic element       '/&
        '    Subroutine called with IEP = ', I10)

end subroutine user_element_dynamic

subroutine user_element_fieldvariables(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, n_int_properties, int_element_properties, &                        ! Input variables
    element_coords, length_coord_array, &                                                                ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                        ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                                ! Input variables
    n_field_variables,field_variable_names, &                                                            ! Field variable definition
    nodal_fieldvariables)                                                                                ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: n_int_properties                                       ! # integer valued properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! No. projected field variables

    integer, intent( in )         :: int_element_properties(n_int_properties)               ! Integer valued element properties


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

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Element stiffness (ROW,COLUMN)



    if ( element_identifier == 1001 ) then              ! Basic fully integrated 3D linear elastic element

        call fieldvars_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &         ! Input variables
            n_properties, element_properties, n_int_properties, int_element_properties, &               ! Input variables
            element_coords, length_coord_array, &                                                       ! Input variables
            dof_increment, dof_total, length_dof_array,  &                                              ! Input variables
            n_state_variables, initial_state_variables,updated_state_variables, &                       ! Input variables
            n_field_variables,field_variable_names, &                                                   ! Field variable definition
            nodal_fieldvariables)      ! Output variables

  
    else  if ( element_identifier == 0 ) then              ! Stub for a new element

        call new_user_element_fieldvariables(lmn, element_identifier, n_nodes, node_property_list, &         ! Input variables
            n_properties, element_properties, n_int_properties, int_element_properties, &               ! Input variables
            element_coords, length_coord_array, &                                                       ! Input variables
            dof_increment, dof_total, length_dof_array,  &                                              ! Input variables
            n_state_variables, initial_state_variables,updated_state_variables, &                       ! Input variables
            n_field_variables,field_variable_names, &                                                   ! Field variable definition
            nodal_fieldvariables)      ! Output variables
 
    else

        write (IOW, 99001) element_identifier
    stop

99001 format ( // ' **** ERROR DETECTED IN SUBROUTINE user_element_fieldvariables ****'/  &
        '   Invalid element type was specified '/, &
        '   Current element types are: '/  &
        '     IEP=1001     Basic fully integrated 3D linear elastic element       '/&
        '    Subroutine called with IEP = ', I10)



end if



end subroutine user_element_fieldvariables


subroutine user_element_lumped_mass(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    density,n_properties, element_properties, n_int_properties, int_element_properties, &                        ! Input variables
    element_coords, length_coord_array, &                                                                ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                     ! Input variables
    n_state_variables, initial_state_variables, &                                                     ! Input variables
    updated_state_variables,element_lumped_mass)      ! Output variables
    use Types
    use ParamIO
    use Globals, only: TIME, DTIME                  ! Total analysis time and time increment
    use Mesh, only : node
    use Element_Utilities, only : N3 => shape_functions_3D
    use Element_Utilities, only : dNdxi3 => shape_function_derivatives_3D
    use Element_Utilities, only : dNdx3 => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi3 => integrationpoints_3D, w3 => integrationweights_3D
    use Element_Utilities, only : dxdxi3 => jacobian_3D
    use Element_Utilities, only : N2 => shape_functions_2D
    use Element_Utilities, only : dNdxi2 => shape_function_derivatives_2D
    use Element_Utilities, only : dNdx2 => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi2 => integrationpoints_2D, w2 => integrationweights_2D
    use Element_Utilities, only : dxdxi2 => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions


    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: n_int_properties                                       ! # integer valued properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    integer, intent( in )         :: int_element_properties(n_int_properties)               ! Integer valued element properties

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element
    real( prec ), intent( in )    :: density                                                ! Density value for element
    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
              
    real( prec ), intent( inout )  :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_lumped_mass(length_dof_array)                     ! Element residual force (ROW)


    real (prec) :: x2D(2,length_coord_array/2)
    real (prec) :: x3D(3,length_coord_array/3)
    real (prec) :: determinant

    logical :: twoD, threeD
    integer :: j,k,kint,n_points,irow

    !     Element lumped mass routine - currently coded to compute lumped mass for standard element types with uniform density

    
    element_lumped_mass = 0.d0

    twoD = .false.
    threeD = .false.
    do j = 1,n_nodes
        if (node_property_list(j)%n_coords==2) twoD = .true.
        if (node_property_list(j)%n_coords==3) threeD = .true.
    end do
          
    if (twoD.and.threeD) return               ! Skip nonstandard element
    if (.not.twoD .and. .not.threeD) return   ! Skip nonstandard element

    if (twoD) then
        if (n_nodes>9) return                                !Nonstandard element
        x2D = reshape(element_coords,(/2,length_coord_array/2/))
        n_points = 0
        if (n_nodes == 3) n_points = 4
        if (n_nodes == 4) n_points = 4
        if (n_nodes == 6) n_points = 7
        if (n_nodes == 8) n_points = 9
        if (n_nodes == 9) n_points = 9
        if (n_points==0) return                             ! Nonstandard element
        call initialize_integration_points(n_points, n_nodes, xi2, w2)
             
        do kint = 1,n_points
            call calculate_shapefunctions(xi2(1:2,kint),n_nodes,N2,dNdxi2)
            dxdxi2 = matmul(x2D(1:2,1:n_nodes),dNdxi2(1:n_nodes,1:2))
            determinant = dxdxi2(1,1)*dxdxi2(2,2) - dxdxi2(2,1)*dxdxi2(1,2)
            !             Lumped mass matrix computed using row sum of consistent mass matrix
            irow = 0
            do j = 1,n_nodes
                do k = 1,node_property_list(j)%n_dof
                    irow = irow + 1
                    element_lumped_mass(irow) = element_lumped_mass(irow) + N2(j)*sum(N2(1:n_nodes))*determinant*w2(kint)
                end do
            end do
        end do
    else if (threeD) then
        x3D = reshape(element_coords,(/3,length_coord_array/3/))
        n_points = 0
        if (n_nodes == 4) n_points = 4
        if (n_nodes == 10) n_points = 5
        if (n_nodes == 8) n_points = 27
        if (n_nodes == 20) n_points = 27
        if (n_points==0) return                                 ! Nonstandard element
        call initialize_integration_points(n_points, n_nodes, xi3, w3)
        do kint = 1,n_points
            call calculate_shapefunctions(xi3(1:3,kint),n_nodes,N3,dNdxi3)
            dxdxi3 = matmul(x3D(1:3,1:n_nodes),dNdxi3(1:n_nodes,1:3))
            determinant =   dxdxi3(1,1)*dxdxi3(2,2)*dxdxi3(3,3)  &
                - dxdxi3(1,1)*dxdxi3(2,3)*dxdxi3(3,2)  &
                - dxdxi3(1,2)*dxdxi3(2,1)*dxdxi3(3,3)  &
                + dxdxi3(1,2)*dxdxi3(2,3)*dxdxi3(3,1)  &
                + dxdxi3(1,3)*dxdxi3(2,1)*dxdxi3(3,2)  &
                - dxdxi3(1,3)*dxdxi3(2,2)*dxdxi3(3,1)
            !             Lumped projection matrix computed  using row sum method
            irow = 0
            do j = 1,n_nodes
                do k = 1,node_property_list(j)%n_dof
                    irow = irow + 1
                    element_lumped_mass(irow) = element_lumped_mass(irow) + N3(j)*sum(N3(1:n_nodes))*determinant*w3(kint)
                end do
            end do
        end do
    endif

    element_lumped_mass = element_lumped_mass*density

end subroutine user_element_lumped_mass



