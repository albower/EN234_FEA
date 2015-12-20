!
!     Subroutines to assemble element level matrices for FEA analysis

!=========================== subroutine user_element_stiffness ===================
!subroutine bbar_hyperelastic(lmn, element_identifier, n_nodes, node_property_list, &          ! Input variables
!    n_properties, element_properties,element_coords, dof_increment, dof_total,  &                  ! Input variables
!    n_state_variables, initial_state_variables, &                                                  ! Input variables
!    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables
!    use Types
!    use ParamIO
!    use Globals, only: TIME, DTIME                  ! Total analysis time and time increment
!    use Mesh, only : node
!    use User_Subroutine_Storage
!
!    use Element_Utilities, only : N => shape_functions_3D
!    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
!    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
!    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
!    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
!    use Element_Utilities, only : dxdxi => jacobian_3D
!    use Element_Utilities, only : initialize_integration_points
!    use Element_Utilities, only : calculate_shapefunctions
!    use Element_Utilities, only : invert_small
!
!    implicit none
!
!    integer, intent( in )         :: lmn                                                    ! Element number
!    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
!    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
!    integer, intent( in )         :: n_properties                                           ! # properties for the element
!    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
!
!    type (node), intent( in )     :: node_property_list(length_node_array)                  ! Data structure describing storage for nodal variables - see below
!    !  type node
!    !      sequence
!    !      integer :: flag                          ! Integer identifier
!    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
!    !      integer :: n_coords                      ! Total no. coordinates for the node
!    !      integer :: dof_index                     ! Index of first DOF in dof array
!    !      integer :: n_dof                         ! Total no. of DOF for node
!    !   end type node
!    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element
!
!    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
!    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
!    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment
!
!    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
!    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
!
!    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
!    real( prec ), intent( inout )  :: updated_state_variables(n_state_variables)             ! State variables at end of time step
!    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
!    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
!
!
!    real (prec)  ::  stress(6)                                ! Stress vector contains [s11, s22, s33, s12, s13, s23] Kirchhoff stress
!    real (prec)  ::  D(6,6)                                   ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
!    real (prec)  ::  B(3,3)                    ! strain = B*(dof_total+dof_increment)
!    real (prec)  ::  Bbar(6,length_dof_array)                 ! locking free strain = Bbar*(dof_total+dof_increment)
!    real (prec)  ::  Bstar(9,length_dof_array)                ! locking free strain = Bbar*(dof_total+dof_increment)
!    real (prec)  ::  G(6,9)                                   ! DB/dF * F^T
!    real (prec)  ::  dxidx(3,3), determinant                  ! Jacobian inverse and determinant
!    real (prec)  ::  x(3,length_coord_array/3)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node
!    real (prec)  ::  F(3,3)                                   ! Def gradient
!    real (prec)  ::  Fbar(3,3)                                ! Locking free deformation gradient
!    real (prec)  ::  Finv(3,3)                                ! Inverse of def gradient
!    real (prec)  ::  J                                        ! det (F)
!    real (prec)  ::  Jbar                                     ! Vol averaged det(F)
!    real (prec)  ::  dNdxvec(3*n_nodes)                ! dN/dx stored as a vector
!    real (prec)  ::  Pbar(length_dof_array,length_dof_array)  ! d^2 Jbar/ du_i^a du_k^b  matrix
!    real (prec)  ::  Pvec(length_dof_array)
!    real (prec)  ::  Pmat(length_dof_array,length_dof_array)
!    real (prec)  ::  P(length_dof_array,length_dof_array)     !
!    real (prec)  ::  S(3,length_dof_array/3)
!    real (prec)  ::  Svec(length_dof_array)
!    real (prec)  ::  Smat(length_dof_array,length_dof_array)
!    real (prec)  ::  el_vol
!    real (prec)  ::  Fscale
!    real (prec)  ::  press
!
!    integer :: n_points        ! No. integration points
!    integer :: kint            ! Integration point index
!    integer :: ie
!    integer :: i
!
!
!    fail = .false.
!
!    x = reshape(element_coords,(/3,length_coord_array/3/))
!
!    if (n_nodes == 4) n_points = 1
!    if (n_nodes == 10) n_points = 4
!    if (n_nodes == 8) n_points = 8
!    if (n_nodes == 20) n_points = 27
!
!    call initialize_integration_points(n_points, n_nodes, xi, w)
!
!    element_residual = 0.d0
!    element_stiffness = 0.d0
!
!
!    !   Volume averaged element variables
!    Jbar = 0.d0
!    dNbardx = 0.d0
!    Pbar = 0.d0
!    el_vol = 0.d0
!    do kint = 1,n_points
!
!        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!        call invert_small(dxdxi,dxidx,determinant)
!        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        do i = 1,3
!            ie = 3*(n_nodes-1)+i
!            F(i,1:3) = matmul(dof_increment(i:ie:3)+dof_total(i:ie:3),dNdx(1:n_nodes,1:3))
!            F(i,i) = F(i,i) + 1.d0
!        end do
!        call invert_small(F,Finv,J)
!        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),Finv)
!        dNbardx = dNbardx + J*dNdx*w(kint)*determinant
!        Jbar = Jbar + J*w(kint)*determinant
!        el_vol = el_vol + w(kint)*determinant
!        dNdxvec(1:3*n_nodes) = reshape(transpose(dNdx(1:n_nodes,1:3)),(/3*n_nodes/))
!        do i = 1,n_nodes
!            Pvec = reshape(spread(transpose(dNdx(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
!            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
!        end do
!        P = spread(dNdxvec,dim=1,ncopies=3*n_nodes)*spread(dNdxvec,dim=2,ncopies=3*n_nodes)
!        Pbar = Pbar +  J*(P-Pmat*transpose(Pmat))*w(kint)*determinant
!    end do
!    dNbardx = dNbardx/Jbar
!    dNdxvec(1:3*n_nodes) = reshape(transpose(dNbardx(1:n_nodes,1:3)),(/3*n_nodes/))
!    Pbar = Pbar/Jbar - spread(dNdxvec,dim=1,ncopies=3*n_nodes)*spread(dNdxvec,dim=2,ncopies=3*n_nodes)
!    Jbar = Jbar/el_vol
!
!
!    do kint = 1,n_points
!        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!        call invert_small(dxdxi,dxidx,determinant)
!        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        do i = 1,3
!            ie = 3*(n_nodes-1)+i
!            F(i,1:3) = matmul(dof_increment(i:ie:3)+dof_total(i:ie:3),dNdx(1:n_nodes,1:3))
!            F(i,i) = F(i,i) + 1.d0
!        end do
!        call invert_small(F,Finv,J)
!        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),Finv)
!        Fscale = (Jbar/J)**(1.d0/3.d0)
!        Fbar = F*Fscale
!
!        call user_material_hyperelastic(element_properties,n_properties,Fbar,Jbar,stress,D)
!
!        press = sum(stress(1:3))/3.d0
!
!        Bbar(1:6,1:length_dof_array) = 0.d0
!        Bbar(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
!        Bbar(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
!        Bbar(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
!        Bbar(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
!        Bbar(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
!        Bbar(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
!        Bbar(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
!        Bbar(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
!        Bbar(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
!
!        S = reshape(matmul(transpose(Bbar),stress),(/3,length_dof_array/3/))
!
!        do i = 1,3
!            Bbar(i,1:3*n_nodes-2:3) = Bbar(i,1:3*n_nodes-2:3) +(dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3.d0
!            Bbar(i,2:3*n_nodes-1:3) = Bbar(i,2:3*n_nodes-1:3)+(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3.d0
!            Bbar(i,3:3*n_nodes:3)   = Bbar(i,3:3*n_nodes:3)  +(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3.d0
!        end do
!
!        do i = 1,n_nodes
!            Pvec = reshape(spread(transpose(dNdx(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
!            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
!            Svec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
!            Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
!        end do
!
!        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(Bbar),stress)*w(kint)*determinant
!
!        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
!            + matmul(transpose(Bbar(1:6,1:3*n_nodes)),matmul(D,Bbar(1:6,1:3*n_nodes)))*w(kint)*determinant  &
!            - Pmat*transpose(Smat)*w(kint)*determinant + press*(Pbar + Pmat*transpose(Pmat))*w(kint)*determinant
!
!
!    end do
!
!    return
!end subroutine bbar_hyperelastic




subroutine bbar_hyperelastic_fieldvariables(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                                    ! Input variables
    dof_increment, dof_total,length_dof_array,  &                                                             ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                                     ! Input variables
    n_field_variables,field_variable_names, &                                                                 ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
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

    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

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

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step

    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Element stiffness (ROW,COLUMN)

    real (prec)  ::  stress(6)                                ! Stress vector contains [s11, s22, s33, s12, s13, s23] Kirchhoff stress
    real (prec)  ::  D(6,6)
    real (prec)  ::  dxidx(3,3), determinant                  ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  F(3,3)                                   ! Def gradient
    real (prec)  ::  Finv(3,3)                                ! Def gradient inverse
    real (prec)  ::  Fbar(3,3)                                ! Locking free deformation gradient
    real (prec)  ::  J                                        ! det (F)
    real (prec)  ::  Jbar                                     ! Vol averaged det(F)
    real (prec)  ::  el_vol
    real (prec)  ::  Fscale

    integer :: n_points        ! No. integration points
    integer :: kint            ! Integration point index
    integer :: ie
    integer :: i
    integer :: k

    logical :: strcmp

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0

    !   Volume averaged Jacobian
    Jbar = 0.d0
    el_vol = 0.d0
    do kint = 1,n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        do i = 1,3
            ie = 3*(n_nodes-1)+i
            F(i,1:3) = matmul(dof_increment(i:ie:3)+dof_total(i:ie:3),dNdx(1:n_nodes,1:3))
            F(i,i) = F(i,i) + 1.d0
        end do
        call invert_small(F,Finv,J)
        Jbar = Jbar + J*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
    end do
    Jbar = Jbar/el_vol

!   Find stress and project it
    do kint = 1,n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        do i = 1,3
            ie = 3*(n_nodes-1)+i
            F(i,1:3) = matmul(dof_increment(i:ie:3)+dof_total(i:ie:3),dNdx(1:n_nodes,1:3))
            F(i,i) = F(i,i) + 1.d0
        end do
        call invert_small(F,Finv,J)
        Fscale = (Jbar/J)**(1.d0/3.d0)
        Fbar = F*Fscale

        call user_material_hyperelastic(element_properties,n_properties,Fbar,Jbar,stress,D)
        stress = stress/Jbar             !  Cauchy stress

        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            end if
         end do


    end do



    return

end subroutine bbar_hyperelastic_fieldvariables

!subroutine user_material_hyperelastic(element_properties,n_properties,F,J,stress,D)
!
!   use Types
!   use ParamIO
!
!   implicit none
!
!   integer, intent(in)           :: n_properties
!   real (prec), intent(in)       :: element_properties(n_properties)
!   real (prec), intent(in)       :: F(3,3)
!   real (prec), intent(in)       :: J
!   real (prec), intent(out)      :: stress(6)
!   real (prec), intent(out)      :: D(6,6)
!
!   real (prec) :: B(3,3)
!
!   real (prec) :: mu
!   real (prec) :: K
!   real (prec) :: ss
!   real (prec) :: trB
!   real (prec) :: M1(6,6)
!   real (prec) :: M2(6,6)
!   real (prec) :: M3(6,6)
!
!   integer :: i
!
!   mu = element_properties(1)
!   K  = element_properties(2)
!
!   B = matmul(F,transpose(F))
!
!   ss = J**(-2.d0/3.d0)
!   do i = 1,3
!     stress(i) = mu*B(i,i)*ss
!   end do
!   trB = sum(stress(1:3))/3.d0
!   stress(1:3) = stress(1:3) - trB + K*J*(J-1.d0)
!   stress(4) = mu*B(1,2)*ss
!   stress(5) = mu*B(1,3)*ss
!   stress(6) = mu*B(2,3)*ss
!
!   M1 = 0.d0
!   M1(1,1) = 2.d0*B(1,1)
!   M1(2,2) = 2.d0*B(2,2)
!   M1(3,3) = 2.d0*B(3,3)
!   M1(1,4) = B(2,1)
!   M1(4,1) = B(2,1)
!   M1(1,5) = B(1,3)
!   M1(5,1) = B(1,3)
!   M1(2,4) = B(1,2)
!   M1(4,2) = B(1,2)
!   M1(6,2) = B(2,3)
!   M1(2,6) = B(2,3)
!   M1(5,3) = B(1,3)
!   M1(3,5) = B(1,3)
!   M1(3,6) = B(2,3)
!   M1(6,3) = B(2,3)
!   M1(4,4) = 0.5d0*(B(1,1)+B(2,2))
!   M1(5,5) = 0.5d0*(B(1,1)+B(3,3))
!   M1(6,6) = 0.5d0*(B(2,2)+B(3,3))
!   M1(4,5) = 0.5d0*B(3,2)
!   M1(5,4) = 0.5d0*B(3,2)
!   M1(4,6) = 0.5d0*B(1,3)
!   M1(6,4) = 0.5d0*B(1,3)
!   M1(5,6) = 0.5d0*B(1,2)
!   M1(6,5) = 0.5d0*B(1,2)
!
!   M2 = 0.d0
!   M2(1,1) = 2.d0*B(1,1)
!   M2(2,2) = 2.d0*B(2,2)
!   M2(3,3) = 2.d0*B(3,3)
!   M2(1,2) = B(1,1) + B(2,2)
!   M2(2,1) = M2(1,2)
!   M2(1,3) = B(1,1) + B(3,3)
!   M2(3,1) = M2(1,3)
!   M2(2,3) = B(2,2) + B(3,3)
!   M2(3,2) = M2(2,3)
!
!   M3 = 0.d0
!   M3(1:3,1:3) = 1.d0
!
!   D = mu*ss*M1 - (2.d0*mu*ss/3.d0) * M2 + (2.d0*trB/3.d0 + K*J*(2.d0*J-1.d0))*M3
!
!   return
!
!end subroutine user_material_hyperelastic


subroutine bbar_hyperelastic(lmn, element_identifier, n_nodes, node_property_list, &          ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                     ! Input variables
    dof_increment, dof_total, length_dof_array, &                                             ! Input variables
    n_state_variables, initial_state_variables, &                                             ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables
    use Types
    use ParamIO
    use Globals, only: TIME, DTIME                  ! Total analysis time and time increment
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

    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Length of coord array
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

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


    real (prec)  ::  stress(6)                                ! Stress vector contains [s11, s22, s33, s12, s13, s23] Kirchhoff stress
    real (prec)  ::  D(6,6)                                   ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(3,3)                    ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bbar(6,length_dof_array)                 ! locking free strain = Bbar*(dof_total+dof_increment)
    real (prec)  ::  Bstar(9,length_dof_array)                ! locking free strain = Bbar*(dof_total+dof_increment)
    real (prec)  ::  G(6,9)                                   ! DB/dF * F^T
    real (prec)  ::  dxidx(3,3), determinant                  ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  F(3,3)                                   ! Def gradient
    real (prec)  ::  Fbar(3,3)                                ! Locking free deformation gradient
    real (prec)  ::  Finv(3,3)                                ! Inverse of def gradient
    real (prec)  ::  J                                        ! det (F)
    real (prec)  ::  Jbar                                     ! Vol averaged det(F)
    real (prec)  ::  dNdxvec(3*n_nodes)                ! dN/dx stored as a vector
    real (prec)  ::  Pbar(length_dof_array,length_dof_array)  ! d^2 Jbar/ du_i^a du_k^b  matrix
    real (prec)  ::  Pvec(length_dof_array)
    real (prec)  ::  Pmat(length_dof_array,length_dof_array)
    real (prec)  ::  P(length_dof_array,length_dof_array)     !
    real (prec)  ::  S(3,length_dof_array/3)
    real (prec)  ::  Svec(length_dof_array)
    real (prec)  ::  Smat(length_dof_array,length_dof_array)
    real (prec)  ::  el_vol
    real (prec)  ::  Fscale
    real (prec)  ::  press

    integer :: n_points        ! No. integration points
    integer :: kint            ! Integration point index
    integer :: ie
    integer :: i


    fail = .false.

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0


    !   Volume averaged element variables
    Jbar = 0.d0
    dNbardx = 0.d0
    Pbar = 0.d0
    el_vol = 0.d0
    do kint = 1,n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        do i = 1,3
            ie = 3*(n_nodes-1)+i
            F(i,1:3) = matmul(dof_increment(i:ie:3)+dof_total(i:ie:3),dNdx(1:n_nodes,1:3))
            F(i,i) = F(i,i) + 1.d0
        end do
        call invert_small(F,Finv,J)
        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),Finv)
        dNbardx = dNbardx + J*dNdx*w(kint)*determinant
        Jbar = Jbar + J*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
        dNdxvec(1:3*n_nodes) = reshape(transpose(dNdx(1:n_nodes,1:3)),(/3*n_nodes/))
        do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdx(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
        end do
        P = spread(dNdxvec,dim=1,ncopies=3*n_nodes)*spread(dNdxvec,dim=2,ncopies=3*n_nodes)
        Pbar = Pbar +  J*(P-Pmat*transpose(Pmat))*w(kint)*determinant
    end do
    dNbardx = dNbardx/Jbar
    dNdxvec(1:3*n_nodes) = reshape(transpose(dNbardx(1:n_nodes,1:3)),(/3*n_nodes/))
    Pbar = Pbar/Jbar - spread(dNdxvec,dim=1,ncopies=3*n_nodes)*spread(dNdxvec,dim=2,ncopies=3*n_nodes)
    Jbar = Jbar/el_vol


    do kint = 1,n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        do i = 1,3
            ie = 3*(n_nodes-1)+i
            F(i,1:3) = matmul(dof_increment(i:ie:3)+dof_total(i:ie:3),dNdx(1:n_nodes,1:3))
            F(i,i) = F(i,i) + 1.d0
        end do
        call invert_small(F,Finv,J)
        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),Finv)
        Fscale = (Jbar/J)**(1.d0/3.d0)
        Fbar = F*Fscale

        call user_material_hyperelastic(element_properties,n_properties,Fbar,Jbar,stress,D)

        press = sum(stress(1:3))/3.d0

        Bbar(1:6,1:length_dof_array) = 0.d0
        Bbar(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        Bbar(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        Bbar(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        Bbar(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        Bbar(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        Bbar(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        Bbar(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        Bbar(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        Bbar(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        S = reshape(matmul(transpose(Bbar),stress),(/3,length_dof_array/3/))

        do i = 1,3
            Bbar(i,1:3*n_nodes-2:3) = Bbar(i,1:3*n_nodes-2:3) +(dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3.d0
            Bbar(i,2:3*n_nodes-1:3) = Bbar(i,2:3*n_nodes-1:3)+(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3.d0
            Bbar(i,3:3*n_nodes:3)   = Bbar(i,3:3*n_nodes:3)  +(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3.d0
        end do

        Bstar(1:9,1:length_dof_array) = 0.d0
        Bstar(1:3,1:3*n_nodes) = Bbar(1:3,1:3*n_nodes)
        Bstar(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        Bstar(5,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        Bstar(6,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        Bstar(7,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        Bstar(8,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        Bstar(9,3:3*n_nodes:3) =   dNdx(1:n_nodes,2)

        do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdx(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
            Svec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
        end do

        B = matmul(Fbar,transpose(Fbar))
        G = 0.d0
        G(1,1) = B(1,1)
        G(2,2) = B(2,2)
        G(3,3) = B(3,3)
        G(4,1) = B(1,2)
        G(5,1) = B(1,3)
        G(4,2) = B(1,2)
        G(6,2) = B(2,3)
        G(5,3) = B(1,3)
        G(6,3) = B(2,3)
        G(1,4) = B(1,2)
        G(4,4) = B(2,2)
        G(5,4) = B(2,3)
        G(2,5) = B(1,2)
        G(4,5) = B(1,1)
        G(6,5) = B(1,3)
        G(1,6) = B(1,3)
        G(4,6) = B(2,3)
        G(5,6) = B(3,3)
        G(3,7) = B(1,3)
        G(5,7) = B(1,1)
        G(6,7) = B(1,2)
        G(2,8) = B(2,3)
        G(4,8) = B(1,3)
        G(6,8) = B(3,3)
        G(3,9) = B(1,3)
        G(5,9) = B(1,2)
        G(6,9) = B(2,2)
 !       G(1:3,1:9) = 2.d0*G(1:3,1:9)
        G = 2.d0*G


        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(Bbar),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(Bbar(1:6,1:3*n_nodes)),matmul(D,matmul(G,Bstar(1:9,1:3*n_nodes))))*w(kint)*determinant  &
            - Pmat*transpose(Smat)*w(kint)*determinant + press*(Pbar + Pmat*transpose(Pmat))*w(kint)*determinant



    end do

    return
end subroutine bbar_hyperelastic

!
subroutine user_material_hyperelastic(element_properties,n_properties,F,J,stress,D)

   use Types
   use ParamIO
   use Element_Utilities, only : invert_small

   implicit none

   integer, intent(in)           :: n_properties
   real (prec), intent(in)       :: element_properties(n_properties)
   real (prec), intent(in)       :: F(3,3)
   real (prec), intent(in)       :: J
   real (prec), intent(out)      :: stress(6)
   real (prec), intent(out)      :: D(6,6)

   real (prec) :: B(3,3)
   real (prec) :: Binv(3,3)
   real (prec) :: Bvec(6)
   real (prec) :: Binvvec(6)
   real (prec) :: eyevec(6)
   real (prec) :: mu
   real (prec) :: K
   real (prec) :: ss
   real (prec) :: trB
   real (prec) :: M1(6,6)
   real (prec) :: M2(6,6)
   real (prec) :: M3(6,6)

   integer :: i

   mu = element_properties(1)
   K  = element_properties(2)

   B = matmul(F,transpose(F))
   call invert_small(B,Binv,ss)

   ss = J**(-2.d0/3.d0)
   do i = 1,3
     stress(i) = mu*B(i,i)*ss
   end do
   trB = sum(stress(1:3))/3.d0
   stress(1:3) = stress(1:3) - trB + K*J*(J-1.d0)
   stress(4) = mu*B(1,2)*ss
   stress(5) = mu*B(1,3)*ss
   stress(6) = mu*B(2,3)*ss
   D = 0.d0
   D(1,1) = 1.d0
   D(2,2) = 1.d0
   D(3,3) = 1.d0
   D(4,4) = 0.5d0
   D(5,5) = 0.5d0
   D(6,6) = 0.5d0
   D = D*mu*ss

   eyevec(1:3) = 1.d0
   eyevec(4:6) = 0.d0
   Bvec(1) = B(1,1)
   Bvec(2) = B(2,2)
   Bvec(3) = B(3,3)
   Bvec(4) = B(1,2)
   Bvec(5) = B(1,3)
   Bvec(6) = B(2,3)
   Binvvec(1) = Binv(1,1)
   Binvvec(2) = Binv(2,2)
   Binvvec(3) = Binv(3,3)
   Binvvec(4) = Binv(1,2)
   Binvvec(5) = Binv(1,3)
   Binvvec(6) = Binv(2,3)

   trB = sum(Bvec(1:3))/3.d0

   D = D + (ss*mu/3.d0)*( trB*spread(eyevec,dim=2,ncopies=6)*spread(Binvvec,dim=1,ncopies=6) &
                        - spread(eyevec,dim=2,ncopies=6)*spread(eyevec,dim=1,ncopies=6)  &
                        - spread(Bvec,dim=2,ncopies=6)*spread(Binvvec,dim=1,ncopies=6) )

   D = D + K*J*(J-0.5d0)*spread(eyevec,dim=2,ncopies=6)*spread(Binvvec,dim=1,ncopies=6)


   return

end subroutine user_material_hyperelastic
