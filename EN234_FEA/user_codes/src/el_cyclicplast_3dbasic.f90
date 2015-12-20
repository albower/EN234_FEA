!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_cyclicplast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
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

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint
    integer      :: n_state_vars_per_intpt, iof

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    n_state_vars_per_intpt = n_state_variables/n_points
  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      
        iof = (kint-1)*n_state_vars_per_intpt + 1

        call nl_kinematic_hardening(n_properties, element_properties, n_state_variables, initial_state_variables(iof), &
           strain, dstrain, stress, D,  updated_state_variables(iof), fail)

        if (fail) return

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

    end do
  

    return
end subroutine el_cyclicplast_3dbasic




!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_cyclicplast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                  !
    dof_increment, dof_total,length_dof_array,  &                                                          ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                                  ! Input variables
    n_field_variables,field_variable_names, &                                                              ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
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
    real( prec ), intent( in )   :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k
    integer      :: n_state_vars_per_intpt, iof

    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: p, ep                          ! Pressure and Mises stress
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     EPROP(1)         Young's modulus
    !     EPROP(2)         Poisson's ratio
    !     EPROP(3)         Thermal expansion coefficient


    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	n_state_vars_per_intpt = n_state_variables/n_points

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)

        iof = n_state_vars_per_intpt*(kint-1)
        sdev(1:6) = updated_state_variables(iof+1:iof+6)
        p = updated_state_variables(iof+7)
        stress = sdev
        stress(1:3) = stress(1:3) + p

        ep = updated_state_variables(iof+14)

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
            else if (strcmp(field_variable_names(k),'ep',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + ep*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_cyclicplast_3dbasic

subroutine nl_kinematic_hardening(n_properties, element_properties, n_state_variables, initial_state_variables, &
           strain, dstrain, stress, D,  updated_state_variables, fail)
    use Types
    use ParamIO

    implicit none

    integer, intent( in )         :: n_properties
    integer, intent( in )         :: n_state_variables

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine

    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step

    real (prec), intent( in )     ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec), intent( out )   ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec), intent( out )    :: D(6,6)

    real (prec) :: E,xnu,Y0,c,gam                                      ! Material properties
    real (prec) :: d11,d12,d44
    real (prec) :: tol                                                 ! Tolerance for Newton iterations
    real (prec) :: err                                                 ! Error in Newton iterations

    real (prec) :: sdev0(6)
    real (prec) :: press0
    real (prec) :: alpha0(6)
    real (prec) :: devol
    real (prec) :: dedev(6)
    real (prec) :: sdevstar(6)                                         ! Deviatoric stress predictor
    real (prec) :: sdev(6)                                             ! Updated deviatoric stress
    real (prec) :: press
    real (prec) :: sequiv
    real (prec) :: dep
    real (prec) :: c1,f1,df1dep,f,dfdep
    real (prec) :: eta,lam,detadep,dlamdep,beta

    integer :: nit, maxit

    E = element_properties(1)
    xnu = element_properties(2)
    Y0 = element_properties(3)
    c = element_properties(4)
    gam = element_properties(5)

    tol = 1.d-9*Y0
    maxit = 30

    sdev0 = initial_state_variables(1:6)
    press0 = initial_state_variables(7)
    alpha0 =  initial_state_variables(8:13)

    devol = sum(dstrain(1:3))
    dedev(1:3) = dstrain(1:3) - devol/3.d0
    dedev(4:6) = 0.5d0*dstrain(4:6)

    sdevstar = sdev0 + E*dedev/(1.d0+xnu)
    press = press0 + E*devol/(3.d0*(1.d0-2.d0*xnu))
    sequiv = dsqrt(1.5d0*( dot_product(sdevstar(1:3)-alpha0(1:3),sdevstar(1:3)-alpha0(1:3)) &
                    + 2.d0*dot_product(sdevstar(4:6)-alpha0(4:6),sdevstar(4:6)-alpha0(4:6)) ) )

    if (sequiv<Y0) then
       stress = sdevstar
       stress(1:3) = stress(1:3) + press
       updated_state_variables(1:6) = sdevstar(1:6)
       updated_state_variables(7) = press
       updated_state_variables(8:14) = initial_state_variables(8:14)
       updated_state_variables(15) = 0.d0
       d44 = 0.5D0*E/(1+xnu)
       d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
       d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
       D = 0.d0
       D(1:3,1:3) = d12
       D(1,1) = d11
       D(2,2) = d11
       D(3,3) = d11
       D(4,4) = d44
       D(5,5) = d44
       D(6,6) = d44
       return
    endif

    dep = 0.d0
    err = 1.d0
    nit = 0
    do while (err>tol)
      c1 = 1.d0/(1.d0+gam*dep)
      f1 = dsqrt(1.5d0*(   dot_product(sdevstar(1:3)-c1*alpha0(1:3),sdevstar(1:3)-c1*alpha0(1:3)) &
                    + 2.d0*dot_product(sdevstar(4:6)-c1*alpha0(4:6),sdevstar(4:6)-c1*alpha0(4:6)) ) )
      df1dep = 1.5d0*gam*c1*c1*c1*( dot_product(sdevstar(1:3)-c1*alpha0(1:3),alpha0(1:3))  &
                           +2.d0*dot_product(sdevstar(4:6)-c1*alpha0(4:6),alpha0(4:6)) )/f1
      f = Y0 + (1.5d0*E/(1.d0+xnu) + c*c1)*dep - f1
      dfdep = 1.5d0*E/(1.d0+xnu) + c*c1 - c*gam*dep*c1*c1 - df1dep
      dep = dep - f/dfdep
      err = f*f
      nit = nit + 1
      if (nit>maxit) then
         fail = .true.
         write(6,*) ' Newton iterations in nl_kinematic_hardening failed to converge '
         write(6,*) ' Forcing cutback '
         return
      endif
    end do

    eta = 1.5d0*E*dep/(   (1.d0+xnu)*( Y0 + (gam*Y0+c)*dep )   )
    lam = 1.d0/(   1.d0  + (1.5d0*E/(1.d0+xnu) - eta*c)*dep/Y0    )

    detadep = 1.5d0*E/( (1.d0+xnu)*(Y0 + (gam*Y0+c)*dep) ) - eta*(gam*Y0+c)/(Y0 + (gam*Y0+c)*dep)
    dlamdep = -lam*lam*( ( 1.5d0*E/(1.d0+xnu) - eta*c )/Y0  - detadep*c*dep/Y0 )
    beta = 1.d0/(Y0 + (1.5d0*E/(1.d0+xnu) + c*c1)*dep)

    sdev = lam*(sdevstar + eta*alpha0)
    stress = sdev
    stress(1:3) = stress(1:3) + press

    d44 = 0.5D0*E*lam/(1.d0+xnu)
    d11 = E*lam/(1.d0 + xnu)/1.5d0 + E/( 3.d0*(1.d0 - 2.d0*xnu) )
    d12 = -E*lam/(1.d0 + xnu)/3.d0 + E/( 3.d0*(1.d0 - 2.d0*xnu) )
    D = 0.d0
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

    D = D + 1.5d0*E/(1.d0+xnu)*beta/dfdep*( &
        dlamdep*spread(sdevstar,dim=2,ncopies=6)*spread(sdevstar,dim=1,ncopies=6) &
        - (dlamdep*eta + lam*detadep)*c1*spread(alpha0,dim=2,ncopies=6)*spread(alpha0,dim=1,ncopies=6) &
        + (dlamdep*eta + lam*detadep)*spread(alpha0,dim=2,ncopies=6)*spread(sdevstar,dim=1,ncopies=6) &
        -  dlamdep*c1*spread(sdevstar,dim=2,ncopies=6)*spread(alpha0,dim=1,ncopies=6) )

    updated_state_variables(1:6) = sdev
    updated_state_variables(7) = press
    updated_state_variables(8:13) = (alpha0 + (c*dep/Y0)*sdev )/(1.d0 + (gam+c/Y0)*dep)
    updated_state_variables(14) = initial_state_variables(14) + dep
    updated_state_variables(15) = dep



end subroutine nl_kinematic_hardening

