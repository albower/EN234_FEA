!
!     Subroutines to assemble element level matrices for FEA analysis


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine gurson_dynamic(lmn, element_identifier, n_nodes, node_property_list, &               ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                              ! Input variables
    n_state_variables, initial_state_variables, &                                               ! Input variables
    updated_state_variables,element_residual,element_deleted)                                   ! Output variables
    use Types
    use Globals, only : TIME, DTIME
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
    use Element_Utilities, only : rotatesymvec
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
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

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine

    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    logical, intent( inout )  :: element_deleted                                             ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint
    integer      :: n_state_vars_per_intpt, iof


    real (prec)  ::  stress0(6)                               ! Initial stress vector contains [s11, s22, s33, s12, s13, s23] Kirchhoff stress
    real (prec)  ::  stress1(6)                               ! Stress vector at end of increment
    real (prec)  ::  Bbar(6,length_dof_array)                 ! locking free strain = Bbar*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant                  ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  dF(3,3)                                  ! Def gradient increment
    real (prec)  ::  Fmid(3,3)                                ! Mid point def gradient
    real (prec)  ::  Fmidinv(3,3)                             ! Inverse of mid point def gradient
    real (prec)  ::  Jmid
    real (prec)  ::  Jmidbar                                  ! Vol averaged det(Fmid)
    real (prec)  ::  dLbar(3,3)                               ! Vol averaged increment in dL_ij
    real (prec)  ::  dLkkbar                                  ! Trace of dLbar
    real (prec)  ::  deps(3,3)                                ! Strain increment
    real (prec)  ::  depskk                                   ! Trace of strain increment
    real (prec)  ::  dW(3,3)                                  ! Spin increment
    real (prec)  ::  dR(3,3)                                  ! Rotation increment
    real (prec)  ::  el_vol
    real (prec)  ::  temp33a(3,3)                             ! Workspace array
    real (prec)  ::  temp33b(3,3)                             ! Workspace array
    real (prec)  ::  temp33c(3,3)                             ! Workspace array
    real (prec)  ::  dummy                                    ! Dummy variable
    real (prec)  ::  Vf                                       ! Void volume fraction (used for element deletion)

    integer :: ie
    integer :: i


    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    n_state_vars_per_intpt = n_state_variables/n_points

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0

    !   Volume averaged element variables
    Jmidbar = 0.d0
    dLbar = 0.d0
    dNbardx = 0.d0
    el_vol = 0.d0
    do kint = 1,n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        do i = 1,3
            ie = 3*(n_nodes-1)+i
            Fmid(i,1:3) = matmul(dof_total(i:ie:3)+0.5d0*dof_increment(i:ie:3),dNdx(1:n_nodes,1:3))
            Fmid(i,i) = Fmid(i,i) + 1.d0
            dF(i,1:3) = matmul(dof_increment(i:ie:3),dNdx(1:n_nodes,1:3))
        end do
        call invert_small(Fmid,Fmidinv,Jmid)
        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),Fmidinv)
        dNbardx = dNbardx + Jmid*dNdx*w(kint)*determinant
        dLbar = dLbar + Jmid*matmul(dF,Fmidinv)*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
        Jmidbar = Jmidbar + Jmid*w(kint)*determinant
    end do
    dNbardx = dNbardx/Jmidbar
    dLkkbar = (dLbar(1,1) + dLbar(2,2) + dLbar(3,3))/Jmidbar
    element_deleted = .true.
    do kint = 1,n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        do i = 1,3
            ie = 3*(n_nodes-1)+i
            Fmid(i,1:3) = matmul(dof_total(i:ie:3)+0.5d0*dof_increment(i:ie:3),dNdx(1:n_nodes,1:3))
            Fmid(i,i) = Fmid(i,i) + 1.d0
            dF(i,1:3) = matmul(dof_increment(i:ie:3),dNdx(1:n_nodes,1:3))
        end do
        call invert_small(Fmid,Fmidinv,Jmid)
        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),Fmidinv)

        deps = matmul(dF,Fmidinv)
        dW = 0.5d0*(deps - transpose(deps))
        deps = 0.5d0*(deps + transpose(deps))
        depskk = deps(1,1) +deps(2,2) + deps(3,3)
        deps(1,1) = deps(1,1) + (dLkkbar - depskk)/3.d0
        deps(2,2) = deps(2,2) + (dLkkbar - depskk)/3.d0
        deps(3,3) = deps(3,3) + (dLkkbar - depskk)/3.d0

        temp33a = eye3_D + 0.5d0*dW
        temp33b = eye3_D - 0.5d0*dW
        call invert_small(temp33b,temp33c,dummy)
        dR = matmul(temp33c,temp33a)

        iof = n_state_vars_per_intpt*(kint-1) + 1
        stress0 = initial_state_variables(iof:iof+5)
        stress1 = rotatesymvec(stress0,dR)
        stress0 = stress1

        call gurson(element_properties,n_properties,initial_state_variables(iof), &
                       updated_state_variables(iof),n_state_vars_per_intpt, &
                                                  deps,stress0,stress1)

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

        do i = 1,3
            Bbar(i,1:3*n_nodes-2:3) = Bbar(i,1:3*n_nodes-2:3) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3.d0
            Bbar(i,2:3*n_nodes-1:3) = Bbar(i,2:3*n_nodes-1:3) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3.d0
            Bbar(i,3:3*n_nodes:3)   = Bbar(i,3:3*n_nodes:3)   + (dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3.d0
        end do

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(Bbar),stress1)*w(kint)*determinant

        Vf = updated_state_variables(iof+7)
        if (Vf< element_properties(13)) element_deleted = .false.


    end do

    return
end subroutine gurson_dynamic


subroutine dynamic_gurson_fieldvariables(lmn, element_identifier, n_nodes, node_property_list, &        ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                       ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
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
    integer, intent( in )         :: length_dof_array                                       !
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

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

    real (prec)  ::  stress(6)                                ! Stress vector contains [s11, s22, s33, s12, s13, s23] Kirchhoff stress
    real (prec)  ::  dxidx(3,3), determinant                  ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  Vf                                       ! Void volume fraction


    integer :: n_points        ! No. integration points
    integer :: kint            ! Integration point index
    integer :: n_state_vars_per_intpt
    integer :: iof

    integer :: k

    logical :: strcmp

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    n_state_vars_per_intpt = n_state_variables/n_points

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0



!   Find stress and project it
    do kint = 1,n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        iof = n_state_vars_per_intpt*(kint-1)
        stress = updated_state_variables(iof+1:iof+6)
        Vf = updated_state_variables(iof+8)
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
            else if (strcmp(field_variable_names(k),'Vf',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + Vf*N(1:n_nodes)*determinant*w(kint)
            end if
         end do


    end do

    return

end subroutine dynamic_gurson_fieldvariables


subroutine gurson(element_properties,n_properties,initial_state_variables,updated_state_variables,n_state_variables, &
                                                                                                      deps,stress0,stress1)
    use Types
    use ParamIO
    use Globals, only : TIME, DTIME

    implicit none

    integer, intent( in )       :: n_properties
    integer, intent( in )       :: n_state_variables

    real (prec), intent( in )   :: element_properties(n_properties)
    real (prec), intent( in )   :: initial_state_variables(n_state_variables)
    real (prec), intent( in )   :: deps(3,3)
    real (prec), intent( in )   :: stress0(6)

    real (prec), intent( out )  :: stress1(6)
    real (prec), intent( out )  :: updated_state_variables(n_state_variables)

    real (prec) ::  dedev(6)            ! Deviatoric strain increment
    real (prec) ::  devol               ! Volumetric strain increment
    real (prec) ::  sdev(6)             ! Deviatoric stress
    real (prec) ::  sdevstar(6)         ! Elastic predictor for deviatoric stress          !
    real (prec) ::  pn

    real (prec) :: E,xnu,Y,edot0,m,q1,q2,q3,fn,sn,en,fc,ff
    real (prec) :: ffbar
    real (prec) :: Vf, ematrix, dematrix
    real (prec) :: fstar

    real (prec) :: sestar, pstar
    real (prec) :: se, p, dsedee, dpdev
    real (prec) :: dee, dev
    real (prec) :: phi, dphidse, dphidp, d2phidse2, d2phidp2, d2phidpdse
    real (prec) :: dphidee, dphidev
    real (prec) :: c1, dc1dse, dc1dp, dc1dee, dc1dev

    real (prec) :: f1, f2, df1dee, df1dev, df2dee, df2dev, det
    real (prec) :: error, tol

    integer :: nit, maxit


    E     = element_properties(1)
    xnu   = element_properties(2)
    Y     = element_properties(3)
    edot0 = element_properties(4)
    m     = element_properties(5)
    q1    = element_properties(6)
    q2    = element_properties(7)
    q3    = element_properties(8)
    fn    = element_properties(9)
    sn    = element_properties(10)
    en    = element_properties(11)
    fc    = element_properties(12)
    ff    = element_properties(13)
    ffbar = (q1 + dsqrt(q1*q1-q3))/q3

    Vf = initial_state_variables(8)
    ematrix = initial_state_variables(7)
    dee = 0.d0 !initial_state_variables(9)
    dev = 0.d0 !initial_state_variables(10)

    devol = deps(1,1) + deps(2,2) + deps(3,3)
    dedev(1) = deps(1,1)
    dedev(2) = deps(2,2)
    dedev(3) = deps(3,3)
    dedev(4) = deps(1,2)
    dedev(5) = deps(1,3)
    dedev(6) = deps(2,3)
    dedev(1:3) = dedev(1:3) - devol/3.d0

    pn = sum(stress0(1:3))/3.d0
    pstar = pn + E*devol/(1.d0-2.d0*xnu)/3.d0
    sdev = stress0
    sdev(1:3) = sdev(1:3) - pn
    sdevstar = sdev + E/(1.d0+xnu)*dedev

    sestar = dsqrt(1.5d0*( dot_product(sdevstar(1:3),sdevstar(1:3)) + 2.d0*dot_product(sdevstar(4:6),sdevstar(4:6)) ) )

    if (sestar + pstar*pstar<1.d-08*Y) then
        stress1 = 0.d0
        updated_state_variables = initial_state_variables
        return
    endif

    fstar = Vf
    if (fstar>fc) then
       if (Vf>ff) then
          fstar = ffbar
       else
          fstar = fc + (ffbar-fc)*(Vf-fc)/(ff-fc)
       endif
    endif

    nit = 0
    maxit = 30
    error = 1.d0
    tol = 1.d-08/Y
    dsedee = -1.5d0*E/(1.d0+xnu)
    dpdev = -E/(1.d0-2.d0*xnu)/3.d0
    do while (error>tol)
      se = sestar + dee*dsedee
      p =  pstar + dev*dpdev
      phi = se*se/(Y*Y) + 2.d0*q1*fstar*cosh(1.5d0*q2*p/Y) - (1.d0+q3*fstar*fstar)
      if (phi<1.d-08) then
         dee = 0.d0
         dev = 0.d0
         se = sestar
         p = pstar
         phi = se*se/(Y*Y) + 2.d0*q1*fstar*cosh(1.5d0*q2*p/Y) - (1.d0+q3*fstar*fstar)
         if (phi<1.d-08) exit
      endif
      phi = sqrt(phi)
      dphidse = se/(phi*Y*Y)
      dphidp = 1.5d0*q1*q2*fstar*sinh(1.5d0*q2*p/Y)/(phi*Y)
      d2phidse2 = 1.d0/(phi*Y*Y) - dphidse*se/(phi*phi*Y*Y)
      d2phidp2 = 2.25d0*q1*q2*q2*fstar*cosh(1.5d0*q2*p/Y)/(phi*Y*Y) &
               - 1.5d0*q1*q2*fstar*sinh(1.5d0*q2*p/Y)*dphidp/(phi*phi*Y)
      d2phidpdse = -se*dphidp/(phi*phi*Y*Y)

      dphidee = dphidse*dsedee
      dphidev = dphidp*dpdev

      c1 = dsqrt(dphidse*dphidse + dphidp*dphidp/4.5d0)
      dc1dse = (dphidse*d2phidse2 + dphidp*d2phidpdse/4.5d0)/c1
      dc1dp = (dphidp*d2phidp2/4.5d0 + dphidse*d2phidpdse)/c1

      dc1dee = dc1dse*dsedee
      dc1dev = dc1dp*dpdev

      f1 = ( c1*dee/(DTIME*edot0) )  - dphidse*phi**m
      f2 = ( c1*dev/(DTIME*edot0) )  - dphidp*phi**m
      df1dee = (dc1dee*dee/(DTIME*edot0) + c1/(DTIME*edot0)) &
             - m*phi**(m-1.d0)*dphidee*dphidse - d2phidse2*dsedee*phi**m
      df1dev = (dc1dev*dee/(DTIME*edot0)) &
            - m*phi**(m-1.d0)*dphidev*dphidse - d2phidpdse*dpdev*phi**m
      df2dee = dc1dee*dev/(DTIME*edot0) &
              - m*phi**(m-1.d0)*dphidee*dphidp - d2phidpdse*dsedee*phi**m
      df2dev = ( dc1dev*dev/(DTIME*edot0) + c1/(DTIME*edot0) ) &
             - m*phi**(m-1.d0)*dphidev*dphidp - d2phidp2*dpdev*phi**m

      det = df1dee*df2dev - df1dev*df2dee

      dee = dee - (f1*df2dev - f2*df1dev)/det
      dev = dev - (f2*df1dee - f1*df2dee)/det

      error = f1*f1 + f2*f2
      nit = nit + 1
      if (nit>maxit) then
         write(IOW,*) ' Newton iterations in gurson failed to converge '
         write(IOW,*) ' Analysis terminated '
         stop
      endif
    end do

    stress1 = sdevstar + dee*dsedee*sdevstar/sestar
    p = pstar + dev*dpdev
    stress1(1:3) = stress1(1:3) + p


    dematrix = 0.d0
    if (phi>1.d-08) dematrix = edot0*DTIME/(1.d0-Vf) * (phi**m/c1) * (dphidse*se + dphidp*p/3.d0)
    Vf = 1.d0 + (Vf-1.d0)*dexp(-dev) + dematrix*fn/(sn*dsqrt(TWOPI_D)) * dexp( -0.5d0*(ematrix-en)*(ematrix-en)/(sn*sn) )
    ematrix = ematrix + dematrix

    updated_state_variables(1:6) = stress1
    updated_state_variables(7) = ematrix
    updated_state_variables(8) = Vf
    updated_state_variables(9) = dee
    updated_state_variables(10) = dev

end subroutine gurson

