!
!     Subroutines to assemble element level matrices for FEA analysis

!==========================SUBROUTINE continuum_element_static_3D ==============================
subroutine continuum_element_static_2D(current_step_number,lmn, element_identifier, &
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
    real (prec)  ::  stress1(3,3)                             ! Stress vector at end of increment
    real (prec)  ::  press                                    ! Hydrostatic stress
    real (prec)  ::  Bbar(4,length_dof_array)                 ! locking free strain = Bbar*(dof_total+dof_increment)
    real (prec)  ::  Bstar(5,length_dof_array)                ! Locking free F = Bstar*(dof_total+dof_increment)
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
    real (prec)  ::  dW(3,3)                                  ! Spin increment components in global basis
    real (prec)  ::  dR(3,3)                                  ! Rotation increment components in global basis
    real (prec)  ::  strain0(3,3)                             ! Strain at start of increment; updated for rigid rotation
    real (prec)  ::  el_vol
    real (prec)  ::  temp33a(3,3)                             ! Workspace array
    real (prec)  ::  temp33b(3,3)                             ! Workspace array
    real (prec)  ::  temp33c(3,3)                             ! Workspace array
    real (prec)  ::  depsdf(4,5)                              ! Derivative of strain increment with respect to def gradient
    real (prec)  ::  dWdf1(4,5)                               ! part of Derivative of spin stress increment with respect to def gradient
    real (prec)  ::  dWdf2(4,5)                               ! part of derivative of spin stress increment with respect to def grad

    real (prec)  ::  H(3,3),Lam(3,3)                          ! Workspace arrays used to compute stiffness
    real (prec)  ::  SH(3,3),SLam(3,3)                        ! Workspace arrays used to compute stiffness
    real (prec)  ::  dummy                                    ! Dummy variable
    real (prec)  ::  Dcorot(4,4)                              ! Corotational material tangent
    real (prec)  ::  Dspin(4,5)                               ! Spin tangent
    real (prec)  ::  Dtherm(4)                                ! Thermal tangent
!
    real (prec)  ::  dNdxvec(2*n_nodes)                       ! dN/dx stored as a vector
    real (prec)  ::  Pbar(length_dof_array,length_dof_array)  ! d^2 Jbar/ du_i^a du_k^b  matrix
    real (prec)  ::  Pvec(length_dof_array)
    real (prec)  ::  Pmat(length_dof_array,length_dof_array)
    real (prec)  ::  P(length_dof_array,length_dof_array)     !
    real (prec)  ::  S(2,length_dof_array/2)
    real (prec)  ::  Svec(length_dof_array)
    real (prec)  ::  Smat(length_dof_array,length_dof_array)

    real (prec)  ::  stress(4)                                ! Stress (in global basis) stored as vector
    real (prec)  ::  strain(4)                                ! Strain (global basis) stored as vector
    real (prec)  ::  strainInc(4)                             ! Strain increment (in global basis) stored as vector
    real (prec)  ::  SSE,  SPD, SCD                           ! Energy dissipation variables
    real (prec)  ::  charLength                               ! Characteristic element length
    real (prec)  ::  rpl                                      ! Rate of plastic heat generation

    real (prec)  ::  drplde(4)                                ! Derivative of plastic work wrt strain (unused)
    real (prec)  ::  drpldt
    real (prec)  ::  abq_TIME(2)
    real (prec)  ::  predef(1)                                ! Dummy argument for predefined field
    real (prec)  ::  dpred(1)                                 ! Dummy argument for predefined field increment
    real (prec)  ::  dummy_state(1)                           ! Dummy state variable array

    integer :: ie
    integer :: i
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

    if (n_nodes == 3) n_points = 1
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
    dNbardx = 0.d0
    Pbar = 0.d0
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
        dNbardx = dNbardx + J1*dNdx*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
        J0bar = J0bar + det33(F0)*w(kint)*determinant
        J1bar = J1bar + J1*w(kint)*determinant
        dNdxvec(1:2*n_nodes) = reshape(transpose(dNdx(1:n_nodes,1:2)),(/2*n_nodes/))
        do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdx(i:i,1:2)),dim=2,ncopies=n_nodes),(/2*n_nodes/))
            Pmat(2*i-1:2*i,1:2*n_nodes) = spread(Pvec,dim=1,ncopies=2)
        end do
        P = spread(dNdxvec,dim=1,ncopies=2*n_nodes)*spread(dNdxvec,dim=2,ncopies=2*n_nodes)
        Pbar = Pbar +  J1*(P-Pmat*transpose(Pmat))*w(kint)*determinant
    end do
    J0bar = J0bar/el_vol
    J1bar = J1bar/el_vol
    charLength = el_vol**(1.d0/2.d0)
    dNbardx = dNbardx/(J1bar*el_vol)
    dNdxvec(1:2*n_nodes) = reshape(transpose(dNbardx(1:n_nodes,1:2)),(/2*n_nodes/))
    Pbar = Pbar/(J1bar*el_vol) - spread(dNdxvec,dim=1,ncopies=2*n_nodes)*spread(dNdxvec,dim=2,ncopies=2*n_nodes)

    do kint = 1,n_points
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

        if (n_state_vars_per_intpt>11)   updated_state_variables(iof+11:iof+n_state_vars_per_intpt-1) = &
                                         initial_state_variables(iof+11:iof+n_state_vars_per_intpt-1)
        SSE = initial_state_variables(iof+8)
        SPD = initial_state_variables(iof+9)
        SCD = initial_state_variables(iof+10)

        abq_time(1:2) = TIME

        if (n_state_vars_per_intpt>11) then
           Call UMAT(stress,updated_state_variables(iof+11:iof+n_state_vars_per_intpt-1),Dcorot,SSE,SPD,SCD, &
                  RPL,Dtherm,DRPLDE,DRPLDT, &
                  strain,strainInc,abq_time,DTIME,BTEMP,BTINC,predef,dpred,material_name, &
                  3,1,4,n_state_vars_per_intpt-11,material_properties,n_properties,xintpt,dR,PNEWDT, &
                  charLength,F0,F1,lmn,kint,0,0,1,current_step_number)
        else
           Call UMAT(stress,dummy_state,Dcorot,SSE,SPD,SCD, &
                  RPL,Dtherm,DRPLDE,DRPLDT, &
                  strain,strainInc,abq_time,DTIME,BTEMP,BTINC,predef,dpred,material_name, &
                  3,1,4,n_state_vars_per_intpt-11,material_properties,n_properties,xintpt,dR,PNEWDT, &
                  charLength,F0,F1,lmn,kint,0,0,1,current_step_number)
        endif
!       Derivative of strain with respect to def gradient, multiplied by F1bar^T
        H = transpose(matmul(F1,Fmidinv))
        Lam = eye3_D - 0.5d0*matmul(dF,Fmidinv)
        depsdf(1,1:5) = [Lam(1,1)*H(1,1),Lam(1,2)*H(1,2),Lam(1,3)*H(1,3), &
                         Lam(1,1)*H(1,2),Lam(1,2)*H(1,1)]
        depsdf(2,1:5) = [Lam(2,1)*H(2,1),Lam(2,2)*H(2,2),Lam(2,3)*H(2,3), &
                         Lam(2,1)*H(2,2),Lam(2,2)*H(2,1)]
        depsdf(3,1:5) = [Lam(3,1)*H(3,1),Lam(3,2)*H(3,2),Lam(3,3)*H(3,3), &
                         Lam(3,1)*H(3,2),Lam(3,2)*H(3,1)]
        depsdf(4,1:5) = [Lam(1,1)*H(2,1)+Lam(2,1)*H(1,1),Lam(1,2)*H(2,2)+Lam(2,2)*H(1,2),Lam(1,3)*H(2,3)+Lam(2,3)*H(1,3), &
                         Lam(1,1)*H(2,2)+Lam(2,1)*H(1,2),Lam(1,2)*H(2,1)+Lam(2,2)*H(1,1)]


        stress1 = reshape([stress(1),stress(4),0.d0, &
                           stress(4),stress(2),0.d0, &
                           0.d0,0.d0,stress(3)], &
                           shape(stress1))

        SH = matmul(stress1,H)
        SLam = matmul(stress1,Lam)
        dWdf1(1,1:5) = [Lam(1,1)*SH(1,1),Lam(1,2)*SH(1,2),Lam(1,3)*SH(1,3), &
                        Lam(1,1)*SH(1,2),Lam(1,2)*SH(1,1)]
        dWdf1(2,1:5) = [Lam(2,1)*SH(2,1),Lam(2,2)*SH(2,2),Lam(2,3)*SH(2,3), &
                        Lam(2,1)*SH(2,2),Lam(2,2)*SH(2,1)]
        dWdf1(3,1:5) = [Lam(3,1)*SH(3,1),Lam(3,2)*SH(3,2),Lam(3,3)*SH(3,3), &
                        Lam(3,1)*SH(3,2),Lam(3,2)*SH(3,1)]
        dWdf1(4,1:5) = [Lam(1,1)*SH(2,1)+Lam(2,1)*SH(1,1),Lam(1,2)*SH(2,2)+Lam(2,2)*SH(1,2),Lam(1,3)*SH(2,3)+Lam(2,3)*SH(1,3), &
                        Lam(1,1)*SH(2,2)+Lam(2,1)*SH(1,2),Lam(1,2)*SH(2,1)+Lam(2,2)*SH(1,1)]

        dWdf2(1,1:5) = [SLam(1,1)*H(1,1),SLam(1,2)*H(1,2),SLam(1,3)*H(1,3), &
                        SLam(1,1)*H(1,2),SLam(1,2)*H(1,1)]
        dWdf2(2,1:5) = [SLam(2,1)*H(2,1),SLam(2,2)*H(2,2),SLam(2,3)*H(2,3), &
                        SLam(2,1)*H(2,2),SLam(2,2)*H(2,1)]
        dWdf2(3,1:5) = [SLam(3,1)*H(3,1),SLam(3,2)*H(3,2),SLam(3,3)*H(3,3), &
                        SLam(3,1)*H(3,2),SLam(3,2)*H(3,1)]
        dWdf2(4,1:5) = [SLam(1,1)*H(2,1)+SLam(2,1)*H(1,1),SLam(1,2)*H(2,2)+SLam(2,2)*H(1,2),SLam(1,3)*H(2,3)+SLam(2,3)*H(1,3), &
                        SLam(1,1)*H(2,2)+SLam(2,1)*H(1,2),SLam(1,2)*H(2,1)+SLam(2,2)*H(1,1)]

        Dspin = 0.5d0*(dWdf1-dWdf2)

        updated_state_variables(iof:iof+3) = stress(1:4)
        updated_state_variables(iof+4:iof+7) = strain(1:4)+strainInc(1:4)
        updated_state_variables(iof+8) = SSE
        updated_state_variables(iof+9) = SPD
        updated_state_variables(iof+10) = SCD

        Bbar(1:4,1:length_dof_array) = 0.d0
        Bbar(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        Bbar(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
        Bbar(4,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        Bbar(4,2:2*n_nodes:2) = dNdx(1:n_nodes,1)

        S = reshape(matmul(transpose(Bbar),stress),(/2,length_dof_array/2/))

        do i = 1,2
            Bbar(i,1:2*n_nodes-1:2) = Bbar(i,1:2*n_nodes-1:2) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/2.d0
            Bbar(i,2:2*n_nodes:2) = Bbar(i,2:2*n_nodes:2) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/2.d0
        end do


        Bstar(1:5,1:length_dof_array) = 0.d0
        Bstar(1:3,1:2*n_nodes) = Bbar(1:3,1:2*n_nodes)
        Bstar(4,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        Bstar(5,2:2*n_nodes:2) = dNdx(1:n_nodes,1)

        do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdx(i:i,1:2)),dim=2,ncopies=n_nodes),(/2*n_nodes/))
            Pmat(2*i-1:2*i,1:2*n_nodes) = spread(Pvec,dim=1,ncopies=2)
            Svec = reshape(spread(S(1:2,i:i),dim=2,ncopies=n_nodes),(/2*n_nodes/))
            Smat(2*i-1:2*i,1:2*n_nodes) = spread(Svec,dim=1,ncopies=2)
        end do

        press = sum(stress(1:2))/2.d0

        element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(Bbar),stress)*w(kint)*determinant*J1bar

        element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
                 + matmul(transpose(Bbar(1:4,1:2*n_nodes)),matmul(Dcorot,matmul(depsdf(1:4,1:5),Bstar(1:5,1:2*n_nodes)))) &
                                                                                                      *w(kint)*determinant*J1bar &
                 + matmul(transpose(Bbar(1:4,1:2*n_nodes)),matmul(Dspin,Bstar(1:5,1:2*n_nodes)))*w(kint)*determinant*J1bar &
                 - Pmat*transpose(Smat)*w(kint)*determinant*J1bar + press*(Pbar + Pmat*transpose(Pmat))*w(kint)*determinant*J1bar

    end do

    return
end subroutine continuum_element_static_2D


!
!==========================SUBROUTINE continuum_element_static_3D ==============================
subroutine continuum_element_static_3D(current_step_number,lmn, element_identifier, &
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
    real (prec)  ::  stress1(3,3)                             ! Stress vector at end of increment
    real (prec)  ::  press                                    ! Hydrostatic stress
    real (prec)  ::  Bbar(6,length_dof_array)                 ! locking free strain = Bbar*(dof_total+dof_increment)
    real (prec)  ::  Bstar(9,length_dof_array)                ! Locking free F = Bstar*(dof_total+dof_increment)
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
    real (prec)  ::  dW(3,3)                                  ! Spin increment components in global basis
    real (prec)  ::  dR(3,3)                                  ! Rotation increment components in global basis
    real (prec)  ::  strain0(3,3)                             ! Strain at start of increment; updated for rigid rotation
    real (prec)  ::  el_vol
    real (prec)  ::  temp33a(3,3)                             ! Workspace array
    real (prec)  ::  temp33b(3,3)                             ! Workspace array
    real (prec)  ::  temp33c(3,3)                             ! Workspace array
    real (prec)  ::  depsdf(6,9)                              ! Derivative of strain increment with respect to def gradient
    real (prec)  ::  dWdf1(6,9)                               ! part of Derivative of spin stress increment with respect to def gradient
    real (prec)  ::  dWdf2(6,9)                               ! part of derivative of spin stress increment with respect to def grad

    real (prec)  ::  H(3,3),Lam(3,3)                          ! Workspace arrays used to compute stiffness
    real (prec)  ::  SH(3,3),SLam(3,3)                        ! Workspace arrays used to compute stiffness
    real (prec)  ::  dummy                                    ! Dummy variable
    real (prec)  ::  Dcorot(6,6)                              ! Corotational material tangent
    real (prec)  ::  Dspin(6,9)                               ! Spin tangent
    real (prec)  ::  Dtherm(6)                                ! Thermal tangent
!
    real (prec)  ::  dNdxvec(3*n_nodes)                       ! dN/dx stored as a vector
    real (prec)  ::  Pbar(length_dof_array,length_dof_array)  ! d^2 Jbar/ du_i^a du_k^b  matrix
    real (prec)  ::  Pvec(length_dof_array)
    real (prec)  ::  Pmat(length_dof_array,length_dof_array)
    real (prec)  ::  P(length_dof_array,length_dof_array)     !
    real (prec)  ::  S(3,length_dof_array/3)
    real (prec)  ::  Svec(length_dof_array)
    real (prec)  ::  Smat(length_dof_array,length_dof_array)

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
    real (prec)  ::  dummy_state(1)                           ! Dummy state variable array
    integer :: ie
    integer :: i
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
    dNbardx = 0.d0
    Pbar = 0.d0
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
        F1 = F0 + dF
        call invert_small(F1,F1inv,J1)
        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),F1inv)
        dNbardx = dNbardx + J1*dNdx*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
        J0bar = J0bar + det33(F0)*w(kint)*determinant
        J1bar = J1bar + J1*w(kint)*determinant
        dNdxvec(1:3*n_nodes) = reshape(transpose(dNdx(1:n_nodes,1:3)),(/3*n_nodes/))
        do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdx(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
        end do
        P = spread(dNdxvec,dim=1,ncopies=3*n_nodes)*spread(dNdxvec,dim=2,ncopies=3*n_nodes)
        Pbar = Pbar +  J1*(P-Pmat*transpose(Pmat))*w(kint)*determinant
    end do
    J0bar = J0bar/el_vol
    J1bar = J1bar/el_vol
    charLength = el_vol**(1.d0/3.d0)
    dNbardx = dNbardx/(J1bar*el_vol)
    dNdxvec(1:3*n_nodes) = reshape(transpose(dNbardx(1:n_nodes,1:3)),(/3*n_nodes/))
    Pbar = Pbar/(J1bar*el_vol) - spread(dNdxvec,dim=1,ncopies=3*n_nodes)*spread(dNdxvec,dim=2,ncopies=3*n_nodes)

    do kint = 1,n_points
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


        if (n_state_vars_per_intpt>15)   then
          updated_state_variables(iof+15:iof+n_state_vars_per_intpt-1) = &
                                         initial_state_variables(iof+15:iof+n_state_vars_per_intpt-1)
        endif
        SSE = initial_state_variables(iof+12)
        SPD = initial_state_variables(iof+13)
        SCD = initial_state_variables(iof+14)

        abq_time(1:2) = TIME

        if (n_state_vars_per_intpt>15) then
           Call UMAT(stress,updated_state_variables(iof+15:iof+n_state_vars_per_intpt-1),Dcorot,SSE,SPD,SCD, &
                  RPL,Dcorot,DRPLDE,DRPLDT, &
                  strain,strainInc,abq_time,DTIME,BTEMP,BTINC,predef,dpred,material_name, &
                  3,3,6,n_state_vars_per_intpt-15,material_properties,n_properties,xintpt,dR,PNEWDT, &
                  charLength,F0,F1,lmn,kint,0,0,1,current_step_number)
        else
           Call UMAT(stress,dummy_state,Dcorot,SSE,SPD,SCD, &
                  RPL,Dcorot,DRPLDE,DRPLDT, &
                  strain,strainInc,abq_time,DTIME,BTEMP,BTINC,predef,dpred,material_name, &
                  3,3,6,n_state_vars_per_intpt-15,material_properties,n_properties,xintpt,dR,PNEWDT, &
                  charLength,F0,F1,lmn,kint,0,0,1,current_step_number)
        endif
!       Derivative of strain with respect to def gradient, multiplied by F1bar^T
        H = transpose(matmul(F1,Fmidinv))
        Lam = eye3_D - 0.5d0*matmul(dF,Fmidinv)
        depsdf(1,1:9) = [Lam(1,1)*H(1,1),Lam(1,2)*H(1,2),Lam(1,3)*H(1,3), &
                         Lam(1,1)*H(1,2),Lam(1,2)*H(1,1),Lam(1,1)*H(1,3),Lam(1,3)*H(1,1),Lam(1,2)*H(1,3),Lam(1,3)*H(1,2)]
        depsdf(2,1:9) = [Lam(2,1)*H(2,1),Lam(2,2)*H(2,2),Lam(2,3)*H(2,3), &
                         Lam(2,1)*H(2,2),Lam(2,2)*H(2,1),Lam(2,1)*H(2,3),Lam(2,3)*H(2,1),Lam(2,2)*H(2,3),Lam(2,3)*H(2,2)]
        depsdf(3,1:9) = [Lam(3,1)*H(3,1),Lam(3,2)*H(3,2),Lam(3,3)*H(3,3), &
                         Lam(3,1)*H(3,2),Lam(3,2)*H(3,1),Lam(3,1)*H(3,3),Lam(3,3)*H(3,1),Lam(3,2)*H(3,3),Lam(3,3)*H(3,2)]
        depsdf(4,1:9) = [Lam(1,1)*H(2,1)+Lam(2,1)*H(1,1),Lam(1,2)*H(2,2)+Lam(2,2)*H(1,2),Lam(1,3)*H(2,3)+Lam(2,3)*H(1,3), &
                         Lam(1,1)*H(2,2)+Lam(2,1)*H(1,2),Lam(1,2)*H(2,1)+Lam(2,2)*H(1,1),Lam(1,1)*H(2,3)+Lam(2,1)*H(1,3), &
                         Lam(1,3)*H(2,1)+Lam(2,3)*H(1,1),Lam(1,2)*H(2,3)+Lam(2,2)*H(1,3),Lam(1,3)*H(2,2)+Lam(2,3)*H(1,2)]
        depsdf(5,1:9) = [Lam(1,1)*H(3,1)+Lam(3,1)*H(1,1),Lam(1,2)*H(3,2)+Lam(3,2)*H(1,2),Lam(1,3)*H(3,3)+Lam(3,3)*H(1,3), &
                         Lam(1,1)*H(3,2)+Lam(3,1)*H(1,2),Lam(1,2)*H(3,1)+Lam(3,2)*H(1,1),Lam(1,1)*H(3,3)+Lam(3,1)*H(1,3), &
                         Lam(1,3)*H(3,1)+Lam(3,3)*H(1,1),Lam(1,2)*H(3,3)+Lam(3,2)*H(1,3),Lam(1,3)*H(3,2)+Lam(3,3)*H(1,2)]
        depsdf(6,1:9) = [Lam(2,1)*H(3,1)+Lam(3,1)*H(2,1),Lam(2,2)*H(3,2)+Lam(3,2)*H(2,2),Lam(2,3)*H(3,3)+Lam(3,3)*H(2,3), &
                         Lam(2,1)*H(3,2)+Lam(3,1)*H(2,2),Lam(2,2)*H(3,1)+Lam(3,2)*H(2,1),Lam(2,1)*H(3,3)+Lam(3,1)*H(2,3), &
                         Lam(2,3)*H(3,1)+Lam(3,3)*H(2,1),Lam(2,2)*H(3,3)+Lam(3,2)*H(2,3),Lam(2,3)*H(3,2)+Lam(3,3)*H(2,2)]


        stress1 = reshape([stress(1),stress(4),stress(5), &
                           stress(4),stress(2),stress(6), &
                           stress(5),stress(6),stress(3)], &
                           shape(stress1))

        SH = matmul(stress1,H)
        SLam = matmul(stress1,Lam)

        dWdf1(1,1:9) = [Lam(1,1)*SH(1,1),Lam(1,2)*SH(1,2),Lam(1,3)*SH(1,3), &
                        Lam(1,1)*SH(1,2),Lam(1,2)*SH(1,1),Lam(1,1)*SH(1,3),Lam(1,3)*SH(1,1),Lam(1,2)*SH(1,3),Lam(1,3)*SH(1,2)]
        dWdf1(2,1:9) = [Lam(2,1)*SH(2,1),Lam(2,2)*SH(2,2),Lam(2,3)*SH(2,3), &
                        Lam(2,1)*SH(2,2),Lam(2,2)*SH(2,1),Lam(2,1)*SH(2,3),Lam(2,3)*SH(2,1),Lam(2,2)*SH(2,3),Lam(2,3)*SH(2,2)]
        dWdf1(3,1:9) = [Lam(3,1)*SH(3,1),Lam(3,2)*SH(3,2),Lam(3,3)*SH(3,3), &
                        Lam(3,1)*SH(3,2),Lam(3,2)*SH(3,1),Lam(3,1)*SH(3,3),Lam(3,3)*SH(3,1),Lam(3,2)*SH(3,3),Lam(3,3)*SH(3,2)]
        dWdf1(4,1:9) = [Lam(1,1)*SH(2,1)+Lam(2,1)*SH(1,1),Lam(1,2)*SH(2,2)+Lam(2,2)*SH(1,2),Lam(1,3)*SH(2,3)+Lam(2,3)*SH(1,3), &
                        Lam(1,1)*SH(2,2)+Lam(2,1)*SH(1,2),Lam(1,2)*SH(2,1)+Lam(2,2)*SH(1,1),Lam(1,1)*SH(2,3)+Lam(2,1)*SH(1,3), &
                        Lam(1,3)*SH(2,1)+Lam(2,3)*SH(1,1),Lam(1,2)*SH(2,3)+Lam(2,2)*SH(1,3),Lam(1,3)*SH(2,2)+Lam(2,3)*SH(1,2)]
        dWdf1(5,1:9) = [Lam(1,1)*SH(3,1)+Lam(3,1)*SH(1,1),Lam(1,2)*SH(3,2)+Lam(3,2)*SH(1,2),Lam(1,3)*SH(3,3)+Lam(3,3)*SH(1,3), &
                        Lam(1,1)*SH(3,2)+Lam(3,1)*SH(1,2),Lam(1,2)*SH(3,1)+Lam(3,2)*SH(1,1),Lam(1,1)*SH(3,3)+Lam(3,1)*SH(1,3), &
                        Lam(1,3)*SH(3,1)+Lam(3,3)*SH(1,1),Lam(1,2)*SH(3,3)+Lam(3,2)*SH(1,3),Lam(1,3)*SH(3,2)+Lam(3,3)*SH(1,2)]
        dWdf1(6,1:9) = [Lam(2,1)*SH(3,1)+Lam(3,1)*SH(2,1),Lam(2,2)*SH(3,2)+Lam(3,2)*SH(2,2),Lam(2,3)*SH(3,3)+Lam(3,3)*SH(2,3), &
                        Lam(2,1)*SH(3,2)+Lam(3,1)*SH(2,2),Lam(2,2)*SH(3,1)+Lam(3,2)*SH(2,1),Lam(2,1)*SH(3,3)+Lam(3,1)*SH(2,3), &
                        Lam(2,3)*SH(3,1)+Lam(3,3)*SH(2,1),Lam(2,2)*SH(3,3)+Lam(3,2)*SH(2,3),Lam(2,3)*SH(3,2)+Lam(3,3)*SH(2,2)]

        dWdf2(1,1:9) = [SLam(1,1)*H(1,1),SLam(1,2)*H(1,2),SLam(1,3)*H(1,3), &
                        SLam(1,1)*H(1,2),SLam(1,2)*H(1,1),SLam(1,1)*H(1,3),SLam(1,3)*H(1,1),SLam(1,2)*H(1,3),SLam(1,3)*H(1,2)]
        dWdf2(2,1:9) = [SLam(2,1)*H(2,1),SLam(2,2)*H(2,2),SLam(2,3)*H(2,3), &
                        SLam(2,1)*H(2,2),SLam(2,2)*H(2,1),SLam(2,1)*H(2,3),SLam(2,3)*H(2,1),SLam(2,2)*H(2,3),SLam(2,3)*H(2,2)]
        dWdf2(3,1:9) = [SLam(3,1)*H(3,1),SLam(3,2)*H(3,2),SLam(3,3)*H(3,3), &
                        SLam(3,1)*H(3,2),SLam(3,2)*H(3,1),SLam(3,1)*H(3,3),SLam(3,3)*H(3,1),SLam(3,2)*H(3,3),SLam(3,3)*H(3,2)]
        dWdf2(4,1:9) = [SLam(1,1)*H(2,1)+SLam(2,1)*H(1,1),SLam(1,2)*H(2,2)+SLam(2,2)*H(1,2),SLam(1,3)*H(2,3)+SLam(2,3)*H(1,3), &
                        SLam(1,1)*H(2,2)+SLam(2,1)*H(1,2),SLam(1,2)*H(2,1)+SLam(2,2)*H(1,1),SLam(1,1)*H(2,3)+SLam(2,1)*H(1,3), &
                        SLam(1,3)*H(2,1)+SLam(2,3)*H(1,1),SLam(1,2)*H(2,3)+SLam(2,2)*H(1,3),SLam(1,3)*H(2,2)+SLam(2,3)*H(1,2)]
        dWdf2(5,1:9) = [SLam(1,1)*H(3,1)+SLam(3,1)*H(1,1),SLam(1,2)*H(3,2)+SLam(3,2)*H(1,2),SLam(1,3)*H(3,3)+SLam(3,3)*H(1,3), &
                        SLam(1,1)*H(3,2)+SLam(3,1)*H(1,2),SLam(1,2)*H(3,1)+SLam(3,2)*H(1,1),SLam(1,1)*H(3,3)+SLam(3,1)*H(1,3), &
                        SLam(1,3)*H(3,1)+SLam(3,3)*H(1,1),SLam(1,2)*H(3,3)+SLam(3,2)*H(1,3),SLam(1,3)*H(3,2)+SLam(3,3)*H(1,2)]
        dWdf2(6,1:9) = [SLam(2,1)*H(3,1)+SLam(3,1)*H(2,1),SLam(2,2)*H(3,2)+SLam(3,2)*H(2,2),SLam(2,3)*H(3,3)+SLam(3,3)*H(2,3), &
                        SLam(2,1)*H(3,2)+SLam(3,1)*H(2,2),SLam(2,2)*H(3,1)+SLam(3,2)*H(2,1),SLam(2,1)*H(3,3)+SLam(3,1)*H(2,3), &
                        SLam(2,3)*H(3,1)+SLam(3,3)*H(2,1),SLam(2,2)*H(3,3)+SLam(3,2)*H(2,3),SLam(2,3)*H(3,2)+SLam(3,3)*H(2,2)]

        Dspin = 0.5d0*(dWdf1-dWdf2)

        updated_state_variables(iof:iof+5) = stress(1:6)
        updated_state_variables(iof+6:iof+11) = strain(1:6)+strainInc(1:6)
        updated_state_variables(iof+12) = SSE
        updated_state_variables(iof+13) = SPD
        updated_state_variables(iof+14) = SCD

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
            Bbar(i,1:3*n_nodes-2:3) = Bbar(i,1:3*n_nodes-2:3) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3.d0
            Bbar(i,2:3*n_nodes-1:3) = Bbar(i,2:3*n_nodes-1:3) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3.d0
            Bbar(i,3:3*n_nodes:3)   = Bbar(i,3:3*n_nodes:3)   + (dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3.d0
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

        press = sum(stress(1:3))/3.d0

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(Bbar),stress)*w(kint)*determinant*J1bar

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
                 + matmul(transpose(Bbar(1:6,1:3*n_nodes)),matmul(Dcorot,matmul(depsdf(1:6,1:9),Bstar(1:9,1:3*n_nodes)))) &
                                                                                                      *w(kint)*determinant*J1bar &
                 + matmul(transpose(Bbar(1:6,1:3*n_nodes)),matmul(Dspin,Bstar(1:9,1:3*n_nodes)))*w(kint)*determinant*J1bar &
                 - Pmat*transpose(Smat)*w(kint)*determinant*J1bar + press*(Pbar + Pmat*transpose(Pmat))*w(kint)*determinant*J1bar


    end do

    return
end subroutine continuum_element_static_3D



!==========================SUBROUTINE continuum_element_dynamic ==============================
subroutine continuum_element_dynamic_3D(lmn, element_identifier, n_nodes, node_property_list, &               ! Input variables
    density,n_properties, material_properties, material_name, &               ! Input variables
    element_coords, length_coord_array, &                                                       ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                              ! Input variables
    n_state_variables, initial_state_variables, &                                               ! Input variables
    updated_state_variables,element_residual,element_deleted)                                   ! Output variables
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

    real( prec ), intent( in )    :: density                                                ! Density of element
    real( prec ), intent( in )    :: material_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine

    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    logical, intent( inout )  :: element_deleted                                            ! Set to .true. to delete element

    character (len=80), intent( in ) :: material_name                                       ! Name of the material

    ! Local Variables
    integer      :: n_points,kint
    integer      :: n_state_vars_per_intpt, iof


    real (prec)  ::  stress0(3,3)                             ! Initial stress vector contains [s11, s22, s33, s12, s13, s23] Kirchhoff stress
    real (prec)  ::  stress1(3,3)                             ! Stress vector at end of increment
    real (prec)  ::  Bbar(6,length_dof_array)                 ! locking free strain = Bbar*(dof_total+dof_increment)
    real (prec)  ::  xintpt(3)                                ! Coords of integration point
    real (prec)  ::  dxidx(3,3), determinant                  ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  F0(3,3)                                  ! Def grad at start of increment, global basis
    real (prec)  ::  F1(3,3)                                  ! Def grad at end of increment, global basis
    real (prec)  ::  dF(3,3)                                  ! Def gradient increment, global basis
    real (prec)  ::  Fmid(3,3)                                ! Mid point def gradient, global basis
    real (prec)  ::  Fmidinv(3,3)                             ! Inverse of mid point def gradient
    
    
    real (prec)  ::  Jmid                                     ! det(Fmid)
    real (prec)  ::  J0bar                                    ! Vol averaged det(F0)
    real (prec)  ::  J1bar                                    ! Vol averaged det(F1)
    real (prec)  ::  Jmidbar                                  ! Vol averaged det(Fmid)
    real (prec)  ::  U0(3,3), V0(3,3), R0(3,3)                ! Left and right polar decomposition of F0
    real (prec)  ::  U1(3,3), V1(3,3), R1(3,3)                ! Left and right polar decomposition of F1
    real (prec)  ::  dLbar(3,3)                               ! Vol averaged increment in dL_ij
    real (prec)  ::  dLkkbar                                  ! Trace of dLbar
    real (prec)  ::  deps(3,3)                                ! Strain increment in global basis
    real (prec)  ::  depskk                                   ! Trace of strain increment
    real (prec)  ::  dW(3,3)                                  ! Spin increment components in global basis
    real (prec)  ::  dR(3,3)                                  ! Rotation increment components in global basis
    real (prec)  ::  el_vol
    
    
    
    real (prec)  ::  dummy(1)                                 ! Dummy variable
!
    real (prec)  ::  defgradOld(9)                            ! Deformation gradient at start of step in global basis
    real (prec)  ::  defgradNew(9)                            ! Deformation gradient at end of step in global basis
    real (prec)  ::  stretchOld(9)                            ! Stretch U at start of step (= V in corotational components)
    real (prec)  ::  stretchNew(9)                            ! Stretch U at end of step (= V in corotational components)
    real (prec)  ::  stressOld(6)                             ! Stress at start of increment (in co-rotational basis)
    real (prec)  ::  stressNew(6)                             ! Stress at end of increment (in co-rotational basis)
    real (prec)  ::  strainInc(6)                             ! Strain increment (in co-rotational basis)
    real (prec)  ::  relSpinInc(9)                            ! Relative spin increment dW-dR in global basis
    real (prec)  ::  enerInternOld(1)
    real (prec)  ::  enerInternNew(1)
    real (prec)  ::  enerInelasOld(1)
    real (prec)  ::  enerInelasNew(1)
    real (prec)  ::  charLength(1)                               ! Characteristic element length
    real (prec)  ::  density_vumat(1)                          ! Density vector for vumat subroutine argument
    real (prec)  ::  tempold(1)
    real (prec)  ::  tempnew(1)
    real (prec)  ::  fieldold(1)
    real (prec)  ::  fieldnew(1)


    integer :: ie
    integer :: i

!   The state variable storage for 3D continuum elements is
!   s11,s22,s33,s12,s13,s23,
!   e11,e22,e33,e12,e13,e23,
!   specific internal energy, specific plastic dissipation, specific creep dissipation,
!   user defined variables
!
!   Vumat does not use the strain variables or the creep dissipation, but they are needed for a UMAT.

    element_deleted = .false.                               ! Element deletion is suppressed

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    n_state_vars_per_intpt = n_state_variables/n_points

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0

    !   Volume averaged element variables
    J0bar = 0.d0
    J1bar = 0.d0
    Jmidbar = 0.d0
    dLbar = 0.d0
    dNbardx = 0.d0
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
        Fmid = F0 + 0.5d0*dF
        F1 = F0 + dF
        call invert_small(Fmid,Fmidinv,Jmid)
        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),Fmidinv)
        dNbardx = dNbardx + Jmid*dNdx*w(kint)*determinant
        dLbar = dLbar + Jmid*matmul(dF,Fmidinv)*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
        Jmidbar = Jmidbar + Jmid*w(kint)*determinant
        J0bar = J0bar + det33(F0)*w(kint)*determinant
        J1bar = J1bar + det33(F1)*w(kint)*determinant
    end do
    J0bar = J0bar/el_vol
    J1bar = J1bar/el_vol
    charLength(1) = el_vol**(1.d0/3.d0)
    dNbardx = dNbardx/Jmidbar
    dLkkbar = (dLbar(1,1) + dLbar(2,2) + dLbar(3,3))/Jmidbar


    do kint = 1,n_points
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
        F1 = F0 + dF
        F0 = F0*(J0bar/det33(F0))**(1.d0/3.d0)
        F1 = F1*(J1bar/det33(F1))**(1.d0/3.d0)
        Fmid = 0.5d0*(F1+F0)
!
!       Stretch and rotation (global basis)
        call polardecomp(F0,V0,U0,R0)
        call polardecomp(F1,V1,U1,R1)

!       Spatial shape function derivatives
        call invert_small(Fmid,Fmidinv,Jmid)
        dNdx(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),Fmidinv)

!       Strain increment
        deps = matmul(dF,Fmidinv)
        dW = 0.5d0*(deps - transpose(deps))
        deps = 0.5d0*(deps + transpose(deps))
        depskk = deps(1,1) + deps(2,2) + deps(3,3)
        deps(1,1) = deps(1,1) + (dLkkbar - depskk)/3.d0
        deps(2,2) = deps(2,2) + (dLkkbar - depskk)/3.d0
        deps(3,3) = deps(3,3) + (dLkkbar - depskk)/3.d0

!      Rotation increment
        dR = R1-R0

!       Assemble the vumat variables
        defgradOld(1:9) = [F0(1,1),F0(2,2),F0(3,3),&        ! F is in the global basis
                           F0(1,2),F0(2,3),F0(3,1),&
                           F0(2,1),F0(3,2),F0(1,3)]
        defgradNew(1:9) = [F1(1,1),F1(2,2),F1(3,3),&
                           F1(1,2),F1(2,3),F1(3,1),&
                           F1(2,1),F1(3,2),F1(1,3)]
        stretchOld(1:9) = [U0(1,1),U0(2,2),U0(3,3),&        ! This is V in co-rotational basis (same as U in global basis)
                           U0(1,2),U0(2,3),U0(3,1),&
                           U0(2,1),U0(3,2),U0(1,3)]
        stretchNew(1:9) = [U1(1,1),U1(2,2),U1(3,3),&
                           U1(1,2),U1(2,3),U1(3,1),&
                           U1(2,1),U1(3,2),U1(1,3)]
        deps = matmul(transpose(R0),matmul(deps,R0))           !    Co-rotational strain increment components
        strainInc(1:6) = [deps(1,1),deps(2,2),deps(3,3), &
                          deps(1,2),deps(2,3),deps(3,1)]

        iof = n_state_vars_per_intpt*(kint-1) + 1
        stress0 = reshape([initial_state_variables(iof),initial_state_variables(iof+3),initial_state_variables(iof+4), &
                           initial_state_variables(iof+3),initial_state_variables(iof+1),initial_state_variables(iof+5), &
                           initial_state_variables(iof+4),initial_state_variables(iof+5),initial_state_variables(iof+2)], &
                           shape(stress0))
        stress0 = matmul(transpose(R0),matmul(stress0,R0))    !     Co-rotational stress components
        stressOld(1:6) = [stress0(1,1),stress0(2,2),stress0(3,3), &
                          stress0(1,2),stress0(2,3),stress0(3,1)]
        relSpinInc(1:9) = [dW(1,1)-dR(1,1),dW(2,2)-dR(2,2),dW(3,3)-dR(3,3),&   ! dW-dR in global basis
                           dW(1,2)-dR(1,2),dW(2,3)-dR(2,3),dW(3,1)-dR(3,1),&
                           dW(2,1)-dR(2,1),dW(3,2)-dR(3,2),dW(1,3)-dR(1,3)]

        enerInternOld = initial_state_variables(iof+12)
        enerInelasOld = initial_state_variables(iof+13)

        density_vumat(1) = density
        tempold(1) = BTEMP
        tempnew(1) = BTEMP+BTINC
        call vumat(1, 3, 3, n_state_variables-15 ,0, n_properties, 0, &
       TIME+DTIME, TIME+DTIME , DTIME, material_name, xintpt, charLength, &
       material_properties, density_vumat, strainInc, relSpinInc, &
       tempold, stretchOld, defgradOld, fieldold, &
       stressOld, initial_state_variables(iof+15:iof+n_state_variables-1), enerInternOld, enerInelasOld, &
       tempnew, stretchNew, defgradNew, fieldnew, &
       stressNew, updated_state_variables(iof+15:iof+n_state_variables-1), enerInternNew, enerInelasNew )

        stress1 = reshape([stressNew(1),stressNew(4),stressNew(6), &  ! ABAQUS VUMAT storage scheme
                           stressNew(4),stressNew(2),stressNew(5), &
                           stressNew(6),stressNew(5),stressNew(3)], &
                           shape(stress1))
        stress1 = matmul(R1,matmul(stress1,transpose(R1)))    !  Stresses in global basis, UMAT storage convention
        stressNew(1:6) = [stress1(1,1),stress1(2,2),stress1(3,3), &
                          stress1(1,2),stress1(1,3),stress1(2,3)]

        updated_state_variables(iof:iof+5) = [stress1(1,1),stress1(2,2),stress1(3,3),stress1(1,2),stress1(1,3),stress1(2,3)]
        updated_state_variables(iof+12) = enerInternNew(1)
        updated_state_variables(iof+13) = enerInelasNew(1)


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

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(Bbar),stressNew)*w(kint)*determinant*J1bar



    end do

    return
end subroutine continuum_element_dynamic_3D


subroutine continuum_element_fieldvariables_2D(lmn, element_identifier, n_nodes, node_property_list, &        ! Input variables
    density,n_properties, material_properties, material_name, &               ! Input variables
    element_coords, length_coord_array, &                                                               ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                       ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
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

    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! No. projected field variables

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


    character (len=80), intent( in ) :: material_name                                       ! Name of the material
    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: density                                                ! Density of element
    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: material_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step

    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Element stiffness (ROW,COLUMN)

    real (prec)  ::  stress(4),strain(4)                      ! Stress vector contains [s11, s22, s33, s12, s13, s23] Cauchy stress
    real (prec)  ::  dxidx(2,2), determinant                  ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node


    integer :: n_points        ! No. integration points
    integer :: kint            ! Integration point index
    integer :: n_state_vars_per_intpt
    integer :: iof
    integer :: k,i,itst
    integer :: nvar

    real (prec) :: svar

    logical :: strcmp

    character (len=100) fvar_string

    x = reshape(element_coords,(/2,length_coord_array/2/))

    if (n_nodes == 3) n_points = 1
    if (n_nodes == 4) n_points = 4
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 8) n_points = 9

    n_state_vars_per_intpt = n_state_variables/n_points

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0



!   Find stress and project it
    do kint = 1,n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        iof = n_state_vars_per_intpt*(kint-1)
        stress = updated_state_variables(iof+1:iof+4)
        strain = updated_state_variables(iof+5:iof+8)
         do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'U',1) ) then
               fvar_string = field_variable_names(k)
               do i = 2, len(fvar_string)
                        !     Check for anything other than a number in the string
                  itst = ichar(fvar_string(i:i))
                  if ( itst<48 .or. itst>57 ) exit
               end do
               if (i>2) then
                  read(fvar_string(2:i-1), *) nvar
                  if (nvar<n_state_vars_per_intpt-11) then
                   svar = updated_state_variables(iof+11+nvar)
                   nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + svar*N(1:n_nodes)*determinant*w(kint)
                  endif
               endif
            end if
         end do


    end do

    return

end subroutine continuum_element_fieldvariables_2D

subroutine continuum_element_fieldvariables_3D(lmn, element_identifier, n_nodes, node_property_list, &        ! Input variables
    density,n_properties, material_properties, material_name, &               ! Input variables
    element_coords, length_coord_array, &                                                               ! Input variables
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
    integer, intent( in )         :: length_coord_array                                     ! # coordinate variables
    integer, intent( in )         :: length_dof_array                                       ! Total # DOFs
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! No. projected field variables

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


    character (len=80), intent( in ) :: material_name                                       ! Name of the material
    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: density                                                ! Density of element
    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: material_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Element stiffness (ROW,COLUMN)

    real (prec)  ::  stress(6),strain(6)                      ! Stress vector contains [s11, s22, s33, s12, s13, s23] Cauchy stress
    real (prec)  ::  dxidx(3,3), determinant                  ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)                ! Re-shaped coordinate array x(i,a) is ith coord of ath node


    integer :: n_points        ! No. integration points
    integer :: kint            ! Integration point index
    integer :: n_state_vars_per_intpt
    integer :: iof
    integer :: k,i,itst
    integer :: nvar

    real (prec) :: svar

    logical :: strcmp

    character (len=100) fvar_string

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
        strain = updated_state_variables(iof+7:iof+12)
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
            else if (strcmp(field_variable_names(k),'E11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'U',1) ) then
               fvar_string = field_variable_names(k)
               do i = 2, len(fvar_string)
                        !     Check for anything other than a number in the string
                  itst = ichar(fvar_string(i:i))
                  if ( itst<48 .or. itst>57 ) exit
               end do
               if (i>2) then
                  read(fvar_string(2:i-1), *) nvar
                  if (nvar<n_state_vars_per_intpt-15) then
                   svar = updated_state_variables(iof+15+nvar)
                   nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + svar*N(1:n_nodes)*determinant*w(kint)
                  endif
               endif
            end if
         end do


    end do

    return

end subroutine continuum_element_fieldvariables_3D

subroutine continuum_element_lumped_mass(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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

end subroutine continuum_element_lumped_mass

