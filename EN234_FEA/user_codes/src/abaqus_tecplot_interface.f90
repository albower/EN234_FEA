!
!    Interface to allow EN234FEA to plot integration point data defined in ABAQUS UEL or VUEL
!    It is used to compute a least-squared fit of nodal values of user defined variables to data at integration points
!    The nodal values are printed to a data file that can be read by tecplot
!    The subroutine must compute the nodal_variables matrix
!
!
SUBROUTINE EN234FEA_ABAQUS_STATE_PROJECTION(staticstep,n_plot_variables,plot_variable_names, &                                                           ! Field variable definition
    nodal_variables,SVARS,ENERGY,NENERGY,NDOFEL,NSVARS, &
    PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME, &
    KSTEP,KINC,JELEM,PREDEF,NPREDF,JPROPS,NJPROP,PERIOD)

    use Types
    use ParamIO
    implicit none

    logical, intent( in )     :: staticstep           ! Set to TRUE if this is a static analysis (UEL) False means dynamics (VUEL)

    integer, intent( in )     :: n_plot_variables     ! No. nodal variables to be defined (requested by user in input file)
    integer, intent( in )     :: NDOFEL               ! Total # DOF for the element
    integer, intent( in )     :: NENERGY              ! No. energy variables (8 for UEL, 12 for VUEL)
    integer, intent( in )     :: NSVARS               ! Total # state variables for the element
    integer, intent( in )     :: NPROPS               ! Total # real valued properties for the element
    integer, intent( in )     :: MCRD                 ! No. coordinates for the element (same convention as ABAQUS UEL)
    integer, intent( in )     :: NNODE                ! No. nodes for the element
    integer, intent( in )     :: JTYPE                ! Element type specified by Un or VUn in input file
    integer, intent( in )     :: KSTEP                ! Step number
    integer, intent( in )     :: KINC                 ! Increment number
    integer, intent( in )     :: JELEM                ! Element number
    integer, intent( in )     :: NPREDF               ! No. predefined field variables
    integer, intent( in )     :: NJPROP              ! No. integer valued element properties

    integer, intent( in )     :: JPROPS(NJPROP )      ! No. integer valued properties

    real (prec), intent( out )   :: nodal_variables(n_plot_variables,NNODE)        ! Nodal projection of plotted variables
    real (prec), intent( in ) :: SVARS(NSVARS)        ! Element state variables
    real (prec), intent( in ) :: PROPS(NPROPS)        ! Real valued element props
    real (prec), intent( in ) :: ENERGY(NENERGY)      ! Energy variables
    real (prec), intent( in ) :: COORDS(MCRD,NNODE)   ! Element coordinates
    real (prec), intent( in ) :: U(NDOFEL)            ! Nodal DOF at end of step
    real (prec), intent( in ) :: DU(NDOFEL)           ! Increment in nodal DOF
    real (prec), intent( in ) :: V(NDOFEL)            ! Velocity
    real (prec), intent( in ) :: A(NDOFEL)            ! Acceleration
    real (prec), intent( in ) :: TIME(2)              ! Time
    real (prec), intent( in ) :: PREDEF(2,NPREDF,NNODE) ! Predefined field variables
    real (prec), intent( in ) :: DTIME                ! Time increment
    real (prec), intent( in ) :: PERIOD               ! Total analysis time

    character (len=100), intent(in) :: plot_variable_names(n_plot_variables)

    ! Local Variables

    integer      :: n_points,kint,k
    !

    double precision  ::  xi_2D(2,9)                          ! Integration points
    double precision  ::  w_2D(9)                             ! Integration weights
    double precision  ::  N_2D(9)                             ! Shape functions
    double precision  ::  dNdxi_2D(9,2)                       ! Shape function derivatives
    double precision  ::  dxdxi_2D(2,2)                        ! Derivative of position wrt normalized coords

    double precision  ::  xi(3,64)                          ! Integration points
    double precision  ::  w(64)                             ! Integration weights
    double precision  ::  N(20)                             ! Shape functions
    double precision  ::  dNdxi(20,3)                       ! Shape function derivatives
    double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
    !
    double precision  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    double precision  ::  sdev(6), p, smises                ! Deviatoric stress, hydrostatic stress, mises stress
    double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant

    logical      :: strcmp

    !
    !       This subroutine is used by EN234FEA to compute nodal values of element state variables
    !       defined at integration points of a UEL.
    !
    !       It is not used by ABAQUS.
    !
    !       Variables that must be computed in this routine
    !       nodal_variables(i,a)       must contain the integral over the element of the ith nodal state variable
    !                                  multiplied by the shape function for the ath node, i.e. int_V state(i)*N^a dV
    !

    if (MCRD == 2 ) then
    if (NNODE == 3) n_points = 1
    if (NNODE == 4) n_points = 4
    if (NNODE == 6) n_points = 4
    if (NNODE == 8) n_points = 9
    if (NNODE == 9) n_points = 9

    call abq_UEL_2D_integrationpoints(n_points, NNODE, xi_2D, w_2D)

    nodal_variables = 0.d0


    !     --  Loop over integration points
    do kint = 1, n_points
        call abq_UEL_2D_shapefunctions(xi_2D(1:2,kint),NNODE,N_2D,dNdxi_2D)
        dxdxi_2D = matmul(coords(1:2,1:NNODE),dNdxi_2D(1:NNODE,1:2))
        determinant = dxdxi_2D(1,1)*dxdxi_2D(2,2)-dxdxi_2D(2,1)*dxdxi_2D(1,2)

        if (NSVARS>=n_points*4) then   ! Store stress at each integration point (if space was allocated to do so)
            stress(1:4) = SVARS(4*kint-3:4*kint)
        else
            write(6,*) ' Error in subroutine EN234FEA_UEL_STATE_PROJECTION'
            write(6,*) ' SVARS vector does not contain enough variables to extract stresses at all integration points '
            write(6,*) ' NSVARS ',NSVARS,' n_points',n_points
            stop
        endif

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*(4)*sdev(4) )*dsqrt(1.5d0)
        do k = 1,n_plot_variables
            if (strcmp(plot_variable_names(k),'S11',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(1)*N_2D(1:NNODE)*determinant*w_2D(kint)
            else if (strcmp(plot_variable_names(k),'S22',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(2)*N_2D(1:NNODE)*determinant*w_2D(kint)
            else if (strcmp(plot_variable_names(k),'S33',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(3)*N_2D(1:NNODE)*determinant*w_2D(kint)
            else if (strcmp(plot_variable_names(k),'S12',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(4)*N_2D(1:NNODE)*determinant*w_2D(kint)
            else if (strcmp(plot_variable_names(k),'SMISES',6) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + smises*N_2D(1:NNODE)*determinant*w_2D(kint)
            endif
        end do

    end do

    else if (MCRD==3) then

    if (NNODE == 4) n_points = 1
    if (NNODE == 10) n_points = 4
    if (NNODE == 8) n_points = 8
    if (NNODE == 20) n_points = 27

    call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)

    nodal_variables = 0.d0


    !     --  Loop over integration points
    do kint = 1, n_points
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)

        if (NSVARS>=n_points*6) then   ! Store stress at each integration point (if space was allocated to do so)
            stress(1:6) = SVARS(6*kint-5:6*kint)
        else
            write(6,*) ' Error in subroutine EN234FEA_UEL_STATE_PROJECTION'
            write(6,*) ' SVARS vector does not contain enough variables to extract stresses at all integration points '
            write(6,*) ' NSVARS ',NSVARS,' n_points',n_points
            stop
        endif

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        do k = 1,n_plot_variables
            if (strcmp(plot_variable_names(k),'S11',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(1)*N(1:NNODE)*determinant*w(kint)
            else if (strcmp(plot_variable_names(k),'S22',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(2)*N(1:NNODE)*determinant*w(kint)
            else if (strcmp(plot_variable_names(k),'S33',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(3)*N(1:NNODE)*determinant*w(kint)
            else if (strcmp(plot_variable_names(k),'S12',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(4)*N(1:NNODE)*determinant*w(kint)
            else if (strcmp(plot_variable_names(k),'S13',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(5)*N(1:NNODE)*determinant*w(kint)
            else if (strcmp(plot_variable_names(k),'S23',3) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + stress(6)*N(1:NNODE)*determinant*w(kint)
            else if (strcmp(plot_variable_names(k),'SMISES',6) ) then
                nodal_variables(k,1:NNODE) = nodal_variables(k,1:NNODE) + smises*N(1:NNODE)*determinant*w(kint)
            endif
        end do

    end do
    end if

!
End Subroutine EN234FEA_ABAQUS_STATE_PROJECTION
