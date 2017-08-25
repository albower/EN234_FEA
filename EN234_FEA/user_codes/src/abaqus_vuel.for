!
!    ABAQUS format VUEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Explicit
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          EN234FEA_VUEL_STATE_PROJECTION          -  used internally by EN234FEA to compute nodal values of integration point variables
!          abq_VUEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_VUEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_VUEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_VUEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!          abq_VUEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_VUEL_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE VUEL(nblock,rhs,amass,dtimeStable,svars,nsvars,
     1  energy,nnode,ndofel,props,nprops,jprops,njprops,
     2  coords,mcrd,u,du,v,a,
     3  jtype,jElem,
     4  time,period,dtimeCur,dtimePrev,kstep,kinc,
     5  lflags,
     6  dMassScaleFactor,
     7  predef,npredef,
     8  jdltyp, adlmag)
    !
      include 'vaba_param.inc'
    ! operational code keys

      parameter ( jMassCalc = 1,   ! lflags(iOpCode)=jMassCalc => calculate mass matrix
     1   jIntForceAndDtStable = 2,  ! lflags(iOpCode) = jIntForceAndDtStable => calculate internal force and stable time step
     2   jExternForce = 3)          ! lflags(iOpCode)=jExternForce => calculate distribted forces
    ! flag indices
      parameter (iProcedure = 1,    ! Index of procedure flag in lflags array
     1   iNlgeom = 2,               ! Index of nonlinear geometry flag in lflags array
     2   iOpCode = 3,               ! Index of operation code in lflags array
     3   nFlags = 3)                ! Dimension of lflags array
    ! energy array indices
      parameter ( iElPd = 1,        ! Plastic Dissipation
     1   iElCd = 2,                 ! Creep Dissipation
     2   iElIe = 3,                 ! Internal Energy
     3   iElTs = 4,                 ! Transverse shear energy
     4   iElDd = 5,                 ! Material Damping Dissipation
     5   iElBv = 6,
     6   iElDe = 7,
     7   iElHe = 8,
     8   iElKe = 9,
     9   iElTh = 10,
     1   iElDmd = 11,
     2   iElDc = 12,
     3   nElEnergy = 12)
    ! predefined variables indices
      parameter ( iPredValueNew = 1,
     2   iPredValueOld = 2,
     3   nPred = 2)
    ! time indices
      parameter (iStepTime = 1,
     1   iTotalTime = 2,
     2   nTime = 2)

      dimension rhs(nblock,ndofel),amass(nblock,ndofel,ndofel),
     1   dtimeStable(nblock),
     2   svars(nblock,nsvars),energy(nblock,nElEnergy),
     3   props(nprops),jprops(njprops),
     4   jElem(nblock),time(nTime),lflags(nFlags),
     5   coords(nblock,nnode,mcrd),
     6   u(nblock,ndofel), du(nblock,ndofel),
     7   v(nblock,ndofel), a(nblock, ndofel),
     8   dMassScaleFactor(nblock),
     9   predef(nblock,nnode,npredef,nPred),
     1   adlmag(nblock)
    !
    !       Variables that must be computed in this routine
    !       RHS(n,i)                   Right hand side vector for nth element block.  In EN234_FEA the dimensions are always RHS(1,ndofel)
    !       amass(n,i,j)               Mass matrix for nth element block
    !                                  The mass matrix for translational DOF in ABAQUS/Explicit must be diagonal.
    !                                  Rotational DOF may be coupled (EN234_FEA assumes everything is diagonal)
    !       SVARS(n,1:NSVARS)          Element state variables for nth element block.  Must be updated in this routine
    !       ENERGY(n,1:nElEnergy)      Energy measures for nth element block.  See parameter list for definitions


    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element (# nodes x # dof per node)
    !       NSVARS                     Total # element state variables
    !       NNODE                      No. nodes on an element
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(N,i)                ith coordinate of Nth node on element (Note this order differs from UEL)
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U(i,m)                     mth dof for ith element at end of step
    !       DU(i,m)                    increment in mth dof for ith element
    !       V(i,m)                     mth velocity (rate of change of dof, using midpoint differencing) for ith element
    !       A(i,m)                     mth acceleration (generalized force/lumped mass) for ith element
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       JDLTYP                     Integer n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(i)                  Distributed load magnitude for ith element block
    !       PREDEF(i,1:NNODE,1:npredf,1:2)  Predefined fields for nodes in ith element block.
    !       PREDEF(....,1)             Value of predefined field at end of step
    !       PREDEF(...,2)              Value of predefined field at start of step
    !       PREDEF(i,k,1,1:2)          Value of temperature at kth node
    !       PREDEF(i,k,2:npedf,1:2)    Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable - see parameter list for definitions

    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      integer      :: i,j,n_points,kint, nfacenodes, ipoin
      integer      :: face_node_list(8)                       ! List of nodes on an element face
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt spatial coords
    !
    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(3,8)                  ! Coords of nodes on an element face
      double precision  ::  xi2(2,9)                          ! Area integration points
      double precision  ::  N2(9)                             ! 2D shape functions
      double precision  ::  dNdxi2(9,2)                       ! 2D shape function derivatives
      double precision  ::  norm(3)                           ! Normal to an element face
      double precision  ::  dxdxi2(3,2)                       ! Derivative of spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  D(6,6)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(6,60)                           ! strain = B*(dof_total)
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      double precision  ::  E, xnu, D44, D11, D12             ! Material properties

    !
    !     Example ABAQUS VUEL implementing 3D linear elastic elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio
    !     PROPS(3)         Mass Density

      if (ndofel<3*NNODE) then
        write(6,*) ' Error in abaqus VUEL '
        write(6,*) ' Variable ndofel must exceed 3*NNODE'
        write(6,*) ' ndofel = ',ndofel,' NNODE = ',NNODE
      endif

      RHS(1:nblock,1:ndofel) = 0.d0

      ENERGY(1:nblock,1:nElEnergy) = 0.d0

      AMASS(1:nblock,1:ndofel,1:ndofel) = 0.d0


      if (lflags(iOpCode)==jMassCalc) then   !       Compute a lumped mass matrix
        n_points = 0
        if (NNODE == 4) n_points = 4
        if (NNODE == 10) n_points = 5
        if (NNODE == 8) n_points = 27
        if (NNODE == 20) n_points = 27
        if (NNODE==0) then                                 ! Nonstandard element
            write(6,*) ' Error in abaqus VUEL '
            write(6,*) ' Mass matrix has not been coded for',
     1       ' element with ',NNODE,' nodes'
            stop
        endif
        call abq_VUEL_3D_integrationpoints(n_points, NNODE, xi, w)

        do iblock = 1,nblock

            do kint = 1,n_points
                call abq_VUEL_3D_shapefunctions(xi(1:3,kint),
     1                                              NNODE,N,dNdxi)
                dxdxi = matmul(transpose(coords(iblock,1:NNODE,1:3)),
     1                                              dNdxi(1:NNODE,1:3))
                determinant =   dxdxi(1,1)*dxdxi(2,2)*dxdxi(3,3)
     1                         - dxdxi(1,1)*dxdxi(2,3)*dxdxi(3,2)
     2                         - dxdxi(1,2)*dxdxi(2,1)*dxdxi(3,3)
     3                         + dxdxi(1,2)*dxdxi(2,3)*dxdxi(3,1)
     4                         + dxdxi(1,3)*dxdxi(2,1)*dxdxi(3,2)
     5                         - dxdxi(1,3)*dxdxi(2,2)*dxdxi(3,1)
                !             Lumped mass matrix computed  using row sum method
                irow = 0
                do j = 1,NNODE
                    irow = 3*j-2
                    amass(iblock,irow,irow) = amass(iblock,irow,irow)
     1                       + N(j)*sum(N(1:NNODE))*determinant*w(kint)
                    amass(iblock,irow+1,irow+1) =
     1                                  amass(iblock,irow+1,irow+1)
     2                       + N(j)*sum(N(1:NNODE))*determinant*w(kint)
                    amass(iblock,irow+2,irow+2) =
     1                            amass(iblock,irow+2,irow+2)
     2                       + N(j)*sum(N(1:NNODE))*determinant*w(kint)
                end do
            end do

        end do

        amass = amass*props(3)


      else if (lflags(iOpCode)==jIntForceAndDtStable) then ! Compute internal force vector and stable time increment

        dtimeStable = dtimecur              ! This leaves time increment unchanged

        if (NNODE == 4) n_points = 1               ! Linear tet
        if (NNODE == 10) n_points = 4              ! Quadratic tet
        if (NNODE == 8) n_points = 8               ! Linear Hex
        if (NNODE == 20) n_points = 27             ! Quadratic hex

        call abq_VUEL_3D_integrationpoints(n_points, NNODE, xi, w)

        D = 0.d0
        E = PROPS(1)
        xnu = PROPS(2)
        d44 = 0.5D0*E/(1+xnu)
        d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
        d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
        D(1:3,1:3) = d12
        D(1,1) = d11
        D(2,2) = d11
        D(3,3) = d11
        D(4,4) = d44
        D(5,5) = d44
        D(6,6) = d44

        do iblock = 1,nblock
            !     --  Loop over integration points
            do kint = 1, n_points
                call abq_VUEL_3D_shapefunctions(xi(1:3,kint),
     1                                               NNODE,N,dNdxi)
                dxdxi = matmul(transpose(coords(iblock,1:NNODE,1:3)),
     1                                            dNdxi(1:NNODE,1:3))
                call abq_VUEL_invert3d(dxdxi,dxidx,determinant)
                dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
                B = 0.d0
                B(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
                B(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
                B(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
                B(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
                B(4,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
                B(5,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
                B(5,3:3*NNODE:3)   = dNdx(1:NNODE,1)
                B(6,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
                B(6,3:3*NNODE:3)   = dNdx(1:NNODE,2)

                strain = matmul(B(1:6,1:3*NNODE),U(iblock,1:3*NNODE))

                stress = matmul(D,strain)

                RHS(iblock,1:3*NNODE) = RHS(iblock,1:3*NNODE)
     1           - matmul(transpose(B(1:6,1:3*NNODE)),stress(1:6))*
     2                                              w(kint)*determinant

                ENERGY(iblock,2) = ENERGY(iblock,2)
     1           + 0.5D0*dot_product(stress,strain)*w(kint)*determinant           ! Store the elastic strain energy

                if (NSVARS>=n_points*6) then   ! Store stress at each integration point (if space was allocated to do so)
                    SVARS(iblock,6*kint-5:6*kint) = stress(1:6)
                endif
            end do
        end do


      else if (lflags(iOpCode)==jExternForce) then
        !
        !   Apply distributed loads
        !
        !   Distributed loads are specified in the input file using the Un option in the input file.
        !   n specifies the face number, following the ABAQUS convention
        !


        do iblock = 1,nblock

            call abq_VUEL_facenodes_3D(NNODE,iabs(JDLTYP),
     1                                     face_node_list,nfacenodes)

            do i = 1,nfacenodes
                face_coords(1:3,i) =
     1                             coords(iblock,face_node_list(i),1:3)
            end do

            if (nfacenodes == 3) n_points = 3
            if (nfacenodes == 6) n_points = 4
            if (nfacenodes == 4) n_points = 4
            if (nfacenodes == 8) n_points = 9

            call abq_VUEL_2D_integrationpoints(n_points,
     1                                             nfacenodes, xi2, w)

            do kint = 1,n_points
                call abq_VUEL_2D_shapefunctions(xi2(1:2,kint),
     1                                            nfacenodes,N2,dNdxi2)
                dxdxi2 = matmul(face_coords(1:3,1:nfacenodes),
     1                                    dNdxi2(1:nfacenodes,1:2))
                norm(1) = (dxdxi2(2,1)*dxdxi2(3,2))
     1                                       -(dxdxi2(2,2)*dxdxi2(3,1))
                norm(2) = (dxdxi2(1,1)*dxdxi2(3,2))
     1                                       -(dxdxi2(1,2)*dxdxi2(3,1))
                norm(3) = (dxdxi2(1,1)*dxdxi2(2,2))
     1                                        -(dxdxi2(1,2)*dxdxi2(2,1))

                do i = 1,nfacenodes
                    ipoin = 3*face_node_list(i)-2
                    RHS(iblock,ipoin:ipoin+2) =
     1                            RHS(iblock,ipoin:ipoin+2)
     2                        - N2(i)*adlmag(iblock)*norm(1:3)*w(kint)      ! Note determinant is already in normal
                end do
            end do
        end do

      else
        write(6,*) ' Error in ABAQUS subroutine VUEL '
        write(6,*) ' Unrecognized LFLAGS option ',lflags(iOpCode)
        stop

      endif

      return

      END SUBROUTINE VUEL


      subroutine abq_VUEL_3D_integrationpoints(n_points, n_nodes, xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(3,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision x1D(4), w1D(4)

    !         Defines integration points and weights for 3D continuum elements

      if (n_nodes  == 4.or.n_nodes ==10 ) then   ! Tetrahedral elements
        if (n_points == 1) then
            xi(1,1) = 0.25D0
            xi(2,1) = 0.25D0
            xi(3,1) = 0.25D0
            w(1) = 1.D0/6.D0
        else if (n_points == 4) then
            xi(1,1) = 0.58541020
            xi(2,1) = 0.13819660
            xi(3,1) = xi(2,1)
            xi(1,2) = xi(2,1)
            xi(2,2) = xi(1,1)
            xi(3,2) = xi(2,1)
            xi(1,3) = xi(2,1)
            xi(2,3) = xi(2,1)
            xi(3,3) = xi(1,1)
            xi(1,4) = xi(2,1)
            xi(2,4) = xi(2,1)
            xi(3,4) = xi(2,1)
            w(1:4) = 1.D0/24.D0
        else if (n_points == 5) then
            xi(1,1) = 0.25d0
            xi(2,1) = 0.25d0
            xi(3,1) = 0.25d0
            xi(1,2) = 0.5d0
            xi(2,2) = 1.d0/6.d0
            xi(3,2) = 1.d0/6.d0
            xi(1,3) = 1.d0/6.d0
            xi(2,3) = 0.5d0
            xi(3,3) = 1.d0/6.d0
            xi(1,4) = 1.d0/6.d0
            xi(2,4) = 1.d0/6.d0
            xi(3,4) = 0.5d0
            xi(1,5) = 1.d0/6.d0
            xi(2,5) = 1.d0/6.d0
            xi(3,5) = 1.d0/6.d0
            w(1) = -4.d0/30.d0
            w(2:5) = 3.d0/40.d0
        else
            write(6,*) ' Incorrect number of integration points for',
     1                                    ' tetrahedral element '
            write(6, *) ' called with ',n_points
            stop
        endif
      else if ( n_nodes == 8 .or. n_nodes == 20 ) then   ! 8 or 20 noded hexahedral elements
        if (n_points == 1) then
            xi(1,1) = 0.D0
            xi(2,1) = 0.D0
            xi(3,1) = 0.D0
            w(1) = 8.D0
        else if (n_points == 8) then
            x1D(1) = -0.5773502692
            x1D(2) =  0.5773502692
            do k = 1,2
                do j = 1,2
                    do i = 1,2
                        n = 4*(k-1) + 2*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                    end do
                end do
            end do
            w(1:8) = 1.D0
        else if (n_points == 27) then
            x1D(1) = -0.7745966692
            x1D(2) = 0.
            x1D(3) = 0.7745966692
            w1D(1) = 0.5555555555D0
            w1D(2) = 0.888888888D0
            w1D(3) = 0.55555555555D0
            do k = 1,3
                do j = 1,3
                    do i = 1,3
                        n = 9*(k-1) + 3*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                        w(n) = w1D(i)*w1D(j)*w1D(k)
                    end do
                end do
            end do
        else if (n_points == 64) then
            x1D(1) = .8611363115940526D+00
            x1D(2) = .3399810435848563D+00
            x1D(3) = -.3399810435848563D+00
            x1D(4) = -.8611363115940526D+00
            w1D(1) = .3478548451374538D+00
            w1D(2) = .6521451548625461D+00
            w1D(3) = .6521451548625461D+00
            w1D(4) = .3478548451374538D+00
            do k = 1,4
                do j = 1,4
                    do i = 1,4
                        n = 16*(k-1) + 4*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                        w(n) = w1D(i)*w1D(j)*w1D(k)
                    end do
                end do
            end do
        endif
      endif

      return

      end subroutine abq_VUEL_3D_integrationpoints

      subroutine abq_VUEL_2D_integrationpoints(n_points, n_nodes, xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(2,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision :: cn,w1,w2,w11,w12,w22

    !         Defines integration points and weights for 2D continuum elements

      if ( n_points==1 ) then
        if ( n_nodes==4 .or. n_nodes==9 ) then    !     ---   4 or 9 noded quad
            xi(1, 1) = 0.D0
            xi(2, 1) = 0.D0
            w(1) = 4.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then !     ---   3 or 6 noded triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = 1.D0/3.D0
            w(1) = 1.D0/2.D0
        end if
      else if ( n_points==3 ) then
        xi(1, 1) = 0.5D0
        xi(2, 1) = 0.5D0
        w(1) = 1.D0/6.D0
        xi(1, 2) = 0.D0
        xi(2, 2) = 0.5D0
        w(2) = w(1)
        xi(1, 3) = 0.5D0
        xi(2, 3) = 0.D0
        w(3) = w(1)
      else if ( n_points==4 ) then
        if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
            !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
            !     43
            !     12
            cn = 0.5773502691896260D0
            xi(1, 1) = -cn
            xi(1, 2) = cn
            xi(1, 3) = cn
            xi(1, 4) = -cn
            xi(2, 1) = -cn
            xi(2, 2) = -cn
            xi(2, 3) = cn
            xi(2, 4) = cn
            w(1) = 1.D0
            w(2) = 1.D0
            w(3) = 1.D0
            w(4) = 1.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then
            !     xi integration points for triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = xi(1, 1)
            w(1) = -27.D0/96.D0
            xi(1, 2) = 0.6D0
            xi(2, 2) = 0.2D0
            w(2) = 25.D0/96.D0
            xi(1, 3) = 0.2D0
            xi(2, 3) = 0.6D0
            w(3) = w(2)
            xi(1, 4) = 0.2D0
            xi(2, 4) = 0.2D0
            w(4) = w(2)
        end if

      else if ( n_points==7 ) then
        ! Quintic integration for triangle
        xi(1,1) = 1.d0/3.d0
        xi(2,1) = xi(1,1)
        w(1) = 0.1125d0
        xi(1,2) = 0.0597158717d0
        xi(2,2) = 0.4701420641d0
        w(2) = 0.0661970763d0
        xi(1,3) = xi(2,2)
        xi(2,3) = xi(1,2)
        w(3) = w(2)
        xi(1,4) = xi(2,2)
        xi(2,4) = xi(2,2)
        w(4) = w(2)
        xi(1,5) = 0.7974269853d0
        xi(2,5) = 0.1012865073d0
        w(5) = 0.0629695902d0
        xi(1,6) = xi(2,5)
        xi(2,6) = xi(1,5)
        w(6) = w(5)
        xi(1,7) = xi(2,5)
        xi(2,7) = xi(2,5)
        w(7) = w(5)
      else if ( n_points==9 ) then
        !     3X3 GAUSS INTEGRATION POINTS
        !     789
        !     456
        !     123
        cn = 0.7745966692414830D0
        xi(1, 1) = -cn
        xi(1, 2) = 0.D0
        xi(1, 3) = cn
        xi(1, 4) = -cn
        xi(1, 5) = 0.D0
        xi(1, 6) = cn
        xi(1, 7) = -cn
        xi(1, 8) = 0.D0
        xi(1, 9) = cn
        xi(2, 1) = -cn
        xi(2, 2) = -cn
        xi(2, 3) = -cn
        xi(2, 4) = 0.D0
        xi(2, 5) = 0.D0
        xi(2, 6) = 0.D0
        xi(2, 7) = cn
        xi(2, 8) = cn
        xi(2, 9) = cn
        w1 = 0.5555555555555560D0
        w2 = 0.8888888888888890D0
        w11 = w1*w1
        w12 = w1*w2
        w22 = w2*w2
        w(1) = w11
        w(2) = w12
        w(3) = w11
        w(4) = w12
        w(5) = w22
        w(6) = w12
        w(7) = w11
        w(8) = w12
        w(9) = w11
      end if

      return

      end subroutine abq_VUEL_2D_integrationpoints


      subroutine abq_VUEL_3D_shapefunctions(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes

      double precision, intent(in) :: xi(3)
      double precision, intent(out) :: f(20)
      double precision, intent(out) :: df(20,3)
      double precision xi4

!   Defines shape functions for 3D continuum elements

      if (n_nodes == 4) then
        f(1) = xi(1)
        f(2) = xi(2)
        f(3) = xi(3)
        f(4) = 1.-xi(1)-xi(2)-xi(3)
        df(1,1) = 1.
        df(2,2) = 1.
        df(3,3) = 1.
        df(4,1) = -1.
        df(4,2) = -1.
        df(4,3) = -1.
      else if (n_nodes == 10) then
        xi4 = 1.D0-xi(1)-xi(2)-xi(3)
        f(1) = (2.*xi(1)-1.)*xi(1)
        f(2) = (2.*xi(2)-1.)*xi(2)
        f(3) = (2.*xi(3)-1.)*xi(3)
        f(4) = (2.*xi4-1.)*xi4
        f(5) = 4.*xi(1)*xi(2)
        f(6) = 4.*xi(2)*xi(3)
        f(7) = 4.*xi(3)*xi(1)
        f(8) = 4.*xi(1)*xi4
        f(9) = 4.*xi(2)*xi4
        f(10) = 4.*xi(3)*xi4
        df(1,1) = (4.*xi(1)-1.)
        df(2,2) = (4.*xi(2)-1.)
        df(3,3) = (4.*xi(3)-1.)
        df(4,1) = -(4.*xi4-1.)
        df(4,2) = -(4.*xi4-1.)
        df(4,3) = -(4.*xi4-1.)
        df(5,1) = 4.*xi(2)
        df(5,2) = 4.*xi(1)
        df(6,2) = 4.*xi(3)
        df(6,3) = 4.*xi(2)
        df(7,1) = 4.*xi(3)
        df(7,3) = 4.*xi(1)
        df(8,1) = 4.*(xi4-xi(1))
        df(8,2) = -4.*xi(1)
        df(8,3) = -4.*xi(1)
        df(9,1) = -4.*xi(2)
        df(9,2) = 4.*(xi4-xi(2))
        df(9,3) = -4.*xi(2)
        df(10,1) = -4.*xi(3)*xi4
        df(10,2) = -4.*xi(3)
        df(10,3) = 4.*(xi4-xi(3))
      else if (n_nodes == 8) then
        f(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))/8.
        f(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))/8.
        f(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))/8.
        f(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))/8.
        f(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))/8.
        f(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))/8.
        f(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))/8.
        f(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))/8.
        df(1,1) = -(1.-xi(2))*(1.-xi(3))/8.
        df(1,2) = -(1.-xi(1))*(1.-xi(3))/8.
        df(1,3) = -(1.-xi(1))*(1.-xi(2))/8.
        df(2,1) = (1.-xi(2))*(1.-xi(3))/8.
        df(2,2) = -(1.+xi(1))*(1.-xi(3))/8.
        df(2,3) = -(1.+xi(1))*(1.-xi(2))/8.
        df(3,1) = (1.+xi(2))*(1.-xi(3))/8.
        df(3,2) = (1.+xi(1))*(1.-xi(3))/8.
        df(3,3) = -(1.+xi(1))*(1.+xi(2))/8.
        df(4,1) = -(1.+xi(2))*(1.-xi(3))/8.
        df(4,2) = (1.-xi(1))*(1.-xi(3))/8.
        df(4,3) = -(1.-xi(1))*(1.+xi(2))/8.
        df(5,1) = -(1.-xi(2))*(1.+xi(3))/8.
        df(5,2) = -(1.-xi(1))*(1.+xi(3))/8.
        df(5,3) = (1.-xi(1))*(1.-xi(2))/8.
        df(6,1) = (1.-xi(2))*(1.+xi(3))/8.
        df(6,2) = -(1.+xi(1))*(1.+xi(3))/8.
        df(6,3) = (1.+xi(1))*(1.-xi(2))/8.
        df(7,1) = (1.+xi(2))*(1.+xi(3))/8.
        df(7,2) = (1.+xi(1))*(1.+xi(3))/8.
        df(7,3) = (1.+xi(1))*(1.+xi(2))/8.
        df(8,1) = -(1.+xi(2))*(1.+xi(3))/8.
        df(8,2) = (1.-xi(1))*(1.+xi(3))/8.
        df(8,3) = (1.-xi(1))*(1.+xi(2))/8.
      else if (n_nodes == 20) then
        f(1)=(1.-xi(1))*(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)/8.
        f(2)=(1.+xi(1))*(1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)/8.
        f(3)=(1.+xi(1))*(1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)/8.
        f(4)=(1.-xi(1))*(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)/8.
        f(5)=(1.-xi(1))*(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)/8.
        f(6)=(1.+xi(1))*(1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)/8.
        f(7)=(1.+xi(1))*(1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)/8.
        f(8)=(1.-xi(1))*(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)/8.
        f(9) = (1.-xi(1)**2.)*(1.-xi(2))*(1.-xi(3))/4.
        f(10) = (1.+xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
        f(11) = (1.-xi(1)**2.)*(1.+xi(2))*(1.-xi(3))/4.
        f(12) = (1.-xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
        f(13) = (1.-xi(1)**2.)*(1.-xi(2))*(1.+xi(3))/4.
        f(14) = (1.+xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
        f(15) = (1.-xi(1)**2.)*(1.+xi(2))*(1.+xi(3))/4.
        f(16) = (1.-xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
        f(17) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
        f(18) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
        f(19) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
        f(20) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
        df(1,1) = (-(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(1,2) = (-(1.-xi(1))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(1,3) = (-(1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.

        df(2,1) = ((1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(2,2) = (-(1.+xi(1))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     1          -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(2,3) = (-(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-xi(3)-2.)
     1           -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.

        df(3,1) = ((1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(3,2) = ((1.+xi(1))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(3,3) = (-(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-xi(3)-2.)
     1           -(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.

        df(4,1) = (-(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(4,2) = ((1.-xi(1))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     1            +(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(4,3) = (-(1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(5,1) = (-(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(5,2) = (-(1.-xi(1))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(5,3) = ((1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,1) = ((1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,2) = (-(1.+xi(1))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     1           -(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,3) = ((1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(7,1) = ((1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(7,2) = ((1.+xi(1))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(7,3) = ((1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,1) = (-(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,2) = ((1.-xi(1))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,3) = ((1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(9,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
        df(9,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
        df(9,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
        df(10,1)  = (1.-xi(2)**2.)*(1.-xi(3))/4.
        df(10,2)  = -2.*xi(2)*(1.+xi(1))*(1.-xi(3))/4.
        df(10,3)  = -(1.-xi(2)**2.)*(1.+xi(1))/4.
        df(11,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
        df(11,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
        df(11,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
        df(12,1)  = -(1.-xi(2)**2.)*(1.-xi(3))/4.
        df(12,2)  = -2.*xi(2)*(1.-xi(1))*(1.-xi(3))/4.
        df(12,3)  = -(1.-xi(2)**2.)*(1.-xi(1))/4.
        df(13,1)  = -2.*xi(1)*(1.-xi(2))*(1.+xi(3))/4.
        df(13,2)  = -(1.-xi(1)**2.)*(1.+xi(3))/4.
        df(13,3)  = (1.-xi(1)**2.)*(1.-xi(2))/4.
        df(14,1)  = (1.-xi(2)**2.)*(1.+xi(3))/4.
        df(14,2)  = -2.*xi(2)*(1.+xi(1))*(1.+xi(3))/4.
        df(14,3)  = (1.-xi(2)**2.)*(1.+xi(1))/4.
        df(15,1)  = 2.*xi(1)*(1.+xi(2))*(1.+xi(3))/4.
        df(15,2)  = (1.-xi(1)**2.)*(1.+xi(3))/4.
        df(15,3)  = (1.-xi(1)**2.)*(1.+xi(2))/4.
        df(16,1)  = -(1.-xi(2)**2.)*(1.+xi(3))/4.
        df(16,2)  = -2.*xi(2)*(1.-xi(1))*(1.+xi(3))/4.
        df(16,3)  = (1.-xi(2)**2.)*(1.-xi(1))/4.
        df(17,1) = -(1.-xi(2))*(1.-xi(3)**2.)/4.
        df(17,2) = -(1.-xi(1))*(1.-xi(3)**2.)/4.
        df(17,3) = -xi(3)*(1.-xi(1))*(1.-xi(2))/2.
        df(18,1) = (1.-xi(2))*(1.-xi(3)**2.)/4.
        df(18,2) = -(1.+xi(1))*(1.-xi(3)**2.)/4.
        df(18,3) = -xi(3)*(1.+xi(1))*(1.-xi(2))/2.
        df(19,1) = (1.+xi(2))*(1.-xi(3)**2.)/4.
        df(19,2) = (1.+xi(1))*(1.-xi(3)**2.)/4.
        df(19,3) = -xi(3)*(1.+xi(1))*(1.+xi(2))/2.
        df(20,1) = -(1.+xi(2))*(1.-xi(3)**2.)/4.
        df(20,2) = (1.-xi(1))*(1.-xi(3)**2.)/4.
        df(20,3) = -xi(3)*(1.-xi(1))*(1.+xi(2))/2.
      endif

      end subroutine abq_VUEL_3D_shapefunctions



      subroutine abq_VUEL_2D_shapefunctions(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes

      double precision, intent(in) :: xi(2)
      double precision, intent(out) :: f(*)
      double precision, intent(out) :: df(9,2)
      double precision g1, g2, g3, dg1, dg2, dg3
      double precision h1, h2, h3, dh1, dh2, dh3
      double precision z,dzdp, dzdq

      if ( n_nodes==3 ) then        !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
        f(1) = xi(1)
        f(2) = xi(2)
        f(3) = 1.D0 - xi(1) - xi(2)
        df(1, 1) = 1.D0
        df(1, 2) = 0.D0
        df(2, 1) = 0.D0
        df(2, 2) = 1.D0
        df(3, 1) = -1.D0
        df(3, 2) = -1.D0
      else if ( n_nodes==4 ) then
        !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
        !     43
        !     12
        g1 = 0.5D0*(1.D0 - xi(1))
        g2 = 0.5D0*(1.D0 + xi(1))
        h1 = 0.5D0*(1.D0 - xi(2))
        h2 = 0.5D0*(1.D0 + xi(2))
        f(1) = g1*h1
        f(2) = g2*h1
        f(3) = g2*h2
        f(4) = g1*h2
        dg1 = -0.5D0
        dg2 = 0.5D0
        dh1 = -0.5D0
        dh2 = 0.5D0
        df(1, 1) = dg1*h1
        df(2, 1) = dg2*h1
        df(3, 1) = dg2*h2
        df(4, 1) = dg1*h2
        df(1, 2) = g1*dh1
        df(2, 2) = g2*dh1
        df(3, 2) = g2*dh2
        df(4, 2) = g1*dh2

      else if ( n_nodes==6 ) then

        !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
        !          3

        !       6      5

        !     1    4     2

        !     P = L1
        !     Q = L2
        !     Z = 1 - P - Q = L3

        z = 1.D0 - xi(1) - xi(2)
        f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
        f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
        f(3) = (2.D0*z - 1.D0)*z
        f(4) = 4.D0*xi(1)*xi(2)
        f(5) = 4.D0*xi(2)*z
        f(6) = 4.D0*xi(1)*z
        dzdp = -1.D0
        dzdq = -1.D0
        df(1, 1) = 4.D0*xi(1) - 1.D0
        df(2, 1) = 0.D0
        df(3, 1) = 4.D0*z*dzdp - dzdp
        df(4, 1) = 4.D0*xi(2)
        df(5, 1) = 4.D0*xi(2)*dzdp
        df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
        df(1, 2) = 0.D0
        df(2, 2) = 4.D0*xi(2) - 1.D0
        df(3, 2) = 4.D0*z*dzdq - dzdq
        df(4, 2) = 4.D0*xi(1)
        df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
        df(6, 2) = 4.D0*xi(1)*dzdq

      else if ( n_nodes==8 ) then
        !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
        f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
        f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
        f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
        f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
        f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
        f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
        f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
        f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
        df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
        df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
        df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
        df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
        df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
        df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
        df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
        df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
        df(5,1) = -xi(1)*(1.-xi(2));
        df(5,2) = -0.5*(1.-xi(1)*xi(1));
        df(6,1) = 0.5*(1.-xi(2)*xi(2));
        df(6,2) = -(1.+xi(1))*xi(2);
        df(7,1) = -xi(1)*(1.+xi(2));
        df(7,2) = 0.5*(1.-xi(1)*xi(1));
        df(8,1) = -0.5*(1.-xi(2)*xi(2));
        df(8,2) = -(1.-xi(1))*xi(2);
      else if ( n_nodes==9 ) then
        !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
        !     789
        !     456
        !     123
        g1 = -.5D0*xi(1)*(1.D0 - xi(1))
        g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
        g3 = .5D0*xi(1)*(1.D0 + xi(1))
        h1 = -.5D0*xi(2)*(1.D0 - xi(2))
        h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
        h3 = .5D0*xi(2)*(1.D0 + xi(2))
        dg1 = xi(1) - 0.5d0
        dg2 = -2.d0*xi(1)
        dg3 = xi(1) + 0.5d0
        dh1 = xi(2)-0.5d0
        dh2 = -2.d0*xi(2)
        dh3 = xi(2) + 0.5d0
        f(1) = g1*h1
        f(2) = g2*h1
        f(3) = g3*h1
        f(4) = g1*h2
        f(5) = g2*h2
        f(6) = g3*h2
        f(7) = g1*h3
        f(8) = g2*h3
        f(9) = g3*h3
        df(1,1) = dg1*h1
        df(1,2) = g1*dh1
        df(2,1) = dg2*h1
        df(2,2) = g2*dh1
        df(3,1) = dg3*h1
        df(3,2) = g3*dh1
        df(4,1) = dg1*h2
        df(4,2) = g1*dh2
        df(5,1) = dg2*h2
        df(5,2) = g2*dh2
        df(6,1) = dg3*h2
        df(6,2) = g3*dh2
        df(7,1) = dg1*h3
        df(7,2) = g1*dh3
        df(8,1) = dg2*h3
        df(8,2) = g2*dh3
        df(9,1) = dg3*h3
        df(9,2) = g3*dh3
      end if

      end subroutine abq_VUEL_2D_shapefunctions

      subroutine abq_VUEL_invert3d(A,A_inverse,determinant)

      double precision, intent(in) :: A(3,3)
      double precision, intent(out) :: A_inverse(3,3)
      double precision, intent(out) :: determinant

      double precision COFACTOR(3,3)

    !   Compute inverse and determinant of 3x3 matrix

      determinant =   A(1,1)*A(2,2)*A(3,3)
     1   - A(1,1)*A(2,3)*A(3,2)
     2   - A(1,2)*A(2,1)*A(3,3)
     3   + A(1,2)*A(2,3)*A(3,1)
     4   + A(1,3)*A(2,1)*A(3,2)
     5   - A(1,3)*A(2,2)*A(3,1)

      IF (determinant==0.d0) THEN
        write(6,*) ' Error in subroutine abq_UEL_inver3d'
        write(6,*) ' A 3x3 matrix has a zero determinant'
        stop
      endif
      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      A_inverse = transpose(COFACTOR) / determinant


      end subroutine abq_VUEL_invert3d

      subroutine abq_VUEL_facenodes_3D(nelnodes,face,list,nfacenodes)
      implicit none

      integer, intent (in)      :: nelnodes
      integer, intent (in)      :: face
      integer, intent (out)     :: list(*)
      integer, intent (out)     :: nfacenodes

    !
    !        Subroutine to return list of nodes on an element face for standard 3D solid elements
    !

      if (nelnodes == 4) then
        nfacenodes = 3
        if   (face == 1) list(1:3) = [1,2,3]
        if (face == 2) list(1:3) = [1,4,2]
        if (face == 3) list(1:3) = [2,4,3]
        if (face == 4) list(1:3) = [3,4,1]
      else if (nelnodes ==6) then
        nfacenodes = 3
        if (face==1) list(1:3) = [1,2,3]
        if (face==2) list(1:3) = [6,5,4]
        if (face==3) list(1:4) = [1,2,5,4]
        if (face==4) list(1:4) = [2,3,6,5]
        if (face==5) list(1:4) = [4,6,3,1]
        if (face>2) nfacenodes = 4
      else if (nelnodes == 10) then
        nfacenodes = 6
        if   (face == 1) list(1:6) = [1,2,3,5,6,7]
        if (face == 2) list(1:6) = [1,4,2,8,9,5]
        if (face == 3) list(1:6) = [2,4,3,9,10,6]
        if (face == 4) list(1:6) = [3,4,1,10,8,7]
      else if (nelnodes == 8) then
        nfacenodes = 4
        if (face==1) list(1:4) = [1,2,3,4]
        if (face==2) list(1:4) = [5,8,7,6]
        if (face==3) list(1:4) = [1,5,6,2]
        if (face==4) list(1:4) = [2,6,7,3]
        if (face==5) list(1:4) = [3,7,8,4]
        if (face==6) list(1:4) = [4,8,5,1]
      else if (nelnodes ==15) then
        nfacenodes = 6
        if (face==1) list(1:6) = [1,2,3,7,8,9]
        if (face==2) list(1:6) = [6,5,4,11,10,12]
        if (face==3) list(1:8) = [1,2,5,4,7,14,10,13]
        if (face==4) list(1:8) = [2,3,6,5,8,15,11,14]
        if (face==5) list(1:8) = [4,6,3,1,12,15,9,13]
        if (face>2) nfacenodes = 8
      else  if (nelnodes == 20) then
        nfacenodes = 8
        if (face == 1) list(1:8) = [1,2,3,4,9,10,11,12]
        if (face == 2) list(1:8) = [5,8,7,6,16,15,14,13]
        if (face == 3) list(1:8) = [1,5,6,2,17,13,18,9]
        if (face == 4) list(1:8) = [2,6,7,3,18,14,19,10]
        if (face == 5) list(1:8) = [3,7,8,4,19,15,6,11]
        if (face == 6) list(1:8) = [4,8,5,1,20,16,17,12]
      endif

      end subroutine abq_VUEL_facenodes_3D


