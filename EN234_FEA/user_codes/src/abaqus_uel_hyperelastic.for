!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D hyperelastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
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
      double precision  ::  B(9,60)                           ! strain = B*(dof_total)
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      
      double precision  ::  mu, Kmod, G11, G22, G33, G44      ! Material properties
      double precision  ::  Gmat(6,6), F(3,3), Finv(3,3), strsq(3,3)
      double precision  ::  Ctemp(3,3), Cinvtemp(3,3), H(6,9)
      double precision  ::  Y(60,60), Sigsq(3,3), HB(6,60)
      double precision  ::  q(9), C(6), Cbar(6), Cinv(6)
      double precision  ::  P(6), Cstar(6), Cbarstar(6), Sigma(6)
      double precision  ::  Jdet, Qsca, Cdet, quan


      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 8               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex

      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)

      if (MLVARX<3*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 3*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0

      Gmat = 0.d0
      mu = PROPS(1)
      Kmod = PROPS(2)
      G11 = PROPS(3)
      G22 = PROPS(4)
      G33 = PROPS(5)
      G44 = PROPS(6)
      Gmat(1,1) = G11
      Gmat(2,2) = G22
      Gmat(3,3) = G33
      Gmat(4,4) = G44
      Gmat(5,5) = G44
      Gmat(6,6) = G44

      ENERGY(1:8) = 0.d0

      do kint = 1, n_points
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
        B = 0.d0
        B(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
        B(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
        B(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
        B(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
        B(5,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
        B(6,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
        B(7,3:3*NNODE:3)   = dNdx(1:NNODE,1)
        B(8,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
        B(9,3:3*NNODE:3)   = dNdx(1:NNODE,2)
        
        F = 0.d0
        do i = 1, 3
            F(i,1:3) = matmul(U(i:3*NNODE-3+i:3),dNdx(1:NNODE,1:3))
            F(i,i) = F(i,i) + 1
        end do
        call abq_UEL_invert3d(F,Finv,Jdet)
        
        call abq_uel_hyper_generateH(F,H)
        
        Ctemp = matmul(transpose(F),F)
        call abq_uel_hyper_sqtovec(Ctemp,C)
        Cbar(1:6) = C(1:6)/(Jdet**(2/3))
        call abq_UEL_invert3d(Ctemp,Cinvtemp,Cdet)
        call abq_uel_hyper_sqtovec(Cinvtemp,Cinv)
        Cstar(1:6) = C(1:6)
        Cstar(4:6) = Cstar(4:6)*2
        Cbarstar(1:6) = Cstar(1:6)/(Jdet**(2/3))
        
        Qsca = 0.d0
        do i = 1, 3
            Qsca = Qsca + Gmat(i,i)*(Cbarstar(i)-1)*(Cbarstar(i)-1)/4
        end do
        do i = 4, 6
            Qsca = Qsca + Gmat(i,i)*Cbarstar(i)*Cbarstar(i)/4
        end do
        
        quan = 0.d0
        do i = 1, 3
            quan = quan + Gmat(i,i)*Cstar(i)*(Cbarstar(i)-1)
        end do
        do i = 4, 6
            quan = quan + Gmat(i,i)*Cstar(i)*Cbarstar(i)
        end do
        do i = 1, 3
            P(i) = (Gmat(i,i)*(Cbarstar(i)-1) - quan*Cinv(i)/3)
     1             /2/(Jdet**(2/3))
        end do
        do i = 4, 6
            P(i) = (Gmat(i,i)*Cbarstar(i) - quan*Cinv(i)/3)
     1             /2/(Jdet**(2/3))
        end do
        Sigma(1:6) = mu*exp(Qsca)*P(1:6)
     1               +Kmod*Jdet*(Jdet-1)*Cinv(1:6)
        call abq_uel_hyper_q(Sigma,F,q)
        
        call abq_uel_hyper_TanMatrix(mu,Kmod,Qsca,Jdet,
     1           Gmat,Cstar,Cbarstar,Cinv,P,D)
        
        call abq_uel_hyper_Y(NNODE,Sigma,dNdx,Y)
        
        call abq_uel_hyper_vectosq(Sigma,Sigsq)
        strsq = matmul(F,matmul(Sigsq,transpose(F)))/Jdet
        call abq_uel_hyper_sqtovec(strsq,stress)
        
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
     1   - matmul(transpose(B(1:9,1:3*NNODE)),q(1:9))*
     2                                   w(kint)*determinant

        HB = matmul(H,B(1:9,1:3*NNODE))
        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1        + (matmul(transpose(HB(1:6,1:3*NNODE)),
     2        matmul(D,HB(1:6,1:3*NNODE))) + Y(1:3*NNODE,1:3*NNODE))
     3        *w(kint)*determinant

        ENERGY(2) = ENERGY(2) + mu/2*(exp(Qsca)-1)
     1              + K/2*(Jdet-1)*(Jdet-1)*w(kint)*determinant

        if (NSVARS>=n_points*6) then   ! Store stress at each integration point (if space was allocated to do so)
            SVARS(6*kint-5:6*kint) = stress(1:6)
        endif
      end do

      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)

      return

      END SUBROUTINE UEL


      subroutine abq_UEL_3D_integrationpoints(n_points, n_nodes, xi, w)

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
            write(6,*) 'Incorrect # of int pts for tetrahedral element '
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

      end subroutine abq_UEL_3D_integrationpoints


      subroutine abq_UEL_3D_shapefunctions(xi,n_nodes,f,df)

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


      end subroutine abq_UEL_3D_shapefunctions


      subroutine abq_uel_hyper_generateH(F,H)
      
      double precision, intent(in) :: F(3,3)
      double precision, intent(out) :: H(6,9)
      
      H = 0.d0
      H(1,1) = F(1,1)
      H(4,1) = F(1,2)
      H(5,1) = F(1,3)
      
      H(2,2) = F(2,2)
      H(4,2) = F(2,1)
      H(6,2) = F(2,3)
      
      H(3,3) = F(3,3)
      H(5,3) = F(3,1)
      H(6,3) = F(3,2)
      
      H(2,4) = F(1,2)
      H(4,4) = F(1,1)
      H(6,4) = F(1,3)
      
      H(1,5) = F(2,1)
      H(4,5) = F(2,2)
      H(5,5) = F(2,3)
      
      H(3,6) = F(1,3)
      H(5,6) = F(1,1)
      H(6,6) = F(1,2)
      
      H(1,7) = F(3,1)
      H(4,7) = F(3,2)
      H(5,7) = F(3,3)
      
      H(3,8) = F(2,3)
      H(5,8) = F(2,1)
      H(6,8) = F(2,2)
      
      H(2,9) = F(3,2)
      H(4,9) = F(3,1)
      H(6,9) = F(3,3)
      
      end subroutine abq_uel_hyper_generateH


      subroutine abq_UEL_invert3d(A,A_inverse,determinant)

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


      end subroutine abq_UEL_invert3d

      
      subroutine abq_facenodes_3D(nelnodes,face,list,nfacenodes)

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

      end subroutine abq_facenodes_3D
      
      
      subroutine abq_uel_hyper_TanMatrix(mu,Kmod,Qsca,Jdet,
     1           Gmat,Cstar,Cbarstar,Cinv,P,D)
      
      double precision, intent(in)  ::  mu, Kmod, Jdet, Qsca, P(6)
      double precision, intent(in)  ::  Cstar(6), Cbarstar(6), Cinv(6)
      double precision, intent(in)  ::  Gmat(6,6)
      double precision, intent(out)  ::  D(6,6)
      double precision  ::  Omega(6,6), M1(6,6), M2(6,6), M3(6,6)
      double precision  ::  V1(6), V2(6), Iden(6)
      
      Iden = 0.d0
      Iden(1:3) = 1
      
      Omega(1,1) = Cinv(1)*Cinv(1)
      Omega(1,2) = Cinv(4)*Cinv(4)
      Omega(1,3) = Cinv(5)*Cinv(5)
      Omega(1,4) = Cinv(1)*Cinv(4)
      Omega(1,5) = Cinv(1)*Cinv(5)
      Omega(1,6) = Cinv(4)*Cinv(5)
      
      Omega(2,1) = Cinv(4)*Cinv(4)
      Omega(2,2) = Cinv(2)*Cinv(2)
      Omega(2,3) = Cinv(6)*Cinv(6)
      Omega(2,4) = Cinv(4)*Cinv(2)
      Omega(2,5) = Cinv(4)*Cinv(6)
      Omega(2,6) = Cinv(2)*Cinv(6)
      
      Omega(3,1) = Cinv(5)*Cinv(5)
      Omega(3,2) = Cinv(6)*Cinv(6)
      Omega(3,3) = Cinv(3)*Cinv(3)
      Omega(3,4) = Cinv(5)*Cinv(6)
      Omega(3,5) = Cinv(5)*Cinv(3)
      Omega(3,6) = Cinv(6)*Cinv(3)
      
      Omega(4,1) = Cinv(1)*Cinv(4)
      Omega(4,2) = Cinv(4)*Cinv(2)
      Omega(4,3) = Cinv(5)*Cinv(6)
      Omega(4,4) = (Cinv(1)*Cinv(2)+Cinv(4)*Cinv(4))/2
      Omega(4,5) = (Cinv(1)*Cinv(6)+Cinv(5)*Cinv(4))/2
      Omega(4,6) = (Cinv(4)*Cinv(6)+Cinv(5)*Cinv(2))/2
      
      Omega(5,1) = Cinv(1)*Cinv(5)
      Omega(5,2) = Cinv(4)*Cinv(6)
      Omega(5,3) = Cinv(5)*Cinv(3)
      Omega(5,4) = (Cinv(1)*Cinv(6)+Cinv(5)*Cinv(4))/2
      Omega(5,5) = (Cinv(1)*Cinv(3)+Cinv(5)*Cinv(5))/2
      Omega(5,6) = (Cinv(4)*Cinv(3)+Cinv(5)*Cinv(6))/2
      
      Omega(6,1) = Cinv(4)*Cinv(5)
      Omega(6,2) = Cinv(2)*Cinv(6)
      Omega(6,3) = Cinv(6)*Cinv(3)
      Omega(6,4) = (Cinv(4)*Cinv(6)+Cinv(5)*Cinv(2))/2
      Omega(6,5) = (Cinv(4)*Cinv(3)+Cinv(5)*Cinv(6))/2
      Omega(6,6) = (Cinv(2)*Cinv(3)+Cinv(6)*Cinv(6))/2
      
      Omega = -Omega
      
      M1 = spread(Cstar,dim=2,ncopies=6)*spread(Cinv,dim=1,ncopies=6)
      M2 = spread(Cinv,dim=2,ncopies=6)*
     1     spread(matmul(Gmat,Cstar),dim=1,ncopies=6)
      M3 = Gmat - (matmul(Gmat,M1)+M2)/3
      M1 = spread(Cinv,dim=2,ncopies=6)*spread(Cinv,dim=1,ncopies=6)
      M3 = M3 + dot_product(Cstar,matmul(Gmat,Cstar))/9*M1
      M3 = M3/(Jdet**(4/3))
      V1 = P - Cinv/3
      M1 = spread(P,dim=2,ncopies=6)*spread(V1,dim=1,ncopies=6)
      M3 = M3 + 2*M1
      V1 = matmul(Gmat,Cbarstar-Iden)
      M1 = spread(Cinv,dim=2,ncopies=6)*spread(V1,dim=1,ncopies=6)
      M2 = dot_product(Cstar,V1)*Omega
      M3 = M3 - (M1+M2)/3/(Jdet**(2/3))
      M3 = M3*mu*exp(Qsca)
      M1 = spread(Cinv,dim=2,ncopies=6)*spread(Cinv,dim=1,ncopies=6)
      D = M3 + Kmod*Jdet*((2*Jdet-1)*M1+2*(Jdet-1)*Omega)
      
      end subroutine abq_uel_hyper_TanMatrix
      
      
      subroutine abq_uel_hyper_q(Sigma,F,q)
      double precision, intent(in)  ::  Sigma(6), F(3,3)
      double precision, intent(out)  ::  q(9)
      double precision  ::  Sigsq(3,3), qsq(3,3)
      
      call abq_uel_hyper_vectosq(Sigma,Sigsq)
      
      qsq = matmul(Sigsq,transpose(F))
      
      q(1) = qsq(1,1)
      q(2) = qsq(2,2)
      q(3) = qsq(3,3)
      q(4) = qsq(2,1)
      q(5) = qsq(1,2)
      q(6) = qsq(3,1)
      q(7) = qsq(1,3)
      q(8) = qsq(3,2)
      q(9) = qsq(2,3)
      
      end subroutine abq_uel_hyper_q
      
      
      subroutine abq_uel_hyper_Y(nn,Sigma,dNdx,Y)
      
      integer, intent(in) :: nn
      double precision, intent(in)  ::  Sigma(6), dNdx(20,3)
      double precision, intent(out)  ::  Y(60,60)
      integer :: i, j
      double precision  ::  Sigsq(3,3), kij
      
      call abq_uel_hyper_vectosq(Sigma,Sigsq)
      
      do i = 1, nn
          do j = 1, nn
              kij = dot_product(dNdx(i,1:3),matmul(Sigsq,dNdx(j,1:3)))  
              Y(3*i-2,3*j-2) = kij
              Y(3*i-1,3*j-1) = kij
              Y(3*i,3*j) = kij
          end do
      end do
      
      end subroutine abq_uel_hyper_Y
      
      
      subroutine abq_uel_hyper_sqtovec(A,V)
      
      double precision, intent(in)  ::  A(3,3)
      double precision, intent(out)  ::  V(6)
      
      V(1) = A(1,1)
      V(2) = A(2,2)
      V(3) = A(3,3)
      V(4) = A(1,2)
      V(5) = A(1,3)
      V(6) = A(2,3)
      
      end subroutine abq_uel_hyper_sqtovec
      
      
      subroutine abq_uel_hyper_vectosq(V,A)
      
      double precision, intent(in)  ::  V(6)
      double precision, intent(out)  ::  A(3,3)
      
      A(1,1) = V(1)
      A(1,2) = V(4)
      A(1,3) = V(5)
      A(2,1) = V(4)
      A(2,2) = V(2)
      A(2,3) = V(6)
      A(3,1) = V(5)
      A(3,2) = V(6)
      A(3,3) = V(3)
      
      end subroutine abq_uel_hyper_vectosq