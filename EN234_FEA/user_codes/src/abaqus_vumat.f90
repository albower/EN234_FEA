!
!    ABAQUS format user material subroutine for explicit dynamics
!
!

     subroutine vumat(nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
       stepTime, totalTime, dt, cmname, coordMp, charLength, &
       props, density, strainInc, relSpinInc, &
       tempOld, stretchOld, defgradOld, fieldOld, &
       stressOld, stateOld, enerInternOld, enerInelasOld, &
       tempNew, stretchNew, defgradNew, fieldNew, &
       stressNew, stateNew, enerInternNew, enerInelasNew )
!
!      include 'vaba_param.inc'
!
      implicit double precision (a-h,o-z)

      dimension props(nprops), density(nblock), coordMp(nblock,*), &
       charLength(nblock), strainInc(nblock,ndir+nshr), &
       relSpinInc(nblock,nshr), tempOld(nblock), &
       stretchOld(nblock,ndir+nshr), &
       defgradOld(nblock,ndir+nshr+nshr), &
       fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr), &
       stateOld(nblock,nstatev), enerInternOld(nblock), &
       enerInelasOld(nblock), tempNew(nblock), &
       stretchNew(nblock,ndir+nshr), &
       defgradNew(nblock,ndir+nshr+nshr), &
       fieldNew(nblock,nfieldv), &
       stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev), &
       enerInternNew(nblock), enerInelasNew(nblock)

       character*80 cmname
!
!      Local variables
!
       integer iblock

       double precision D(6,6)
       double precision E,xnu,d11,d12,d44
!
!      Conventions for storing tensors:
!
!      Deformation gradients are provided as components in the global basis
!      (ABAQUS also allows the user to define a fixed local coord system defined with *ORIENTATION)
!
!      Stresses and stretches are stored as components in a basis that 'rotates with the material'
!      These are defined as follows:
!          Let e_i be the initial basis vectors (the global ijk directions, or user-specified basis vectors)
!          Let R be the rotation tensor, defined through the polar decomposition of deformation gradient F=RU=VR
!          Then define rotated basis vectors m_i = R e_i.
!          The co-rotational components of a tensor are defined as A_ij = m_i . A . m_j
!          The components A_ij can also be interpreted as the global components of the tensor  R^T A R
!
!
!      The components of symmetric tensors (stress,strainInc and stretch) are stored as vectors in the following order
!                      For 3D problems (ndir=3,nshr=3) [11,22,33,12,23,13]
!                      For 2D problems (ndir=3,nshr=1) [11,22,33,12]
!      The components of unsymmetric tensors (defgrad and relSpinInc) are stored as vectors in the following order
!                      For 3D problems (ndir=3,nshr=3) [11,22,33,12,23,31,21,32,13]
!                      For 2D problems (ndir=3,nshr=1) [11,22,33,12,21]
!
!      The stresses are Cauchy (true) stress
!
!      nblock                   No. integration points in data block (data are provided in 'blocks' of integration points)
!                               EN234FEA always has only 1 int pt in each block
!      ndir                     No. direct tensor components
!      nshr                     No. shear tensor components
!      nstatev                  No. user defined state variables (declared in input file)
!      nprops                   No. material properties
!      lanneal                  Annealing flag - if set to 1, then stresses are zeroed, and state vars should be reset to initial state
!      stepTime                 Elapsed time in this step at end of increment
!      totalTime                Total elapsed time at end of time increment
!      dt                       time step
!      cmname                   Material name
!      coordMp(n,i)             ith coord of nth integration point
!      charLength(i)            Characteristic length of element (in EN234FEA this is (element volume)^1/3)
!      props(k)                 kth material property
!      density                  density
!      strainInc(n,i)           ith strain increment component for nth int pt.  Strain components are in a
!                               basis that rotates with material, stored in order [de11, de22, de33, de12, de23, de13]
!                               NOTE THAT SHEAR STRAIN COMPONENTS ARE NOT DOUBLED IN THE VECTOR
!      relSpinInc(n,i)          ith component of relative spin increment tensor.   The global components of relSpinInc are
!                               dW_ij - dR_ij   where dW is the skew part of displacement gradient increment, and
!                               dR_ij is the rotation increment.
!      tempOld                  Temperature at start of increment
!      stretchOld(n,i)          ith component of right stretch V in co-rotational basis (equivalent to U in global basis) at start of increment
!      defgradOld(n,i)          ith component of deformation gradient in fixed global basis at start of increment
!      fieldOld(n,i)            ith user-defined field variable at nth integration point at start of increment
!      stressOld(n,i)           ith component of Cauchy stress in co-rotational basis (equivalent to R^T sigma R in fixed basis)
!      stateOld(n,i)            ith user defined material state variable at start of increment
!      enerInternOld(n)         specific internal energy at nth integration point, start of increment
!      enerInelasOld(n)         specific irreversible dissipation at nth integration point, start of increment
!      tempNew                  Temperature at end of increment
!      stretchNew(n,i)          ith component of right stretch V in co-rotational basis (equivalent to U in global basis) at start of increment
!      defgradNew(n,i)          ith component of deformation gradient in fixed global basis at end of increment
!      fieldNew(n,i)            ith user-defined field variable at nth integration point at end of increment
!      stressNew(n,i)           ith component of Cauchy stress in co-rotational basis (equivale nt to R^T sigma R in fixed basis)
!      stateNew(n,i)            ith user defined material state variable at end of increment
!      enerInternNew(n)         specific internal energy at nth integration point, end of increment
!      enerInelasNew(n)         specific irreversible dissipation at nth integration point, end of increment
!
!       Coded for 3D problems; small stretch assumption (works approximately for finite rotations)
!
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
        D(4,4) = 2.d0*d44    ! Important - shear strain components are not doubled in a VUMAT
        D(5,5) = 2.d0*d44
        D(6,6) = 2.d0*d44

        do iblock = 1,nblock
            stressNew(iblock,1:6) = stressOld(iblock,1:6) + matmul(D,strainInc(iblock,1:6))
        end do


!
End Subroutine vumat
