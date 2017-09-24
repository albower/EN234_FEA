program en234fea
  use Types
  use ParamIO
  use Globals
  use Controlparameters
  implicit none
 
!
!   Homework Assignments 2017
!   You can test the ABAQUS UEL, UMAT and VUMAT codes that you develop for EN234
!   by uncommenting the relevant input files listed below.
!
!   You will also need to write some code for each assignment:
!           A user material for a static analysis should be in a subroutine called UMAT( )
!           A user material for a dynamic analysis should be in a subroutine called VUMAT( )
!           A user element for a static analysis should be in a subroutine called UEL( )
!           A user element for a dynamic analysis should be in a subroutine called VUEL( )
!
!     ABAQUS has a number of other user subroutines, but those are not duplicated in EN234FEA
!

!   Demo codes - these provide examples of coding and testing ABAQUS user elements in EN234FEA
!
!   Small strain linear elasticity - the UEL is in Abaqus_uel_3d.for

!   Use the full path if you are using Intel Parallel studio (you will need to change the path)
!   infil = 'C:/Users/Bower/Source/Repos/EN234_FEA/EN234_FEA/input_files/Abaqus_uel_linear_elastic_3d.in'
!   Eclipse can handle the relative path
!   infil = './input_files/Abaqus_uel_linear_elastic_3d.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = 'C:/Users/Bower/Source/Repos/EN234_FEA/EN234_FEA/input_files/Abaqus_uel_linear_elastic_3d.in'
!   outfil = './Output_files/Abaqus_uel_linear_elastic_3d.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!   Linear elastic plate with a central hole using an ABAQUS UEL
!   infil = './input_files/Abaqus_uel_holeplate_3d.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_holeplate_3d.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!   Simple 1 element demonstration of an ABAQUS VUEL
!   The source code for the user element is in abaqus_vuel.for
!   infil = './input_files/Abaqus_vuel_linear_elastic_3d.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_vuel_linear_elastic_3d.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
!

!   Runs an explicit dynamic simulation of a 3D plate with a central hole with and ABAQUS VUEL
!   This simulation will take a few minutes to run (running in release mode will speed it up)
  
!
!   To run with intel parallel studio you have to put in the full file path
!   infil = 'C:/Users/Bower/Source/Repos/EN234_FEA/EN234_FEA/input_files/Abaqus_vuel_holeplate_3d.in'
!   Eclipse can handle the relative path
!   infil = './input_files/Abaqus_uel_holeplate_3d.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   Fila path for parallel studio
!   outfil = 'C:/Users/Bower/Source/Repos/EN234_FEA/EN234_FEA/Output_files/Abaqus_vuel_holeplate_3d.out'
!   outfil = './output_files/Abaqus_uel_holeplate_3d.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


   
!  Tests an ABAQUS format UMAT subroutine (in abaqus_umat_elastic.for) with 2 8 noded quadrilateral elements
!   infil = './input_files/Abaqus_umat_linear_elastic_3d.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_umat_linear_elastic_3d.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!  Tests the UMAT on the hole in a plate problem.
!   infil = './input_files/Abaqus_umat_holeplate_3d.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_umat_holeplate_3d.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!   infil = './input_files/Abaqus_vumat_linear_elastic_3d.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_vumat_linear_elastic_3d.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!   Homework 3: develop and test an ABAQUS user element implementing 2D linear elasticity with full integration

!   Simple test of a 2D plane element
!   infil = './input_files/Abaqus_uel_linear_elastic_2d.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_linear_elastic_2d.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!  Solve hole-in-a-plate problem with 4 noded quadrilateral elements
!   infil = './input_files/Abaqus_uel_holeplate_2d_quad4.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_holeplate_2d_quad4.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!  Solve hole-in-a-plate problem with 8 noded quads
!   infil = './input_files/Abaqus_uel_holeplate_2d_quad8.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_holeplate_2d_quad8.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!  Solve hole-in-a-plate problem with 3 noded triangles
!   infil = './input_files/Abaqus_uel_holeplate_2d_tri3.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_holeplate_2d_tri3.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!   infil = './input_files/Abaqus_uel_holeplate_2d_tri6.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_holeplate_2d_tri6.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!  HW5  Cantilever beam to test incompatible mode elements

!   infil = './input_files/Abaqus_uel_cantilever.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_cantilever.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!  HW6  Porous elasticity UMAT

!   infil = './input_files/Abaqus_umat_porous_elastic.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_umat_porous_elastic.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

 !  HW7 Hyperelastic user element
!   infil = './input_files/Abaqus_uel_hyperelastic.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_hyperelastic.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!   Hyperelastic umat
!  infil = './input_files/Abaqus_umat_hyperelastic2.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/Abaqus_umat_hyperelastic.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!   HW8 - phase field modeling with elasticity
!   Single element test
!   infil = './input_files/Abaqus_uel_phasefield_1el.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_phasefield_1el.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!   infil = './input_files/Abaqus_uel_phasefield_coarse.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_phasefield_coarse.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!   infil = './input_files/Abaqus_uel_phasefield_fine.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_phasefield_fine.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!   Homework 9 - McCormick model with 1 element
!   infil = './input_files/Abaqus_vumat_McCormick.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_vumat_McCormick.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!   Homework 10 - Continuum beam element solution to end loaded cantilever beam
!   infil = './input_files/Abaqus_uel_continuum_beam.in'
!   open (unit = IOR, file = infil, status = 'old', ERR=500)
!   outfil = './Output_files/Abaqus_uel_continuum_beam.out'
!   open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


  call read_input_file
  
   if (printinitialmesh) call print_initial_mesh

  if (checkstiffness) call check_stiffness(checkstiffness_elementno)
  if (checktangent) call check_tangent(checktangent_materialno)

  if (staticstep) then
      call compute_static_step
      if (checkstiffness) call check_stiffness(checkstiffness_elementno)
      if (checktangent) call check_tangent(checktangent_materialno)
  endif
  
  if (explicitdynamicstep) call explicit_dynamic_step
  
  write(6,*) ' Program completed successfully '

  stop
  
  500 write(6,*) ' Error opening input or output file '
  write(6,*) infil

  
  

  
end program en234fea
