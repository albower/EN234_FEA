program en234fea
  use Types
  use ParamIO
  use Globals
  use Controlparameters
  implicit none
 
  character (len=80) :: VS_root_folder
  character (len=80) :: Eclipse_root_folder

  VS_root_folder = 'H:/Repos/EN234_FEA/EN234_FEA/'   ! This should work with Intel Studio on the remote desktop if you follow the instructions for cloning your EN234FEA fork 
  Eclipse_root_folder = './'   !  This should work with Eclipse

  VS_root_folder = 'C:/Users/Bower/Source/Repos/EN234_FEA/EN234_FEA/'
  root_directory = VS_root_folder
    
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
   infil = 'input_files/Abaqus_uel_linear_elastic_3d.in'
   outfil = 'Output_files/Abaqus_uel_linear_elastic_3d.out'

!   Linear elastic plate with a central hole using an ABAQUS UEL
!   infil = 'input_files/Abaqus_uel_holeplate_3d.in'
!   outfil = 'Output_files/Abaqus_uel_holeplate_3d.out'

!   Simple 1 element demonstration of an ABAQUS VUEL
!   The source code for the user element is in abaqus_vuel.for
!   infil = 'input_files/Abaqus_vuel_linear_elastic_3d.in'
!   outfil = 'Output_files/Abaqus_vuel_linear_elastic_3d.out'
!

!   Runs an explicit dynamic simulation of a 3D plate with a central hole with and ABAQUS VUEL
!   This simulation will take a few minutes to run (running in release mode will speed it up)
!
!   infil = 'input_files/Abaqus_uel_holeplate_3d.in'
!   outfil = 'output_files/Abaqus_uel_holeplate_3d.out'
   
!  Tests an ABAQUS format UMAT subroutine (in abaqus_umat_elastic.for) with two 8 noded quadrilateral elements
!   infil = 'input_files/Abaqus_umat_linear_elastic_3d.in'
!   outfil = 'Output_files/Abaqus_umat_linear_elastic_3d.out'

!  Tests the UMAT on the hole in a plate problem.
!   infil = 'input_files/Abaqus_umat_holeplate_3d.in'
!   outfil = 'Output_files/Abaqus_umat_holeplate_3d.out'

!   Tests the VUMAT on a hole in a plate problem
!   infil = 'input_files/Abaqus_vumat_linear_elastic_3d.in'
!   outfil = 'Output_files/Abaqus_vumat_linear_elastic_3d.out'


!   Homework 3: develop and test an ABAQUS user element implementing 2D linear elasticity with full integration

!   Simple test of a 2D plane element
!   infil = 'input_files/Abaqus_uel_linear_elastic_2d.in'
!   outfil = 'Output_files/Abaqus_uel_linear_elastic_2d.out'

!  Solve hole-in-a-plate problem with 4 noded quadrilateral elements
!   infil = 'input_files/Abaqus_uel_holeplate_2d_quad4.in'
!   outfil = 'Output_files/Abaqus_uel_holeplate_2d_quad4.out'

!  Solve hole-in-a-plate problem with 8 noded quads
!   infil = 'input_files/Abaqus_uel_holeplate_2d_quad8.in'
!   outfil = 'Output_files/Abaqus_uel_holeplate_2d_quad8.out'

!  Solve hole-in-a-plate problem with 3 noded triangles
!   infil = 'input_files/Abaqus_uel_holeplate_2d_tri3.in'
!   outfil = 'Output_files/Abaqus_uel_holeplate_2d_tri3.out'

!   infil = 'input_files/Abaqus_uel_holeplate_2d_tri6.in'
!   outfil = 'Output_files/Abaqus_uel_holeplate_2d_tri6.out'

!  HW5  Cantilever beam to test incompatible mode elements

!   infil = 'input_files/Abaqus_uel_cantilever.in'
!   outfil = 'Output_files/Abaqus_uel_cantilever.out'

!  HW6  Porous elasticity UMAT

!   infil = 'input_files/Abaqus_umat_porous_elastic.in'
!   outfil = 'Output_files/Abaqus_umat_porous_elastic.out'

!  HW7 Hyperelastic user element
!   infil = 'input_files/Abaqus_uel_hyperelastic.in'
!   outfil = 'Output_files/Abaqus_uel_hyperelastic.out'

!   Hyperelastic umat
!  infil = 'input_files/Abaqus_umat_hyperelastic2.in'
!  outfil = 'Output_files/Abaqus_umat_hyperelastic.out'

!   HW8 - phase field modeling with elasticity
!   Single element test
!   infil = 'input_files/Abaqus_uel_phasefield_1el.in'
!   outfil = 'Output_files/Abaqus_uel_phasefield_1el.out'

!   infil = 'input_files/Abaqus_uel_phasefield_coarse.in'
!   outfil = 'Output_files/Abaqus_uel_phasefield_coarse.out'

!   infil = 'input_files/Abaqus_uel_phasefield_fine.in'
!   outfil = 'Output_files/Abaqus_uel_phasefield_fine.out'


!   Homework 9 - McCormick model with 1 element
!   infil = 'input_files/Abaqus_vumat_McCormick.in'
!   outfil = 'Output_files/Abaqus_vumat_McCormick.out'



!   Homework 10 - Continuum beam element solution to end loaded cantilever beam
!   infil = 'input_files/Abaqus_uel_continuum_beam.in'
!   outfil = 'Output_files/Abaqus_uel_continuum_beam.out'

   infil = trim(root_directory)//trim(infil)
   outfil = trim(root_directory)//trim(outfil)
   open (unit = IOR, file = trim(infil), status = 'old', ERR=500)
   open (UNIT = IOW, FILE = trim(outfil), STATUS = 'unknown', ERR=500)
   
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
 

  
  

  
end program en234fea
