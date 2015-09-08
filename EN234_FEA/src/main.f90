program en234fea
  use Types
  use ParamIO
  use Globals
  use Controlparameters

  implicit none

!  Demo codes - basic 3D linear elasticity
!
!  infil = './input_files/linear_elastic_3d.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/linear_elastic_3d.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

  infil = './input_files/Holeplate_3d.in'
  open (unit = IOR, file = infil, status = 'old', ERR=500)
  outfil = './Output_files/Holeplate_3d.out'
  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!  infil = './input_files/linear_elastic_3d_dynamic.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/linear_elastic_3d_dynamic.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

! This simulation takes a few minutes - be patient!
!  infil = './input_files/holeplate_3d_dynamic.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/holeplate_3d_dynamic.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!  None of the files below will work until you write the codes that will use them!
!
!  Homework 3
!  Basic 2 element test (one or two elements)
!  infil = './input_files/linear_elastic_2d.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/linear_elastic_2d.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
! Homework 3, Basic 2D linear elasticity with different element types.
!  infil = './input_files/holeplate_2d_tri3.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/holeplate_2d_tri3.out'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  infil = './input_files/holeplate_2d_tri6.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/holeplate_2d_tri6.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
!  infil = './input_files/holeplate_2d_quad4.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/holeplate_2d_quad4.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
!  infil = './input_files/holeplate_2d_quad8.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/holeplate_2d_quad8.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
!

!
!  Homework 4, crack tip elements and the J integral
!  infil = './input_files/crack_tri6.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/crack_tri6.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

!  Homework 5, small-strain B bar element - test with same files as in HW3, but
!  try approaching incompressible limit by making Poisson's ratio close to 0.5

!
! Homework 6: small-strain Armstrong-Frederick kinematic hardening model
!  infil = './input_files/cyclic_plastic_3d.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/cyclic_plastic_3d.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

! Homework 7, stretch a hyperelastic bar, check stiffness.
!  infil = './input_files/Hyperelastic_bar_stretch.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/hyperelastic_bar_stretch.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
!!
!!  Homework 7, stretch and rotate a hyperelastic bar
!  infil = './input_files/Hyperelastic_stretch_rotate.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/hyperelastic_stretch_rotate.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
!
! Homework 7, stretch a hyperelastic plate with a central hole
!  infil = './input_files/Holeplate_hyperelastic.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/Holeplate_hyperelastic.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!!  Homework 8, solve the 2D Cahn-Hilliard equation
!  infil = './input_files/cahn_hilliard_2d_fine.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/cahn_hilliard_2d_fine.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!!  Homework 9, Dynamic fracture with explicit dynamics, finite strain Gurson model.
!  infil = './input_files/notch_fracture_dynamic.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/notch_fracture_dynamic.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

  call read_input_file
  
   if (printinitialmesh) call print_initial_mesh

  if (checkstiffness) call check_stiffness(checkstiffness_elementno)

  if (staticstep) then
      call compute_static_step
      if (checkstiffness) call check_stiffness(checkstiffness_elementno)
  endif
  
  if (explicitdynamicstep) call explicit_dynamic_step
  
  write(6,*) ' Program completed successfully '

  stop
  
  500 write(6,*) ' Error opening input or output file '
  

  
  

  
end program en234fea
