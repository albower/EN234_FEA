program en234fea
  use Types
  use ParamIO
  use Globals
  use Controlparameters
  use Element_Utilities, only : eigenvecs33
  use Element_Utilities, only : sqrtM33
  implicit none

  real (prec) :: matrix(3,3)
  real (prec) :: eigenvalues(3)
  real (prec) :: eigenvectors(3,3)

!  Uncomment the lines below to get input and output file names from the console
!  This does not work well in eclipse (and does not work at all when debugging in eclipse)
!  100 write (6, 99001, advance = 'no')
!  read (5, '(A100)') infil
!  open (unit = IOR, file = infil, status = 'old', ERR = 100)
!  200 write (6, 99002, ADVANCE = 'no')
!  read (5, '(A100)') outfil
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR = 200)
!  99001 format ('Please specify the input file name:  ')
!  99002 format ('Please specify the log file name:    ')
!

   matrix = reshape((/1.5d0,0.5d0,0.d0,0.5d0,1.5d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))
   eigenvectors = sqrtM33(matrix)

   write(6,*) eigenvectors(1,1:3)
   write(6,*) eigenvectors(2,1:3)
   write(6,*) eigenvectors(3,1:3)


!  infil = './input_files/hyperelastic_3d.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/hyperelastic_3d.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

  stop

  call read_input_file
  
  if (printinitialmesh) call print_initial_mesh
  
  if (checkstiffness) call check_stiffness(checkstiffness_elementno)

  if (staticstep) call compute_static_step
  
  if (explicitdynamicstep) call explicit_dynamic_step
  
  write(6,*) ' Program completed successfully '

  stop
  
  500 write(6,*) ' Error opening input or output file '
  

  
  

  
end program en234fea
