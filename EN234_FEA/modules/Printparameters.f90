 module Printparameters
 
   use Types
 
   integer, save :: state_print_steps                    ! No. steps between state prints
   integer, save :: user_print_steps                     ! No. steps between user prints
   integer, save :: n_field_variables                    ! No. field variables to print (must be less than 18
   integer, save :: n_user_print_files                   ! No. user print files
   integer, save :: n_user_print_parameters              ! No. user print parameters
 
   logical, save :: userprint
   logical, save :: stateprint
   logical, save :: print_dof
   logical, save :: print_displacedmesh
   logical, save :: combinezones
  
   real (prec), save :: displacementscalefactor
 
   character (len=100), save :: state_print_filename
   integer, save :: state_print_unit

   character (len=100), save :: initial_mesh_print_filename
   integer, save :: initial_mesh_print_unit

   
   logical, save, allocatable :: zone_print_flag(:)
   integer, save, allocatable :: user_print_units(:)
   integer, save, allocatable :: zone_dimension(:)
   integer, save, allocatable :: zone_ndof(:)
   character (len=100), save, allocatable :: user_print_filenames(:)
   character (len=100), save, allocatable :: field_variable_names(:)  
   real (prec), save, allocatable :: user_print_parameters(:)

 
 
 
 
 end module
