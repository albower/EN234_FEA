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
   logical, save :: use_lumped_projection_matrix
  
   real (prec), save :: displacementscalefactor
 
   character (len=100), save :: state_print_filename
   integer, save :: state_print_unit

   character (len=100), save :: initial_mesh_print_filename
   integer, save :: initial_mesh_print_unit

   
   logical, save, allocatable :: zone_print_flag(:)                       ! Set to .true. to print a zone
   integer, save, allocatable :: user_print_units(:)                      ! Unit numbers for user print files
   integer, save, allocatable :: zone_dimension(:)                        ! Specifies whether a zone is 2D or 3D
   integer, save, allocatable :: zone_ndof(:)                             ! N DOF for nodes in a zone
   character (len=100), save, allocatable :: user_print_filenames(:)      ! Names of user print files
   character (len=100), save, allocatable :: field_variable_names(:)      ! Names of field variables in a STATE PRINT
   real (prec), save, allocatable :: user_print_parameters(:)             ! List of user print parameters

 
 
 
 
 end module
