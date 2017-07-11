module Staticstepparameters

   use Types

   real (prec), save :: timestep_min
   real (prec), save :: timestep_max
   real (prec), save :: timestep_initial
   real (prec), save :: max_total_time
   real (prec), save :: newtonraphson_tolerance
   real (prec), save :: abq_PNEWDT                 ! ABAQUS UEL specified time step
   
   integer, save :: max_no_steps
   integer, save :: solvertype
   integer, save :: max_newton_iterations
   
   integer, save :: current_step_number

   logical, save :: nonlinear  
 

end module
