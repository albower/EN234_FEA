module Dynamicstepparameters

   use Types

   real (prec), save :: dynamic_timestep
   real (prec), save :: abq_minStableTimeStep
   real (prec), save :: abq_prevTimeStep
   real (prec), save :: total_dynamic_time
   integer, save :: no_dynamic_steps
   integer, save :: current_dynamic_step
 
 
end module
