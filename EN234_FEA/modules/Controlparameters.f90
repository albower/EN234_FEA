module Controlparameters

   use Types

   logical, save :: printinitialmesh                ! Activates an initial mesh print
   logical, save :: checkstiffness                  ! Activates a check to ensure residual is consistent with element stiffness
   logical, save :: staticstep                      ! Activates a static analysis
   logical, save :: explicitdynamicstep             ! Activates an explicit dynamic analysis
   logical, save :: userstep                        ! Activates a user-controlled time stepping analysis
    
   integer, save :: checkstiffness_elementno        ! Specifies the element flag for stiffness check
  
  
   
end module
