     subroutine user_prescribedforce(node_number,idof,parameter_list,nparam,force_value)
     use Types
     ! use Globals          ! Time and time increment can be accessed through this module
     ! use Mesh             ! Mesh and solution data can be accessed through this module
     
     integer, intent(in)     :: node_number
     integer, intent(in)     :: idof
     integer, intent(in)     :: nparam
     
     real(prec), intent(in)  :: parameter_list(nparam)
     real(prec), intent(out) :: force_value
     
     ! User subroutine to compute value of a prescribed degree of freedom
     ! Subroutine must return a value for force_value
     
    force_value = 0.d0
     
     
     end subroutine user_prescribedforce