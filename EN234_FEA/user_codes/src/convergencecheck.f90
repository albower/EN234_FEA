  subroutine convergencecheck(iteration_number,converged)
  use Types
  use ParamIO
  use Mesh, only : nodal_force_norm, unbalanced_force_norm, correction_norm
  use Staticstepparameters, only : newtonraphson_tolerance
  
  integer, intent(in)    :: iteration_number
  logical, intent(out)   :: converged
  
   write (IOW, 99042) iteration_number,correction_norm,nodal_force_norm,unbalanced_force_norm,  &
                   unbalanced_force_norm/(nodal_force_norm+correction_norm), newtonraphson_tolerance

99042 format ( // '  Newton Raphson iteration   ', I5/&
                  '  Correction norm            ', D15.6/      &
                  '  Generalized Force norm     ', D15.6/      &
                  '  Out of balance force norm  ', D15.6/      &
                  '  Ratio                      ', D15.6/     &
                  '  Tolerance                  ', D15.6)
   converged = .false.
   if (nodal_force_norm+correction_norm>0.d0) then
     if (unbalanced_force_norm/(nodal_force_norm+correction_norm) < newtonraphson_tolerance) converged = .true.
   else
     if (unbalanced_force_norm==0.d0) converged = .true.
   endif
  end subroutine convergencecheck
