     subroutine user_initialdof(parameter_list,n_parameters,node_number,coords,n_coords,dof_total,dof_increment,n_dof)
     use Types
     implicit none
     
     integer, intent(in)     :: n_parameters
     integer, intent(in)     :: n_coords
     integer, intent(in)     :: n_dof
     integer, intent(in)     :: node_number
     
     real(prec), intent(in)  :: parameter_list(n_parameters)
     real(prec), intent(in)  :: coords(n_coords)
!
     real(prec), intent(out) :: dof_total(n_dof)
     real(prec), intent(out) :: dof_increment(n_dof)

     real(prec) :: A, k

     ! User subroutine to set initial value of a DOF
     ! Subroutine must return values for dof_total and dof_increment
     
     dof_total = 0.d0
     dof_increment = 0.d0
     
     A = parameter_list(1)
     k = parameter_list(2)

     dof_total(2) = A*sin(k*coords(1))*sin(k*coords(2))
     
     
     end subroutine user_initialdof
