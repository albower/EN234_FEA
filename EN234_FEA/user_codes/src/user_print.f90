subroutine user_print(n_steps)
  use Types
  use ParamIO
  use Globals, only : TIME, DTIME
  use Mesh
  use Printparameters
  implicit none
  
  integer, intent(in) :: n_steps
  
  integer :: n, butler_volmer_flag, n1, n2, conc_print_steps
  
  real (prec), save :: total_capacity
!
!  Use this file to process or print time histories of the solution, or to print a non-standard mesh.
!
    conc_print_steps = int(user_print_parameters(1))
!   Write the chemical potential profile and concentration profile to a file
    if (mod(n_steps,conc_print_steps)==0) then
      write(user_print_units(1),'(A)') 'VARIABLES = X,MU,C,DMU,DC'
      write(user_print_units(1),'(A10,E10.4,A1)') ' ZONE, T="',TIME+DTIME,'"'
      butler_volmer_flag = 0  
      if (node_list(1)%n_coords==0) butler_volmer_flag=1   ! If the first node has no coords it is the voltage node
      do n = 1+butler_volmer_flag,n_nodes
         n1 = 2*n-1-butler_volmer_flag
         n2 = 2*n-butler_volmer_flag
         write(user_print_units(1),'(5(1x,D12.5))')   coords(n-butler_volmer_flag),dof_increment(n1)+dof_total(n1),&
                                                    dof_increment(n2)+dof_total(n2), &
                                                    dof_increment(n1),dof_increment(n2)  
      end do
    endif

!   Write the voltage and current history to a file
    if (TIME<1.d-12) then
      write(user_print_units(2),'(A)') 'VARIABLES = TIME,V,I,CAPACITY'
      total_capacity = 0.d0
    ENDIF
    total_capacity = total_capacity - rforce(1)*DTIME
    write(user_print_units(2),'(4(1x,D12.5))') TIME+DTIME,dof_increment(1)+dof_total(1),rforce(1),total_capacity


end subroutine user_print