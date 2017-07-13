subroutine compute_static_step
    use Types
    use ParamIO
    use Globals
    use Controlparameters, only : abaqusformat
    use Mesh, only : dof_increment, dof_total, initial_state_variables, updated_state_variables
    use Mesh, only : n_elements, n_nodes
    use Boundaryconditions, only : n_constraints
    use Staticstepparameters
    use Stiffness       

    implicit none

    logical :: fail
    logical :: continue_timesteps
    logical :: converged
    logical :: activatestateprint
    logical :: activateuserprint
    integer :: status
    
    integer :: iteration
    integer :: i
    
    real (prec) :: new_time_increment

    DTIME = timestep_initial
    TIME = 0.d0
    continue_timesteps = .true.
 
    current_step_number = 1
    if (solvertype==1) then                        ! Direct solver

        if (allocated(node_numbers)) deallocate(node_numbers)
        if (allocated(node_order_index)) deallocate(node_order_index)
        allocate(node_numbers(n_nodes+n_constraints), stat = status)
        allocate(node_order_index(n_nodes+n_constraints), stat = status)

        call generate_node_numbers(1,n_elements,.true.,node_numbers,node_order_index)
        call compute_profile
        call allocate_direct_stiffness

        do while (continue_timesteps)

            call assemble_direct_stiffness(fail)

            if (fail) then                          ! Force a timestep cutback if stiffness computation fails
                converged=.false.
                iteration = 0
                call compute_static_time_increment(iteration,converged,continue_timesteps, &
                    activatestateprint,activateuserprint,new_time_increment)
                DTIME = new_time_increment
                dof_increment = 0.d0
                cycle
            endif
            call apply_direct_boundaryconditions
            call solve_direct

            converged = .true.

            do iteration = 1,max_newton_iterations
                call assemble_direct_stiffness(fail)

                if (fail) exit                         ! Force a cutback if stiffness computation fails

                call apply_direct_boundaryconditions
                call convergencecheck(iteration,converged)
                if (converged) exit
                call solve_direct

            end do

            call compute_static_time_increment(iteration,converged,continue_timesteps, &
                activatestateprint,activateuserprint,new_time_increment)
         
            if (.not.converged) then
                dof_increment = 0.d0
                DTIME = new_time_increment
                cycle
            endif

            if (activatestateprint)  call print_state
            if (activateuserprint) call user_print(current_step_number)

            !        Update solution and continue
            current_step_number = current_step_number + 1
            dof_total = dof_total + dof_increment
            initial_state_variables = updated_state_variables
            TIME = TIME + DTIME
            DTIME = new_time_increment
            write(6,*) ' Step ',current_step_number
        end do
      
    else if (solvertype==2) then                        ! Conjugate gradient solver
     
        call initialize_cg_solver

        do while (continue_timesteps)
            call assemble_conjugate_gradient_stiffness(fail)
            if (fail) then                          ! Force a timestep cutback if stiffness computation fails
                converged=.false.
                iteration = 0
                call compute_static_time_increment(iteration,converged,continue_timesteps, &
                    activatestateprint,activateuserprint,new_time_increment)
                dof_increment = 0.d0
                cycle
            endif
            call apply_cojugategradient_boundaryconditions
            call solve_conjugategradient(fail)
            if (fail) then                          ! Force a timestep cutback if CG iteration fails
                converged=.false.
                iteration = 0
                call compute_static_time_increment(iteration,converged,continue_timesteps, &
                    activatestateprint,activateuserprint,new_time_increment)
                DTIME = new_time_increment
                dof_increment = 0.d0
                cycle
            endif
            converged = .true.
         
            if (nonlinear) then                          ! Nonlinear problem - activate Newton iterations
                do iteration = 1,max_newton_iterations
                    call assemble_conjugate_gradient_stiffness(fail)
                    if (fail) then                          ! Force a cutback if stiffness computation fails
                        converged = .false.
                        exit
                    endif
                    call apply_cojugategradient_boundaryconditions
                    call convergencecheck(iteration,converged)
                    if (converged) exit
                    call solve_conjugategradient(fail)
                end do
            endif

            call compute_static_time_increment(iteration,converged,continue_timesteps, &
                activatestateprint,activateuserprint,new_time_increment)
         
            if (.not.converged) then
                DTIME = new_time_increment
                dof_increment = 0.d0
                cycle
            endif

            if (abq_PNEWDT<1.d0) then    ! Timestep cutback forced by ABAQUS UEL
               DTIME = abq_PNEWDT*DTIME
               dof_increment = 0.d0
               cycle
            endif
         
            if (activatestateprint) call print_state
            if (activateuserprint) call user_print(current_step_number)

            !        Update solution and continue
            current_step_number = current_step_number + 1
            dof_total = dof_total + dof_increment
            initial_state_variables = updated_state_variables
            TIME = TIME + DTIME
            DTIME = new_time_increment
        end do
     
    endif

end subroutine compute_static_step

subroutine compute_static_time_increment(iteration,converged,continue_timesteps, &
    activatestateprint,activateuserprint,new_time_increment)
    use Types
    use ParamIO
    use Globals, only : TIME,DTIME
    use Controlparameters, only : abaqusformat
    use Staticstepparameters
    use Printparameters
    implicit none
   
    integer, intent (in)         :: iteration
    logical, intent (in)         :: converged
   
   
    logical, intent (out)        :: continue_timesteps
    logical, intent (out)        :: activatestateprint
    logical, intent (out)        :: activateuserprint
   
    real (prec), intent (out)    :: new_time_increment
   
    continue_timesteps = .true.
    activatestateprint = .false.
    activateuserprint  = .false.
    new_time_increment = DTIME

    if (abaqusformat) then
       if (abq_PNEWDT<1)  then
          new_time_increment = abq_PNEWDT*DTIME

        write(IOW,'(//A)')        ' +++ A timestep cutback was specified in an ABAQUS UEL or UMAT +++ '
        write(IOW,'(A35,G13.5)')  '     Timestep has been reduced to: ',new_time_increment

          return
       endif
    endif

    if (.not.converged) then
        new_time_increment = DTIME/2.D0
        if (new_time_increment<timestep_min) then
            write(IOW,'(A)') ' Time step has been reduced to the minimum allowable value without convergence '
            write(IOW,'(A)') ' Analysis has been terminated '
            write(6,'(A)') ' timestep cut below min value '
            write(6,'(A)') ' Analysis terminated'
            stop
        endif
        write(IOW,'(//A)')        ' +++ Newton-Raphson iterations have not converged +++ '
        write(IOW,'(A35,G13.5)')  '     Timestep has been reduced to: ',new_time_increment
        return
    endif

    if (iteration<max_newton_iterations/5) new_time_increment = 1.25D0*DTIME

    if (abaqusformat.and.abq_PNEWDT<2.d0)  new_time_increment = abq_PNEWDT*DTIME

    if (new_time_increment>timestep_max) new_time_increment = timestep_max
 
    if (TIME+DTIME+new_time_increment >= max_total_time) then
        new_time_increment = max_total_time-TIME-DTIME
        if (new_time_increment<=0.d0) continue_timesteps = .false.
    endif
    
    if (current_step_number == max_no_steps) continue_timesteps = .false.

    if (stateprint) then
        if (current_step_number==1) activatestateprint = .true.
        if (mod(current_step_number,state_print_steps)==0) activatestateprint = .true.
        if (.not.continue_timesteps) activatestateprint = .true.
    endif
    if (userprint) then
        if (current_step_number==1) activateuserprint = .true.
        if (mod(current_step_number,user_print_steps)==0) activateuserprint = .true.
        if (.not.continue_timesteps) activateuserprint = .true.
    endif 
    
    write(IOW,'(//A16,I5,A27)') '    Step number ',current_step_number,' has completed successfully'
    write(IOW,'(A36,G12.5)')    '    Steps remaining:                ',max_no_steps-current_step_number
    write(IOW,'(A36,G12.5)')    '    Elapsed time:                   ',TIME+DTIME
    write(IOW,'(A36,G12.5)')    '    Time remaining:                 ',max_total_time-(TIME+DTIME)
    write(IOW,'(A36,G12.5)')    '    Time step has been adjusted to: ',new_time_increment
    write(IOW,'(A36,G12.5)')    '    Max time step:                  ',timestep_max
    write(IOW,'(A36,G12.5)')    '    Min time step:                  ',timestep_min
   
end subroutine compute_static_time_increment
