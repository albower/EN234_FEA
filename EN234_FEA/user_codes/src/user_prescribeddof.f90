     subroutine user_prescribeddof(node_number,idof,parameter_list,nparam,dofvalue,ignoredof)
     use Types
     ! use Globals         ! Time and time increment can be accessed through this module
      use Mesh, only : node,node_list             ! Mesh and solution data can be accessed through this module
      use Mesh, only : coords
      use Mesh, only : dof_total
      use Staticstepparameters, only : current_step_number
      use Element_Utilities, only : cross_product
      use Mesh, only : dof_increment

      implicit none
     
     integer, intent(in)     :: node_number
     integer, intent(in)     :: idof
     integer, intent(in)     :: nparam
     
     real(prec), intent(in)  :: parameter_list(nparam)
     real(prec), intent(out) :: dofvalue
     
     logical, intent (inout) :: ignoredof          ! Set this flag to .true. to skip enforcing a DOF

     integer :: tensile_dir
     integer :: n_extension_steps
     integer :: n_rotation_steps
     integer :: rotation_axis
     integer :: n_base_nodes
     integer :: n_end_nodes
     integer :: base_fixed_node
     integer :: base_constrained_node
     integer :: base_constrained_dir
     integer :: i

     real (prec) :: extension
     real (prec) :: rotation_angle
     real (prec) :: angle

     real (prec) :: refcoords(3)
     real (prec) :: currentcoords(3)
     real (prec) :: axis(3)
     real (prec) :: newcoords(3)


     ! User subroutine to compute value of a prescribed degree of freedom
     ! Subroutine must return a value for dofvalue
     
     ! Test code to apply a uniaxial tensile stretch followed by a rigid body rotation to a bar
     !

     tensile_dir = int(parameter_list(1))
     extension = parameter_list(2)
     n_extension_steps = int(parameter_list(3))
     rotation_angle = parameter_list(4)
     rotation_axis = int(parameter_list(5))
     n_rotation_steps = int(parameter_list(6))
     n_base_nodes = int(parameter_list(7))
     n_end_nodes = int(parameter_list(8))
     base_fixed_node = int(parameter_list(9))
     base_constrained_node = int(parameter_list(10))
     base_constrained_dir = int(parameter_list(11))

     ignoredof = .true.


     if (current_step_number<=n_extension_steps) then
        do i = 12,12+n_base_nodes-1
           if (node_number==int(parameter_list(i)) .and. idof==tensile_dir) then
              dofvalue = 0.d0
              ignoredof = .false.
           endif
        end do
        do i = 12+n_base_nodes,12+n_base_nodes+n_end_nodes-1
           if (node_number == int(parameter_list(i)) .and. idof==tensile_dir) then
              dofvalue = extension*current_step_number/n_extension_steps
              ignoredof = .false.
           endif
        end do
        if (node_number == base_fixed_node) then
           dofvalue = 0.d0
           ignoredof = .false.
        endif
        if (node_number == base_constrained_node) then
           if (idof==base_constrained_dir) then
               dofvalue = 0.d0
               ignoredof = .false.
           endif
        endif

     else if (current_step_number<=n_extension_steps+n_rotation_steps) then
        refcoords(1:3) = coords(node_list(node_number)%coord_index:node_list(node_number)%coord_index+2)
        currentcoords(1:3) = refcoords(1:3) + dof_total(node_list(node_number)%dof_index:node_list(node_number)%dof_index+2)

        angle = PI_D*rotation_angle/(n_rotation_steps*180.d0)
        axis = 0.d0
        axis(rotation_axis) = 1.d0
        newcoords = dcos(angle)*currentcoords + (1.d0-dcos(angle))*dot_product(axis,currentcoords)*axis &
                                              + dsin(angle)*cross_product(axis,currentcoords)
        dofvalue = newcoords(idof)-refcoords(idof)
        ignoredof = .false.
     else
        dofvalue = 0.d0
        ignoredof = .false.

     endif

     return
     
     
     end subroutine user_prescribeddof
