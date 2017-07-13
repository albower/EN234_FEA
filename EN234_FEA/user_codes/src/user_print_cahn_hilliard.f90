subroutine user_print(n_steps)
  use Types
  use ParamIO
  use Globals, only : TIME, DTIME
  use Mesh
  use Printparameters
  use Staticstepparameters, only: current_step_number
  use User_Subroutine_Storage, only : length_state_variable_array ! Max no. state variables on any element

  implicit none
  
  integer, intent(in) :: n_steps                                 ! Current step number
  
  integer ::  lmn
  integer ::  status
  integer ::  n,k
  integer ::  n_state_vars_per_intpt                                         ! No. state variables per integration point
  real (prec) ::   vol_averaged_strain(6)                                    ! Volume averaged strain in an element
  real (prec), allocatable ::   vol_averaged_state_variables(:)              ! Volume averaged state variables in an element
  real (prec) :: arclen

!   User subroutine for 2017 HW problems
!   The parameter is the HW problem number


  if (n_user_print_parameters==0) then
       write(IOW,*) ' Error in subroutine user_print '
       write(IOW,*) ' At least one parameter must be provided (specifying the HW problem number '
       stop
  endif

  write(6,*) ' User Print '
  write(6,*) ' User print parameters ',user_print_parameters(1)
  if (user_print_parameters(1)==8) then
      print_dof = .true.
      print_displacedmesh = .true.
      displacementscalefactor = 1.d0
      n_field_variables = 4
      zone_dimension(1) = 2
      zone_ndof(1) = 4

      state_print_unit = user_print_units(1)

      allocate(field_variable_names(n_field_variables), stat=status)

      field_variable_names(1) = 'S11'
      field_variable_names(2) = 'S22'
      field_variable_names(3) = 'S33'
      field_variable_names(4) = 'S12'

      combinezones = .true.

      !  Suppress nodes 5-8 on all elements

      do lmn = 1,n_elements
          element_list(lmn)%n_nodes = 4
      end do

      call print_state

      do lmn = 1,n_elements
          element_list(lmn)%n_nodes = 8
      end do

      deallocate(field_variable_names)

  else if (user_print_parameters(1)==9) then
  
      allocate(vol_averaged_state_variables(length_state_variable_array), stat=status)

      if (status/=0) then
          write(IOW,*) ' Error in subroutine user_print'
          write(IOW,*) ' Unable to allocate memory for state variables '
          stop
      endif

      lmn = 1     ! The element number

      call compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_variables,length_state_variable_array, &
          n_state_vars_per_intpt)

      if (TIME<1.d-12) then
          if (n_state_vars_per_intpt<6) then
              write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23'
          else
              write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23'
          endif
      endif

      if (n_state_vars_per_intpt<6) then
          write(user_print_units(1),'(7(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6)
      else
          write(user_print_units(1),'(13(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6),vol_averaged_state_variables(1:6)
      endif

      deallocate(vol_averaged_state_variables)
  else if (user_print_parameters(1)==10)  then
     write(6,*) ' Printing '
     if (current_step_number==1) then
       write(user_print_units(1),'(A)') 'VARIABLES = X,Y,Ux,Uy,theta'
       write(user_print_units(2),'(A)') 'VARIABLES = S,Fx,Fy,Mz'
     endif
     write(user_print_units(1),'(A10,E12.4,A1)') ' ZONE, T="',TIME+DTIME,'"'
     write(user_print_units(2),'(A10,E12.4,A1)') ' ZONE, T="',TIME+DTIME,'"'

     do n = 1,n_nodes
        write(6,'(5(1x,D17.6))') coords(2*n-1),coords(2*n), &
        (dof_total(3*n-k)+dof_increment(3*n-k), k=2,0,-1)
        write(user_print_units(1),'(5(1x,D17.6))') coords(2*n-1),coords(2*n), &
        (dof_total(3*n-k)+dof_increment(3*n-k), k=2,0,-1)
     end do
     arclen = 0.d0
     do n = 1,n_nodes-1
        arclen = arclen + 0.5d0*dsqrt((coords(2*(n+1)-1)-coords(2*n-1))**2. +  (coords(2*(n+1))-coords(2*n))**2.)
        write(user_print_units(2),'(4(1x,D17.6))') arclen, &
        (updated_state_variables(3*n-k), k=2,0,-1)
        arclen = arclen + 0.5d0*dsqrt((coords(2*(n+1)-1)-coords(2*n-1))**2. +  (coords(2*(n+1))-coords(2*n))**2.)
     end do

  else
      write(IOW,*) ' No user prints are defined for HW problem ',user_print_parameters(1)
  endif



end subroutine user_print


subroutine compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_vars,length_output_array, &
                                                                                                       n_state_vars_per_intpt)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent ( in )      :: lmn                                          ! Element number
    integer, intent ( in )      :: length_output_array

    real (prec), intent( out )  ::  vol_averaged_strain(6)
    real (prec), intent( out )  ::  vol_averaged_state_vars(length_output_array)

    integer, intent( out )      :: n_state_vars_per_intpt

    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element


    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i
    integer      :: n_coords, n_dof
    integer      :: iof
    integer      :: status

    real (prec)  ::  el_vol
    real (prec), allocatable  ::  B(:,:)               ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  dstrain(6)                        ! Strain increment vector
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,length_dof_array), stat=status)

    if (status/=0) then
       write(IOW,*) ' Error in subroutine compute_volume_average_3D'
       write(IOW,*) ' Unable to allocate memory for element variables '
       stop
    endif
    !
    ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    vol_averaged_strain = 0.d0
    vol_averaged_state_vars = 0.d0
    el_vol = 0.d0
    n_state_vars_per_intpt = n_state_variables/n_points

    if (n_state_vars_per_intpt>size(vol_averaged_state_vars)) then
       write(IOW,*) ' Error detected in subroutine compute_element_volume_average_3d '
       write(IOW,*) ' The element contains ',n_state_vars_per_intpt
       write(IOW,*) ' but the array storing averaged state variables has length ',size(vol_averaged_state_vars)
       stop
    endif
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)

        iof = n_state_vars_per_intpt*(kint-1)+1
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        strain = matmul(B(1:6,1:3*n_nodes),dof_total(1:3*n_nodes))
        dstrain = matmul(B(1:6,1:3*n_nodes),dof_increment(1:3*n_nodes))

        vol_averaged_strain(1:6) = vol_averaged_strain(1:6) + (strain(1:6)+dstrain(1:6))*w(kint)*determinant

        if (n_state_vars_per_intpt>0) then
           vol_averaged_state_vars(1:n_state_vars_per_intpt) = vol_averaged_state_vars(1:n_state_vars_per_intpt) &
                              + updated_state_variables(iof:iof+n_state_vars_per_intpt-1)*w(kint)*determinant
        endif

        el_vol = el_vol + w(kint)*determinant

    end do

    vol_averaged_strain = vol_averaged_strain/el_vol
    vol_averaged_state_vars = vol_averaged_state_vars/el_vol

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)

    return




end subroutine compute_element_volume_average_3D
