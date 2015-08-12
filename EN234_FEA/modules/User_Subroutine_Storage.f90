module User_Subroutine_Storage

! Specifies array dimensions required to store user subroutine data
! These variables are computed automatically while reading input file

integer, save :: length_node_array          ! Max nodes on a user element 
integer, save :: length_coord_array         ! max #coordinates * Max # nodes 
integer, save :: length_dof_array           ! max # DOF * max # nodes 


end module
