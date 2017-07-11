
subroutine traction_boundarycondition_static(flag,ndims,ndof,nfacenodes,element_coords,length_coord_array,&
                             element_dof_increment,element_dof_total,length_dof_array,traction,ntract,&
                             element_stiffness,element_residual)               ! Output variables
   use Types
   use ParamIO
   use Element_Utilities, only : N1 => shape_functions_1D
   use Element_Utilities, only : dNdxi1 => shape_function_derivatives_1D
   use Element_Utilities, only : dNdx1 => shape_function_spatial_derivatives_1D
   use Element_Utilities, only : xi1 => integrationpoints_1D, w1 => integrationweights_1D
   use Element_Utilities, only : N2 => shape_functions_2D
   use Element_Utilities, only : dNdxi2 => shape_function_derivatives_2D
   use Element_Utilities, only : dNdx2 => shape_function_spatial_derivatives_2D
   use Element_Utilities, only : xi2 => integrationpoints_2D, w2 => integrationweights_2D
   use Element_Utilities, only : initialize_integration_points
   use Element_Utilities, only : calculate_shapefunctions
   implicit none

   integer, intent( in )         :: flag                                      ! Flag specifying traction type. 1=direct value; 2=history+direct; 3=normal+history
   integer, intent( in )         :: ndims                                     ! No. coords for nodes
   integer, intent( in )         :: ndof                                      ! No. DOFs for nodes
   integer, intent( in )         :: nfacenodes                                ! No. nodes on face
   integer, intent( in )         :: length_coord_array                        ! Total # coords
   integer, intent( in )         :: length_dof_array                          ! Total # DOF
   integer, intent( in )         :: ntract                                    ! # traction components
   real( prec ), intent( in )    :: element_coords(length_coord_array)        ! List of coords
   real( prec ), intent( in )    :: element_dof_total(length_dof_array)       ! List of DOFs (not currently used, provided for extension to finite deformations)
   real( prec ), intent( in )    :: element_dof_increment(length_dof_array)   ! List of DOF increments (not used)
   real( prec ), intent( in )    :: traction(ntract)                          ! Traction value

   real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
   real( prec ), intent( out )   :: element_residual(length_dof_array)      ! Element residual force (ROW)
  
  ! Local Variables

   integer :: n_points,k,kint
   real (prec) :: dxdxi1(2),dxdxi2(3,2), x1(2,length_coord_array/2), x2(3,length_coord_array/3), norm(3)
   real (prec) :: determinant

   ! Routine to compute contribution to element residual from distributed force applied to a solid element face
   


   element_residual = 0.D0
   element_stiffness = 0.d0

   if (ndims==2) then                               ! 2D elements
     x1 = reshape(element_coords,(/2,length_coord_array/2/)) 
     n_points = nfacenodes
     call initialize_integration_points(n_points, nfacenodes, xi1, w1)

     do kint = 1, n_points
       call calculate_shapefunctions(xi1(kint:kint,1),nfacenodes,N1,dNdxi1)      
       dxdxi1 = matmul(x1(1:2,1:nfacenodes),dNdxi1(1:nfacenodes,1))
       determinant = dsqrt(dot_product(dxdxi1,dxdxi1))
       
       if (flag<3) then  ! Traction given directly
          do k = 1,ndof
            element_residual(k:ndof*(nfacenodes-1)+k:ndof) = &
                              element_residual(k:ndof*(nfacenodes-1)+k:ndof) + N1(1:nfacenodes)*traction(k)*w1(kint)*determinant
          end do
        else      ! Applied normal to element face
          element_residual(1:2*nfacenodes-1:2) = &
                              element_residual(1:2*nfacenodes-1:2) + N1(1:nfacenodes)*traction(1)*dxdxi1(2)*w1(kint)
          element_residual(2:2*nfacenodes:2) = element_residual(2:2*nfacenodes:2) - N1(1:nfacenodes)*traction(2)*dxdxi1(1)*w1(kint)
        endif
     end do
   else if (ndims==3) then                               ! 3D elements
     x2 = reshape(element_coords,(/3,length_coord_array/3/)) 
     if (nfacenodes == 3) n_points = 3
     if (nfacenodes == 6) n_points = 4
     if (nfacenodes == 4) n_points = 4
     if (nfacenodes == 8) n_points = 9
     call initialize_integration_points(n_points, nfacenodes, xi2, w2)   
     
     do kint = 1,n_points
       call calculate_shapefunctions(xi2(1:2,kint),nfacenodes,N2,dNdxi2)
       dxdxi2 = matmul(x2(1:3,1:nfacenodes),dNdxi2(1:nfacenodes,1:2))
       norm(1) = (dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
       norm(2) = (dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
       norm(3) = (dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))
       determinant = dsqrt(dot_product(norm,norm))
       if (flag<3) then  ! Traction given directly
          do k = 1,ndof
            element_residual(k:ndof*(nfacenodes-1)+k:ndof) = &
                         element_residual(k:ndof*(nfacenodes-1)+k:ndof) + N2(1:nfacenodes)*traction(k)*w2(kint)*determinant   ! Traction given directly
          end do
        else      ! Applied normal to element face
          element_residual(1:3*nfacenodes-2:3) = &
                          element_residual(1:3*nfacenodes-2:3) - N2(1:nfacenodes)*traction(1)*norm(1)*w2(kint)      ! Note determinant is already in normal
          element_residual(2:3*nfacenodes-1:3) = &
                          element_residual(2:3*nfacenodes-1:3) - N2(1:nfacenodes)*traction(1)*norm(2)*w2(kint)
          element_residual(3:3*nfacenodes:3) = &
                           element_residual(3:3*nfacenodes:3) - N2(1:nfacenodes)*traction(1)*norm(3)*w2(kint)
        endif
     end do
   endif

   return
end subroutine traction_boundarycondition_static

subroutine traction_boundarycondition_dynamic(flag,ndims,ndof,nfacenodes,element_coords,length_coord_array,&
                             element_dof_increment,element_dof_total,length_dof_array,traction,ntract,&
                             element_residual)               ! Output variables
   use Types
   use ParamIO
   use Element_Utilities, only : N1 => shape_functions_1D
   use Element_Utilities, only : dNdxi1 => shape_function_derivatives_1D
   use Element_Utilities, only : dNdx1 => shape_function_spatial_derivatives_1D
   use Element_Utilities, only : xi1 => integrationpoints_1D, w1 => integrationweights_1D
   use Element_Utilities, only : N2 => shape_functions_2D
   use Element_Utilities, only : dNdxi2 => shape_function_derivatives_2D
   use Element_Utilities, only : dNdx2 => shape_function_spatial_derivatives_2D
   use Element_Utilities, only : xi2 => integrationpoints_2D, w2 => integrationweights_2D
   use Element_Utilities, only : initialize_integration_points
   use Element_Utilities, only : calculate_shapefunctions
   implicit none

   integer, intent( in )         :: flag                     ! Flag specifying traction type. 1=direct value; 2=history+direct; 3=normal+history                        
   integer, intent( in )         :: ndims                    ! No. coords for nodes
   integer, intent( in )         :: ndof                     ! No. DOFs for nodes
   integer, intent( in )         :: nfacenodes               ! No. nodes on face
   integer, intent( in )         :: length_coord_array       ! Total # coords
   integer, intent( in )         :: length_dof_array         ! Total # DOF
   integer, intent( in )         :: ntract                   ! # traction components
   real( prec ), intent( in )    :: element_coords(length_coord_array)        ! List of coords
   real( prec ), intent( in )    :: element_dof_total(length_dof_array)     ! List of DOFs (not currently used, provided for extension to finite deformations)
   real( prec ), intent( in )    :: element_dof_increment(length_dof_array) ! List of DOF increments (not used)
   real( prec ), intent( in )    :: traction(ntract)              ! Traction value

   real( prec ), intent( out )   :: element_residual(length_dof_array)      ! Element residual force (ROW)
  
  ! Local Variables

   integer :: n_points,k,kint
   real (prec) :: dxdxi1(2),dxdxi2(3,2), x1(2,length_coord_array/2), x2(3,length_coord_array/3), norm(3)
   real (prec) :: determinant

   ! Routine to compute contribution to element residual from distributed force applied to a solid element face
   
   element_residual = 0.D0


   if (ndims==2) then                               ! 2D elements
     x1 = reshape(element_coords,(/2,length_coord_array/2/))
     n_points = nfacenodes
     call initialize_integration_points(n_points, nfacenodes, xi1, w1)

     do kint = 1, n_points
       call calculate_shapefunctions(xi1(kint:kint,1),nfacenodes,N1,dNdxi1)      
       dxdxi1 = matmul(x1(1:2,1:nfacenodes),dNdxi1(1:nfacenodes,1))
       determinant = dsqrt(dot_product(dxdxi1,dxdxi1))
       
       if (flag<3) then  ! Traction given directly
          do k = 1,ndof
            element_residual(k:ndof*(nfacenodes-1)+k:ndof) = &
                        element_residual(k:ndof*(nfacenodes-1)+k:ndof) + N1(1:nfacenodes)*traction(k)*w1(kint)*determinant
          end do
        else      ! Applied normal to element face
          element_residual(1:2*nfacenodes-1:2) = &
                        element_residual(1:2*nfacenodes-1:2) + N1(1:nfacenodes)*traction(1)*dxdxi1(2)*w1(kint)
          element_residual(2:2*nfacenodes:2) = element_residual(2:2*nfacenodes:2) - N1(1:nfacenodes)*traction(2)*dxdxi1(1)*w1(kint)
        endif
     end do
   else if (ndims==3) then                               ! 3D elements
     x2 = reshape(element_coords,(/3,length_coord_array/3/))
     if (nfacenodes == 3) n_points = 3
     if (nfacenodes == 6) n_points = 4
     if (nfacenodes == 4) n_points = 4
     if (nfacenodes == 8) n_points = 9

     call initialize_integration_points(n_points, nfacenodes, xi2, w2)   
     do kint = 1,n_points

       call calculate_shapefunctions(xi2(1:2,kint),nfacenodes,N2,dNdxi2)
       dxdxi2 = matmul(x2(1:3,1:nfacenodes),dNdxi2(1:nfacenodes,1:2))
       norm(1) = (dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
       norm(2) = (dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
       norm(3) = (dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))
       determinant = dsqrt(dot_product(norm,norm))
       if (flag<3) then  ! Traction given directly
          do k = 1,ndof
            element_residual(k:ndof*(nfacenodes-1)+k:ndof) = &
                      element_residual(k:ndof*(nfacenodes-1)+k:ndof) + N2(1:nfacenodes)*traction(k)*w2(kint)*determinant   ! Traction given directly
          end do
        else      ! Applied normal to element face
          element_residual(1:3*nfacenodes-2:3) = &
                          element_residual(1:3*nfacenodes-2:3) + N2(1:nfacenodes)*traction(1)*norm(1)*w2(kint)      ! Note determinant is already in normal
          element_residual(2:3*nfacenodes-1:3) = &
                           element_residual(2:3*nfacenodes-1:3) - N2(1:nfacenodes)*traction(2)*norm(2)*w2(kint)
          element_residual(3:3*nfacenodes:3) = &
                            element_residual(3:3*nfacenodes:3) - N2(1:nfacenodes)*traction(3)*norm(3)*w2(kint)
        endif
     end do
   endif

   return
end subroutine traction_boundarycondition_dynamic
