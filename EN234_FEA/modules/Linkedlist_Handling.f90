module Linkedlist_Handling

use Types

!  Data structures and subroutines to store integers in a linked list array of bins
!  These are used to store sparse 2D integer arrays.  

integer, parameter :: databinsize=10

type Integer_Linked_List                  ! Type for a data bin
  sequence
  integer :: n_entries                    ! No. data entries in the bin
  integer :: data_list(databinsize)       ! The data bin
  integer :: next                         ! Index of next bin in the linked list if current bin is full; zero if this is the last bin
end type Integer_Linked_List

contains

   subroutine addlistdata(start_index,linked_list_array,linked_list_end,datavalue)     
      use Types
      use ParamIO
      implicit none
      integer, intent(in)    :: start_index                                            ! The index of the first bin in the linked list
      integer, intent(in)    :: datavalue
      type (Integer_Linked_List), intent(inout) :: linked_list_array(:)
      integer, intent(inout) :: linked_list_end                                        ! The last filled bin in the linked list array
      
      integer :: currentindex
      integer :: ndata

! Add the integer datavalue to a linked list of data bins

      currentindex = start_index
      do while (.true.)
        ndata = linked_list_array(currentindex)%n_entries
        if (ndata<databinsize) then           ! Add to current bin if space remains
           linked_list_array(currentindex)%n_entries=linked_list_array(currentindex)%n_entries+1
           linked_list_array(currentindex)%data_list(linked_list_array(currentindex)%n_entries) = datavalue
           exit
        else
           if (linked_list_array(currentindex)%next>0) then          ! Go to next bin if it has been created
             currentindex = linked_list_array(currentindex)%next
           else                                                      ! Otherwise create a new bin at the end of the linked list
             linked_list_end = linked_list_end + 1
             if (linked_list_end>size(linked_list_array)) then
                 write(IOW,*) ' *** Error detected in subroutine addlistdata ***'
                 write(IOW,*) ' Insufficient memory was found for an integer linked list '
                 write(IOW,*) ' Its current size is ',size(linked_list_array)
                 stop
             endif
             linked_list_array(currentindex)%next = linked_list_end
             currentindex = linked_list_end
             linked_list_array(currentindex)%n_entries=0
             linked_list_array(currentindex)%next = 0
           endif
        endif     
      end do
    end subroutine addlistdata
    
    subroutine adddistinctlistdata(start_index,linked_list_array,linked_list_end,datavalue) 
      use Types
      use ParamIO
      implicit none
      
      integer, intent(in)    :: start_index
      integer, intent(in)    :: datavalue
      type (Integer_Linked_List), intent(inout) :: linked_list_array(:)
      integer, intent(inout) :: linked_list_end
      
      integer :: i
      integer :: currentindex
      integer :: ndata

!     Add datavalue to linked list only if it has not previously been stored
    
      currentindex = start_index
      do while (.true.)
        ndata = linked_list_array(currentindex)%n_entries
        do i = 1,ndata                     ! Check if value is already present
           if (linked_list_array(currentindex)%data_list(i) == datavalue) return
        end do
        if (ndata<databinsize) then           ! Add to current bin if space remains
           linked_list_array(currentindex)%n_entries=linked_list_array(currentindex)%n_entries+1
           linked_list_array(currentindex)%data_list(linked_list_array(currentindex)%n_entries) = datavalue
           exit
        else
           if (linked_list_array(currentindex)%next>0) then          ! Go to next bin if it has been created
             currentindex = linked_list_array(currentindex)%next
           else           ! Create a new bin at the end of the linked list
             linked_list_end = linked_list_end + 1
             if (linked_list_end>size(linked_list_array)) then
                 write(IOW,*) ' *** Error detected in subroutine adddistinctlistdata ***'
                 write(IOW,*) ' Insufficient memory was found for an integer linked list '
                 write(IOW,*) ' Its current size is ',size(linked_list_array)
                 stop
             endif
             linked_list_array(currentindex)%next = linked_list_end
             currentindex = linked_list_end
             linked_list_array(currentindex)%n_entries=0
             linked_list_array(currentindex)%next = 0

           endif
        endif     
      end do
    end subroutine adddistinctlistdata
    
    subroutine extractlinkedlistdata(start_index,linked_list_array,datalist,n_datavalues)
      use Types
      use ParamIO
      implicit none
      integer, intent(in)    :: start_index
      integer, intent(out)   :: datalist(:)
      integer, intent(out)   :: n_datavalues
      type (Integer_Linked_List), intent(in) :: linked_list_array(:)

      integer :: currentindex
      integer :: ndata
      integer :: i

!     Extract all members in the linked list starting at index start_index.

      currentindex = start_index
      n_datavalues = 0
      do while (.true.)
        ndata = linked_list_array(currentindex)%n_entries
        do i = 1,ndata                     ! Extract data from the current bin
          n_datavalues = n_datavalues  + 1
          if (n_datavalues>size(datalist)) then
             write(IOW,'(A)') ' *** Error in subroutine extractlinkedlistdata *** '
             write(IOW,'(A)') ' Insufficient memory to store extracted data '
             write(IOW,*) datalist
             write(IOW,*) ' start index ',start_index
             write(IOW,*) ' Currentindex ',currentindex
             stop
          endif
          datalist(n_datavalues) = linked_list_array(currentindex)%data_list(i)
        end do
        if (linked_list_array(currentindex)%next ==0) exit       ! Exit if this is the last bin
        currentindex = linked_list_array(currentindex)%next      ! Proceed to the next bin
      end do


    end subroutine extractlinkedlistdata
end module
