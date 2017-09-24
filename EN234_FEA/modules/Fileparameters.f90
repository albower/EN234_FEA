 module Fileparameters
  use Types
  
  integer, save :: n_total_files
  integer, save, allocatable :: file_units(:)
  character ( len = 100 ), save, allocatable :: file_list(:)
  
  contains
  
  subroutine files(fname,rw, iun)
  use Types
  use ParamIO
  implicit none

  integer, intent( in )                   :: rw
  integer, intent( out )                  :: iun
  character ( len = 100 ), intent( in )   :: fname

  ! Local Variables
  integer :: k
  logical :: strcmp
  logical :: ex
  character (len=220) :: callstr
  character (len=200) :: fullpath
    
  !     Look for file name in list of open files
  
  do k = 1, n_total_files
    if ( strcmp(fname,file_list(k), 100) ) then
      iun = file_units(k)
      return
    end if
  end do

  !     rw determines how the file will be opened
  !     rw = 1 : File is opened for writing in ASCII
  !     rw = -1 : File is opened for reading in ASCII
  !     (after reading, the file must be closed by
  !     the reading subroutine)
  
  !     File has not been opened - store the file name and open it.
  
  !     Move existing file with name
  !     FNAME to FNAME~

  fullpath = trim(root_directory)//trim(fname)
  inquire(FILE=fullpath,EXIST=ex)
  if (rw==1 .and. ex  ) then
    if ( IOPSYS==0 ) then
      callstr = 'mv ' // trim(fullpath) // ' ' //trim(fullpath)//'~'
    else if ( IOPSYS==1 ) then
       callstr = 'move ' // trim(fullpath) // ' ' //trim(fullpath)//'~'
    end if
    call system(trim(callstr))
  endif
  !     Open the file and add it to the list of open files

  !     Set unit number.  The largest open unit number is always in file_units(n_total_files)

  if ( n_total_files==0 ) then
    iun = 11
  else
    iun = file_units(n_total_files) + 1
  end if

  n_total_files = n_total_files + 1
    
  file_units(n_total_files) = iun
  file_list(n_total_files) = trim(fname)
  if (rw==1) then
    continue
    open (UNIT = iun, FILE = trim(fullpath), STATUS = 'NEW')   ! File opened for writing
  else if (rw==-1) then
    open (UNIT = iun, file = trim(fullpath), STATUS = 'OLD')   ! File opened for reading
  endif

end subroutine files

  
 end module
