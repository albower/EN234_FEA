! ParamIO.f90                                                 

module ParamIO
  implicit none

  ! --  Unit to use for log file
  integer, parameter :: IOW = 1
  ! --  Unit to use for intput file
  integer, parameter :: IOR = 2
  ! --  Parameter specifying operating system type.  0 for UNIX, 1 for Windows
  integer, parameter :: IOPSYS = 1
  integer, parameter :: IWT = 1
  
  logical, parameter :: echo = .false.           ! Set to .true. to echo lines of input file
  
  integer :: iblnk                              ! Set to 1 for blank line
  integer :: nstr                               ! No. strings parsed
  character ( len = 100 ) :: strin
  character ( len = 100 ) :: strpar(100)
  integer :: ityp(100), lenstr(100)

  character ( len = 100 ) :: infil, outfil
  character ( len = 100 ) :: root_directory
  

end module ParamIO
