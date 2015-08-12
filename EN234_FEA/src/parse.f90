!  Routines to parse input file 
  
  !======================SUBROUTINE PARSE =======================
subroutine parse(strin, strpar, length, nstr, lenstr, ityp, iblnk)
    use Types
    implicit none

    integer, intent( in )                    :: length
    integer, intent( out )                   :: nstr
    integer, intent( out )                   :: iblnk
    integer, intent( out )                   :: lenstr(*)
    integer, intent( out )                   :: ityp(*)
    character ( len = 100 ), intent( inout )   :: strin
    character ( len = 100 ), intent( inout )  :: strpar(length)


    ! Local Variables
    integer :: i, itst, itst2, jp, k, n, nd, ne, nnum, np, npm, nz,iasc
    character ( len = 1 ) blank

    blank = ' '
    !  logical :: strcmp

    !     Parser for free format input

    !     STRIN     Character string of length LEN
    !     STRIN should contain characters; integers or real #s separated
    !     by commas.  The routine will separate each field into the parsed
    !     string STRPAR(N).  Blanks are removed, blank lines or lines starting
    !     with * are ignored completely. STRPAR(N) contains strings left justified
    !     and padded on the right with blanks.
    !     STRPAR(N) Nth string extracted from STRIN
    !     LENGTH    Length of input character string
    !     NSTR      No. strings extracted from STRIN
    !     LENSTR(N) Length of Nth string (total no. non blank characters)
    !     ITYP(N)   Type of Nth string. 0 for integer; 1 for real #; 2 for character string
    !     IBLNK     Set to 1 if blank line or comment encountered; zero otherwise.

    !     Fill all output strings with blanks
    do k = 1, length
        do i = 1, 100
            call adchar(strpar(k), 100, ' ', i)
        end do
        lenstr(k) = 0
    end do

    !     Check for blank line or comment line

    iblnk = 1
    do k = 1, length
        if ( ichar(strin(k:k))==ichar('%') ) return            ! First nonblank is a comment
        if ( .not.(ichar(strin(k:k))==ichar(' ')) ) then
            iblnk = 0
            exit
        endif
    end do
    if ( iblnk==1 ) return

    nstr = 1
    !     Loop over input string.
    do i = 1, length
    
        if ( .not.(ichar(strin(i:i))==ichar(' ')) ) then  !     Ignore blanks
            if ( strin(i:i)==',' ) then   !     Check for comma
                nstr = nstr + 1
            else if (ichar(strin(i:i))==ichar('%') ) then  ! Comment
                exit
            else   !     If character is anything other than blank or comma, put it in the output string
                jp = lenstr(nstr) + 1
                call adchar(strpar(nstr), 100, strin(i:i), jp)
                lenstr(nstr) = jp
            end if
        end if
    end do

    !     Determine the type of each output string

    do n = 1, nstr
        ityp(n) = 0
        nd = 0
        ne = 0
        nz = 0
        npm = 0
        np = 0
        nnum = 0
        do i = 1, lenstr(n)
            !     Check for anything other than a number in the string
            itst = iasc(strpar(n), 100, i)
            if ( itst>=48 .and. itst<=57 ) then
                nnum = nnum + 1
              !     Count the number of Ds or ds
            else if ( itst==68 .or. itst==100 ) then
                nd = nd + 1
                !     Check if character following D or d is +or-
                itst2 = iasc(strpar(n), 100, i + 1)
                if ( itst2==43 .or. itst2==45 ) npm = npm - 1
              !     Count the number of Es or es
            else if ( itst==69 .or. itst==101 ) then
                ne = ne + 1
                !     Check if character following E or e is +or-
                itst2 = iasc(strpar(n), 100, i + 1)
                if ( itst2==43 .or. itst2==45 ) npm = npm - 1
              !     Count the number of + or -
            else if ( itst==43 .or. itst==45 ) then
                npm = npm + 1
              !     Count the number of .
            else if ( itst==46 ) then
                np = np + 1
            else
                nz = nz + 1
            end if
        end do
        !     If there was a single decimal point we might have a floating point #
        if ( np==1 ) ityp(n) = 1
        !     Too many .s; +; -; Ee; Dd to be a number
        if ( np>1 .or. nd>1 .or. ne>1 .or. npm>1 ) ityp(n) = 2
        !     Non numerical character in string - must be character string
        if ( nz>0 ) ityp(n) = 2
        !     No numbers in string -- treated as character string
        if ( nnum==0 ) ityp(n) = 2
    end do

end subroutine parse
!===========================SUBROUTINE ADCHAR =========================
subroutine adchar(strin, length, character, ipoin)
    use Types
    implicit none
  
    integer                              :: length
    integer                              :: ipoin
    character ( len = 1 )                :: character
    character ( len = 100 )               :: strin

    !     Add character CHAR to position IPOIN in string STRIN
    strin(ipoin:ipoin) = character

end subroutine adchar
!
!===========================FUNCTION IASC =============================
function iasc(strin, length, ipoin)
    use Types
    implicit none

    integer, intent( in )                  :: length
    integer, intent( in )                  :: ipoin
    character ( len = 100 ), intent( in )   :: strin
    integer                                :: iasc
  
    ! Determine ASCII character # for character at position IPOIN in string STRIN
  
    iasc = ichar(strin(ipoin:ipoin))
  
end function iasc
!==========================FUNCTION STRCMP ====================
function strcmp(wd1, wd2, length)
    use Types
    implicit none

    integer, intent( in )                  :: length
    character ( len = 100 )                :: wd1
    character (len=length)                 :: wd2
    logical                                :: strcmp

    ! Local Variables
    integer :: i, ia, ib
  
    !     Determine whether character strings WD1 and WD2 match
    strcmp = .false.
    do i = 1, length
        ia = ichar(wd1(i:i))
        ib = ichar(wd2(i:i))
        if ( ia/=ib .and. ia + 32/=ib .and. ia/=ib + 32 ) return
    end do
    strcmp = .true.

end function strcmp
