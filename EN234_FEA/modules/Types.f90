! Types.f90
!
! Module to provide default machine precision types, common
! mathematical constants and some enumerated logical types.
!
! Usage:
!          data_type( [KIND =] < precision > ) :: variable_name
!
! where < precision > is one of the precisions listed below.
!
! Precisions:
!   
!   Integers:
!   ---------  
!     huge  : 4-byte integer
!     long  : 2-byte integer
!     short : 1-byte integer
!
!   Reals:
!   ------
!     single : single precision real
!     double : double precision real
!     
!   Complex:
!   --------
!     single_complex : single precision complex
!     double_complex : double precision complex
!
!   Logical:
!   --------
!     default_logical : default logical type
!
! The following constants are also supplied:
!
!   Single Precision : PI, HALFPI, TWOPI, SQRT2
!   Double Precision : PI_D, HALFPI_D, TWOPI_D, SQRT2_D 
!

module Types
  implicit none

  ! Integers
  integer, parameter :: huge = selected_int_kind( 9 )
  integer, parameter :: long = selected_int_kind( 4 )
  integer, parameter :: short = selected_int_kind( 2 )

  ! Single and Double Precision Reals
  integer, parameter :: single = kind( 1.0 ) 
  integer, parameter :: double = kind( 1.D0 )

  ! Single and Double Precision Complex
  integer, parameter :: single_complex = kind( ( 1.0, 1.0 ) ) 
  integer, parameter :: double_complex = kind( ( 1.D0, 1.D0 ) )

  ! Logical 
  integer, parameter :: default_logical = kind( .true. )

  ! Frequently Used Mathematical Constants
  real( single ), parameter :: PI = &
       3.141592653589793238462643383279502884197_single
  real( single ), parameter :: HALFPI = PI / 2.0
  real( single ), parameter :: TWOPI = 2.0 * PI
  real( single ), parameter :: SQRT2 = &
       1.41421356237309504880168872420969807856967_single
  real( double ), parameter :: PI_D = &
       3.141592653589793238462643383279502884197_double
  real( double ), parameter :: HALFPI_D = PI_D / 2.D0
  real( double ), parameter :: TWOPI_D = 2.D0 * PI_D
  real( double ), parameter :: SQRT2_D = &
       1.41421356237309504880168872420969807856967_double

  ! Enumerated Logical Data Types
  logical( default_logical ), parameter :: yes = .true.
  logical( default_logical ), parameter :: no = .false.
  logical( default_logical ), parameter :: successful = .true.
  logical( default_logical ), parameter :: failed = .false.

  ! Useful matrices
  real ( single ), dimension(2,2), parameter :: eye2 = reshape((/1.0,0.0,0.0,1.0/),(/2,2/))
  real ( single ), dimension(3,3), parameter :: eye3 = reshape((/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/),(/3,3/))


  real ( double ), dimension(2,2), parameter :: eye2_D = reshape((/1.d0,0.d0,0.d0,1.d0/),(/2,2/))
  real ( double ), dimension(3,3), parameter :: eye3_D = reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))

  ! Set default precision for reals
  integer, parameter :: prec = double

end module Types


