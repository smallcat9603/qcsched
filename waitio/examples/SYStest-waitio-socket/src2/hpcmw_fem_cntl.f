!C
!C***
!C*** HPCmw_fem_cntl
!C***
!C
      module hpcmw_fem_cntl
      use hpcmw_util

      integer(kind=kint ), dimension(100) :: hpcmwIarray
      real   (kind=kreal), dimension(100) :: hpcmwRarray

      real   (kind=kreal), parameter :: O8th= 0.125d0

      end module hpcmw_fem_cntl
