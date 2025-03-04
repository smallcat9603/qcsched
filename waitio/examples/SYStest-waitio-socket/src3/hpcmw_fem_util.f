!C
!C***
!C*** HPCmw_fem_util
!C***
!C
      module hpcmw_fem_util
      use hpcmw_util

      real(kind=kreal), dimension(2,2,8) :: PNQ, PNE, PNT
      real(kind=kreal), dimension(2)     :: WEI, POS
      integer(kind=kint), dimension(100) :: NCOL1, NCOL2

      real(kind=kreal), dimension(2,2,2,8) :: SHAPE
      real(kind=kreal), dimension(2,2,2,8) :: PNX, PNY, PNZ
      real(kind=kreal), dimension(2,2,2  ) :: DETJ

      end module hpcmw_fem_util
