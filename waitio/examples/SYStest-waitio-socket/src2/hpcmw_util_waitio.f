!C
!C***
!C*** HPCmw_util
!C***
!C
!C    MPI settings
!C    FILE names
!C
      module hpcmw_util
      include 'mpif.h'
      include 'precision.inc'

      integer :: PEsmpTOT, FTflag
      integer :: PETOT, my_rank, SOLVER_COMM, WAITIO_SOLVER_COMM, errno
      integer :: my_wrank

      character(len= 4 ) ::  penum, penum_left
      character(len= 9 ) ::  fname, rname

      end module hpcmw_util
