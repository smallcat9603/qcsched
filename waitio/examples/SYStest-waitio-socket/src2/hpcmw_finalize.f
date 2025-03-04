!C
!C***
!C*** HPCMW_FINALIZE
!C***
!C
      subroutine HPCMW_FINALIZE
      use hpcmw_fem_cntl
      implicit REAL*8 (A-H,O-Z)

      call MPI_FINALIZE (ierr)
      if (my_rank.eq.0) stop ' * normal termination'

      return
      end
