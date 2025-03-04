!C
!C***
!C*** HPCMW_INIT
!C***
!C
!C    INIT. HPCMW-FEM process's
!C
      subroutine HPCMW_INIT
      use hpcmw_fem_cntl
      implicit REAL*8 (A-H,O-Z)

      call MPI_INIT      (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr )
      call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr )

      hpcmwRarray= 0.d0
      hpcmwIarray= 0

      return
      end
