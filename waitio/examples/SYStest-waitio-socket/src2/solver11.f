      module SOLVER11

      contains
        subroutine SOLVE11 (hpcmwIarray, hpcmwRarray, ITERkk) 

        use hpcmw_solver_matrix
        use hpcmw_solver_cntl
        use hpcmw_fem_mesh

        use solver_CG

        implicit REAL*8 (A-H,O-Z)

        real(kind=kreal), dimension(3,3) :: ALU
        real(kind=kreal), dimension(3)   :: PW  

        integer :: ERROR, ICFLAG
        character(len=char_length) :: BUF

        data ICFLAG/0/

        integer(kind=kint ), dimension(:) :: hpcmwIarray
        real   (kind=kreal), dimension(:) :: hpcmwRarray


!C
!C +------------+
!C | PARAMETERs |
!C +------------+
!C===
        ITER      = hpcmwIarray(1)
        METHOD    = hpcmwIarray(2)
        PRECOND   = hpcmwIarray(3)
        NSET      = hpcmwIarray(4)
        iterPREmax= hpcmwIarray(5)

        RESID     = hpcmwRarray(1)
        SIGMA_DIAG= hpcmwRarray(2)

        if (iterPREmax.lt.1) iterPREmax= 1
        if (iterPREmax.gt.4) iterPREmax= 4
!C===

!C
!C +------------------+
!C | ITERATIVE solver |
!C +------------------+
!C===
        call CG                                                         &
     &   ( N, NP, NPL,                                                  &
     &     D, AL, indexL, itemL,                                        &
     &     B, X,  RESID, ITER, ERROR,                                   &
     &           my_rank, NEIBPETOT, NEIBPE,                            &
     &           NOD_STACK_IMPORT,                                      &
     &           NOD_STACK_EXPORT, NOD_EXPORT)
      ITERactual= ITER
!C===

      end subroutine SOLVE11
      end module SOLVER11
