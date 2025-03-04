!C
!C***
!C*** HPCmw_solver_cntl
!C***
!C
      module hpcmw_solver_cntl
      use hpcmw_util

      integer(kind=kint) :: METHOD, PRECOND, NSET, iterPREmax
      integer(kind=kint) :: ITER, ITERactual

      real   (kind=kreal) :: RESID, SIGMA_DIAG, SIGMA

      end module hpcmw_solver_cntl
