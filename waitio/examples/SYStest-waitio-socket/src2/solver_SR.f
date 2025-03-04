!C
!C*** 
!C*** module solver_SR
!C***
!C
      module solver_SR
      contains
!C
!C*** SOLVER_SEND_RECV
!C
      subroutine SOLVER_SEND_RECV                                       &
     &                ( N0, N, NEIBPETOT, NEIBPE,                       &
     &                                        STACK_IMPORT,             &
     &                                        STACK_EXPORT, NOD_EXPORT, &
     &                  WS, X, my_rank)

      implicit REAL*8 (A-H,O-Z)
      include  'mpif.h'
      include  'precision.inc'

      integer(kind=kint )                , intent(in)   ::  N, N0
      integer(kind=kint )                , intent(in)   ::  NEIBPETOT
      integer(kind=kint ), pointer :: NEIBPE      (:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:), NOD_EXPORT(:)
      real   (kind=kreal), dimension(N), intent(inout):: WS
      real   (kind=kreal), dimension(N), intent(inout):: X
      integer                            , intent(in)   :: my_rank

      integer(kind=kint ), dimension(:,:), save,allocatable :: sta1
      integer(kind=kint ), dimension(:  ), save,allocatable :: req1

      integer, save :: NFLAG
      data NFLAG/0/

!C
!C-- INIT.
      if (NFLAG.eq.0) then
        allocate (sta1(MPI_STATUS_SIZE,2*NEIBPETOT))
        allocate (req1(2*NEIBPETOT))
        NFLAG= 1
      endif
!C
!C-- SEND
      do neib= 1, NEIBPETOT
        istart= STACK_EXPORT(neib-1)
        inum  = STACK_EXPORT(neib  ) - istart
!$omp parallel do private (k,ii)
        do k= istart+1, istart+inum
              ii= NOD_EXPORT(k)
           WS(k)= X(ii)
        enddo
        call MPI_Isend (WS(istart+1), inum,  MPI_DOUBLE_PRECISION,      &
     &                  NEIBPE(neib), 0, MPI_COMM_WORLD, req1(neib),    &
     &                  ierr)
      enddo

!C
!C-- RECV
      do neib= 1, NEIBPETOT
        istart= STACK_IMPORT(neib-1)
        inum  = STACK_IMPORT(neib  ) - istart
        call MPI_Irecv (X(istart+N0+1),inum,MPI_DOUBLE_PRECISION,       &
     &                  NEIBPE(neib), 0, MPI_COMM_WORLD,                &
     &                  req1(neib+NEIBPETOT), ierr)
      enddo
      call MPI_Waitall (2*NEIBPETOT, req1, sta1, ierr)

      end subroutine SOLVER_SEND_RECV
      end module     solver_SR



