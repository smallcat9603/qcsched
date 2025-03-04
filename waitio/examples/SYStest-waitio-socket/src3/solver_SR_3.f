!C
!C*** 
!C*** module solver_SR_3
!C***
!C
      module solver_SR_3
      contains
!C
!C*** SOLVER_SEND_RECV
!C
      subroutine  SOLVER_SEND_RECV_3                                    &
     &                ( N0, N, NEIBPETOT, NEIBPE,                       &
     &                                        STACK_IMPORT, NOD_IMPORT, &
     &                                        STACK_EXPORT, NOD_EXPORT, &
     &                  WS, WR, X, SOLVER_COMM,my_rank)

      implicit REAL*8 (A-H,O-Z)
      include  'mpif.h'
      include  'precision.inc'

      integer(kind=kint )                , intent(in)   ::  N, N0
      integer(kind=kint )                , intent(in)   ::  NEIBPETOT
      integer(kind=kint ), pointer :: NEIBPE      (:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:), NOD_IMPORT(:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:), NOD_EXPORT(:)
      real   (kind=kreal), dimension(3*N), intent(inout):: WS
      real   (kind=kreal), dimension(3*N), intent(inout):: WR
      real   (kind=kreal), dimension(3*N), intent(inout):: X
      integer                            , intent(in)   ::SOLVER_COMM
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
!$omp parallel do private (ii)
        do k= istart+1, istart+inum
               ii   = 3*NOD_EXPORT(k)
           WS(3*k-2)= X(ii-2)
           WS(3*k-1)= X(ii-1)
           WS(3*k  )= X(ii  )
        enddo
!$omp end parallel do
        call MPI_ISEND (WS(3*istart+1), 3*inum,MPI_DOUBLE_PRECISION,    &
     &                  NEIBPE(neib), 0, MPI_COMM_WORLD, req1(neib),    &
     &                  ierr)
      enddo

!C
!C-- RECEIVE
      do neib= 1, NEIBPETOT
        istart= STACK_IMPORT(neib-1)
        inum  = STACK_IMPORT(neib  ) - istart
        call MPI_IRECV (X(3*(istart+N0)+1),3*inum,MPI_DOUBLE_PRECISION, &
     &                  NEIBPE(neib), 0, MPI_COMM_WORLD,                &
     &                  req1(neib+NEIBPETOT), ierr)
      enddo
      call MPI_WAITALL (2*NEIBPETOT, req1, sta1, ierr)

      end subroutine solver_send_recv_3

!C
!C*** SOLVER_SEND_RECV_PC
!C
      subroutine  SOLVER_SEND_RECV_3_PC                                 &
     &                ( N0, N, NEIBPETOT, NEIBPE,                       &
     &                                        STACK_IMPORT, NOD_IMPORT, &
     &                                        STACK_EXPORT, NOD_EXPORT, &
     &                  WS, WR, X, SOLVER_COMM, my_rank, REQS)

      implicit REAL*8 (A-H,O-Z)
      include  'mpif.h'
      include  'precision.inc'

      integer(kind=kint )                , intent(in)   ::  N, N0
      integer(kind=kint )                , intent(in)   ::  NEIBPETOT
      integer(kind=kint ), pointer :: NEIBPE      (:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:), NOD_IMPORT(:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:), NOD_EXPORT(:)
      real   (kind=kreal), dimension(3*N), intent(inout):: WS, WR
      real   (kind=kreal), dimension(3*N), intent(inout):: X
      integer                            , intent(in)   ::SOLVER_COMM
      integer                            , intent(in)   :: my_rank

      integer(kind=kint ), dimension(:) , intent(inout) :: REQS

      integer, save :: NFLAG
      data NFLAG/0/

!C
!C-- SEND/RECV
      do neib= 1, NEIBPETOT
        istart= STACK_EXPORT(neib-1)
        inum  = STACK_EXPORT(neib  ) - istart
!$omp parallel do private (ii)
        do k= istart+1, istart+inum
               ii   = 3*NOD_EXPORT(k)
           WS(3*k-2)= X(ii-2)
           WS(3*k-1)= X(ii-1)
           WS(3*k  )= X(ii  )
        enddo
      enddo

      call MPI_Startall(2*NEIBPETOT, REQS(1), ierr)
      call MPI_Waitall (2*NEIBPETOT, REQS(1), MPI_STATUSES_IGNORE, ierr) 

      end subroutine solver_send_recv_3_PC
      end module     solver_SR_3



