!C
!C*** 
!C*** module solver_CG
!C***
!C
      module solver_CG
      contains
!C
!C*** CG
!C
!C    CG solves the linear system Ax = b 
!C    using the Conjugate Gradient iterative method with the following
!C    preconditioners for SMP nodes:
!C
      subroutine CG                                                     &
     &                 (N, NP, NPL,                                     &
     &                  D, AL, INL, IAL,                                &
     &                  B,  X, RESID,  ITER, ERROR,                     &
     &                  my_rank, NEIBPETOT, NEIBPE,                     &
     &                  STACK_IMPORT,                                   &
     &                  STACK_EXPORT, NOD_EXPORT)

      use  solver_SR

      implicit REAL*8(A-H,O-Z)
      include  'precision.inc'
      include  'mpif.h'

      integer(kind=kint ), intent(in):: N, NP, NPL, my_rank
      integer(kind=kint ), intent(in):: NEIBPETOT

      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID

      real(kind=kreal), dimension(NP) , intent(inout):: B, X
      real(kind=kreal), dimension(NPL), intent(inout):: AL
      real(kind=kreal), dimension(NP) , intent(inout):: D

      integer(kind=kint ), dimension(0:NP ),intent(in) :: INL
      integer(kind=kint ), dimension(  NPL),intent(in) :: IAL

      real(kind=kreal), dimension(:,:), allocatable :: WW
      real(kind=kreal), dimension(:)  , allocatable :: WS

      integer(kind=kint ), pointer :: NEIBPE(:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:), NOD_EXPORT(:)

      integer(kind=kint), parameter ::  R= 1
      integer(kind=kint), parameter ::  Z= 2
      integer(kind=kint), parameter ::  Q= 2
      integer(kind=kint), parameter ::  P= 3
      integer(kind=kint), parameter :: DD= 4

      integer(kind=kint ) :: MAXIT, IFLAG

      real   (kind=kreal) :: TOL
      data IFLAG/0/

!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      ERROR= 0

      allocate (WS(NP), WW(NP+128,4))

      WW= 0.d0
      COMMtime= 0.d0
      COMPtime= 0.d0

      MAXIT  = ITER
       TOL   = RESID           

!$omp parallel do
      do i= 1, N
        WW(i,1)= 0.d0
        WW(i,2)= 0.d0
        WW(i,3)= 0.d0
        WW(i,4)= 0.d0
        WS(i  )= 0.d0
        X (i  )= 0.d0
      enddo

      do i= 1+N, NP
        WW(i,1)= 0.d0
        WW(i,2)= 0.d0
        WW(i,3)= 0.d0
        WW(i,4)= 0.d0
        WS(i  )= 0.d0
        X (i  )= 0.d0
      enddo
!C===

!C
!C +-----------------------+
!C | {r0}= {b} - [A]{xini} |
!C +-----------------------+
!C===

!$omp parallel do private(j)
      do j= 1, N
        WW(j,R )= B(j)
        WW(j,DD)= 1.d0/D(j)
      enddo


      BNRM20= 0.d0
!$omp parallel do private(i) reduction (+:BNRM20)
      do i= 1, N
        BNRM20= BNRM20+B(i)**2
      enddo

      call MPI_Allreduce (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
      if (BNRM2.eq.0.d0) BNRM2= 1.d0
      ITER = 0
!C===

      S1_time= MPI_WTIME()
!      call start_collection ("solver")
      do iter= 1, MAXIT
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!C
!C************************************************* Conjugate Gradient Iteration

!C
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===

!$omp parallel do private (i)
      do i= 1, N
        WW(i,Z)= WW(i,DD)*WW(i,R)
      enddo
!C===
      
!C
!C +---------------+
!C | {RHO}= {r}{z} |
!C +---------------+
!C===
      RHO0= 0.d0

!$omp parallel do private(i) reduction(+: RHO0)
      do i= 1, N
        RHO0= RHO0 + WW(i,R)*WW(i,Z)
      enddo

      call MPI_Allreduce (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,           &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
!C===

!C
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
      if ( ITER.eq.1 ) then

!$omp parallel do
        do i= 1, N
          WW(i,P)= WW(i,Z)
        enddo
       else
         BETA= RHO / RHO1

!$omp parallel do
         do i= 1, N
           WW(i,P)= WW(i,Z) + BETA*WW(i,P)
         enddo
      endif
!C===

!C
!C +-------------+
!C | {q}= [A]{p} |
!C +-------------+
!C===        

      call SOLVER_SEND_RECV                                             &
     &   ( N, NP, NEIBPETOT, NEIBPE, STACK_IMPORT,                      &
     &     STACK_EXPORT, NOD_EXPORT, WS, WW(1,P), my_rank )

!$omp parallel do private (j,k)
      do j= 1, N
        WW(j,Q)= D(j)*WW(j,P)
      enddo

!$omp parallel do private (j,k)
      do j= 1, N
      do k= INL(j-1)+1, INL(j)
        WW(j,Q)= WW(j,Q) + AL(k)*WW(IAL(k),P)
      enddo
      enddo
!C===

!C
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
      C10= 0.d0
!$omp parallel do private(i) reduction(+:C10)
      do i= 1, N
        C10= C10 + WW(i,P)*WW(i,Q)
      enddo

      call MPI_Allreduce (C10, C1, 1, MPI_DOUBLE_PRECISION,             &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
      ALPHA= RHO / C1
!C===

!C
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===

!$omp parallel do private(i) shared (ALPHA)
      do i= 1, N
         X(i)  = X (i)   + ALPHA * WW(i,P)
        WW(i,R)= WW(i,R) - ALPHA * WW(i,Q)
      enddo

      DNRM20= 0.d0
!$omp parallel do private(i) reduction(+:DNRM20)
      do i= 1, N
        DNRM20= DNRM20 + WW(i,R)**2 
      enddo

      call MPI_Allreduce (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,       &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
      RESID= dsqrt(DNRM2/BNRM2)

!C##### ITERATION HISTORY
!        if (my_rank.eq.0) write (*,1000) ITER, RESID
 1000   format (i5, 1pe16.6)
! 1010   format (1pe16.6)
!C#####

        if ( RESID.le.TOL   ) exit
        if ( ITER .eq.MAXIT ) ERROR= -300

        RHO1 = RHO                                                             
      enddo
!C===

!C
!C-- INTERFACE data EXCHANGE
   30 continue
!      call stop_collection ("solver")
      E1_time= MPI_WTIME()

      if (ERROR.eq.-300) ITER= MAXIT
      if (my_rank.eq.0) write (*, 1000) ITER, RESID

      if (my_rank.eq.0) then
        COMPtime= (E1_time - S1_time)/dfloat(ITER)
        write (*,'(a,2(1pe16.6))') '*solver ', E1_time - S1_time
      endif

      call SOLVER_SEND_RECV                                             &
     &   ( N, NP, NEIBPETOT, NEIBPE, STACK_IMPORT,                      &
     &     STACK_EXPORT, NOD_EXPORT, WS, X, my_rank)

      end subroutine        CG
      end module     solver_CG
