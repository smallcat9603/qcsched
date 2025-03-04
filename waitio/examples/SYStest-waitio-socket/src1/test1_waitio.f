      implicit REAL*8(A-H,O-Z)
      real(kind=8) :: Interval, ttt, Source
      integer(kind=4) :: dest
      integer(kind=4 ), dimension(:,:), save,allocatable :: sta1
      integer(kind=4 ), dimension(:,:), save,allocatable :: req1

      include 'mpif.h'
      include 'waitio_mpif.h'

      call MPI_Init(ierr)
      open (11,file='s1.dat',status='unknown')
      read (11,*) ITERmax, Period, Interval, Val
      write (*,'(a,i8)')      '##ITERmax ', ITERmax
      write (*,'(a,1pe16.6)') '##Period  ',  Period
      write (*,'(a,1pe16.6)') '##Interval',  Interval
      write (*,'(a,1pe16.6//)') '##Val     ',  Val
      close (11)

      allocate (sta1(WAITIO_STATUS_SIZE,2*2))
      allocate (req1(WAITIO_REQUEST_SIZE,2*2))
      call WAITIO_CREATE_UNIVERSE_PBHEAD (WAITIO_SOLVER_COMM, ierr)

      open (12,file='s1.out',status='unknown')

      pi  = 4.d0*datan(1.d0)
      coef= pi/180.d0 
      time= 0.d0
      S0time= mpi_wtime ()

      do 
         S1time= mpi_wtime ()

         do 
            do iter= 1, ITERmax
               a= 1.d0
            enddo
            E0time= mpi_wtime (ierr)
            if ((E0time-S1time).gt.Interval) exit
         enddo

         ttt= E0time - S0time
         Source= Val * dsin(ttt*coef)
!     rewind (12)
         write (12,'(2(1pe16.6))') ttt, Source
         write (*,'(2(1pe16.6))') ttt, Source
         do dest=1, 2
            call WAITIO_MPI_Isend (ttt, 1,  
     &           WAITIO_MPI_DOUBLE_PRECISION,      
     &          dest, 0, WAITIO_SOLVER_COMM, req1(1,2*(dest-1)+1), ierr)
            if(ierr.ne.0) goto 100
            call WAITIO_MPI_Isend (Source, 1,  
     &           WAITIO_MPI_DOUBLE_PRECISION,      
     &          dest, 0, WAITIO_SOLVER_COMM, req1(1,2*(dest-1)+2), ierr)
            if(ierr.ne.0) goto 100
         enddo
         call WAITIO_MPI_Waitall (4, req1, sta1, ierr)
         if(ierr.ne.0) goto 100

         if(ttt.gt.Period) exit
      enddo
 100  continue
      call MPI_Finalize()
      stop
      end
      

