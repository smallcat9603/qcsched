      implicit REAL*8(A-H,O-Z)
      real(kind=8) :: Interval
      include 'mpif.h'

      call MPI_Init(ierr)
      open (11,file='s1.dat',status='unknown')
      read (11,*) ITERmax, Period, Interval, Val
      write (*,'(a,i8)')      '##ITERmax ', ITERmax
      write (*,'(a,1pe16.6)') '##Period  ',  Period
      write (*,'(a,1pe16.6)') '##Interval',  Interval
      write (*,'(a,1pe16.6//)') '##Val     ',  Val
      close (11)

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
!      rewind (12)
        write (12,'(2(1pe16.6))') ttt, Source
      enddo

      call MPI_Finalize()
      stop
      end
      

