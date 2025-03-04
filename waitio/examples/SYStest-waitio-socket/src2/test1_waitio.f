      program SOLVER33_TEST_SMP

      use solver11
      use hpcmw_all
      implicit REAL*8(A-H,O-Z)
      include 'waitio_mpif.h'
      integer, dimension(:), allocatable :: OLDtoNEWpe
      integer ::  dest
      integer(kind=kint ), dimension(:,:), save,allocatable :: sta1
      integer(kind=kint ), dimension(:,:), save,allocatable :: req1

      allocate (sta1(WAITIO_STATUS_SIZE,1*6))
      allocate (req1(WAITIO_REQUEST_SIZE,1*6))

!C
!C +-------+
!C | INIT. |
!C +-------+
!C=== 
      call HPCMW_INIT
      call INPUT_CNTL
      allocate (OLDtoNEWpe(PETOT))

      call INPUT_GRID (OLDtoNEWpe, ITERkk)
      ITERkk= 0
      call MAT_CON0 (ITERkk)
      call MAT_CON1 

      allocate (X(NP), X0(NP))

      do j= 1, NP
        X (j)= 0.d0
      enddo

      ITERmaxTIME= hpcmwIarray(10)
      DelT       = hpcmwRarray(10)
!C===
!c      open (12,file='s1.out',status='unknown')
      open (22,file='s2.out',status='unknown')

      TIME= 0.d0
      call MPI_Barrier (MPI_COMM_WORLD, ierr)

      S0_time= MPI_WTIME()
!C
!C>>>**************************************************************** Unstady Loop - START
      do ITERkk= 1, ITERmaxTIME

        TIME= TIME + DelT
        if (my_rank.eq.0) then 
          write (*,'(/a,i8,1pe16.6)') '*** STEP', ITERkk, TIME

          do 
!c            read (12,*,end=100) ttt, Source
             dest=0
             call WAITIO_MPI_Irecv (ttt, 1,  
     &            WAITIO_MPI_DOUBLE_PRECISION,      
     &            dest, 0, WAITIO_SOLVER_COMM, req1(1,1), ierr)
             call WAITIO_MPI_Irecv (Source, 1,  
     &            WAITIO_MPI_DOUBLE_PRECISION,      
     &            dest, 0, WAITIO_SOLVER_COMM, req1(1,2), ierr)

             call WAITIO_MPI_Waitall (2, req1, sta1, ierr)
             if (ttt.ge.TIME) goto 100
          enddo
 100      continue
!c          write (*, '(/a,2(1pe16.6))') '<< get from wrank=0 ', ttt, 
!c     &            Source
        endif

        do j= 1, NP
          X0(j)= X(j)
        enddo

        S1_time= MPI_WTIME()

        call MAT_ASS_MAIN_1
        call MAT_ASS_BC 

        if (my_rank.eq.0) then
          ii= OLDtoNEW(1)
          B(ii)= B(ii) + Source
        endif

        call SOLVE11 (hpcmwIarray, hpcmwRarray, ITERkk)
        E1_time= MPI_WTIME()

        if (my_rank.eq.PETOT-1) then
          ii3= OLDtoNEW(N)          
          VAL= X(ii3)
        endif

        call MPI_Bcast (VAL, 1, MPI_DOUBLE_PRECISION,                   &
     &                  PETOT-1, MPI_COMM_WORLD, ierr)

        if (my_rank.eq.0) then
          ii1= OLDtoNEW(1)
          ii2= OLDtoNEW(2)
          write (*, '(6(1pe16.6))') TIME, ttt, Source, X(ii1), 
     &                             X(ii2), VAL
          write (22,'(6(1pe16.6))') TIME, ttt, Source, X(ii1), 
     &                             X(ii2), VAL
          dest= 2
          call WAITIO_MPI_Isend (TIME, 1,  
     &         WAITIO_MPI_DOUBLE_PRECISION,      
     &         dest, 0, WAITIO_SOLVER_COMM, req1(1,1), ierr)
          call WAITIO_MPI_Isend (ttt, 1,  
     &         WAITIO_MPI_DOUBLE_PRECISION,      
     &         dest, 0, WAITIO_SOLVER_COMM, req1(1,2), ierr)
          call WAITIO_MPI_Isend (Source, 1,  
     &         WAITIO_MPI_DOUBLE_PRECISION,      
     &         dest, 0, WAITIO_SOLVER_COMM, req1(1,3), ierr)
          call WAITIO_MPI_Isend (X(ii1), 1,  
     &         WAITIO_MPI_DOUBLE_PRECISION,      
     &         dest, 0, WAITIO_SOLVER_COMM, req1(1,4), ierr)
          call WAITIO_MPI_Isend (X(ii2), 1,  
     &         WAITIO_MPI_DOUBLE_PRECISION,      
     &         dest, 0, WAITIO_SOLVER_COMM, req1(1,5), ierr)
          call WAITIO_MPI_Isend (VAL, 1,  
     &         WAITIO_MPI_DOUBLE_PRECISION,      
     &         dest, 0, WAITIO_SOLVER_COMM, req1(1,6), ierr)

          call WAITIO_MPI_Waitall (6, req1, sta1, ierr)
        endif

      enddo
!C>>>**************************************************************** Unstady Loop - END

      E0_time= MPI_WTIME()
      COMPtime= E0_time-S0_time
      call MPI_Allreduce (COMPtime, CTsum, 1, MPI_DOUBLE_PRECISION,     &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce (COMPtime, CTmax, 1, MPI_DOUBLE_PRECISION,     &
     &                    MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce (COMPtime, CTmin, 1, MPI_DOUBLE_PRECISION,     &
     &                    MPI_MIN, MPI_COMM_WORLD, ierr)

      call MPI_COMM_SIZE (MPI_COMM_WORLD, nPETOT, ierr )

      if (my_rank.eq.0) then
      write (*,'(//a,3(1pe16.6))') '*elapsed', CTmax, CTmin,
     &                                         CTsum/dfloat(nPETOT)
      endif

      close (22)

      call HPCMW_FINALIZE
      end program SOLVER33_TEST_SMP

!C
!C***
!C*** ERROR_EXIT
!C***
!C
      subroutine ERROR_EXIT (IFLAG, my_rank)
      return
      end
      

