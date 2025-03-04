      program SOLVER33_TEST_SMP

      use solver11
      use hpcmw_all
      implicit REAL*8(A-H,O-Z)
      integer, dimension(:), allocatable :: OLDtoNEWpe
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
      open (12,file='s1.out',status='unknown')
      open (13,file='s2.out',status='unknown')
      open (22,file='s3.out',status='unknown')

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
            read (13,*,end=100) ttt, tt1, s0, s1, s2, TBOU
            if (ttt.ge.TIME) goto 100
          enddo
 100      continue
          do 
            read (12,*,end=200) ttx, Source
            if (ttx.ge.TIME) goto 200
          enddo
 200      continue
        endif

        call MPI_Bcast (ttt,  1, MPI_DOUBLE_PRECISION,                  &
     &                  0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast (TBOU, 1, MPI_DOUBLE_PRECISION,                  &
     &                  0, MPI_COMM_WORLD, ierr)

        do j= 1, NP
          X0(j)= X(j)
        enddo

        S1_time= MPI_WTIME()

        call MAT_ASS_MAIN_1
        call MAT_ASS_BC (TBOU)

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
          write (*, '(7(1pe16.6))') TIME, ttt, Source, X(ii1), 
     &                             X(ii2), VAL, TBOU
          write (22,'(7(1pe16.6))') TIME, ttt, Source, X(ii1), 
     &                             X(ii2), VAL, TBOU
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
      

