!C
!C***
!C*** MAT_CON1
!C***
!C
      subroutine MAT_CON1
      use hpcmw_all
      implicit REAL*8 (A-H,O-Z)
      integer, dimension(0:1024) :: IW3
      integer, dimension(:,:), allocatable :: IW2

      call INIT_ORDERING
      call PREP_SMP
      call FINAL_ORDERING
      call ONEdim_ORDERING

      if (PETOT.ne.1) then
        do k= 1, NOD_STACK_EXPORT(NEIBPETOT)
          i= NOD_EXPORT(k)
          in= OLDtoNEW(i)
          NOD_EXPORT_NEW(k)= in
        enddo

        do k= 1, NOD_STACK_EXPORT(NEIBPETOT)
          NOD_EXPORT(k)= NOD_EXPORT_NEW(k)
        enddo
      endif

      contains
!C
!C***
!C*** INIT_ORDERING
!C***
!C
      subroutine INIT_ORDERING

      if (allocated(IW)) deallocate (IW)
      allocate (IW(NP), IVECT(0:NP))
      allocate (OLDtoNEW(NP), NEWtoOLD(NP))

      IVECT   = 0
      IW      = 0

      do i= 1, NP
        OLDtoNEW(i)= i
        NEWtoOLD(i)= i
      enddo

      NLmax= -NP
      NUmax= -NP

!      call RCM (N,NP,NL,INL,IAL,NEWtoOLD,OLDtoNEW)

      do i= 1, NP
        NLmax= max(NLmax, INL(i))
      enddo

      NCOLORtot= NHYP
!      if (my_rank.eq.0) write (*,'(a,i8)') 'color number: ', NCOLORtot

      end subroutine INIT_ORDERING
!C
!C***
!C*** PREP_SMP
!C***
!C
      subroutine PREP_SMP
      end subroutine PREP_SMP

!C
!C***
!C*** FINAL_ORDERING
!C***
!C
      subroutine FINAL_ORDERING

      N2= 256
      NU= 26
      NL= 26

      INL= 0
      IAL= 0

*voption novec
      do icel= 1, ICELTOT
        in1= OLDtoNEW(ICELNOD(icel,1))
        in2= OLDtoNEW(ICELNOD(icel,2))
        in3= OLDtoNEW(ICELNOD(icel,3))
        in4= OLDtoNEW(ICELNOD(icel,4))
        in5= OLDtoNEW(ICELNOD(icel,5))
        in6= OLDtoNEW(ICELNOD(icel,6))
        in7= OLDtoNEW(ICELNOD(icel,7))
        in8= OLDtoNEW(ICELNOD(icel,8))

        call FIND_TS_NODE (in1,in2)
        call FIND_TS_NODE (in1,in3)
        call FIND_TS_NODE (in1,in4)
        call FIND_TS_NODE (in1,in5)
        call FIND_TS_NODE (in1,in6)
        call FIND_TS_NODE (in1,in7)
        call FIND_TS_NODE (in1,in8)

        call FIND_TS_NODE (in2,in1)
        call FIND_TS_NODE (in2,in3)
        call FIND_TS_NODE (in2,in4)
        call FIND_TS_NODE (in2,in5)
        call FIND_TS_NODE (in2,in6)
        call FIND_TS_NODE (in2,in7)
        call FIND_TS_NODE (in2,in8)

        call FIND_TS_NODE (in3,in1)
        call FIND_TS_NODE (in3,in2)
        call FIND_TS_NODE (in3,in4)
        call FIND_TS_NODE (in3,in5)
        call FIND_TS_NODE (in3,in6)
        call FIND_TS_NODE (in3,in7)
        call FIND_TS_NODE (in3,in8)

        call FIND_TS_NODE (in4,in1)
        call FIND_TS_NODE (in4,in2)
        call FIND_TS_NODE (in4,in3)
        call FIND_TS_NODE (in4,in5)
        call FIND_TS_NODE (in4,in6)
        call FIND_TS_NODE (in4,in7)
        call FIND_TS_NODE (in4,in8)

        call FIND_TS_NODE (in5,in1)
        call FIND_TS_NODE (in5,in2)
        call FIND_TS_NODE (in5,in3)
        call FIND_TS_NODE (in5,in4)
        call FIND_TS_NODE (in5,in6)
        call FIND_TS_NODE (in5,in7)
        call FIND_TS_NODE (in5,in8)

        call FIND_TS_NODE (in6,in1)
        call FIND_TS_NODE (in6,in2)
        call FIND_TS_NODE (in6,in3)
        call FIND_TS_NODE (in6,in4)
        call FIND_TS_NODE (in6,in5)
        call FIND_TS_NODE (in6,in7)
        call FIND_TS_NODE (in6,in8)

        call FIND_TS_NODE (in7,in1)
        call FIND_TS_NODE (in7,in2)
        call FIND_TS_NODE (in7,in3)
        call FIND_TS_NODE (in7,in4)
        call FIND_TS_NODE (in7,in5)
        call FIND_TS_NODE (in7,in6)
        call FIND_TS_NODE (in7,in8)

        call FIND_TS_NODE (in8,in1)
        call FIND_TS_NODE (in8,in2)
        call FIND_TS_NODE (in8,in3)
        call FIND_TS_NODE (in8,in4)
        call FIND_TS_NODE (in8,in5)
        call FIND_TS_NODE (in8,in6)
        call FIND_TS_NODE (in8,in7)
      enddo

      do in= 1, NP
        NN= INL(in)
        do k= 1, NN
          if (k.gt.26) write (*,*) "###ERR011", my_rank
          NCOL1(k)= IAL(k,in)
        enddo
          call mSORT (NCOL1, NCOL2, NN)
        do k= NN, 1, -1
          if (NN-k+1.gt.26) write (*,*) "###ERR012", my_rank
          IAL(NN-k+1,in)= NCOL1(NCOL2(k))
        enddo        
      enddo
!C===

      end subroutine FINAL_ORDERING

!C
!C***
!C*** ONEdim_ORDERING
!C***
!C
      subroutine ONEdim_ORDERING
      
      if (allocated(indexL)) deallocate (indexL)
      allocate (indexL(0:NP))      

      indexL(0)= 0

!$omp parallel do private(j)
      do j= 1, N
        indexL(j)= 0
      enddo

      do j= N+1, NP
        indexL(j)= 0
      enddo

      do i= 1, NP
        indexL(   i)= indexL(   i-1) + INL(i)
      enddo

      NPL= indexL(NP)
      if (allocated(itemL)) deallocate (itemL)
      allocate (itemL(NPL))

!$omp parallel do private(j,jsL,jeL,k)
      do j = 1, N
        jsL= indexL(j-1)+1
        jeL= indexL(j)
        do k= jsL, jeL
          itemL(k)= 0
        enddo
      enddo

      do i= 1, NP
        do k= 1, INL(i)
                kk = k + indexL(i-1)
          if (k.gt.26) write (*,*) "###ERR021", my_rank
          itemL(kk)=        IAL(k,i)
        enddo
      enddo

      deallocate (INL, IAL)
      end subroutine ONEdim_ORDERING

        subroutine FIND_TS_NODE (ip1,ip2)

          do kk= 1, INL(ip1)
            if (ip2.eq.IAL(kk,ip1)) return
          enddo
          icou= INL(ip1) + 1
          if (icou.gt.26) write (*,*) "###ERR022", my_rank
          IAL(icou,ip1)= ip2
          INL(     ip1)= icou
          return

        end subroutine FIND_TS_NODE

      end subroutine MAT_CON1
