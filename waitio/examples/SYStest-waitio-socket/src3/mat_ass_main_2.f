!C
!C***
!C*** MAT_ASS_MAIN_2
!C***
!C
      subroutine MAT_ASS_MAIN_2 
      use hpcmw_all
      implicit REAL*8 (A-H,O-Z)
      integer(kind=kint), dimension(8) :: nodLOCAL
      integer(kind=kint) , dimension(:,:,:), allocatable :: NtoELEM
      real   (kind=kreal), dimension(:,:,:), allocatable :: Emat

      Time_1= 0.d0
      Time_2= 0.d0

      S1_time= MPI_WTIME()
!C
!C +------------------+
!C | ELEMENT Coloring |
!C +------------------+
!C===
      allocate (ELMCOLORindex(0:NP))
      allocate (ELMCOLORitem (ICELTOT))
      if (allocated (IWKX)) deallocate (IWKX)
      allocate (IWKX(0:NP,3))

      IWKX= 0
      icou= 0
      do icol= 1, NP
        do i= 1, NP
          IWKX(i,1)= 0
        enddo
        do icel= 1, ICELTOT
          if (IWKX(icel,2).eq.0) then
            in1= ICELNOD(icel,1)
            in2= ICELNOD(icel,2)
            in3= ICELNOD(icel,3)
            in4= ICELNOD(icel,4)
            in5= ICELNOD(icel,5)
            in6= ICELNOD(icel,6)
            in7= ICELNOD(icel,7)
            in8= ICELNOD(icel,8)

            ip1= IWKX(in1,1)
            ip2= IWKX(in2,1)
            ip3= IWKX(in3,1)
            ip4= IWKX(in4,1)
            ip5= IWKX(in5,1)
            ip6= IWKX(in6,1)
            ip7= IWKX(in7,1)
            ip8= IWKX(in8,1)

            isum= ip1 + ip2 + ip3 + ip4 + ip5 + ip6 + ip7 + ip8
            if (isum.eq.0) then 
              icou= icou + 1
              IWKX(icol,3)= icou
              IWKX(icel,2)= icol
              ELMCOLORitem(icou)= icel

              IWKX(in1,1)= 1
              IWKX(in2,1)= 1
              IWKX(in3,1)= 1
              IWKX(in4,1)= 1
              IWKX(in5,1)= 1
              IWKX(in6,1)= 1
              IWKX(in7,1)= 1
              IWKX(in8,1)= 1
              if (icou.eq.ICELTOT) goto 100            
            endif
          endif
        enddo
      enddo

 100  continue
      ELMCOLORtot= icol
      IWKX(0          ,3)= 0
      IWKX(ELMCOLORtot,3)= ICELTOT

!      do icol= 0, ELMCOLORtot
!        ELMCOLORindex(icol)= IWKX(icol,3)
!      enddo

      do icol= 1, ELMCOLORtot
        nn1= IWKX(icol,3) - IWKX(icol-1,3) 
        nn2= nn1/PEsmpTOT
        nnr= nn1 - nn2*PEsmpTOT
        do ip= 1, PEsmpTOT
          IWKX(PEsmpTOT*(icol-1)+ip,2)= nn2
        enddo
        do ip= 1, nnr
          IWKX(PEsmpTOT*(icol-1)+ip,2)= nn2 + 1
        enddo
      enddo

      ELMCOLORindex(0)= 0
      do icol= 1, ELMCOLORtot
        do ip= 1, PEsmpTOT
          is1= PEsmpTOT*(icol-1) + ip      
          is0= PEsmpTOT*(icol-1) + ip - 1
          ELMCOLORindex(is1)= ELMCOLORindex(is0) + IWKX(is1,2)      
        enddo
      enddo

!      write (*,'(a,2i8)') '### Number of Element Colors', 
!     &                     my_rank, ELMCOLORtot
!      deallocate (IWKX)
!C===

!C
!C-- MATRIX init.

      WEI(1)= +1.0000000000D+00
      WEI(2)= +1.0000000000D+00

      POS(1)= -0.5773502692D+00
      POS(2)= +0.5773502692D+00

      if (allocated(AL)) deallocate (AL)
      if (allocated(D )) deallocate (D )
      if (allocated(B )) deallocate (B )
      if (allocated(X )) deallocate (X )

      allocate (AL(NPL))
      allocate (D(NP), B(NP), X(NP))

!$omp parallel do private(ip,j, jS, jE)
        do j= 1, N
          X(j)= 0.d0
          B(j)= 0.d0
          D(j)= 0.d0
        enddo


!$omp parallel do private(j,jsL,jeL,k)
        do j= 1, N
          jsL= indexL(j-1)+1
          jeL= indexL(j)
          do k= jsL, jeL
            AL(k)= 0.d0
          enddo
        enddo

        do j= 1+N, NP
          X(j)= 0.d0
          B(j)= 0.d0
          D(j)= 0.d0
          jsL= indexL(j-1)+1
          jeL= indexL(j)
          do k= jsL, jeL
            AL(k)= 0.d0
          enddo
        enddo

!C
!C-- INIT.
!C     PNQ   - 1st-order derivative of shape function by QSI
!C     PNE   - 1st-order derivative of shape function by ETA
!C     PNT   - 1st-order derivative of shape function by ZET
!C
      QNOD= 0.125d0

      do kp= 1, 2
      do jp= 1, 2
      do ip= 1, 2
        QP1= 1.d0 + POS(ip)
        QM1= 1.d0 - POS(ip)
        EP1= 1.d0 + POS(jp)
        EM1= 1.d0 - POS(jp)
        TP1= 1.d0 + POS(kp)
        TM1= 1.d0 - POS(kp)
        SHAPE(ip,jp,kp,1)= O8th * QM1 * EM1 * TM1
        SHAPE(ip,jp,kp,2)= O8th * QP1 * EM1 * TM1
        SHAPE(ip,jp,kp,3)= O8th * QP1 * EP1 * TM1
        SHAPE(ip,jp,kp,4)= O8th * QM1 * EP1 * TM1
        SHAPE(ip,jp,kp,5)= O8th * QM1 * EM1 * TP1
        SHAPE(ip,jp,kp,6)= O8th * QP1 * EM1 * TP1
        SHAPE(ip,jp,kp,7)= O8th * QP1 * EP1 * TP1
        SHAPE(ip,jp,kp,8)= O8th * QM1 * EP1 * TP1
        PNQ(jp,kp,1)= - O8th * EM1 * TM1
        PNQ(jp,kp,2)= + O8th * EM1 * TM1
        PNQ(jp,kp,3)= + O8th * EP1 * TM1
        PNQ(jp,kp,4)= - O8th * EP1 * TM1
        PNQ(jp,kp,5)= - O8th * EM1 * TP1
        PNQ(jp,kp,6)= + O8th * EM1 * TP1
        PNQ(jp,kp,7)= + O8th * EP1 * TP1
        PNQ(jp,kp,8)= - O8th * EP1 * TP1
        PNE(ip,kp,1)= - O8th * QM1 * TM1
        PNE(ip,kp,2)= - O8th * QP1 * TM1
        PNE(ip,kp,3)= + O8th * QP1 * TM1
        PNE(ip,kp,4)= + O8th * QM1 * TM1
        PNE(ip,kp,5)= - O8th * QM1 * TP1
        PNE(ip,kp,6)= - O8th * QP1 * TP1
        PNE(ip,kp,7)= + O8th * QP1 * TP1
        PNE(ip,kp,8)= + O8th * QM1 * TP1
        PNT(ip,jp,1)= - O8th * QM1 * EM1
        PNT(ip,jp,2)= - O8th * QP1 * EM1
        PNT(ip,jp,3)= - O8th * QP1 * EP1
        PNT(ip,jp,4)= - O8th * QM1 * EP1
        PNT(ip,jp,5)= + O8th * QM1 * EM1
        PNT(ip,jp,6)= + O8th * QP1 * EM1
        PNT(ip,jp,7)= + O8th * QP1 * EP1
        PNT(ip,jp,8)= + O8th * QM1 * EP1
      enddo
      enddo
      enddo

      S2_time= MPI_WTIME()

!      call start_collection ("mat_ass_all")
      do icol= 1, ELMCOLORtot
!$omp parallel do private (icel0,icel,in1,in2,in3,in4,in5,in6,in7,in8)  &
!$omp&            private (in10,in20,in30,in40,in50,in60,in70,in80)     &
!$omp&            private (nodLOCAL,ie,je,ip,jp,kk,iiS,iiE,iDlu,k)      &
!$omp&            private (PNXi,PNYi,PNZi,PNXj,PNYj,PNZj,AMAT)          &
!$omp&            private (ipn,jpn,kpn,coef)                            &
!$omp&            private (ib0,IBLOCKx,iS0,iE0,iSk,iEk,icou,iip)        &
!$omp&            private (NtoELEM,EMAT)                                &
!$omp&            private (DETJ,PNX,PNY,PNZ,X1,X2,X3,X4,X5,X6,X7,X8)    &
!$omp&            private (Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)                     &
!$omp&            private (Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8)
      do iip= 1, PEsmpTOT
        if (.not.allocated(NtoELEM)) allocate (NtoELEM(8,8,NBLOCK))
        if (.not.allocated(EMAT   )) allocate (EMAT   (8,8,NBLOCK))

        icou= 0
        iS0=  ELMCOLORindex(PEsmpTOT*(icol-1)+iip-1)
        iE0=  ELMCOLORindex(PEsmpTOT*(icol-1)+iip  )
        IBLOCKx= (iE0-iS0-1)/NBLOCK + 1
      do ib0= 1, IBLOCKx
        iSk= iS0 + (ib0-1)*NBLOCK + 1
        iEk= iSk + NBLOCK - 1
        if (iEk.ge.iE0) iEk= iE0 

!C
!C-- Address of Global Matrix
      do icel0= iSk, iEk
        icel= ELMCOLORitem(icel0)

        in1= ICELNOD(icel,1)
        in2= ICELNOD(icel,2)
        in3= ICELNOD(icel,3)
        in4= ICELNOD(icel,4)
        in5= ICELNOD(icel,5)
        in6= ICELNOD(icel,6)
        in7= ICELNOD(icel,7)
        in8= ICELNOD(icel,8)

        nodLOCAL(1)= in1
        nodLOCAL(2)= in2
        nodLOCAL(3)= in3
        nodLOCAL(4)= in4
        nodLOCAL(5)= in5
        nodLOCAL(6)= in6
        nodLOCAL(7)= in7
        nodLOCAL(8)= in8

        do je= 1, 8
          jp = nodLOCAL(je)
        do ie= 1, 8
          ip = nodLOCAL(ie)

          kk = 0

          if (ip.eq.jp) then
            kk= ip
           else
            iiS= indexL(ip-1) + 1
            iiE= indexL(ip  )
            do k= iiS, iiE
              if ( itemL(k).eq.jp ) then
                kk  = k
                exit
              endif
            enddo
          endif

          NtoELEM(ie,je,icel0+1-iSk)= kk          
        enddo
        enddo
      enddo
!C
!C-- MATRIX Assembly: LOCAL
      do icel0= iSk, iEk
        icel= ELMCOLORitem(icel0)

        in1= ICELNOD(icel,1)
        in2= ICELNOD(icel,2)
        in3= ICELNOD(icel,3)
        in4= ICELNOD(icel,4)
        in5= ICELNOD(icel,5)
        in6= ICELNOD(icel,6)
        in7= ICELNOD(icel,7)
        in8= ICELNOD(icel,8)

        nodLOCAL(1)= in1
        nodLOCAL(2)= in2
        nodLOCAL(3)= in3
        nodLOCAL(4)= in4
        nodLOCAL(5)= in5
        nodLOCAL(6)= in6
        nodLOCAL(7)= in7
        nodLOCAL(8)= in8

        X1= 0.d0
        Y1= 0.d0
        Z1= 0.d0

        X2= DX
        Y2= 0.d0
        Z2= 0.d0

        X3= DX
        Y3= DY
        Z3= 0.d0

        X4= 0.d0
        Y4= DY
        Z4= 0.d0

        X5= 0.d0
        Y5= 0.d0
        Z5= DZ

        X6= DX
        Y6= 0.d0
        Z6= DZ

        X7= DX
        Y7= DY
        Z7= DZ

        X8= 0.d0
        Y8= DY
        Z8= DZ

        call JACOBI (DETJ, PNQ, PNE, PNT, PNX, PNY, PNZ,                &
     &               X1, X2, X3, X4, X5, X6, X7, X8,                    &
     &               Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8,                    &
     &               Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8 )

        do je= 1, 8
        do ie= 1, 8
          AMAT= 0.d0
          do kpn= 1, 2
          do jpn= 1, 2
          do ipn= 1, 2
            coef= dabs(DETJ(ipn,jpn,kpn))*WEI(ipn)*WEI(jpn)*WEI(kpn)

            PNXi= PNX(ipn,jpn,kpn,ie)
            PNYi= PNY(ipn,jpn,kpn,ie)
            PNZi= PNZ(ipn,jpn,kpn,ie)

            PNXj= PNX(ipn,jpn,kpn,je)
            PNYj= PNY(ipn,jpn,kpn,je)
            PNZj= PNZ(ipn,jpn,kpn,je)

            AMAT= AMAT + (PNXi*PNXj+PNYi*PNYj+PNZi*PNZj)*coef
          enddo
          enddo
          enddo
          EMAT(ie,je,icel0+1-iSK)= AMAT
        enddo
        enddo
      enddo

!C
!C-- MATRIX Assembly: GLOBAL
      do icel0= iSk, iEk
        do je= 1, 8
        do ie= 1, 8
          kk= NtoELEM(ie,je,icel0+1-iSk)
          if (ie.ne.je) then
            AL(kk)= AL(kk) + EMAT(ie,je,icel0+1-iSK)
           else
            D (kk)= D (kk) + EMAT(ie,je,icel0+1-iSK)
            B (kk)= B (kk) + QNOD
          endif
        enddo
        enddo
      enddo

      enddo
      enddo
      enddo

      S3_time= MPI_WTIME()
!      call stop_collection ("mat_ass_all")

      if (my_rank.eq.0) write (*,'(a,1pe16.6)') '*e-color', 
     &                                           S2_time-S1_time
      if (my_rank.eq.0) write (*,'(a,1pe16.6)') '*mat-ass', 
     &                                           S3_time-S2_time

      return
      end
