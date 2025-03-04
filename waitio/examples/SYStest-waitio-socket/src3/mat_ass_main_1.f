!C
!C***
!C*** MAT_ASS_MAIN_1
!C***
!C
      subroutine MAT_ASS_MAIN_1
      use hpcmw_all
      implicit REAL*8 (A-H,O-Z)
      integer(kind=kint), dimension(8) :: nodLOCAL
      integer, save :: NFLAG
      data NFLAG/0/


      if (NFLAG.ne.2) then
!C
!C-- MATRIX init.

      WEI(1)= +1.0000000000D+00
      WEI(2)= +1.0000000000D+00

      POS(1)= -0.5773502692D+00
      POS(2)= +0.5773502692D+00

      if (allocated(AL)) deallocate (AL)
      if (allocated(D )) deallocate (D )
      if (allocated(B )) deallocate (B )
      if (allocated(ALUG)) deallocate (ALUG)

      allocate (AL(NPL))
      allocate (D(NP), B(NP))

        NFLAG= 1
      endif

!$omp parallel do private(ip,j, jS, jE)
        do j= 1, N
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
          B(j)= 0.d0
          D(j)= 0.d0
          jsL= indexL(j-1)+1
          jeL= indexL(j)
          do k= jsL, jeL
            AL(k)= 0.d0
          enddo
        enddo

      if (NFLAG.ne.2) then
!C
!C-- INIT.
!C     PNQ   - 1st-order derivative of shape function by QSI
!C     PNE   - 1st-order derivative of shape function by ETA
!C     PNT   - 1st-order derivative of shape function by ZET
!C
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

        NFLAG= 2
      endif

      DelT= hpcmwRarray(10)
      COND= hpcmwRarray(11)
      RHOC= hpcmwRarray(12)*hpcmwRarray(13)

      QNOD= 0.125d0*DX*DY*DZ*DelT
      RNOD= 0.125d0*DX*DY*DZ*RHOC

      do icel= 1, ICELTOT
        in1= OLDtoNEW(ICELNOD(icel,1))
        in2= OLDtoNEW(ICELNOD(icel,2))
        in3= OLDtoNEW(ICELNOD(icel,3))
        in4= OLDtoNEW(ICELNOD(icel,4))
        in5= OLDtoNEW(ICELNOD(icel,5))
        in6= OLDtoNEW(ICELNOD(icel,6))
        in7= OLDtoNEW(ICELNOD(icel,7))
        in8= OLDtoNEW(ICELNOD(icel,8))
!C
!C== JACOBIAN & INVERSE JACOBIAN
        nodLOCAL(1)= in1
        nodLOCAL(2)= in2
        nodLOCAL(3)= in3
        nodLOCAL(4)= in4
        nodLOCAL(5)= in5
        nodLOCAL(6)= in6
        nodLOCAL(7)= in7
        nodLOCAL(8)= in8
!C
!C== CONSTRUCT the GLOBAL MATRIX
        do ie= 1, 8
          ip = nodLOCAL(ie)
        do je= 1, 8
          jp = nodLOCAL(je)

          IDlu= 0
          if (ip.eq.jp) IDlu= 0

          kk= 0
            iiS= indexL(ip-1) + 1
            iiE= indexL(ip  )
            do k= iiS, iiE
              if ( itemL(k).eq.jp) then
                kk= k
                IDlu= -1
              endif
            enddo

          PNXi= 0.d0
          PNYi= 0.d0
          PNZi= 0.d0
          PNXj= 0.d0
          PNYj= 0.d0
          PNZj= 0.d0

          
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

            AMAT= AMAT + (PNXi*PNXj+PNYi*PNYj+PNZi*PNZj)*coef*DelT
          enddo
          enddo
          enddo

          if (IDlu.ne.0) then
            AL(kk)= AL(kk) + AMAT*COND
           else
            D(ip)= D(ip) + AMAT*COND + RNOD
            B(ip)= B(ip) + RNOD*X0(ip)
          endif

        enddo
        enddo
      enddo


      return
      end
