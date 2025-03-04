!C    
!C***
!C*** RCM
!C***
!C
      subroutine RCM  (N, NP, NLU, INLU, IALU, 
     &                 NEWtoOLD, OLDtoNEW)

      implicit REAL*8(A-H,O-Z)
      integer, dimension(NP)    :: INLU, NEWtoOLD, OLDtoNEW
      integer, dimension(NLU,NP):: IALU

      integer, dimension(:)  , allocatable :: INLUw, IW
      integer, dimension(:,:), allocatable :: IALUw

!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      allocate (IW(NP))
 
      IW = 0

       INmin= NP
      NODmin= 0

      do i= 1, N
        if (INLU(i).lt.INmin) then
          INmin = INLU(i)
          NODmin= i
        endif
      enddo
 200  continue


      if (NODmin.eq.0) NODmin= 1

      IW(NODmin)= 1

      NEWtoOLD(1     )= NODmin
      OLDtoNEW(NODmin)= 1

      icol= 1
!C===

!C
!C +----------------+
!C | CM-redordering |
!C +----------------+
!C===
      icouG= 1
      do icol= 2, N
        do i= 1, N
          if (IW(i).eq.icol-1) then
            do k= 1, INLU(i)
              in= IALU(k,i)
              if (IW(in).eq.0) then
                IW(in)= icol
                icouG= icouG + 1
              endif
            enddo
          endif
        enddo
        if (icouG.eq.N) exit
      enddo

      NCOLORtot= icol
      icoug    = 0
      do ic= 1, NCOLORtot
        do i= 1, N
          if (IW(i).eq.ic) then
            icoug= icoug + 1
            NEWtoOLD(icoug)= i
            OLDtoNEW(i    )= icoug
          endif
        enddo
      enddo
!C===

!C
!C +-----------------+
!C | MATRIX transfer |
!C +-----------------+
!C===
      allocate (INLUw(NP), IALUw(NLU,NP))
      do j= 1, NLU
        do i= 1, NP
          IW(i) = IALU(j,NEWtoOLD(i))
        enddo
        do i= 1, NP
          IALU(j,i) = IW(i)
        enddo
      enddo

      do i= 1, NP
        IW(i) = INLU(NEWtoOLD(i))
      enddo

      do i= 1, NP
        INLUw(i) = IW(i)
      enddo

      do j= 1, NLU
        do i= 1, NP
          if (IALU(j,i).eq.0) then
            IALUw(j,i) = 0
           else
            IALUw(j,i) = OLDtoNEW(IALU(j,i))
          endif
        enddo
      enddo

      INLU= 0
      IALU= 0

      do i= 1, N
        INLU(i)= INLUw(i)
        do j= 1, INLUw(i)
          IALU(j,i)= IALUw(j,i)
        enddo
      enddo
!C===
      deallocate (IW, INLUw, IALUw)

      return
      end



