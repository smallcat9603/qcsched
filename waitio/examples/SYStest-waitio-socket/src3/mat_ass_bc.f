!C
!C***
!C*** MAT_ASS_BC
!C***
!C
      subroutine MAT_ASS_BC (TBOU)
      use hpcmw_all
      implicit REAL*8 (A-H,O-Z)

      if (allocated(IWKX)) deallocate (IWKX)
      allocate (IWKX(NP,2))
      IWKX= 0

!C
!C== Z=Zmax

!$omp parallel do private(in)
      do in= 1, NP
        IWKX(in,1)= 0
      enddo

      ib0= 4
!$omp parallel do private(ib,in)
      do ib= BC_STACK(ib0-1)+1, BC_STACK(ib0)
        in= OLDtoNEW(BC_NOD(ib))
        IWKX(in,1)= 1
      enddo

!$omp parallel do private(in,iS,iE,k)
      do in= 1, N
        iS= indexL(in-1) + 1
        iE= indexL(in  )
        do k= iS, iE
          if (IWKX(itemL(k),1).eq.1.and.IWKX(in,1).eq.0) then
            B(itemL(k))= B(itemL(k)) - AL(k)*TBOU
            AL(k)= 0.d0
          endif
        enddo
      enddo

!$omp parallel do private(in,iS,iE,k)
      do in= 1, N
        if (IWKX(in,1).eq.1) then
          D(in)= 1.d0
          B(in)= TBOU
          
          iS= indexL(in-1) + 1
          iE= indexL(in  )
          do k= iS, iE
            AL(k)= 0.d0
          enddo
        endif
      enddo

      deallocate (IWKX)

      return
      end
