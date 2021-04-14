        subroutine REGAVERAGE(
     & coarse
     & ,icoarselo0,icoarselo1
     & ,icoarsehi0,icoarsehi1
     & ,fine
     & ,ifinelo0,ifinelo1
     & ,ifinehi0,ifinehi1
     & ,icoarboxlo0,icoarboxlo1
     & ,icoarboxhi0,icoarboxhi1
     & ,irefboxlo0,irefboxlo1
     & ,irefboxhi0,irefboxhi1
     & ,numfinepercoar
     & ,reftocoar
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1)
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1)
      integer icoarboxlo0,icoarboxlo1
      integer icoarboxhi0,icoarboxhi1
      integer irefboxlo0,irefboxlo1
      integer irefboxhi0,irefboxhi1
      integer numfinepercoar
      integer reftocoar
        integer iic,jjc
        integer iie,jje
        integer iif,jjf
        REAL*8 weight
        weight = (1.0d0) / numfinepercoar
      do jjc = icoarboxlo1,icoarboxhi1
      do iic = icoarboxlo0,icoarboxhi0
      do jje = irefboxlo1,irefboxhi1
      do iie = irefboxlo0,irefboxhi0
        iif = reftocoar*iic + iie
        jjf = reftocoar*jjc + jje
        coarse(iic,jjc) = coarse(iic,jjc)
     & + weight*fine(iif,jjf)
      enddo
      enddo
      enddo
      enddo
        ch_flops=ch_flops+(icoarboxhi0- icoarboxlo0+1)*(icoarboxhi1- ico
     &arboxlo1+1)*(irefboxhi0- irefboxlo0+1)*(irefboxhi1- irefboxlo1+1)*
     &(2+2)
        return
        end
