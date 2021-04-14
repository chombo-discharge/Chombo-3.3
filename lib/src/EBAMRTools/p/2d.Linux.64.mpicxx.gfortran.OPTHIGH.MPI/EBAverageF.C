#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine EBAVERAGE(
     &           coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,icoarboxlo0,icoarboxlo1
     &           ,icoarboxhi0,icoarboxhi1
     &           ,refrat
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1)
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer icoarboxlo0,icoarboxlo1
      integer icoarboxhi0,icoarboxhi1
      integer refrat
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer icc, jcc
      integer iff, jff
      integer ii, jj
      real_t refscale, coarsesum, frefine
      refscale = one
      frefine = refrat
      do ii = 1, CH_SPACEDIM
         refscale = refscale/frefine
      enddo
      
      do jcc = icoarboxlo1,icoarboxhi1
      do icc = icoarboxlo0,icoarboxhi0

      coarsesum = zero
      
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0

      
      iff = icc*refrat + ii
      jff = jcc*refrat + jj
      coarsesum = coarsesum + fine(iff,jff)
      
      enddo
      enddo
      coarse(icc,jcc) = coarsesum * refscale
      
      enddo
      enddo
      return
      end
      subroutine EBAVERARZ(
     &           coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,icoarboxlo0,icoarboxlo1
     &           ,icoarboxhi0,icoarboxhi1
     &           ,refrat
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           ,dxcoar
     &           ,dxfine
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1)
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer icoarboxlo0,icoarboxlo1
      integer icoarboxhi0,icoarboxhi1
      integer refrat
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      REAL_T dxcoar
      REAL_T dxfine
      integer icc, jcc
      integer iff, jff
      integer ii, jj
      real_t coarsesum, coarsevol, finevol, sumfinevol
      
      do jcc = icoarboxlo1,icoarboxhi1
      do icc = icoarboxlo0,icoarboxhi0

      coarsesum = zero
      sumfinevol = zero
      
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0

      
      iff = icc*refrat + ii
      jff = jcc*refrat + jj
      finevol = (iff + half)*dxfine*dxfine*dxfine
      sumfinevol = sumfinevol + finevol
      coarsesum = coarsesum + finevol*fine(iff,jff)
      
      enddo
      enddo
      coarsevol = (icc + half)*dxcoar*dxcoar*dxcoar
      coarse(icc,jcc) = coarsesum/coarsevol
      
      enddo
      enddo
      return
      end
      subroutine EBAVERAGEFACE(
     &           coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,icoarboxlo0,icoarboxlo1
     &           ,icoarboxhi0,icoarboxhi1
     &           ,refrat
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           ,idir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1)
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer icoarboxlo0,icoarboxlo1
      integer icoarboxhi0,icoarboxhi1
      integer refrat
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer idir
      integer idoloop, jdoloop
      integer icc, jcc
      integer iff, jff
      integer ii, jj
      real_t refscale, coarsesum, frefine
      
      idoloop = 1-CHF_ID(idir,0)
      jdoloop = 1-CHF_ID(idir,1)
      refscale = one
      frefine = refrat
      do ii = 1, (CH_SPACEDIM - 1)
         refscale = refscale/frefine
      enddo
      
      do jcc = icoarboxlo1,icoarboxhi1
      do icc = icoarboxlo0,icoarboxhi0

      coarsesum = zero
#if (CH_SPACEDIM == 3)
      do kk = 0,(refrat-1)*kdoloop
#endif
         do jj = 0,(refrat-1)*jdoloop
            do ii = 0,(refrat-1)*idoloop
               
               iff = icc*refrat + ii
               jff = jcc*refrat + jj
               coarsesum = coarsesum + fine(iff,jff)
            enddo
         enddo
#if (CH_SPACEDIM == 3)
      enddo
#endif
      coarse(icc,jcc) = coarsesum * refscale
      
      enddo
      enddo
      return
      end
      subroutine EBCOARSEN(
     &           coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,h2laplfine
     &           ,ih2laplfinelo0,ih2laplfinelo1
     &           ,ih2laplfinehi0,ih2laplfinehi1
     &           ,icoarboxlo0,icoarboxlo1
     &           ,icoarboxhi0,icoarboxhi1
     &           ,refrat
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1)
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer ih2laplfinelo0,ih2laplfinelo1
      integer ih2laplfinehi0,ih2laplfinehi1
      REAL_T h2laplfine(
     &           ih2laplfinelo0:ih2laplfinehi0,
     &           ih2laplfinelo1:ih2laplfinehi1)
      integer icoarboxlo0,icoarboxlo1
      integer icoarboxhi0,icoarboxhi1
      integer refrat
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer icc, jcc
      integer iff, jff
      integer ii, jj
      real_t refscale,coarsesum,frefine,coarseave,laplsum, laplave
      refscale = one
      frefine = refrat
      do ii = 1, CH_SPACEDIM
         refscale = refscale/2
      enddo
      
      do jcc = icoarboxlo1,icoarboxhi1
      do icc = icoarboxlo0,icoarboxhi0

      laplsum   = zero
      coarsesum = zero
      
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0

      
      iff = icc*refrat + ii
      jff = jcc*refrat + jj
      coarsesum = coarsesum +       fine(iff,jff)
      laplsum   = laplsum   + h2laplfine(iff,jff)
      
      enddo
      enddo
      laplave   =   laplsum * refscale
      coarseave = coarsesum * refscale
      coarse(icc,jcc) = coarseave - laplave/eight
      
      enddo
      enddo
      return
      end
      subroutine H2LAPL1DADDITIVE(
     &           opphidir
     &           ,iopphidirlo0,iopphidirlo1
     &           ,iopphidirhi0,iopphidirhi1
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,idir
     &           ,iloboxlo0,iloboxlo1
     &           ,iloboxhi0,iloboxhi1
     &           ,haslo
     &           ,ihiboxlo0,ihiboxlo1
     &           ,ihiboxhi0,ihiboxhi1
     &           ,hashi
     &           ,icenterboxlo0,icenterboxlo1
     &           ,icenterboxhi0,icenterboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer iopphidirlo0,iopphidirlo1
      integer iopphidirhi0,iopphidirhi1
      REAL_T opphidir(
     &           iopphidirlo0:iopphidirhi0,
     &           iopphidirlo1:iopphidirhi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      integer idir
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer haslo
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer hashi
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer i,j
      integer ioff,joff
      
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

      opphidir(i,j) = opphidir(i,j) +
     $     (    phi(i+ioff,j+joff)
     &     -two*phi(i     ,j     )
     &     +    phi(i-ioff,j-joff))
      
      enddo
      enddo
      if (haslo .eq. 1) then
         
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

         opphidir(i,j) = opphidir(i,j) +
     $        (    phi(i       ,j       )
     &        -two*phi(i+  ioff,j+  joff)
     &        +    phi(i+2*ioff,j+2*joff))
         
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
         
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

         opphidir(i,j) = opphidir(i,j) +
     $        (    phi(i       ,j       )
     &        -two*phi(i-  ioff,j-  joff)
     &        +    phi(i-2*ioff,j-2*joff))
         
      enddo
      enddo
      endif
      return
      end
