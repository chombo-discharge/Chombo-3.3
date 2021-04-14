#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine NODEINTERPMG_GETWEIGHTS(
     &           nref
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           ,wtcrnr
     &           ,iwtcrnrlo0,iwtcrnrlo1
     &           ,iwtcrnrhi0,iwtcrnrhi1
     &           ,nwtcrnrcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nref
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer nwtcrnrcomp
      integer iwtcrnrlo0,iwtcrnrlo1
      integer iwtcrnrhi0,iwtcrnrhi1
      REAL_T wtcrnr(
     &           iwtcrnrlo0:iwtcrnrhi0,
     &           iwtcrnrlo1:iwtcrnrhi1,
     &           0:nwtcrnrcomp-1)
      integer iref0,iref1
      integer ib0,ib1
      integer nvwt
      integer ibmax0,ibmax1
      REAL_T refinv, wt
      REAL_T fraci0,fraci1
      REAL_T wti0,wti1
      refinv = one / nref
      nvwt = 0
      
      do iref1 = ibreflo1,ibrefhi1
      do iref0 = ibreflo0,ibrefhi0

         
         call maxb(iref0, ibmax0) 
         call maxb(iref1, ibmax1) 
         
         fraci0 = iref0 * refinv 
         fraci1 = iref1 * refinv 
         
         do ib0 = 0, ibmax0 
         do ib1 = 0, ibmax1 
            wt = one
            
            call wtside(ib0, fraci0, wti0)
            wt = wt * wti0 
            call wtside(ib1, fraci1, wti1)
            wt = wt * wti1 
            wtcrnr( ib0,ib1, nvwt ) = wt
         
         end do 
         end do 
         nvwt = nvwt + 1
      
      enddo
      enddo
      return
      end
      subroutine NODEINTERPMG(
     &           fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,nfinecomp
     &           ,crse
     &           ,icrselo0,icrselo1
     &           ,icrsehi0,icrsehi1
     &           ,ncrsecomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,nref
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           ,wtcrnr
     &           ,iwtcrnrlo0,iwtcrnrlo1
     &           ,iwtcrnrhi0,iwtcrnrhi1
     &           ,nwtcrnrcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           0:nfinecomp-1)
      integer ncrsecomp
      integer icrselo0,icrselo1
      integer icrsehi0,icrsehi1
      REAL_T crse(
     &           icrselo0:icrsehi0,
     &           icrselo1:icrsehi1,
     &           0:ncrsecomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer nref
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer nwtcrnrcomp
      integer iwtcrnrlo0,iwtcrnrlo1
      integer iwtcrnrhi0,iwtcrnrhi1
      REAL_T wtcrnr(
     &           iwtcrnrlo0:iwtcrnrhi0,
     &           iwtcrnrlo1:iwtcrnrhi1,
     &           0:nwtcrnrcomp-1)
      integer iref0,iref1, icrse0,icrse1
      integer ifine0,ifine1, ib0,ib1;
      integer var, ncomp, nvwt
      integer ibmax0,ibmax1
      integer icmin0,icmin1
      integer icmax0,icmax1;
      REAL_T csum, refinv
      ncomp = nfinecomp
      if (ncomp .ne. ncrsecomp) then
         print *, 'nodeinterpmg: fine and crse incompatible'
         call MAYDAY_ERROR()
      endif
      refinv = one / nref
      
      icmin0 = iregionlo0 
      icmin1 = iregionlo1 
      nvwt = 0
      
      do iref1 = ibreflo1,ibrefhi1
      do iref0 = ibreflo0,ibrefhi0

         
         call maxb(iref0, ibmax0) 
         call maxb(iref1, ibmax1) 
         
         icmax0 = iregionhi0 + (1-ibmax0) 
         icmax1 = iregionhi1 + (1-ibmax1) 
         
         do icrse0 = icmin0, icmax0 
         do icrse1 = icmin1, icmax1 
            
            ifine0 = nref*icrse0 + iref0 
            ifine1 = nref*icrse1 + iref1 
            do var = 0, ncomp-1
               csum = 0
               
               do ib0 = 0, ibmax0 
               do ib1 = 0, ibmax1 
                  csum = csum + wtcrnr( ib0,ib1, nvwt ) *
     &              crse( icrse0+ib0,icrse1+ib1, var)
               
               end do 
               end do 
               fine ( ifine0,ifine1, var ) = csum +
     &              fine ( ifine0,ifine1, var )
            end do
         
         end do 
         end do 
         nvwt = nvwt + 1
      
      enddo
      enddo
      return
      end
      subroutine WTSIDE(
     &           i
     &           ,frac
     &           ,wt
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer i
      REAL_T frac
      REAL_T wt
      if (i .eq. 0) then
         wt = one - frac
      else
         wt = frac
      endif
      return
      end
      subroutine MAXB(
     &           iref
     &           ,ibmax
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iref
      integer ibmax
      if (iref .eq. 0) then
         ibmax = 0
      else
         ibmax = 1
      endif
      return
      end
