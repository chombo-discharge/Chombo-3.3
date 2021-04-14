#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine DOTPRODUCT(
     &           dotprodout
     &           ,afab
     &           ,iafablo0,iafablo1
     &           ,iafabhi0,iafabhi1
     &           ,nafabcomp
     &           ,bfab
     &           ,ibfablo0,ibfablo1
     &           ,ibfabhi0,ibfabhi1
     &           ,nbfabcomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T dotprodout
      integer nafabcomp
      integer iafablo0,iafablo1
      integer iafabhi0,iafabhi1
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           0:nafabcomp-1)
      integer nbfabcomp
      integer ibfablo0,ibfablo1
      integer ibfabhi0,ibfabhi1
      REAL_T bfab(
     &           ibfablo0:ibfabhi0,
     &           ibfablo1:ibfabhi1,
     &           0:nbfabcomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer startcomp
      integer endcomp
      integer i0,i1
      integer nv
      dotprodout = 0
      do nv=startcomp,endcomp,1
         
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

         dotprodout = dotprodout +
     &        afab(i0,i1,nv)*
     &        bfab(i0,i1,nv)
         
      enddo
      enddo
      enddo
      return
      end
