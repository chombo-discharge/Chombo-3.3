#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine MASKVALUE(
     &           mask
     &           ,imasklo0,imasklo1
     &           ,imaskhi0,imaskhi1
     &           ,test
     &           ,itestlo0,itestlo1
     &           ,itesthi0,itesthi1
     &           ,ibxlo0,ibxlo1
     &           ,ibxhi0,ibxhi1
     &           ,val
     &           ,onoff
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      REAL_T mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1)
      integer itestlo0,itestlo1
      integer itesthi0,itesthi1
      integer test(
     &           itestlo0:itesthi0,
     &           itestlo1:itesthi1)
      integer ibxlo0,ibxlo1
      integer ibxhi0,ibxhi1
      integer val
      integer onoff
      integer i0,i1
      if (onoff .eq. 1) then
         
      do i1 = ibxlo1,ibxhi1
      do i0 = ibxlo0,ibxhi0

         if (test(i0,i1) .eq. val) then
            mask(i0,i1) = one
         else
            mask(i0,i1) = 0
         endif
         
      enddo
      enddo
      else
         
      do i1 = ibxlo1,ibxhi1
      do i0 = ibxlo0,ibxhi0

         if (test(i0,i1) .eq. val) then
            mask(i0,i1) = 0
         else
            mask(i0,i1) = one
         endif
         
      enddo
      enddo
      endif
      return
      end
