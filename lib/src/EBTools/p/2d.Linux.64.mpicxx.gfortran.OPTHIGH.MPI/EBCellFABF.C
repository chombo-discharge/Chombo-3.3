#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine SETCOVERED(
     &           vel
     &           ,ivello0,ivello1
     &           ,ivelhi0,ivelhi1
     &           ,m
     &           ,mask
     &           ,imasklo0,imasklo1
     &           ,imaskhi0,imaskhi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1)
      REAL_T m
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      BYTE_T mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      if ( m .eq. 0) then
       
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

         vel(i,j)= vel(i,j)*mask(i,j)
       
      enddo
      enddo
       ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)
      else
       
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
 
         vel(i,j)= vel(i,j)*mask(i,j)
     &                     + (1-mask(i,j))*m
       
      enddo
      enddo
       ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*4
      endif
      return
      end
