      subroutine HALVEINT(
     & arr
     & ,iarrlo0,iarrlo1
     & ,iarrhi0,iarrhi1
     & ,narrcomp
     & ,ibxlo0,ibxlo1
     & ,ibxhi0,ibxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer narrcomp
      integer iarrlo0,iarrlo1
      integer iarrhi0,iarrhi1
      integer arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1,
     & 0:narrcomp-1)
      integer ibxlo0,ibxlo1
      integer ibxhi0,ibxhi1
      integer i0,i1
      integer var
      do var = 0, narrcomp-1
      do i1 = ibxlo1,ibxhi1
      do i0 = ibxlo0,ibxhi0
            arr(i0,i1, var) = arr(i0,i1, var) / 2
      enddo
      enddo
      enddo
      return
      end
