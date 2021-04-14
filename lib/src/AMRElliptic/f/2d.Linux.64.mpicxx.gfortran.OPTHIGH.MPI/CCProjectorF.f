      subroutine CCPEXTRAPTODOMFACE(
     & facevel
     & ,ifacevello0,ifacevello1
     & ,ifacevelhi0,ifacevelhi1
     & ,facedir
     & ,ishift
     & ,ifaceboxlo0,ifaceboxlo1
     & ,ifaceboxhi0,ifaceboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacevello0,ifacevello1
      integer ifacevelhi0,ifacevelhi1
      REAL*8 facevel(
     & ifacevello0:ifacevelhi0,
     & ifacevello1:ifacevelhi1)
      integer facedir
      integer ishift
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      integer i,j
      integer ioff,joff
      ioff = ishift*chf_id(0,facedir)
      joff = ishift*chf_id(1,facedir)
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      facevel(i,j) =
     & ((2.0d0)*facevel(i+ ioff,j+ joff)
     & - facevel(i+2*ioff,j+2*joff))
      enddo
      enddo
      ch_flops=ch_flops+(ifaceboxhi0- ifaceboxlo0+1)*(ifaceboxhi1- iface
     &boxlo1+1)
      return
      end
      subroutine CCPAVECELLTOFACE(
     & facevel
     & ,ifacevello0,ifacevello1
     & ,ifacevelhi0,ifacevelhi1
     & ,cellvel
     & ,icellvello0,icellvello1
     & ,icellvelhi0,icellvelhi1
     & ,facedir
     & ,ifaceboxlo0,ifaceboxlo1
     & ,ifaceboxhi0,ifaceboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacevello0,ifacevello1
      integer ifacevelhi0,ifacevelhi1
      REAL*8 facevel(
     & ifacevello0:ifacevelhi0,
     & ifacevello1:ifacevelhi1)
      integer icellvello0,icellvello1
      integer icellvelhi0,icellvelhi1
      REAL*8 cellvel(
     & icellvello0:icellvelhi0,
     & icellvello1:icellvelhi1)
      integer facedir
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      integer i,j
      integer ioff,joff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      facevel(i,j) =
     & ( cellvel(i ,j )
     & + cellvel(i-ioff,j-joff)
     & )*(0.500d0)
      enddo
      enddo
      ch_flops=ch_flops+(ifaceboxhi0- ifaceboxlo0+1)*(ifaceboxhi1- iface
     &boxlo1+1)*2
      return
      end
      subroutine CCPAVEFACETOCELL(
     & cellgrad
     & ,icellgradlo0,icellgradlo1
     & ,icellgradhi0,icellgradhi1
     & ,facegrad
     & ,ifacegradlo0,ifacegradlo1
     & ,ifacegradhi0,ifacegradhi1
     & ,facedir
     & ,icellboxlo0,icellboxlo1
     & ,icellboxhi0,icellboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer icellgradlo0,icellgradlo1
      integer icellgradhi0,icellgradhi1
      REAL*8 cellgrad(
     & icellgradlo0:icellgradhi0,
     & icellgradlo1:icellgradhi1)
      integer ifacegradlo0,ifacegradlo1
      integer ifacegradhi0,ifacegradhi1
      REAL*8 facegrad(
     & ifacegradlo0:ifacegradhi0,
     & ifacegradlo1:ifacegradhi1)
      integer facedir
      integer icellboxlo0,icellboxlo1
      integer icellboxhi0,icellboxhi1
      integer i, j
      integer ioff, joff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      do j = icellboxlo1,icellboxhi1
      do i = icellboxlo0,icellboxhi0
      cellgrad(i,j) =
     & (facegrad(i+ioff,j+joff)
     & +facegrad(i ,j ))*(0.500d0)
      enddo
      enddo
      ch_flops=ch_flops+(icellboxhi0- icellboxlo0+1)*(icellboxhi1- icell
     &boxlo1+1)*2
      return
      end
      subroutine CCPCELLGRADFROMFACEDATA(
     & cellgrad
     & ,icellgradlo0,icellgradlo1
     & ,icellgradhi0,icellgradhi1
     & ,facedata
     & ,ifacedatalo0,ifacedatalo1
     & ,ifacedatahi0,ifacedatahi1
     & ,facedir
     & ,dx
     & ,icellboxlo0,icellboxlo1
     & ,icellboxhi0,icellboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer icellgradlo0,icellgradlo1
      integer icellgradhi0,icellgradhi1
      REAL*8 cellgrad(
     & icellgradlo0:icellgradhi0,
     & icellgradlo1:icellgradhi1)
      integer ifacedatalo0,ifacedatalo1
      integer ifacedatahi0,ifacedatahi1
      REAL*8 facedata(
     & ifacedatalo0:ifacedatahi0,
     & ifacedatalo1:ifacedatahi1)
      integer facedir
      REAL*8 dx
      integer icellboxlo0,icellboxlo1
      integer icellboxhi0,icellboxhi1
      integer i, j
      integer ioff, joff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      do j = icellboxlo1,icellboxhi1
      do i = icellboxlo0,icellboxhi0
      cellgrad(i,j) =
     & (facedata(i+ioff,j+joff)
     & -facedata(i ,j ))/dx
      enddo
      enddo
      ch_flops=ch_flops+(icellboxhi0- icellboxlo0+1)*(icellboxhi1- icell
     &boxlo1+1)*2
      return
      end
