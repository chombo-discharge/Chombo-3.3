      subroutine NODEOPLAP(
     & lofphi
     & ,ilofphilo0,ilofphilo1
     & ,ilofphihi0,ilofphihi1
     & ,nlofphicomp
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlofphicomp
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL*8 lofphi(
     & ilofphilo0:ilofphihi0,
     & ilofphilo1:ilofphihi1,
     & 0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 dxinv2, lphi
      integer var, ncomp
      integer i, j
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = (1.0d0) / (dx*dx)
      do var = 0, ncomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
            lphi = ( (phi(i+1 , j , var)
     & - phi(i , j , var) )
     & - (phi(i , j , var)
     & - phi(i-1 , j , var) ) ) * dxinv2
     & + ( (phi(i , j+1, var)
     & - phi(i , j , var) )
     & - (phi(i , j , var)
     & - phi(i , j-1, var) ) ) * dxinv2
            lofphi(i, j, var) = lphi
      enddo
      enddo
      end do
      return
      end
      subroutine NODEOPLAPPOINT(
     & lofphi
     & ,ilofphilo0,ilofphilo1
     & ,ilofphihi0,ilofphihi1
     & ,nlofphicomp
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,pt
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlofphicomp
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL*8 lofphi(
     & ilofphilo0:ilofphihi0,
     & ilofphilo1:ilofphihi1,
     & 0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer pt(0:1)
      REAL*8 dx
      REAL*8 dxinv2, lphi
      integer var, ncomp
      integer i, j
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = (1.0d0) / (dx*dx)
      i = pt(0)
      j = pt(1)
      do var = 0, ncomp-1
         lphi = ( (phi(i+1 , j , var)
     & - phi(i , j , var) )
     & - (phi(i , j , var)
     & - phi(i-1 , j , var) ) ) * dxinv2
     & + ( (phi(i , j+1, var)
     & - phi(i , j , var) )
     & - (phi(i , j , var)
     & - phi(i , j-1, var) ) ) * dxinv2
         lofphi(i, j, var) = lphi
      end do
      return
      end
      subroutine NODEGRAD(
     & grdphi
     & ,igrdphilo0,igrdphilo1
     & ,igrdphihi0,igrdphihi1
     & ,ngrdphicomp
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ngrdphicomp
      integer igrdphilo0,igrdphilo1
      integer igrdphihi0,igrdphihi1
      REAL*8 grdphi(
     & igrdphilo0:igrdphihi0,
     & igrdphilo1:igrdphihi1,
     & 0:ngrdphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 dxinvh
      integer var, ncomp, gbase
      integer i, j
      ncomp = nphicomp
      if (2 * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = (0.500d0) / dx
      do var = 0, ncomp-1
         gbase = 2 * var
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
            grdphi(i, j, gbase) =
     & ( phi(i+1 , j , var)
     & - phi(i-1 , j , var) ) * dxinvh
            grdphi(i, j, gbase + 1) =
     & ( phi(i , j+1 , var)
     & - phi(i , j-1 , var) ) * dxinvh
      enddo
      enddo
      end do
      return
      end
      subroutine NODEGRADPOINT(
     & grdphi
     & ,igrdphilo0,igrdphilo1
     & ,igrdphihi0,igrdphihi1
     & ,ngrdphicomp
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,pt
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ngrdphicomp
      integer igrdphilo0,igrdphilo1
      integer igrdphihi0,igrdphihi1
      REAL*8 grdphi(
     & igrdphilo0:igrdphihi0,
     & igrdphilo1:igrdphihi1,
     & 0:ngrdphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer pt(0:1)
      REAL*8 dx
      REAL*8 dxinvh
      integer var, ncomp, gbase
      integer i, j
      ncomp = nphicomp
      if (2 * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = (0.500d0) / dx
      i = pt(0)
      j = pt(1)
      do var = 0, ncomp-1
         gbase = 2 * var
            grdphi(i, j, gbase) =
     & ( phi(i+1 , j , var)
     & - phi(i-1 , j , var) ) * dxinvh
            grdphi(i, j, gbase + 1) =
     & ( phi(i , j+1 , var)
     & - phi(i , j-1 , var) ) * dxinvh
      end do
      return
      end
      subroutine NODEGSRBLEVELLAP(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,rhs
     & ,irhslo0,irhslo1
     & ,irhshi0,irhshi1
     & ,nrhscomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,dx
     & ,redBlack
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & 0:nrhscomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      integer redBlack
      REAL*8 lambda
      REAL*8 dxinv2, lphi
      integer i, j
      integer imin, imax, var, ncomp, indtot
      dxinv2 = (1.0d0)/(dx*dx)
      lambda = (dx*dx) / ((2.0d0)*2)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAY_ERROR()
      endif
      do var = 0, ncomp - 1
            do j = iregionlo1, iregionhi1
               imin = iregionlo0
               indtot = imin + j
               imin = imin + mod(indtot + redBlack, 2)
               imax = iregionhi0
               do i = imin, imax, 2
            lphi = ( (phi(i+1 , j , var)
     & + phi(i-1 , j , var) ) ) * dxinv2
     & + ( (phi(i , j+1, var)
     & + phi(i , j-1, var) ) ) * dxinv2
                  phi(i,j, var) =
     & lambda * (lphi - rhs(i,j, var))
               enddo
            enddo
      enddo
      return
      end
