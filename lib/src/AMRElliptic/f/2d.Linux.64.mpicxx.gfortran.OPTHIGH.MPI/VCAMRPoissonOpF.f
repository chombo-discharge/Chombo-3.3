      subroutine GSRBHELMHOLTZVC2D(
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
     & ,alpha
     & ,aCoef
     & ,iaCoeflo0,iaCoeflo1
     & ,iaCoefhi0,iaCoefhi1
     & ,naCoefcomp
     & ,beta
     & ,bCoef0
     & ,ibCoef0lo0,ibCoef0lo1
     & ,ibCoef0hi0,ibCoef0hi1
     & ,nbCoef0comp
     & ,bCoef1
     & ,ibCoef1lo0,ibCoef1lo1
     & ,ibCoef1hi0,ibCoef1hi1
     & ,nbCoef1comp
     & ,lambda
     & ,ilambdalo0,ilambdalo1
     & ,ilambdahi0,ilambdahi1
     & ,nlambdacomp
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
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0,iaCoeflo1
      integer iaCoefhi0,iaCoefhi1
      REAL*8 aCoef(
     & iaCoeflo0:iaCoefhi0,
     & iaCoeflo1:iaCoefhi1,
     & 0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1
      integer ibCoef0hi0,ibCoef0hi1
      REAL*8 bCoef0(
     & ibCoef0lo0:ibCoef0hi0,
     & ibCoef0lo1:ibCoef0hi1,
     & 0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1
      integer ibCoef1hi0,ibCoef1hi1
      REAL*8 bCoef1(
     & ibCoef1lo0:ibCoef1hi0,
     & ibCoef1lo1:ibCoef1hi1,
     & 0:nbCoef1comp-1)
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1
      integer ilambdahi0,ilambdahi1
      REAL*8 lambda(
     & ilambdalo0:ilambdahi0,
     & ilambdalo1:ilambdahi1,
     & 0:nlambdacomp-1)
      integer redBlack
      REAL*8 dxinv,lofphi
      integer n,ncomp,indtot,imin,imax
      integer i,j
      ncomp = nphicomp
      if (ncomp .ne. nphicomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef1comp) then
         call MAYDAYERROR()
      endif
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp - 1
          do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
              lofphi =
     & alpha * aCoef(i,j,n) * phi(i,j,n)
     & - beta *
     & (
     & bCoef0(i+1,j ,n)
     & * (phi(i+1,j ,n)-phi(i ,j ,n))
     & - bCoef0(i ,j ,n)
     & * (phi(i ,j ,n)-phi(i-1,j ,n))
     & + bCoef1(i ,j+1,n)
     & * (phi(i ,j+1,n)-phi(i ,j ,n))
     & - bCoef1(i ,j ,n)
     & * (phi(i ,j ,n)-phi(i ,j-1,n))
     & ) * dxinv
              phi(i,j,n) = phi(i,j,n)
     & - lambda(i,j,n) * (lofphi - rhs(i,j,n))
            enddo
          enddo
      enddo
      return
      end
      subroutine VCCOMPUTEOP2D(
     & lofphi
     & ,ilofphilo0,ilofphilo1
     & ,ilofphihi0,ilofphihi1
     & ,nlofphicomp
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,alpha
     & ,aCoef
     & ,iaCoeflo0,iaCoeflo1
     & ,iaCoefhi0,iaCoefhi1
     & ,naCoefcomp
     & ,beta
     & ,bCoef0
     & ,ibCoef0lo0,ibCoef0lo1
     & ,ibCoef0hi0,ibCoef0hi1
     & ,nbCoef0comp
     & ,bCoef1
     & ,ibCoef1lo0,ibCoef1lo1
     & ,ibCoef1hi0,ibCoef1hi1
     & ,nbCoef1comp
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
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0,iaCoeflo1
      integer iaCoefhi0,iaCoefhi1
      REAL*8 aCoef(
     & iaCoeflo0:iaCoefhi0,
     & iaCoeflo1:iaCoefhi1,
     & 0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1
      integer ibCoef0hi0,ibCoef0hi1
      REAL*8 bCoef0(
     & ibCoef0lo0:ibCoef0hi0,
     & ibCoef0lo1:ibCoef0hi1,
     & 0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1
      integer ibCoef1hi0,ibCoef1hi1
      REAL*8 bCoef1(
     & ibCoef1lo0:ibCoef1hi0,
     & ibCoef1lo1:ibCoef1hi1,
     & 0:nbCoef1comp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 dxinv
      integer n,ncomp
      integer i,j
      ncomp = nphicomp
      if (ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef1comp) then
         call MAYDAYERROR()
      endif
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          lofphi(i,j,n) =
     & alpha * aCoef(i,j,n) * phi(i,j,n)
     & - beta *
     & (
     & bCoef0(i+1,j ,n)
     & * (phi(i+1,j ,n) - phi(i ,j ,n))
     & - bCoef0(i ,j ,n)
     & * (phi(i ,j ,n) - phi(i-1,j ,n))
     & + bCoef1(i ,j+1,n)
     & * (phi(i ,j+1,n) - phi(i ,j ,n))
     & - bCoef1(i ,j ,n)
     & * (phi(i ,j ,n) - phi(i ,j-1,n))
     & ) * dxinv
      enddo
      enddo
      enddo
      return
      end
      subroutine VCCOMPUTERES2D(
     & res
     & ,ireslo0,ireslo1
     & ,ireshi0,ireshi1
     & ,nrescomp
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,rhs
     & ,irhslo0,irhslo1
     & ,irhshi0,irhshi1
     & ,nrhscomp
     & ,alpha
     & ,aCoef
     & ,iaCoeflo0,iaCoeflo1
     & ,iaCoefhi0,iaCoefhi1
     & ,naCoefcomp
     & ,beta
     & ,bCoef0
     & ,ibCoef0lo0,ibCoef0lo1
     & ,ibCoef0hi0,ibCoef0hi1
     & ,nbCoef0comp
     & ,bCoef1
     & ,ibCoef1lo0,ibCoef1lo1
     & ,ibCoef1hi0,ibCoef1hi1
     & ,nbCoef1comp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0,ireslo1
      integer ireshi0,ireshi1
      REAL*8 res(
     & ireslo0:ireshi0,
     & ireslo1:ireshi1,
     & 0:nrescomp-1)
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
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0,iaCoeflo1
      integer iaCoefhi0,iaCoefhi1
      REAL*8 aCoef(
     & iaCoeflo0:iaCoefhi0,
     & iaCoeflo1:iaCoefhi1,
     & 0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1
      integer ibCoef0hi0,ibCoef0hi1
      REAL*8 bCoef0(
     & ibCoef0lo0:ibCoef0hi0,
     & ibCoef0lo1:ibCoef0hi1,
     & 0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1
      integer ibCoef1hi0,ibCoef1hi1
      REAL*8 bCoef1(
     & ibCoef1lo0:ibCoef1hi0,
     & ibCoef1lo1:ibCoef1hi1,
     & 0:nbCoef1comp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 dxinv
      integer n,ncomp
      integer i,j
      ncomp = nphicomp
      if (ncomp .ne. nrescomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef1comp) then
         call MAYDAYERROR()
      endif
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          res(i,j,n) =
     & rhs(i,j,n)
     & - (alpha * aCoef(i,j,n) * phi(i,j,n)
     & - beta *
     & (
     & bCoef0(i+1,j ,n)
     & * (phi(i+1,j ,n) - phi(i ,j ,n))
     & - bCoef0(i ,j ,n)
     & * (phi(i ,j ,n) - phi(i-1,j ,n))
     & + bCoef1(i ,j+1,n)
     & * (phi(i ,j+1,n) - phi(i ,j ,n))
     & - bCoef1(i ,j ,n)
     & * (phi(i ,j ,n) - phi(i ,j-1,n))
     & ) * dxinv
     & )
      enddo
      enddo
      enddo
      return
      end
      subroutine RESTRICTRESVC2D(
     & res
     & ,ireslo0,ireslo1
     & ,ireshi0,ireshi1
     & ,nrescomp
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,rhs
     & ,irhslo0,irhslo1
     & ,irhshi0,irhshi1
     & ,nrhscomp
     & ,alpha
     & ,aCoef
     & ,iaCoeflo0,iaCoeflo1
     & ,iaCoefhi0,iaCoefhi1
     & ,naCoefcomp
     & ,beta
     & ,bCoef0
     & ,ibCoef0lo0,ibCoef0lo1
     & ,ibCoef0hi0,ibCoef0hi1
     & ,nbCoef0comp
     & ,bCoef1
     & ,ibCoef1lo0,ibCoef1lo1
     & ,ibCoef1hi0,ibCoef1hi1
     & ,nbCoef1comp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0,ireslo1
      integer ireshi0,ireshi1
      REAL*8 res(
     & ireslo0:ireshi0,
     & ireslo1:ireshi1,
     & 0:nrescomp-1)
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
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0,iaCoeflo1
      integer iaCoefhi0,iaCoefhi1
      REAL*8 aCoef(
     & iaCoeflo0:iaCoefhi0,
     & iaCoeflo1:iaCoefhi1,
     & 0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1
      integer ibCoef0hi0,ibCoef0hi1
      REAL*8 bCoef0(
     & ibCoef0lo0:ibCoef0hi0,
     & ibCoef0lo1:ibCoef0hi1,
     & 0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1
      integer ibCoef1hi0,ibCoef1hi1
      REAL*8 bCoef1(
     & ibCoef1lo0:ibCoef1hi0,
     & ibCoef1lo1:ibCoef1hi1,
     & 0:nbCoef1comp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 denom,dxinv,lofphi
      integer n,ncomp
      integer i,j
      integer ii,jj
      ncomp = nphicomp
      dxinv = (1.0d0) / (dx*dx)
      denom = 2 *2
      do n = 0, ncomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          ii = i/2
          jj = j/2
          lofphi =
     & alpha * aCoef(i,j,n) * phi(i,j,n)
     & - beta *
     & (
     & bCoef0(i+1,j ,n)
     & * (phi(i+1,j ,n)-phi(i ,j ,n))
     & - bCoef0(i ,j ,n)
     & * (phi(i ,j ,n)-phi(i-1,j ,n))
     & + bCoef1(i ,j+1,n)
     & * (phi(i ,j+1,n)-phi(i ,j ,n))
     & - bCoef1(i ,j ,n)
     & * (phi(i ,j ,n)-phi(i ,j-1,n))
     & ) * dxinv
          res(ii,jj,n) = res(ii,jj,n)
     & + (rhs(i,j,n) - lofphi) / denom
      enddo
      enddo
      enddo
      return
      end
      subroutine SUMFACES(
     & lhs
     & ,ilhslo0,ilhslo1
     & ,ilhshi0,ilhshi1
     & ,nlhscomp
     & ,beta
     & ,bCoefs
     & ,ibCoefslo0,ibCoefslo1
     & ,ibCoefshi0,ibCoefshi1
     & ,nbCoefscomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dir
     & ,scale
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nlhscomp
      integer ilhslo0,ilhslo1
      integer ilhshi0,ilhshi1
      REAL*8 lhs(
     & ilhslo0:ilhshi0,
     & ilhslo1:ilhshi1,
     & 0:nlhscomp-1)
      REAL*8 beta
      integer nbCoefscomp
      integer ibCoefslo0,ibCoefslo1
      integer ibCoefshi0,ibCoefshi1
      REAL*8 bCoefs(
     & ibCoefslo0:ibCoefshi0,
     & ibCoefslo1:ibCoefshi1,
     & 0:nbCoefscomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      REAL*8 scale
      REAL*8 sumVal
      integer i,j
      integer ii,jj
      integer n
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      do n = 0, nlhscomp-1
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          sumVal = bCoefs(i+ii,j+jj,n)
     & + bCoefs(i ,j ,n)
          lhs(i,j,n) = lhs(i,j,n) + scale * beta * sumVal
      enddo
      enddo
      enddo
      return
      end
