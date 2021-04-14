#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine OPERATORLAP4(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1
     &           ,ilofphihi0,ilofphihi1
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlofphicomp
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1,
     &           0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      REAL_T dxinv, lap
      integer n,ncomp
      integer i,j
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
        
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          lap = ( 
     &   sixteen*phi(i-1,j  ,n) - phi(i-2,j  ,n)
     & + sixteen*phi(i+1,j  ,n) - phi(i+2,j  ,n)
     & + sixteen*phi(i  ,j-1,n) - phi(i  ,j-2,n)
     & + sixteen*phi(i  ,j+1,n) - phi(i  ,j+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,j,n) )
     &       * twelfth * dxinv
          lofphi(i,j,n) = alpha*phi(i,j,n)+beta*lap
        
      enddo
      enddo
      enddo
      return
      end
      subroutine OPERATORLAPRES4(
     &           r
     &           ,irlo0,irlo1
     &           ,irhi0,irhi1
     &           ,nrcomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrcomp
      integer irlo0,irlo1
      integer irhi0,irhi1
      REAL_T r(
     &           irlo0:irhi0,
     &           irlo1:irhi1,
     &           0:nrcomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      REAL_T dxinv, lap, lhs
      integer n,ncomp
      integer i,j
      ncomp = nphicomp
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
         
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          lap = ( 
     &   sixteen*phi(i-1,j  ,n) - phi(i-2,j  ,n)
     & + sixteen*phi(i+1,j  ,n) - phi(i+2,j  ,n)
     & + sixteen*phi(i  ,j-1,n) - phi(i  ,j-2,n)
     & + sixteen*phi(i  ,j+1,n) - phi(i  ,j+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,j,n) )
     &       * twelfth * dxinv
          lhs = alpha*phi(i,j,n) + beta*lap
          r(i,j,n) = rhs(i,j,n) - lhs
         
      enddo
      enddo
      enddo
      return
      end
      subroutine RESTRICTRES4(
     &           res
     &           ,ireslo0,ireslo1
     &           ,ireshi0,ireshi1
     &           ,nrescomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,alpha
     &           ,beta
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0,ireslo1
      integer ireshi0,ireshi1
      REAL_T res(
     &           ireslo0:ireshi0,
     &           ireslo1:ireshi1,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      REAL_T alpha
      REAL_T beta
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T dx
      REAL_T denom,dxinv,lofphi,lap
      integer n,ncomp
      integer i,j
      integer ii,jj
      ncomp = nphicomp
      dxinv = one / (dx*dx)
      denom = D_TERM(2, *2, *2)
      do n = 0, ncomp-1
        
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          
          ii = i/2 
          jj = j/2 
          lap = ( 
     &   sixteen*phi(i-1,j  ,n) - phi(i-2,j  ,n)
     & + sixteen*phi(i+1,j  ,n) - phi(i+2,j  ,n)
     & + sixteen*phi(i  ,j-1,n) - phi(i  ,j-2,n)
     & + sixteen*phi(i  ,j+1,n) - phi(i  ,j+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,j,n) )
     &       * twelfth * dxinv
          lofphi = alpha*phi(i,j,n) + beta*lap
          res(ii,jj,n) = res(ii,jj,n)
     &                            + (rhs(i,j,n) - lofphi) / denom
        
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBLAPLACIAN4(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           ,tmp
     &           ,itmplo0,itmplo1
     &           ,itmphi0,itmphi1
     &           ,ntmpcomp
     &           ,redBlack
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T dx
      integer ntmpcomp
      integer itmplo0,itmplo1
      integer itmphi0,itmphi1
      REAL_T tmp(
     &           itmplo0:itmphi0,
     &           itmplo1:itmphi1,
     &           0:ntmpcomp-1)
      integer redBlack
      REAL_T dx2t, thD
      integer i,j
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = twelve*dx*dx
      thD  = thirtieth/CH_SPACEDIM
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               tmp(i,j,n) = thD*( 
     &           sixteen*phi(i+1,j,n) - phi(i+2,j,n)
     &         + sixteen*phi(i-1,j,n) - phi(i-2,j,n)
     &         + sixteen*phi(i,j+1,n) - phi(i,j+2,n)
     &         + sixteen*phi(i,j-1,n) - phi(i,j-2,n)
     &         - dx2t*rhs(i,j,n) )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else if (redBlack .eq. black) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,j,n) = thD*( 
     &           sixteen*tmp(i+1,j,n) - tmp(i+2,j,n)
     &         + sixteen*tmp(i-1,j,n) - tmp(i-2,j,n)
     &         + sixteen*tmp(i,j+1,n) - tmp(i,j+2,n)
     &         + sixteen*tmp(i,j-1,n) - tmp(i,j-2,n)
     &         - dx2t*rhs(i,j,n) )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,j,n) = tmp(i,j,n)
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else
           call MAYDAYERROR()
        end if
      enddo
      return
      end
      subroutine SORLAPLACIAN4(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T dx
      REAL_T dx2t, thD, tmp, omega
      integer i,j
      integer n,ncomp
      dx2t = twelve*dx*dx
      thD  = thirtieth/CH_SPACEDIM
      omega = 0.47
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         tmp = thD*( 
     &        sixteen*phi(i+1,j,n) - phi(i+2,j,n)
     &        + sixteen*phi(i-1,j,n) - phi(i-2,j,n)
     &        + sixteen*phi(i,j+1,n) - phi(i,j+2,n)
     &        + sixteen*phi(i,j-1,n) - phi(i,j-2,n)
     &        - dx2t*rhs(i,j,n) )
         phi(i,j,n) = omega*tmp
     &        + (one-omega)*phi(i,j,n)
         
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBHELMHOLTZ4(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           ,tmp
     &           ,itmplo0,itmplo1
     &           ,itmphi0,itmphi1
     &           ,ntmpcomp
     &           ,alpha
     &           ,beta
     &           ,redBlack
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T dx
      integer ntmpcomp
      integer itmplo0,itmplo1
      integer itmphi0,itmphi1
      REAL_T tmp(
     &           itmplo0:itmphi0,
     &           itmplo1:itmphi1,
     &           0:ntmpcomp-1)
      REAL_T alpha
      REAL_T beta
      integer redBlack
      REAL_T dx2t, lambda, lap, dxinv, helm
      integer i,j
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = twelve*dx*dx
      dxinv = one/(dx*dx)
      lambda = one/(alpha - beta*thirty*CH_SPACEDIM*twelfth*dxinv)
      lambda = lambda*(0.60)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
          lap = ( 
     &   sixteen*phi(i-1,j  ,n) - phi(i-2,j  ,n)
     & + sixteen*phi(i+1,j  ,n) - phi(i+2,j  ,n)
     & + sixteen*phi(i  ,j-1,n) - phi(i  ,j-2,n)
     & + sixteen*phi(i  ,j+1,n) - phi(i  ,j+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,j,n) )
     &       * twelfth * dxinv
          helm = alpha*phi(i,j,n) + beta*lap
          tmp(i,j,n) = phi(i,j,n) +
     &      lambda*( rhs(i,j,n) - helm )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else if (redBlack .eq. black) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               lap = ( 
     &           sixteen*tmp(i+1,j,n) - tmp(i+2,j,n)
     &         + sixteen*tmp(i-1,j,n) - tmp(i-2,j,n)
     &         + sixteen*tmp(i,j+1,n) - tmp(i,j+2,n)
     &         + sixteen*tmp(i,j-1,n) - tmp(i,j-2,n)
     &                     -(thirty*CH_SPACEDIM)*tmp(i,j,n) )
     &       * twelfth * dxinv
               helm = alpha*tmp(i,j,n) + beta*lap
               phi(i,j,n) = tmp(i,j,n) +
     &              lambda*( rhs(i,j,n) - helm )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,j,n) = tmp(i,j,n)
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else
           call MAYDAYERROR()
        end if
      enddo
      return
      end
      subroutine NEWGETFLUX4(
     &           flux
     &           ,ifluxlo0,ifluxlo1
     &           ,ifluxhi0,ifluxhi1
     &           ,nfluxcomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,beta_dx
     &           ,a_idir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfluxcomp
      integer ifluxlo0,ifluxlo1
      integer ifluxhi0,ifluxhi1
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           ifluxlo1:ifluxhi1,
     &           0:nfluxcomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T beta_dx
      integer a_idir
      INTEGER ncomp,n
      integer ii, jj
      integer i , j 
      ncomp = nphicomp
      
      ii = CHF_ID(a_idir, 0)
      jj = CHF_ID(a_idir, 1)
      do n = 0, ncomp-1
          
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          flux(i,j,n) = beta_dx * twelfth *
     &        ( fifteen*phi(i,j,n)
     &           + phi(i-2*ii,j-2*jj,n)
     &           - phi(i+ii,j+jj,n)
     &           - fifteen*phi(i-ii,j-jj,n) )
          
      enddo
      enddo
      enddo
      return
      end
      subroutine PROLONGLINEAR(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,ncoarsecomp
     &           ,ifineBoxlo0,ifineBoxlo1
     &           ,ifineBoxhi0,ifineBoxhi1
     &           ,icrseBoxlo0,icrseBoxlo1
     &           ,icrseBoxhi0,icrseBoxhi1
     &           ,r
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           0:ncoarsecomp-1)
      integer ifineBoxlo0,ifineBoxlo1
      integer ifineBoxhi0,ifineBoxhi1
      integer icrseBoxlo0,icrseBoxlo1
      integer icrseBoxhi0,icrseBoxhi1
      integer r
      INTEGER ncomp, n
      integer i ,j 
      integer ic,jc
      ncomp = nphicomp
      do n = 0, ncomp-1
          
      do j = ifineBoxlo1,ifineBoxhi1
      do i = ifineBoxlo0,ifineBoxhi0

          
           ic = i/r
           jc = j/r
           phi(i,j,n) =  phi(i,j,n) +
     &          coarse(ic,jc,n)
          
           if (ic.ne.icrseBoxhi0 .and.
     &         (ic*r.lt.i .or. ic.eq.icrseBoxlo0)) then
              phi(i,j,n) =  phi(i,j,n) +
     &             (coarse(ic+1,jc,n)
     &              - coarse(ic,jc,n))/r*(i+half-ic*r-half*r)
           else
              phi(i,j,n) =  phi(i,j,n) +
     &             (- coarse(ic-1,jc,n)
     &              + coarse(ic,jc,n))/r*(i+half-ic*r-half*r)
           endif
           if (jc.ne.icrseBoxhi1 .and.
     &         (jc*r.lt.j .or. jc.eq.icrseBoxlo1)) then
              phi(i,j,n) =  phi(i,j,n) +
     &             (coarse(ic,jc+1,n)
     &              - coarse(ic,jc,n))/r*(j+half-jc*r-half*r)
           else
              phi(i,j,n) =  phi(i,j,n) +
     &             (- coarse(ic,jc-1,n)
     &              + coarse(ic,jc,n))/r*(j+half-jc*r-half*r)
           endif
          
      enddo
      enddo
      enddo
      return
      end
