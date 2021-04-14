#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine NODEOPLAP(
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
      REAL_T dxinv2, lphi
      integer var, ncomp
      integer  i, j
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = one / (dx*dx)
      do var = 0, ncomp-1
         
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

            
            lphi = ( (phi(i+1 , j  , var)
     &              - phi(i   , j  , var) )
     &            -  (phi(i   , j  , var)
     &              - phi(i-1 , j  , var) ) ) * dxinv2 
     &         +   ( (phi(i   , j+1, var)
     &              - phi(i   , j  , var) )
     &            -  (phi(i   , j  , var)
     &              - phi(i   , j-1, var) ) ) * dxinv2 
            lofphi(i, j, var) =  lphi
         
      enddo
      enddo
      end do
      return
      end
      subroutine NODEOPLAPPOINT(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1
     &           ,ilofphihi0,ilofphihi1
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,pt
     &           ,dx
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
      integer pt(0:1)
      REAL_T dx
      REAL_T dxinv2, lphi
      integer var, ncomp
      integer  i, j
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = one / (dx*dx)
      
      i = pt(0) 
      j = pt(1) 
      do var = 0, ncomp-1
         
         lphi = ( (phi(i+1 , j  , var)
     &           - phi(i   , j  , var) )
     &         -  (phi(i   , j  , var)
     &           - phi(i-1 , j  , var) ) ) * dxinv2 
     &      +   ( (phi(i   , j+1, var)
     &           - phi(i   , j  , var) )
     &         -  (phi(i   , j  , var)
     &           - phi(i   , j-1, var) ) ) * dxinv2 
         lofphi(i, j, var) =  lphi
      end do
      return
      end
      subroutine NODEGRAD(
     &           grdphi
     &           ,igrdphilo0,igrdphilo1
     &           ,igrdphihi0,igrdphihi1
     &           ,ngrdphicomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ngrdphicomp
      integer igrdphilo0,igrdphilo1
      integer igrdphihi0,igrdphihi1
      REAL_T grdphi(
     &           igrdphilo0:igrdphihi0,
     &           igrdphilo1:igrdphihi1,
     &           0:ngrdphicomp-1)
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
      REAL_T dxinvh
      integer var, ncomp, gbase
      integer  i, j
      ncomp = nphicomp
      if (CH_SPACEDIM * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = half / dx
      do var = 0, ncomp-1
         gbase = CH_SPACEDIM * var
         
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

            
            grdphi(i, j, gbase) =
     &           ( phi(i+1 , j   , var)
     &           - phi(i-1 , j   , var) ) * dxinvh 
            grdphi(i, j, gbase + 1) =
     &           ( phi(i   , j+1 , var)
     &           - phi(i   , j-1 , var) ) * dxinvh 
         
      enddo
      enddo
      end do
      return
      end
      subroutine NODEGRADPOINT(
     &           grdphi
     &           ,igrdphilo0,igrdphilo1
     &           ,igrdphihi0,igrdphihi1
     &           ,ngrdphicomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,pt
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ngrdphicomp
      integer igrdphilo0,igrdphilo1
      integer igrdphihi0,igrdphihi1
      REAL_T grdphi(
     &           igrdphilo0:igrdphihi0,
     &           igrdphilo1:igrdphihi1,
     &           0:ngrdphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer pt(0:1)
      REAL_T dx
      REAL_T dxinvh
      integer var, ncomp, gbase
      integer  i, j
      ncomp = nphicomp
      if (CH_SPACEDIM * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = half / dx
      
      i = pt(0) 
      j = pt(1) 
      do var = 0, ncomp-1
         gbase = CH_SPACEDIM * var
            
            grdphi(i, j, gbase) =
     &           ( phi(i+1 , j   , var)
     &           - phi(i-1 , j   , var) ) * dxinvh 
            grdphi(i, j, gbase + 1) =
     &           ( phi(i   , j+1 , var)
     &           - phi(i   , j-1 , var) ) * dxinvh 
      end do
      return
      end
      subroutine NODEGSRBLEVELLAP(
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
      integer redBlack
      REAL_T lambda
      REAL_T dxinv2, lphi
      integer i, j
      integer imin, imax, var, ncomp, indtot
      dxinv2 = one/(dx*dx)
      lambda = (dx*dx) / (two*CH_SPACEDIM)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAY_ERROR()
      endif
      do var = 0, ncomp - 1
#if CH_SPACEDIM>=3
         do k = iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM>=2
            do j = iregionlo1, iregionhi1
#endif
               imin = iregionlo0
               indtot = imin  + j 
               imin = imin + mod(indtot + redBlack, 2)
               imax = iregionhi0
               do i = imin, imax, 2
#ifndef NDEBUG
                  if (mod(i +j, 2) .ne. redBlack) then
                     print *, 'NODEGSRBLEVELLAP:  computing ',
     &                    i,  j, 
     &                    ' at pass ', redBlack
                  endif
#endif
            
            lphi = ( (phi(i+1 , j  , var)
     &              + phi(i-1 , j  , var) ) ) * dxinv2 
     &         +   ( (phi(i   , j+1, var)
     &              + phi(i   , j-1, var) ) ) * dxinv2 
                  phi(i,j, var) =
     &                 lambda * (lphi - rhs(i,j, var))
               enddo
#if CH_SPACEDIM>=2
            enddo
#endif
#if CH_SPACEDIM>=3
         enddo
#endif
      enddo
      return
      end
