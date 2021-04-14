#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
        subroutine SETDIRICHLETFACEFLUX(
     &           faceFlux
     &           ,ifaceFluxlo0,ifaceFluxlo1
     &           ,ifaceFluxhi0,ifaceFluxhi1
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,value
     &           ,dx
     &           ,idir
     &           ,iside
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifaceFluxlo0,ifaceFluxlo1
      integer ifaceFluxhi0,ifaceFluxhi1
      REAL_T faceFlux(
     &           ifaceFluxlo0:ifaceFluxhi0,
     &           ifaceFluxlo1:ifaceFluxhi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      REAL_T value
      REAL_T dx(0:1)
      integer idir
      integer iside
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t ihdx
        ihdx = two/dx(idir)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          faceFlux(i,j) = iside * ihdx * (phi(i,j) - value)
        
      enddo
      enddo
        return
        end
        subroutine SETHODIRICHLETFACEFLUX(
     &           faceFlux
     &           ,ifaceFluxlo0,ifaceFluxlo1
     &           ,ifaceFluxhi0,ifaceFluxhi1
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,value
     &           ,dx
     &           ,idir
     &           ,iside
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ifaceFluxlo0,ifaceFluxlo1
      integer ifaceFluxhi0,ifaceFluxhi1
      REAL_T faceFlux(
     &           ifaceFluxlo0:ifaceFluxhi0,
     &           ifaceFluxlo1:ifaceFluxhi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      REAL_T value
      REAL_T dx(0:1)
      integer idir
      integer iside
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j, ioff,joff
        real_t dxinv
        ioff = iside*chf_id(0,idir)
                  joff = iside*chf_id(1,idir)
        dxinv = one/dx(idir)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          faceFlux(i,j) =
     &       iside * dxinv * (three*phi(i,j)
     &                  - one/three*phi(i+ioff,j+joff)
     &                  - eight/three*value)
        
      enddo
      enddo
        return
        end
