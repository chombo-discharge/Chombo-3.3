#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine NODEAVERAGE(
     &           coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,nfinecomp
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,ref_ratio
     &           ,weight
     &           ,iweightlo0,iweightlo1
     &           ,iweighthi0,iweighthi1
     &           ,nweightcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           0:nfinecomp-1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer ref_ratio
      integer nweightcomp
      integer iweightlo0,iweightlo1
      integer iweighthi0,iweighthi1
      REAL_T weight(
     &           iweightlo0:iweighthi0,
     &           iweightlo1:iweighthi1,
     &           0:nweightcomp-1)
      integer var
      integer icrse0,icrse1
      integer ifine0,ifine1
      integer ii0,ii1
      REAL_T csum
      do var = 0, ncoarsecomp - 1
         
      do icrse1 = iblo1,ibhi1
      do icrse0 = iblo0,ibhi0

            csum = 0
            
      do ii1 = iweightlo1,iweighthi1
      do ii0 = iweightlo0,iweighthi0

               
               ifine0 = icrse0*ref_ratio + ii0 
               ifine1 = icrse1*ref_ratio + ii1 
               csum = csum + weight( ii0,ii1, 0) *
     &                 fine( ifine0,ifine1, var )
            
      enddo
      enddo
            coarse( icrse0,icrse1, var ) = csum
         
      enddo
      enddo
      end do
      return
      end
      subroutine NODEAVERAGEPOINT(
     &           coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,nfinecomp
     &           ,pcrse
     &           ,ref_ratio
     &           ,weight
     &           ,iweightlo0,iweightlo1
     &           ,iweighthi0,iweighthi1
     &           ,nweightcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           0:nfinecomp-1)
      integer pcrse(0:1)
      integer ref_ratio
      integer nweightcomp
      integer iweightlo0,iweightlo1
      integer iweighthi0,iweighthi1
      REAL_T weight(
     &           iweightlo0:iweighthi0,
     &           iweightlo1:iweighthi1,
     &           0:nweightcomp-1)
      integer var
      integer ifine0,ifine1
      integer ii0,ii1
      REAL_T csum, weightpt, finept
      do var = 0, ncoarsecomp - 1
         csum = 0
         
      do ii1 = iweightlo1,iweighthi1
      do ii0 = iweightlo0,iweighthi0

            
            ifine0 = pcrse(0)*ref_ratio + ii0 
            ifine1 = pcrse(1)*ref_ratio + ii1 
            weightpt = weight( ii0,ii1, 0)
            finept = fine( ifine0,ifine1, var )
            csum = csum + weightpt*finept
         
      enddo
      enddo
         coarse(pcrse(0),pcrse(1), var ) = csum
      end do
      return
      end
      subroutine NODEAVERAGE_GETWEIGHTS(
     &           weight
     &           ,iweightlo0,iweightlo1
     &           ,iweighthi0,iweighthi1
     &           ,nweightcomp
     &           ,ref_ratio
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nweightcomp
      integer iweightlo0,iweightlo1
      integer iweighthi0,iweighthi1
      REAL_T weight(
     &           iweightlo0:iweighthi0,
     &           iweightlo1:iweighthi1,
     &           0:nweightcomp-1)
      integer ref_ratio
      integer ext, nxtrm
      integer ii0,ii1
      REAL_T ref_scale
      ext = ref_ratio / 2
      ref_scale = one / (ref_ratio**CH_SPACEDIM)
      
      do ii1 = iweightlo1,iweighthi1
      do ii0 = iweightlo0,iweighthi0

         nxtrm = 0
         
         if (iabs(ii0) .eq. ext) nxtrm = nxtrm + 1 
         if (iabs(ii1) .eq. ext) nxtrm = nxtrm + 1 
         weight( ii0,ii1, 0) = ref_scale * half**nxtrm
      
      enddo
      enddo
      return
      end
