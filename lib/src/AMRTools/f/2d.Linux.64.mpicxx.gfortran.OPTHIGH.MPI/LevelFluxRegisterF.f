      subroutine INCREMENTFINE(
     & fine
     & ,ifinelo0,ifinelo1
     & ,ifinehi0,ifinehi1
     & ,nfinecomp
     & ,cFine
     & ,icFinelo0,icFinelo1
     & ,icFinehi0,icFinehi1
     & ,ncFinecomp
     & ,ifineBoxlo0,ifineBoxlo1
     & ,ifineBoxhi0,ifineBoxhi1
     & ,nRef
     & ,scale
     & ,srcStart
     & ,destStart
     & ,ncomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & 0:nfinecomp-1)
      integer ncFinecomp
      integer icFinelo0,icFinelo1
      integer icFinehi0,icFinehi1
      REAL*8 cFine(
     & icFinelo0:icFinehi0,
     & icFinelo1:icFinehi1,
     & 0:ncFinecomp-1)
      integer ifineBoxlo0,ifineBoxlo1
      integer ifineBoxhi0,ifineBoxhi1
      integer nRef(0:1)
      REAL*8 scale
      integer srcStart
      integer destStart
      integer ncomp
      integer i0,i1
      integer ii0,ii1
      integer var, srcComp, destComp
      do var=0, ncomp-1
         srcComp = srcStart + var
         destComp = destStart + var
      do i1 = ifineBoxlo1,ifineBoxhi1
      do i0 = ifineBoxlo0,ifineBoxhi0
            ii0=i0/nRef(0)
            ii1=i1/nRef(1)
            cFine(ii0,ii1,destComp) =
     & cFine(ii0,ii1, destComp) +
     & scale * fine(i0,i1, srcComp)
      enddo
      enddo
      enddo
      return
      end
