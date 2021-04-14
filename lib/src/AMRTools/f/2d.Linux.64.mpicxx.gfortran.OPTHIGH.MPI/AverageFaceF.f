      subroutine AVERAGEFACE(
     & coarse
     & ,icoarselo0,icoarselo1
     & ,icoarsehi0,icoarsehi1
     & ,ncoarsecomp
     & ,fine
     & ,ifinelo0,ifinelo1
     & ,ifinehi0,ifinehi1
     & ,nfinecomp
     & ,icrseBoxlo0,icrseBoxlo1
     & ,icrseBoxhi0,icrseBoxhi1
     & ,dir
     & ,nRef
     & ,refFactor
     & ,irefBoxlo0,irefBoxlo1
     & ,irefBoxhi0,irefBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & 0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & 0:nfinecomp-1)
      integer icrseBoxlo0,icrseBoxlo1
      integer icrseBoxhi0,icrseBoxhi1
      integer dir
      integer nRef
      integer refFactor
      integer irefBoxlo0,irefBoxlo1
      integer irefBoxhi0,irefBoxhi1
      integer ic0,ic1
      integer ifine0,ifine1
      integer var
      integer ii0,ii1
      REAL*8 crseSum, ref_scale
      ref_scale = ((1.0d0)/refFactor)**(2 -1)
      do var=0, ncoarsecomp-1
      do ic1 = icrseBoxlo1,icrseBoxhi1
      do ic0 = icrseBoxlo0,icrseBoxhi0
         crseSum = 0
      do ii1 = irefBoxlo1,irefBoxhi1
      do ii0 = irefBoxlo0,irefBoxhi0
         ifine0=ic0*nRef+ii0
         ifine1=ic1*nRef+ii1
            crseSum = crseSum + fine(ifine0,ifine1,var)
      enddo
      enddo
            coarse(ic0,ic1,var) = ref_scale*crseSum
      enddo
      enddo
       enddo
       return
       end
      subroutine AVERAGEFACEHARMONIC(
     & coarse
     & ,icoarselo0,icoarselo1
     & ,icoarsehi0,icoarsehi1
     & ,ncoarsecomp
     & ,fine
     & ,ifinelo0,ifinelo1
     & ,ifinehi0,ifinehi1
     & ,nfinecomp
     & ,icrseBoxlo0,icrseBoxlo1
     & ,icrseBoxhi0,icrseBoxhi1
     & ,dir
     & ,nRef
     & ,refFactor
     & ,irefBoxlo0,irefBoxlo1
     & ,irefBoxhi0,irefBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & 0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & 0:nfinecomp-1)
      integer icrseBoxlo0,icrseBoxlo1
      integer icrseBoxhi0,icrseBoxhi1
      integer dir
      integer nRef
      integer refFactor
      integer irefBoxlo0,irefBoxlo1
      integer irefBoxhi0,irefBoxhi1
      integer ic0,ic1
      integer ifine0,ifine1
      integer var
      integer ii0,ii1
      REAL*8 crseSum, ref_scale
      ref_scale = ((1.0d0)/refFactor)**(2 -1)
      do var=0, ncoarsecomp-1
      do ic1 = icrseBoxlo1,icrseBoxhi1
      do ic0 = icrseBoxlo0,icrseBoxhi0
         crseSum = 0
      do ii1 = irefBoxlo1,irefBoxhi1
      do ii0 = irefBoxlo0,irefBoxhi0
         ifine0=ic0*nRef+ii0
         ifine1=ic1*nRef+ii1
            crseSum = crseSum + (1.0d0)/fine(ifine0,ifine1,var)
      enddo
      enddo
            coarse(ic0,ic1,var) = (1.0d0)/(ref_scale*crseSum)
      enddo
      enddo
       enddo
       return
       end
