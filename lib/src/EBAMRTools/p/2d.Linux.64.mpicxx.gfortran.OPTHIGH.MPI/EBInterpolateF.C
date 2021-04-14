#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine EBINTERPCONSTANT(
     &           fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,refratio
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer refratio
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer icc, jcc
      integer iff, jff
      integer ii, jj
      
      do jcc = iblo1,ibhi1
      do icc = iblo0,ibhi0

      
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0

      
      iff = icc*refratio + ii
      jff = jcc*refratio + jj
      fine(iff,jff) = coarse(icc,jcc)
      
      enddo
      enddo
      
      enddo
      enddo
      return
      end
      subroutine EBCENTRALSLOPE(
     &           slope
     &           ,islopelo0,islopelo1
     &           ,islopehi0,islopehi1
     &           ,state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer islopelo0,islopelo1
      integer islopehi0,islopehi1
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1)
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer dir
      integer i,ii, j,jj
      
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      
      do j = iblo1,ibhi1
      do i = iblo0,ibhi0

      slope(i,j) = half *(
     &     state(i+ii, j+jj) -
     &     state(i-ii, j-jj))
      
      enddo
      enddo
      return
      end
      subroutine EBHISIDESLOPE(
     &           slope
     &           ,islopelo0,islopelo1
     &           ,islopehi0,islopehi1
     &           ,state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer islopelo0,islopelo1
      integer islopehi0,islopehi1
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1)
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer dir
      integer i,ii, j,jj
      
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      
      do j = iblo1,ibhi1
      do i = iblo0,ibhi0

      slope(i,j) =
     &     state(i+ii,j+jj) -
     &     state(i   ,j   )
      
      enddo
      enddo
      return
      end
      subroutine EBLOSIDESLOPE(
     &           slope
     &           ,islopelo0,islopelo1
     &           ,islopehi0,islopehi1
     &           ,state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer islopelo0,islopelo1
      integer islopehi0,islopehi1
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1)
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer dir
      integer i,ii, j,jj
      
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      
      do j = iblo1,ibhi1
      do i = iblo0,ibhi0

      slope(i,j) =
     &     state(i   ,j   ) -
     &     state(i-ii,j-jj)
      
      enddo
      enddo
      return
      end
      subroutine EBMAXMINMOD(
     &           mmslope
     &           ,immslopelo0,immslopelo1
     &           ,immslopehi0,immslopehi1
     &           ,loslope
     &           ,iloslopelo0,iloslopelo1
     &           ,iloslopehi0,iloslopehi1
     &           ,hislope
     &           ,ihislopelo0,ihislopelo1
     &           ,ihislopehi0,ihislopehi1
     &           ,islopeboxlo0,islopeboxlo1
     &           ,islopeboxhi0,islopeboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer immslopelo0,immslopelo1
      integer immslopehi0,immslopehi1
      REAL_T mmslope(
     &           immslopelo0:immslopehi0,
     &           immslopelo1:immslopehi1)
      integer iloslopelo0,iloslopelo1
      integer iloslopehi0,iloslopehi1
      REAL_T loslope(
     &           iloslopelo0:iloslopehi0,
     &           iloslopelo1:iloslopehi1)
      integer ihislopelo0,ihislopelo1
      integer ihislopehi0,ihislopehi1
      REAL_T hislope(
     &           ihislopelo0:ihislopehi0,
     &           ihislopelo1:ihislopehi1)
      integer islopeboxlo0,islopeboxlo1
      integer islopeboxhi0,islopeboxhi1
      integer i,j
      REAL_T deltal, deltar, mono, rsign, finslope
      
      do j = islopeboxlo1,islopeboxhi1
      do i = islopeboxlo0,islopeboxhi0

      deltal = loslope(i,j)
      deltar = hislope(i,j)
      mono = deltal*deltar
      if(mono .gt. zero) then
         rsign = sign(one, deltal + deltar)
         finslope = rsign*(min(abs(deltal), abs(deltar)))
      else
         finslope = zero
      endif
      mmslope(i,j) = finslope
      
      enddo
      enddo
      return
      end
      subroutine EBINTERPLINEAR(
     &           fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,slope
     &           ,islopelo0,islopelo1
     &           ,islopehi0,islopehi1
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,dir
     &           ,refratio
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer islopelo0,islopelo1
      integer islopehi0,islopehi1
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer dir
      integer refratio
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer icc, jcc
      integer iff, jff
      integer ii, jj
      integer id
      REAL_T dxf
      
      do jcc = iblo1,ibhi1
      do icc = iblo0,ibhi0

      
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0

      
      iff = icc*refratio + ii
      jff = jcc*refratio + jj
      
      if(dir .eq. 0) then
         id = ii
      else if(dir .eq. 1) then
         id = jj
      endif
      dxf = -half +((id+half) / refratio)
      fine(iff,jff) =
     &     fine(iff,jff) +
     &     dxf*slope(icc, jcc)
      
      enddo
      enddo
      
      enddo
      enddo
      return
      end
      subroutine EBINTERPSMOOTHERLINEAR(
     &           fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,refratio
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer refratio
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      return
      end
      subroutine EBINTERPQUADRATIC(
     &           fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,refratio
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer refratio
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      return
      end
      subroutine EBINTERPQUADSHIFT(
     &           fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,refratio
     &           ,ishift
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer refratio
      integer ishift(0:1)
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      return
      end
