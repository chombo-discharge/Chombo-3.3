      subroutine AVEFACESCALTOFACEVECT(
     & facevect
     & ,ifacevectlo0,ifacevectlo1
     & ,ifacevecthi0,ifacevecthi1
     & ,nfacevectcomp
     & ,facescal
     & ,ifacescallo0,ifacescallo1
     & ,ifacescalhi0,ifacescalhi1
     & ,facedir
     & ,vectdir
     & ,idcalcfacelo0,idcalcfacelo1
     & ,idcalcfacehi0,idcalcfacehi1
     & ,ioffboxlo0,ioffboxlo1
     & ,ioffboxhi0,ioffboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfacevectcomp
      integer ifacevectlo0,ifacevectlo1
      integer ifacevecthi0,ifacevecthi1
      REAL*8 facevect(
     & ifacevectlo0:ifacevecthi0,
     & ifacevectlo1:ifacevecthi1,
     & 0:nfacevectcomp-1)
      integer ifacescallo0,ifacescallo1
      integer ifacescalhi0,ifacescalhi1
      REAL*8 facescal(
     & ifacescallo0:ifacescalhi0,
     & ifacescallo1:ifacescalhi1)
      integer facedir
      integer vectdir
      integer idcalcfacelo0,idcalcfacelo1
      integer idcalcfacehi0,idcalcfacehi1
      integer ioffboxlo0,ioffboxlo1
      integer ioffboxhi0,ioffboxhi1
      integer i0,i1
      integer ioff0,ioff1
      integer numpts
      REAL*8 wt, tot
      if (facedir .eq. vectdir) then
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0
            facevect(i0,i1, vectdir) = facescal(i0,i1)
      enddo
      enddo
      else
         numpts = 2**2
         wt = (1.0d0) / (numpts * (1.0d0))
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0
            tot = 0
      do ioff1 = ioffboxlo1,ioffboxhi1
      do ioff0 = ioffboxlo0,ioffboxhi0
               tot = tot + facescal(i0 +ioff0,i1 +ioff1)
      enddo
      enddo
            facevect(i0,i1, vectdir) = wt * tot
      enddo
      enddo
      endif
      return
      end
      subroutine AVESCALTOFACE(
     & facescal
     & ,ifacescallo0,ifacescallo1
     & ,ifacescalhi0,ifacescalhi1
     & ,cellscal
     & ,icellscallo0,icellscallo1
     & ,icellscalhi0,icellscalhi1
     & ,idir
     & ,idcalcfacelo0,idcalcfacelo1
     & ,idcalcfacehi0,idcalcfacehi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacescallo0,ifacescallo1
      integer ifacescalhi0,ifacescalhi1
      REAL*8 facescal(
     & ifacescallo0:ifacescalhi0,
     & ifacescallo1:ifacescalhi1)
      integer icellscallo0,icellscallo1
      integer icellscalhi0,icellscalhi1
      REAL*8 cellscal(
     & icellscallo0:icellscalhi0,
     & icellscallo1:icellscalhi1)
      integer idir
      integer idcalcfacelo0,idcalcfacelo1
      integer idcalcfacehi0,idcalcfacehi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0
         facescal(i0,i1) = (0.500d0) * (
     & cellscal(i0 +f2cLo0,i1 +f2cLo1) +
     & cellscal(i0 +f2cHi0,i1 +f2cHi1) )
      enddo
      enddo
      return
      end
      subroutine AVECELLTOFACE(
     & facevel
     & ,ifacevello0,ifacevello1
     & ,ifacevelhi0,ifacevelhi1
     & ,cellvel
     & ,icellvello0,icellvello1
     & ,icellvelhi0,icellvelhi1
     & ,ncellvelcomp
     & ,idir
     & ,idcalcfacelo0,idcalcfacelo1
     & ,idcalcfacehi0,idcalcfacehi1
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
      integer ncellvelcomp
      integer icellvello0,icellvello1
      integer icellvelhi0,icellvelhi1
      REAL*8 cellvel(
     & icellvello0:icellvelhi0,
     & icellvello1:icellvelhi1,
     & 0:ncellvelcomp-1)
      integer idir
      integer idcalcfacelo0,idcalcfacelo1
      integer idcalcfacehi0,idcalcfacehi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      f2cLo0= -1*CHF_ID(0, idir)
      f2cLo1= -1*CHF_ID(1, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      f2cHi1= 0*CHF_ID(1, idir)
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0
         facevel(i0,i1) = (0.500d0) * (
     & cellvel(i0 +f2cLo0,i1 +f2cLo1, idir) +
     & cellvel(i0 +f2cHi0,i1 +f2cHi1, idir) )
      enddo
      enddo
      return
      end
      subroutine AVEFACETOCELL(
     & cellvel
     & ,icellvello0,icellvello1
     & ,icellvelhi0,icellvelhi1
     & ,ncellvelcomp
     & ,facevel
     & ,ifacevello0,ifacevello1
     & ,ifacevelhi0,ifacevelhi1
     & ,idir
     & ,idcalccelllo0,idcalccelllo1
     & ,idcalccellhi0,idcalccellhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ncellvelcomp
      integer icellvello0,icellvello1
      integer icellvelhi0,icellvelhi1
      REAL*8 cellvel(
     & icellvello0:icellvelhi0,
     & icellvello1:icellvelhi1,
     & 0:ncellvelcomp-1)
      integer ifacevello0,ifacevello1
      integer ifacevelhi0,ifacevelhi1
      REAL*8 facevel(
     & ifacevello0:ifacevelhi0,
     & ifacevello1:ifacevelhi1)
      integer idir
      integer idcalccelllo0,idcalccelllo1
      integer idcalccellhi0,idcalccellhi1
      integer i0,i1
      integer c2fLo0,c2fLo1
      integer c2fHi0,c2fHi1
      c2fLo0= 0*CHF_ID(0, idir)
      c2fLo1= 0*CHF_ID(1, idir)
      c2fHi0= 1*CHF_ID(0, idir)
      c2fHi1= 1*CHF_ID(1, idir)
      do i1 = idcalccelllo1,idcalccellhi1
      do i0 = idcalccelllo0,idcalccellhi0
         cellvel(i0,i1, idir) = (0.500d0) * (
     & facevel(i0 +c2fLo0,i1 +c2fLo1) +
     & facevel(i0 +c2fHi0,i1 +c2fHi1) )
      enddo
      enddo
      return
      end
      subroutine MAGNITUDEF(
     & magdata
     & ,imagdatalo0,imagdatalo1
     & ,imagdatahi0,imagdatahi1
     & ,data
     & ,idatalo0,idatalo1
     & ,idatahi0,idatahi1
     & ,ndatacomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imagdatalo0,imagdatalo1
      integer imagdatahi0,imagdatahi1
      REAL*8 magdata(
     & imagdatalo0:imagdatahi0,
     & imagdatalo1:imagdatahi1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1,
     & 0:ndatacomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer iv
      REAL*8 cur,sum
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         sum = 0
         do iv = 0, ndatacomp-1
            cur = data(i0,i1, iv)
            sum = sum + cur*cur
         enddo
         magdata(i0,i1) = sqrt(sum)
      enddo
      enddo
      return
      end
        subroutine GETRELGRADF(
     & du
     & ,idulo0,idulo1
     & ,iduhi0,iduhi1
     & ,u
     & ,iulo0,iulo1
     & ,iuhi0,iuhi1
     & ,idir
     & ,iloBoxlo0,iloBoxlo1
     & ,iloBoxhi0,iloBoxhi1
     & ,hasLo
     & ,ihiBoxlo0,ihiBoxlo1
     & ,ihiBoxhi0,ihiBoxhi1
     & ,hasHi
     & ,icenterBoxlo0,icenterBoxlo1
     & ,icenterBoxhi0,icenterBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idulo0,idulo1
      integer iduhi0,iduhi1
      REAL*8 du(
     & idulo0:iduhi0,
     & idulo1:iduhi1)
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL*8 u(
     & iulo0:iuhi0,
     & iulo1:iuhi1)
      integer idir
      integer iloBoxlo0,iloBoxlo1
      integer iloBoxhi0,iloBoxhi1
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1
      integer ihiBoxhi0,ihiBoxhi1
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1
      integer icenterBoxhi0,icenterBoxhi1
        integer i0,i1
        integer ioff0,ioff1
        REAL*8 diff,aver
      ioff0=CHF_ID(0, idir)
      ioff1=CHF_ID(1, idir)
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0
          diff = (0.500d0)*(u(i0 +ioff0,i1 +ioff1)
     & -u(i0 -ioff0,i1 -ioff1))
          aver = (0.500d0)*(u(i0 +ioff0,i1 +ioff1)
     & +u(i0 -ioff0,i1 -ioff1))
          du(i0,i1) = diff / aver
      enddo
      enddo
        if (hasLo .eq. 1) then
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0
            diff = u(i0 +ioff0,i1 +ioff1) - u(i0,i1)
            aver = (0.500d0)*(u(i0 +ioff0,i1 +ioff1) + u(i0,i1))
            du(i0,i1) = diff / aver
      enddo
      enddo
        endif
        if (hasHi .eq. 1) then
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0
            diff = u(i0,i1) - u(i0 -ioff0,i1 -ioff1)
            aver = (0.500d0)*(u(i0,i1) + u(i0 -ioff0,i1 -ioff1))
            du(i0,i1) = diff / aver
      enddo
      enddo
        endif
        return
        end
      subroutine POSTNORMALSOURCE(
     & dWminus
     & ,idWminuslo0,idWminuslo1
     & ,idWminushi0,idWminushi1
     & ,ndWminuscomp
     & ,dWplus
     & ,idWpluslo0,idWpluslo1
     & ,idWplushi0,idWplushi1
     & ,ndWpluscomp
     & ,W
     & ,iWlo0,iWlo1
     & ,iWhi0,iWhi1
     & ,nWcomp
     & ,advVel
     & ,iadvVello0,iadvVello1
     & ,iadvVelhi0,iadvVelhi1
     & ,dt
     & ,dx
     & ,idir
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndWminuscomp
      integer idWminuslo0,idWminuslo1
      integer idWminushi0,idWminushi1
      REAL*8 dWminus(
     & idWminuslo0:idWminushi0,
     & idWminuslo1:idWminushi1,
     & 0:ndWminuscomp-1)
      integer ndWpluscomp
      integer idWpluslo0,idWpluslo1
      integer idWplushi0,idWplushi1
      REAL*8 dWplus(
     & idWpluslo0:idWplushi0,
     & idWpluslo1:idWplushi1,
     & 0:ndWpluscomp-1)
      integer nWcomp
      integer iWlo0,iWlo1
      integer iWhi0,iWhi1
      REAL*8 W(
     & iWlo0:iWhi0,
     & iWlo1:iWhi1,
     & 0:nWcomp-1)
      integer iadvVello0,iadvVello1
      integer iadvVelhi0,iadvVelhi1
      REAL*8 advVel(
     & iadvVello0:iadvVelhi0,
     & iadvVello1:iadvVelhi1)
      REAL*8 dt
      REAL*8 dx
      integer idir
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer c2fLo0,c2fLo1
      integer c2fHi0,c2fHi1
      integer n
      REAL*8 dudx
      c2fLo0= 0*CHF_ID(0, idir)
      c2fLo1= 0*CHF_ID(1, idir)
      c2fHi0= 1*CHF_ID(0, idir)
      c2fHi1= 1*CHF_ID(1, idir)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         dudx = advVel(i0 +c2fHi0,i1 +c2fHi1)
     & - advVel(i0 +c2fLo0,i1 +c2fLo1)
         dudx = dt*dudx / dx * (0.500d0)
         do n=0, nWcomp-1
            dWplus(i0,i1, n) =
     & dWPlus(i0,i1, n) - dudx * W(i0,i1,n)
            dWMinus(i0,i1, n) =
     & dWMinus(i0,i1, n) - dudx * W(i0,i1,n)
         end do
      enddo
      enddo
      return
      end
      subroutine RIEMANNF(
     & Wgdnv
     & ,iWgdnvlo0,iWgdnvlo1
     & ,iWgdnvhi0,iWgdnvhi1
     & ,nWgdnvcomp
     & ,WLeft
     & ,iWLeftlo0,iWLeftlo1
     & ,iWLefthi0,iWLefthi1
     & ,nWLeftcomp
     & ,WRight
     & ,iWRightlo0,iWRightlo1
     & ,iWRighthi0,iWRighthi1
     & ,nWRightcomp
     & ,advVel
     & ,iadvVello0,iadvVello1
     & ,iadvVelhi0,iadvVelhi1
     & ,idir
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nWgdnvcomp
      integer iWgdnvlo0,iWgdnvlo1
      integer iWgdnvhi0,iWgdnvhi1
      REAL*8 Wgdnv(
     & iWgdnvlo0:iWgdnvhi0,
     & iWgdnvlo1:iWgdnvhi1,
     & 0:nWgdnvcomp-1)
      integer nWLeftcomp
      integer iWLeftlo0,iWLeftlo1
      integer iWLefthi0,iWLefthi1
      REAL*8 WLeft(
     & iWLeftlo0:iWLefthi0,
     & iWLeftlo1:iWLefthi1,
     & 0:nWLeftcomp-1)
      integer nWRightcomp
      integer iWRightlo0,iWRightlo1
      integer iWRighthi0,iWRighthi1
      REAL*8 WRight(
     & iWRightlo0:iWRighthi0,
     & iWRightlo1:iWRighthi1,
     & 0:nWRightcomp-1)
      integer iadvVello0,iadvVello1
      integer iadvVelhi0,iadvVelhi1
      REAL*8 advVel(
     & iadvVello0:iadvVelhi0,
     & iadvVello1:iadvVelhi1)
      integer idir
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer n
      REAL*8 sl,sr
      REAL*8 so
      REAL*8 ustar
      do n=0, nWLeftcomp-1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            sl = WLeft(i0,i1, n)
            sr = WRight(i0,i1, n)
            ustar = advVel(i0,i1)
            if (ustar .gt. 0) then
               so = sl
            else
               so = sr
            endif
            if (abs(ustar) .lt. 1.0d-9) then
               so = (0.500d0)*(sl+sr)
            endif
            Wgdnv(i0,i1, n) = so
      enddo
      enddo
      end do
      return
      end
      subroutine QUASILINEARUPDATE(
     & AdWdx
     & ,iAdWdxlo0,iAdWdxlo1
     & ,iAdWdxhi0,iAdWdxhi1
     & ,nAdWdxcomp
     & ,WHalf
     & ,iWHalflo0,iWHalflo1
     & ,iWHalfhi0,iWHalfhi1
     & ,nWHalfcomp
     & ,cellVel
     & ,icellVello0,icellVello1
     & ,icellVelhi0,icellVelhi1
     & ,scale
     & ,idir
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nAdWdxcomp
      integer iAdWdxlo0,iAdWdxlo1
      integer iAdWdxhi0,iAdWdxhi1
      REAL*8 AdWdx(
     & iAdWdxlo0:iAdWdxhi0,
     & iAdWdxlo1:iAdWdxhi1,
     & 0:nAdWdxcomp-1)
      integer nWHalfcomp
      integer iWHalflo0,iWHalflo1
      integer iWHalfhi0,iWHalfhi1
      REAL*8 WHalf(
     & iWHalflo0:iWHalfhi0,
     & iWHalflo1:iWHalfhi1,
     & 0:nWHalfcomp-1)
      integer icellVello0,icellVello1
      integer icellVelhi0,icellVelhi1
      REAL*8 cellVel(
     & icellVello0:icellVelhi0,
     & icellVello1:icellVelhi1)
      REAL*8 scale
      integer idir
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer c2fLo0,c2fLo1
      integer c2fHi0,c2fHi1
      integer n
      c2fLo0= 0*CHF_ID(0, idir)
      c2fLo1= 0*CHF_ID(1, idir)
      c2fHi0= 1*CHF_ID(0, idir)
      c2fHi1= 1*CHF_ID(1, idir)
      do n=0, nAdWdxcomp-1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            AdWdx(i0,i1, n) =
     & scale * cellVel(i0,i1) *
     & (WHalf(i0 +c2fHi0,i1 +c2fHi1, n) -
     & WHalf(i0 +c2fLo0,i1 +c2fLo1, n))
      enddo
      enddo
      end do
      return
      end
      subroutine AVECELLVECTOFACEVEC(
     & facevec
     & ,ifaceveclo0,ifaceveclo1
     & ,ifacevechi0,ifacevechi1
     & ,nfaceveccomp
     & ,cellvec
     & ,icellveclo0,icellveclo1
     & ,icellvechi0,icellvechi1
     & ,ncellveccomp
     & ,facedir
     & ,idcalcfacelo0,idcalcfacelo1
     & ,idcalcfacehi0,idcalcfacehi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfaceveccomp
      integer ifaceveclo0,ifaceveclo1
      integer ifacevechi0,ifacevechi1
      REAL*8 facevec(
     & ifaceveclo0:ifacevechi0,
     & ifaceveclo1:ifacevechi1,
     & 0:nfaceveccomp-1)
      integer ncellveccomp
      integer icellveclo0,icellveclo1
      integer icellvechi0,icellvechi1
      REAL*8 cellvec(
     & icellveclo0:icellvechi0,
     & icellveclo1:icellvechi1,
     & 0:ncellveccomp-1)
      integer facedir
      integer idcalcfacelo0,idcalcfacelo1
      integer idcalcfacehi0,idcalcfacehi1
      integer i0,i1
      integer f2cLo0,f2cLo1
      integer f2cHi0,f2cHi1
      integer idir
      f2cLo0= -1*CHF_ID(0, facedir)
      f2cLo1= -1*CHF_ID(1, facedir)
      f2cHi0= 0*CHF_ID(0, facedir)
      f2cHi1= 0*CHF_ID(1, facedir)
      do idir = 0, 2 -1
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0
            facevec(i0,i1, idir) = (0.500d0) * (
     & cellvec(i0 +f2cLo0,i1 +f2cLo1, idir) +
     & cellvec(i0 +f2cHi0,i1 +f2cHi1, idir) )
      enddo
      enddo
      enddo
      return
      end
