      subroutine CONDUCTIVITYGSRB(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,rhs
     & ,irhslo0,irhslo1
     & ,irhshi0,irhshi1
     & ,relcoef
     & ,irelcoeflo0,irelcoeflo1
     & ,irelcoefhi0,irelcoefhi1
     & ,acoef
     & ,iacoeflo0,iacoeflo1
     & ,iacoefhi0,iacoefhi1
     & ,b0
     & ,ib0lo0,ib0lo1
     & ,ib0hi0,ib0hi1
     & ,b1
     & ,ib1lo0,ib1lo1
     & ,ib1hi0,ib1hi1
     & ,b2
     & ,ib2lo0,ib2lo1
     & ,ib2hi0,ib2hi1
     & ,alpha
     & ,beta
     & ,dx
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,redBlack
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1)
      integer irelcoeflo0,irelcoeflo1
      integer irelcoefhi0,irelcoefhi1
      REAL*8 relcoef(
     & irelcoeflo0:irelcoefhi0,
     & irelcoeflo1:irelcoefhi1)
      integer iacoeflo0,iacoeflo1
      integer iacoefhi0,iacoefhi1
      REAL*8 acoef(
     & iacoeflo0:iacoefhi0,
     & iacoeflo1:iacoefhi1)
      integer ib0lo0,ib0lo1
      integer ib0hi0,ib0hi1
      REAL*8 b0(
     & ib0lo0:ib0hi0,
     & ib0lo1:ib0hi1)
      integer ib1lo0,ib1lo1
      integer ib1hi0,ib1hi1
      REAL*8 b1(
     & ib1lo0:ib1hi0,
     & ib1lo1:ib1hi1)
      integer ib2lo0,ib2lo1
      integer ib2hi0,ib2hi1
      REAL*8 b2(
     & ib2lo0:ib2hi0,
     & ib2lo1:ib2hi1)
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dx
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer redBlack
      integer i,j
      REAL*8 laplphi, dx0
      integer imin,imax,indtot
      dx0 = beta/(dx * dx)
         do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               laplphi =
     & (b0(i+1,j )*(phi(i+1,j ) - phi(i ,j ))
     & - b0(i ,j )*(phi(i ,j ) - phi(i-1,j )))*dx0
     & +(b1(i ,j+1)*(phi(i ,j+1) - phi(i ,j ))
     & - b1(i ,j )*(phi(i ,j ) - phi(i ,j-1)))*dx0
               laplphi = laplphi + alpha * acoef(i,j) * phi(i,j)
               phi(i,j) = phi(i,j) + relcoef(i,j)*(rhs(i,j) - laplphi)
            enddo
         enddo
      return
      end
        subroutine EBCOREGAPPLYDOMAINFLUX(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,faceflux
     & ,ifacefluxlo0,ifacefluxlo1
     & ,ifacefluxhi0,ifacefluxhi1
     & ,bc
     & ,ibclo0,ibclo1
     & ,ibchi0,ibchi1
     & ,dx
     & ,side
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
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer ifacefluxlo0,ifacefluxlo1
      integer ifacefluxhi0,ifacefluxhi1
      REAL*8 faceflux(
     & ifacefluxlo0:ifacefluxhi0,
     & ifacefluxlo1:ifacefluxhi1)
      integer ibclo0,ibclo1
      integer ibchi0,ibchi1
      REAL*8 bc(
     & ibclo0:ibchi0,
     & ibclo1:ibchi1)
      REAL*8 dx
      integer side
      integer idir
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j, ioff,joff
        REAL*8 scaledflux, tol
        tol = 1.0d-15
        ioff = chf_id(0,idir)
        joff = chf_id(1,idir)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        if (side.eq.1) then
           if (ABS(bc(i-ioff,j-joff)) .GT. tol) then
              scaledflux = faceflux(i-ioff,j-joff)/bc(i-ioff,j-joff)
              phi(i,j) = phi(i-ioff,j-joff) + scaledflux*dx
           else
              phi(i,j) = phi(i-ioff,j-joff)
           endif
        else
           if(ABS(bc(i+ioff,j+joff)) .GT. tol) then
              scaledflux = faceflux(i+ioff,j+joff)/bc(i+ioff,j+joff)
              phi(i,j) = phi(i+ioff,j+joff) - scaledflux*dx
           else
              phi(i,j) = phi(i+ioff,j+joff)
           endif
        endif
      enddo
      enddo
        return
        end
      subroutine APPLYOPEBCONDNOBCS(
     & opphidir
     & ,iopphidirlo0,iopphidirlo1
     & ,iopphidirhi0,iopphidirhi1
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,aco
     & ,iacolo0,iacolo1
     & ,iacohi0,iacohi1
     & ,b0
     & ,ib0lo0,ib0lo1
     & ,ib0hi0,ib0hi1
     & ,b1
     & ,ib1lo0,ib1lo1
     & ,ib1hi0,ib1hi1
     & ,b2
     & ,ib2lo0,ib2lo1
     & ,ib2hi0,ib2hi1
     & ,dx
     & ,alpha
     & ,beta
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iopphidirlo0,iopphidirlo1
      integer iopphidirhi0,iopphidirhi1
      REAL*8 opphidir(
     & iopphidirlo0:iopphidirhi0,
     & iopphidirlo1:iopphidirhi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer iacolo0,iacolo1
      integer iacohi0,iacohi1
      REAL*8 aco(
     & iacolo0:iacohi0,
     & iacolo1:iacohi1)
      integer ib0lo0,ib0lo1
      integer ib0hi0,ib0hi1
      REAL*8 b0(
     & ib0lo0:ib0hi0,
     & ib0lo1:ib0hi1)
      integer ib1lo0,ib1lo1
      integer ib1hi0,ib1hi1
      REAL*8 b1(
     & ib1lo0:ib1hi0,
     & ib1lo1:ib1hi1)
      integer ib2lo0,ib2lo1
      integer ib2hi0,ib2hi1
      REAL*8 b2(
     & ib2lo0:ib2hi0,
     & ib2lo1:ib2hi1)
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      REAL*8 laplphi, dxinvsq
      dxinvsq = (1.0d0)/(dx * dx)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      laplphi =
     & (b0(i+1,j )*(phi(i+1,j ) - phi(i ,j ))
     & - b0(i ,j )*(phi(i ,j ) - phi(i-1,j )))*dxinvsq
     & +(b1(i ,j+1)*(phi(i ,j+1) - phi(i ,j ))
     & - b1(i ,j )*(phi(i ,j ) - phi(i ,j-1)))*dxinvsq
      opphidir(i,j) = alpha*aco(i,j)*phi(i,j) + beta*laplphi
      enddo
      enddo
      ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(4 + 6
     &*2)
      return
      end
      subroutine GSCOLOREBCONDOP(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,rhs
     & ,irhslo0,irhslo1
     & ,irhshi0,irhshi1
     & ,relco
     & ,irelcolo0,irelcolo1
     & ,irelcohi0,irelcohi1
     & ,aco
     & ,iacolo0,iacolo1
     & ,iacohi0,iacohi1
     & ,b0
     & ,ib0lo0,ib0lo1
     & ,ib0hi0,ib0hi1
     & ,b1
     & ,ib1lo0,ib1lo1
     & ,ib1hi0,ib1hi1
     & ,b2
     & ,ib2lo0,ib2lo1
     & ,ib2hi0,ib2hi1
     & ,dx
     & ,alpha
     & ,beta
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1)
      integer irelcolo0,irelcolo1
      integer irelcohi0,irelcohi1
      REAL*8 relco(
     & irelcolo0:irelcohi0,
     & irelcolo1:irelcohi1)
      integer iacolo0,iacolo1
      integer iacohi0,iacohi1
      REAL*8 aco(
     & iacolo0:iacohi0,
     & iacolo1:iacohi1)
      integer ib0lo0,ib0lo1
      integer ib0hi0,ib0hi1
      REAL*8 b0(
     & ib0lo0:ib0hi0,
     & ib0lo1:ib0hi1)
      integer ib1lo0,ib1lo1
      integer ib1hi0,ib1hi1
      REAL*8 b1(
     & ib1lo0:ib1hi0,
     & ib1lo1:ib1hi1)
      integer ib2lo0,ib2lo1
      integer ib2hi0,ib2hi1
      REAL*8 b2(
     & ib2lo0:ib2hi0,
     & ib2lo1:ib2hi1)
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j,ncolors
      REAL*8 laplphi, dxinvsq, opphipt, rhspt, relcopt, phipt
      dxinvsq = (1.0d0)/(dx * dx)
      do j = iboxlo1,iboxhi1,2
      do i = iboxlo0,iboxhi0,2
      laplphi =
     & (b0(i+1,j )*(phi(i+1,j ) - phi(i ,j ))
     & - b0(i ,j )*(phi(i ,j ) - phi(i-1,j )))*dxinvsq
     & +(b1(i ,j+1)*(phi(i ,j+1) - phi(i ,j ))
     & - b1(i ,j )*(phi(i ,j ) - phi(i ,j-1)))*dxinvsq
      opphipt = alpha*aco(i,j)*phi(i,j) + beta*laplphi
      phipt = phi(i,j)
      rhspt = rhs(i,j)
      relcopt = relco(i,j)
      phi(i,j) = phipt + relcopt*(rhspt - opphipt)
      enddo
      enddo
      ncolors = 2 *2
      ch_flops=ch_flops+((iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(4 + 
     &6*2))/ncolors
      return
      end
      subroutine CONDUCTIVITYINPLACE(
     & opphidir
     & ,iopphidirlo0,iopphidirlo1
     & ,iopphidirhi0,iopphidirhi1
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,b0
     & ,ib0lo0,ib0lo1
     & ,ib0hi0,ib0hi1
     & ,b1
     & ,ib1lo0,ib1lo1
     & ,ib1hi0,ib1hi1
     & ,b2
     & ,ib2lo0,ib2lo1
     & ,ib2hi0,ib2hi1
     & ,beta
     & ,dx
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iopphidirlo0,iopphidirlo1
      integer iopphidirhi0,iopphidirhi1
      REAL*8 opphidir(
     & iopphidirlo0:iopphidirhi0,
     & iopphidirlo1:iopphidirhi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer ib0lo0,ib0lo1
      integer ib0hi0,ib0hi1
      REAL*8 b0(
     & ib0lo0:ib0hi0,
     & ib0lo1:ib0hi1)
      integer ib1lo0,ib1lo1
      integer ib1hi0,ib1hi1
      REAL*8 b1(
     & ib1lo0:ib1hi0,
     & ib1lo1:ib1hi1)
      integer ib2lo0,ib2lo1
      integer ib2hi0,ib2hi1
      REAL*8 b2(
     & ib2lo0:ib2hi0,
     & ib2lo1:ib2hi1)
      REAL*8 beta
      REAL*8 dx
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      REAL*8 laplphi, dx0
      dx0 = beta/(dx * dx)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      laplphi =
     & (b0(i+1,j )*(phi(i+1,j ) - phi(i ,j ))
     & - b0(i ,j )*(phi(i ,j ) - phi(i-1,j )))*dx0
     & +(b1(i ,j+1)*(phi(i ,j+1) - phi(i ,j ))
     & - b1(i ,j )*(phi(i ,j ) - phi(i ,j-1)))*dx0
      opphidir(i,j) = opphidir(i,j) + laplphi
      enddo
      enddo
      return
      end
      subroutine INCRAPPLYEBCO(
     & lhs
     & ,ilhslo0,ilhslo1
     & ,ilhshi0,ilhshi1
     & ,interiorflux
     & ,iinteriorfluxlo0,iinteriorfluxlo1
     & ,iinteriorfluxhi0,iinteriorfluxhi1
     & ,domainfluxlo
     & ,idomainfluxlolo0,idomainfluxlolo1
     & ,idomainfluxlohi0,idomainfluxlohi1
     & ,domainfluxhi
     & ,idomainfluxhilo0,idomainfluxhilo1
     & ,idomainfluxhihi0,idomainfluxhihi1
     & ,beta
     & ,dx
     & ,iloboxlo0,iloboxlo1
     & ,iloboxhi0,iloboxhi1
     & ,ihiboxlo0,ihiboxlo1
     & ,ihiboxhi0,ihiboxhi1
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
     & ,haslo
     & ,hashi
     & ,facedir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ilhslo0,ilhslo1
      integer ilhshi0,ilhshi1
      REAL*8 lhs(
     & ilhslo0:ilhshi0,
     & ilhslo1:ilhshi1)
      integer iinteriorfluxlo0,iinteriorfluxlo1
      integer iinteriorfluxhi0,iinteriorfluxhi1
      REAL*8 interiorflux(
     & iinteriorfluxlo0:iinteriorfluxhi0,
     & iinteriorfluxlo1:iinteriorfluxhi1)
      integer idomainfluxlolo0,idomainfluxlolo1
      integer idomainfluxlohi0,idomainfluxlohi1
      REAL*8 domainfluxlo(
     & idomainfluxlolo0:idomainfluxlohi0,
     & idomainfluxlolo1:idomainfluxlohi1)
      integer idomainfluxhilo0,idomainfluxhilo1
      integer idomainfluxhihi0,idomainfluxhihi1
      REAL*8 domainfluxhi(
     & idomainfluxhilo0:idomainfluxhihi0,
     & idomainfluxhilo1:idomainfluxhihi1)
      REAL*8 beta
      REAL*8 dx
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer haslo
      integer hashi
      integer facedir
      integer ii,i,jj,j
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      lhs(i,j) = lhs(i,j)
     $ +beta*
     $ (interiorflux(i+ii,j+jj)
     $ -interiorflux(i ,j ))/dx
      enddo
      enddo
      if(haslo .eq. 1) then
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         lhs(i,j) = lhs(i,j)
     $ + beta*
     $ (interiorflux(i+ii,j+jj)
     $ -domainfluxlo(i ,j ))/dx
      enddo
      enddo
      endif
      if(hashi .eq. 1) then
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         lhs(i,j) = lhs(i,j)
     $ + beta*
     $ (domainfluxhi(i+ii,j+jj)
     $ -interiorflux(i ,j ))/dx
      enddo
      enddo
      endif
      return
      end
      subroutine DECRINVRELCOEFEBCO(
     & relcoef
     & ,irelcoeflo0,irelcoeflo1
     & ,irelcoefhi0,irelcoefhi1
     & ,bcoef
     & ,ibcoeflo0,ibcoeflo1
     & ,ibcoefhi0,ibcoefhi1
     & ,beta
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dx
     & ,idir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer irelcoeflo0,irelcoeflo1
      integer irelcoefhi0,irelcoefhi1
      REAL*8 relcoef(
     & irelcoeflo0:irelcoefhi0,
     & irelcoeflo1:irelcoefhi1)
      integer ibcoeflo0,ibcoeflo1
      integer ibcoefhi0,ibcoefhi1
      REAL*8 bcoef(
     & ibcoeflo0:ibcoefhi0,
     & ibcoeflo1:ibcoefhi1)
      REAL*8 beta
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 dx
      integer idir
      integer i,j
      integer ii,jj
      ii = chf_id(idir, 0)
      jj = chf_id(idir, 1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      relcoef(i,j) = relcoef(i,j)
     $ - beta*(
     $ bcoef(i+ii,j+jj) +
     $ bcoef(i ,j ))/(dx*dx)
      enddo
      enddo
      return
      end
      subroutine INVERTLAMBDAEBCO(
     & lambda
     & ,ilambdalo0,ilambdalo1
     & ,ilambdahi0,ilambdahi1
     & ,safety
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilambdalo0,ilambdalo1
      integer ilambdahi0,ilambdahi1
      REAL*8 lambda(
     & ilambdalo0:ilambdahi0,
     & ilambdalo1:ilambdahi1)
      REAL*8 safety
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         lambda(i,j) = safety/(lambda(i,j))
      enddo
      enddo
      return
      end
      subroutine GETFLUXEBCO(
     & flux
     & ,ifluxlo0,ifluxlo1
     & ,ifluxhi0,ifluxhi1
     & ,bcoef
     & ,ibcoeflo0,ibcoeflo1
     & ,ibcoefhi0,ibcoefhi1
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,iopphiboxlo0,iopphiboxlo1
     & ,iopphiboxhi0,iopphiboxhi1
     & ,dx
     & ,idir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifluxlo0,ifluxlo1
      integer ifluxhi0,ifluxhi1
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1)
      integer ibcoeflo0,ibcoeflo1
      integer ibcoefhi0,ibcoefhi1
      REAL*8 bcoef(
     & ibcoeflo0:ibcoefhi0,
     & ibcoeflo1:ibcoefhi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer iopphiboxlo0,iopphiboxlo1
      integer iopphiboxhi0,iopphiboxhi1
      REAL*8 dx
      integer idir
      integer i,j
      integer ioff,joff
      REAL*8 oneoverdx
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      do j = iopphiboxlo1,iopphiboxhi1
      do i = iopphiboxlo0,iopphiboxhi0
      oneoverdx = bcoef(i,j)/dx
      flux(i,j) =
     $ oneoverdx*(
     $ phi(i ,j ) -
     $ phi(i-ioff,j-joff) )
      enddo
      enddo
      return
      end
      subroutine GSRBEBCO(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,lph
     & ,ilphlo0,ilphlo1
     & ,ilphhi0,ilphhi1
     & ,rhs
     & ,irhslo0,irhslo1
     & ,irhshi0,irhshi1
     & ,lam
     & ,ilamlo0,ilamlo1
     & ,ilamhi0,ilamhi1
     & ,icoloredboxlo0,icoloredboxlo1
     & ,icoloredboxhi0,icoloredboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer ilphlo0,ilphlo1
      integer ilphhi0,ilphhi1
      REAL*8 lph(
     & ilphlo0:ilphhi0,
     & ilphlo1:ilphhi1)
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1)
      integer ilamlo0,ilamlo1
      integer ilamhi0,ilamhi1
      REAL*8 lam(
     & ilamlo0:ilamhi0,
     & ilamlo1:ilamhi1)
      integer icoloredboxlo0,icoloredboxlo1
      integer icoloredboxhi0,icoloredboxhi1
      integer i,j
      REAL*8 phio, lamo, rhso, lpho
      do j = icoloredboxlo1,icoloredboxhi1,2
      do i = icoloredboxlo0,icoloredboxhi0,2
         phio = phi(i,j)
         lamo = lam(i,j)
         rhso = rhs(i,j)
         lpho = lph(i,j)
         phi(i,j) =
     $ phi(i,j) +
     $ lam(i,j)*(
     $ rhs(i,j) -
     $ lph(i,j))
      enddo
      enddo
      return
      end
