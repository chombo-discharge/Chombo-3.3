#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine NWOEBVTOPOINTLPH(
     &           lphi
     &           ,iv
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,acofab
     &           ,iacofablo0,iacofablo1
     &           ,iacofabhi0,iacofabhi1
     &           ,eta0fab
     &           ,ieta0fablo0,ieta0fablo1
     &           ,ieta0fabhi0,ieta0fabhi1
     &           ,eta1fab
     &           ,ieta1fablo0,ieta1fablo1
     &           ,ieta1fabhi0,ieta1fabhi1
     &           ,eta2fab
     &           ,ieta2fablo0,ieta2fablo1
     &           ,ieta2fabhi0,ieta2fabhi1
     &           ,lam0fab
     &           ,ilam0fablo0,ilam0fablo1
     &           ,ilam0fabhi0,ilam0fabhi1
     &           ,lam1fab
     &           ,ilam1fablo0,ilam1fablo1
     &           ,ilam1fabhi0,ilam1fabhi1
     &           ,lam2fab
     &           ,ilam2fablo0,ilam2fablo1
     &           ,ilam2fabhi0,ilam2fabhi1
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      REAL_T lphi(0:1)
      integer iv(0:1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer iacofablo0,iacofablo1
      integer iacofabhi0,iacofabhi1
      REAL_T acofab(
     &           iacofablo0:iacofabhi0,
     &           iacofablo1:iacofabhi1)
      integer ieta0fablo0,ieta0fablo1
      integer ieta0fabhi0,ieta0fabhi1
      REAL_T eta0fab(
     &           ieta0fablo0:ieta0fabhi0,
     &           ieta0fablo1:ieta0fabhi1)
      integer ieta1fablo0,ieta1fablo1
      integer ieta1fabhi0,ieta1fabhi1
      REAL_T eta1fab(
     &           ieta1fablo0:ieta1fabhi0,
     &           ieta1fablo1:ieta1fabhi1)
      integer ieta2fablo0,ieta2fablo1
      integer ieta2fabhi0,ieta2fabhi1
      REAL_T eta2fab(
     &           ieta2fablo0:ieta2fabhi0,
     &           ieta2fablo1:ieta2fabhi1)
      integer ilam0fablo0,ilam0fablo1
      integer ilam0fabhi0,ilam0fabhi1
      REAL_T lam0fab(
     &           ilam0fablo0:ilam0fabhi0,
     &           ilam0fablo1:ilam0fabhi1)
      integer ilam1fablo0,ilam1fablo1
      integer ilam1fabhi0,ilam1fabhi1
      REAL_T lam1fab(
     &           ilam1fablo0:ilam1fabhi0,
     &           ilam1fablo1:ilam1fabhi1)
      integer ilam2fablo0,ilam2fablo1
      integer ilam2fabhi0,ilam2fabhi1
      REAL_T lam2fab(
     &           ilam2fablo0:ilam2fabhi0,
     &           ilam2fablo1:ilam2fabhi1)
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      real_t divf(0:CH_SPACEDIM-1)
      real_t aphi(0:CH_SPACEDIM-1)
      real_t  etaL(0:CH_SPACEDIM-1)
      real_t  etaH(0:CH_SPACEDIM-1)
      real_t  lamL(0:CH_SPACEDIM-1)
      real_t  lamH(0:CH_SPACEDIM-1)
      real_t  fluxL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  fluxH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  divuL(0:CH_SPACEDIM-1)
      real_t  divuH(0:CH_SPACEDIM-1)
      integer i,j, facedir, derivdir, veldir
      integer iif,jjf
      integer iid,jjd
      
      i = iv(0)
      j = iv(1)
      
      etaL(0) =eta0fab(i   ,j   )
      etaL(1) =eta1fab(i   ,j   )
                                  
      etaH(0) =eta0fab(i+1 ,j   )
      etaH(1) =eta1fab(i   ,j+1 )
      
      lamL(0) =lam0fab(i   ,j   )
      lamL(1) =lam1fab(i   ,j   )
                                  
      lamH(0) =lam0fab(i+1 ,j   )
      lamH(1) =lam1fab(i   ,j+1 )
      do facedir = 0, CH_SPACEDIM-1
         
         iif = chf_id(facedir, 0)
         jjf = chf_id(facedir, 1)
         do derivdir = 0, CH_SPACEDIM-1
            
            iid = chf_id(derivdir, 0)
            jjd = chf_id(derivdir, 1)
            do veldir = 0, CH_SPACEDIM-1
               if(facedir .eq. derivdir) then
                  gphiH(veldir, derivdir ,facedir) = (phi(i+iid,j+jjd,veldir) - phi(i    ,j    ,veldir))/dx
                  gphiL(veldir, derivdir ,facedir) = (phi(i    ,j    ,veldir) - phi(i-iid,j-jjd,veldir))/dx
               else
                  gphiH(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phi(i+iid+iif,j+jjd+jjf,veldir) - phi(i-iid+iif,j-jjd+jjf,veldir) + 
     $                 phi(i+iid    ,j+jjd    ,veldir) - phi(i-iid    ,j-jjd    ,veldir)  )
                  gphiL(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phi(i+iid-iif,j+jjd-jjf,veldir) - phi(i-iid-iif,j-jjd-jjf,veldir) + 
     $                 phi(i+iid    ,j+jjd    ,veldir) - phi(i-iid    ,j-jjd    ,veldir)  )
               endif
            enddo
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         divuL(facedir) = zero
         divuH(facedir) = zero
         do veldir = 0, CH_SPACEDIM-1
            divuL(facedir)= divuL(facedir) + gphiL(veldir, veldir, facedir)
            divuH(facedir)= divuH(facedir) + gphiH(veldir, veldir, facedir)
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         do veldir = 0, CH_SPACEDIM-1
            fluxL(veldir, facedir) = etaL(facedir)*(gphiL(facedir, veldir, facedir) + gphiL(veldir, facedir, facedir))
            fluxH(veldir, facedir) = etaH(facedir)*(gphiH(facedir, veldir, facedir) + gphiH(veldir, facedir, facedir))
            if(veldir .eq. facedir) then
               fluxL(veldir, facedir) = fluxL(veldir, facedir) + lamL(facedir)*divuL(facedir)
               fluxH(veldir, facedir) = fluxH(veldir, facedir) + lamH(facedir)*divuH(facedir)
            endif
         enddo
      enddo
      do veldir = 0, CH_SPACEDIM-1
         divf(veldir) = zero
         do facedir = 0, CH_SPACEDIM-1
            divf(veldir) = divf(veldir)
     $           + (fluxH(veldir, facedir)-fluxL(veldir, facedir))/dx
         enddo
      enddo
      
      aphi(0) = acofab(i,j)*phi(i,j, 0)
      aphi(1) = acofab(i,j)*phi(i,j, 1)
      
      lphi(0) =  alpha*aphi(0) + beta*divf(0)
      lphi(1) =  alpha*aphi(1) + beta*divf(1)
      ch_flops=ch_flops + 4 + 6*CH_SPACEDIM + 12*(CH_SPACEDIM*CH_SPACEDIM)
      return
      end
      subroutine GSRBNWOEBVTOP(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,lambda
     &           ,ilambdalo0,ilambdalo1
     &           ,ilambdahi0,ilambdahi1
     &           ,nlambdacomp
     &           ,acofab
     &           ,iacofablo0,iacofablo1
     &           ,iacofabhi0,iacofabhi1
     &           ,eta0fab
     &           ,ieta0fablo0,ieta0fablo1
     &           ,ieta0fabhi0,ieta0fabhi1
     &           ,eta1fab
     &           ,ieta1fablo0,ieta1fablo1
     &           ,ieta1fabhi0,ieta1fabhi1
     &           ,eta2fab
     &           ,ieta2fablo0,ieta2fablo1
     &           ,ieta2fabhi0,ieta2fabhi1
     &           ,lam0fab
     &           ,ilam0fablo0,ilam0fablo1
     &           ,ilam0fabhi0,ilam0fabhi1
     &           ,lam1fab
     &           ,ilam1fablo0,ilam1fablo1
     &           ,ilam1fabhi0,ilam1fabhi1
     &           ,lam2fab
     &           ,ilam2fablo0,ilam2fablo1
     &           ,ilam2fabhi0,ilam2fabhi1
     &           ,dx
     &           ,alpha
     &           ,beta
     &           ,icoloredboxlo0,icoloredboxlo1
     &           ,icoloredboxhi0,icoloredboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


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
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1
      integer ilambdahi0,ilambdahi1
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1,
     &           0:nlambdacomp-1)
      integer iacofablo0,iacofablo1
      integer iacofabhi0,iacofabhi1
      REAL_T acofab(
     &           iacofablo0:iacofabhi0,
     &           iacofablo1:iacofabhi1)
      integer ieta0fablo0,ieta0fablo1
      integer ieta0fabhi0,ieta0fabhi1
      REAL_T eta0fab(
     &           ieta0fablo0:ieta0fabhi0,
     &           ieta0fablo1:ieta0fabhi1)
      integer ieta1fablo0,ieta1fablo1
      integer ieta1fabhi0,ieta1fabhi1
      REAL_T eta1fab(
     &           ieta1fablo0:ieta1fabhi0,
     &           ieta1fablo1:ieta1fabhi1)
      integer ieta2fablo0,ieta2fablo1
      integer ieta2fabhi0,ieta2fabhi1
      REAL_T eta2fab(
     &           ieta2fablo0:ieta2fabhi0,
     &           ieta2fablo1:ieta2fabhi1)
      integer ilam0fablo0,ilam0fablo1
      integer ilam0fabhi0,ilam0fabhi1
      REAL_T lam0fab(
     &           ilam0fablo0:ilam0fabhi0,
     &           ilam0fablo1:ilam0fabhi1)
      integer ilam1fablo0,ilam1fablo1
      integer ilam1fabhi0,ilam1fabhi1
      REAL_T lam1fab(
     &           ilam1fablo0:ilam1fabhi0,
     &           ilam1fablo1:ilam1fabhi1)
      integer ilam2fablo0,ilam2fablo1
      integer ilam2fabhi0,ilam2fabhi1
      REAL_T lam2fab(
     &           ilam2fablo0:ilam2fabhi0,
     &           ilam2fablo1:ilam2fabhi1)
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      integer icoloredboxlo0,icoloredboxlo1
      integer icoloredboxhi0,icoloredboxhi1
      REAL_T lphi(0:CH_SPACEDIM-1)
      real_t divf(0:CH_SPACEDIM-1)
      real_t aphi(0:CH_SPACEDIM-1)
      real_t  etaL(0:CH_SPACEDIM-1)
      real_t  etaH(0:CH_SPACEDIM-1)
      real_t  lamL(0:CH_SPACEDIM-1)
      real_t  lamH(0:CH_SPACEDIM-1)
      real_t  fluxL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  fluxH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  divuL(0:CH_SPACEDIM-1)
      real_t  divuH(0:CH_SPACEDIM-1)
      integer i,j, facedir, derivdir, veldir
      integer iif,jjf
      integer iid,jjd, ncolors
      
      do j = icoloredboxlo1,icoloredboxhi1,2
      do i = icoloredboxlo0,icoloredboxhi0,2

      
      etaL(0) =eta0fab(i   ,j   )
      etaL(1) =eta1fab(i   ,j   )
                                  
      etaH(0) =eta0fab(i+1 ,j   )
      etaH(1) =eta1fab(i   ,j+1 )
      
      lamL(0) =lam0fab(i   ,j   )
      lamL(1) =lam1fab(i   ,j   )
                                  
      lamH(0) =lam0fab(i+1 ,j   )
      lamH(1) =lam1fab(i   ,j+1 )
      do facedir = 0, CH_SPACEDIM-1
         
         iif = chf_id(facedir, 0)
         jjf = chf_id(facedir, 1)
         do derivdir = 0, CH_SPACEDIM-1
            
            iid = chf_id(derivdir, 0)
            jjd = chf_id(derivdir, 1)
            do veldir = 0, CH_SPACEDIM-1
               if(facedir .eq. derivdir) then
                  gphiH(veldir, derivdir ,facedir) = (phi(i+iid,j+jjd,veldir) - phi(i    ,j    ,veldir))/dx
                  gphiL(veldir, derivdir ,facedir) = (phi(i    ,j    ,veldir) - phi(i-iid,j-jjd,veldir))/dx
               else
                  gphiH(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phi(i+iid+iif,j+jjd+jjf,veldir) - phi(i-iid+iif,j-jjd+jjf,veldir) + 
     $                 phi(i+iid    ,j+jjd    ,veldir) - phi(i-iid    ,j-jjd    ,veldir)  )
                  gphiL(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phi(i+iid-iif,j+jjd-jjf,veldir) - phi(i-iid-iif,j-jjd-jjf,veldir) + 
     $                 phi(i+iid    ,j+jjd    ,veldir) - phi(i-iid    ,j-jjd    ,veldir)  )
               endif
            enddo
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         divuL(facedir) = zero
         divuH(facedir) = zero
         do veldir = 0, CH_SPACEDIM-1
            divuL(facedir)= divuL(facedir) + gphiL(veldir, veldir, facedir)
            divuH(facedir)= divuH(facedir) + gphiH(veldir, veldir, facedir)
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         do veldir = 0, CH_SPACEDIM-1
            fluxL(veldir, facedir) = etaL(facedir)*(gphiL(facedir, veldir, facedir) + gphiL(veldir, facedir, facedir))
            fluxH(veldir, facedir) = etaH(facedir)*(gphiH(facedir, veldir, facedir) + gphiH(veldir, facedir, facedir))
            if(veldir .eq. facedir) then
               fluxL(veldir, facedir) = fluxL(veldir, facedir) + lamL(facedir)*divuL(facedir)
               fluxH(veldir, facedir) = fluxH(veldir, facedir) + lamH(facedir)*divuH(facedir)
            endif
         enddo
      enddo
      do veldir = 0, CH_SPACEDIM-1
         divf(veldir) = zero
         do facedir = 0, CH_SPACEDIM-1
            divf(veldir) = divf(veldir)
     $           + (fluxH(veldir, facedir)-fluxL(veldir, facedir))/dx
         enddo
      enddo
      
      aphi(0) = acofab(i,j)*phi(i,j, 0)
      aphi(1) = acofab(i,j)*phi(i,j, 1)
      
      lphi(0) =  alpha*aphi(0) + beta*divf(0)
      lphi(1) =  alpha*aphi(1) + beta*divf(1)
      
      phi(i,j, 0) =  phi(i,j,0) + lambda(i,j,0)*(rhs(i,j,0) - lphi(0))
      phi(i,j, 1) =  phi(i,j,1) + lambda(i,j,1)*(rhs(i,j,1) - lphi(1))
      
      enddo
      enddo
      ncolors = 2 *2
      ch_flops=ch_flops+(icoloredboxhi0- icoloredboxlo0+1)*(icoloredboxhi1- icoloredboxlo1+1)*(4 + 9*CH_SPACEDIM + 12*(CH_SPACEDIM*CH_SPACEDIM))/ncolors
      return
      end
      subroutine CELLGRADEBVTOP(
     &           grad
     &           ,igradlo0,igradlo1
     &           ,igradhi0,igradhi1
     &           ,vel
     &           ,ivello0,ivello1
     &           ,ivelhi0,ivelhi1
     &           ,dx
     &           ,iloboxlo0,iloboxlo1
     &           ,iloboxhi0,iloboxhi1
     &           ,ihiboxlo0,ihiboxlo1
     &           ,ihiboxhi0,ihiboxhi1
     &           ,icenterboxlo0,icenterboxlo1
     &           ,icenterboxhi0,icenterboxhi1
     &           ,haslo
     &           ,hashi
     &           ,divdir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradlo0,igradlo1
      integer igradhi0,igradhi1
      REAL_T grad(
     &           igradlo0:igradhi0,
     &           igradlo1:igradhi1)
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1)
      REAL_T dx
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer haslo
      integer hashi
      integer divdir
      integer ii,i,jj,j
      
      ii = chf_id(divdir, 0)
      jj = chf_id(divdir, 1)
      
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

      grad(i,j) =
     $     (    vel(i+ii,j+jj)
     $     -    vel(i-ii,j-jj) )/(two*dx)
      
      enddo
      enddo
      if(haslo.eq.1) then
         
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

         grad(i,j) =
     $        (    vel(i+ii,j+jj)
     $        -    vel(i   ,j   ) )/(dx)
         
      enddo
      enddo
      endif
      if(hashi.eq.1) then
         
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

         grad(i,j) =
     $        (    vel(i   ,j   )
     $        -    vel(i-ii,j-jj) )/(dx)
         
      enddo
      enddo
      endif
      return
      end
      subroutine INCRAPPLYEBVTOP(
     &           lhs
     &           ,ilhslo0,ilhslo1
     &           ,ilhshi0,ilhshi1
     &           ,interiorflux
     &           ,iinteriorfluxlo0,iinteriorfluxlo1
     &           ,iinteriorfluxhi0,iinteriorfluxhi1
     &           ,domainfluxlo
     &           ,idomainfluxlolo0,idomainfluxlolo1
     &           ,idomainfluxlohi0,idomainfluxlohi1
     &           ,domainfluxhi
     &           ,idomainfluxhilo0,idomainfluxhilo1
     &           ,idomainfluxhihi0,idomainfluxhihi1
     &           ,beta
     &           ,dx
     &           ,iloboxlo0,iloboxlo1
     &           ,iloboxhi0,iloboxhi1
     &           ,ihiboxlo0,ihiboxlo1
     &           ,ihiboxhi0,ihiboxhi1
     &           ,icenterboxlo0,icenterboxlo1
     &           ,icenterboxhi0,icenterboxhi1
     &           ,haslo
     &           ,hashi
     &           ,facedir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ilhslo0,ilhslo1
      integer ilhshi0,ilhshi1
      REAL_T lhs(
     &           ilhslo0:ilhshi0,
     &           ilhslo1:ilhshi1)
      integer iinteriorfluxlo0,iinteriorfluxlo1
      integer iinteriorfluxhi0,iinteriorfluxhi1
      REAL_T interiorflux(
     &           iinteriorfluxlo0:iinteriorfluxhi0,
     &           iinteriorfluxlo1:iinteriorfluxhi1)
      integer idomainfluxlolo0,idomainfluxlolo1
      integer idomainfluxlohi0,idomainfluxlohi1
      REAL_T domainfluxlo(
     &           idomainfluxlolo0:idomainfluxlohi0,
     &           idomainfluxlolo1:idomainfluxlohi1)
      integer idomainfluxhilo0,idomainfluxhilo1
      integer idomainfluxhihi0,idomainfluxhihi1
      REAL_T domainfluxhi(
     &           idomainfluxhilo0:idomainfluxhihi0,
     &           idomainfluxhilo1:idomainfluxhihi1)
      REAL_T beta
      REAL_T dx
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
     $     +beta*
     $     (interiorflux(i+ii,j+jj)
     $     -interiorflux(i   ,j   ))/dx
      
      enddo
      enddo
      if(haslo .eq. 1) then
         
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

         lhs(i,j) = lhs(i,j)
     $        + beta*
     $        (interiorflux(i+ii,j+jj)
     $        -domainfluxlo(i   ,j   ))/dx
         
      enddo
      enddo
      endif
      if(hashi .eq. 1) then
         
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

         lhs(i,j) = lhs(i,j)
     $        + beta*
     $        (domainfluxhi(i+ii,j+jj)
     $        -interiorflux(i   ,j   ))/dx
         
      enddo
      enddo
      endif
      return
      end
      subroutine NORMALGRADVISCDIRCH(
     &           gradvelface
     &           ,igradvelfacelo0,igradvelfacelo1
     &           ,igradvelfacehi0,igradvelfacehi1
     &           ,velcomp
     &           ,ivelcomplo0,ivelcomplo1
     &           ,ivelcomphi0,ivelcomphi1
     &           ,velvalu
     &           ,ivelvalulo0,ivelvalulo1
     &           ,ivelvaluhi0,ivelvaluhi1
     &           ,ifaceboxlo0,ifaceboxlo1
     &           ,ifaceboxhi0,ifaceboxhi1
     &           ,dx
     &           ,iside
     &           ,divdir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradvelfacelo0,igradvelfacelo1
      integer igradvelfacehi0,igradvelfacehi1
      REAL_T gradvelface(
     &           igradvelfacelo0:igradvelfacehi0,
     &           igradvelfacelo1:igradvelfacehi1)
      integer ivelcomplo0,ivelcomplo1
      integer ivelcomphi0,ivelcomphi1
      REAL_T velcomp(
     &           ivelcomplo0:ivelcomphi0,
     &           ivelcomplo1:ivelcomphi1)
      integer ivelvalulo0,ivelvalulo1
      integer ivelvaluhi0,ivelvaluhi1
      REAL_T velvalu(
     &           ivelvalulo0:ivelvaluhi0,
     &           ivelvalulo1:ivelvaluhi1)
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      REAL_T dx
      integer iside
      integer divdir
      integer i,id,j,jd
      real_t val0, val1, val2, denom
      
      id = chf_id(divdir, 0)
      jd = chf_id(divdir, 1)
      if(iside.eq.-1) then
         
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

         val0 = velvalu(i   ,j   )
         val1 = velcomp(i   ,j   )
         val2 = velcomp(i+id,j+jd)
         denom = half*dx
         gradvelface(i,j) = (val1-val0)/denom
         
      enddo
      enddo
      else
         
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

         val0 = velvalu(i     ,j     )
         val1 = velcomp(i-  id,j-  jd)
         val2 = velcomp(i-2*id,j-2*jd)
         denom = -half*dx
         gradvelface(i,j) = (val1-val0)/denom
         
      enddo
      enddo
      endif
      return
      end
      subroutine SLIPWALLGRAD(
     &           gradvelface
     &           ,igradvelfacelo0,igradvelfacelo1
     &           ,igradvelfacehi0,igradvelfacehi1
     &           ,velcomp
     &           ,ivelcomplo0,ivelcomplo1
     &           ,ivelcomphi0,ivelcomphi1
     &           ,ifaceboxlo0,ifaceboxlo1
     &           ,ifaceboxhi0,ifaceboxhi1
     &           ,dx
     &           ,veldir
     &           ,divdir
     &           ,iside
     &           ,iloboxlo0,iloboxlo1
     &           ,iloboxhi0,iloboxhi1
     &           ,ihiboxlo0,ihiboxlo1
     &           ,ihiboxhi0,ihiboxhi1
     &           ,icenterboxlo0,icenterboxlo1
     &           ,icenterboxhi0,icenterboxhi1
     &           ,haslo
     &           ,hashi
     &           ,facedir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradvelfacelo0,igradvelfacelo1
      integer igradvelfacehi0,igradvelfacehi1
      REAL_T gradvelface(
     &           igradvelfacelo0:igradvelfacehi0,
     &           igradvelfacelo1:igradvelfacehi1)
      integer ivelcomplo0,ivelcomplo1
      integer ivelcomphi0,ivelcomphi1
      REAL_T velcomp(
     &           ivelcomplo0:ivelcomphi0,
     &           ivelcomplo1:ivelcomphi1)
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      REAL_T dx
      integer veldir
      integer divdir
      integer iside
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer haslo
      integer hashi
      integer facedir
      integer i,id,iv,j,jd,jv
      real_t  val1, val2
      
      id = chf_id(facedir, 0)
      jd = chf_id(facedir, 1)
      
      iv = chf_id(divdir, 0)
      jv = chf_id(divdir, 1)
      if(facedir.eq.veldir) then
         if(iside.eq.-1) then
            
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

            val1 = velcomp(i   ,j   )
            gradvelface(i,j) = two*val1/dx
            
      enddo
      enddo
         else
            
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

            val1 = velcomp(i-  id,j-  jd)
            gradvelface(i,j) = -two*val1/dx
            
      enddo
      enddo
         endif
      else
         if(iside.eq.-1) then
            
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

            val1 = velcomp(i+iv,j+jv)
            val2 = velcomp(i-iv,j-jv)
            gradvelface(i,j) = (val1-val2)/(two*dx)
            
      enddo
      enddo
            if(haslo.eq.1) then
               
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

               val1 = velcomp(i+iv,j+jv)
               val2 = velcomp(i   ,j   )
               gradvelface(i,j) = (val1-val2)/(dx)
               
      enddo
      enddo
            endif
            if(hashi.eq.1) then
               
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

               val1 = velcomp(i   ,j   )
               val2 = velcomp(i-iv,j-jv)
               gradvelface(i,j) = (val1-val2)/(dx)
               
      enddo
      enddo
            endif
         else
         
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

            val1 = velcomp(i-id+iv,j-jd+jv)
            val2 = velcomp(i-id-iv,j-jd-jv)
            gradvelface(i,j) = (val1-val2)/(two*dx)
         
      enddo
      enddo
         if(haslo.eq.1) then
            
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

            val1 = velcomp(i-id+iv,j-jd+jv)
            val2 = velcomp(i-id   ,j-jd   )
            gradvelface(i,j) = (val1-val2)/(dx)
            
      enddo
      enddo
         endif
         if(hashi.eq.1) then
            
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

            val1 = velcomp(i-id   ,j-jd   )
            val2 = velcomp(i-id-iv,j-jd-jv)
            gradvelface(i,j) = (val1-val2)/(dx)
            
      enddo
      enddo
         endif
         endif
      endif
      return
      end
      subroutine VELDOTSIGMA(
     &           veldotsig
     &           ,iveldotsiglo0,iveldotsiglo1
     &           ,iveldotsighi0,iveldotsighi1
     &           ,nveldotsigcomp
     &           ,vel
     &           ,ivello0,ivello1
     &           ,ivelhi0,ivelhi1
     &           ,nvelcomp
     &           ,sig
     &           ,isiglo0,isiglo1
     &           ,isighi0,isighi1
     &           ,nsigcomp
     &           ,ifaceboxlo0,ifaceboxlo1
     &           ,ifaceboxhi0,ifaceboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nveldotsigcomp
      integer iveldotsiglo0,iveldotsiglo1
      integer iveldotsighi0,iveldotsighi1
      REAL_T veldotsig(
     &           iveldotsiglo0:iveldotsighi0,
     &           iveldotsiglo1:iveldotsighi1,
     &           0:nveldotsigcomp-1)
      integer nvelcomp
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1,
     &           0:nvelcomp-1)
      integer nsigcomp
      integer isiglo0,isiglo1
      integer isighi0,isighi1
      REAL_T sig(
     &           isiglo0:isighi0,
     &           isiglo1:isighi1,
     &           0:nsigcomp-1)
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      integer i,j
      integer icomp
      real_t sum
      
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

      sum = zero
      do icomp = 0, CH_SPACEDIM-1
         sum = sum + vel(i,j, icomp)*sig(i,j, icomp)
      enddo
      veldotsig(i,j, 0) = sum
      
      enddo
      enddo
      return
      end
      subroutine INCRPOINTDOTPROD(
     &           sigmadotgradu
     &           ,isigmadotgradulo0,isigmadotgradulo1
     &           ,isigmadotgraduhi0,isigmadotgraduhi1
     &           ,gradu
     &           ,igradulo0,igradulo1
     &           ,igraduhi0,igraduhi1
     &           ,sigma
     &           ,isigmalo0,isigmalo1
     &           ,isigmahi0,isigmahi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isigmadotgradulo0,isigmadotgradulo1
      integer isigmadotgraduhi0,isigmadotgraduhi1
      REAL_T sigmadotgradu(
     &           isigmadotgradulo0:isigmadotgraduhi0,
     &           isigmadotgradulo1:isigmadotgraduhi1)
      integer igradulo0,igradulo1
      integer igraduhi0,igraduhi1
      REAL_T gradu(
     &           igradulo0:igraduhi0,
     &           igradulo1:igraduhi1)
      integer isigmalo0,isigmalo1
      integer isigmahi0,isigmahi1
      REAL_T sigma(
     &           isigmalo0:isigmahi0,
     &           isigmalo1:isigmahi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      real_t incr
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

      incr =  gradu(i,j)*sigma(i,j)
      sigmadotgradu(i,j) = sigmadotgradu(i,j) + incr
      
      enddo
      enddo
      return
      end
      subroutine INCRCCSIGMALAMBDA(
     &           sigma
     &           ,isigmalo0,isigmalo1
     &           ,isigmahi0,isigmahi1
     &           ,gradu
     &           ,igradulo0,igradulo1
     &           ,igraduhi0,igraduhi1
     &           ,ngraducomp
     &           ,lambda
     &           ,ilambdalo0,ilambdalo1
     &           ,ilambdahi0,ilambdahi1
     &           ,igridlo0,igridlo1
     &           ,igridhi0,igridhi1
     &           ,diagcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isigmalo0,isigmalo1
      integer isigmahi0,isigmahi1
      REAL_T sigma(
     &           isigmalo0:isigmahi0,
     &           isigmalo1:isigmahi1)
      integer ngraducomp
      integer igradulo0,igradulo1
      integer igraduhi0,igraduhi1
      REAL_T gradu(
     &           igradulo0:igraduhi0,
     &           igradulo1:igraduhi1,
     &           0:ngraducomp-1)
      integer ilambdalo0,ilambdalo1
      integer ilambdahi0,ilambdahi1
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1)
      integer igridlo0,igridlo1
      integer igridhi0,igridhi1
      integer diagcomp
      integer i,j
      
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

        sigma(i,j) = sigma(i,j)
     $     + lambda(i,j)*gradu(i,j, diagcomp)
      
      enddo
      enddo
      return
      end
      subroutine INCRCCSIGMAETA(
     &           sigma
     &           ,isigmalo0,isigmalo1
     &           ,isigmahi0,isigmahi1
     &           ,gradu
     &           ,igradulo0,igradulo1
     &           ,igraduhi0,igraduhi1
     &           ,ngraducomp
     &           ,eta
     &           ,ietalo0,ietalo1
     &           ,ietahi0,ietahi1
     &           ,igridlo0,igridlo1
     &           ,igridhi0,igridhi1
     &           ,gradcomp
     &           ,trancomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isigmalo0,isigmalo1
      integer isigmahi0,isigmahi1
      REAL_T sigma(
     &           isigmalo0:isigmahi0,
     &           isigmalo1:isigmahi1)
      integer ngraducomp
      integer igradulo0,igradulo1
      integer igraduhi0,igraduhi1
      REAL_T gradu(
     &           igradulo0:igraduhi0,
     &           igradulo1:igraduhi1,
     &           0:ngraducomp-1)
      integer ietalo0,ietalo1
      integer ietahi0,ietahi1
      REAL_T eta(
     &           ietalo0:ietahi0,
     &           ietalo1:ietahi1)
      integer igridlo0,igridlo1
      integer igridhi0,igridhi1
      integer gradcomp
      integer trancomp
      integer i,j
      
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

      sigma(i,j) = sigma(i,j)
     $     +  eta(i,j)*
     $     (gradu(i,j, gradcomp)
     $     +gradu(i,j, trancomp))
      
      enddo
      enddo
      return
      end
      subroutine INCRFACETOCELLAVERAGE(
     &           cellcoef
     &           ,icellcoeflo0,icellcoeflo1
     &           ,icellcoefhi0,icellcoefhi1
     &           ,facecoef
     &           ,ifacecoeflo0,ifacecoeflo1
     &           ,ifacecoefhi0,ifacecoefhi1
     &           ,facedir
     &           ,igridlo0,igridlo1
     &           ,igridhi0,igridhi1
     &           ,nfacesper
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer icellcoeflo0,icellcoeflo1
      integer icellcoefhi0,icellcoefhi1
      REAL_T cellcoef(
     &           icellcoeflo0:icellcoefhi0,
     &           icellcoeflo1:icellcoefhi1)
      integer ifacecoeflo0,ifacecoeflo1
      integer ifacecoefhi0,ifacecoefhi1
      REAL_T facecoef(
     &           ifacecoeflo0:ifacecoefhi0,
     &           ifacecoeflo1:ifacecoefhi1)
      integer facedir
      integer igridlo0,igridlo1
      integer igridhi0,igridhi1
      REAL_T nfacesper
      integer ii,i,jj,j
      
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

        cellcoef(i,j) = cellcoef(i,j)   +
     $     (facecoef(i+ii,j+jj)
     $     +facecoef(i   ,j   ))/nfacesper
      
      enddo
      enddo
      return
      end
