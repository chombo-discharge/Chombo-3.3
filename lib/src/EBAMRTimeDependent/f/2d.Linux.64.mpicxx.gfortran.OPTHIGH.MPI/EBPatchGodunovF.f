      subroutine MINFLAT(
     & flattening
     & ,iflatteninglo0,iflatteninglo1
     & ,iflatteninghi0,iflatteninghi1
     & ,zetadir
     & ,izetadirlo0,izetadirlo1
     & ,izetadirhi0,izetadirhi1
     & ,nzetadircomp
     & ,du
     & ,idulo0,idulo1
     & ,iduhi0,iduhi1
     & ,nducomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iflatteninglo0,iflatteninglo1
      integer iflatteninghi0,iflatteninghi1
      REAL*8 flattening(
     & iflatteninglo0:iflatteninghi0,
     & iflatteninglo1:iflatteninghi1)
      integer nzetadircomp
      integer izetadirlo0,izetadirlo1
      integer izetadirhi0,izetadirhi1
      REAL*8 zetadir(
     & izetadirlo0:izetadirhi0,
     & izetadirlo1:izetadirhi1,
     & 0:nzetadircomp-1)
      integer nducomp
      integer idulo0,idulo1
      integer iduhi0,iduhi1
      REAL*8 du(
     & idulo0:iduhi0,
     & idulo1:iduhi1,
     & 0:nducomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      integer iv
      REAL*8 sumdu,minflattot,minzetadir
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      sumdu = (0.0d0)
      do iv = 0,nducomp - 1
         sumdu = sumdu + du(i,j,iv)
      enddo
      if (sumdu .lt. (0.0d0)) then
         minflattot = zetadir(i,j,0)
         do iv = 1,nducomp - 1
            minzetadir = zetadir(i,j,iv)
            minflattot = min(minflattot,minzetadir)
         enddo
         flattening(i,j) = minflattot
      else
         flattening(i,j) = (1.0d0)
      endif
      enddo
      enddo
      return
      end
      subroutine GETDPTWO(
     & delta2p
     & ,idelta2plo0,idelta2plo1
     & ,idelta2phi0,idelta2phi1
     & ,delta1p
     & ,idelta1plo0,idelta1plo1
     & ,idelta1phi0,idelta1phi1
     & ,idir
     & ,iloboxlo0,iloboxlo1
     & ,iloboxhi0,iloboxhi1
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1
     & ,ihiboxhi0,ihiboxhi1
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idelta2plo0,idelta2plo1
      integer idelta2phi0,idelta2phi1
      REAL*8 delta2p(
     & idelta2plo0:idelta2phi0,
     & idelta2plo1:idelta2phi1)
      integer idelta1plo0,idelta1plo1
      integer idelta1phi0,idelta1phi1
      REAL*8 delta1p(
     & idelta1plo0:idelta1phi0,
     & idelta1plo1:idelta1phi1)
      integer idir
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer haslo
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer hashi
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer i ,j
      integer ioff,joff
      REAL*8 dp1hi, dp1lo
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      dp1hi = delta1p(i+ioff,j+joff)
      dp1lo = delta1p(i-ioff,j-joff)
      delta2p(i,j) = dp1hi + dp1lo
      enddo
      enddo
      if (haslo .eq. 1) then
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         dp1hi = delta1p(i+ioff,j+joff)
         dp1lo = delta1p(i ,j )
         delta2p(i,j) = dp1hi + dp1lo
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         dp1hi = delta1p(i ,j )
         dp1lo = delta1p(i-ioff,j-joff)
         delta2p(i,j) = dp1hi + dp1lo
      enddo
      enddo
      endif
      return
      end
      subroutine GETFLAT(
     & zetatwiddle
     & ,izetatwiddlelo0,izetatwiddlelo1
     & ,izetatwiddlehi0,izetatwiddlehi1
     & ,delta1p
     & ,idelta1plo0,idelta1plo1
     & ,idelta1phi0,idelta1phi1
     & ,delta2p
     & ,idelta2plo0,idelta2plo1
     & ,idelta2phi0,idelta2phi1
     & ,bulkmin
     & ,ibulkminlo0,ibulkminlo1
     & ,ibulkminhi0,ibulkminhi1
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer izetatwiddlelo0,izetatwiddlelo1
      integer izetatwiddlehi0,izetatwiddlehi1
      REAL*8 zetatwiddle(
     & izetatwiddlelo0:izetatwiddlehi0,
     & izetatwiddlelo1:izetatwiddlehi1)
      integer idelta1plo0,idelta1plo1
      integer idelta1phi0,idelta1phi1
      REAL*8 delta1p(
     & idelta1plo0:idelta1phi0,
     & idelta1plo1:idelta1phi1)
      integer idelta2plo0,idelta2plo1
      integer idelta2phi0,idelta2phi1
      REAL*8 delta2p(
     & idelta2plo0:idelta2phi0,
     & idelta2plo1:idelta2phi1)
      integer ibulkminlo0,ibulkminlo1
      integer ibulkminhi0,ibulkminhi1
      REAL*8 bulkmin(
     & ibulkminlo0:ibulkminhi0,
     & ibulkminlo1:ibulkminhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      REAL*8 d,r0,r1,ratio,strength
      REAL*8 d1pvof, d2pvof, smallp
      data d /0.33d0/
      data r0 /0.75d0/
      data r1 /0.85d0/
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      strength = abs(delta1p(i,j)/bulkmin(i,j))
      smallp = 1.0e-7
      if (strength .ge. d) then
         d1pvof = abs(delta1p(i,j))
         d2pvof = max(abs(delta2p(i,j)),smallp)
         ratio = d1pvof/d2pvof
         if (ratio .le. r0) then
            zetatwiddle(i,j) = (1.0d0)
         else if (ratio .ge. r1) then
            zetatwiddle(i,j) = (0.0d0)
         else
            zetatwiddle(i,j) = (1.0d0) - (ratio - r0)/(r1 - r0)
         endif
      else
         zetatwiddle(i,j) = (1.0d0)
      endif
      enddo
      enddo
      return
      end
      subroutine GETGRAD(
     & du
     & ,idulo0,idulo1
     & ,iduhi0,iduhi1
     & ,u
     & ,iulo0,iulo1
     & ,iuhi0,iuhi1
     & ,idir
     & ,iloboxlo0,iloboxlo1
     & ,iloboxhi0,iloboxhi1
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1
     & ,ihiboxhi0,ihiboxhi1
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
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
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer haslo
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer hashi
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer i ,j
      integer ioff,joff
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      du(i,j) = (0.500d0)*
     & ( u(i+ioff,j+joff)
     & - u(i-ioff,j-joff))
      enddo
      enddo
      if (haslo .eq. 1) then
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         du(i,j) =
     & ( u(i+ioff,j+joff)
     & - u(i ,j ))
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         du(i,j) =
     & ( u(i ,j )
     & - u(i-ioff,j-joff))
      enddo
      enddo
      endif
      return
      end
      subroutine GETRELATIVEGRAD(
     & du
     & ,idulo0,idulo1
     & ,iduhi0,iduhi1
     & ,u
     & ,iulo0,iulo1
     & ,iuhi0,iuhi1
     & ,idir
     & ,iloboxlo0,iloboxlo1
     & ,iloboxhi0,iloboxhi1
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1
     & ,ihiboxhi0,ihiboxhi1
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
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
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer haslo
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer hashi
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer i ,j
      integer ioff,joff
      REAL*8 diff, ave
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      diff = (0.500d0)*
     & ( u(i+ioff,j+joff)
     & - u(i-ioff,j-joff))
      ave = (0.500d0)*
     & ( u(i+ioff,j+joff)
     & + u(i-ioff,j-joff))
      du(i,j) = diff/ave
      enddo
      enddo
      if (haslo .eq. 1) then
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         diff =
     & ( u(i+ioff,j+joff)
     & - u(i ,j ))
         ave =(0.500d0)*
     & ( u(i+ioff,j+joff)
     & + u(i ,j ))
         du(i,j) = diff/ave
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         diff =
     & ( u(i ,j )
     & - u(i-ioff,j-joff))
         ave = (0.500d0)*
     & ( u(i ,j )
     & + u(i-ioff,j-joff))
         du(i,j) = diff/ave
      enddo
      enddo
      endif
      return
      end
      subroutine MAGNITUDE(
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
      integer i,j
      integer iv
      REAL*8 cur,sum
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      sum = (0.0d0)
      do iv = 0,ndatacomp-1
         cur = data(i,j,iv)
         sum = sum + cur*cur
      enddo
      magdata(i,j) = sqrt(sum)
      enddo
      enddo
      return
      end
      subroutine MIN3PTS(
     & mindata
     & ,imindatalo0,imindatalo1
     & ,imindatahi0,imindatahi1
     & ,data
     & ,idatalo0,idatalo1
     & ,idatahi0,idatahi1
     & ,idir
     & ,iloboxlo0,iloboxlo1
     & ,iloboxhi0,iloboxhi1
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1
     & ,ihiboxhi0,ihiboxhi1
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer imindatalo0,imindatalo1
      integer imindatahi0,imindatahi1
      REAL*8 mindata(
     & imindatalo0:imindatahi0,
     & imindatalo1:imindatahi1)
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1)
      integer idir
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer haslo
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer hashi
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer i ,j
      integer ioff,joff
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      mindata(i,j) = min(
     & data(i ,j ),
     & data(i+ioff,j+joff),
     & data(i-ioff,j-joff))
      enddo
      enddo
      if (haslo .ne. 0) then
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         mindata(i,j) = min(
     & data(i ,j ),
     & data(i+ioff,j+joff))
      enddo
      enddo
      endif
      if (hashi .ne. 0) then
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         mindata(i,j) = min(
     & data(i ,j ),
     & data(i-ioff,j-joff))
      enddo
      enddo
      endif
      return
      end
      subroutine SECONDSLOPEDIFFS(
     & deltawc
     & ,ideltawclo0,ideltawclo1
     & ,ideltawchi0,ideltawchi1
     & ,ndeltawccomp
     & ,deltawl
     & ,ideltawllo0,ideltawllo1
     & ,ideltawlhi0,ideltawlhi1
     & ,ndeltawlcomp
     & ,deltawr
     & ,ideltawrlo0,ideltawrlo1
     & ,ideltawrhi0,ideltawrhi1
     & ,ndeltawrcomp
     & ,w
     & ,iwlo0,iwlo1
     & ,iwhi0,iwhi1
     & ,nwcomp
     & ,numslopes
     & ,idir
     & ,iloboxlo0,iloboxlo1
     & ,iloboxhi0,iloboxhi1
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1
     & ,ihiboxhi0,ihiboxhi1
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndeltawccomp
      integer ideltawclo0,ideltawclo1
      integer ideltawchi0,ideltawchi1
      REAL*8 deltawc(
     & ideltawclo0:ideltawchi0,
     & ideltawclo1:ideltawchi1,
     & 0:ndeltawccomp-1)
      integer ndeltawlcomp
      integer ideltawllo0,ideltawllo1
      integer ideltawlhi0,ideltawlhi1
      REAL*8 deltawl(
     & ideltawllo0:ideltawlhi0,
     & ideltawllo1:ideltawlhi1,
     & 0:ndeltawlcomp-1)
      integer ndeltawrcomp
      integer ideltawrlo0,ideltawrlo1
      integer ideltawrhi0,ideltawrhi1
      REAL*8 deltawr(
     & ideltawrlo0:ideltawrhi0,
     & ideltawrlo1:ideltawrhi1,
     & 0:ndeltawrcomp-1)
      integer nwcomp
      integer iwlo0,iwlo1
      integer iwhi0,iwhi1
      REAL*8 w(
     & iwlo0:iwhi0,
     & iwlo1:iwhi1,
     & 0:nwcomp-1)
      integer numslopes
      integer idir
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer haslo
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer hashi
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer i ,j ,lvar
      integer ioff,joff
      REAL*8 dwr,dwl
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      do lvar = 0, numslopes - 1
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
         dwr = w(i+ioff,j+joff,lvar)
     & - w(i ,j ,lvar)
         dwl = w(i ,j ,lvar)
     & - w(i-ioff,j-joff,lvar)
         deltawr(i,j,lvar) = dwr
         deltawl(i,j,lvar) = dwl
         deltawc(i,j,lvar) = (0.500d0)*(dwr + dwl)
      enddo
      enddo
         if (haslo .ne. 0) then
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            dwr = w(i+ioff,j+joff,lvar)
     & - w(i ,j ,lvar)
            deltawc(i,j,lvar) = dwr
            deltawl(i,j,lvar) = dwr
            deltawr(i,j,lvar) = dwr
      enddo
      enddo
         endif
         if (hashi .ne. 0) then
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            dwl = w(i ,j ,lvar)
     & - w(i-ioff,j-joff,lvar)
            deltawc(i,j,lvar) = dwl
            deltawl(i,j,lvar) = dwl
            deltawr(i,j,lvar) = dwl
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine FORTHSLOPEDIFFS(
     & delta4wc
     & ,idelta4wclo0,idelta4wclo1
     & ,idelta4wchi0,idelta4wchi1
     & ,ndelta4wccomp
     & ,w
     & ,iwlo0,iwlo1
     & ,iwhi0,iwhi1
     & ,nwcomp
     & ,delta2w
     & ,idelta2wlo0,idelta2wlo1
     & ,idelta2whi0,idelta2whi1
     & ,ndelta2wcomp
     & ,numslopes
     & ,idir
     & ,iloboxlo0,iloboxlo1
     & ,iloboxhi0,iloboxhi1
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1
     & ,ihiboxhi0,ihiboxhi1
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndelta4wccomp
      integer idelta4wclo0,idelta4wclo1
      integer idelta4wchi0,idelta4wchi1
      REAL*8 delta4wc(
     & idelta4wclo0:idelta4wchi0,
     & idelta4wclo1:idelta4wchi1,
     & 0:ndelta4wccomp-1)
      integer nwcomp
      integer iwlo0,iwlo1
      integer iwhi0,iwhi1
      REAL*8 w(
     & iwlo0:iwhi0,
     & iwlo1:iwhi1,
     & 0:nwcomp-1)
      integer ndelta2wcomp
      integer idelta2wlo0,idelta2wlo1
      integer idelta2whi0,idelta2whi1
      REAL*8 delta2w(
     & idelta2wlo0:idelta2whi0,
     & idelta2wlo1:idelta2whi1,
     & 0:ndelta2wcomp-1)
      integer numslopes
      integer idir
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer haslo
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer hashi
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer i ,j ,lvar
      integer ioff,joff
      REAL*8 dwr,dwl, vall, slol, valr, slor
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      do lvar = 0, numslopes - 1
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
         valr = w(i+ioff,j+joff,lvar)
         slor = delta2w(i+ioff,j+joff,lvar)
         vall = w(i-ioff,j-joff,lvar)
         slol = delta2w(i-ioff,j-joff,lvar)
         dwl = vall + (0.250d0)*slol
         dwr = valr - (0.250d0)*slor
         delta4wc(i,j,lvar) = (2.000d0 / 3.000d0)*(dwr - dwl)
      enddo
      enddo
         if (haslo .ne. 0) then
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            delta4wc(i,j,lvar) = delta2w(i,j,lvar)
      enddo
      enddo
         endif
         if (hashi .ne. 0) then
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            delta4wc(i,j,lvar) = delta2w(i,j,lvar)
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine APPLYFLAT(
     & dw
     & ,idwlo0,idwlo1
     & ,idwhi0,idwhi1
     & ,ndwcomp
     & ,flattening
     & ,iflatteninglo0,iflatteninglo1
     & ,iflatteninghi0,iflatteninghi1
     & ,numslopes
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndwcomp
      integer idwlo0,idwlo1
      integer idwhi0,idwhi1
      REAL*8 dw(
     & idwlo0:idwhi0,
     & idwlo1:idwhi1,
     & 0:ndwcomp-1)
      integer iflatteninglo0,iflatteninglo1
      integer iflatteninghi0,iflatteninghi1
      REAL*8 flattening(
     & iflatteninglo0:iflatteninghi0,
     & iflatteninglo1:iflatteninghi1)
      integer numslopes
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j,lvar
      do lvar = 0, numslopes - 1
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         dw(i,j,lvar) = flattening(i,j)
     & * dw(i,j,lvar)
      enddo
      enddo
      enddo
      return
      end
      subroutine INCSOURCE(
     & prim
     & ,iprimlo0,iprimlo1
     & ,iprimhi0,iprimhi1
     & ,nprimcomp
     & ,source
     & ,isourcelo0,isourcelo1
     & ,isourcehi0,isourcehi1
     & ,nsourcecomp
     & ,scale
     & ,idcalclo0,idcalclo1
     & ,idcalchi0,idcalchi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nprimcomp
      integer iprimlo0,iprimlo1
      integer iprimhi0,iprimhi1
      REAL*8 prim(
     & iprimlo0:iprimhi0,
     & iprimlo1:iprimhi1,
     & 0:nprimcomp-1)
      integer nsourcecomp
      integer isourcelo0,isourcelo1
      integer isourcehi0,isourcehi1
      REAL*8 source(
     & isourcelo0:isourcehi0,
     & isourcelo1:isourcehi1,
     & 0:nsourcecomp-1)
      REAL*8 scale
      integer idcalclo0,idcalclo1
      integer idcalchi0,idcalchi1
      integer i,j, iv
      REAL*8 increment
      do iv = 0,nprimcomp - 1
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         increment = scale*source(i,j, iv)
         prim(i,j, iv) = prim(i,j, iv) + increment
      enddo
      enddo
      enddo
      return
      end
      subroutine VLLIMITER(
     & slopeprim
     & ,islopeprimlo0,islopeprimlo1
     & ,islopeprimhi0,islopeprimhi1
     & ,nslopeprimcomp
     & ,slopeleft
     & ,islopeleftlo0,islopeleftlo1
     & ,islopelefthi0,islopelefthi1
     & ,nslopeleftcomp
     & ,sloperigh
     & ,isloperighlo0,isloperighlo1
     & ,isloperighhi0,isloperighhi1
     & ,nsloperighcomp
     & ,idcalclo0,idcalclo1
     & ,idcalchi0,idcalchi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nslopeprimcomp
      integer islopeprimlo0,islopeprimlo1
      integer islopeprimhi0,islopeprimhi1
      REAL*8 slopeprim(
     & islopeprimlo0:islopeprimhi0,
     & islopeprimlo1:islopeprimhi1,
     & 0:nslopeprimcomp-1)
      integer nslopeleftcomp
      integer islopeleftlo0,islopeleftlo1
      integer islopelefthi0,islopelefthi1
      REAL*8 slopeleft(
     & islopeleftlo0:islopelefthi0,
     & islopeleftlo1:islopelefthi1,
     & 0:nslopeleftcomp-1)
      integer nsloperighcomp
      integer isloperighlo0,isloperighlo1
      integer isloperighhi0,isloperighhi1
      REAL*8 sloperigh(
     & isloperighlo0:isloperighhi0,
     & isloperighlo1:isloperighhi1,
     & 0:nsloperighcomp-1)
      integer idcalclo0,idcalclo1
      integer idcalchi0,idcalchi1
      integer i,j, iv
      REAL*8 dql, dqr, dqlim
      do iv = 0,nslopeprimcomp - 1
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         dql = slopeleft(i,j, iv)
         dqr = sloperigh(i,j, iv)
         dqlim = slopeprim(i,j, iv)
         call pointvllimiter(dqlim, dql, dqr)
         slopeprim(i,j,iv) = dqlim
      enddo
      enddo
      enddo
      return
      end
      subroutine POINTVLLIMITER(
     & dqlim
     & ,dql
     & ,dqr
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 dqlim
      REAL*8 dql
      REAL*8 dqr
      REAL*8 dqc
      dqc = dqlim
      dqlim = min((2.0d0)*abs(dql),(2.0d0)*abs(dqr))
      dqlim = min(dqlim, abs(dqc))
      if (dql*dqr .lt. (0.0d0)) then
         dqlim = (0.0d0)
      else
         dqlim = dqlim*sign((1.0d0), dql)
         ch_flops=ch_flops+2
      endif
      ch_flops=ch_flops+8
      return
      end
        subroutine DIVUEDGE(
     & divu
     & ,idivulo0,idivulo1
     & ,idivuhi0,idivuhi1
     & ,facedir
     & ,iloboxlo0,iloboxlo1
     & ,iloboxhi0,iloboxhi1
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1
     & ,ihiboxhi0,ihiboxhi1
     & ,hashi
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivulo0,idivulo1
      integer idivuhi0,idivuhi1
      REAL*8 divu(
     & idivulo0:idivuhi0,
     & idivulo1:idivuhi1)
      integer facedir
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer haslo
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer hashi
      integer i,j,ioff,joff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      if (haslo .eq. 1) then
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         divu(i,j) = divu(i+ioff,j+joff)
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         divu(i,j) = divu(i-ioff,j-joff)
      enddo
      enddo
      endif
      return
      end
        subroutine AVEFLUXTOFACE(
     & faceflux
     & ,ifacefluxlo0,ifacefluxlo1
     & ,ifacefluxhi0,ifacefluxhi1
     & ,ccflux
     & ,iccfluxlo0,iccfluxlo1
     & ,iccfluxhi0,iccfluxhi1
     & ,facedir
     & ,ifaceboxlo0,ifaceboxlo1
     & ,ifaceboxhi0,ifaceboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacefluxlo0,ifacefluxlo1
      integer ifacefluxhi0,ifacefluxhi1
      REAL*8 faceflux(
     & ifacefluxlo0:ifacefluxhi0,
     & ifacefluxlo1:ifacefluxhi1)
      integer iccfluxlo0,iccfluxlo1
      integer iccfluxhi0,iccfluxhi1
      REAL*8 ccflux(
     & iccfluxlo0:iccfluxhi0,
     & iccfluxlo1:iccfluxhi1)
      integer facedir
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      integer i,j,ioff,joff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      faceflux(i,j) = (0.500d0)*(
     & ccflux(i ,j ) +
     & ccflux(i - ioff,j - joff))
      enddo
      enddo
      return
      end
        subroutine DIVUONED(
     & divu
     & ,idivulo0,idivulo1
     & ,idivuhi0,idivuhi1
     & ,velnorm
     & ,ivelnormlo0,ivelnormlo1
     & ,ivelnormhi0,ivelnormhi1
     & ,facedir
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivulo0,idivulo1
      integer idivuhi0,idivuhi1
      REAL*8 divu(
     & idivulo0:idivuhi0,
     & idivulo1:idivuhi1)
      integer ivelnormlo0,ivelnormlo1
      integer ivelnormhi0,ivelnormhi1
      REAL*8 velnorm(
     & ivelnormlo0:ivelnormhi0,
     & ivelnormlo1:ivelnormhi1)
      integer facedir
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer i,j,ioff,joff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      divu(i,j) =
     & velnorm(i ,j ) -
     & velnorm(i - ioff,j - joff)
      enddo
      enddo
      return
      end
        subroutine DIVUTRAN(
     & divu
     & ,idivulo0,idivulo1
     & ,idivuhi0,idivuhi1
     & ,slopevel
     & ,islopevello0,islopevello1
     & ,islopevelhi0,islopevelhi1
     & ,facedir
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivulo0,idivulo1
      integer idivuhi0,idivuhi1
      REAL*8 divu(
     & idivulo0:idivuhi0,
     & idivulo1:idivuhi1)
      integer islopevello0,islopevello1
      integer islopevelhi0,islopevelhi1
      REAL*8 slopevel(
     & islopevello0:islopevelhi0,
     & islopevello1:islopevelhi1)
      integer facedir
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      integer i,j,ioff,joff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      divu(i,j) = divu(i,j)+ (0.500d0)*(
     & slopevel(i ,j ) +
     & slopevel(i-ioff,j-joff))
      enddo
      enddo
      return
      end
        subroutine ARTVISC(
     & f
     & ,iflo0,iflo1
     & ,ifhi0,ifhi1
     & ,nfcomp
     & ,u
     & ,iulo0,iulo1
     & ,iuhi0,iuhi1
     & ,nucomp
     & ,divu
     & ,idivulo0,idivulo1
     & ,idivuhi0,idivuhi1
     & ,coeff
     & ,idir
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,numcons
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfcomp
      integer iflo0,iflo1
      integer ifhi0,ifhi1
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1,
     & 0:nfcomp-1)
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL*8 u(
     & iulo0:iuhi0,
     & iulo1:iuhi1,
     & 0:nucomp-1)
      integer idivulo0,idivulo1
      integer idivuhi0,idivuhi1
      REAL*8 divu(
     & idivulo0:idivuhi0,
     & idivulo1:idivuhi1)
      REAL*8 coeff
      integer idir
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer numcons
      REAL*8 dx
        integer i , j
        integer ioff, joff
        integer iv
        REAL*8 fc,dv,s1,s2
        ioff = chf_id(0,idir)
        joff = chf_id(1,idir)
        do iv = 0,numcons - 1
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            fc = f (i ,j ,iv)
            dv = divu(i ,j )
            s1 = u (i ,j ,iv)
            s2 = u (i-ioff,j-joff,iv)
            f(i,j,iv) = fc - coeff*max(-dv, 0.d0)*(s1-s2)
      enddo
      enddo
        enddo
        return
        end
      subroutine UPDATE(
     & idcalclo0,idcalclo1
     & ,idcalchi0,idcalchi1
     & ,state
     & ,istatelo0,istatelo1
     & ,istatehi0,istatehi1
     & ,nstatecomp
     & ,flux
     & ,ifluxlo0,ifluxlo1
     & ,ifluxhi0,ifluxhi1
     & ,nfluxcomp
     & ,facedir
     & ,nconserved
     & ,dtbydx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idcalclo0,idcalclo1
      integer idcalchi0,idcalchi1
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL*8 state(
     & istatelo0:istatehi0,
     & istatelo1:istatehi1,
     & 0:nstatecomp-1)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1
      integer ifluxhi0,ifluxhi1
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & 0:nfluxcomp-1)
      integer facedir
      integer nconserved
      REAL*8 dtbydx
      integer i, j
      integer ioff, joff
      integer spacedim,iv
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      spacedim = 2
      do iv = 0,nconserved - 1
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         state(i,j,iv) = state(i,j,iv) -
     & dtbydx *
     & ( flux(i+ioff,j+joff,iv)
     & - flux(i ,j ,iv))
      enddo
      enddo
      enddo
      return
      end
      subroutine DIVERGEF(
     & idcalclo0,idcalclo1
     & ,idcalchi0,idcalchi1
     & ,divf
     & ,idivflo0,idivflo1
     & ,idivfhi0,idivfhi1
     & ,ndivfcomp
     & ,flux
     & ,ifluxlo0,ifluxlo1
     & ,ifluxhi0,ifluxhi1
     & ,nfluxcomp
     & ,facedir
     & ,nconserved
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idcalclo0,idcalclo1
      integer idcalchi0,idcalchi1
      integer ndivfcomp
      integer idivflo0,idivflo1
      integer idivfhi0,idivfhi1
      REAL*8 divf(
     & idivflo0:idivfhi0,
     & idivflo1:idivfhi1,
     & 0:ndivfcomp-1)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1
      integer ifluxhi0,ifluxhi1
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & 0:nfluxcomp-1)
      integer facedir
      integer nconserved
      REAL*8 dx
      integer i, j
      integer ioff, joff
      integer spacedim,iv
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      spacedim = 2
      do iv = 0,nconserved - 1
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         divf(i,j,iv) = divf(i,j,iv) +
     & (flux(i+ioff,j+joff,iv)
     & -flux(i ,j ,iv))/dx
      enddo
      enddo
      enddo
      return
      end
      subroutine REGUPDATE(
     & idcalclo0,idcalclo1
     & ,idcalchi0,idcalchi1
     & ,consstate
     & ,iconsstatelo0,iconsstatelo1
     & ,iconsstatehi0,iconsstatehi1
     & ,nconsstatecomp
     & ,divf
     & ,idivflo0,idivflo1
     & ,idivfhi0,idivfhi1
     & ,ndivfcomp
     & ,nconserved
     & ,dt
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idcalclo0,idcalclo1
      integer idcalchi0,idcalchi1
      integer nconsstatecomp
      integer iconsstatelo0,iconsstatelo1
      integer iconsstatehi0,iconsstatehi1
      REAL*8 consstate(
     & iconsstatelo0:iconsstatehi0,
     & iconsstatelo1:iconsstatehi1,
     & 0:nconsstatecomp-1)
      integer ndivfcomp
      integer idivflo0,idivflo1
      integer idivfhi0,idivfhi1
      REAL*8 divf(
     & idivflo0:idivfhi0,
     & idivflo1:idivfhi1,
     & 0:ndivfcomp-1)
      integer nconserved
      REAL*8 dt
      integer i, j
      integer iv
      do iv = 0,nconserved - 1
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         consstate(i,j,iv) = consstate(i,j,iv)
     & - dt*divf(i,j,iv)
      enddo
      enddo
      enddo
      return
      end
