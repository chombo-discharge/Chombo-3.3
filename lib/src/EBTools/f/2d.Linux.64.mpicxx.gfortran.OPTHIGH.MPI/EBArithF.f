      subroutine VOLWGTSUM(
     & src
     & ,isrclo0,isrclo1
     & ,isrchi0,isrchi1
     & ,nsrccomp
     & ,volfrac
     & ,ivolfraclo0,ivolfraclo1
     & ,ivolfrachi0,ivolfrachi1
     & ,nvolfraccomp
     & ,norm
     & ,volume
     & ,comp
     & ,pval
     & ,idoreg
     & ,idoirr
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nsrccomp
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL*8 src(
     & isrclo0:isrchi0,
     & isrclo1:isrchi1,
     & 0:nsrccomp-1)
      integer nvolfraccomp
      integer ivolfraclo0,ivolfraclo1
      integer ivolfrachi0,ivolfrachi1
      REAL*8 volfrac(
     & ivolfraclo0:ivolfrachi0,
     & ivolfraclo1:ivolfrachi1,
     & 0:nvolfraccomp-1)
      REAL*8 norm
      REAL*8 volume
      integer comp
      integer pval
      integer idoreg
      integer idoirr
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer i,j
      integer ncompsrc,counter
      REAL*8 eps, rabs, vfrac
      logical isreg, iscov, isirr, usethis, useirr, usereg
      ncompsrc = nsrccomp
      if(pval .lt. 0) then
         call MAYDAY_ERROR
      endif
      if(ncompsrc .le. comp) then
         call MAYDAY_ERROR()
      endif
      eps = 1.0e-9
      usereg = (idoreg .eq. 1)
      useirr = (idoirr .eq. 1)
      counter=0
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
      vfrac = volfrac(i,j,0)
      iscov = (vfrac .lt. eps)
      isreg = (abs((1.0d0)-vfrac) .lt. eps)
      isirr = (.not. isreg).and.(.not. iscov)
      usethis = ((isreg .and. usereg).or.(isirr .and. useirr))
      if(usethis) then
         counter=counter+1
         if(pval .eq. 0) then
            rabs = abs(src(i,j,comp))
            norm = max(norm, rabs)
         elseif(pval .eq. 1) then
            rabs = abs(src(i,j,comp))
            norm = norm + vfrac*rabs
         elseif(pval .eq. 2) then
            rabs = src(i,j,comp)*src(i,j,comp)
            norm = norm + vfrac*rabs
         else
            call MAYDAY_ERROR()
         endif
         volume = volume + vfrac
      endif
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1) + counter*(pval+1)
      return
      end
      subroutine ADDTWOFAB(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,src
     & ,isrclo0,isrclo1
     & ,isrchi0,isrchi1
     & ,nsrccomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,srccomp
     & ,destcomp
     & ,numcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL*8 src(
     & isrclo0:isrchi0,
     & isrclo1:isrchi1,
     & 0:nsrccomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer srccomp
      integer destcomp
      integer numcomp
      integer i,j
      integer nsrc, ndest,ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     & ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp+numcomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,ndest) =
     & dst(i,j,ndest) +
     & src(i,j,nsrc)
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*numcomp
      return
      end
      subroutine SCALEADDTWOFAB(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,src
     & ,isrclo0,isrclo1
     & ,isrchi0,isrchi1
     & ,nsrccomp
     & ,scale
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,srccomp
     & ,destcomp
     & ,numcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL*8 src(
     & isrclo0:isrchi0,
     & isrclo1:isrchi1,
     & 0:nsrccomp-1)
      REAL*8 scale
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer srccomp
      integer destcomp
      integer numcomp
      integer i,j
      integer nsrc, ndest, ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     & ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,ndest) =
     & dst(i,j,ndest) +
     & src(i,j,nsrc) * scale
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*2*numcomp
      return
      end
      subroutine AXBYFAB(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,nxcomp
     & ,y
     & ,iylo0,iylo1
     & ,iyhi0,iyhi1
     & ,nycomp
     & ,a
     & ,b
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,srccomp
     & ,destcomp
     & ,numcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1,
     & 0:nxcomp-1)
      integer nycomp
      integer iylo0,iylo1
      integer iyhi0,iyhi1
      REAL*8 y(
     & iylo0:iyhi0,
     & iylo1:iyhi1,
     & 0:nycomp-1)
      REAL*8 a
      REAL*8 b
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer srccomp
      integer destcomp
      integer numcomp
      integer i,j
      integer nsrc, ndest, ncompsrc, ncompdst
      ncompsrc = nxcomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     & ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp-1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,ndest) =
     & x(i,j,nsrc) * a
     & + y(i,j,nsrc) * b
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*3*numcomp
      return
      end
      subroutine AXBYFABCOMP(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,nxcomp
     & ,y
     & ,iylo0,iylo1
     & ,iyhi0,iyhi1
     & ,nycomp
     & ,a
     & ,b
     & ,destcomp
     & ,xcomp
     & ,ycomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1,
     & 0:nxcomp-1)
      integer nycomp
      integer iylo0,iylo1
      integer iyhi0,iyhi1
      REAL*8 y(
     & iylo0:iyhi0,
     & iylo1:iyhi1,
     & 0:nycomp-1)
      REAL*8 a
      REAL*8 b
      integer destcomp
      integer xcomp
      integer ycomp
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer i,j
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
      dst(i,j,destcomp) =
     & x(i,j,xcomp) * a
     & + y(i,j,ycomp) * b
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*3
      return
      end
      subroutine SUBTRACTTWOFAB(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,src
     & ,isrclo0,isrclo1
     & ,isrchi0,isrchi1
     & ,nsrccomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,srccomp
     & ,destcomp
     & ,numcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL*8 src(
     & isrclo0:isrchi0,
     & isrclo1:isrchi1,
     & 0:nsrccomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer srccomp
      integer destcomp
      integer numcomp
      integer i,j
      integer nsrc, ndest,ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     & ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp - 1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,ndest) =
     & dst(i,j,ndest) -
     & src(i,j,nsrc)
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*numcomp
      return
      end
      subroutine MULTIPLYTWOFAB(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,src
     & ,isrclo0,isrclo1
     & ,isrchi0,isrchi1
     & ,nsrccomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,srccomp
     & ,destcomp
     & ,numcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL*8 src(
     & isrclo0:isrchi0,
     & isrclo1:isrchi1,
     & 0:nsrccomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer srccomp
      integer destcomp
      integer numcomp
      integer i,j
      integer nsrc, ndest,ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     & ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp - 1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,ndest) =
     & dst(i,j,ndest)*
     & src(i,j,nsrc)
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*numcomp
      return
      end
      subroutine DIVIDETWOFAB(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,src
     & ,isrclo0,isrclo1
     & ,isrchi0,isrchi1
     & ,nsrccomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,srccomp
     & ,destcomp
     & ,numcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL*8 src(
     & isrclo0:isrchi0,
     & isrclo1:isrchi1,
     & 0:nsrccomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer srccomp
      integer destcomp
      integer numcomp
      integer i,j
      integer nsrc, ndest,ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     & ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp - 1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,ndest) =
     & dst(i,j,ndest)/
     & src(i,j,nsrc)
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*numcomp
      return
      end
      subroutine SUBTRACTFABR(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,src
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      REAL*8 src
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer i,j
      integer n, ncompdst
      ncompdst = ndstcomp
      do n = 0, ncompdst - 1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,n) =
     & dst(i,j,n) - src
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*ncompdst
      return
      end
      subroutine ADDFABR(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,src
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      REAL*8 src
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer i,j
      integer n, ncompdst
      ncompdst = ndstcomp
      do n = 0, ncompdst - 1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,n) =
     & dst(i,j,n) + src
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*ncompdst
      return
      end
      subroutine MULTIPLYFABR(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,src
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      REAL*8 src
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer i,j
      integer n,ncompdst
      ncompdst = ndstcomp
      do n = 0, ncompdst - 1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,n) =
     & dst(i,j,n) * src
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*ncompdst
      return
      end
      subroutine DIVIDEFABR(
     & dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,ndstcomp
     & ,src
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1,
     & 0:ndstcomp-1)
      REAL*8 src
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer i,j
      integer n,ncompdst
      ncompdst = ndstcomp
      do n = 0, ncompdst - 1
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         dst(i,j,n) =
     & dst(i,j,n)/src
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*ncompdst
      return
      end
      subroutine EBDOTPRODUCT(
     & dotprodout
     & ,afab
     & ,iafablo0,iafablo1
     & ,iafabhi0,iafabhi1
     & ,nafabcomp
     & ,bfab
     & ,ibfablo0,ibfablo1
     & ,ibfabhi0,ibfabhi1
     & ,nbfabcomp
     & ,volfrac
     & ,ivolfraclo0,ivolfraclo1
     & ,ivolfrachi0,ivolfrachi1
     & ,nvolfraccomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,icomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 dotprodout
      integer nafabcomp
      integer iafablo0,iafablo1
      integer iafabhi0,iafabhi1
      REAL*8 afab(
     & iafablo0:iafabhi0,
     & iafablo1:iafabhi1,
     & 0:nafabcomp-1)
      integer nbfabcomp
      integer ibfablo0,ibfablo1
      integer ibfabhi0,ibfabhi1
      REAL*8 bfab(
     & ibfablo0:ibfabhi0,
     & ibfablo1:ibfabhi1,
     & 0:nbfabcomp-1)
      integer nvolfraccomp
      integer ivolfraclo0,ivolfraclo1
      integer ivolfrachi0,ivolfrachi1
      REAL*8 volfrac(
     & ivolfraclo0:ivolfrachi0,
     & ivolfraclo1:ivolfrachi1,
     & 0:nvolfraccomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer icomp
      integer i,j, ncomp, counter
      REAL*8 eps
      if(icomp .lt. 0) then
         call MAYDAY_ERROR()
      endif
      ncomp = nafabcomp
      if(ncomp .le. icomp) then
         call MAYDAY_ERROR()
      endif
      ncomp = nbfabcomp
      if(ncomp .le. icomp) then
         call MAYDAY_ERROR()
      endif
      dotprodout = (0.0d0)
      eps = 1.0e-10
      counter=0
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
      if(volfrac(i,j,0) .gt. eps) then
         counter=counter+1
         dotprodout = dotprodout +
     & afab(i,j,icomp)*bfab(i,j,icomp)
      endif
      enddo
      enddo
      ch_flops=ch_flops+counter*2
      return
      end
      subroutine MAXFAB(
     & aval
     & ,afab
     & ,iafablo0,iafablo1
     & ,iafabhi0,iafabhi1
     & ,nafabcomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,acomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 aval
      integer nafabcomp
      integer iafablo0,iafablo1
      integer iafabhi0,iafabhi1
      REAL*8 afab(
     & iafablo0:iafabhi0,
     & iafablo1:iafabhi1,
     & 0:nafabcomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer acomp
      integer i,j
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
        aval = max(aval, afab(i,j,acomp))
      enddo
      enddo
      return
      end
      subroutine MINFAB(
     & aval
     & ,afab
     & ,iafablo0,iafablo1
     & ,iafabhi0,iafabhi1
     & ,nafabcomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,acomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 aval
      integer nafabcomp
      integer iafablo0,iafablo1
      integer iafabhi0,iafabhi1
      REAL*8 afab(
     & iafablo0:iafabhi0,
     & iafablo1:iafabhi1,
     & 0:nafabcomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer acomp
      integer i,j
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
        aval = min(aval, afab(i,j,acomp))
      enddo
      enddo
      return
      end
