      subroutine REGCORRECTTANVEL(
     & vel
     & ,ivello0,ivello1
     & ,ivelhi0,ivelhi1
     & ,grad
     & ,igradlo0,igradlo1
     & ,igradhi0,igradhi1
     & ,iinteriorboxlo0,iinteriorboxlo1
     & ,iinteriorboxhi0,iinteriorboxhi1
     & ,veldir
     & ,graddir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1)
      integer igradlo0,igradlo1
      integer igradhi0,igradhi1
      REAL*8 grad(
     & igradlo0:igradhi0,
     & igradlo1:igradhi1)
      integer iinteriorboxlo0,iinteriorboxlo1
      integer iinteriorboxhi0,iinteriorboxhi1
      integer veldir
      integer graddir
      integer ioffvel , joffvel
      integer ioffgrad, joffgrad
      integer i,j
      REAL*8 correction, factor
      ioffvel = chf_id(0,veldir)
      joffvel = chf_id(1,veldir)
      ioffgrad = chf_id(0,graddir)
      joffgrad = chf_id(1,graddir)
      factor = (1.0d0)/(4.0d0)
      do j = iinteriorboxlo1,iinteriorboxhi1
      do i = iinteriorboxlo0,iinteriorboxhi0
      correction = factor*
     $ ( grad(i ,j )
     $ + grad(i+ioffgrad-ioffvel,j+joffgrad-joffvel)
     $ + grad(i -ioffvel,j -joffvel)
     $ + grad(i+ioffgrad ,j+joffgrad ))
      vel(i,j) = vel(i,j) - correction
      enddo
      enddo
      ch_flops=ch_flops+(iinteriorboxhi0- iinteriorboxlo0+1)*(iinteriorb
     &oxhi1- iinteriorboxlo1+1)*5+1
      return
      end
      subroutine REGCORRECTTANVELVARDENS(
     & vel
     & ,ivello0,ivello1
     & ,ivelhi0,ivelhi1
     & ,grad
     & ,igradlo0,igradlo1
     & ,igradhi0,igradhi1
     & ,oneoverdens
     & ,ioneoverdenslo0,ioneoverdenslo1
     & ,ioneoverdenshi0,ioneoverdenshi1
     & ,iinteriorboxlo0,iinteriorboxlo1
     & ,iinteriorboxhi0,iinteriorboxhi1
     & ,veldir
     & ,graddir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1)
      integer igradlo0,igradlo1
      integer igradhi0,igradhi1
      REAL*8 grad(
     & igradlo0:igradhi0,
     & igradlo1:igradhi1)
      integer ioneoverdenslo0,ioneoverdenslo1
      integer ioneoverdenshi0,ioneoverdenshi1
      REAL*8 oneoverdens(
     & ioneoverdenslo0:ioneoverdenshi0,
     & ioneoverdenslo1:ioneoverdenshi1)
      integer iinteriorboxlo0,iinteriorboxlo1
      integer iinteriorboxhi0,iinteriorboxhi1
      integer veldir
      integer graddir
      integer ioffvel , joffvel
      integer ioffgrad, joffgrad
      integer i,j
      REAL*8 correction, factor
      ioffvel = chf_id(0,veldir)
      joffvel = chf_id(1,veldir)
      ioffgrad = chf_id(0,graddir)
      joffgrad = chf_id(1,graddir)
      factor = (1.0d0)/(4.0d0)
      do j = iinteriorboxlo1,iinteriorboxhi1
      do i = iinteriorboxlo0,iinteriorboxhi0
      correction = factor*
     $ ( grad(i ,j )
     $ + grad(i+ioffgrad-ioffvel,j+joffgrad-joffvel)
     $ + grad(i -ioffvel,j -joffvel)
     $ + grad(i+ioffgrad ,j+joffgrad ))
      vel(i,j) = vel(i,j) - correction*oneoverdens(i,j)
      enddo
      enddo
      ch_flops=ch_flops+(iinteriorboxhi0- iinteriorboxlo0+1)*(iinteriorb
     &oxhi1- iinteriorboxlo1+1)*5+1
      return
      end
      subroutine MACDIVERGEF(
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
      ch_flops=ch_flops+(idcalchi0- idcalclo0+1)*(idcalchi1- idcalclo1+1
     &)*3
      return
      end
      subroutine MACGRADPHI(
     & gradphi
     & ,igradphilo0,igradphilo1
     & ,igradphihi0,igradphihi1
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,facedir
     & ,dx
     & ,ifaceboxlo0,ifaceboxlo1
     & ,ifaceboxhi0,ifaceboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer igradphilo0,igradphilo1
      integer igradphihi0,igradphihi1
      REAL*8 gradphi(
     & igradphilo0:igradphihi0,
     & igradphilo1:igradphihi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer facedir
      REAL*8 dx
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      integer i,j
      integer ioff,joff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      gradphi(i,j) =
     & ( phi(i ,j )
     & - phi(i-ioff,j-joff)
     & )/dx
      enddo
      enddo
      ch_flops=ch_flops+(ifaceboxhi0- ifaceboxlo0+1)*(ifaceboxhi1- iface
     &boxlo1+1)*2
      return
      end
