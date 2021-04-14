        subroutine REGCONSTANTINTERP(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,ivalidboxlo0,ivalidboxlo1
     & ,ivalidboxhi0,ivalidboxhi1
     & ,dir
     & ,nghost
     & ,side
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
      integer ivalidboxlo0,ivalidboxlo1
      integer ivalidboxhi0,ivalidboxhi1
      integer dir
      integer nghost
      integer side
        integer ioff,joff
        integer ighost, jghost
        integer ivalid, jvalid
        integer whichghost
        if((side.ne.1).and.(side.ne.-1)) then
           call MAYDAY_ERROR()
        endif
        ioff = chf_id(0,dir)*side
        joff = chf_id(1,dir)*side
      do jvalid = ivalidboxlo1,ivalidboxhi1
      do ivalid = ivalidboxlo0,ivalidboxhi0
        do whichghost = 1, nghost
           ighost = ivalid + whichghost*ioff
           jghost = jvalid + whichghost*joff
           phi(ighost,jghost) = phi(ivalid,jvalid)
        enddo
      enddo
      enddo
        return
        end
