#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine TRAPWEIGHTS(
     &           wt
     &           ,iwtlo0,iwtlo1
     &           ,iwthi0,iwthi1
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dv
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iwtlo0,iwtlo1
      integer iwthi0,iwthi1
      REAL_T wt(
     &           iwtlo0:iwthi0,
     &           iwtlo1:iwthi1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T dv
      integer i0,i1
      
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

         wt(i0,i1) = dv
         
         if ((i0 .eq. iregionlo0) .or.
     &        (i0 .eq. iregionhi0)) then
             wt(i0,i1) = wt(i0,i1) * half
         endif 
         if ((i1 .eq. iregionlo1) .or.
     &        (i1 .eq. iregionhi1)) then
             wt(i0,i1) = wt(i0,i1) * half
         endif 
      
      enddo
      enddo
      return
      end
      subroutine NODEDOTPRODUCT(
     &           dotprodout
     &           ,afab
     &           ,iafablo0,iafablo1
     &           ,iafabhi0,iafabhi1
     &           ,nafabcomp
     &           ,bfab
     &           ,ibfablo0,ibfablo1
     &           ,ibfabhi0,ibfabhi1
     &           ,nbfabcomp
     &           ,wfab
     &           ,iwfablo0,iwfablo1
     &           ,iwfabhi0,iwfabhi1
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T dotprodout
      integer nafabcomp
      integer iafablo0,iafablo1
      integer iafabhi0,iafabhi1
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           0:nafabcomp-1)
      integer nbfabcomp
      integer ibfablo0,ibfablo1
      integer ibfabhi0,ibfabhi1
      REAL_T bfab(
     &           ibfablo0:ibfabhi0,
     &           ibfablo1:ibfabhi1,
     &           0:nbfabcomp-1)
      integer iwfablo0,iwfablo1
      integer iwfabhi0,iwfabhi1
      REAL_T wfab(
     &           iwfablo0:iwfabhi0,
     &           iwfablo1:iwfabhi1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer startcomp
      integer endcomp
      integer i0,i1
      integer nv
      dotprodout = 0
      do nv=startcomp,endcomp,1
         
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

            dotprodout = dotprodout +
     &        afab(i0,i1,nv)*
     &        bfab(i0,i1,nv)*
     &        wfab(i0,i1)
         
      enddo
      enddo
      enddo
      return
      end
      subroutine NODENORM(
     &           normout
     &           ,p
     &           ,afab
     &           ,iafablo0,iafablo1
     &           ,iafabhi0,iafabhi1
     &           ,nafabcomp
     &           ,wfab
     &           ,iwfablo0,iwfablo1
     &           ,iwfabhi0,iwfabhi1
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T normout
      integer p
      integer nafabcomp
      integer iafablo0,iafablo1
      integer iafabhi0,iafabhi1
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           0:nafabcomp-1)
      integer iwfablo0,iwfablo1
      integer iwfabhi0,iwfabhi1
      REAL_T wfab(
     &           iwfablo0:iwfabhi0,
     &           iwfablo1:iwfabhi1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer startcomp
      integer endcomp
      integer i0,i1
      integer nv
      REAL_T pwrinv
      normout = 0
      if (p .eq. 1) then
         do nv = startcomp, endcomp, 1
            
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

               normout = normout +
     &           wfab(i0,i1) *
     &           abs(afab(i0,i1, nv))
            
      enddo
      enddo
         enddo
      elseif (p .eq. 2) then
         do nv = startcomp, endcomp, 1
            
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

            normout = normout +
     &           wfab(i0,i1) *
     &           afab(i0,i1, nv) *
     &           afab(i0,i1, nv)
            
      enddo
      enddo
         enddo
         normout = sqrt(normout)
      else
         do nv = startcomp, endcomp, 1
            
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

               normout = normout +
     &           wfab(i0,i1) *
     &           (afab(i0,i1, nv)**p)
            
      enddo
      enddo
         enddo
         pwrinv = one / p
         normout = normout**pwrinv
      endif
      return
      end
      subroutine NODEINTEGRAL(
     &           ans
     &           ,afab
     &           ,iafablo0,iafablo1
     &           ,iafabhi0,iafabhi1
     &           ,nafabcomp
     &           ,wfab
     &           ,iwfablo0,iwfablo1
     &           ,iwfabhi0,iwfabhi1
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T ans
      integer nafabcomp
      integer iafablo0,iafablo1
      integer iafabhi0,iafabhi1
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           0:nafabcomp-1)
      integer iwfablo0,iwfablo1
      integer iwfabhi0,iwfabhi1
      REAL_T wfab(
     &           iwfablo0:iwfabhi0,
     &           iwfablo1:iwfabhi1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer startcomp
      integer endcomp
      integer i0,i1
      integer nv
      ans = 0
      do nv = startcomp, endcomp, 1
         
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

            ans = ans +
     &        wfab(i0,i1) * afab(i0,i1, nv)
         
      enddo
      enddo
      enddo
      return
      end
      subroutine NODEMAXNORM(
     &           normout
     &           ,afab
     &           ,iafablo0,iafablo1
     &           ,iafabhi0,iafabhi1
     &           ,nafabcomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T normout
      integer nafabcomp
      integer iafablo0,iafablo1
      integer iafabhi0,iafabhi1
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           0:nafabcomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer startcomp
      integer endcomp
      integer i0,i1
      integer nv
      REAL_T this
      normout = 0
      do nv = startcomp, endcomp, 1
         
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

            this = abs(afab(i0,i1, nv))
            if (this .gt. normout) then
               normout = this
            endif
         
      enddo
      enddo
      enddo
      return
      end
