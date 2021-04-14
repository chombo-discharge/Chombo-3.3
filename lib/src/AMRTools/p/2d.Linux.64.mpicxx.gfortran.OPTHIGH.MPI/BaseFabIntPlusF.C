#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine BASEFABINTPLUS(
     &           sum
     &           ,isumlo0,isumlo1
     &           ,isumhi0,isumhi1
     &           ,nsumcomp
     &           ,piece
     &           ,ipiecelo0,ipiecelo1
     &           ,ipiecehi0,ipiecehi1
     &           ,npiececomp
     &           ,ibxlo0,ibxlo1
     &           ,ibxhi0,ibxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nsumcomp
      integer isumlo0,isumlo1
      integer isumhi0,isumhi1
      integer sum(
     &           isumlo0:isumhi0,
     &           isumlo1:isumhi1,
     &           0:nsumcomp-1)
      integer npiececomp
      integer ipiecelo0,ipiecelo1
      integer ipiecehi0,ipiecehi1
      integer piece(
     &           ipiecelo0:ipiecehi0,
     &           ipiecelo1:ipiecehi1,
     &           0:npiececomp-1)
      integer ibxlo0,ibxlo1
      integer ibxhi0,ibxhi1
      integer i0,i1
      integer var
      do var = 0, nsumcomp-1
         
      do i1 = ibxlo1,ibxhi1
      do i0 = ibxlo0,ibxhi0

            sum(i0,i1, var) = sum(i0,i1, var) +
     &        piece(i0,i1, var)
         
      enddo
      enddo
      enddo
      return
      end
