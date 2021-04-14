#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine INTERPCONSTANT(
     &           fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,nfinecomp
     &           ,coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,ncoarsecomp
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,ref_ratio
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           0:nfinecomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           0:ncoarsecomp-1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer ref_ratio
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer var
      integer ic0,ic1
      integer if0,if1
      integer ii0,ii1
      do var = 0, ncoarsecomp - 1
         
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

            
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               if1 = ic1*ref_ratio + ii1
               fine(if0,if1,var) = coarse(ic0,ic1,var)
            
      enddo
      enddo
         
      enddo
      enddo
      end do
      return
      end
      subroutine INTERPCENTRALSLOPE(
     &           slope
     &           ,islopelo0,islopelo1
     &           ,islopehi0,islopehi1
     &           ,nslopecomp
     &           ,state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,nstatecomp
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopecomp
      integer islopelo0,islopelo1
      integer islopehi0,islopehi1
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1,
     &           0:nslopecomp-1)
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           0:nstatecomp-1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer dir
      integer i0,i1
      integer ii0,ii1
      integer var
      
      ii0=CHF_ID(0, dir)

      ii1=CHF_ID(1, dir)

      do var = 0, nstatecomp - 1
         
      do i1 = iblo1,ibhi1
      do i0 = iblo0,ibhi0

          slope (i0,i1,var) = half * (
     &        state (i0+ii0,i1+ii1,var) -
     &        state (i0-ii0,i1-ii1,var) )
          
      enddo
      enddo
       end do
      return
      end
      subroutine INTERPHISIDESLOPE(
     &           slope
     &           ,islopelo0,islopelo1
     &           ,islopehi0,islopehi1
     &           ,nslopecomp
     &           ,state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,nstatecomp
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopecomp
      integer islopelo0,islopelo1
      integer islopehi0,islopehi1
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1,
     &           0:nslopecomp-1)
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           0:nstatecomp-1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer dir
      integer i0,i1
      integer ii0,ii1
      integer var
      
      ii0=CHF_ID(0, dir)

      ii1=CHF_ID(1, dir)

      do var = 0, nstatecomp - 1
         
      do i1 = iblo1,ibhi1
      do i0 = iblo0,ibhi0

          slope (i0,i1,var) =
     &          state ( i0+ii0,i1+ii1, var)
     &        - state ( i0,i1, var)
          
      enddo
      enddo
       enddo
      return
      end
      subroutine INTERPLOSIDESLOPE(
     &           slope
     &           ,islopelo0,islopelo1
     &           ,islopehi0,islopehi1
     &           ,nslopecomp
     &           ,state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,nstatecomp
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopecomp
      integer islopelo0,islopelo1
      integer islopehi0,islopehi1
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1,
     &           0:nslopecomp-1)
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           0:nstatecomp-1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer dir
      integer i0,i1, var
      integer ii0,ii1
      
      ii0=CHF_ID(0, dir)

      ii1=CHF_ID(1, dir)

      do var = 0, nstatecomp - 1
         
      do i1 = iblo1,ibhi1
      do i0 = iblo0,ibhi0

         slope (i0,i1,var) =
     &        state (  i 0, i 1, var) -
     &        state (  i0-ii0, i1-ii1, var)
          
      enddo
      enddo
       end do
      return
      end
      subroutine INTERPLIMIT(
     &           islope
     &           ,iislopelo0,iislopelo1
     &           ,iislopehi0,iislopehi1
     &           ,nislopecomp
     &           ,jslope
     &           ,ijslopelo0,ijslopelo1
     &           ,ijslopehi0,ijslopehi1
     &           ,njslopecomp
     &           ,kslope
     &           ,ikslopelo0,ikslopelo1
     &           ,ikslopehi0,ikslopehi1
     &           ,nkslopecomp
     &           ,state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,nstatecomp
     &           ,ibcoarselo0,ibcoarselo1
     &           ,ibcoarsehi0,ibcoarsehi1
     &           ,ibnlo0,ibnlo1
     &           ,ibnhi0,ibnhi1
     &           ,iphysdomainlo0,iphysdomainlo1
     &           ,iphysdomainhi0,iphysdomainhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nislopecomp
      integer iislopelo0,iislopelo1
      integer iislopehi0,iislopehi1
      REAL_T islope(
     &           iislopelo0:iislopehi0,
     &           iislopelo1:iislopehi1,
     &           0:nislopecomp-1)
      integer njslopecomp
      integer ijslopelo0,ijslopelo1
      integer ijslopehi0,ijslopehi1
      REAL_T jslope(
     &           ijslopelo0:ijslopehi0,
     &           ijslopelo1:ijslopehi1,
     &           0:njslopecomp-1)
      integer nkslopecomp
      integer ikslopelo0,ikslopelo1
      integer ikslopehi0,ikslopehi1
      REAL_T kslope(
     &           ikslopelo0:ikslopehi0,
     &           ikslopelo1:ikslopehi1,
     &           0:nkslopecomp-1)
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           0:nstatecomp-1)
      integer ibcoarselo0,ibcoarselo1
      integer ibcoarsehi0,ibcoarsehi1
      integer ibnlo0,ibnlo1
      integer ibnhi0,ibnhi1
      integer iphysdomainlo0,iphysdomainlo1
      integer iphysdomainhi0,iphysdomainhi1
      integer  i0, i1, var
      integer  ii0, ii1
      integer  in0, in1
      REAL_T statemax, statemin, deltasum,  eta
      do var = 0, nislopecomp - 1
         
      do i1 = ibcoarselo1,ibcoarsehi1
      do i0 = ibcoarselo0,ibcoarsehi0

             statemax = state ( i0,i1, var )
             statemin = state ( i0,i1, var )
             
      do ii1 = ibnlo1,ibnhi1
      do ii0 = ibnlo0,ibnhi0

             
                 in0 = i0 + ii0
                 in1 = i1 + ii1
                 if (
                 
     &                in0 .ge. istatelo0 .and.
     &                in0 .le. istatehi0 
     &                .and.
     &                in1 .ge. istatelo1 .and.
     &                in1 .le. istatehi1 
     &                ) then
                    statemax = max ( statemax, state(in0,in1,var))
                    statemin = min ( statemin, state(in0,in1,var))
                 endif
             
      enddo
      enddo
             deltasum = half * (
                
     &            abs ( islope ( i0,i1, var ) )
     &            +
     &            abs ( jslope ( i0,i1, var ) )
     &            )
#if CH_SPACEDIM > 3
                call MAYDAY_ERROR()
#endif
              eta = min(statemax - state(i0,i1,var),
     &                  state(i0,i1,var) - statemin)
              if( eta .le. 1.e-9*abs(statemax) ) then
                 eta = zero
              else
              if (deltasum .gt. eta) then
                eta = eta/deltasum
              else
                eta = one
              endif
              endif
              
              islope ( i0,i1, var ) =
     &             eta * islope ( i0,i1, var ) 
              jslope ( i0,i1, var ) =
     &             eta * jslope ( i0,i1, var ) 
#if CH_SPACEDIM > 3
              call MAYDAY_ERROR()
#endif
         
      enddo
      enddo
      end do
      return
      end
      subroutine INTERPLINEAR(
     &           fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,nfinecomp
     &           ,slope
     &           ,islopelo0,islopelo1
     &           ,islopehi0,islopehi1
     &           ,nslopecomp
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,dir
     &           ,ref_ratio
     &           ,ibreflo0,ibreflo1
     &           ,ibrefhi0,ibrefhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           0:nfinecomp-1)
      integer nslopecomp
      integer islopelo0,islopelo1
      integer islopehi0,islopehi1
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1,
     &           0:nslopecomp-1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer dir
      integer ref_ratio
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer ic0,ic1
      integer if0,if1
      integer ii0,ii1
      integer var, id
      REAL_T dxf
      do var = 0, nfinecomp - 1
          
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

              
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  if1 = ic1*ref_ratio + ii1
              
                  if (dir .eq. 0) then
                      id = ii0
                  else if (dir .eq. 1) then
                      id = ii1
                  endif
                  dxf = -half + ( (id+half) / ref_ratio )
                  fine( if0,if1,var) =
     &                 fine( if0,if1,var) +
     &                 dxf * slope (   ic 0,  ic 1, var )
              
      enddo
      enddo
          
      enddo
      enddo
      end do
      return
      end
      subroutine INTERPHOMO_OLD(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,x1
     &           ,dxCrse
     &           ,idir
     &           ,ihilo
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
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T x1
      REAL_T dxCrse
      integer idir
      integer ihilo
      REAL_T x2, denom, idenom, x, xsquared, m1, m2
      REAL_T q1, q2
      REAL_T pa, pb, a, b
      INTEGER ncomp,  n
      INTEGER ii0,ii1
      INTEGER i0,i1
      x2 = half*(three*x1+dxCrse)
      denom = one-((x1+x2)/x1)
      idenom = one/(denom)
      x = two*x1
      xsquared = x*x
      m1 = one/(x1*x1)
      m2 = one/(x1*(x1-x2))
      q1 = one/(x1-x2)
      q2 = x1+x2
      ihilo = ihilo*(-1)
      ncomp = nphicomp
      
      ii0= ihilo*CHF_ID(0, idir)

      ii1= ihilo*CHF_ID(1, idir)

      do n = 0, ncomp-1
          
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

          pa=phi(i0+2*ii0,i1+2*ii1,n)
          pb=phi(i0+ii0,i1+ii1,n)
          a=((pb-pa)*m1 - (pb)*m2)*idenom
          b=(pb)*q1 - a*q2
          phi(i0,i1,n) = a*xsquared + b*x + pa
          
      enddo
      enddo
      enddo
      return
      end
      subroutine INTERPHOMO(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,x1
     &           ,dxCrse
     &           ,idir
     &           ,ihilo
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
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T x1
      REAL_T dxCrse
      integer idir
      integer ihilo
      REAL_T c1, c2
      REAL_T pa, pb
      INTEGER ncomp,  n
      INTEGER ii0,ii1
      INTEGER i0,i1
      c1=two*(dxCrse-x1)/(dxCrse+x1)
      c2=   -(dxCrse-x1)/(dxCrse+three*x1)
      ihilo = ihilo*(-1)
      ncomp = nphicomp
      
      ii0= ihilo*CHF_ID(0, idir)

      ii1= ihilo*CHF_ID(1, idir)

      do n = 0, ncomp-1
          
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

          pa=phi(i0+ii0,i1+ii1,n)
          pb=phi(i0+2*ii0,i1+2*ii1,n)
          phi(i0,i1,n) = c1*pa + c2*pb
          
      enddo
      enddo
      enddo
      return
      end
      subroutine INTERPHOMOLINEAR(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,x1
     &           ,dxCrse
     &           ,idir
     &           ,ihilo
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
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL_T x1
      REAL_T dxCrse
      integer idir
      integer ihilo
      INTEGER ncomp,  n
      INTEGER ii0,ii1
      INTEGER i0,i1
      REAL_T pa, factor
      ihilo = ihilo*(-1)
      ncomp = nphicomp
      
      ii0= ihilo*CHF_ID(0, idir)

      ii1= ihilo*CHF_ID(1, idir)

      factor = one - two*x1/(x1+dxCrse)
      do n = 0, ncomp-1
          
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

          pa=phi(i0+ii0,i1+ii1,n)
          phi(i0,i1,n) = factor*pa
          
      enddo
      enddo
      enddo
      return
      end
