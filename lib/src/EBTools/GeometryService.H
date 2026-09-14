#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#ifndef _GEOMETRYSERVICE_H_
#define _GEOMETRYSERVICE_H_

#include <cmath>
#include <cstdlib>

#include "REAL.H"
#include "LoHiSide.H"
#include "RealVect.H"
#include "Box.H"
#include "IntVect.H"

#include "EBISBox.H"
#include "VolIndex.H"
#include "FaceIndex.H"
#include "IrregNode.H"
#include "BaseIVFAB.H"
#include "DisjointBoxLayout.H"
#include "NamespaceHeader.H"

///
/**
   This is a pure base class that EBIndexSpace uses
   to interface with geometry-generation classes.
 */
class GeometryService
{
public:

  enum InOut
  {
    Regular,
    Covered,
    Irregular
  };

  GeometryService();

  virtual ~GeometryService();

  ///
  /**
      Return true if every cell in region is regular at the
      refinement described by dx.
  */
  virtual bool isRegular(const Box&           a_region,
                         const ProblemDomain& a_domain,
                         const RealVect&      a_origin,
                         const Real&          a_dx) const = 0;

  ///
  /**
      Return true if every cell in region is covered at the
      refinement described by dx.
  */
  virtual bool isCovered(const Box&           a_region,
                         const ProblemDomain& a_domain,
                         const RealVect&      a_origin,
                         const Real&          a_dx) const = 0;

  virtual bool isIrregular(const Box&           a_region,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const ;
  ///
  /**
     irregGraph contains a complete list of irregular vofs. \\
     validRegion is the the region of the graph \\
     domain is the domain of compuation \\
     regIrregCovered = 1 for regular, 0 for irregular, -1 for covered.
     regIrregCovered must contain the valid region grown by 1. \\
     In chombospeak, \\
     ghostregion = (grow(validRegion, 1) & domain);
  */
  /*
  virtual void fillGraph(BaseFab<int>&        a_regIrregCovered,
                         Vector<IrregNode>&   a_nodes,
                         const Box&           a_validRegion,
                         const Box&           a_ghostRegion,
                         const ProblemDomain& a_domain,
                         const RealVect&      a_origin,
                         const Real&          a_dx) const = 0;
  */
  
  ///
  /**
    This overloaded method will be called when the GeometryService
    was asked to make the grids with  makeGrids()
  */
  virtual void fillGraph(BaseFab<int>&        a_regIrregCovered,
                         Vector<IrregNode>&   a_nodes,
                         const Box&           a_validRegion,
                         const Box&           a_ghostRegion,
                         const ProblemDomain& a_domain,
                         const RealVect&      a_origin,
                         const Real&          a_dx,
                         const DataIndex&     a_di ) const = 0;
  // {
  //   return fillGraph( a_regIrregCovered, a_nodes, a_validRegion, a_ghostRegion,
  //                    a_domain, a_origin, a_dx );
  // }


  virtual bool canGenerateMultiCells() const;

  virtual InOut InsideOutside(const Box&           a_region,
                              const ProblemDomain& a_domain,
                              const RealVect&      a_origin,
                              const Real&          a_dx) const ;

      
  ///return true if the geometry service generates higher order moments
  /**
     Most of the geometry generating objects we have do *not* generate higher order
     moments (they only generate the standard stuff (volume fractions, area fractions and centroids)).
     WrappedGShop, does generate the higher order stuff so that one can return true;
  */
  virtual bool generatesHigherOrderMoments() const
    {
      return false; 
    }


  ///
  /**
    This overloaded method will be called when the GeometryService
    was asked to make the grids with  makeGrids()
  */

  virtual InOut InsideOutside(const Box&           a_region,
                              const ProblemDomain& a_domain,
                              const RealVect&      a_origin,
                              const Real&          a_dx,
                              const DataIndex&     a_di ) const
  {
    return InsideOutside( a_region, a_domain, a_origin, a_dx );
  }

  ///handy functions to do rectangle intersections in real space
  static bool intersection(const RealVect& a_lo1, const RealVect& a_hi1,
                           const RealVect& a_lo2, const RealVect& a_hi2);

  static bool intersection(const Box&           a_region,
                           const RealVect&      a_origin,
                           const Real&          a_dx,
                           const RealVect&      a_lower,
                           const RealVect&      a_upper);


  /** Perform any adjustments needed after box layout is created during EBISLevel construction
      This function is used because subclasses may have distributed data that they will
      want to adjust based on the new layout. The default implementation is empty.
  */
  virtual void postMakeBoxLayout(const DisjointBoxLayout& a_dbl,
                                 const RealVect& a_dx);


  ///
  /**
     A GeometryService has three options for implementing this function
     1) do nothing, allow the empty base implementation to remain in place as a null-op
     2) take the makeGrids call as a directive: Here are the grids EBIndexSpace is wanting to use, configure
         yourself accordingly to make this efficient for you.  
     3) discard the DisjointBoxLayout EBIndexSpace would like and insert your own implementation of layout
        EBIndexSpace will faithfully use a_grids returned from this function, including it's load balance. 
  */
  /** How many numbers a geometry service keeps per cut cell so that the cell can be refined later.

      Zero means the service cannot refine a cell it has already made, which is the default and
      what every existing service does.  A service that answers with a positive count is promising
      that getSurfaces will hand those numbers back and that fillRefinedGraph will accept them.
  */
  virtual int numSurfaceComponents() const;

  /** Hand back what was kept for each cut cell of a box, so that the index space can carry it.

      The values are numSurfaceComponents() per cell, in the order of a_cells.  A service that
      keeps nothing returns nothing, which is the default.
      \param a_cells   Cut cells of the box, in the order their values appear.
      \param a_values  numSurfaceComponents() values per entry of a_cells.
      Asked for by region rather than by box index, because the level that asks may have grown
      since the service made it and its boxes are then no longer the service's own.
      \param a_region  Cells wanted.
      \param a_dx      Grid spacing of the level the cells belong to.
  */
  virtual void getSurfaces(Vector<IntVect>& a_cells,
                           Vector<Real>&    a_values,
                           const Box&       a_region,
                           const Real&      a_dx) const;

  /** Spacing of the coarser level a region has to be cut from, or a negative number if the
      service cannot cut it.

      Which level a region is cut from is a property of the geometry and not of the rank asking,
      since two cells cut from different ancestors need not agree on the face between them.
      \param a_ghostRegion Region wanted, including its ghost cells, on its own level.
      \param a_dx          Grid spacing of that region.
  */
  virtual Real refinedFillParentDx(const Box& a_ghostRegion, const Real& a_dx) const;

  /** Fill the graph and the nodes of a box on a level that was never made, by cutting the cells
      of a coarser level that was.

      The parent cells are handed in rather than looked up, because they live on the index space
      rather than on the service and because which rank holds them is the index space's business.
      \param a_regIrregCovered Regular/Covered/Irregular cells.
      \param a_nodes           Nodes.
      \param a_validRegion     Box wanted.
      \param a_ghostRegion     Box wanted, including its ghost cells.
      \param a_domain          Domain of the level being filled.
      \param a_probLo          Lower-left corner of the computational domain.
      \param a_dx              Grid spacing of the level being filled.
      \param a_parents         Values kept for the cut cells of the coarser level.
      \param a_parentDx        Grid spacing of that coarser level.
      \return False if the box could not be cut, in which case nothing is written.
  */
  virtual bool fillRefinedGraph(BaseFab<int>&           a_regIrregCovered,
                                Vector<IrregNode>&      a_nodes,
                                const Box&              a_validRegion,
                                const Box&              a_ghostRegion,
                                const ProblemDomain&    a_domain,
                                const RealVect&         a_probLo,
                                const Real&             a_dx,
                                const BaseIVFAB<Real>&  a_parents,
                                const Real&             a_parentDx) const;

  virtual void makeGrids( const ProblemDomain&      a_domain,
                          DisjointBoxLayout&        a_grids,
                          const int&                a_maxGridSize,
                          const int&                a_maxIrregGridSize )
  { // default implementation of makeGrids is to ignore it.  GeometryService is ignorant of BoxLayouts
    // or load balances and just answers geometric queries
  }

protected:

};
#include "NamespaceFooter.H"
#endif
