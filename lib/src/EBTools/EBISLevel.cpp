#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "parstream.H"
#include "memtrack.H"
#include "CH_Attach.H"
#include "IndexTM.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "LayoutIterator.H"
#include "BRMeshRefine.H"
#include "AMRIO.H"

#include "EBCFCopy.H"
#include "EBISLevel.H"
#include "EBIndexSpace.H"
#include "EBGraphFactory.H"
#include "EBDataFactory.H"
#include "BaseIVFactory.H"
#include "LoadBalance.H"
#include "EBISLayout.H"
#include "VoFIterator.H"
#include "IrregNode.H"
#include "AllRegularService.H"
#include "PolyGeom.H"
#include "EBLevelDataOps.H"
#include "FaceIterator.H"
#include "NamespaceHeader.H"


EBIndexSpace* Chombo_EBIS::s_instance = NULL;
bool          Chombo_EBIS::s_aliased  = false;
int EBISLevel::s_ebislGhost = 6;
bool EBISLevel::s_distributedData = false;
EBIndexSpace* Chombo_EBIS::instance()
{
  if ((!s_aliased) && (s_instance == NULL))
  {
    s_instance = new EBIndexSpace();
  }

  return  s_instance;
}
////
void 
EBISLevel::
checkGraph() const
{
#ifndef NDEBUG
  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & graph  = m_graph[dit()];
      const Box     & grid   = m_grids[dit()];
      IntVectSet ivsgrid(grid);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          for(FaceIterator faceit(ivsgrid, graph, idir, FaceStop::SurroundingNoBoundary); faceit.ok(); ++faceit)
            {
              const FaceIndex& face = faceit();
              for(SideIterator sit; sit.ok(); ++sit)
                {
                  const IntVect& iv = face.gridIndex(sit());
                  if(graph.getRegion().contains(iv) && graph.isCovered(iv))
                    {
                      pout() << "cell " << iv << "is both covered and part of a face" << endl;
                      MayDay::Error("inconsistent graph description");
                    }
                }
            }
        }
    }
#endif
}
////
void Chombo_EBIS::alias(const EBIndexSpace* a_input)
{
  s_instance = (EBIndexSpace*)(a_input);
  s_aliased  = true;
}

Real EBISLevel::s_tolerance = 1.0e-12;
bool EBISLevel::s_verbose   = false;
bool EBISLevel::s_recursive = false;

long long EBISLevel::numVoFsOnProc() const
{
  CH_TIME("EBISLevel::numVoFsOnProc");
  long long retval = 0;
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& ebgraph = m_graph[dit()];
      long long numVoFsBox = ebgraph.numVoFs(m_grids.get(dit()));
      retval += numVoFsBox;
    }

  return retval;
}

Real EBISLevel::totalVolFracOnProc() const
{
  CH_TIME("EBISLevel::totalVolFracOnProc");
  Real retval = 0.0;

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& ebgraph = m_graph[dit()];
      const EBData&  ebdata  = m_data[dit()];

      EBISBox ebisbox;
      ebisbox.define(ebgraph,ebdata,dit());

      const Box& gridBox = m_grids.get(dit());

      Real totalVolFrac = 0.0;

      if (ebgraph.isRegular(gridBox))
      {
        totalVolFrac = gridBox.numPts();
      }
      else if (!ebgraph.isCovered(gridBox))
      {
        IntVectSet ivs(gridBox);

        for (VoFIterator vofit(ivs,ebgraph); vofit.ok(); ++vofit)
        {
          totalVolFrac += ebisbox.volFrac(vofit());
        }
      }

      retval += totalVolFrac;
    }

  return retval;
}

///
bool isPowerOfTwo (int x)
{
  while (((x % 2) == 0) && (x > 1)) /* While x is even and > 1 */
    {
      x /= 2;
    }
 return (x == 1);
}
///
void EBISLevel::makeLoads(Vector<unsigned long long>&       a_loads,
                          Vector<Box>&                      a_boxes,
                          const Box&                        a_region,
                          const ProblemDomain&              a_domain,
                          const GeometryService&            a_geoserver,
                          const RealVect&                   a_origin,
                          const Real&                       a_dx,
                          const int                         a_ncellmax)
{
  CH_TIME("EBISLevel::makeLoads");
#ifdef CH_MPI
  if(EBIndexSpace::s_useMemoryLoadBalance)
    {
      pout() << "using memory for load balance" << endl;
      for (int i = 0; i < a_boxes.size(); i++)
        {
          if (a_boxes[i].ixType() ==  IndexType::TheNodeType())
            {
              a_boxes[i].convert(IndexType::TheCellType());
            }
        }
      Vector<int> procs;
      LoadBalance(procs, a_boxes);
      DisjointBoxLayout dbl(a_boxes, procs);
      DataIterator dit = dbl.dataIterator();

      dit.enablePeak();
      dit.clearPeak();
      EBGraphFactory graphfact(a_domain);
      EBDataFactory   datafact;
      LevelData<EBGraph> graph(dbl, 1, IntVect::Unit, graphfact);
      LayoutData<Vector<IrregNode> > allNodes(dbl);
      defineGraphFromGeo(graph, allNodes, a_geoserver, dbl,      
                         a_domain,  a_origin, a_dx);

      LevelData<EBData>   data(dbl, 1, IntVect::Zero,  datafact);
      for (dit.reset(); dit.ok(); ++dit)
        {
          data[dit()].defineVoFData( graph[dit()], dbl.get(dit()));
          data[dit()].defineFaceData(graph[dit()], dbl.get(dit()));
        }
      dit.disablePeak();
      dit.mergePeak();
      a_loads = dit.getPeak();
    }
  else
#endif
    {
      pout() << "using old ebis load balance" << endl;
      for (int i = 0; i < a_boxes.size(); i++)
        {
          if (a_boxes[i].ixType() ==  IndexType::TheNodeType())
            {
              a_boxes[i].convert(IndexType::TheCellType());
              a_loads[i] = 8;
            }
          else
            {
              a_loads[i] = 1;
            }
        }
    }
}
///
void EBISLevel::makeBoxes(Vector<Box>&               a_boxes,
                          Vector<unsigned long long>&              a_loads,
                          const Box&                 a_region,
                          const ProblemDomain&       a_domain,
                          const GeometryService&     a_geoserver,
                          const RealVect&            a_origin,
                          const Real&                a_dx,
                          const int                  a_ncellmax)
{
  bool allPowersOfTwo = true;
  for(int idir = 0; idir < CH_SPACEDIM; idir++)
    {
      bool powerOfTwoThisDir = isPowerOfTwo(a_domain.size(idir));
      allPowersOfTwo = allPowersOfTwo && powerOfTwoThisDir;
    }
  if(allPowersOfTwo && s_recursive) //the recursive makeboxes really only likes powers of two
    {
      pout() << "EBISLevel::makeBoxes -- doing recursive" << endl;
      std::list<Box> boxes;
      makeBoxes(boxes, a_region, a_domain, a_geoserver, a_origin, a_dx, a_ncellmax);
      a_boxes.resize(boxes.size());
      a_loads.resize(boxes.size());
      std::list<Box>::iterator it = boxes.begin();
      for (int i = 0; i < a_boxes.size(); ++i, ++it)
        {
          a_boxes[i]=*it;
        }
      //mortonOrdering(a_boxes);
    }
  else
    {
      if (EBIndexSpace::s_MFSingleBox)
        {
          pout() << "EBISLevel::makeBoxes -- doing single box" << endl;
          a_boxes.resize(1);
          a_loads.resize(1);
          a_boxes[0] = a_region;
          a_loads[0] = 1;
          return;
        }
      pout() << "EBISLevel::makeBoxes -- doing domain split" << endl;
      domainSplit(a_domain, a_boxes, a_ncellmax, 1);
      //mortonOrdering(a_boxes);
      a_loads.resize(a_boxes.size(), 1);
    }
  makeLoads(a_loads, a_boxes, a_region, a_domain, a_geoserver, a_origin, a_dx, a_ncellmax);
}

void EBISLevel::makeBoxes(std::list<Box>&        a_boxes,
                          const Box&             a_region,
                          const ProblemDomain&   a_domain,
                          const GeometryService& a_geoserver,
                          const RealVect&        a_origin,
                          const Real&            a_dx,
                          const int              a_ncellmax)
{
  int longdir;
  int length = a_region.longside(longdir);

  if (length > a_ncellmax)
    {
      int n = length/2;
      //CH_assert(n*2==length);
      Box low(a_region), high;
      high = low.chop(longdir, a_region.smallEnd(longdir)+n);
      makeBoxes(a_boxes, low,  a_domain,
                a_geoserver, a_origin, a_dx, a_ncellmax);
      makeBoxes(a_boxes, high, a_domain,
                a_geoserver, a_origin, a_dx, a_ncellmax);
    }
  else
    {
      if (a_geoserver.InsideOutside(a_region, a_domain, a_origin, a_dx) == GeometryService::Irregular)
        {
          Box n = a_region;
          n.convert(IndexType::TheNodeType());
          a_boxes.push_back(n);
        }
      else
        {
          a_boxes.push_back(a_region);
        }
    }
}

#ifdef CH_USE_HDF5
EBISLevel::EBISLevel(HDF5Handle& a_handle)
{
  CH_TIME("EBISLevel::EBISLevel_hdf5");
  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;
  m_hasSurface  = false;


  HDF5HeaderData header;
  header.readFromFile(a_handle);
  m_origin = header.m_realvect["EBIS_origin"];
  m_domain = header.m_box     ["EBIS_domain"];
  m_dx =     header.m_real    ["EBIS_dx"]    ;
  m_tolerance = m_dx*1E-4;
  m_level = 0;

  //read in the grids
  Vector<Box> boxes;
  read(a_handle,boxes);
  Vector<int> procAssign;
  LoadBalance(procAssign, boxes);
  m_grids.define(boxes, procAssign);//this should use m_domain for periodic...
  EBGraphFactory graphfact(m_domain);
  m_graph.define(m_grids, 1, IntVect::Zero, graphfact);

  //read the graph  in from the file
  std::string graphName("EBIS_graph");
  int eekflag = read(a_handle, m_graph, graphName, m_grids, Interval(), false);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }

  //need a ghosted layout so that the data can be defined properly
  LevelData<EBGraph> ghostGraph(m_grids, 1, IntVect::Unit, graphfact);
  Interval interv(0,0);
  m_graph.copyTo(interv, ghostGraph, interv);

  //now the data for the graph
  EBDataFactory dataFact;
  m_data.define(m_grids, 1, IntVect::Zero, dataFact);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_data[dit()].defineVoFData(ghostGraph[dit()], m_grids.get(dit()));
      m_data[dit()].defineFaceData(ghostGraph[dit()], m_grids.get(dit()));
    }
  //read the data  in from the file
  std::string  dataName("EBIS_data");

  eekflag = read(a_handle, m_data ,  dataName, m_grids, Interval(), false);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
#if 0
  pout() << "EBISLevel::EBISLevel 1 - m_grids - m_dx: " << m_dx << endl;
  pout() << "--------" << endl;
  pout() << m_grids.boxArray().size() << endl;
  pout() << "--------" << endl;
  pout() << endl;
#endif
}

void EBISLevel::write(HDF5Handle& a_handle) const
{
  CH_TIME("EBISLevel::write");
  HDF5HeaderData header;
  //this naming stuff kinda depends on the fact
  //that we are only outputting the finest level.
  //we could be slick and incorporate the
  //level number in there if we wanted.
  header.m_int["num_levels"] = 1;
  header.m_int["num_components"] = 1;
  header.m_string["component_0"] = "phi0";
  header.m_realvect["EBIS_origin"] = m_origin;
  header.m_box     ["EBIS_domain"] = m_domain.domainBox();
  header.m_real    ["EBIS_dx"]     = m_dx;
  header.writeToFile(a_handle);
  //write the grids to the file
  CH_XD::write(a_handle, m_grids);

  std::string graphName("EBIS_graph");
  std::string  dataName("EBIS_data");
  int eekflag = CH_XD::write(a_handle, m_graph, graphName);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }
  eekflag = CH_XD::write(a_handle, m_data ,  dataName);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
}
#endif
///
void
EBISLevel::defineGraphFromGeo(LevelData<EBGraph>             & a_graph,
                              LayoutData<Vector<IrregNode> > & a_allNodes,
                              const GeometryService          & a_geoserver,
                              const DisjointBoxLayout        & a_grids,
                              const ProblemDomain            & a_domain,
                              const RealVect                 & a_origin,
                              const Real                     & a_dx)
{
  CH_TIME("EBISLevel::defineGraphFromGeo");
  //define the graph stuff
  const DataIterator dit = a_grids.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) 
    {
      const DataIndex din = dit[mybox];
      
      const Box region = grow(a_grids[din],1) & a_domain;
      const Box ghostRegion = grow(region, 1) & a_domain;

      EBGraph& ebgraph = a_graph[din];
      GeometryService::InOut inout;

      inout = a_geoserver.InsideOutside(region, a_domain, a_origin, a_dx, din);
  
      if (inout == GeometryService::Regular)
        {
          ebgraph.setToAllRegular();
        }
      else if (inout == GeometryService::Covered)
        {
          ebgraph.setToAllCovered();
        }
      else
        {
          BaseFab<int>       regIrregCovered(ghostRegion, 1);
          Vector<IrregNode>&  nodes = a_allNodes[din];

          // if (!a_distributedData)
          //   {
          //     a_geoserver.fillGraph(regIrregCovered, nodes, region,
          //                           ghostRegion, a_domain,
          //                           a_origin, a_dx);
          //   }
          // else
          //   {
          a_geoserver.fillGraph(regIrregCovered, nodes, region,
                                ghostRegion, a_domain,
                                a_origin, a_dx, din);
          // }
          ebgraph.buildGraph(regIrregCovered, nodes, region, a_domain);
          
        }
    }
}

void EBISLevel::simplifyGraphFromGeo(LevelData<EBGraph>             & a_graph,
				     const GeometryService          & a_geoserver,
				     const DisjointBoxLayout        & a_grids,
				     const ProblemDomain            & a_domain,
				     const RealVect                 & a_origin,
				     const Real                     & a_dx)
{
  CH_TIME("EBISLevel::simplifyGraphFromGeo");

  //define the graph stuff
  const DataIterator dit = a_grids.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++)
    {
      const DataIndex& din = dit[mybox];
    
      const Box region = grow(a_grids[din],1) & a_domain;

      EBGraph& ebgraph = a_graph[din];
      GeometryService::InOut inout = a_geoserver.InsideOutside(region, a_domain, a_origin, a_dx, din);

      if (inout == GeometryService::Regular)
        {
          ebgraph.setToAllRegular();
        }
      else if (inout == GeometryService::Covered)
        {
          ebgraph.setToAllCovered();
        }
    }
}

EBISLevel::EBISLevel(const ProblemDomain   & a_domain,
                     const RealVect        & a_origin,
                     const Real            & a_dx,
                     const GeometryService & a_geoserver,
                     int                     a_nCellMax,
                     const bool            & a_fixRegularNextToMultiValued)
{
  // this is the method called by EBIndexSpace::buildFirstLevel
  CH_TIME("EBISLevel::EBISLevel_geoserver_domain");
  pout() << "Entering EBISLevel::EBISLevel called by EBIndexSpace::buildFirstLevel..." << endl;

  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;
  m_hasSurface  = false;

  m_domain = a_domain;
  m_dx = a_dx;
  m_tolerance = a_dx*1E-4;
  m_origin = a_origin;

  m_level = 0;

  m_geoserver = &a_geoserver;

  Vector<Box> vbox;
  Vector<unsigned long long> irregCount;

  if(!s_distributedData){
    {
      CH_TIME("EBISLevel::EBISLevel_makeboxes");
      makeBoxes(vbox,
		irregCount,
		a_domain.domainBox(),
		a_domain,
		a_geoserver,
		a_origin,
		a_dx,
		a_nCellMax);
    }

    // pout()<<vbox<<"\n\n";
    //load balance the boxes
    Vector<int> procAssign;
    //UnLongLongLoadBalance(procAssign, irregCount, vbox);
    basicLoadBalance(procAssign, vbox.size());
    //   pout()<<irregCount<<std::endl;
    //   pout()<<procAssign<<std::endl;
    pout() << "before defining grids" << endl;
    m_grids.define(vbox, procAssign,a_domain);//this should use a_domain for periodic
    pout() << "after defining grids" << endl;
  }
  else
    {
      CH_TIME("EBISLevel::EBISLevel_makegrids");
      // permit the geometry service to construct a layout, or accept an already defined layout from EBIndexSpace

      (const_cast<GeometryService*>(&a_geoserver))->makeGrids(a_domain, m_grids, a_nCellMax, 15);
    }

  RealVect dx2D;
  for (int i = 0; i < SpaceDim; i++)
    {
      dx2D[i]=a_dx;
    }

  (const_cast<GeometryService*>(&a_geoserver))->postMakeBoxLayout(m_grids,dx2D);
  LayoutData<Vector<IrregNode> > allNodes(m_grids);

  EBGraphFactory graphfact(a_domain);
  m_graph.define(m_grids, 1, IntVect::Unit, graphfact);

  defineGraphFromGeo(m_graph, allNodes, a_geoserver, m_grids,
                     m_domain,m_origin, m_dx);

  checkGraph();

  EBDataFactory dataFact;
  m_data.define(m_grids, 1, IntVect::Zero, dataFact);

  const DataIterator& dit = m_grids.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++)
    {
      const DataIndex& din = dit[mybox];
      
      m_data[din].define(m_graph[din], allNodes[din], m_grids[din]);
    }

  if (a_geoserver.canGenerateMultiCells())
    {
      if (a_fixRegularNextToMultiValued)
        {
          fixRegularNextToMultiValued();
        }
    }
  pout() << "Exiting EBISLevel::EBISLevel called by EBIndexSpace::buildFirstLevel..." << endl;
}

//now fix the multivalued next to regular thing for the graph and the data
//the oldgraph/newgraph thing is necessary because the graphs are
//reference counted and they have to be kept consistent with the data
int EBISLevel::numSurfaceComponents() const
{
  return (m_geoserver == NULL) ? 0 : m_geoserver->numSurfaceComponents();
}

void EBISLevel::defineSurfaces()
{
  CH_TIME("EBISLevel::defineSurfaces");

  if (m_hasSurface)
    {
      return;
    }

  m_hasSurface = true;

  const int ncomp = this->numSurfaceComponents();

  if (ncomp <= 0)
    {
      // The geometry service keeps nothing, which is every service but the polyhedral one, and
      // which is what makes a level built by any of them unable to be refined.
      return;
    }

  LayoutData<IntVectSet>        sets(m_grids);
  LayoutData<Vector<IntVect> >  cells(m_grids);
  LayoutData<Vector<Real> >     values(m_grids);

  DataIterator dit = m_grids.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex din = dit();

      m_geoserver->getSurfaces(cells[din], values[din], m_grids[din], m_dx);

      IntVectSet ivs;

      for (int n = 0; n < cells[din].size(); n++)
        {
          ivs |= cells[din][n];
        }

      sets[din] = ivs;
    }

  // BaseIVFactory reads the graph through an EBISLayout, and this level's own grids with no
  // ghost cells is the layout whose graph is exactly m_graph.  It costs a copy of the level
  // while this runs; the alternative is a factory that takes a LevelData<EBGraph> directly,
  // which is worth adding if this ever shows up in a profile.
  EBISLayout ebisl;
  this->fillEBISLayout(ebisl, m_grids, 0);

  BaseIVFactory<Real> factory(ebisl, sets);

  m_surface.define(m_grids, ncomp, IntVect::Zero, factory);

  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex din = dit();

      BaseIVFAB<Real>&       fab      = m_surface[din];
      const Vector<IntVect>& theCells = cells[din];
      const Vector<Real>&    theVals  = values[din];

      for (int n = 0; n < theCells.size(); n++)
        {
          // the generator mandates single-valued cut cells, so each one is a single volume
          const VolIndex vof(theCells[n], 0);

          for (int comp = 0; comp < ncomp; comp++)
            {
              fab(vof, comp) = theVals[n*ncomp + comp];
            }
        }
    }
}

void EBISLevel::extendTo(const Vector<Box>& a_newBoxes, EBISLevel& a_coarser)
{
  CH_TIME("EBISLevel::extendTo");

  CH_assert(a_newBoxes.size() > 0);
  CH_assert(m_geoserver != NULL);

  // The parents have to be on the coarser level rather than on the geometry service, since the
  // service keeps them on layouts of its own and a box being cut wants the parents of its ghost
  // cells too -- which generally sit on another rank.
  a_coarser.defineSurfaces();

  const int ncomp = a_coarser.numSurfaceComponents();

  CH_assert(ncomp > 0);

  const Interval interv(0, ncomp - 1);

  Vector<Box> newBoxes = a_newBoxes;
  Vector<int> newRanks;

  LoadBalance(newRanks, newBoxes);

  DisjointBoxLayout newGrids(newBoxes, newRanks, m_domain);

  // The parents of a box reach one coarse cell past its own coarsening, since the graph is built
  // over the box grown by one and the data asks the graph for one further still.  Two is taken
  // rather than one so that the parent graph reaches past the parent data, which is what defining
  // face data needs.
  Vector<Box> parentBoxes(newBoxes.size());

  for (int i = 0; i < newBoxes.size(); i++)
    {
      CH_assert(newBoxes[i] == refine(coarsen(newBoxes[i], 2), 2));

      parentBoxes[i] = coarsen(newBoxes[i], 2);
    }

  DisjointBoxLayout parentGrids(parentBoxes, newRanks, a_coarser.m_domain);

  const int     parentGhost     = 2;
  const IntVect parentGhostVect = parentGhost*IntVect::Unit;

  EBISLayout parentLayout;
  a_coarser.fillEBISLayout(parentLayout, parentGrids, parentGhost);

  // A surface was kept for every cell the generator called irregular, so the graph is what says
  // where they are -- which is also what lays the destination out to receive them.
  LayoutData<IntVectSet> parentSets(parentGrids);

  DataIterator pdit = parentGrids.dataIterator();

  for (pdit.begin(); pdit.ok(); ++pdit)
    {
      const Box grown = grow(parentGrids[pdit()], parentGhostVect) & a_coarser.m_domain.domainBox();

      parentSets[pdit()] = parentLayout[pdit()].getEBGraph().getIrregCells(grown);
    }

  BaseIVFactory<Real> parentFact(parentLayout, parentSets);

  LevelData<BaseIVFAB<Real> > parents(parentGrids, ncomp, parentGhostVect, parentFact);

  {
    Copier copier(a_coarser.m_grids, parentGrids, a_coarser.m_domain, parentGhostVect);

    a_coarser.m_surface.copyTo(interv, parents, interv, copier);
  }

  // cut the new boxes
  EBGraphFactory graphfact(m_domain);
  EBDataFactory  datafact;

  LevelData<EBGraph> newGraph(newGrids, 1, IntVect::Unit, graphfact);
  LevelData<EBData>  newData (newGrids, 1, IntVect::Zero, datafact);

  DataIterator ndit = newGrids.dataIterator();

  for (ndit.begin(); ndit.ok(); ++ndit)
    {
      const DataIndex din = ndit();

      const Box region      = grow(newGrids[din], 1) & m_domain.domainBox();
      const Box ghostRegion = grow(region,        1) & m_domain.domainBox();

      BaseFab<int>      regIrregCovered(ghostRegion, 1);
      Vector<IrregNode> nodes;

      const bool cut = m_geoserver->fillRefinedGraph(regIrregCovered, nodes, region, ghostRegion,
                                                     m_domain, m_origin, m_dx,
                                                     parents[din], a_coarser.m_dx);

      if (!cut)
        {
          MayDay::Error("EBISLevel::extendTo - the geometry service could not cut a box it was asked for");
        }

      newGraph[din].buildGraph(regIrregCovered, nodes, region, m_domain);
      newData [din].define(newGraph[din], nodes, newGrids[din]);
    }

  // the level is what it was, plus the new boxes
  Vector<Box> allBoxes;
  Vector<int> allRanks;

  for (LayoutIterator lit = m_grids.layoutIterator(); lit.ok(); ++lit)
    {
      allBoxes.push_back(m_grids[lit()]);
      allRanks.push_back(m_grids.procID(lit()));
    }

  for (int i = 0; i < newBoxes.size(); i++)
    {
      allBoxes.push_back(newBoxes[i]);
      allRanks.push_back(newRanks[i]);
    }

  DisjointBoxLayout allGrids(allBoxes, allRanks, m_domain);

  const Interval one(0, 0);

  // EBData is laid out against its graph and keeps a reference to it, so the two have to be
  // rebuilt together: defining the data against a temporary graph and then replacing this
  // level's graph underneath it leaves the data describing a graph the level no longer holds.
  // That is the "oldgraph/newgraph" hazard this file warns about above.
  LevelData<EBGraph> oldGraph;
  oldGraph.define(m_graph, graphfact);

  LevelData<EBData> oldData;
  oldData.define(m_data, datafact);

  m_grids = allGrids;

  m_graph.define(m_grids, 1, IntVect::Unit, graphfact);

  oldGraph.copyTo(one, m_graph, one);
  newGraph.copyTo(one, m_graph, one);

  m_data.define(m_grids, 1, IntVect::Zero, datafact);

  DataIterator adit = m_grids.dataIterator();

  for (adit.begin(); adit.ok(); ++adit)
    {
      m_data[adit()].defineVoFData (m_graph[adit()], m_grids[adit()]);
      m_data[adit()].defineFaceData(m_graph[adit()], m_grids[adit()]);
    }

  oldData.copyTo(one, m_data, one);
  newData.copyTo(one, m_data, one);


  // what was kept for the old cells no longer covers the level, and the cache is laid out against
  // grids that have changed
  m_surface.clear();
  m_hasSurface = false;
  m_cache.clear();
}

void EBISLevel::attachFinerNodesFrom(EBISLevel& a_finer, const TreeIntVectSet& a_cells)
{
  CH_TIME("EBISLevel::attachFinerNodesFrom");

  DisjointBoxLayout fineFromCoar;
  refine(fineFromCoar, m_grids, 2);
  fineFromCoar.close();

  EBGraphFactory fineFact(a_finer.m_domain);

  LevelData<EBGraph> sharedFineGraph(fineFromCoar, 1, IntVect::Zero, fineFact);

  const Interval one(0, 0);

  a_finer.m_graph.copyTo(one, sharedFineGraph, one);

  const IntVectSet wanted(a_cells);

  long long disagreed = 0;

  DataIterator dit = m_grids.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex din = dit();

      IntVectSet cells = wanted;
      cells &= m_grids[din];

      if (cells.isEmpty())
        {
          continue;
        }

      disagreed += m_graph[din].attachFinerNodes(sharedFineGraph[din], cells);
    }

  disagreed = EBLevelDataOps::parallelSum(disagreed);

  if (disagreed > 0)
    {
      pout() << "    " << disagreed
             << " extended cells hold a different number of vofs than the level under them" << endl;
    }
}

void EBISLevel::fixRegularNextToMultiValued()
{
  CH_TIME("EBISLevel::fixRegularNextToMultiValued");

  EBGraphFactory graphfact(m_domain);
  LayoutData<IntVectSet> vofsToChange(m_grids);
  LevelData<EBGraph> oldGhostGraph(m_grids, 1, 2*IntVect::Unit, graphfact);
  Interval interv(0,0);

  if(s_distributedData){
    this->simplifyGraphFromGeo(oldGhostGraph, *m_geoserver, m_grids, m_domain, m_origin, m_dx);
  }

  m_graph.copyTo(interv, oldGhostGraph, interv);

  const DataIterator dit = m_grids.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++)
    {
      const DataIndex& din = dit[mybox];
      
      CH_TIME("EBISLevel::fixRegularNextToMultiValued_loop1");
      m_graph[din].getRegNextToMultiValued(vofsToChange[din],
                                             oldGhostGraph[din]);

      //only this box's own cells are converted here.  a ghost cell belongs to the box that holds
      //it and is converted there, and the graph step already leaves it alone since it lies
      //outside the graph this box owns.  the moments beside it reach one cell out, and the
      //multivalued face this step then looks across would reach a second
      vofsToChange[din] &= m_grids[din];

      m_graph[din].addFullIrregularVoFs(vofsToChange[ din],
                                          oldGhostGraph[din]);
    }

  EBDataFactory datafact;
  LevelData<EBGraph> newGhostGraph(m_grids, 1, 2*IntVect::Unit, graphfact);
  LevelData<EBData>  newGhostData(m_grids, 1, IntVect::Unit, datafact);

  if(s_distributedData){ 
    simplifyGraphFromGeo(newGhostGraph, *m_geoserver, m_grids, m_domain, m_origin, m_dx);
  }

  m_graph.copyTo(interv, newGhostGraph, interv);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++)
    {
      const DataIndex& din = dit[mybox]; 

      CH_TIME("EBISLevel::fixRegularNextToMultiValued_loop2");

      const Box localBox = grow(m_grids[din],1) & m_domain;

      newGhostData[din].defineVoFData(oldGhostGraph[din],  localBox);
      newGhostData[din].defineFaceData(oldGhostGraph[din], localBox);
    }

  m_data.copyTo(interv,  newGhostData,  interv);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++)
    {
      const DataIndex& din = dit[mybox]; 

      CH_TIME("EBISLevel::fixRegularNextToMultiValued_loop3");

      m_data[din].addFullIrregularVoFs(vofsToChange[din],
				       newGhostGraph[din],
				       newGhostData[din].getVolData(),
				       oldGhostGraph[din]);
    }
}

//checks to see the vofs are in the correct cells.
//checks to see that the faces are over the correct cells
//checks that volume fractions, area fractions are positive
//bail out with MayDay::Error if anything fails
void EBISLevel::sanityCheck(const EBIndexSpace* const a_ebisPtr)
{
#if 0
  pout() << "EBISLevel::sanityCheck" << endl;
  CH_TIME("EBISLevel::sanityCheck");
  EBISLayout ghostLayout;

  a_ebisPtr->fillEBISLayout(ghostLayout, m_grids, m_domain, 1);

  Real maxcentval = 0.5+ s_tolerance;
  for (DataIterator dit = m_grids.dataIterator();  dit.ok(); ++dit)
    {
      const Box& thisBox = m_grids.get(dit());
      const EBISBox& ebisBox = ghostLayout[dit()];
      for (BoxIterator bit(thisBox); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          Vector<VolIndex> vofs = ebisBox.getVoFs(iv);
          for (int ivof = 0; ivof < vofs.size(); ivof++)
            {
              const VolIndex& vof = vofs[ivof];
              if (vof.gridIndex() != iv)
                {
                  pout() << "EBISLevel::sanityCheck: Error" << endl;
                  pout() << "VoF at Intvect = " << iv
                         << "has grid index = " << vof.gridIndex() << endl;
                  MayDay::Error("EBISLevel::sanityCheck Error 1");
                }
              if (vof.cellIndex() < 0)
                {
                  pout() << "EBISLevel::sanityCheck: Error" << endl;
                  pout() << "VoF at Intvect = " << iv
                         << "has negative cell index = " << vof.cellIndex() << endl;
                  MayDay::Error("EBISLevel::sanityCheck Error 2");
                }
              Real volFrac = ebisBox.volFrac(vof);
              if (volFrac < 0.0)
                {
                  pout() << "EBISLevel::sanityCheck: Error" << endl;
                  pout() << "VoF at Intvect = " << iv
                         << "has invalid volume fraction = " << volFrac << endl;
                  MayDay::Error("EBISLevel::sanityCheck Error 5");
                }
              RealVect volCentroid = ebisBox.centroid(vof);
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  Real volcentdir = volCentroid[idir];
                  if (volFrac > s_tolerance)
                    {
                      if (volcentdir > maxcentval || volcentdir < -maxcentval)
                        {
                          pout() << "EBISLevel::sanityCheck: Error" << endl;
                          pout() << "VoF at Intvect = " << iv
                                 << " has invalid vol centroid = " << volcentdir
                                 << " at direction "<< idir << endl;
                          MayDay::Error("EBISLevel::sanityCheck Error 51");
                        }
                    }
                }
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      Vector<FaceIndex> faces = ebisBox.getFaces(vof, idir, sit());
                      IntVect iv2 = iv + sign(sit())*BASISV(idir);
                      ////check for regular next to covered and multivalued next to regular
                      if (m_domain.contains(iv2))
                        {
                          if (ebisBox.isRegular(iv))
                            {
                              if (ebisBox.isCovered(iv2))
                                {
                                  pout() << iv << " is regular and " <<  iv2 << " is covered" << endl;
                                  MayDay::Error("EBISLevel::sanityCheck error 420 ");
                                }
                              else
                                {
                                  Vector<VolIndex> otherVoFs = ebisBox.getVoFs(iv2);
                                  if (otherVoFs.size() > 1)
                                    {
                                      pout() << iv << " is regular and " <<  iv2 << " is multivalued" << endl;
                                      MayDay::Error("EBISLevel::sanityCheck error 420.2 ");
                                    }
                                }
                            }
                        }
                      IntVect ivlo, ivhi;
                      if (sit() == Side::Lo)
                        {
                          ivlo = iv2;
                          ivhi = iv;
                        }
                      else
                        {
                          ivlo = iv;
                          ivhi = iv2;
                        }
                      for (int iface = 0; iface < faces.size(); iface++)
                        {
                          const FaceIndex& face = faces[iface];
                          if (face.gridIndex(Side::Lo) != ivlo)
                            {
                              pout() << "EBISLevel::sanityCheck: Error" << endl;
                              pout() << "face at IntVects = " << ivlo << "  " << ivhi
                                     << "has low IntVect  = " << face.gridIndex(Side::Lo)
                                     << endl;
                              MayDay::Error("EBISLevel::sanityCheck Error 3");
                            }
                          if (face.gridIndex(Side::Hi) != ivhi)
                            {
                              pout() << "EBISLevel::sanityCheck: Error" << endl;
                              pout() << "face at IntVects = " << ivlo << "  " << ivhi
                                     << "has high IntVect = " << face.gridIndex(Side::Hi)
                                     << endl;
                              MayDay::Error("EBISLevel::sanityCheck Error 4");
                            }
                          Real areaFrac = ebisBox.areaFrac(face);
                          if (areaFrac  < 0.0)
                            {
                              pout() << "EBISLevel::sanityCheck: Error" << endl;
                              pout() << "VoF at Intvect = " << iv
                                     << "has invalid area fraction = " << areaFrac << endl;
                              MayDay::Error("EBISLevel::sanityCheck Error 51");
                            }
                          if (areaFrac  >  s_tolerance)
                            {
                              RealVect faceCentroid = ebisBox.centroid(face);
                              for (int idir = 0; idir < SpaceDim; idir++)
                                {
                                  Real facecentdir = faceCentroid[idir];
                                  if (facecentdir > maxcentval || facecentdir < -maxcentval)
                                    {
                                      pout() << "EBISLevel::sanityCheck: Error" << endl;
                                      pout() << "VoF at Intvect = " << iv
                                             << " has invalid face centroid = " << facecentdir
                                             << " at direction "<< idir << endl;
                                      MayDay::Error("EBISLevel::sanityCheck Error 51");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif
}

EBISLevel::EBISLevel()
{
  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;
  m_hasSurface  = false;

  m_level = 0;

}

//steps to coarsen an ebislevel:
//1. coarsen vofs
//1a. make a layout over refine(mydbl,2) and copy
//    fine layout into it
//1b.do connectivity bizbaz to make my vofs, volfrac, vof->fineVofs
//2. make faces doing connectivity jive
//2.23 make coarse geometric stuff (centroids and all that) from fine
//3. make coarse layout from coarsen(finelayout). and copy my data into it
//   to make fine->coarserVoF
// (so finer ebislevel does change in this function)
void EBISLevel::dumpDebug(const string& a_string)
{
  if (m_domain.domainBox() == EBGraphImplem::s_doDebug)
    {
      pout() << a_string <<  ": ";
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          if (m_grids[dit()].contains(EBGraphImplem::s_ivDebug))
            {
              pout() << "EBIS1: " << EBGraphImplem::s_ivDebug;
              if (m_graph[dit()].isRegular(EBGraphImplem::s_ivDebug))
                {
                  pout() << " is regular" << endl;
                }
              else if (m_graph[dit()].isCovered(EBGraphImplem::s_ivDebug))
                {
                  pout() << " is covered" << endl;
                }
              else
                {
                  pout() << " is irregular" << endl;
                }
            }
        }

    }
}

void EBISLevel::coarsenFrom(EBISLevel& a_fineEBIS, bool a_fixRegularNextToMultiValued)
{
  CH_TIME("EBISLevel::coarsenFrom");

//  pout() << "before coarsenVoFs " << endl;
  //create coarsened vofs from fine.
  //the fine graph and data on the refinement of these grids, and the coarse graph on these grids,
  //are each read by both coarsening steps.  build them once with the widest ghost region either
  //step needs -- three for the fine graph, two for the fine data, one for the coarse graph -- so
  //that the level is exchanged once per object rather than once per step
  DisjointBoxLayout fineFromCoarDBL;
  refine(fineFromCoarDBL, m_grids, 2);
  fineFromCoarDBL.close();

  // Which of this level's boxes the finer level actually reaches. Coarsening can only speak for
  // those; the rest keep whatever this level was built with.
  IntVectSet fineCoverage;
  {
    const Vector<Box>& fineBoxes = a_fineEBIS.m_grids.boxArray();

    for (int ibox = 0; ibox < fineBoxes.size(); ibox++)
      {
        fineCoverage |= coarsen(fineBoxes[ibox], 2);
      }
  }

  // Coarsening a cell reads what lies under the cells around it as well as under its own, so a
  // cell comes up through coarsening when the finer level reaches a cell past it.  The outermost
  // ring of what the finer level was carried over is therefore not coarsened: it is what the
  // cells inside it read, and it keeps what it was generated with.  Deciding this a cell at a
  // time rather than a box at a time is what keeps that ring one cell thick.  A ring a box thick
  // would have to be paid for by carrying the finer level a whole box further than it is wanted,
  // which for a surface one or two boxes thick is most of it again.
  //
  // The ring is found a box at a time, from the boxes of the finer level rather than from a set
  // over the domain.  Eroding a set is the same thing said globally, and says it far too dearly:
  // what the finer level does not reach is nearly the whole domain, and growing a set that size
  // costs more than everything else here put together.
  Vector<Box> coveredBoxes;
  {
    const Vector<Box>& fineBoxes = a_fineEBIS.m_grids.boxArray();

    for (int ibox = 0; ibox < fineBoxes.size(); ibox++)
      {
        coveredBoxes.push_back(coarsen(fineBoxes[ibox], 2));
      }
  }

  LayoutData<IntVectSet> coarsenCells(m_grids);
  LayoutData<BaseFab<bool> > reached(m_grids);
  {
    const DataIterator& maskDit = m_grids.dataIterator();

    for (int mybox = 0; mybox < maskDit.size(); mybox++)
      {
        const DataIndex din = maskDit[mybox];

        Box grown = m_grids[din];
        grown.grow(1);
        grown &= m_domain;

        reached[din].resize(grown, 1);
        reached[din].setVal(false);

        for (int ibox = 0; ibox < coveredBoxes.size(); ibox++)
          {
            const Box overlap = coveredBoxes[ibox] & grown;

            if (!overlap.isEmpty())
              {
                reached[din].setVal(true, overlap, 0, 1);
              }
          }

        coarsenCells[din] = IntVectSet(DenseIntVectSet(m_grids[din], false));

        for (BoxIterator bit(m_grids[din]); bit.ok(); ++bit)
          {
            bool all = reached[din](bit(), 0);

            for (int idir = 0; idir < SpaceDim && all; idir++)
              {
                for (SideIterator sit; sit.ok(); ++sit)
                  {
                    const IntVect iv = bit() + sign(sit()) * BASISV(idir);

                    if (m_domain.contains(iv) && !reached[din](iv, 0))
                      {
                        all = false;
                      }
                  }
              }

            if (all)
              {
                coarsenCells[din] |= bit();
              }
          }
      }
  }

  // What each rank holds of the mask, so that a cell can ask whether the cell across a box
  // boundary came up through coarsening without anyone holding a set over the whole domain
  LevelData<FArrayBox> coarsenedMask(m_grids, 1, IntVect::Unit);
  {
    const DataIterator& maskDit = m_grids.dataIterator();

    for (int mybox = 0; mybox < maskDit.size(); mybox++)
      {
        const DataIndex din = maskDit[mybox];

        coarsenedMask[din].setVal(0.0);

        for (IVSIterator ivsIt(coarsenCells[din]); ivsIt.ok(); ++ivsIt)
          {
            coarsenedMask[din](ivsIt(), 0) = 1.0;
          }
      }

    coarsenedMask.exchange();
  }

  EBGraphFactory ebgraphfactfine(a_fineEBIS.m_domain);
  EBGraphFactory ebgraphfactcoar(m_domain);
  EBDataFactory  ebdatafactshared;

  LevelData<EBGraph> sharedFineGraph(fineFromCoarDBL, 1, 3*IntVect::Unit, ebgraphfactfine);
  LevelData<EBGraph> sharedCoarGraph(m_grids,         1,   IntVect::Unit, ebgraphfactcoar);
  LevelData<EBData>  sharedFineData (fineFromCoarDBL, 1, 2*IntVect::Unit, ebdatafactshared);

  {
    Interval sharedInterv(0,0);

    if(s_distributedData){
      simplifyGraphFromGeo(sharedFineGraph, *m_geoserver, fineFromCoarDBL, m_domain, m_origin, m_dx);
    }
    a_fineEBIS.m_graph.copyTo(sharedInterv, sharedFineGraph, sharedInterv);

    const DataIterator& sharedDit = m_grids.dataIterator();

    const int sharedNbox = sharedDit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < sharedNbox; mybox++)
      {
        const DataIndex din = sharedDit[mybox];

        Box localBox = grow(fineFromCoarDBL.get(din), 2);
        localBox &= a_fineEBIS.m_domain;
        sharedFineData[din].defineVoFData(sharedFineGraph[din], localBox);
        sharedFineData[din].defineFaceData(sharedFineGraph[din], localBox);
      }

    a_fineEBIS.m_data.copyTo(sharedInterv, sharedFineData, sharedInterv);
  }

  // The cells on the edge of what is about to be coarsened, and how many VoFs they hold now. The
  // cells outside were generated with arcs naming these as they stand; if coarsening changes one
  // of them, those arcs are stale and nothing goes back to rewrite them.
  LayoutData<IntVectSet> edgeCells(m_grids);
  LayoutData<std::vector<long long>> edgeBefore(m_grids);

  {
    const DataIterator& edgeDit = m_grids.dataIterator();

    for (int mybox = 0; mybox < edgeDit.size(); mybox++)
      {
        const DataIndex din = edgeDit[mybox];

        for (IVSIterator ivsIt(coarsenCells[din]); ivsIt.ok(); ++ivsIt)
          {
            const Box neighbourhood = grow(Box(ivsIt(), ivsIt()), 1) & reached[din].box();

            for (BoxIterator nit(neighbourhood); nit.ok(); ++nit)
              {
                if (!reached[din](nit(), 0))
                  {
                    edgeCells[din] |= ivsIt();

                    break;
                  }
              }
          }

        for (IVSIterator ivsIt(edgeCells[din]); ivsIt.ok(); ++ivsIt)
          {
            edgeBefore[din].push_back(m_graph[din].numVoFs(ivsIt()));
          }
      }
  }

  // Coarsening a cell asks its neighbours what lies under them, whether those were coarsened or
  // not: the face between two coarse cells is decided by testing whether the fine vofs under one
  // reach the fine vofs under the other.  A cell in a box that was not coarsened was generated
  // instead, and holds no such record.  Write it here, for the cells the fine level reaches,
  // without touching the cells themselves.  The record is what the neighbour is read for; the
  // cell keeps the vofs, arcs and moments it was generated with.
  {
    long long disagreed = 0;

    const DataIterator& recordDit = m_grids.dataIterator();

    for (int mybox = 0; mybox < recordDit.size(); mybox++)
      {
        const DataIndex din = recordDit[mybox];

        IntVectSet reachable(DenseIntVectSet(m_grids[din], false));

        for (BoxIterator bit(m_grids[din]); bit.ok(); ++bit)
          {
            if (reached[din](bit(), 0))
              {
                reachable |= bit();
              }
          }

        reachable -= coarsenCells[din];

        if (reachable.isEmpty())
          {
            continue;
          }

        disagreed += m_graph[din].attachFinerNodes(sharedFineGraph[din], reachable);
      }

    disagreed = EBLevelDataOps::parallelSum(disagreed);

    if (disagreed > 0)
      {
        pout() << "    " << disagreed
               << " generated cells hold a different number of vofs than the fine level under them" << endl;
      }
  }

  coarsenVoFs(a_fineEBIS, sharedFineGraph, sharedFineData, sharedCoarGraph, coarsenCells);

//  pout() << "before coarsenFacess " << endl;
  //overallMemoryUsage();
  //create coarse faces from fine
  coarsenFaces(a_fineEBIS, sharedFineGraph, sharedFineData, sharedCoarGraph, coarsenCells);

  // The two sides of the edge of what was coarsened read the faces they share differently.  Make
  // them agree while the coarse ghost graph still describes the graph the data sits on
  reconcileSeam(coarsenedMask);
  //overallMemoryUsage();
  //fix the regular next to the multivalued cells
  //to be full irregular cells
  //  dumpDebug(string("EBIS::before FRNTM"));
  if (a_fixRegularNextToMultiValued)
    {
      fixRegularNextToMultiValued();
    }
  //  dumpDebug(string("EBIS::after FRNTM"));

  //overallMemoryUsage();
  // fix the fine->coarseVoF thing.
//  pout() << "before fix fine to coarse " << endl;
  fixFineToCoarse(a_fineEBIS);
  checkGraph();

  // Where a level is only coarsened in part, the cells that were not coarsened were generated
  // with arcs naming their neighbours as they stood. Coarsening may since have changed one of
  // those neighbours -- given it another VoF, or turned it from whole to cut -- and nothing goes
  // back to rewrite the arcs that point at it. Say so here rather than letting a stencil walk a
  // graph that does not join up.
  {
    long long stale = 0;

    const DataIterator& edgeDit = m_grids.dataIterator();

    for (int mybox = 0; mybox < edgeDit.size(); mybox++)
      {
        const DataIndex din = edgeDit[mybox];

        int which = 0;

        for (IVSIterator ivsIt(edgeCells[din]); ivsIt.ok(); ++ivsIt, ++which)
          {
            if (m_graph[din].numVoFs(ivsIt()) != edgeBefore[din][which])
              {
                stale++;
              }
          }
      }

    stale = EBLevelDataOps::parallelSum(stale);

    if (stale > 0)
      {
        pout() << "    " << stale << " cells on the edge of what was coarsened changed underneath their neighbours"
               << endl;

        MayDay::Error("EBISLevel::coarsenFrom - coarsening changed cells on the edge of the region it covers, and the "
                      "cells outside still hold arcs describing them as they were. The region coarsened has to reach "
                      "past every cell coarsening changes.");
      }
  }
}

void EBISLevel::coarsenVoFs(EBISLevel&          a_fineEBIS,
                            LevelData<EBGraph>& a_fineGraph,
                            LevelData<EBData>&  a_fineData,
                            LevelData<EBGraph>& a_coarGraph,
                            const LayoutData<IntVectSet>& a_coarsenCells)
{
  CH_TIME("EBISLevel::coarsenVoFs");

  //so that i can do the vofs and faces of this level.
  DisjointBoxLayout fineFromCoarDBL;
  refine(fineFromCoarDBL, m_grids, 2);
  fineFromCoarDBL.close();

  //the fine graph is filled by the caller and shared with coarsenFaces
  LevelData<EBGraph>& fineFromCoarEBGraph = a_fineGraph;

  Interval interv(0,0);

  const DataIterator dit = m_grids.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) 
    {
      const DataIndex din = dit[mybox];

      if (a_coarsenCells[din].isEmpty())
        {
          continue;
        }

      const EBGraph& fineEBGraph = fineFromCoarEBGraph[din];
      const Box& coarRegion      = m_grids[din];
      EBGraph& coarEBGraph = m_graph[din];

      if (a_coarsenCells[din].contains(coarRegion))
        {
          coarEBGraph.coarsenVoFs(fineEBGraph, coarRegion);
        }
      else
        {
          coarEBGraph.coarsenVoFs(fineEBGraph, a_coarsenCells[din]);
        }
    }

  //the coarse ghost graph is filled here and reused by coarsenFaces, which needs the same layout,
  //the same ghost region and the same m_graph that this step has just written
  LevelData<EBGraph>& coarGhostEBGraph = a_coarGraph;

  if(s_distributedData){
    simplifyGraphFromGeo(coarGhostEBGraph, *m_geoserver, m_grids, m_domain, m_origin, m_dx);
  }

  m_graph.copyTo(interv, coarGhostEBGraph, interv);

  //dumpDebug(string("EBIS::coarsenVoFs"));

  //the fine data is defined and filled by the caller over the wider region coarsenFaces needs
  LevelData<EBData>& fineFromCoarEBData = a_fineData;

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) 
    {
      const DataIndex din = dit[mybox];

      if (a_coarsenCells[din].isEmpty())
        {
          continue;
        }

      const EBGraph& fineEBGraph =  fineFromCoarEBGraph[din];
      const EBData& fineEBData = fineFromCoarEBData[din];
      const EBGraph& coarEBGraph = coarGhostEBGraph[din];

      if (a_coarsenCells[din].contains(m_grids[din]))
        {
          m_data[din].coarsenVoFs(fineEBData, fineEBGraph, coarEBGraph, m_grids[din]);

          continue;
        }

      //the cells the finer level does not reach keep what they were generated with, but the data
      //has to be laid out again on the graph coarsening has just changed.  hold their moments
      //aside, lay it out, put them back, and let coarsening fill the rest
      IntVectSet keep = coarEBGraph.getIrregCells(m_grids[din]);
      keep -= a_coarsenCells[din];

      Vector<VolIndex> keptVoFs;
      Vector<VolData>  keptData;

      for (VoFIterator vofit(keep, coarEBGraph); vofit.ok(); ++vofit)
        {
          keptVoFs.push_back(vofit());
          keptData.push_back(m_data[din].getVolData()(vofit(), 0));
        }

      m_data[din].defineVoFData(coarEBGraph, m_grids[din]);

      for (int ikept = 0; ikept < keptVoFs.size(); ikept++)
        {
          m_data[din].getVolData()(keptVoFs[ikept], 0) = keptData[ikept];
        }

      IntVectSet fill = coarEBGraph.getIrregCells(m_grids[din]);
      fill &= a_coarsenCells[din];

      m_data[din].coarsenVoFs(fineEBData, fineEBGraph, coarEBGraph, fill);
    }
}

void EBISLevel::reconcileSeam(const LevelData<FArrayBox>& a_coarsenedMask)
{
  CH_TIME("EBISLevel::reconcileSeam");

  // What each cell holds for each of its own faces, so that a cell can be told what the cell
  // across a box boundary holds for the face they share.  Sending it as cell data rather than as
  // face data is what makes the answer unambiguous: a face belongs to two cells and both boxes
  // store it, so exchanging the faces themselves delivers whichever of the two the copy reached
  // last, which for half of them is the value we already had
  const int numSides = 2 * SpaceDim;

  LevelData<FArrayBox> faceView(m_grids, numSides, IntVect::Unit);

  const DataIterator& dit = m_grids.dataIterator();

  for (int mybox = 0; mybox < dit.size(); mybox++)
    {
      const DataIndex din = dit[mybox];

      const EBGraph& graph = m_graph[din];
      const EBData&  data  = m_data[din];

      FArrayBox& view = faceView[din];

      view.setVal(0.0);

      for (BoxIterator bit(m_grids[din]); bit.ok(); ++bit)
        {
          const IntVect iv = bit();

          if (graph.isCovered(iv))
            {
              continue;
            }

          if (graph.isRegular(iv))
            {
              for (int iside = 0; iside < numSides; iside++)
                {
                  view(iv, iside) = 1.0;
                }

              continue;
            }

          const Vector<VolIndex> vofs = graph.getVoFs(iv);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  Real area = 0.0;

                  for (int ivof = 0; ivof < vofs.size(); ivof++)
                    {
                      const Vector<FaceIndex> faces = graph.getFaces(vofs[ivof], idir, sit());

                      for (int iface = 0; iface < faces.size(); iface++)
                        {
                          area += data.areaFrac(faces[iface]);
                        }
                    }

                  view(iv, 2 * idir + ((sit() == Side::Hi) ? 1 : 0)) = area;
                }
            }
        }
    }

  faceView.exchange();

  long long ambiguous = 0;

  for (int mybox = 0; mybox < dit.size(); mybox++)
    {
      const DataIndex din = dit[mybox];

      const EBGraph&   graph = m_graph[din];
      const FArrayBox& view  = faceView[din];

      IntVectSet seam = graph.getIrregCells(m_grids[din]);

      for (IVSIterator trim(seam); trim.ok(); ++trim)
        {
          if (a_coarsenedMask[din](trim(), 0) > 0.5)
            {
              seam -= trim();
            }
        }

      for (IVSIterator ivsIt(seam); ivsIt.ok(); ++ivsIt)
        {
          const IntVect iv = ivsIt();

          bool touchesCoarsened = false;

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  if (a_coarsenedMask[din](iv + sign(sit()) * BASISV(idir), 0) > 0.5)
                    {
                      touchesCoarsened = true;
                    }
                }
            }

          if (!touchesCoarsened)
            {
              continue;
            }

          const Vector<VolIndex> vofs = graph.getVoFs(iv);

          // A cell holding more than one vof, or facing more than one across the seam, leaves no
          // way to say which of them the shared face belongs to.  Say how many rather than guess
          bool single = (vofs.size() == 1);

          for (int idir = 0; idir < SpaceDim && single; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  if (a_coarsenedMask[din](iv + sign(sit()) * BASISV(idir), 0) < 0.5)
                    {
                      continue;
                    }

                  if (graph.getFaces(vofs[0], idir, sit()).size() > 1)
                    {
                      single = false;
                    }
                }
            }

          if (!single)
            {
              ambiguous++;

              continue;
            }

          const VolIndex& vof = vofs[0];

          RealVect apertureVector = RealVect::Zero;

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  const IntVect other = iv + sign(sit()) * BASISV(idir);

                  const Vector<FaceIndex> faces = graph.getFaces(vof, idir, sit());

                  Real area = 0.0;

                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      FaceData& here = m_data[din].getFaceData(idir)(faces[iface], 0);

                      // take what the coarsened side holds, and keep it, so that whichever side a
                      // reader takes the face from it reads the same number
                      if (a_coarsenedMask[din](other, 0) > 0.5)
                        {
                          here.m_areaFrac = view(other, 2 * idir + ((sit() == Side::Hi) ? 0 : 1));
                        }

                      area += here.m_areaFrac;
                    }

                  apertureVector[idir] += (sit() == Side::Hi) ? area : -area;
                }
            }

          VolData& volData = m_data[din].getVolData()(vof, 0);

          const Real area = apertureVector.vectorLength();

          volData.m_averageFace.m_bndryArea = area;
          volData.m_averageFace.m_normal    = (area > 0.0) ? (apertureVector / area) : RealVect::Zero;
        }
    }

  ambiguous = EBLevelDataOps::parallelSum(ambiguous);

  if (ambiguous > 0)
    {
      pout() << "    " << ambiguous
             << " cells on the edge of what was coarsened face more than one vof and were left alone" << endl;
    }
}

void EBISLevel::fixFineToCoarse(EBISLevel& a_fineEBIS)
{
  CH_TIME("EBISLevel::fixFineToCoarse");
  // make a coarse layout from the fine layout so that we
  //can fix the fine->coarseVoF thing.
  DisjointBoxLayout coarFromFineDBL;

  coarsen(coarFromFineDBL, a_fineEBIS.m_graph.getBoxes(), 2);
  coarFromFineDBL.close();
  EBGraphFactory ebgraphfact(m_domain);
  LevelData<EBGraph> coarFromFineEBGraph(coarFromFineDBL,1, IntVect::Zero, ebgraphfact);
  Interval interv(0,0);

  if(s_distributedData){ // Won't work because of coarsening
    simplifyGraphFromGeo(coarFromFineEBGraph, *m_geoserver, coarFromFineDBL, a_fineEBIS.m_domain, m_origin, a_fineEBIS.m_dx);
  }
  m_graph.copyTo(interv, coarFromFineEBGraph, interv);

  for (DataIterator dit = a_fineEBIS.m_grids.dataIterator(); dit.ok(); ++dit)
    {
      EBGraph& fineEBGraph       =  a_fineEBIS.m_graph[dit()];
      const EBGraph& coarEBGraph = coarFromFineEBGraph[dit()];

      coarEBGraph.fixFineToCoarse(fineEBGraph);
    }

}

void EBISLevel::coarsenFaces(EBISLevel&          a_fineEBIS,
                             LevelData<EBGraph>& a_fineGraph,
                             LevelData<EBData>&  a_fineData,
                             LevelData<EBGraph>& a_coarGraph,
                             const LayoutData<IntVectSet>& a_coarsenCells)
{
  CH_TIME("EBISLevel::coarsenFaces");
  //now make a fine ebislayout with two ghost cell
  //on the same mapping as m_dbl
  //so that i can do the vofs and faces of this level.
  //this one will have the fine from coarse stuff fixed
  DisjointBoxLayout fineFromCoarDBL;
  refine(fineFromCoarDBL, m_grids, 2);
  fineFromCoarDBL.close();

  EBGraphFactory ebgraphfactcoar(m_domain);

  //both graphs were filled by coarsenVoFs.  the fine one is untouched since then and already
  //carries the three ghost cells this step needs, and the coarse one reflects the m_graph that
  //coarsenVoFs left behind, which is exactly what this step wants
  LevelData<EBGraph>& fineEBGraphGhostLD = a_fineGraph;
  LevelData<EBGraph>& coarEBGraphGhostLD = a_coarGraph;

  Interval interv(0,0);

  const DataIterator dit = m_grids.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) 
    {
      const DataIndex din = dit[mybox];  

      if (a_coarsenCells[din].isEmpty())
        {
          continue;
        }

      const EBGraph& fineEBGraphGhost = fineEBGraphGhostLD[din];
      const EBGraph& coarEBGraphGhost = coarEBGraphGhostLD[din];
      EBGraph& coarEBGraph = m_graph[din];

      if (a_coarsenCells[din].contains(m_grids[din]))
        {
          coarEBGraph.coarsenFaces(coarEBGraphGhost, fineEBGraphGhost);
        }
      else
        {
          coarEBGraph.coarsenFaces(coarEBGraphGhost, fineEBGraphGhost, a_coarsenCells[din]);
        }
    }
  //redefine coarebghostgraphld so i can use the faces for the ebdata
  coarEBGraphGhostLD.define(m_grids, 1,  IntVect::Unit, ebgraphfactcoar);

  if(s_distributedData){  // Should work because everything is on this level. 
    simplifyGraphFromGeo(coarEBGraphGhostLD, *m_geoserver, m_grids, m_domain, m_origin, m_dx);
  }
  m_graph.copyTo(interv, coarEBGraphGhostLD, interv);

  //the fine data was defined and filled by the caller over the region this step needs
  LevelData<EBData>& fineEBDataGhostLD = a_fineData;

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) 
    {
      const DataIndex din = dit[mybox];      

      if (a_coarsenCells[din].isEmpty())
        {
          continue;
        }

      const EBData&   fineEBData      = fineEBDataGhostLD[din];
      const EBGraph& fineEBGraphGhost = fineEBGraphGhostLD[din];
      const EBGraph& coarEBGraphGhost = coarEBGraphGhostLD[din];

      EBData& coarEBData   = m_data[din];

      if (a_coarsenCells[din].contains(m_grids[din]))
        {
          coarEBData.coarsenFaces(fineEBData, fineEBGraphGhost, coarEBGraphGhost, m_grids.get(din));

          continue;
        }

      //as in coarsenVoFs: the faces of the cells the finer level does not reach are held aside
      //while the data is laid out again.  only faces between two such cells are kept, since an
      //arc to a cell that was coarsened may have gone or arrived and coarsening fills those
      IntVectSet keep = coarEBGraphGhost.getIrregCells(m_grids[din]);
      keep -= a_coarsenCells[din];

      Vector<FaceIndex> keptFaces[SpaceDim];
      Vector<FaceData>  keptData[SpaceDim];

      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          for (FaceIterator faceit(keep, coarEBGraphGhost, faceDir, FaceStop::SurroundingWithBoundary);
               faceit.ok(); ++faceit)
            {
              const FaceIndex& face = faceit();

              if (a_coarsenCells[din].contains(face.gridIndex(Side::Lo)) ||
                  a_coarsenCells[din].contains(face.gridIndex(Side::Hi)))
                {
                  continue;
                }

              keptFaces[faceDir].push_back(face);
              keptData[faceDir].push_back(coarEBData.getFaceData(faceDir)(face, 0));
            }
        }

      coarEBData.defineFaceData(coarEBGraphGhost, m_grids.get(din));

      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          for (int ikept = 0; ikept < keptFaces[faceDir].size(); ikept++)
            {
              coarEBData.getFaceData(faceDir)(keptFaces[faceDir][ikept], 0) = keptData[faceDir][ikept];
            }
        }

      IntVectSet fill = coarEBGraphGhost.getIrregCells(m_grids[din]);
      fill &= a_coarsenCells[din];

      coarEBData.coarsenFaces(fineEBData, fineEBGraphGhost, coarEBGraphGhost, fill);
    }

}

EBISLevel::EBISLevel(EBISLevel             & a_fineEBIS,
                     const GeometryService & a_geoserver,
                     int                     a_nCellMax,
                     const bool            & a_fixRegularNextToMultiValued)
{ // method used by EBIndexSpace::buildNextLevel
  CH_TIME("EBISLevel::EBISLevel_fineEBIS");

  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;
  m_hasSurface  = false;

  m_domain = coarsen(a_fineEBIS.m_domain,2);
  m_dx = 2.*a_fineEBIS.m_dx;
  m_tolerance = 2.*a_fineEBIS.m_tolerance;
  m_origin = a_fineEBIS.m_origin;

  m_geoserver = &a_geoserver;

  m_level = a_fineEBIS.m_level + 1;

 
//  pout() << "before make boxes" << endl;
  if(!s_distributedData){
    Vector<Box> vbox;
    Vector<unsigned long long> irregCount;
    {
      CH_TIME("EBISLevel::EBISLevel_fineEBIS_makeboxes 2");
      makeBoxes(vbox, irregCount, m_domain.domainBox(), m_domain, a_geoserver,
		m_origin, m_dx, a_nCellMax);
    }

    //pout()<<vbox<<"\n\n";
    //load balance the boxes
    Vector<int> procAssign;
    //UnLongLongLoadBalance(procAssign, irregCount, vbox);
    basicLoadBalance(procAssign, vbox.size());
    //pout()<<procAssign<<std::endl;
    //define the layout.  this includes the domain and box stuff
    m_grids.define(vbox, procAssign);//this should use m_domain for periodic
  }
  else {
    (const_cast<GeometryService*>(&a_geoserver))->makeGrids(m_domain, m_grids, a_nCellMax, 15);
  }
  

  EBGraphFactory ebgraphfact(m_domain);
//  pout() << "before defining grids" << endl;
  m_graph.define(m_grids, 1, IntVect::Zero, ebgraphfact);
  
  EBDataFactory ebdatafact;
  m_data.define(m_grids, 1, IntVect::Zero, ebdatafact);
  
  coarsenFrom(a_fineEBIS, a_fixRegularNextToMultiValued);
#if 0
  pout() << "EBISLevel::EBISLevel 4 - m_grids - m_dx: " << m_dx << endl;
  pout() << "--------" << endl;
  pout() << m_grids.boxArray().size() << endl;
  pout() << "--------" << endl;
  pout() << endl;
#endif
}

EBISLevel::~EBISLevel()
{
}

void EBISLevel::fillEBISLayout(EBISLayout&              a_ebisLayout,
                               const DisjointBoxLayout& a_grids,
                               const int&               a_nghost) const
{
  CH_assert(a_nghost >= 0);
  
  //a_ebisLayout.define(m_domain, a_grids, a_nghost, m_graph, m_data);
  //return; // caching disabled for now.... ugh.  bvs

  EBISLayout& l = m_cache[a_grids];
  if (!l.isDefined() || (a_nghost > l.getGhost()))
    {
      CH_TIME("ebisllevel::fillebislayout cache miss");
      //int thisghost = Max(s_ebislGhost, a_nghost);
      int thisghost = a_nghost;
      l.define(m_domain, a_grids, thisghost, m_graph, m_data);
      m_cacheMisses++;
      m_cacheStale++;
      //pout()<<"a_nghost:"<<a_nghost;
    }
  else
    {
      CH_TIME("cache_hit");
      m_cacheHits++;
    }
  a_ebisLayout = l;// refcount is at least 2 now.
  if (m_cacheStale == 1)
    {
      refreshCache();
      m_cacheStale = 0;
    }
  //int isize = m_cache.size();
  //  pout()<<"fillebisl::m_level:"<<m_level<<" m_cache.size():"<<isize<<" m_cacheHits:"<<m_cacheHits<<" m_cacheMisses:"<<m_cacheMisses<<"\n";
}

void EBISLevel::dumpCache() const
{
  pout()<<std::endl;
  dmap::iterator d = m_cache.begin();
  while (d != m_cache.end())
    {
      d++;
          
    }

} 
void EBISLevel::refreshCache() const
{
  dmap::iterator d = m_cache.begin();
  while (d != m_cache.end())
    {
      if(d->second.refCount() ==1)
        {
          m_cache.erase(d++);
        }
      else
        {
          d++;
        }
    }

 // int s=m_cache.size();
 // pout()<<" m_level:"<<m_level<<" m_cache.size():"<<s<<" m_cacheHits:"<<m_cacheHits<<" m_cacheMisses:"<<m_cacheMisses<<"\n";
}

void EBISLevel::clearMultiBoundaries()
{
  CH_TIME("EBISLevel::clearMultiBoundaries");
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      EBData& data = m_data[dit];
      data.clearMultiBoundaries();
    }
}

void EBISLevel::setBoundaryPhase(int phase)
{
  CH_TIME("EBISLevel::setBoundaryPhase");
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      EBData& data = m_data[dit];
      data.setBoundaryPhase(phase);
    }
}

//throughout this routine A refers to objects belonging to this
//ebindexspace,  B refers to the other one.
void setOtherVoFSingleValued(VolData&        a_vol,
                             const VolIndex& a_sourceVoF,
                             int&            a_otherPhase,
                             EBData&         a_otherFluidData,
                             bool&           a_sourceIvContainedInOtherFluid,
                             bool&           a_skipCell)
{
  VolIndex vother = a_sourceVoF;
  a_skipCell = false;
  bool fixBoundaryData = false;
  VolData refVol;
  //if a cell is full and has a unit boundary area, it means that the irregular
  //boundary lives on a cell boundary.  We have to link to a cell over.
  // In this case, in the presence of multi valued cells, to quote Yeats:
  //Mere anarchy is loosed upon the world,
  //The blood-dimmed tide is loosed, and everywhere
  //The ceremony of innocence is drowned.
  if ((a_vol.m_volFrac == 1) && (a_vol.m_averageFace.m_bndryArea == 1))
    {
      a_skipCell = true;
      //figure out which cell we are moving to. Then we use cellIndex = 0
      //because the cells involved are single valued.
      int whichWayToGo = -23;
      int sign = 1;
      bool found = false;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real normDir = a_vol.m_averageFace.m_normal[idir];
          if (Abs(normDir) == 1)
            {
              found = true;
              whichWayToGo = idir;
              sign = -1;
              if (normDir < 0)
                {
                  sign = 1;
                }
            }
        }
      if (!found)
        {
          MayDay::Error("EBIndexSpace::setOtherVoFSingleValued - Unit normal direction not found");
        }
      IntVect ivOther = sign*BASISV(whichWayToGo);
      ivOther += a_sourceVoF.gridIndex();
      vother = VolIndex(ivOther, 0);
    }
  else if ((a_vol.m_volFrac == 0) &&
           (a_vol.m_averageFace.m_bndryArea == 1) &&
           (a_sourceIvContainedInOtherFluid))
    {
      a_skipCell = true;
      refVol = a_otherFluidData.getVolData()(a_sourceVoF, 0);
      // use -1*normal of the VoF with this IV in the other fluid
      RealVect normal = -refVol.m_averageFace.m_normal;
      //figure out which cell we are moving to. Then we use cellIndex = 0
      //because the cells involved are single valued.
      int whichWayToGo = -23;
      int sign = 1;
      bool found = false;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real normDir = normal[idir];
          if (Abs(normDir) == 1)
            {
              found = true;
              whichWayToGo = idir;
              sign = 1;
              if (normDir < 0)
                {
                  sign = -1;
                }
            }
        }
      if (found)
        {
          IntVect ivOther = sign*BASISV(whichWayToGo);
          ivOther += a_sourceVoF.gridIndex();
          vother = VolIndex(ivOther, 0);
          fixBoundaryData = true;
        }
      else
        {
          MayDay::Warning("EBIndexSpace::setOtherVoFSingleValued - Unit normal direction not found; EB probably intersects cell corner");
        }
    }
  a_vol.m_phaseFaces[0].m_volIndex   = vother;
  a_vol.m_phaseFaces[0].m_bndryPhase = a_otherPhase;
  if (fixBoundaryData)
    {
      Real eps = 1e-6;
      BoundaryData&    boundaryData0 = a_vol.m_phaseFaces[0];
      BoundaryData& avgBoundaryData0 = a_vol.m_averageFace;
      BoundaryData&    boundaryData1 = refVol.m_averageFace;
      if ((boundaryData0.m_normal.vectorLength() < eps) &&
          (boundaryData1.m_normal.vectorLength() > 1-eps) &&
          (boundaryData1.m_normal.vectorLength() < 1+eps))
        {
          // fix boundary at phaseFace
          boundaryData0.m_normal = -1.0*(boundaryData1.m_normal);
          boundaryData0.m_bndryArea = boundaryData1.m_bndryArea;
          boundaryData0.m_bndryCentroid = boundaryData1.m_bndryCentroid;
          // make average face reflect new phaseFace data
          avgBoundaryData0.m_normal = boundaryData0.m_normal;
          avgBoundaryData0.m_bndryArea = boundaryData0.m_bndryArea;
          avgBoundaryData0.m_bndryCentroid = boundaryData0.m_bndryCentroid;
        }
    }
}

void testAndFixBoundaryData(VolData&       a_volData0,
                            const VolData& a_volData1)
{
  Real eps = 1e-6;
  BoundaryData&    boundaryData0 = a_volData0.m_phaseFaces[0];
  BoundaryData& avgBoundaryData0 = a_volData0.m_averageFace;
  const BoundaryData&    boundaryData1 = a_volData1.m_phaseFaces[0];
  //  const BoundaryData& avgBoundaryData1 = a_volData1.m_averageFace;
  if ((boundaryData0.m_bndryArea > 0) &&
      (boundaryData0.m_normal.vectorLength() < eps) &&
      (boundaryData1.m_normal.vectorLength() > 1-eps) &&
      (boundaryData1.m_normal.vectorLength() < 1+eps))
    {
      // fix boundary at phaseFace
      boundaryData0.m_normal = -1.0*(boundaryData1.m_normal);
      boundaryData0.m_bndryArea = boundaryData1.m_bndryArea;
      boundaryData0.m_bndryCentroid = boundaryData1.m_bndryCentroid;
      // make average face reflect new phaseFace data
      avgBoundaryData0.m_normal = boundaryData0.m_normal;
      avgBoundaryData0.m_bndryArea = boundaryData0.m_bndryArea;
      avgBoundaryData0.m_bndryCentroid = boundaryData0.m_bndryCentroid;
    }
}

void EBISLevel::reconcileIrreg(EBISLevel& a_otherPhase)
{
  CH_TIME("EBISLevel::reconcileIrreg");

  //both layouts are made independently so will not work
  //with the same data iterator but they ought to because they
  //are both breaking up the domain the same way.   This is yet
  //another subtle thing thing that will break when the
  //number of fluids is greater than two.

  //define ghosted coarse graph so i can use the faces for the ebdata
  EBGraphFactory ebgraphfactcoarA(m_domain);
  LevelData<EBGraph> ghostedEBGraphLDCoarA(m_grids, 1, IntVect::Unit,
                                           ebgraphfactcoarA);
  Interval interv(0, 0);
  m_graph.copyTo(interv, ghostedEBGraphLDCoarA, interv);

  EBGraphFactory ebgraphfactcoarB(a_otherPhase.m_domain);
  LevelData<EBGraph> ghostedEBGraphLDCoarB(a_otherPhase.m_grids,
                                           1,
                                           IntVect::Unit,
                                           ebgraphfactcoarB);
  a_otherPhase.m_graph.copyTo(interv, ghostedEBGraphLDCoarB, interv);

  DataIterator dita = m_grids.dataIterator();
  DataIterator ditb= a_otherPhase.m_grids.dataIterator();
  dita.begin(); ditb.begin();
  for ( ; dita.ok(); ++dita, ++ditb)
    {
      EBGraph& ebgrapCoarA =              m_graph[dita];
      EBGraph& ebgrapCoarB = a_otherPhase.m_graph[ditb];
      EBData&  ebdataCoarA =              m_data[dita];
      EBData&  ebdataCoarB = a_otherPhase.m_data[ditb];
      Box region = ebgrapCoarA.getRegion();

      // ghosted graphs for filling faces for ebdata
      EBGraph& ghostedEBGraphCoarA = ghostedEBGraphLDCoarA[dita];
      EBGraph& ghostedEBGraphCoarB = ghostedEBGraphLDCoarB[ditb];

      CH_assert(ebgrapCoarB.getRegion() == region);
      IntVectSet seta = ebgrapCoarA.getIrregCells(region);
      IntVectSet setb = ebgrapCoarB.getIrregCells(region);

      IntVectSet abDifference = seta - setb;
      IntVectSet baDifference = setb - seta;
      //this should only happen when one side changed a regular to an
      //irregular.   This means that the other has to change a covered
      //to an irregular
      ebgrapCoarB.addEmptyIrregularVoFs(abDifference);
      ebgrapCoarA.addEmptyIrregularVoFs(baDifference);
      // repeat add for ghosted graphs
      ghostedEBGraphCoarB.addEmptyIrregularVoFs(abDifference);
      ghostedEBGraphCoarA.addEmptyIrregularVoFs(baDifference);
      // use ghosted graphs to add to ebdata
      ebdataCoarB.addEmptyIrregularVoFs(abDifference, ghostedEBGraphCoarB);
      ebdataCoarA.addEmptyIrregularVoFs(baDifference, ghostedEBGraphCoarA);
    } //end loop over grids
}

void EBISLevel::levelStitch(EBISLevel&       a_otherPhase,
                            const EBISLevel* a_finePtrA,
                            const EBISLevel* a_finePtrB)
{
  CH_TIME("EBISLevel::levelStitch");

  //I have no idea what to do if only one of the inputs is null.
  //either both or neither makes sense
  CH_assert(((a_finePtrA != NULL) && (a_finePtrB != NULL)) ||
            ((a_finePtrA == NULL) && (a_finePtrB == NULL)));
  EBISLayout ebislFineA, ebislFineB;
  DisjointBoxLayout dblFineA, dblFineB;
  if (a_finePtrA != NULL)
    {
      int nghost = 0; //should not need any as this is all about connectivity within a coarse cell
      refine(dblFineA,              m_grids, 2);
      refine(dblFineB, a_otherPhase.m_grids, 2);
      a_finePtrA->fillEBISLayout(ebislFineA, dblFineA, nghost);
      a_finePtrB->fillEBISLayout(ebislFineB, dblFineB, nghost);
    }

  //both layouts are made independently so will not work
  //with the same data iterator but they ought to because they
  //are both breaking up the domain the same way.   This is yet
  //another subtle thing thing that will break when the
  //number of fluids is greater than two.
  DataIterator dita = m_grids.dataIterator();
  DataIterator ditb= a_otherPhase.m_grids.dataIterator();
  dita.begin(); ditb.begin();
  for ( ; dita.ok(); ++dita, ++ditb)
    {
      EBGraph& ebgrapCoarA =              m_graph[dita];
      EBGraph& ebgrapCoarB = a_otherPhase.m_graph[ditb];
      EBData&  ebdataCoarA =              m_data[dita];
      EBData&  ebdataCoarB = a_otherPhase.m_data[ditb];
      // Box region = ebgrapCoarA.getRegion();
      Box region = m_grids[dita];

      // CH_assert(ebgrapCoarB.getRegion() == region);
      CH_assert(a_otherPhase.m_grids[ditb] == region);
      IntVectSet seta = ebgrapCoarA.getIrregCells(region);
      IntVectSet setb = ebgrapCoarB.getIrregCells(region);

      //different EBIS's can (correctly) disagree about which cells are multivalued
      IntVectSet setMultiA = ebgrapCoarA.getMultiCells(region);
      IntVectSet setMultiB = ebgrapCoarB.getMultiCells(region);
      IntVectSet setMulti = setMultiA | setMultiB;
      int aphase =              m_phase;
      int bphase = a_otherPhase.m_phase;

      IntVectSet setANoMulti = seta - setMulti;
      IntVectSet setBNoMulti = setb - setMulti;
      {
        //EBIndexSpace does the right thing with all the geometric information
        //in the case of two single valued vofs.   All that needs to be set is
        //the phase on the other side of the irregular face.
        for (IVSIterator ita(setANoMulti); ita.ok(); ++ita)
          {
            VolIndex v(ita(), 0); //0 because we know this whole set is single valued

            VolData& volDatA = ebdataCoarA.getVolData()(v,0);
            volDatA.m_phaseFaces.resize(1);
            //this sets boundary area and normal and centroid to
            //whatever the ebindexspace set it to.
            volDatA.m_phaseFaces[0]=volDatA.m_averageFace;
            bool skip = false;
            bool vContainedInOtherIVS = setBNoMulti.contains(v.gridIndex());
            setOtherVoFSingleValued(volDatA, v, bphase, ebdataCoarB, vContainedInOtherIVS, skip);
            if (skip && setMultiB.contains(volDatA.m_phaseFaces[0].m_volIndex.gridIndex()))
              {
                MayDay::Error("coordinate face boundary also a multi valued one");
              }
          }//end loop over cells  where both phases are singlevalued

        for (IVSIterator itb(setBNoMulti); itb.ok(); ++itb)
          {
            VolIndex v(itb(), 0); //0 because we know this whole set is single valued

            VolData& volDatB = ebdataCoarB.getVolData()(v,0);
            volDatB.m_phaseFaces.resize(1);
            //this sets boundary area and normal and centroid to
            //whatever the ebindexspace set it to.
            volDatB.m_phaseFaces[0]=volDatB.m_averageFace;
            bool skip = false;
            bool vContainedInOtherIVS = setANoMulti.contains(v.gridIndex());
            setOtherVoFSingleValued(volDatB, v, aphase, ebdataCoarA, vContainedInOtherIVS, skip);
            if (skip && setMultiA.contains(volDatB.m_phaseFaces[0].m_volIndex.gridIndex()))
              {
                MayDay::Error("coordinate face boundary also a multi valued one");
              }
          }//end loop over cells where both phases are singlevalued

        // now fix incorrect boundary data, where possible
        // first, for the phase a volData
//         for (IVSIterator ita(setANoMulti); ita.ok(); ++ita)
        // hack to get around mismatching sets
        IntVectSet abIntersectionNoMulti = setANoMulti & setBNoMulti;
        for (IVSIterator ita(abIntersectionNoMulti); ita.ok(); ++ita)
          {
            VolIndex v(ita(), 0);

            VolData& volDatA = ebdataCoarA.getVolData()(v,0);
            VolIndex vother = volDatA.m_phaseFaces[0].m_volIndex;
            // need to fix the normal
            VolData& volDatB = ebdataCoarB.getVolData()(vother,0);
            // only send the first boundary data, since single valued cell
            testAndFixBoundaryData(volDatA, volDatB);
          }
        // next, for the phase b volData
//         for (IVSIterator itb(setBNoMulti); itb.ok(); ++itb)
        // hack to get around mismatching sets
        for (IVSIterator itb(abIntersectionNoMulti); itb.ok(); ++itb)
          {
            VolIndex v(itb(), 0); //0 because we know this whole set is single valued

            VolData& volDatB = ebdataCoarB.getVolData()(v,0);
            VolIndex vother = volDatB.m_phaseFaces[0].m_volIndex;
            // need to fix the normal
            VolData& volDatA = ebdataCoarA.getVolData()(vother,0);
            // only send the first boundary data, since single valued cell
            testAndFixBoundaryData(volDatB, volDatA);
          }
      }
      if (!setMulti.isEmpty())
        {
          //now have to do the union of the each set
          //of multivalued cells.   either phase being multivalued can
          //make either ebindex space get the wrong answer for
          //the geometric information
          const EBISBox& ebisbxFineA = ebislFineA[dita()];
          const EBISBox& ebisbxFineB = ebislFineB[ditb()];

          IVSIterator it(setMulti); //see above derivation.
          for (it.begin(); it.ok(); ++it)
            {
              cellStitch(ebdataCoarA, ebdataCoarB,
                         ebgrapCoarA, ebgrapCoarB,
                         ebisbxFineA, ebisbxFineB,
                         it(), aphase, bphase);
            }//end loop over multivalued cells
        }
    } //end loop over grids
}

void getFinerBoundaryData(Vector<Vector<BoundaryData> > & a_bndryDataFineA,
                          Vector<VolIndex>&               a_otherCellIndex,
                          const VolIndex&                 a_coarVoFA,
                          const EBGraph&                  a_ebgrapCoarA,
                          const EBGraph&                  a_ebgrapCoarB,
                          const EBISBox&                  a_ebisbxFineA)
{
  //divdes fine boundary datas a by which coarse vof b they are connected to.
  // to which they are connected?

  const IntVect ivCoar = a_coarVoFA.gridIndex();

  Vector<VolIndex>          coarVoFsB = a_ebgrapCoarB.getVoFs(ivCoar);
  Vector<Vector<VolIndex> > fineVoFsB(coarVoFsB.size());
  a_bndryDataFineA.resize(coarVoFsB.size());
  a_otherCellIndex.resize(coarVoFsB.size());
  for (int ib = 0; ib < coarVoFsB.size(); ib++)
    {
      fineVoFsB[ib] = a_ebgrapCoarB.refine(coarVoFsB[ib]);
      a_otherCellIndex[ib] = coarVoFsB[ib];
    }

  Vector<VolIndex> allFinerVoFsA = a_ebgrapCoarA.refine(a_coarVoFA);
  for (int ib = 0; ib < coarVoFsB.size(); ib++)
    {
      a_bndryDataFineA[ib].resize(0);
      for (int ifineB = 0; ifineB < fineVoFsB[ib].size(); ifineB++)
        {
          const VolIndex& vofFineB = fineVoFsB[ib][ifineB];
          for (int ifineA = 0; ifineA < allFinerVoFsA.size(); ifineA++)
            {
              const VolIndex& vofFineA = allFinerVoFsA[ifineA];
              if (vofFineA.gridIndex() == vofFineB.gridIndex())
                {
                  if (a_ebisbxFineA.isIrregular(vofFineA.gridIndex()))
                    {
                      const VolData& voldat = a_ebisbxFineA.getEBData().getVolData()(vofFineA, 0);
                      const Vector<BoundaryData>& boundaryDat = voldat.m_phaseFaces;
                      for (int iface = 0; iface < boundaryDat.size(); iface++)
                        {
                          if (boundaryDat[iface].m_volIndex.cellIndex() == vofFineB.cellIndex())
                            {
                              a_bndryDataFineA[ib].push_back(boundaryDat[iface]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void coarsenBoundaryInfo(EBData&                                a_ebdataCoar,
                         const EBISBox&                         a_ebisbxFine,
                         const Vector< Vector<BoundaryData > >& a_finerBoundaryData,
                         const VolIndex&                        a_vof,
                         const int&                             a_otherPhase,
                         const Vector<VolIndex>&                a_otherCellIndex)
{
  //this coarsening factor for boundary area fractions can be thought of this way
  //(1) take fine area and multiply by dxf^{SD-1}
  //(2) add up areas (now we have the real coarse area).
  //(3) divide by dxc^{SD-1} = dxf{SD-1}*2{SD-1}
  Real faceCoarsenFactor = D_TERM(1.0, * 2.0,* 2.0);

  //the point of this routine is to fix the boundary area and normal
  //stuff in this object
  VolData& coarData = a_ebdataCoar.getVolData()(a_vof, 0);
  coarData.m_phaseFaces.resize(a_finerBoundaryData.size());
  for (int ib = 0; ib < a_finerBoundaryData.size(); ib++)
    {
      //initialization is important  in both of these things.
      //since i am taking the average
      RealVect aveNormal  = RealVect::Zero;
      RealVect aveCentro  = RealVect::Zero;
      Real sumArea=0;
      Real aveArea=0;
      int numfaces=0;
      for (int jb = 0; jb < a_finerBoundaryData[ib].size(); jb++)
        {
          const BoundaryData& boundaryDat = a_finerBoundaryData[ib][jb];
          const Real    & fineArea = boundaryDat.m_bndryArea;
          const RealVect& fineNorm = boundaryDat.m_normal;
          const RealVect& fineCent = boundaryDat.m_bndryCentroid;
          sumArea += fineArea;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              aveNormal[idir] += fineArea*fineNorm[idir];
              aveCentro[idir] += fineArea*fineCent[idir];
            }
          numfaces++;
        }

      //we are taking an average-weighted average
      //of each normal and centroid components so we need to divide by the area
      if (sumArea > 1.0e-12)
        {
          aveNormal /= sumArea;
          aveCentro /= sumArea;
        }
      Real sumSq;
      PolyGeom::unifyVector(aveNormal, sumSq);
      //this takes into account the fact that boundary areas are normalized
      //by dx^(SD-1).  See note above.
      aveArea  = sumArea/faceCoarsenFactor;

      coarData.m_phaseFaces[ib].m_normal            = aveNormal;
      coarData.m_phaseFaces[ib].m_bndryCentroid     = aveCentro;
      coarData.m_phaseFaces[ib].m_bndryArea         = aveArea;
      coarData.m_phaseFaces[ib].m_bndryPhase        = a_otherPhase;
      coarData.m_phaseFaces[ib].m_volIndex          = a_otherCellIndex[ib];
    }

}

void EBISLevel::cellStitch(EBData&        a_ebdataCoarA,
                           EBData&        a_ebdataCoarB,
                           const EBGraph& a_ebgrapCoarA,
                           const EBGraph& a_ebgrapCoarB,
                           const EBISBox& a_ebisbxFineA,
                           const EBISBox& a_ebisbxFineB,
                           const IntVect& a_iv,
                           const int&     a_phaseA,
                           const int&     a_phaseB)
{
  //pout() << "entering cellStitch" << endl;
  Vector<VolIndex> coarVoFsA = a_ebgrapCoarA.getVoFs(a_iv);
  Vector<VolIndex> coarVoFsB = a_ebgrapCoarB.getVoFs(a_iv);
  for (int ivof = 0; ivof < coarVoFsA.size(); ivof++)
    {
      const VolIndex& vof = coarVoFsA[ivof];
      // each element in finerVoFs is a list of finer vofs whose
      // opposing number on the other side of the irregular boundary
      // are connected to each other.   This will correspond with
      // the Vector<VolIndex> of the coarse vofs in the *other* phase.
      Vector<Vector< BoundaryData> > finerBndryData;
      Vector<VolIndex>  cellIndexB;
      getFinerBoundaryData(finerBndryData, cellIndexB, vof, a_ebgrapCoarA, a_ebgrapCoarB, a_ebisbxFineA);


      //once the above list is made, we can correctly coarsen the boundary
      //information using only the information in this phase.
      coarsenBoundaryInfo(a_ebdataCoarA, a_ebisbxFineA, finerBndryData, vof, a_phaseB, cellIndexB);
    }

  for (int ivof = 0; ivof < coarVoFsB.size(); ivof++)
    {
      const VolIndex& vof = coarVoFsB[ivof];
      // each element in finerVoFs is a list of finer vofs whose
      // opposing number on the other side of the irregular boundary
      // are connected to each other.   This will correspond with
      // the Vector<VolIndex> of the coarse vofs in the *other* phase.
      Vector<Vector< BoundaryData> > finerBndryData;
      Vector<VolIndex>  cellIndexA;
      getFinerBoundaryData(finerBndryData, cellIndexA, vof, a_ebgrapCoarB, a_ebgrapCoarA, a_ebisbxFineB);


      //once the above list is made, we can correctly coarsen the boundary
      //information using only the information in this phase.
      coarsenBoundaryInfo(a_ebdataCoarB, a_ebisbxFineB, finerBndryData, vof, a_phaseA, cellIndexA);
    }
}

#ifdef CH_USE_HDF5
void EBISLevel::write(HDF5Handle& a_handle,
                      const int& a_levelNumber) const
{
  CH_TIME("EBISLevel::write_with_lev_number");
  HDF5HeaderData header;
  char levelcstr[256];
  sprintf(levelcstr, "%d", a_levelNumber);
  string levelstring(levelcstr);

  header.m_realvect["EBIS_origin"] = m_origin;

  string domstring = string("EBIS_domain_lev_") + levelstring;
  header.m_box[domstring] = m_domain.domainBox();
  
  string dxstring = string("EBIS_dx_lev_") + levelstring;
  header.m_real[dxstring]     = m_dx;

  //  pout() << "writing handle  for level " << a_levelNumber << endl;
  header.writeToFile(a_handle);

  //write the grids to the file
  string boxstring = string("EBIS_boxes_lev_") + levelstring;
  CH_XD::write(a_handle, m_grids, boxstring);

  string graphstring = string("EBIS_graph_lev_") + levelstring;
  string datastring = string("EBIS_data_lev_") + levelstring;

  //  pout() << "writing graph for level " << a_levelNumber << endl;
  int eekflag = CH_XD::write(a_handle, m_graph, graphstring);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }
  //  pout() << "writing data  for level " << a_levelNumber << endl;
  eekflag = CH_XD::write(a_handle, m_data ,  datastring);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
}

EBISLevel::EBISLevel(HDF5Handle& a_handle, 
                     const int& a_levelNumber)
{
  CH_TIME("EBISLevel::EBISLevel_hdf5_with_lev_number");

  m_cacheMisses = 0;
  m_cacheHits   = 0;
  m_cacheStale  = 0;
  m_hasSurface  = false;
  char levelcstr[256];
  sprintf(levelcstr, "%d", a_levelNumber);
  string levelstring(levelcstr);

  HDF5HeaderData header;
  header.readFromFile(a_handle);
  m_origin = header.m_realvect["EBIS_origin"];

  string domstring = string("EBIS_domain_lev_") + levelstring;
  m_domain = header.m_box[domstring];

  string dxstring = string("EBIS_dx_lev_") + levelstring;
  m_dx =     header.m_real[dxstring]    ;

  m_tolerance = m_dx*1E-4;

  //read in the grids
  Vector<Box> boxes;
  string boxstring = string("EBIS_boxes_lev_") + levelstring;
  read(a_handle,boxes, boxstring);
  Vector<int> procAssign;
  LoadBalance(procAssign, boxes);
  m_grids.define(boxes, procAssign);//this should use m_domain for periodic...
  EBGraphFactory graphfact(m_domain);
  m_graph.define(m_grids, 1, IntVect::Zero, graphfact);


  string graphstring = string("EBIS_graph_lev_") + levelstring;
  int eekflag = read(a_handle, m_graph, graphstring, m_grids, Interval(), false);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }

  //need a ghosted layout so that the data can be defined properly
  LevelData<EBGraph> ghostGraph(m_grids, 1, IntVect::Unit, graphfact);
  Interval interv(0,0);
  m_graph.copyTo(interv, ghostGraph, interv);

  //now the data for the graph
  EBDataFactory dataFact;
  m_data.define(m_grids, 1, IntVect::Zero, dataFact);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_data[dit()].defineVoFData(ghostGraph[dit()], m_grids.get(dit()));
      m_data[dit()].defineFaceData(ghostGraph[dit()], m_grids.get(dit()));
    }
  //read the data  in from the file
  string  datastring  = string("EBIS_data_lev_") + levelstring;
  eekflag = read(a_handle, m_data ,  datastring, m_grids, Interval(), false);
  if (eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
#if 0
  pout() << "EBISLevel::EBISLevel 5 - m_grids - m_dx: " << m_dx << endl;
  pout() << "--------" << endl;
  pout() << m_grids.boxArray().size() << endl;
  pout() << "--------" << endl;
  pout() << endl;
#endif
}
#endif

const ProblemDomain& EBISLevel::getDomain() const
{
  return m_domain;
}

const Real& EBISLevel::getDX() const
{
  return m_dx;
}

const RealVect& EBISLevel::getOrigin() const
{
  return m_origin;
}

const DisjointBoxLayout& EBISLevel::getGrids() const
{
  return m_grids;
}

// All boxes containing an irregular cell
DisjointBoxLayout EBISLevel::getIrregGrids(const ProblemDomain& a_domain) const
{
  DisjointBoxLayout irregGrids;
  DataIterator dit =   m_graph.dataIterator();
  Vector<Box> localBoxes;

  for (dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if (graph.hasIrregular())
        {
          localBoxes.push_back(b);

        }
    }
  Vector<Vector<Box> > allBoxes;

  gather(allBoxes, localBoxes, uniqueProc(SerialTask::compute));

  broadcast(allBoxes, uniqueProc(SerialTask::compute));

  Vector<Box> boxes;
  for (int i = 0; i < allBoxes.size(); i++)
    {
      boxes.append(allBoxes[i]);
    }
  Vector<int> procs;
  LoadBalance(procs, boxes);
  irregGrids = DisjointBoxLayout(boxes, procs, a_domain);

  return irregGrids;
}

// All boxes with only irregular or regular cells
DisjointBoxLayout EBISLevel::getFlowGrids(const ProblemDomain& a_domain) const
{
  DisjointBoxLayout flowGrids;
  DataIterator dit =   m_graph.dataIterator();
  Vector<Box> localBoxes;

  for (dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if (graph.hasIrregular() || graph.isAllRegular())
        {
          localBoxes.push_back(b);

        }
    }
  Vector<Vector<Box> > allBoxes;

  gather(allBoxes, localBoxes, uniqueProc(SerialTask::compute));

  broadcast(allBoxes, uniqueProc(SerialTask::compute));
  Vector<Box> boxes;
  for (int i = 0; i < allBoxes.size(); i++)
    {
      boxes.append(allBoxes[i]);
    }
  Vector<int> procs;
  LoadBalance(procs, boxes);

  flowGrids = DisjointBoxLayout(boxes, procs, a_domain);
  return flowGrids;
}

// All boxes with an irregular or covered cell
DisjointBoxLayout EBISLevel::getCoveredGrids(const ProblemDomain& a_domain) const
{
  DisjointBoxLayout coveredGrids;
  DataIterator dit = m_graph.dataIterator();
  Vector<Box> localBoxes;

  for (dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if (graph.hasIrregular() || graph.isAllCovered())
        {
          localBoxes.push_back(b);
        }
    }
  Vector<Vector<Box> > allBoxes;

  gather(allBoxes, localBoxes, uniqueProc(SerialTask::compute));

  broadcast(allBoxes, uniqueProc(SerialTask::compute));

  Vector<Box> boxes;
  for (int i = 0; i < allBoxes.size(); i++)
    {
      boxes.append(allBoxes[i]);
    }
  Vector<int> procs;
  LoadBalance(procs, boxes);
  coveredGrids = DisjointBoxLayout(boxes, procs, a_domain);

  return coveredGrids;
}

IntVectSet EBISLevel::irregCells() const
{
  DataIterator dit =   m_graph.dataIterator();
  IntVectSet rtn;

  for (dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if (graph.hasIrregular())
        {
          rtn |= graph.getIrregCells(b);
        }
    }
  return rtn;
}

void EBISLevel::getGraphSummary(long long & a_irrVoFs,
                                long long & a_arcs,
                                long long & a_multiVoFs,
                                long long & a_zeroVoFs,
                                long long & a_zeroVoFsArcs)
{
  CH_TIME("EBISLevel::getGraphSummary");

  a_irrVoFs      = 0;
  a_arcs         = 0;
  a_multiVoFs    = 0;
  a_zeroVoFs     = 0;
  a_zeroVoFsArcs = 0;

  for (DataIterator dit(m_grids); dit.ok(); ++dit)
  {
    const Box & curBox = m_grids[dit()];

    EBGraph & curGraph = m_graph[dit()];
    const EBData  & curData  = m_data[dit()];

    if (curGraph.m_implem->m_tag == EBGraphImplem::HasIrregular)
    {

      BaseFab<GraphNode> & curGraphFAB  = curGraph.m_implem->m_graph;
      const BaseIVFAB<VolData> & curDataIVFAB = curData.m_implem->m_volData;

      for (BoxIterator bit(curBox); bit.ok(); ++bit)
      {
        const IntVect   & curIV = bit();
        GraphNode & curGraphNode = curGraphFAB(curIV,0);

        if (curGraphNode.isIrregular())
        {
          int vofs = curGraphNode.m_cellList->size();

          a_irrVoFs += vofs;

          if (vofs > 1)
          {
            a_multiVoFs += vofs;
          }

          for (int vof=0; vof < vofs; vof++)
          {
            VolIndex volIndex(curIV,vof);

            Real volFrac = curDataIVFAB(volIndex,0).m_volFrac;

            if (volFrac == 0.0)
            {
              a_zeroVoFs++;
            }

            const GraphNodeImplem & curGraphNodeImplem = (*curGraphNode.m_cellList)[vof];
            int localVoFArcs = 0;

            for (int iside = 0; iside < 2*SpaceDim; iside++)
            {
              int arcs = curGraphNodeImplem.m_arc[iside].size();
              localVoFArcs += arcs;

            }

            if (volFrac == 0.0)
            {
              a_zeroVoFsArcs += localVoFArcs;

              if (vofs == 1 && localVoFArcs == 0)
              {
                delete curGraphNode.m_cellList;
                curGraphNode.m_cellList = 0;
                a_irrVoFs--;
              }
              else
              {
                a_arcs += localVoFArcs;
              }
            }
            else
            {
              a_arcs += localVoFArcs;
            }
          }
        }
      }
    }
  }

  a_irrVoFs      = EBLevelDataOps::parallelSum(a_irrVoFs);
  a_arcs         = EBLevelDataOps::parallelSum(a_arcs);
  a_multiVoFs    = EBLevelDataOps::parallelSum(a_multiVoFs);
  a_zeroVoFs     = EBLevelDataOps::parallelSum(a_zeroVoFs);
  a_zeroVoFsArcs = EBLevelDataOps::parallelSum(a_zeroVoFsArcs);
}

void EBISLevel::printGraphSummary(char const * a_prefix)
{
  CH_TIME("EBISLevel::printGraphSummary");

  long long irrVoFs;
  long long arcs;
  long long multiVoFs;
  long long zeroVoFs;
  long long zeroVoFsArcs;

  getGraphSummary(irrVoFs,arcs,multiVoFs,zeroVoFs,zeroVoFsArcs);

  pout() << a_prefix << "irreg VoFs    : " << irrVoFs      << endl;
  pout() << a_prefix << "arcs          : " << arcs         << endl;
  pout() << a_prefix << "multi VoFs    : " << multiVoFs    << endl;
  pout() << a_prefix << "zero VoFs     : " << zeroVoFs     << endl;
  pout() << a_prefix << "zero VoFs arcs: " << zeroVoFsArcs << endl;
}

#include "NamespaceFooter.H"
