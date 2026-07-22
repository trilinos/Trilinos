// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

// STK_search objects
#include "stk_search/CoarseSearch.hpp"

// Copied from PeriodicBoundarySearch
struct SelectorMapping {
  SelectorMapping(const stk::mesh::Selector & domain, const stk::mesh::Selector & range) : m_domain(domain), m_range(range) {}
  stk::mesh::Selector m_domain;
  stk::mesh::Selector m_range;
  stk::mesh::Selector m_unique_range;
};

typedef double Scalar;
typedef std::vector<SelectorMapping> SelectorMappingVector;
typedef stk::search::IdentProc<stk::mesh::EntityKey> SearchId;
typedef stk::search::Point<Scalar> Point;
typedef stk::search::Sphere<Scalar> Sphere;
typedef std::vector< std::pair<Sphere,SearchId> > SphereIdVector;
typedef std::vector<std::pair<SearchId,SearchId> > SearchPairVector;
typedef std::vector<std::pair<stk::mesh::EntityKey,stk::mesh::EntityKey> > SearchPairSet;

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tCoarseSearch, basic)
{
  using Teuchos::RCP;
  using Teuchos::Tuple;

  panzer_stk::SquareQuadMeshFactory mesh_factory;

  // setup mesh
  /////////////////////////////////////////////
  RCP<panzer_stk::STK_Interface> mesh;
  {
     RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
     pl->set("X Blocks",2);
     pl->set("Y Blocks",1);
     pl->set("X Elements",6);
     pl->set("Y Elements",4);
     mesh_factory.setParameterList(pl);
     mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
  }
  TEST_ASSERT(mesh!=Teuchos::null);

  // create a side 1 and side 2 vector (trivially for now)

  SphereIdVector side_1_vector, side_2_vector;
  double tolerance = 1e-10;
  SearchPairVector results;
  stk::mesh::EntityRank r = mesh->getNodeRank();
  auto myrank = mesh->getBulkData()->parallel_rank();

  // this test is designed for 4 mpi processes
  // each process will own a line of points, side 1 (0 and 1) side 2 (2 and 3)
  // process 0 and 3 will only have matches for 1/2 of their points

  size_t npts = 10;
  size_t start[4] = {0,npts,npts/2,npts/2+npts};

  for (size_t ind=0; ind<npts; ++ind) {
    stk::mesh::EntityId id(ind+myrank*1000);
    stk::mesh::EntityKey key(r,id);
    SearchId search_id(key,myrank);
    Point center(0,start[myrank]+ind,0);
    if (myrank<2) side_1_vector.emplace_back( Sphere(center,tolerance), search_id );
    if (myrank>1) side_2_vector.emplace_back( Sphere(center,tolerance), search_id );
  }

  stk::search::coarse_search(side_1_vector,side_2_vector,stk::search::KDTREE,MPI_COMM_WORLD,results);

  bool fail = true;
  unsigned int matches = 0;

  for ( auto & x : results ) {
    if (myrank<2) {
      if (x.first.proc()==myrank) ++matches;
    } else {
      if (x.second.proc()==myrank) ++matches;
    }
  }

  if ( (myrank == 0) || (myrank == 3) ) {
    if (matches == npts/2) fail = false;
  } else {
    if (matches == npts) fail = false;
  }

  TEUCHOS_ASSERT(fail!=true);

}

}
