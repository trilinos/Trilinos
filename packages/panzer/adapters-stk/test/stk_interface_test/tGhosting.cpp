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
#include "Panzer_STK_SetupUtilities.hpp"

#include "Shards_BasicTopologies.hpp"

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
//#include <stk_rebalance/ZoltanPartition.hpp>

namespace panzer_stk {

inline bool XOR(bool A,bool B)
{ return ! ( (A && B) || ( !A && !B)); }

class LocalIdCompare {
public:
   LocalIdCompare(const Teuchos::RCP<const STK_Interface> & mesh) : mesh_(mesh) {}
   bool operator()(stk::mesh::Entity a,stk::mesh::Entity b) const
   { return mesh_->elementLocalId(a) < mesh_->elementLocalId(b); }

private:
   Teuchos::RCP<const STK_Interface> mesh_;
};

// This test was modified to its current lame state when the
// construction of the local element IDs was automated in the
// STK_Interface. (Independent of order of addition in the mesh

TEUCHOS_UNIT_TEST(tGhosting, get_neighbor_elements)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",1);
   pl->set("X Procs",2);
   pl->set("Y Procs",2);
   pl->set("X Elements",2);
   pl->set("Y Elements",3);

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   TEUCHOS_ASSERT(numprocs==4);

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   mesh->writeToExodus("TEST.exo");

   {
     std::vector<stk::mesh::Entity> neighbors;
     mesh->getNeighborElements(neighbors);

     std::size_t vec[4];
     vec[0] = 8;
     vec[1] = 8;
     vec[2] = 6;
     vec[3] = 6;

     TEST_EQUALITY(neighbors.size(),vec[rank]);
   }

   {
     std::vector<stk::mesh::Entity> neighbors;
     mesh->getNeighborElements("eblock-0_0",neighbors);

     std::size_t vec[4];
     vec[0] = 4;
     vec[1] = 4;
     vec[2] = 3;
     vec[3] = 3;

     TEST_EQUALITY(neighbors.size(),vec[rank]);
   }

   {
     std::vector<stk::mesh::Entity> neighbors;
     mesh->getNeighborElements("eblock-1_0",neighbors);

     std::size_t vec[4];
     vec[0] = 4;
     vec[1] = 4;
     vec[2] = 3;
     vec[3] = 3;

     TEST_EQUALITY(neighbors.size(),vec[rank]);
   }
}

}
