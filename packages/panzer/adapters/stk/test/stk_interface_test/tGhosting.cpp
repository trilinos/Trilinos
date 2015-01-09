// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include <stk_rebalance/ZoltanPartition.hpp>

namespace panzer_stk_classic {

inline bool XOR(bool A,bool B)
{ return ! ( (A && B) || ( !A && !B)); }

class LocalIdCompare {
public:
   LocalIdCompare(const Teuchos::RCP<const STK_Interface> & mesh) : mesh_(mesh) {}
   bool operator()(stk_classic::mesh::Entity * a,stk_classic::mesh::Entity * b) const 
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

   int numprocs = stk_classic::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk_classic::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   TEUCHOS_ASSERT(numprocs==4);

   SquareQuadMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   mesh->writeToExodus("TEST.exo");

   {
     std::vector<stk_classic::mesh::Entity*> neighbors;
     mesh->getNeighborElements(neighbors);

     std::size_t vec[4];
     vec[0] = 8;
     vec[1] = 8;
     vec[2] = 6;
     vec[3] = 6;

     TEST_EQUALITY(neighbors.size(),vec[rank]);
   }

   {
     std::vector<stk_classic::mesh::Entity*> neighbors;
     mesh->getNeighborElements("eblock-0_0",neighbors);

     std::size_t vec[4];
     vec[0] = 4;
     vec[1] = 4;
     vec[2] = 3;
     vec[3] = 3;

     TEST_EQUALITY(neighbors.size(),vec[rank]);
   }

   {
     std::vector<stk_classic::mesh::Entity*> neighbors;
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
