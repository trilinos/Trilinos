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
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ParameterList.hpp>

#include "PanzerCore_config.hpp"

#include "Panzer_ConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"
#include "Panzer_Relations.hpp"

#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopology.hpp>

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {
namespace unit_test {


TEUCHOS_UNIT_TEST(tFullRelations, test)
{

  RCP<Teuchos::ParameterList>  mesh_pl(new Teuchos::ParameterList("Mesh"));
  mesh_pl->set("X Blocks", 1);
  mesh_pl->set("Y Blocks", 1);
  mesh_pl->set("Y Blocks", 1);
  mesh_pl->set("X Elements", 5);
  mesh_pl->set("Y Elements", 4);
  mesh_pl->set("Z Elements", 3);
  mesh_pl->set("X0", 0.0);
  mesh_pl->set("Y0", 0.0);
  mesh_pl->set("Z0", 0.0);
  mesh_pl->set("Xf", 1.0);
  mesh_pl->set("Yf", 1.0);
  mesh_pl->set("Zf", 1.0);

  RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::CubeTetMeshFactory);
  mesh_factory->setParameterList(mesh_pl);
  RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  const Teuchos::RCP<panzer::ConnManager<int,int> >
    conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager<int>(mesh));


  const CellTopologyData & myCellData = *shards::getCellTopologyData<shards::Pyramid< 5> >();
  struct shards::CellTopology my_topo(&myCellData);

  Relations rel(conn_manager);


  for (int i=0;i<=my_topo.getDimension(); ++i) {
    int cnt = my_topo.getSubcellCount(i);
    std::cout << cnt<<std::endl;
  }

  int cnt = my_topo.getSubcellCount(my_topo.getDimension()-1);
  for (int i=0; i<cnt; ++i) {
    const char *key =my_topo.getCellTopologyData(my_topo.getDimension()-1, i)->name;
    int num_nodes = my_topo.getNodeCount(my_topo.getDimension()-1, i);
    std::cout << num_nodes<<" key "<<key<<std::endl;
    for (int j=0;j<num_nodes; ++j)
      std::cout << " root node " << my_topo.getNodeMap(my_topo.getDimension()-1, i, j)<<std::endl;

  }



}

} // end unit test
} // end panzer
