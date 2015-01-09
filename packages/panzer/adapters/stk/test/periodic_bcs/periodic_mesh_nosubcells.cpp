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

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Tuple.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_PeriodicBC_Matcher.hpp"
#include "Panzer_STK_PeriodicBC_Parser.hpp"
#include "Panzer_STK_PeriodicBC_MatchConditions.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_DOFManager.hpp"

#include "Epetra_MpiComm.h"

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"

#include <string>

typedef Intrepid::FieldContainer<double> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk_classic {


  template <typename IntrepidType>
  RCP<const panzer::FieldPattern> buildFieldPattern()
  {
     // build a geometric pattern from a single basis
     RCP<Intrepid::Basis<double,FieldContainer> > basis = rcp(new IntrepidType);
     RCP<const panzer::FieldPattern> pattern = rcp(new panzer::IntrepidFieldPattern(basis));
     return pattern;
  }

  TEUCHOS_UNIT_TEST(periodic_mesh_nosubcells, serial)
  {
    using Teuchos::RCP;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==1);

    // panzer::pauseToAttach();

    panzer_stk_classic::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",1);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",2);
       pl->set("Z Elements",2);
       pl->set("X Procs",1);
       pl->set("Y Procs",1);
       pl->set("Z Procs",1);
       pl->set("X0",0.0);
       pl->set("Y0",0.0);
       pl->set("Z0",0.0);
       pl->set("Xf",6.0);
       pl->set("Yf",3.0);
       pl->set("Zf",4.0);
       mesh_factory.setParameterList(pl);
       mesh_factory.buildSubcells(false);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk_classic::PlaneMatcher side_matcher(1,2);
    panzer_stk_classic::PlaneMatcher front_matcher(0,1);
    mesh->addPeriodicBC(panzer_stk_classic::buildPeriodicBC_Matcher("right","left",side_matcher));
    mesh->addPeriodicBC(panzer_stk_classic::buildPeriodicBC_Matcher("front","back",front_matcher));

    mesh->writeToExodus("file.exo");

    // connection manager
    /////////////////////////////////////////////
    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer> >();

    Teuchos::RCP<panzer::ConnManager<int,panzer::Ordinal64> > connMngr 
          = Teuchos::rcp(new panzer_stk_classic::STKConnManager<panzer::Ordinal64>(mesh));

    Teuchos::RCP<panzer::DOFManager<int,panzer::Ordinal64> > dofManager
          = Teuchos::rcp(new panzer::DOFManager<int,panzer::Ordinal64>(connMngr,MPI_COMM_WORLD));
    dofManager->addField("VAR",fp);
    dofManager->buildGlobalUnknowns();

    std::vector<panzer::Ordinal64> owned;
    dofManager->getOwnedIndices(owned);

    std::size_t unkCount = owned.size();
    std::size_t nodeCount = mesh->getEntityCounts(mesh->getNodeRank());
    out << "Unknown Count = " << unkCount 
        << ", Node Count = "  << nodeCount << std::endl;

    TEST_ASSERT(unkCount < nodeCount);
  }
}
