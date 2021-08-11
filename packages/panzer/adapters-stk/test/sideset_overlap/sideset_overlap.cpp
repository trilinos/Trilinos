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

#include "Panzer_STK_CheckSidesetOverlap.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"

namespace panzer_test {

  TEUCHOS_UNIT_TEST(periodic_bcs, sidesetOverlap_2D)
  {
    // setup mesh
    Teuchos::RCP<panzer_stk::STK_Interface> mesh;
    {
      Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Blocks",2);
      pl->set("Y Blocks",1);
      pl->set("X Elements",2);
      pl->set("Y Elements",2);
      pl->set("X Procs",2);
      pl->set("Y Procs",2);
      panzer_stk::SquareQuadMeshFactory mesh_factory;
      mesh_factory.setParameterList(pl);
      mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // Test is hard coded to four processes to make sure we cover
    // parallel split of sidesets
    TEST_ASSERT(mesh->getComm()->getSize() == 4);

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("left","right",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","bottom",*mesh) );

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("top","bottom",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","bottom",*mesh) );
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, sidesetOverlap_3D)
  {
    // setup mesh
    Teuchos::RCP<panzer_stk::STK_Interface> mesh;
    {
      Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",2);
       pl->set("Z Blocks",2);
       pl->set("X Elements",2);
       pl->set("Y Elements",2);
       pl->set("Z Elements",2);
       pl->set("X Procs",1);
       pl->set("Y Procs",2);
       pl->set("Z Procs",2);
       panzer_stk::CubeHexMeshFactory mesh_factory;
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // Test is hard coded to four processes to make sure we cover
    // parallel split of sidesets
    TEST_ASSERT(mesh->getComm()->getSize() == 4);

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("left","right",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","bottom",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","front",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","back",*mesh) );

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("top","bottom",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","bottom",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","front",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","back",*mesh) );

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("back","front",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("back","left",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("back","right",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("back","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("back","bottom",*mesh) );
  }
}
