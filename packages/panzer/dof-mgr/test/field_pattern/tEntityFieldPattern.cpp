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

#include <string>
#include <iostream>

#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_EdgeFieldPattern.hpp"
#include "Panzer_FaceFieldPattern.hpp"
#include "Panzer_ElemFieldPattern.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

  std::string note = "***   NOTE: UNIT TEST BASED ON SEPT 2010   ***\n"
    "***   INTREPID AND SHARDS Trilinos-dev     ***\n"
    "***   DOXYGEN WEBSITE                      ***\n";

  // Nodal FP
  TEUCHOS_UNIT_TEST(tEntityFieldPattern, node)
  {
    using namespace shards;
    out << note << std::endl;

    auto hexa = CellTopology(getCellTopologyData<Hexahedron<8>>());
    auto tria = CellTopology(getCellTopologyData<Triangle<6>>());
    auto line = CellTopology(getCellTopologyData<Line<2>>());

    NodalFieldPattern fp_hexa (hexa);
    NodalFieldPattern fp_tria (tria);
    NodalFieldPattern fp_line (line);

    TEST_ASSERT (fp_hexa.getDimension()==3);
    TEST_ASSERT (fp_tria.getDimension()==2);
    TEST_ASSERT (fp_line.getDimension()==1);

    // Nodes
    TEST_ASSERT (fp_hexa.getSubcellIndices(0,0).size()==1);
    TEST_ASSERT (fp_tria.getSubcellIndices(0,0).size()==1);
    TEST_ASSERT (fp_line.getSubcellIndices(0,0).size()==1);

    // Edges
    TEST_ASSERT (fp_hexa.getSubcellIndices(1,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(1,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(1,0).size()==0);

    // Faces
    TEST_ASSERT (fp_hexa.getSubcellIndices(2,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(2,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(2,0).size()==0);

    // Volumes
    TEST_ASSERT (fp_hexa.getSubcellIndices(3,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(3,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(3,0).size()==0);
  }

  // Edge FP
  TEUCHOS_UNIT_TEST(tEntityFieldPattern, edge)
  {
    using namespace shards;
    out << note << std::endl;

    auto hexa = CellTopology(getCellTopologyData<Hexahedron<8>>());
    auto tria = CellTopology(getCellTopologyData<Triangle<6>>());
    auto line = CellTopology(getCellTopologyData<Line<2>>());

    EdgeFieldPattern fp_hexa (hexa);
    EdgeFieldPattern fp_tria (tria);
    EdgeFieldPattern fp_line (line);

    TEST_ASSERT (fp_hexa.getDimension()==3);
    TEST_ASSERT (fp_tria.getDimension()==2);
    TEST_ASSERT (fp_line.getDimension()==1);

    // Nodes
    TEST_ASSERT (fp_hexa.getSubcellIndices(0,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(0,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(0,0).size()==0);

    // Edges
    TEST_ASSERT (fp_hexa.getSubcellIndices(1,0).size()==1);
    TEST_ASSERT (fp_tria.getSubcellIndices(1,0).size()==1);
    TEST_ASSERT (fp_line.getSubcellIndices(1,0).size()==1);

    // Faces
    TEST_ASSERT (fp_hexa.getSubcellIndices(2,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(2,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(2,0).size()==0);

    // Volumes
    TEST_ASSERT (fp_hexa.getSubcellIndices(3,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(3,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(3,0).size()==0);
  }

  // Face FP
  TEUCHOS_UNIT_TEST(tEntityFieldPattern, face)
  {
    using namespace shards;
    out << note << std::endl;

    auto hexa = CellTopology(getCellTopologyData<Hexahedron<8>>());
    auto tria = CellTopology(getCellTopologyData<Triangle<6>>());
    auto line = CellTopology(getCellTopologyData<Line<2>>());

    FaceFieldPattern fp_hexa (hexa);
    FaceFieldPattern fp_tria (tria);
    FaceFieldPattern fp_line (line);

    TEST_ASSERT (fp_hexa.getDimension()==3);
    TEST_ASSERT (fp_tria.getDimension()==2);
    TEST_ASSERT (fp_line.getDimension()==1);

    // Nodes
    TEST_ASSERT (fp_hexa.getSubcellIndices(0,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(0,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(0,0).size()==0);

    // Edges
    TEST_ASSERT (fp_hexa.getSubcellIndices(1,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(1,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(1,0).size()==0);

    // Faces
    TEST_ASSERT (fp_hexa.getSubcellIndices(2,0).size()==1);
    TEST_ASSERT (fp_tria.getSubcellIndices(2,0).size()==1);
    TEST_ASSERT (fp_line.getSubcellIndices(2,0).size()==1);

    // Volumes
    TEST_ASSERT (fp_hexa.getSubcellIndices(3,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(3,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(3,0).size()==0);
  }

  // Elem FP
  TEUCHOS_UNIT_TEST(tEntityFieldPattern, elem)
  {
    using namespace shards;
    out << note << std::endl;

    auto hexa = CellTopology(getCellTopologyData<Hexahedron<8>>());
    auto tria = CellTopology(getCellTopologyData<Triangle<6>>());
    auto line = CellTopology(getCellTopologyData<Line<2>>());

    ElemFieldPattern fp_hexa (hexa);
    ElemFieldPattern fp_tria (tria);
    ElemFieldPattern fp_line (line);

    TEST_ASSERT (fp_hexa.getDimension()==3);
    TEST_ASSERT (fp_tria.getDimension()==2);
    TEST_ASSERT (fp_line.getDimension()==1);

    // Nodes
    TEST_ASSERT (fp_hexa.getSubcellIndices(0,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(0,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(0,0).size()==0);

    // Edges
    TEST_ASSERT (fp_hexa.getSubcellIndices(1,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(1,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(1,0).size()==1); // Edge==Element

    // Faces
    TEST_ASSERT (fp_hexa.getSubcellIndices(2,0).size()==0);
    TEST_ASSERT (fp_tria.getSubcellIndices(2,0).size()==1); // Face==Element
    TEST_ASSERT (fp_line.getSubcellIndices(2,0).size()==0);

    // Volumes
    TEST_ASSERT (fp_hexa.getSubcellIndices(3,0).size()==1);
    TEST_ASSERT (fp_tria.getSubcellIndices(3,0).size()==0);
    TEST_ASSERT (fp_line.getSubcellIndices(3,0).size()==0);
  }
}
