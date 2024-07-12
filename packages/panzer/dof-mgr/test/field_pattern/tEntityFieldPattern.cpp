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
