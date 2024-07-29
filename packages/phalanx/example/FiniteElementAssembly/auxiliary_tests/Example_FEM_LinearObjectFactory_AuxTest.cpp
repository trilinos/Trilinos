// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_config.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Mesh.hpp"
#include "LinearObjectFactory.hpp"
#include <limits>

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TEUCHOS_UNIT_TEST(mesh, coord_gids)
{
  using namespace std;
  using namespace Teuchos;
  
  const int nx=2;
  const int ny=3;
  const int nz=4;
  const double lx = 1.0;
  const double ly = 2.0;
  const double lz = 3.0;
  phx_example::Mesh mesh(nx,ny,nz,lx,ly,lz);

  const int num_equations = 2;
  phx_example::LinearObjectFactory lof(mesh.getNumNodes(),
                                       num_equations,
                                       mesh.getGlobalIndices());

  const auto x = lof.createSolutionVector("x");
  TEST_EQUALITY(x.size(), ((nx+1)*(ny+1)*(nz+1)*num_equations));

  auto J = lof.createJacobianMatrix("J");

  // Mirror calls deep_copy on underlying obejcts so host copy is up
  // to date.
  auto graph = Kokkos::create_mirror(J.graph);
  
  const int max_entries_per_row = 27*num_equations;
  const int min_entries_per_row = 8*num_equations;
  for (decltype(graph)::size_type row=0; row < graph.numRows(); ++row) {
    TEST_ASSERT(graph.rowConst(row).length <= max_entries_per_row);
    TEST_ASSERT(graph.rowConst(row).length >= min_entries_per_row);
    std::cout << "row(" << row << ") = ";
    for (int j=0; j < graph.rowConst(row).length; ++j) {
      std::cout << graph.rowConst(row).colidx(j);
      if (j < (graph.rowConst(row).length - 1))
        std::cout << ",";
    }
    std::cout << std::endl;   
  }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
