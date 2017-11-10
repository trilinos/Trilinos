// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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
  for (size_t row=0; row < graph.numRows(); ++row) {
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
