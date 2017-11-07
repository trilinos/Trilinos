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
#include "Workset.hpp"
#include "WorksetBuilder.hpp"
#include <limits>

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TEUCHOS_UNIT_TEST(workset, builder)
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

  // 24 elements / size 3 = 8 worksets
  {
    std::vector<Workset> worksets(9); // will be resized internally
    const int workset_size = 3;
    WorksetBuilder builder;
    builder.buildWorksets(workset_size,mesh,worksets);
    TEST_EQUALITY(worksets.size(),static_cast<size_t>(8));
    TEST_EQUALITY(worksets[0].num_cells_,static_cast<int>(3));
    TEST_EQUALITY(worksets[0].first_cell_global_index_,static_cast<int>(0));
    TEST_EQUALITY(worksets[1].num_cells_,static_cast<int>(3));
    TEST_EQUALITY(worksets[1].first_cell_global_index_,static_cast<int>(3));
    TEST_EQUALITY(worksets[2].num_cells_,static_cast<int>(3));
    TEST_EQUALITY(worksets[2].first_cell_global_index_,static_cast<int>(6));
    TEST_EQUALITY(worksets[3].num_cells_,static_cast<int>(3));
    TEST_EQUALITY(worksets[3].first_cell_global_index_,static_cast<int>(9));
    TEST_EQUALITY(worksets[4].num_cells_,static_cast<int>(3));
    TEST_EQUALITY(worksets[4].first_cell_global_index_,static_cast<int>(12));
    TEST_EQUALITY(worksets[5].num_cells_,static_cast<int>(3));
    TEST_EQUALITY(worksets[5].first_cell_global_index_,static_cast<int>(15));
    TEST_EQUALITY(worksets[6].num_cells_,static_cast<int>(3));
    TEST_EQUALITY(worksets[6].first_cell_global_index_,static_cast<int>(18));
    TEST_EQUALITY(worksets[7].num_cells_,static_cast<int>(3));
    TEST_EQUALITY(worksets[7].first_cell_global_index_,static_cast<int>(21));
  }

  // 24 elements / size 5 = 5 worksets
  {
    std::vector<Workset> worksets(1); // will be resized internally
    const int workset_size = 5;
    WorksetBuilder builder;
    builder.buildWorksets(workset_size,mesh,worksets);
    TEST_EQUALITY(worksets.size(),static_cast<size_t>(5));
    TEST_EQUALITY(worksets[0].num_cells_,static_cast<int>(5));
    TEST_EQUALITY(worksets[0].first_cell_global_index_,static_cast<int>(0));
    TEST_EQUALITY(worksets[1].num_cells_,static_cast<int>(5));
    TEST_EQUALITY(worksets[1].first_cell_global_index_,static_cast<int>(5));
    TEST_EQUALITY(worksets[2].num_cells_,static_cast<int>(5));
    TEST_EQUALITY(worksets[2].first_cell_global_index_,static_cast<int>(10));
    TEST_EQUALITY(worksets[3].num_cells_,static_cast<int>(5));
    TEST_EQUALITY(worksets[3].first_cell_global_index_,static_cast<int>(15));
    TEST_EQUALITY(worksets[4].num_cells_,static_cast<int>(4));
    TEST_EQUALITY(worksets[4].first_cell_global_index_,static_cast<int>(20));
  }

}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
