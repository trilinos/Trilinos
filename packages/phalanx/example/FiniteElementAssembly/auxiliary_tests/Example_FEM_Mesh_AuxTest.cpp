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
#include "Phalanx_KokkosUtilities.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Mesh.hpp"
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
  const int nodes_per_element = 8;
  const int qp_per_element = 8;
  phx_example::Mesh mesh(nx,ny,nz,lx,ly,lz,2);
  out << mesh << std::endl;
  
  // Check GIDS
  const auto& gids = mesh.getGlobalIndices();
  TEST_EQUALITY(static_cast<int>(gids.extent(0)), nx * ny * nz); // num elements
  TEST_EQUALITY(static_cast<int>(gids.extent(1)), nodes_per_element); // linear nodal element
  TEST_EQUALITY(gids(23,0),33);
  TEST_EQUALITY(gids(23,1),53);
  TEST_EQUALITY(gids(23,2),58);
  TEST_EQUALITY(gids(23,3),38);
  TEST_EQUALITY(gids(23,4),34);
  TEST_EQUALITY(gids(23,5),54);
  TEST_EQUALITY(gids(23,6),59);
  TEST_EQUALITY(gids(23,7),39);

  // Check coordinates
  const auto& coords = mesh.getCoordinates();
  const double tol = std::numeric_limits<double>::epsilon() * 100.0;
  const double dx = lx / static_cast<double>(nx);
  const double dy = ly / static_cast<double>(ny);
  const double dz = lz / static_cast<double>(nz);
  TEST_FLOATING_EQUALITY(coords(23,0,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(coords(23,0,1),ly-dy,tol);
  TEST_FLOATING_EQUALITY(coords(23,0,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(coords(23,1,0),lx,tol);
  TEST_FLOATING_EQUALITY(coords(23,1,1),ly-dy,tol);
  TEST_FLOATING_EQUALITY(coords(23,1,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(coords(23,2,0),lx,tol);
  TEST_FLOATING_EQUALITY(coords(23,2,1),ly,tol);
  TEST_FLOATING_EQUALITY(coords(23,2,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(coords(23,3,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(coords(23,3,1),ly,tol);
  TEST_FLOATING_EQUALITY(coords(23,3,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(coords(23,3,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(coords(23,3,1),ly,tol);
  TEST_FLOATING_EQUALITY(coords(23,3,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(coords(23,4,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(coords(23,4,1),ly-dy,tol);
  TEST_FLOATING_EQUALITY(coords(23,4,2),lz,tol);
  
  TEST_FLOATING_EQUALITY(coords(23,5,0),lx,tol);
  TEST_FLOATING_EQUALITY(coords(23,5,1),ly-dy,tol);
  TEST_FLOATING_EQUALITY(coords(23,5,2),lz,tol);
  
  TEST_FLOATING_EQUALITY(coords(23,6,0),lx,tol);
  TEST_FLOATING_EQUALITY(coords(23,6,1),ly,tol);
  TEST_FLOATING_EQUALITY(coords(23,6,2),lz,tol);
  
  TEST_FLOATING_EQUALITY(coords(23,7,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(coords(23,7,1),ly,tol);
  TEST_FLOATING_EQUALITY(coords(23,7,2),lz,tol);

  // Check jacobians via area integration
  const int num_cells = nx * ny * nz;
  const auto& weights = mesh.getWeights();
  const auto& det_jac = mesh.getDetJac();
  for (int cell=0; cell<num_cells; ++cell) {
    double volume = 0.0;

    for (int qp=0; qp < qp_per_element; ++qp) {
      volume += weights(qp) * det_jac(cell,qp);
    }

    TEST_FLOATING_EQUALITY(volume,dx*dy*dz,tol);
  }

  // Test values and gradients by projecting a linear function to the qp
  Kokkos::View<double**> T("T",num_cells,nodes_per_element); // <cell,node>
  for (int cell=0; cell < num_cells; ++cell) {
    for (int node=0; node < nodes_per_element; ++node) {
      const auto& x = coords(cell,node,0);
      const auto& y = coords(cell,node,1);
      const auto& z = coords(cell,node,2);
      T(cell,node) = 2.0 * (lx-x) / lx + 3.0 * (ly-y) / ly + 4.0 * (lz-z) / lz;
    }
  }

  // Perform projection to qp

  Kokkos::View<double**> T_qp("T_qp",num_cells,qp_per_element);
  Kokkos::deep_copy(T_qp,0.0);

  Kokkos::View<double***> DTDX_qp_local("DTDX_qp_local",num_cells,qp_per_element,3);
  Kokkos::deep_copy(DTDX_qp_local,0.0);

  const auto& N = mesh.getBasis();
  const auto& DNDX = mesh.getGradBasisRef();
  const auto& qp_coords = mesh.getQPCoordinates();
  for (int cell=0; cell < num_cells; ++cell) {
    for (int qp=0; qp < qp_per_element; ++qp) {
      for (int basis=0; basis < nodes_per_element; ++basis) {
        T_qp(cell,qp) += T(cell,basis) * N(qp,basis);
        for (int dim=0; dim < 3; ++dim) {
          DTDX_qp_local(cell,qp,dim) += T(cell,basis) * DNDX(qp,basis,dim); // in local basis coords
        }
      }
      //out << "qp_x=" << qp_coords(cell,qp,0) << std::endl;
      TEST_FLOATING_EQUALITY(T_qp(cell,qp), (2.0*(lx-qp_coords(cell,qp,0))/lx)+(3.0*(ly-qp_coords(cell,qp,1))/ly)+(4.0*(lz-qp_coords(cell,qp,2))/lz), tol);
    }
  }
  
  // Test gradient values from reference basis
  Kokkos::View<double***> DTDX_qp("DTDX_qp",num_cells,qp_per_element,3);
  Kokkos::deep_copy(DTDX_qp,0.0);
  const auto& invJac = mesh.getInvJac();
  for (int cell=0; cell < num_cells; ++cell) {
    for (int qp=0; qp < qp_per_element; ++qp) {
      for (int dim1=0; dim1 < 3; ++dim1) {
        for (int dim2=0; dim2 < 3; ++dim2) {
          DTDX_qp(cell,qp,dim1) += DTDX_qp_local(cell,qp,dim1) * invJac(cell,qp,dim1,dim2);
        }
      }
      TEST_FLOATING_EQUALITY(DTDX_qp(cell,qp,0), -2.0/lx, tol);
      TEST_FLOATING_EQUALITY(DTDX_qp(cell,qp,1), -3.0/ly, tol);
      TEST_FLOATING_EQUALITY(DTDX_qp(cell,qp,2), -4.0/lz, tol);
    }
  }

  // Test gradient values from real basis
  Kokkos::deep_copy(DTDX_qp,0.0);
  const auto& grad_basis_real = mesh.getGradBasisReal();
  for (int cell=0; cell < num_cells; ++cell) {
    for (int qp=0; qp < qp_per_element; ++qp) {
      for (int basis=0; basis < nodes_per_element; ++basis) {
        for (int dim=0; dim < 3; ++dim) {
          DTDX_qp(cell,qp,dim) += T(cell,basis) * grad_basis_real(cell,qp,basis,dim);
        }
      }
      TEST_FLOATING_EQUALITY(DTDX_qp(cell,qp,0), -2.0/lx, tol);
      TEST_FLOATING_EQUALITY(DTDX_qp(cell,qp,1), -3.0/ly, tol);
      TEST_FLOATING_EQUALITY(DTDX_qp(cell,qp,2), -4.0/lz, tol);
    }
  }
  
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
