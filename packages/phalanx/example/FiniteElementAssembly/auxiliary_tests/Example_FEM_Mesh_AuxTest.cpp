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
#include "Kokkos_NumericTraits.hpp"
#include <limits>
#include <cmath>
#include <algorithm>

#define PHX_TEST_FLOAT_EQUAL(A,B,tol,failed) if ( (std::fabs(A-B) / ((std::fabs(A) > std::fabs(B) ? std::fabs(A) : std::fabs(B)) + 1e-12) ) > tol) failed += 1;

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
  phx_example::Mesh mesh(nx,ny,nz,lx,ly,lz);
  out << mesh << std::endl;
  
  // Check GIDS
  const auto& gids = mesh.getGlobalIndices();
  auto host_gids = Kokkos::create_mirror_view(gids);
  Kokkos::deep_copy(host_gids,gids);
  TEST_EQUALITY(static_cast<int>(host_gids.extent(0)), nx * ny * nz); // num elements
  TEST_EQUALITY(static_cast<int>(host_gids.extent(1)), nodes_per_element); // linear nodal element
  TEST_EQUALITY(host_gids(23,0),33);
  TEST_EQUALITY(host_gids(23,1),53);
  TEST_EQUALITY(host_gids(23,2),58);
  TEST_EQUALITY(host_gids(23,3),38);
  TEST_EQUALITY(host_gids(23,4),34);
  TEST_EQUALITY(host_gids(23,5),54);
  TEST_EQUALITY(host_gids(23,6),59);
  TEST_EQUALITY(host_gids(23,7),39);

  // Check coordinates
  const auto& coords = mesh.getCoordinates();
  auto host_coords = Kokkos::create_mirror_view(coords);
  Kokkos::deep_copy(host_coords,coords);
  const double tol = std::numeric_limits<double>::epsilon() * 100.0;
  const double dx = lx / static_cast<double>(nx);
  const double dy = ly / static_cast<double>(ny);
  const double dz = lz / static_cast<double>(nz);
  TEST_FLOATING_EQUALITY(host_coords(23,0,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,0,1),ly-dy,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,0,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(host_coords(23,1,0),lx,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,1,1),ly-dy,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,1,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(host_coords(23,2,0),lx,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,2,1),ly,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,2,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(host_coords(23,3,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,3,1),ly,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,3,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(host_coords(23,3,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,3,1),ly,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,3,2),lz-dz,tol);
  
  TEST_FLOATING_EQUALITY(host_coords(23,4,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,4,1),ly-dy,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,4,2),lz,tol);
  
  TEST_FLOATING_EQUALITY(host_coords(23,5,0),lx,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,5,1),ly-dy,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,5,2),lz,tol);
  
  TEST_FLOATING_EQUALITY(host_coords(23,6,0),lx,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,6,1),ly,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,6,2),lz,tol);
  
  TEST_FLOATING_EQUALITY(host_coords(23,7,0),lx-dx,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,7,1),ly,tol);
  TEST_FLOATING_EQUALITY(host_coords(23,7,2),lz,tol);

  // Check jacobians via area integration
  const int num_cells = nx * ny * nz;
  const auto& weights = mesh.getWeights();
  const auto& det_jac = mesh.getDetJac();
  auto host_weights = Kokkos::create_mirror_view(weights);
  auto host_det_jac = Kokkos::create_mirror_view(det_jac);
  Kokkos::deep_copy(host_weights,weights);
  Kokkos::deep_copy(host_det_jac,det_jac);
  for (int cell=0; cell<num_cells; ++cell) {
    double volume = 0.0;

    for (int qp=0; qp < qp_per_element; ++qp) {
      volume += host_weights(qp) * host_det_jac(cell,qp);
    }

    TEST_FLOATING_EQUALITY(volume,dx*dy*dz,tol);
  }

  // Test values and gradients by projecting a linear function to the qp
  Kokkos::View<double**,PHX::MemSpace> T("T",num_cells,nodes_per_element); // <cell,node>
  Kokkos::parallel_for("Fill T values",Kokkos::RangePolicy<PHX::ExecSpace>(0,num_cells),KOKKOS_LAMBDA (const int& cell){
      //for (int cell=0; cell < num_cells; ++cell) {
    for (int node=0; node < nodes_per_element; ++node) {
      const auto& x = coords(cell,node,0);
      const auto& y = coords(cell,node,1);
      const auto& z = coords(cell,node,2);
      T(cell,node) = 2.0 * (lx-x) / lx + 3.0 * (ly-y) / ly + 4.0 * (lz-z) / lz;
    }
  });
  PHX::ExecSpace().fence();

  // Perform projection to qp
  Kokkos::View<double**,PHX::MemSpace> T_qp("T_qp",num_cells,qp_per_element);
  Kokkos::deep_copy(T_qp,0.0);

  Kokkos::View<double***,PHX::MemSpace> DTDX_qp_local("DTDX_qp_local",num_cells,qp_per_element,3);
  Kokkos::deep_copy(DTDX_qp_local,0.0);

  const auto& N = mesh.getBasis();
  const auto& DNDX = mesh.getGradBasisRef();
  const auto& qp_coords = mesh.getQPCoordinates();
  int num_failures = 0;
  Kokkos::parallel_reduce("Compute and check T at qp",Kokkos::RangePolicy<PHX::ExecSpace>(0,num_cells),KOKKOS_LAMBDA (const int& cell,int& lfailures) {
    for (int qp=0; qp < qp_per_element; ++qp) {
      for (int basis=0; basis < nodes_per_element; ++basis) {
        T_qp(cell,qp) += T(cell,basis) * N(qp,basis);
        for (int dim=0; dim < 3; ++dim) {
          DTDX_qp_local(cell,qp,dim) += T(cell,basis) * DNDX(qp,basis,dim); // in local basis coords
        }
      }
      double analytic_solution = (2.0*(lx-qp_coords(cell,qp,0))/lx)+(3.0*(ly-qp_coords(cell,qp,1))/ly)+(4.0*(lz-qp_coords(cell,qp,2))/lz);
      PHX_TEST_FLOAT_EQUAL(T_qp(cell,qp), analytic_solution, tol, lfailures);
    }
  },num_failures);
  TEST_EQUALITY(num_failures,0);

  // Test gradient values from reference basis
  num_failures = 0;
  Kokkos::View<double***> DTDX_qp("DTDX_qp",num_cells,qp_per_element,3);
  Kokkos::deep_copy(DTDX_qp,0.0);
  const auto& invJac = mesh.getInvJac();
  Kokkos::parallel_reduce("Compute and check DTDX at qp (ref)",Kokkos::RangePolicy<PHX::ExecSpace>(0,num_cells),KOKKOS_LAMBDA (const int& cell,int& lfailures) {
    for (int qp=0; qp < qp_per_element; ++qp) {
      for (int dim1=0; dim1 < 3; ++dim1) {
        for (int dim2=0; dim2 < 3; ++dim2) {
          DTDX_qp(cell,qp,dim1) += DTDX_qp_local(cell,qp,dim1) * invJac(cell,qp,dim1,dim2);
        }
      }
      PHX_TEST_FLOAT_EQUAL(DTDX_qp(cell,qp,0), -2.0/lx, tol, lfailures);
      PHX_TEST_FLOAT_EQUAL(DTDX_qp(cell,qp,1), -3.0/ly, tol, lfailures);
      PHX_TEST_FLOAT_EQUAL(DTDX_qp(cell,qp,2), -4.0/lz, tol, lfailures);
    }
  },num_failures);
  TEST_EQUALITY(num_failures,0);

  // Test gradient values from real basis
  num_failures = 0;
  Kokkos::deep_copy(DTDX_qp,0.0);
  const auto& grad_basis_real = mesh.getGradBasisReal();
  Kokkos::parallel_reduce("Compute and check DTDX at qp (real)",Kokkos::RangePolicy<PHX::ExecSpace>(0,num_cells),KOKKOS_LAMBDA (const int& cell,int& lfailures) {
    for (int qp=0; qp < qp_per_element; ++qp) {
      for (int basis=0; basis < nodes_per_element; ++basis) {
        for (int dim=0; dim < 3; ++dim) {
          DTDX_qp(cell,qp,dim) += T(cell,basis) * grad_basis_real(cell,qp,basis,dim);
        }
      }
      PHX_TEST_FLOAT_EQUAL(DTDX_qp(cell,qp,0), -2.0/lx, tol, lfailures);
      PHX_TEST_FLOAT_EQUAL(DTDX_qp(cell,qp,1), -3.0/ly, tol, lfailures);
      PHX_TEST_FLOAT_EQUAL(DTDX_qp(cell,qp,2), -4.0/lz, tol, lfailures);
    }
  },num_failures);
  TEST_EQUALITY(num_failures,0);

}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
