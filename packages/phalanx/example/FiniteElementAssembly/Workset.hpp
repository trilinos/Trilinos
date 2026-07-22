// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_MY_WORKSET_HPP
#define PHX_EXAMPLE_MY_WORKSET_HPP

#include "Phalanx_config.hpp"
#include "Teuchos_RCP.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "Mesh.hpp"

struct Workset {
  int num_cells_;
  int first_cell_global_index_;

  // Cell global indices <cell,node>
  Kokkos::View<int**,PHX::Device> gids_;
  
  // Weights for integration rule <qp>
  Kokkos::View<double*,PHX::Device> weights_;

  // Determinant of Jacobian <cell,qp>
  Kokkos::View<double**,PHX::Device> det_jac_;  
  
  // Basis <qp,basis>
  Kokkos::View<double**,PHX::Device> basis_;

  // Gradient of basis in real space <cell,qp,basis,dim>
  Kokkos::View<double****,PHX::Device> grad_basis_real_;

  // Solution vector (Required only for Device DAG support)
  Kokkos::View<double*,PHX::Device> global_solution_;

  // Global residual vector, must be atomic (Required only for Device DAG support)
  Kokkos::View<double*,PHX::Device,Kokkos::MemoryTraits<Kokkos::Atomic>> global_residual_atomic_;

  // Global Jacobian  matrix (Required only for Device DAG support)
  KokkosSparse::CrsMatrix<double,int,PHX::Device> global_jacobian_;

  // thread team size (Required only for Host DAG support)
  int team_size_;
  
  // vector team size  (Required only for Host DAG support)
  int vector_size_;
};

#endif
