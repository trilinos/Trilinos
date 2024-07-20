// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_MESH_BUILDER_HPP
#define PHX_EXAMPLE_MESH_BUILDER_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

namespace phx_example {

class Mesh {

  using team_t =  Kokkos::TeamPolicy<PHX::exec_space>::member_type;
  
  const int nex_;
  const int ney_;
  const int nez_;
  const int nel_;
  const double lx_;
  const double ly_;
  const double lz_;
  
  // Cell global indices <cell,node>
  Kokkos::View<int**,PHX::Device> gids_;

  // Coordinates <cell,node,dimension>
  Kokkos::View<double***,PHX::Device> coords_;

  // Coordinates <cell,qp,dimension>
  Kokkos::View<double***,PHX::Device> qp_coords_;

  // Quad points <qp,dimension> in local basis cocordinates
  Kokkos::View<double**,PHX::Device> qp_;

  // Weights for integration rule <qp>
  Kokkos::View<double*,PHX::Device> weights_;
  
  // Basis <qp,basis>
  Kokkos::View<double**,PHX::Device> basis_;

  // Gradient of Basis on reference element <qp,basis,dimension>
  Kokkos::View<double***,PHX::Device> grad_basis_ref_;
  
  // Jacobian transform <cell,qp,dim,dim>
  Kokkos::View<double****,PHX::Device> jac_;
  
  // Inverse Jacobian Transform <cell,qp,dim,dim>
  Kokkos::View<double****,PHX::Device> inv_jac_;

  // Determinant of Jacobian <cell,qp>
  Kokkos::View<double**,PHX::Device> det_jac_;  

  // Gradient of basis in real space <cell,qp,basis,dim>
  Kokkos::View<double****,PHX::Device> grad_basis_real_;

public:
  
  struct ComputeJac_Tag{};
  struct ComputeInvJac_Tag{};
  struct ComputeCoords_Tag{};
  struct ComputeGradBasisReal_Tag{};

  KOKKOS_INLINE_FUNCTION
  void operator() (const ComputeJac_Tag& tag, const team_t& team) const;
  
  KOKKOS_INLINE_FUNCTION
  void operator() (const ComputeInvJac_Tag& tag, const team_t& team) const;
  
  KOKKOS_INLINE_FUNCTION
  void operator() (const ComputeCoords_Tag& tag, const team_t& team) const;
  
  KOKKOS_INLINE_FUNCTION
  void operator() (const ComputeGradBasisReal_Tag& tag, const team_t& team) const;

public:

  Mesh(const int num_elements_x,
       const int num_elements_y,
       const int num_elements_z,
       const double length_x,
       const double length_y,
       const double length_z);

  int getNumElements() const;
  int getNumNodes() const; 
  const Kokkos::View<int**,PHX::Device>& getGlobalIndices() const;
  const Kokkos::View<double***,PHX::Device>& getCoordinates() const;
  const Kokkos::View<double***,PHX::Device>& getQPCoordinates() const;
  const Kokkos::View<double*,PHX::Device> getWeights() const;
  const Kokkos::View<double**,PHX::Device> getBasis() const;
  const Kokkos::View<double***,PHX::Device> getGradBasisRef() const;
  const Kokkos::View<double****,PHX::Device> getJac() const;
  const Kokkos::View<double****,PHX::Device> getInvJac() const;
  const Kokkos::View<double**,PHX::Device> getDetJac() const;
  const Kokkos::View<double****,PHX::Device> getGradBasisReal() const;
  
  void print(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const Mesh& b);

} // namespace phx_example
  
#endif
