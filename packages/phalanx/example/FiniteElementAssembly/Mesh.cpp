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

#include "Mesh.hpp"

namespace phx_example {
  
//**********************************************************************
Mesh::Mesh(const int num_elements_x,
           const int num_elements_y,
           const int num_elements_z,
           const double length_x,
           const double length_y,
           const double length_z) :
  nex_(num_elements_x),
  ney_(num_elements_y),
  nez_(num_elements_z),
  nel_(nex_*ney_*nez_),
  lx_(length_x),
  ly_(length_y),
  lz_(length_z),
  gids_("gids",nel_,8),  // Assume nodal linear elements
  coords_("coords",nel_,8,3), // Assume nodal linear elements
  qp_coords_("qp_coords",nel_,8,3), // Assume nodal linear elements
  qp_("qp",8,3),
  weights_("weights",8),
  basis_("basis",8,8),
  grad_basis_ref_("grad_basis_ref",8,8,3),
  jac_("jac",nel_,8,3,3),
  inv_jac_("inv_jac",nel_,8,3,3),
  det_jac_("det_jac",nel_,8),
  grad_basis_real_("grad_basis_real",nel_,8,8,3)
{  
  // Create unstructured data objects even though this is structured
  // underneath
  const double dx = lx_ / static_cast<double>(nex_);
  const double dy = ly_ / static_cast<double>(ney_);
  const double dz = lz_ / static_cast<double>(nez_);

  // Loop over elements and assign a global index
  const int nx_offset = (ney_ + 1) * (nez_ + 1);
  const int ny_offset = nez_ + 1;
  const int nz_offset = 1;
  for (int nx=0; nx < nex_; ++nx) {  
    for (int ny=0; ny < ney_; ++ny) {
      for (int nz=0; nz < nez_; ++nz) {

        const int cell = nx * ney_ * nez_ + ny * nez_ + nz;
        // std::cout << "cell=" << cell << std::endl;
        // std::cout << "offsets(" << nx_offset << "," << ny_offset
        //           << "," << nz_offset << ")" << std::endl;
        
        gids_(cell,0) = nx * nx_offset + ny * ny_offset + nz;
        gids_(cell,1) = gids_(cell,0) + nx_offset;
        gids_(cell,2) = gids_(cell,0) + nx_offset + ny_offset;
        gids_(cell,3) = gids_(cell,0) + ny_offset;
        gids_(cell,4) = gids_(cell,0) + nz_offset;
        gids_(cell,5) = gids_(cell,1) + nz_offset;
        gids_(cell,6) = gids_(cell,2) + nz_offset;
        gids_(cell,7) = gids_(cell,3) + nz_offset;

        coords_(cell,0,0) = nx * dx;
        coords_(cell,0,1) = ny * dy;
        coords_(cell,0,2) = nz * dz;
        coords_(cell,1,0) = nx * dx + dx;
        coords_(cell,1,1) = ny * dy;
        coords_(cell,1,2) = nz * dz;
        coords_(cell,2,0) = nx * dx + dx;
        coords_(cell,2,1) = ny * dy + dy;
        coords_(cell,2,2) = nz * dz;
        coords_(cell,3,0) = nx * dx;
        coords_(cell,3,1) = ny * dy + dy;
        coords_(cell,3,2) = nz * dz;
        coords_(cell,4,0) = nx * dx;
        coords_(cell,4,1) = ny * dy;
        coords_(cell,4,2) = nz * dz + dz;
        coords_(cell,5,0) = nx * dx + dx;
        coords_(cell,5,1) = ny * dy;
        coords_(cell,5,2) = nz * dz + dz;
        coords_(cell,6,0) = nx * dx + dx;
        coords_(cell,6,1) = ny * dy + dy;
        coords_(cell,6,2) = nz * dz + dz;
        coords_(cell,7,0) = nx * dx;
        coords_(cell,7,1) = ny * dy + dy;
        coords_(cell,7,2) = nz * dz + dz;
      } // nz
    } // ny
  } // nx

  const double value = 1.0/std::sqrt(3.0);
  qp_(0,0) = -value;
  qp_(0,1) = -value;
  qp_(0,2) = -value;
  qp_(1,0) = +value;
  qp_(1,1) = -value;
  qp_(1,2) = -value;
  qp_(2,0) = +value;
  qp_(2,1) = +value;
  qp_(2,2) = -value;
  qp_(3,0) = -value;
  qp_(3,1) = +value;
  qp_(3,2) = -value;
  qp_(4,0) = -value;
  qp_(4,1) = -value;
  qp_(4,2) = +value;
  qp_(5,0) = +value;
  qp_(5,1) = -value;
  qp_(5,2) = +value;
  qp_(6,0) = +value;
  qp_(6,1) = +value;
  qp_(6,2) = +value;
  qp_(7,0) = -value;
  qp_(7,1) = +value;
  qp_(7,2) = +value;

  for (int qp=0; qp < static_cast<int>(qp_.extent(0)); ++qp) {
    double chi = qp_(qp,0);
    double eta = qp_(qp,1);
    double mu = qp_(qp,2);
    basis_(qp,0) = (1.-chi) * (1.-eta) * (1.-mu) / 8.0;
    basis_(qp,1) = (1.+chi) * (1.-eta) * (1.-mu) / 8.0;
    basis_(qp,2) = (1.+chi) * (1.+eta) * (1.-mu) / 8.0;
    basis_(qp,3) = (1.-chi) * (1.+eta) * (1.-mu) / 8.0;
    basis_(qp,4) = (1.-chi) * (1.-eta) * (1.+mu) / 8.0;
    basis_(qp,5) = (1.+chi) * (1.-eta) * (1.+mu) / 8.0;
    basis_(qp,6) = (1.+chi) * (1.+eta) * (1.+mu) / 8.0;
    basis_(qp,7) = (1.-chi) * (1.+eta) * (1.+mu) / 8.0;

    grad_basis_ref_(qp,0,0) = - (1.-eta) * (1.-mu) / 8.0;
    grad_basis_ref_(qp,0,1) = - (1.-chi) * (1.-mu) / 8.0;
    grad_basis_ref_(qp,0,2) = - (1.-chi) * (1.-eta) / 8.0;

    grad_basis_ref_(qp,1,0) = + (1.-eta) * (1.-mu) / 8.0;
    grad_basis_ref_(qp,1,1) = - (1.+chi) * (1.-mu) / 8.0;
    grad_basis_ref_(qp,1,2) = - (1.+chi) * (1.-eta) / 8.0;
    
    grad_basis_ref_(qp,2,0) = + (1.+eta) * (1.-mu) / 8.0;
    grad_basis_ref_(qp,2,1) = + (1.+chi) * (1.-mu) / 8.0;
    grad_basis_ref_(qp,2,2) = - (1.+chi) * (1.+eta) / 8.0;

    grad_basis_ref_(qp,3,0) = - (1.+eta) * (1.-mu) / 8.0;
    grad_basis_ref_(qp,3,1) = + (1.-chi) * (1.-mu) / 8.0;
    grad_basis_ref_(qp,3,2) = - (1.-chi) * (1.+eta) / 8.0;

    grad_basis_ref_(qp,4,0) = - (1.-eta) * (1.+mu) / 8.0;
    grad_basis_ref_(qp,4,1) = - (1.-chi) * (1.+mu) / 8.0;
    grad_basis_ref_(qp,4,2) = + (1.-chi) * (1.-eta) / 8.0;

    grad_basis_ref_(qp,5,0) = + (1.-eta) * (1.+mu) / 8.0;
    grad_basis_ref_(qp,5,1) = - (1.+chi) * (1.+mu) / 8.0;
    grad_basis_ref_(qp,5,2) = + (1.+chi) * (1.-eta) / 8.0;

    grad_basis_ref_(qp,6,0) = + (1.+eta) * (1.+mu) / 8.0;
    grad_basis_ref_(qp,6,1) = + (1.+chi) * (1.+mu) / 8.0;
    grad_basis_ref_(qp,6,2) = + (1.+chi) * (1.+eta) / 8.0;

    grad_basis_ref_(qp,7,0) = - (1.+eta) * (1.+mu) / 8.0;
    grad_basis_ref_(qp,7,1) = + (1.-chi) * (1.+mu) / 8.0;
    grad_basis_ref_(qp,7,2) = + (1.-chi) * (1.+eta) / 8.0;
  }

  // weights
  Kokkos::deep_copy(weights_, 1.0);

  // jac
  Kokkos::TeamPolicy<PHX::exec_space> tp(nel_,Kokkos::AUTO());
  Kokkos::deep_copy(jac_,0.0);
  Kokkos::parallel_for(Kokkos::TeamPolicy<ComputeJac_Tag,PHX::exec_space>(nel_,Kokkos::AUTO()),*this);

  // inv_jac, det_jac
  Kokkos::parallel_for(Kokkos::TeamPolicy<ComputeInvJac_Tag,PHX::exec_space>(nel_,Kokkos::AUTO()),*this);

  // qp coords
  Kokkos::deep_copy(qp_coords_,0.0);
  Kokkos::parallel_for(Kokkos::TeamPolicy<ComputeCoords_Tag,PHX::exec_space>(nel_,Kokkos::AUTO()),*this);

  // transform basis gradients to real space
  Kokkos::deep_copy(grad_basis_real_,0.0);
  Kokkos::parallel_for(Kokkos::TeamPolicy<ComputeGradBasisReal_Tag,PHX::exec_space>(nel_,Kokkos::AUTO()),*this);
}

//**********************************************************************
KOKKOS_INLINE_FUNCTION
void Mesh::operator() (const ComputeJac_Tag& , const team_t& team) const
{
  const int cell = team.league_rank();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,qp_.extent(0)), [=] (const int& qp) {
      for (int basis=0; basis < static_cast<int>(basis_.extent(1)); ++basis) {
        for (int i=0; i < 3; ++i) {
          for (int j=0; j < 3; ++j) {
            jac_(cell,qp,i,j) += coords_(cell,basis,i) * grad_basis_ref_(qp,basis,j);
          }
        }
      }
  });
}

//**********************************************************************
KOKKOS_INLINE_FUNCTION
void Mesh::operator() (const ComputeInvJac_Tag& , const team_t& team) const
{
  const int cell = team.league_rank();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,qp_.extent(0)), [=] (const int& qp) {
      inv_jac_(cell,qp,0,0) = jac_(cell,qp,1,1) * jac_(cell,qp,2,2) - jac_(cell,qp,1,2) * jac_(cell,qp,2,1);
      inv_jac_(cell,qp,1,1) = jac_(cell,qp,2,2) * jac_(cell,qp,0,0) - jac_(cell,qp,2,0) * jac_(cell,qp,0,2);
      inv_jac_(cell,qp,2,2) = jac_(cell,qp,0,0) * jac_(cell,qp,1,1) - jac_(cell,qp,0,1) * jac_(cell,qp,1,0);
      inv_jac_(cell,qp,0,1) = jac_(cell,qp,1,2) * jac_(cell,qp,2,0) - jac_(cell,qp,1,0) * jac_(cell,qp,2,2);
      inv_jac_(cell,qp,1,2) = jac_(cell,qp,2,0) * jac_(cell,qp,0,1) - jac_(cell,qp,2,1) * jac_(cell,qp,0,0);
      inv_jac_(cell,qp,2,0) = jac_(cell,qp,0,1) * jac_(cell,qp,1,2) - jac_(cell,qp,0,2) * jac_(cell,qp,1,1);
      inv_jac_(cell,qp,1,0) = jac_(cell,qp,2,1) * jac_(cell,qp,0,2) - jac_(cell,qp,0,1) * jac_(cell,qp,2,2);
      inv_jac_(cell,qp,2,1) = jac_(cell,qp,0,2) * jac_(cell,qp,1,0) - jac_(cell,qp,1,2) * jac_(cell,qp,0,0);
      inv_jac_(cell,qp,0,2) = jac_(cell,qp,1,0) * jac_(cell,qp,1,1) - jac_(cell,qp,2,0) * jac_(cell,qp,1,1);
      
      det_jac_(cell,qp) = jac_(cell,qp,0,0) * inv_jac_(cell,qp,0,0)
        + jac_(cell,qp,0,1) * inv_jac_(cell,qp,1,0)
        + jac_(cell,qp,0,2) * inv_jac_(cell,qp,2,0);
      
      for (int i=0; i < 3; ++i)
        for (int j=0; j < 3; ++j)
          inv_jac_(cell,qp,i,j) = inv_jac_(cell,qp,i,j) / det_jac_(cell,qp);
  });
}
  
//**********************************************************************
KOKKOS_INLINE_FUNCTION
void Mesh::operator() (const ComputeCoords_Tag& , const team_t& team) const
{
  const int cell = team.league_rank();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,qp_.extent(0)), [=] (const int& qp) {
      for (int basis=0; basis < static_cast<int>(basis_.extent(1)); ++basis) {
        qp_coords_(cell,qp,0) += basis_(qp,basis) * coords_(cell,basis,0);
        qp_coords_(cell,qp,1) += basis_(qp,basis) * coords_(cell,basis,1);
        qp_coords_(cell,qp,2) += basis_(qp,basis) * coords_(cell,basis,2);
      }
  });
}
 
//**********************************************************************
KOKKOS_INLINE_FUNCTION
void Mesh::operator() (const ComputeGradBasisReal_Tag& , const team_t& team) const
{
  const int cell = team.league_rank();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,grad_basis_real_.extent(1)), [=] (const int& qp) {
      const int num_basis = static_cast<int>(grad_basis_real_.extent(2)); 
      for (int basis=0; basis < num_basis; ++basis)
        for (int dim1=0; dim1 < 3; ++dim1)
          for (int dim2=0; dim2 < 3; ++dim2)
            grad_basis_real_(cell,qp,basis,dim1) += grad_basis_ref_(qp,basis,dim1) * inv_jac_(cell,qp,dim1,dim2);
  });
}

//**********************************************************************
int Mesh::getNumElements() const
{ return nel_; }

//**********************************************************************
int Mesh::getNumNodes() const
{ return (nex_+1)*(ney_+1)*(nez_+1); }

//**********************************************************************
const Kokkos::View<int**,PHX::Device>& Mesh::getGlobalIndices() const
{ return gids_; }

//**********************************************************************
const Kokkos::View<double***,PHX::Device>& Mesh::getCoordinates() const
{ return coords_; }

//**********************************************************************
const Kokkos::View<double***,PHX::Device>& Mesh::getQPCoordinates() const
{ return qp_coords_; }

//**********************************************************************
const Kokkos::View<double*,PHX::Device> Mesh::getWeights() const
{return weights_;}

  //**********************************************************************
const Kokkos::View<double**,PHX::Device> Mesh::getBasis() const
{return basis_;}

//**********************************************************************
const Kokkos::View<double***,PHX::Device> Mesh::getGradBasisRef() const
{return grad_basis_ref_; }
  
//**********************************************************************
const Kokkos::View<double****,PHX::Device> Mesh::getJac() const
{ return jac_; }
  
//**********************************************************************
const Kokkos::View<double****,PHX::Device> Mesh::getInvJac() const
{ return inv_jac_; }

//**********************************************************************
const Kokkos::View<double**,PHX::Device> Mesh::getDetJac() const
{ return det_jac_; }

//**********************************************************************
  const Kokkos::View<double****,PHX::Device> Mesh::getGradBasisReal() const
{ return grad_basis_real_; }

//**********************************************************************
void Mesh::print(std::ostream& os) const
{
  os << "Mesh:" << std::endl;
  os << "m_num_elements_x=" << nex_ << std::endl; 
  os << "m_num_elements_y=" << ney_ << std::endl; 
  os << "m_num_elements_z=" << nez_ << std::endl; 
  os << "m_length_x=" << lx_ << std::endl;
  os << "m_length_y=" << ly_ << std::endl;
  os << "m_length_z=" << lz_ << std::endl;
  //os << "m_num_equations=" << neq_ << std::endl;
  for (int cell=0; cell < nel_; ++cell) {
    for (int node=0; node < 8; ++node) {
      os << "el(" << cell << "): gid(" << node << ")="
         << gids_(cell,node)
         << ", coords("
         << coords_(cell,node,0) << ","
         << coords_(cell,node,1) << ","
         << coords_(cell,node,2) << ")"
         << std::endl; 
    }
  }
}

//**********************************************************************
std::ostream& operator<<(std::ostream& os, const Mesh& e)
{
  e.print(os);
  return os;
}

//**********************************************************************

} // namespace phx_example
