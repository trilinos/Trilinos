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


#include <cmath>
#include "Element_Linear2D.hpp"

//**********************************************************************
Element_Linear2D::Element_Linear2D(std::vector<size_type> global_node_ids,
				   size_type global_element_index,
				   size_type local_element_index,
				   std::vector<double> x_node_coords,
				   std::vector<double> y_node_coords) :
  m_global_element_index(global_element_index),
  m_local_element_index(local_element_index),
  m_global_node_ids(global_node_ids),
  m_owns_node(4, false)
{
 const int numQPs = this->numQuadraturePoints();
 const int numNodes = this->numNodes(); 
 const int numDims = 2;
 m_coords = Kokkos::View<double**,PHX::Device>("mcoords", numNodes,numDims);
 m_phi = Kokkos::View<double**,PHX::Device>("phi", numQPs, numNodes);
 m_grad_phi = Kokkos::View<double***,PHX::Device>("m_grad_phi", numQPs, numNodes, numDims);
 m_grad_phi_xy = Kokkos::View<double***,PHX::Device>("m_grad_phi_xy", numQPs, numNodes, numDims);
 m_det_jacobian = Kokkos::View<double*,PHX::Device>("m_det_jacobian", numQPs);
 m_weights = Kokkos::View<double*,PHX::Device>("m_weights", numQPs);  
  // Set node coordinatates
  m_coords(0,0) = x_node_coords[0];
  m_coords(0,1) = y_node_coords[0];
  m_coords(1,0) = x_node_coords[1];
  m_coords(1,1) = y_node_coords[1];
  m_coords(2,0) = x_node_coords[2];
  m_coords(2,1) = y_node_coords[2];
  m_coords(3,0) = x_node_coords[3];
  m_coords(3,1) = y_node_coords[3];
  
  // Quadrature rule defines number of quadrature points - 4
  std::vector<double> chi(4);
  std::vector<double> eta(4);

  chi[0] = -1.0 / std::sqrt(3);
  chi[1] =  1.0 / std::sqrt(3);
  chi[2] =  1.0 / std::sqrt(3);
  chi[3] = -1.0 / std::sqrt(3);

  eta[0] = -1.0 / std::sqrt(3);
  eta[1] = -1.0 / std::sqrt(3);
  eta[2] =  1.0 / std::sqrt(3);
  eta[3] =  1.0 / std::sqrt(3);

  m_weights[0] = 1.0;
  m_weights[1] = 1.0;
  m_weights[2] = 1.0;
  m_weights[3] = 1.0;

  for (size_type qp=0; qp < this->numQuadraturePoints(); ++qp) {
    // Phi
  
    evaluatePhi(chi[qp], eta[qp], m_phi, qp);
    
    // Grad Phi in local element coordinates
    evaluateGradPhi(chi[qp], eta[qp], m_grad_phi, qp);
    
    // Determinant of Jacobian and basis function gradients in 
    // real space
    evaluateDetJacobianAndGradients(chi[qp], eta[qp], m_det_jacobian(qp),
				    m_grad_phi, m_grad_phi_xy, qp);
  }
  
  

  for (int qp=0; qp < this->numQuadraturePoints(); ++qp)
    for (int node=0; node < this->numNodes(); ++node)
      for (int dim=0; dim < numDims; ++dim)
	m_grad_phi_xy(qp,node,dim) = 
 	  (1.0 / m_det_jacobian(qp)) * m_grad_phi(qp,node,dim);

}

//**********************************************************************
Element_Linear2D& Element_Linear2D::operator=(const Element_Linear2D& right)
{
  m_global_element_index = right.m_global_element_index;

  m_local_element_index = right.m_local_element_index;

  m_global_node_ids = right.m_global_node_ids;

  m_owns_node = right.m_owns_node;

//  m_coords_mem = right.m_coords_mem;
//  m_phi_mem = m_phi_mem;
//  m_grad_phi_mem = right.m_grad_phi_mem;
//  m_grad_phi_xy_mem = right.m_grad_phi_xy_mem;
//  m_det_jacobian_mem = right.m_det_jacobian_mem;
//  m_weights_mem = right.m_weights_mem;

  m_coords = right.m_coords;
  m_phi = right.m_phi;
  m_grad_phi = right.m_grad_phi;
  m_grad_phi_xy = right.m_grad_phi_xy;
  m_det_jacobian = right.m_det_jacobian;
  m_weights = right.m_weights;

  return *this;
}

//**********************************************************************
Element_Linear2D::size_type Element_Linear2D::numQuadraturePoints() const
{
  return 4;
}

//**********************************************************************
Element_Linear2D::size_type Element_Linear2D::numNodes() const
{
  return 4;
}

//**********************************************************************
const std::vector<Element_Linear2D::size_type>& 
Element_Linear2D::globalNodeIds() const
{
  return m_global_node_ids;
}

//**********************************************************************
Element_Linear2D::size_type Element_Linear2D::
globalNodeId(size_type local_node_index) const
{
  return m_global_node_ids[local_node_index];
}

//**********************************************************************
Element_Linear2D::size_type Element_Linear2D::globalElementIndex() const
{
  return m_global_element_index;
}

//**********************************************************************
Element_Linear2D::size_type Element_Linear2D::localElementIndex() const
{
  return m_local_element_index;
}

//**********************************************************************
bool Element_Linear2D::ownsNode(size_type local_node_index) const
{
  return m_owns_node[local_node_index];
}

//**********************************************************************
void 
Element_Linear2D::setOwnsNode(size_type local_node_index, bool owns_node)
{
  m_owns_node[local_node_index] = owns_node;
}

//**********************************************************************
const Kokkos::View<double**,PHX::Device> 
Element_Linear2D::nodeCoordinates() const
{
  return m_coords;
}

//**********************************************************************
const Kokkos::View<double**,PHX::Device> 
Element_Linear2D::basisFunctions() const
{
  return m_phi;
}

//**********************************************************************
const Kokkos::View<double***,PHX::Device> 
Element_Linear2D::basisFunctionGradientsLocalSpace() const
{
  return m_grad_phi;
}

//**********************************************************************
const Kokkos::View<double***,PHX::Device> 
Element_Linear2D::basisFunctionGradientsRealSpace() const
{
  return m_grad_phi_xy;
}

//**********************************************************************
const Kokkos::View<double*,PHX::Device>
Element_Linear2D::detJacobian() const
{
  return m_det_jacobian;
}

//**********************************************************************
const Kokkos::View<double*,PHX::Device> 
Element_Linear2D::quadratureWeights() const
{
  return m_weights;
}

//**********************************************************************
void Element_Linear2D::evaluatePhi(double chi, double eta, 
		    Kokkos::View<double**,PHX::Device> phi, int qp)
{
  phi(qp,0) = 0.25 * (1.0 - chi) * (1 - eta);
  phi(qp,1) = 0.25 * (1.0 + chi) * (1 - eta);
  phi(qp,2) = 0.25 * (1.0 + chi) * (1 + eta);
  phi(qp,3) = 0.25 * (1.0 - chi) * (1 + eta);
}

//**********************************************************************
void Element_Linear2D::evaluateGradPhi(double chi, double eta,
	   Kokkos::View<double***,PHX::Device>  grad_phi, int qp)
{
  grad_phi(qp,0,0) = -0.25 * (1 - eta);
  grad_phi(qp,0,1) = -0.25 * (1 - chi);

  grad_phi(qp,1,0) =  0.25 * (1 - eta);
  grad_phi(qp,1,1) = -0.25 * (1 + chi);

  grad_phi(qp,2,0) =  0.25 * (1 + eta);
  grad_phi(qp,2,1) =  0.25 * (1 + chi);

  grad_phi(qp,3,0) = -0.25 * (1 + eta);
  grad_phi(qp,3,1) =  0.25 * (1 - chi);
}

//**********************************************************************
void Element_Linear2D::
evaluateDetJacobianAndGradients(double chi, double eta, double& det_jac,
	const Kokkos::View<double***,PHX::Device>  grad_phi,
	Kokkos::View<double***,PHX::Device> grad_phi_xy, int qp)
{
  double
  dx_dchi = 0.25 * ( ( m_coords(1,0) - m_coords(0,0) ) * (1.0 - eta) +
		     ( m_coords(2,0) - m_coords(3,0) ) * (1.0 + eta) 
		   );

  double
  dx_deta = 0.25 * ( ( m_coords(1,1) - m_coords(0,1) ) * (1.0 - eta) +
		     ( m_coords(2,1) - m_coords(3,1) ) * (1.0 + eta) 
		   );

  double
  dy_dchi = 0.25 * ( ( m_coords(3,0) - m_coords(0,0) ) * (1.0 - chi) +
		     ( m_coords(2,0) - m_coords(1,0) ) * (1.0 + chi) 
		   );

  double
  dy_deta = 0.25 * ( ( m_coords(3,1) - m_coords(0,1) ) * (1.0 - chi) +
		     ( m_coords(2,1) - m_coords(1,1) ) * (1.0 + chi) 
		   );

  det_jac = dx_dchi * dy_deta - dx_deta * dy_dchi;

  double inv_det_jac = 1.0 / det_jac;

  for (size_type node = 0; node < this->numNodes(); ++node) {

    grad_phi_xy(qp, node,0) = inv_det_jac * 
      (dy_deta * grad_phi(qp, node, 0) - dy_dchi * grad_phi(qp, node, 1));

    grad_phi_xy(qp, node,1) = inv_det_jac * 
      (-dx_deta * grad_phi(qp, node, 0) + dx_dchi * grad_phi(qp, node, 1));
  }
}

//**********************************************************************
void Element_Linear2D::print(std::ostream& os) const
{
  os << "Element: gid = " << m_global_element_index << ", lid = " 
     << m_local_element_index << std::endl;
  os << "  coords: " << std::endl;
  for (size_type i=0; i < this->numNodes(); ++i)
    os << "    node[" << i << "]: gid = " << m_global_node_ids[i] 
       << "  coords =  (" << m_coords(i,0) << "," << m_coords(i,1) 
       << "), owns = " << m_owns_node[i] << std::endl;

  if (false) {
    os << "\n  m_grad_phi_xy(QP,Node,Dim):" << std::endl;
    for (size_type qp=0; qp < this->numQuadraturePoints(); ++qp)
      for (size_type node=0; node < this->numNodes(); ++node)
	for (size_type dim=0; dim < 2; ++dim) {
	  os << "    m_grad_phi_xy(" << qp << "," << node << "," << dim
	     << ") = " << m_grad_phi_xy(qp,node,dim) << std::endl;
	}
  }
  
}

//**********************************************************************
std::ostream& operator<<(std::ostream& os, const Element_Linear2D& e)
{
  e.print(os);
  return os;
}

//**********************************************************************
