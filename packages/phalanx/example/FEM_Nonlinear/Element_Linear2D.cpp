// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include <cmath>
#include "Element_Linear2D.hpp"

//**********************************************************************
Element_Linear2D::Element_Linear2D(std::vector<unsigned> global_node_ids,
				   std::size_t global_element_index,
				   std::size_t local_element_index,
				   std::vector<double> x_node_coords,
				   std::vector<double> y_node_coords) :
  m_global_element_index(global_element_index),
  m_local_element_index(local_element_index),
  m_global_node_ids(global_node_ids),
  m_owns_node(4, false),
  m_coords_mem(Teuchos::arcp<double>(4*2)),
  m_phi_mem(Teuchos::arcp<double>(4*4)),
  m_grad_phi_mem(Teuchos::arcp<double>(4*4*2)),
  m_grad_phi_xy_mem(Teuchos::arcp<double>(4*4*2)),
  m_det_jacobian_mem(Teuchos::arcp<double>(4)),
  m_weights_mem(Teuchos::arcp<double>(4)),
  m_coords(&(m_coords_mem[0]),4,2),
  m_phi(&(m_phi_mem[0]),4,4),
  m_grad_phi(&(m_grad_phi_mem[0]),4,4,2),
  m_grad_phi_xy(&(m_grad_phi_xy_mem[0]),4,4,2),
  m_det_jacobian(&(m_det_jacobian_mem[0]),4),
  m_weights(&(m_weights_mem[0]),4)
{ 
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

  for (std::size_t qp=0; qp < this->numQuadraturePoints(); ++qp) {
    // Phi
    PHX::Array<double,PHX::NaturalOrder,Node> phi_qp = m_phi.truncate(qp);
    evaluatePhi(chi[qp], eta[qp], phi_qp);
    
    // Grad Phi in local element coordinates
    PHX::Array<double,PHX::NaturalOrder,Node,Dim> grad_phi_qp = 
      m_grad_phi.truncate(qp);
    evaluateGradPhi(chi[qp], eta[qp], grad_phi_qp);
    
    // Determinant of Jacobian and basis function gradients in 
    // real space
    PHX::Array<double,PHX::NaturalOrder,Node,Dim> grad_phi_xy_qp = 
      m_grad_phi_xy.truncate(qp);
    evaluateDetJacobianAndGradients(chi[qp], eta[qp], m_det_jacobian(qp),
				    grad_phi_qp, grad_phi_xy_qp);
  }
  
  
//   for (std::size_t i=0; i < m_grad_phi_xy.size(); ++i)
//     m_grad_phi_xy[i] = 0.0;

//   for (std::size_t qp=0; qp < this->numQuadraturePoints(); ++qp)
//     for (std::size_t node=0; node < this->numNodes(); ++node)
//       for (std::size_t dim=0; dim < 2; ++dim)
// 	m_grad_phi_xy(qp,node,dim) = 
// 	  (1.0 / m_det_jacobian(qp)) * m_grad_phi(qp,node,dim);

}

//**********************************************************************
Element_Linear2D& Element_Linear2D::operator=(const Element_Linear2D& right)
{
  m_global_element_index = right.m_global_element_index;

  m_local_element_index = right.m_local_element_index;

  m_global_node_ids = right.m_global_node_ids;

  m_owns_node = right.m_owns_node;

  m_coords_mem = right.m_coords_mem;
  m_phi_mem = m_phi_mem;
  m_grad_phi_mem = right.m_grad_phi_mem;
  m_grad_phi_xy_mem = right.m_grad_phi_xy_mem;
  m_det_jacobian_mem = right.m_det_jacobian_mem;
  m_weights_mem = right.m_weights_mem;

  m_coords = right.m_coords;
  m_phi = right.m_phi;
  m_grad_phi = right.m_grad_phi;
  m_grad_phi_xy = right.m_grad_phi_xy;
  m_det_jacobian = right.m_det_jacobian;
  m_weights = right.m_weights;

  return *this;
}

//**********************************************************************
std::size_t Element_Linear2D::numQuadraturePoints() const
{
  return 4;
}

//**********************************************************************
std::size_t Element_Linear2D::numNodes() const
{
  return 4;
}

//**********************************************************************
const std::vector<unsigned>& Element_Linear2D::globalNodeIds() const
{
  return m_global_node_ids;
}

//**********************************************************************
unsigned Element_Linear2D::globalNodeId(std::size_t local_node_index) const
{
  return m_global_node_ids[local_node_index];
}

//**********************************************************************
std::size_t Element_Linear2D::globalElementIndex() const
{
  return m_global_element_index;
}

//**********************************************************************
std::size_t Element_Linear2D::localElementIndex() const
{
  return m_local_element_index;
}

//**********************************************************************
bool Element_Linear2D::ownsNode(std::size_t local_node_index) const
{
  return m_owns_node[local_node_index];
}

//**********************************************************************
void 
Element_Linear2D::setOwnsNode(std::size_t local_node_index, bool owns_node)
{
  m_owns_node[local_node_index] = owns_node;
}

//**********************************************************************
const PHX::Array<double,PHX::NaturalOrder,Node,Dim>& 
Element_Linear2D::nodeCoordinates() const
{
  return m_coords;
}

//**********************************************************************
const PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node>& 
Element_Linear2D::basisFunctions() const
{
  return m_phi;
}

//**********************************************************************
const PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node,Dim>& 
Element_Linear2D::basisFunctionGradientsLocalSpace() const
{
  return m_grad_phi;
}

//**********************************************************************
const PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node,Dim>& 
Element_Linear2D::basisFunctionGradientsRealSpace() const
{
  return m_grad_phi_xy;
}

//**********************************************************************
const PHX::Array<double,PHX::NaturalOrder,QuadPoint>& 
Element_Linear2D::detJacobian() const
{
  return m_det_jacobian;
}

//**********************************************************************
const PHX::Array<double,PHX::NaturalOrder,QuadPoint>& 
Element_Linear2D::quadratureWeights() const
{
  return m_weights;
}

//**********************************************************************
void Element_Linear2D::evaluatePhi(double chi, double eta, 
			   PHX::Array<double,PHX::NaturalOrder,Node>& phi)
{
  phi(0) = 0.25 * (1.0 - chi) * (1 - eta);
  phi(1) = 0.25 * (1.0 + chi) * (1 - eta);
  phi(2) = 0.25 * (1.0 + chi) * (1 + eta);
  phi(3) = 0.25 * (1.0 - chi) * (1 + eta);
}

//**********************************************************************
void Element_Linear2D::evaluateGradPhi(double chi, double eta,
		  PHX::Array<double,PHX::NaturalOrder,Node,Dim>& grad_phi)
{
  grad_phi(0,0) = -0.25 * (1 - eta);
  grad_phi(0,1) = -0.25 * (1 - chi);

  grad_phi(1,0) =  0.25 * (1 - eta);
  grad_phi(1,1) = -0.25 * (1 + chi);

  grad_phi(2,0) =  0.25 * (1 + eta);
  grad_phi(2,1) =  0.25 * (1 + chi);

  grad_phi(3,0) = -0.25 * (1 + eta);
  grad_phi(3,1) =  0.25 * (1 - chi);
}

//**********************************************************************
void Element_Linear2D::
evaluateDetJacobianAndGradients(double chi, double eta, double& det_jac,
	       const PHX::Array<double,PHX::NaturalOrder,Node,Dim>& grad_phi,
	       PHX::Array<double,PHX::NaturalOrder,Node,Dim>& grad_phi_xy)
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

  for (std::size_t node = 0; node < this->numNodes(); ++node) {

    grad_phi_xy(node,0) = inv_det_jac * 
      (dy_deta * grad_phi(node, 0) - dy_dchi * grad_phi(node, 1));

    grad_phi_xy(node,1) = inv_det_jac * 
      (-dx_deta * grad_phi(node, 0) + dx_dchi * grad_phi(node, 1));
  }
}

//**********************************************************************
void Element_Linear2D::print(std::ostream& os) const
{
  os << "Element: gid = " << m_global_element_index << ", lid = " 
     << m_local_element_index << std::endl;
  os << "  coords: " << std::endl;
  for (std::size_t i=0; i < this->numNodes(); ++i)
    os << "    node[" << i << "]: gid = " << m_global_node_ids[i] 
       << "  coords =  (" << m_coords(i,0) << "," << m_coords(i,1) 
       << "), owns = " << m_owns_node[i] << std::endl;

  if (false) {
    os << "\n  m_grad_phi_xy(QP,Node,Dim):" << std::endl;
    for (std::size_t qp=0; qp < this->numQuadraturePoints(); ++qp)
      for (std::size_t node=0; node < this->numNodes(); ++node)
	for (std::size_t dim=0; dim < 2; ++dim) {
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
