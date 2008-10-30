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

#ifndef PHX_EXAMPLE_ELEMENT_2D_HPP
#define PHX_EXAMPLE_ELEMENT_2D_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Dimension.hpp"
#include "Phalanx_Array.hpp"


/** \brief 2D Linear Lagrangian Finite Element

  Note: We mix the integration rule with the element.  Not the best
  idea, but it simplifies the example.

  Local Element Coordinates (chi-eta space):

\verbatim
        eta
        ^
        |
        |
  -1,+1 |  +1,+1
    +--------+
    |   |    |
    |   |    |
    |   +----|------> chi
    |        |
    |        |
    +--------+
  -1,-1    +1,-1

\endverbatim

  Local Node/QP Numbering:

  + = node
  X = quadrature point

\verbatim
    3        2
    +--------+
    |  3  2  |
    |  X  X  |
    |        |
    |  X  X  |
    |  0  1  |
    +--------+
    0        1
\endverbatim
*/
class Element_Linear2D {
  
public:

  Element_Linear2D() {}
  
  Element_Linear2D(std::vector<unsigned> global_node_ids,
		   std::size_t global_element_index,
		   std::size_t local_element_index,
		   std::vector<double> x_node_coords,
		   std::vector<double> y_node_coords);
  
  ~Element_Linear2D() {}
  
  Element_Linear2D& operator=(const Element_Linear2D& right);

  std::size_t numQuadraturePoints() const;

  std::size_t numNodes() const;

  //! Returns a vector of global ids for all nodes
  const std::vector<unsigned>& globalNodeIds() const;

  //! Returns the global node id, given the local node index
  unsigned globalNodeId(std::size_t local_node_index) const;

  //! Returns the global index for this element
  std::size_t globalElementIndex() const;

  //! Returns the local processor index for this element
  std::size_t localElementIndex() const;

  //! Returns true if the node is owned by the calling process
  bool ownsNode(std::size_t local_node_index) const;

  //! Set to true if this node is owned by calling process
  void setOwnsNode(std::size_t local_node_index, bool owns_node=true);

  //! Returns nodal coordinates
  const PHX::Array<double,PHX::NaturalOrder,Node,Dim>& 
  nodeCoordinates() const;
  
  //! Returns values of basis functions at the quadrature points
  const PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node>& 
  basisFunctions() const;
  
  //! Returns gradient of basis functions at the quadrature points
  //! in local element coordinate system, d_phi/d_chi and d_phi/d_eta
  const PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node,Dim>& 
  basisFunctionGradientsLocalSpace() const;
  
  //! Returns gradient of basis functions at the quadrature points
  //! in real coordinate system, d_phi/d_x and d_phi/d_y
  const PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node,Dim>& 
  basisFunctionGradientsRealSpace() const;
  
  //! Returns determinant of jacobian transform values at quadrature points
  const PHX::Array<double,PHX::NaturalOrder,QuadPoint>& 
  detJacobian() const;

  //! Returns jacobian transform values at the quadrature points
  const PHX::Array<double,PHX::NaturalOrder,QuadPoint>& 
  quadratureWeights() const;  
  
  void print(std::ostream& os) const;

private:

  void evaluatePhi(double chi, double eta, 
		   PHX::Array<double,PHX::NaturalOrder,Node>& phi);

  void evaluateGradPhi(double chi, double eta, 
		    PHX::Array<double,PHX::NaturalOrder,Node,Dim>& grad_phi);

  void evaluateDetJacobianAndGradients(double chi, double eta, double& det_jac,
	         const PHX::Array<double,PHX::NaturalOrder,Node,Dim>& grad_phi,
	         PHX::Array<double,PHX::NaturalOrder,Node,Dim>& grad_phi_xy);
  
private:
  
  std::size_t m_global_element_index;

  std::size_t m_local_element_index;

  std::vector<unsigned> m_global_node_ids;

  std::vector<bool> m_owns_node;

  Teuchos::ArrayRCP<double> m_coords_mem;
  
  Teuchos::ArrayRCP<double> m_phi_mem;
  
  Teuchos::ArrayRCP<double> m_grad_phi_mem;

  Teuchos::ArrayRCP<double> m_grad_phi_xy_mem;

  Teuchos::ArrayRCP<double> m_det_jacobian_mem;

  Teuchos::ArrayRCP<double> m_weights_mem;

  PHX::Array<double,PHX::NaturalOrder,Node,Dim> m_coords;
  
  PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node> m_phi;

  PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node,Dim> m_grad_phi;

  PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node,Dim> m_grad_phi_xy;

  PHX::Array<double,PHX::NaturalOrder,QuadPoint> m_det_jacobian;

  PHX::Array<double,PHX::NaturalOrder,QuadPoint> m_weights;

};

std::ostream& operator<<(std::ostream& os, const Element_Linear2D& e);

#endif
