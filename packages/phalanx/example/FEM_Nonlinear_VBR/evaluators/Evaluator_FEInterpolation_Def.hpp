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

#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

//**********************************************************************
PHX_EVALUATOR_CTOR(FEInterpolation,p) :
  val_node(p.get<std::string>("Node Variable Name"), 
	   p.get< Teuchos::RCP<PHX::DataLayout> >("Node Data Layout") ),
  val_qp(p.get<std::string>("QP Variable Name"), 
	 p.get< Teuchos::RCP<PHX::DataLayout> >("QP Scalar Data Layout") ),
  val_grad_qp(p.get<std::string>("Gradient QP Variable Name"), 
	      p.get< Teuchos::RCP<PHX::DataLayout> >("QP Vector Data Layout") )
{ 
  this->addDependentField(val_node);
  this->addEvaluatedField(val_qp);
  this->addEvaluatedField(val_grad_qp);
  
  this->setName("FEInterpolation");
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(FEInterpolation,fm)
{
  this->utils.setFieldData(val_node,fm);
  this->utils.setFieldData(val_qp,fm);
  this->utils.setFieldData(val_grad_qp,fm);

  // Get dimensions of MDArray
  typename std::vector< typename PHX::template MDField<ScalarT,PHX::NaturalOrder,Cell,Node>::size_type > dims;
  val_node.dimensions(dims);
  num_nodes = dims[1];

  val_grad_qp.dimensions(dims);
  num_qp = dims[1];
  num_dim = dims[2];
}

//**********************************************************************
PHX_EVALUATE_FIELDS(FEInterpolation,cell_data)
{ 
  using namespace PHX;

  std::vector<Element_Linear2D>::iterator cell_it = cell_data.begin;

  // Loop over number of cells
  for (std::size_t cell = 0; cell < cell_data.num_cells; ++cell) {
    
    const Array<double,NaturalOrder,QuadPoint,Node>& phi = 
      cell_it->basisFunctions();

    const Array<double,NaturalOrder,QuadPoint,Node,Dim>& grad_phi = 
      cell_it->basisFunctionGradientsRealSpace();

    // Loop over quad points of cell
    for (std::size_t qp = 0; qp < num_qp; ++qp) {
      
      val_qp(cell,qp) = 0.0;

      for (std::size_t dim = 0; dim < num_dim; ++dim)
	val_grad_qp(cell,qp,dim) = 0.0;

      // Sum nodal contributions to qp
      for (std::size_t node = 0; node < num_nodes; ++node) {

	val_qp(cell,qp) += phi(qp,node) * val_node(cell,node);
	
	for (std::size_t dim = 0; dim < num_dim; ++dim)
	  val_grad_qp(cell,qp,dim) += 
	    grad_phi(qp,node,dim) * val_node(cell,node);
	
      }
    }
    
    ++cell_it;
 
  }
    
//   std::cout << "FEINterpolation: val_node" << std::endl;
//   val_node.print(std::cout,true);
//   std::cout << "FEINterpolation: val_qp" << std::endl;
//   val_qp.print(std::cout,true);
//   std::cout << "FEINterpolation: val_grad_qp" << std::endl;
//   val_grad_qp.print(std::cout,true);

}

//**********************************************************************
