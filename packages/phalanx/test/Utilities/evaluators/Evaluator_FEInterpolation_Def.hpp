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
	 p.get< Teuchos::RCP<PHX::DataLayout> >("QP Data Layout") ),
  val_grad_qp(p.get<std::string>("Gradient QP Variable Name"), 
	      p.get< Teuchos::RCP<PHX::DataLayout> >("QP Data Layout") )
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
}

//**********************************************************************
PHX_EVALUATE_FIELDS(FEInterpolation,cell_data)
{ 
  
  const int nodes_per_cell = val_node.fieldTag().dataLayout().size();
  const int qp_per_cell = val_qp.fieldTag().dataLayout().size();

  // Loop over number of cells
  for (std::size_t cell = 0; cell < cell_data.size(); ++cell) {
    
    std::vector< std::vector<double> >& phi = 
      cell_data[cell].getBasisFunctions();
    std::vector< std::vector< MyVector<double> > >& grad_phi = 
      cell_data[cell].getBasisFunctionGradients();

    int node_offset = cell * nodes_per_cell;
    int qp_offset = cell * qp_per_cell;
    
    // Loop over quad points of cell
    for (int qp = 0; qp < qp_per_cell; ++qp) {
      
      val_qp[qp_offset + qp] = 0.0;
      val_grad_qp[qp_offset + qp] = MyVector<ScalarT>(0.0, 0.0, 0.0);
      
      // Sum nodal contributions to qp
      for (int node = 0; node < nodes_per_cell; ++node) {
	
	val_qp[qp_offset + qp] += phi[qp][node] * val_node[node_offset + node];

	val_grad_qp[qp_offset + qp] += 
	  grad_phi[qp][node] * val_node[node_offset + node];
      }      
    }
    
  }
    
}

//**********************************************************************
