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


#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

//**********************************************************************
template<typename EvalT, typename Traits>
FEInterpolation<EvalT, Traits>::
FEInterpolation(const Teuchos::ParameterList& p) :
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
template<typename EvalT, typename Traits> 
void FEInterpolation<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(val_node,fm);
  this->utils.setFieldData(val_qp,fm);
  this->utils.setFieldData(val_grad_qp,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void FEInterpolation<EvalT, Traits>::
evaluateFields(typename Traits::EvalData cell_data)
{ 
  
  const int nodes_per_cell = val_node.fieldTag().dataLayout().dimension(1);
  const int qp_per_cell = val_qp.fieldTag().dataLayout().dimension(1);

  std::vector<MyCell>::iterator cell_it = cell_data.begin;

  // Loop over number of cells
  for (std::size_t cell = 0; cell < cell_data.num_cells; ++cell) {
    
    std::vector< std::vector<double> >& phi = 
      cell_it->getBasisFunctions();
    std::vector< std::vector< MyVector<double> > >& grad_phi = 
      cell_it->getBasisFunctionGradients();

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
   
    ++cell_it;
 
  }
    
}

//**********************************************************************
