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
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_TypeStrings.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename Traits>
ScatterResidual<PHX::MyTraits::Residual, Traits>::
ScatterResidual(const Teuchos::ParameterList& p)
{ 
  scatter_operation = Teuchos::rcp(new PHX::Tag<ScalarT>("Scatter",p.get< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout")));

  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Residual Names"));

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  f = p.get< Teuchos::RCP<Epetra_Vector> >("Residual Vector");
  
  val.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<ScalarT,Cell,Node> mdf(names[eq],dl);
    val[eq] = mdf;
    this->addDependentField(val[eq]);
  }

  this->addEvaluatedField(*scatter_operation);

  this->setName("Scatter Residual");
}

// **********************************************************************
template<typename Traits> 
void ScatterResidual<PHX::MyTraits::Residual, Traits>::
postRegistrationSetup(typename Traits::SetupData data,
		      PHX::FieldManager<Traits>& fm)
{
  for (std::size_t eq = 0; eq < val.size(); ++eq)
    this->utils.setFieldData(val[eq],fm);

  num_nodes = val[0].dimension(1);
  num_eq = val.size();
}

// **********************************************************************
template<typename Traits>
void ScatterResidual<PHX::MyTraits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  std::vector<Element_Linear2D>::iterator element = workset.begin;

  std::size_t cell = 0;
  for (; element != workset.end; ++element,++cell) {
    
    for (int node = 0; node < num_nodes; node++) {
      unsigned node_GID = element->globalNodeId(node);
      //int firstDOF = f->Map().LID(node_GID * num_eq);
      int firstDOF = f->Map().LID(node_GID) * num_eq;
      for (std::size_t eq = 0; eq < val.size(); eq++)
	(*f)[firstDOF + eq] += (val[eq])(cell,node);
    }

  }

}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename Traits>
ScatterResidual<PHX::MyTraits::Jacobian, Traits>::
ScatterResidual(const Teuchos::ParameterList& p)
{ 
  scatter_operation = Teuchos::rcp(new PHX::Tag<ScalarT>("Scatter",p.get< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout")));

  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Residual Names"));

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  f = p.get< Teuchos::RCP<Epetra_Vector> >("Residual Vector");
  
  if (p.template isType< Teuchos::RCP<Epetra_VbrMatrix> >("Jacobian Matrix"))
    Jac = p.template get< Teuchos::RCP<Epetra_VbrMatrix> >("Jacobian Matrix");
      
  val.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<ScalarT,Cell,Node> mdf(names[eq],dl);
    val[eq] = mdf;
    this->addDependentField(val[eq]);
  }

  this->addEvaluatedField(*scatter_operation);

  this->setName("Scatter Residual(Jacobian)");
}

// **********************************************************************
template<typename Traits> 
void ScatterResidual<PHX::MyTraits::Jacobian, Traits>::
postRegistrationSetup(typename Traits::SetupData data,
		      PHX::FieldManager<Traits>& fm)
{
  for (std::size_t eq = 0; eq < val.size(); ++eq)
    this->utils.setFieldData(val[eq],fm);

  num_nodes = val[0].dimension(1);
  num_eq = val.size();
}

// **********************************************************************
template<typename Traits>
void ScatterResidual<PHX::MyTraits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  std::vector<Element_Linear2D>::iterator element = workset.begin;

  std::size_t cell = 0;
  for (; element != workset.end; ++element,++cell) {
    
    // Sum element residual and Jacobian into global residual, Jacobian
    // Loop over nodes in element
    int row, col;
    int lrow, lcol;
    for (int node_row = 0; node_row < num_nodes; node_row++) {
      
      int row_dim;
      int num_block_entries;
      int* block_indices;
      Jac->BeginExtractMyBlockRowView(Jac->LRID(element->globalNodeId(node_row)),row_dim, num_block_entries,block_indices);
      std::vector<Epetra_SerialDenseMatrix*> matrices(num_block_entries); 
      for (std::size_t i = 0; i < matrices.size();  ++i)
	Jac->ExtractEntryView(matrices[i]);

      // Loop over equations per node
      for (int eq_row = 0; eq_row < num_eq; eq_row++) {
	
	lrow = num_eq * node_row + eq_row;
	
	// Global row
	row = static_cast<int>( f->Map().LID(element->globalNodeId(node_row)) * num_eq + eq_row);
	
	// Sum residual
	if (f != Teuchos::null)
	  //f->SumIntoGlobalValue(row, 0, val[eq_row](cell,node_row).val());
	  (*f)[row] += val[eq_row](cell,node_row).val();
	
// 	std::cout << "val[" << eq_row << "](" << cell << "," << node_row << ") = " 
// 		  << val[eq_row](cell,node_row).val() << std::endl;


	// Check derivative array is nonzero
	if (val[eq_row](cell,node_row).hasFastAccess()) {
	  
	  // Loop over nodes in element
	  for (int node_col=0; node_col<num_nodes; node_col++){
	    
	    // Loop over equations per node
	    for (int eq_col=0; eq_col<num_eq; eq_col++) {
	      lcol = num_eq * node_col + eq_col;
	      
	      // Global column
	      col = static_cast<int>(Jac->LCID(element->globalNodeId(node_col)) * num_eq + eq_col);
	      
	      // Sum Jacobian
	      
	      Epetra_SerialDenseMatrix* block = 0;
	      for (int i = 0; i < num_block_entries; ++i) {
		if ( block_indices[i] == (Jac->LCID(element->globalNodeId(node_col))) )
		  block = matrices[i];
	      }
	      TEST_FOR_EXCEPTION(block == 0, std::logic_error,"Failed to find block column index for this entry!");

	      (*block)(eq_row,eq_col) += val[eq_row](cell,node_row).fastAccessDx(lcol); 
	      
	    //   Jac->SumIntoGlobalValues(row, 1, 
// 				       &(val[eq_row](cell,node_row).fastAccessDx(lcol)),
// 				       &col);
	      
	    } // column equations
	    
	  } // column nodes
	  
	} // has fast access
	
      } // row equations
      
    } // row node

  } // element

}

// **********************************************************************
// Specialization: Jv
// **********************************************************************

template<typename Traits>
ScatterResidual<PHX::MyTraits::Jv, Traits>::
ScatterResidual(const Teuchos::ParameterList& p)
{ 
  // Not implemented yet!
}

// **********************************************************************
template<typename Traits> 
void ScatterResidual<PHX::MyTraits::Jv, Traits>::
postRegistrationSetup(typename Traits::SetupData data,
		      PHX::FieldManager<Traits>& fm)
{
  // Not implemented yet!
}

// **********************************************************************
template<typename Traits>
void ScatterResidual<PHX::MyTraits::Jv, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Error: GatherSolution not implemented for \"Jv\" evaluation type!");
}

// **********************************************************************
