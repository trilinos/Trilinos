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
    PHX::MDField<ScalarT,PHX::NaturalOrder,Cell,Node> mdf(names[eq],dl);
    val[eq] = mdf;
    this->addDependentField(val[eq]);
  }

  this->addEvaluatedField(*scatter_operation);

  this->setName("Scatter Residual");
}

// **********************************************************************
template<typename Traits> 
void ScatterResidual<PHX::MyTraits::Residual, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& fm)
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
    
    for (std::size_t node = 0; node < num_nodes; node++) {
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
    PHX::MDField<ScalarT,PHX::NaturalOrder,Cell,Node> mdf(names[eq],dl);
    val[eq] = mdf;
    this->addDependentField(val[eq]);
  }

  this->addEvaluatedField(*scatter_operation);

  this->setName("Scatter Residual(Jacobian)");
}

// **********************************************************************
template<typename Traits> 
void ScatterResidual<PHX::MyTraits::Jacobian, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& fm)
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
    unsigned int lrow, lcol;
    for (unsigned int node_row = 0; node_row < num_nodes; node_row++) {
      
      int row_dim;
      int num_block_entries;
      int* block_indices;
      Jac->BeginExtractMyBlockRowView(Jac->LRID(element->globalNodeId(node_row)),row_dim, num_block_entries,block_indices);
      std::vector<Epetra_SerialDenseMatrix*> matrices(num_block_entries); 
      for (std::size_t i = 0; i < matrices.size();  ++i)
	Jac->ExtractEntryView(matrices[i]);

      // Loop over equations per node
      for (unsigned int eq_row = 0; eq_row < num_eq; eq_row++) {
	
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
	  for (unsigned int node_col=0; node_col<num_nodes; node_col++){
	    
	    // Loop over equations per node
	    for (unsigned int eq_col=0; eq_col<num_eq; eq_col++) {
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
