// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename Traits>
ScatterResidual<PHX::MyTraits::Residual, Traits>::
ScatterResidual(const Teuchos::RCP<PHX::FieldTag>& in_scatter_tag,
                const std::string& residual_name,
                const Teuchos::RCP<PHX::DataLayout>& residual_layout,
                const int& in_equation_index,
                const int& in_num_equations,
                const Kokkos::View<double*,PHX::Device>& in_global_residual) :
  scatter_tag(in_scatter_tag),
  residual_contribution(residual_name,residual_layout),
  global_residual_atomic(in_global_residual),
  equation_index(in_equation_index),
  num_equations(in_num_equations)
{ 
  this->addEvaluatedField(*scatter_tag);
  this->addDependentField(residual_contribution);
  this->setName("Scatter Residual: "+residual_name);
}

// **********************************************************************
template<typename Traits>
void ScatterResidual<PHX::MyTraits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  gids = workset.gids_;
  cell_global_offset_index = workset.first_cell_global_index_;
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),*this);
}

// **********************************************************************
template<typename Traits>
PHALANX_HIP_HACK_KOKKOS_FUNCTION
void ScatterResidual<PHX::MyTraits::Residual,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int local_cell = team.league_rank();
  if (team.team_rank() == 0) {
    // Fix gcc 5/6 lambda bug by changing to capture by value (potentially less efficient)
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,residual_contribution.extent(1)), [=] (const int& node) {
      const int residual_index = gids(cell_global_offset_index+local_cell,node) * num_equations + equation_index;
      global_residual_atomic(residual_index) += residual_contribution(local_cell,node);
    });
  }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename Traits>
ScatterResidual<PHX::MyTraits::Jacobian, Traits>::
ScatterResidual(const Teuchos::RCP<PHX::FieldTag>& in_scatter_tag,
                const std::string& residual_name,
                const Teuchos::RCP<PHX::DataLayout>& residual_layout,
                const int& in_equation_index,
                const int& in_num_equations,
                const Kokkos::View<double*,PHX::Device>& in_global_residual,
                const KokkosSparse::CrsMatrix<double,int,PHX::Device>& in_global_jacobian) :
  scatter_tag(in_scatter_tag),
  residual_contribution(residual_name,residual_layout),
  global_residual_atomic(in_global_residual),
  global_jacobian(in_global_jacobian),
  equation_index(in_equation_index),
  num_equations(in_num_equations)
{ 
  this->addEvaluatedField(*scatter_tag);
  this->addDependentField(residual_contribution);
  this->setName("Scatter Jacobian: "+residual_name);
}

// **********************************************************************
template<typename Traits>
void ScatterResidual<PHX::MyTraits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  gids = workset.gids_;
  cell_global_offset_index = workset.first_cell_global_index_;
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),*this);
}

// **********************************************************************
template<typename Traits>
PHALANX_HIP_HACK_KOKKOS_FUNCTION
void ScatterResidual<PHX::MyTraits::Jacobian,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int cell = team.league_rank();
  const int num_nodes = residual_contribution.extent(1);

  if (team.team_rank() == 0) {
    // Fix gcc 5/6 lambda bug by changing to capture by value (potentially less efficient)
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,num_nodes), [=] (const int& node) {
      const int global_row_index = gids(cell_global_offset_index+cell,node) * num_equations + equation_index;
      global_residual_atomic(global_row_index) += residual_contribution(cell,node).val();
    });
  }

  // Fix gcc 5/6 lambda bug by changing to capture by value (potentially less efficient)
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_nodes), [=] (const int& node) {

    const int global_row_index = gids(cell_global_offset_index+cell,node) * num_equations + equation_index;

    // loop over nodes
    for (int col_node=0; col_node < num_nodes; ++col_node) {

      // loop over equations
      // Fix gcc 5/6 lambda bug by changing to capture by value (potentially less efficient)
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,num_equations),[=] (const int& col_eq) {
        const int global_col_index = gids(cell_global_offset_index+cell,col_node) * num_equations + col_eq;
        const int derivative_index = col_node * num_equations + col_eq;
        global_jacobian.sumIntoValues(global_row_index,&global_col_index,1,
                                      &(residual_contribution(cell,node).fastAccessDx(derivative_index)),
                                      false,true);
      });

    } 
  });
}

  /*
  std::size_t cell = 0;
  for (; element != workset.end; ++element,++cell) {
    
    // Sum element residual and Jacobian into global residual, Jacobian
    // Loop over nodes in element
    for (int node = 0; node < num_nodes; node++) {
      
      int node_GID = element->globalNodeId(node);
      int firstDOF = Jac->RowMap().LID(node_GID * num_eq);

      // Loop over equations per node
      for (int eq = 0; eq < num_eq; eq++) {
	
	int row = firstDOF + eq;
	
	// Sum residual
	if (f != Teuchos::null)
	  f->SumIntoMyValue(row, 0, val[eq](cell,node).val());
	

	// Check derivative array is nonzero
	if (val[eq](cell,node).hasFastAccess()) {
	  
	  // Loop over nodes in element
	  int firstcol = -1;
	  for (int node_col=0; node_col<num_nodes; node_col++){
	    firstcol =  Jac->RowMap().LID(element->globalNodeId(node_col) * num_eq);
	    
	    // Loop over equations per node
	    for (int eq_col=0; eq_col<num_eq; eq_col++) {
	      
              int lcol = num_eq * node_col + eq_col;
	      int col = firstcol + eq_col;
	      
	      // Sum Jacobian
	      Jac->SumIntoMyValues(row, 1, 
				   &(val[eq](cell,node).fastAccessDx(lcol)),
				   &col);
	      
	    } // column equations
	    
	  } // column nodes
	  
	} // has fast access
	
      } // row equations
      
    } // row node

  } // element
  */

// **********************************************************************
// Specialization: Jv
// **********************************************************************
/*
template<typename Traits>
ScatterResidual<PHX::MyTraits::Jv, Traits>::
ScatterResidual(const Teuchos::ParameterList& p)
{ 
  scatter_operation = Teuchos::rcp(new PHX::Tag<ScalarT>("Scatter",p.get< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout")));

  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Residual Names"));

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  f = p.get< Teuchos::RCP<Epetra_Vector> >("Residual Vector");
  
  Jac = p.get< Teuchos::RCP<Epetra_CrsMatrix> >("Jacobian Matrix");
  
  val.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<ScalarT,Cell,Node> mdf(names[eq],dl);
    val[eq] = mdf;
    this->addDependentField(val[eq]);
  }

  this->addEvaluatedField(*scatter_operation);

  this->setName("Scatter Residual(Jv)");
}

// **********************************************************************
template<typename Traits> 
void ScatterResidual<PHX::MyTraits::Jv, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  for (std::size_t eq = 0; eq < val.size(); ++eq)
    this->utils.setFieldData(val[eq],fm);

  num_nodes = val[0].dimension(1);
  num_eq = val.size();
}

// **********************************************************************
template<typename Traits>
void ScatterResidual<PHX::MyTraits::Jv, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  std::vector<Element_Linear2D>::iterator element = workset.begin;

  std::size_t cell = 0;
  for (; element != workset.end; ++element,++cell) {
    
    // Sum element residual and Jacobian into global residual, Jacobian
    // Loop over nodes in element
    for (int node = 0; node < num_nodes; node++) {
      
      int node_GID = element->globalNodeId(node);
      int firstDOF = Jac->RowMap().LID(node_GID * num_eq);

      // Loop over equations per node
      for (int eq = 0; eq < num_eq; eq++) {
	
	int row = firstDOF + eq;
	
	// Sum residual
	// 	if (f != Teuchos::null)
	// 	  f->SumIntoMyValue(row, 0, val[eq](cell,node).val());
	
	workset.Jv->SumIntoMyValue(row, 0, 
				     val[eq](cell,node).fastAccessDx(0));
	
	

      } // row equations
      
    } // row node

  } // element

}
*/
// **********************************************************************
