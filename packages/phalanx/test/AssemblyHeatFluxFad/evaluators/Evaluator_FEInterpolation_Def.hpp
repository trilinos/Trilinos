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

//**********************************************************************
template<typename EvalT, typename Traits>
FEInterpolation<EvalT, Traits>::
FEInterpolation(const Teuchos::ParameterList& p) :
  val_node(p.get<std::string>("Node Variable Name"), 
	   p.get< Teuchos::RCP<PHX::DataLayout> >("Node Data Layout") ),
  val_qp(p.get<std::string>("QP Variable Name"), 
	 p.get< Teuchos::RCP<PHX::DataLayout> >("QP Scalar Data Layout") ),
  val_grad_qp(p.get<std::string>("Gradient QP Variable Name"), 
	      p.get< Teuchos::RCP<PHX::DataLayout> >("QP Vector Data Layout") ),
  dummy("Dummy", 
	p.get< Teuchos::RCP<PHX::DataLayout> >("QP Vector Data Layout") )
{ 
  this->addDependentField(val_node);
  this->addEvaluatedField(val_qp);
  this->addEvaluatedField(val_grad_qp);

  // for unit testing
  this->addEvaluatedField(dummy);
  
  this->setName("FEInterpolation");
}

//**********************************************************************
template<typename EvalT, typename Traits> 
void FEInterpolation<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData /* d */,
		      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(val_node,fm);
  this->utils.setFieldData(val_qp,fm);
  this->utils.setFieldData(val_grad_qp,fm);
  this->utils.setFieldData(dummy,fm);

  num_nodes = static_cast<PHX::index_size_type>(val_node.extent(1));
  num_qp = static_cast<PHX::index_size_type>(val_grad_qp.extent(1));
  num_dim = static_cast<PHX::index_size_type>(val_grad_qp.extent(2));
}
//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void FEInterpolation<EvalT, Traits>::operator () (const int i) const
{  
  for (PHX::index_size_type qp = 0; qp < num_qp; ++qp) {
    val_qp(i,qp) = 0.0;

    for (PHX::index_size_type dim = 0; dim < num_dim; ++dim)
      val_grad_qp(i,qp,dim) = 0.0;

    // Sum nodal contributions to qp
    for (PHX::index_size_type node = 0; node < num_nodes; ++node) {
      val_qp(i,qp) += phi(qp, node) * val_node(i,node);
      for (PHX::index_size_type dim = 0; dim < num_dim; ++dim)
	val_grad_qp(i,qp,dim) += grad_phi(qp, node, dim) * val_node(i,node);
    }
  }
}
//**********************************************************************
template<typename EvalT, typename Traits>
void FEInterpolation<EvalT, Traits>::
evaluateFields(typename Traits::EvalData cell_data)
{ 
  std::vector<MyCell>::iterator cell_it = cell_data.begin;
  phi = cell_it->getBasisFunctions();
  grad_phi = cell_it->getBasisFunctionGradients();
  Kokkos::parallel_for(cell_data.num_cells, *this);  
}
//**********************************************************************
#ifdef  PHX_ENABLE_KOKKOS_AMT
template<typename EvalT, typename Traits>
Kokkos::Future<void,PHX::exec_space>
FEInterpolation<EvalT, Traits>::
createTask(Kokkos::TaskScheduler<PHX::exec_space>& policy,
	   const int& work_size,
           const std::vector<Kokkos::Future<void,PHX::exec_space>>& dependent_futures,
	   typename Traits::EvalData workset)
{
  std::vector<MyCell>::iterator cell_it = workset.begin;
  phi = cell_it->getBasisFunctions();
  grad_phi = cell_it->getBasisFunctionGradients();
  auto dep_future = policy.when_all(dependent_futures.size(),dependent_futures.data());
  return policy.host_spawn(PHX::TaskWrap<PHX::exec_space,FEInterpolation<EvalT, Traits>>(work_size,*this),Kokkos::TaskTeam,dep_future);
}
#endif

//**********************************************************************
