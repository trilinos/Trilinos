// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_VP_NONLINEAR_SOURCE_DEF_HPP
#define PHX_EXAMPLE_VP_NONLINEAR_SOURCE_DEF_HPP

//**********************************************************************
template<typename EvalT, typename Traits> NonlinearSource<EvalT, Traits>::
NonlinearSource(const Teuchos::ParameterList& p) :
  source("Nonlinear Source", 
	 p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout")),
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout")),
  temp("Temperature", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout"))
{ 
  this->addEvaluatedField(source);
  this->addDependentField(density);
  this->addDependentField(temp);

  this->setName("NonlinearSource");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData /* d */,
		      PHX::FieldManager<Traits>& vm)
{
  // NOTE: We no longer need manually call setFieldData(). It happens
  // automatically. We set them here to make sure backwards
  // compatiblity still works. We demonstrate both valid calling
  // pathways.
  this->utils.setFieldData(source,vm);
  vm.template getFieldData<EvalT>(density);
  vm.template getFieldData<EvalT>(temp);

  cell_data_size = source.size() / source.dimension(0);
}
//*********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void NonlinearSource<EvalT, Traits>::operator () (const int i) const
{
  for (PHX::index_size_type ip = 0; ip < static_cast<PHX::index_size_type>(density.extent(1)); ++ip)
    source(i,ip) =  density(i,ip) * temp(i,ip) * temp(i,ip);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
evaluateFields(typename Traits::EvalData d)
{ 
 Kokkos::parallel_for(d.num_cells, *this);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
preEvaluate(typename Traits::PreEvalData /* d */)
{ 
  using namespace std;
  cout << "In Source Pre Op" << endl;
}

//**********************************************************************
template< typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
postEvaluate(typename Traits::PostEvalData /* d */)
{ 
  using namespace std;
  cout << "In Source Post Op" << endl;
}

//**********************************************************************

#endif
