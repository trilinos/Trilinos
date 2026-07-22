// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ScatterResidual_Tpetra_Hessian_impl_hpp__
#define __Panzer_ScatterResidual_Tpetra_Hessian_impl_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// Tpetra scatter residual file

namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
ScatterResidual_Tpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>::
ScatterResidual_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & /* indexer */,
                       const Teuchos::ParameterList& p) 
{
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ = 
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  this->addEvaluatedField(*scatterHolder_);

  this->setName(scatterName+" Scatter Residual Tpetra (Hessian)");
}
  
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void
ScatterResidual_Tpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* vm */) 
{
}

template<typename TRAITS,typename LO,typename GO,typename NodeT>
void
ScatterResidual_Tpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData /* d */) 
{
}
  
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void
ScatterResidual_Tpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData /* workset */) 
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                             "ScatterResidual_Tpetra<Hessian> is not yet implemented"); // just in case
}
  
}

#endif // end hessian support

#endif
