// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_FUNCTIONAL_IMPL_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_FUNCTIONAL_IMPL_HPP

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_ResponseBase.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Teuchos_ArrayRCP.hpp"

namespace panzer {

/** This class handles responses with values aggregated
  * on each finite element cell.
  */
template<typename EvalT, typename Traits>
ResponseScatterEvaluator_Functional<EvalT,Traits>::
ResponseScatterEvaluator_Functional(const std::string & name,
                                    const CellData & cd,
                                    const Teuchos::RCP<FunctionalScatterBase> & functionalScatter)
  : responseName_(name)
  , scatterObj_(functionalScatter)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  std::string dummyName = ResponseBase::buildLookupName(name) + " dummy target";

  // build dummy target tag
  RCP<PHX::DataLayout> dl_dummy = rcp(new PHX::MDALayout<panzer::Dummy>(0));
  scatterHolder_ = rcp(new PHX::Tag<ScalarT>(dummyName,dl_dummy));
  this->addEvaluatedField(*scatterHolder_);

  // build dendent field
  RCP<PHX::DataLayout> dl_cell = rcp(new PHX::MDALayout<panzer::Cell>(cd.numCells()));
  cellIntegral_ = PHX::MDField<const ScalarT,panzer::Cell>(name,dl_cell);
  this->addDependentField(cellIntegral_);

  std::string n = "Functional Response Scatter: " + name;
  this->setName(n);
}

template<typename EvalT, typename Traits>
ResponseScatterEvaluator_Functional<EvalT,Traits>::
ResponseScatterEvaluator_Functional(const std::string & integrandName,
                                    const std::string & responseName,
                                    const CellData & cd,
                                    const Teuchos::RCP<FunctionalScatterBase> & functionalScatter)
  : responseName_(responseName)
  , scatterObj_(functionalScatter)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  std::string dummyName = ResponseBase::buildLookupName(responseName) + " dummy target";

  // build dummy target tag
  RCP<PHX::DataLayout> dl_dummy = rcp(new PHX::MDALayout<panzer::Dummy>(0));
  scatterHolder_ = rcp(new PHX::Tag<ScalarT>(dummyName,dl_dummy));
  this->addEvaluatedField(*scatterHolder_);

  // build dendent field
  RCP<PHX::DataLayout> dl_cell = rcp(new PHX::MDALayout<panzer::Cell>(cd.numCells()));
  cellIntegral_ = PHX::MDField<const ScalarT,panzer::Cell>(integrandName,dl_cell);
  this->addDependentField(cellIntegral_);

  std::string n = "Functional Response Scatter: " + responseName;
  this->setName(n);
}

template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_Functional<EvalT,Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  // extract linear object container
  responseObj_ = Teuchos::rcp_dynamic_cast<Response_Functional<EvalT> >(
                                   d.gedc->getDataObject(ResponseBase::buildLookupName(responseName_)),true);
}


template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_Functional<EvalT,Traits>::
evaluateFields(typename Traits::EvalData d)
{
  auto cellIntegral_h = Kokkos::create_mirror_view ( cellIntegral_.get_static_view());
  Kokkos::deep_copy(cellIntegral_h, cellIntegral_.get_static_view());
  for(index_t i=0;i<d.num_cells;i++) {
    responseObj_->value += cellIntegral_h(i);
  }
}

template < >
void ResponseScatterEvaluator_Functional<panzer::Traits::Jacobian,panzer::Traits>::
evaluateFields(panzer::Traits::EvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::SpmdVectorBase;
  using Thyra::ProductVectorBase;

  TEUCHOS_ASSERT(scatterObj_!=Teuchos::null);
  TEUCHOS_ASSERT(responseObj_->getGhostedVector()!=Teuchos::null);

  RCP<ProductVectorBase<double> > prod_dgdx = Thyra::castOrCreateNonconstProductVectorBase(responseObj_->getGhostedVector());

  std::vector<Teuchos::ArrayRCP<double> > local_dgdxs;
  for(int b=0;b<prod_dgdx->productSpace()->numBlocks();b++) {
    // grab local data for inputing
    Teuchos::ArrayRCP<double> local_dgdx;
    RCP<SpmdVectorBase<double> > dgdx = rcp_dynamic_cast<SpmdVectorBase<double> >(prod_dgdx->getNonconstVectorBlock(b));
    dgdx->getNonconstLocalData(ptrFromRef(local_dgdx));

    TEUCHOS_ASSERT(!local_dgdx.is_null());

    local_dgdxs.push_back(local_dgdx);
  }

  scatterObj_->scatterDerivative(cellIntegral_,d,this->wda,local_dgdxs);
}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
template < >
void ResponseScatterEvaluator_Functional<panzer::Traits::Hessian,panzer::Traits>::
evaluateFields(panzer::Traits::EvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::SpmdVectorBase;
  using Thyra::ProductVectorBase;

  TEUCHOS_ASSERT(scatterObj_!=Teuchos::null);
  TEUCHOS_ASSERT(responseObj_->getGhostedVector()!=Teuchos::null);

  RCP<ProductVectorBase<double> > prod_dgdx = Thyra::castOrCreateNonconstProductVectorBase(responseObj_->getGhostedVector());

  std::vector<Teuchos::ArrayRCP<double> > local_dgdxs;
  for(int b=0;b<prod_dgdx->productSpace()->numBlocks();b++) {
    // grab local data for inputing
    Teuchos::ArrayRCP<double> local_dgdx;
    RCP<SpmdVectorBase<double> > dgdx = rcp_dynamic_cast<SpmdVectorBase<double> >(prod_dgdx->getNonconstVectorBlock(b));
    dgdx->getNonconstLocalData(ptrFromRef(local_dgdx));

    TEUCHOS_ASSERT(!local_dgdx.is_null());

    local_dgdxs.push_back(local_dgdx);
  }

  // TEUCHOS_ASSERT(false);
  scatterObj_->scatterHessian(cellIntegral_,d,this->wda,local_dgdxs);
}
#endif

}

#endif
