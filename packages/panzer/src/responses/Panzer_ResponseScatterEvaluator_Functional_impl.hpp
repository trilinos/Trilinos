// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_FUNCTIONAL_IMPL_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_FUNCTIONAL_IMPL_HPP

#include <iostream>
#include <string>

#include "Panzer_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_ResponseBase.hpp"
#include "Panzer_Dimension.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace panzer {

template <typename LO,typename GO>
void FunctionalScatter<LO,GO>::scatterDerivative(const PHX::MDField<panzer::Traits::Jacobian::ScalarT,panzer::Cell> & cellIntegral,
                                                panzer::Traits::EvalData workset, 
                                                Teuchos::ArrayRCP<double> & dgdx) const
{
  std::vector<int> LIDs;
 
  // for convenience pull out some objects from workset
  std::string blockId = workset.block_id;
  const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

  // NOTE: A reordering of these loops will likely improve performance
  //       The "getGIDFieldOffsets may be expensive.  However the
  //       "getElementGIDs" can be cheaper. However the lookup for LIDs
  //       may be more expensive!

  // scatter operation for each cell in workset
  for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
    std::size_t cellLocalId = localCellIds[worksetCellIndex];
    LIDs = globalIndexer_->getElementLIDs(cellLocalId); 

    // loop over basis functions
    for(std::size_t i=0;i<LIDs.size();i++)
      dgdx[LIDs[i]] += cellIntegral(worksetCellIndex).dx(i); // its possible functional is independent of solution value!
  }
}

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
  cellIntegral_ = PHX::MDField<ScalarT,panzer::Cell>(name,dl_cell);
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
  cellIntegral_ = PHX::MDField<ScalarT,panzer::Cell>(integrandName,dl_cell);
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
                                   d.getDataObject(ResponseBase::buildLookupName(responseName_)),true);
}


template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_Functional<EvalT,Traits>::
postRegistrationSetup(typename Traits::SetupData d,
                      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(cellIntegral_,fm);
}

template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_Functional<EvalT,Traits>::
evaluateFields(typename Traits::EvalData d)
{
  for(std::size_t i=0;i<d.num_cells;i++) {
    responseObj_->value += cellIntegral_(i);
  }
}

template < >
void ResponseScatterEvaluator_Functional<panzer::Traits::Jacobian,panzer::Traits>::
evaluateFields(panzer::Traits::EvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::SpmdVectorBase;

  TEUCHOS_ASSERT(scatterObj_!=Teuchos::null);

  // grab local data for inputing
  Teuchos::ArrayRCP<double> local_dgdx;
  RCP<SpmdVectorBase<double> > dgdx = rcp_dynamic_cast<SpmdVectorBase<double> >(responseObj_->getGhostedVector());
  dgdx->getNonconstLocalData(ptrFromRef(local_dgdx));
  TEUCHOS_ASSERT(!local_dgdx.is_null());

  scatterObj_->scatterDerivative(cellIntegral_,d,local_dgdx);
}

}

#endif
