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

#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_FUNCTIONAL_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_FUNCTIONAL_HPP

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Response_Functional.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

class FunctionalScatterBase {
public:
  virtual ~FunctionalScatterBase() {}

  virtual void scatterDerivative(const PHX::MDField<const panzer::Traits::Jacobian::ScalarT,panzer::Cell> & cellIntegral,
                                 panzer::Traits::EvalData workset, 
                                 WorksetDetailsAccessor& wda,
                                 const std::vector<Teuchos::ArrayRCP<double> > & dgdx) const = 0;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  virtual void scatterHessian(const PHX::MDField<const panzer::Traits::Hessian::ScalarT,panzer::Cell> & cellIntegral,
                              panzer::Traits::EvalData workset, 
                              WorksetDetailsAccessor& wda,
                              const std::vector<Teuchos::ArrayRCP<double> > & d2gdx2) const = 0;
#endif
};
 
template <typename LO,typename GO>
class FunctionalScatter : public FunctionalScatterBase {
public:
   FunctionalScatter(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & globalIndexer)
   { 
     if(globalIndexer!=Teuchos::null)
       ugis_.push_back(globalIndexer);
   }

   FunctionalScatter(const std::vector<Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > > & ugis)
     : ugis_(ugis) {}

   void scatterDerivative(const PHX::MDField<const panzer::Traits::Jacobian::ScalarT,panzer::Cell> & cellIntegral,
                         panzer::Traits::EvalData workset, 
                         WorksetDetailsAccessor& wda,
                         const std::vector<Teuchos::ArrayRCP<double> > & dgdx) const;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
   void scatterHessian(const PHX::MDField<const panzer::Traits::Hessian::ScalarT,panzer::Cell> & cellIntegral,
                       panzer::Traits::EvalData workset, 
                       WorksetDetailsAccessor& wda,
                       const std::vector<Teuchos::ArrayRCP<double> > & d2gdx2) const;
#endif

private:
 
   std::vector<Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > > ugis_;
};

/** This class handles responses with values aggregated
  * on each finite element cell.
  */
template<typename EvalT, typename Traits>
class ResponseScatterEvaluator_Functional : public panzer::EvaluatorWithBaseImpl<Traits>,
                                            public PHX::EvaluatorDerived<EvalT, Traits>  { 
public:

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator_Functional(const std::string & name,const CellData & cd,
                                      const Teuchos::RCP<FunctionalScatterBase> & functionalScatter);
  ResponseScatterEvaluator_Functional(const std::string & integrandName,const std::string & responseName,const CellData & cd,
                                      const Teuchos::RCP<FunctionalScatterBase> & functionalScatter);

  void evaluateFields(typename Traits::EvalData d);

  void preEvaluate(typename Traits::PreEvalData d);

private:
  typedef typename EvalT::ScalarT ScalarT;

  std::string responseName_;
  Teuchos::RCP<Response_Functional<EvalT> > responseObj_;

  Teuchos::RCP<PHX::FieldTag> scatterHolder_; // dummy target
  PHX::MDField<const ScalarT,panzer::Cell> cellIntegral_; // holds cell integrals
  Teuchos::RCP<FunctionalScatterBase> scatterObj_;
};

template <typename LO,typename GO>
void FunctionalScatter<LO,GO>::scatterDerivative(const PHX::MDField<const panzer::Traits::Jacobian::ScalarT,panzer::Cell> & cellIntegral,
                                                panzer::Traits::EvalData workset, 
                                                WorksetDetailsAccessor& wda,
                                                const std::vector<Teuchos::ArrayRCP<double> > & dgdx) const 
{
  Kokkos::View<const LO*, PHX::Device> LIDs;
 
  // for convenience pull out some objects from workset
  std::string blockId = wda(workset).block_id;

  std::vector<int> blockOffsets;
  computeBlockOffsets(blockId,ugis_,blockOffsets);

  TEUCHOS_ASSERT(dgdx.size()==ugis_.size());

  // scatter operation for each cell in workset
  const std::vector<std::size_t> & localCellIds = wda(workset).cell_local_ids;
  for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
    std::size_t cellLocalId = localCellIds[worksetCellIndex];

    for(std::size_t b=0;b<ugis_.size();b++) {
      int start = blockOffsets[b];

      LIDs = ugis_[b]->getElementLIDs(cellLocalId); 

      Teuchos::ArrayRCP<double> dgdx_b = dgdx[b];

      // loop over basis functions
      for(std::size_t i=0;i<LIDs.size();i++) {
        dgdx_b[LIDs[i]] += cellIntegral(worksetCellIndex).dx(start+i); // its possible functional is independent of solution value!
      }
    }
  }
}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
template <typename LO,typename GO>
void FunctionalScatter<LO,GO>::scatterHessian(const PHX::MDField<const panzer::Traits::Hessian::ScalarT,panzer::Cell> & cellIntegral,
                                                panzer::Traits::EvalData workset, 
                                                WorksetDetailsAccessor& wda,
                                                const std::vector<Teuchos::ArrayRCP<double> > & d2gdx2) const 
{
  Kokkos::View<const LO*, PHX::Device> LIDs;
 
  // for convenience pull out some objects from workset
  std::string blockId = wda(workset).block_id;

  std::vector<int> blockOffsets;
  computeBlockOffsets(blockId,ugis_,blockOffsets);

  TEUCHOS_ASSERT(d2gdx2.size()==ugis_.size());

  // scatter operation for each cell in workset
  const std::vector<std::size_t> & localCellIds = wda(workset).cell_local_ids;
  for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
    std::size_t cellLocalId = localCellIds[worksetCellIndex];

    for(std::size_t b=0;b<ugis_.size();b++) {
      int start = blockOffsets[b];

      LIDs = ugis_[b]->getElementLIDs(cellLocalId); 

      Teuchos::ArrayRCP<double> d2gdx2_b = d2gdx2[b];

      // loop over basis functions
      for(std::size_t i=0;i<LIDs.size();i++) {
        d2gdx2_b[LIDs[i]] += cellIntegral(worksetCellIndex).dx(start+i).dx(0); // its possible functional is independent of solution value!
      }
    }
  }
}
#endif

}

#endif
