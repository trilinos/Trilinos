// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_FUNCTIONAL_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_FUNCTIONAL_HPP

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Response_Functional.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"

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
   FunctionalScatter(const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer)
   {
     if(globalIndexer!=Teuchos::null)
       ugis_.push_back(globalIndexer);
   }

   FunctionalScatter(const std::vector<Teuchos::RCP<const panzer::GlobalIndexer> > & ugis)
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

   std::vector<Teuchos::RCP<const panzer::GlobalIndexer> > ugis_;
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
  // for convenience pull out some objects from workset
  std::string blockId = wda(workset).block_id;

  std::vector<int> blockOffsets;
  computeBlockOffsets(blockId,ugis_,blockOffsets);

  TEUCHOS_ASSERT(dgdx.size()==ugis_.size());

  auto cellIntegral_h = Kokkos::create_mirror_view(cellIntegral.get_view());
  Kokkos::deep_copy(cellIntegral_h, cellIntegral.get_view());

  const std::vector<std::size_t> & localCellIds = wda(workset).cell_local_ids;

  for(std::size_t b=0;b<ugis_.size();b++) {
    int start = blockOffsets[b];

    auto LIDs = ugis_[b]->getLIDs();
    auto LIDs_h = Kokkos::create_mirror_view(LIDs);
    Kokkos::deep_copy(LIDs_h, LIDs);

    Teuchos::ArrayRCP<double> dgdx_b = dgdx[b];

    // scatter operation for each cell in workset
    for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      // loop over basis functions
      for(std::size_t i=0;i<LIDs_h.extent(1);i++) {
        dgdx_b[LIDs_h(cellLocalId, i)] += cellIntegral_h(worksetCellIndex).dx(start+i); // its possible functional is independent of solution value!
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
  PHX::View<const LO*> LIDs;

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
