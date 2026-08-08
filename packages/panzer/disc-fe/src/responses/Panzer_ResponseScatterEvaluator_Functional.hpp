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

/**
 * \brief Abstract interface for scattering a functional response's derivative back into global sensitivity vectors.
 *
 * ResponseScatterEvaluator_Functional computes a scalar functional
 * integral per cell; a FunctionalScatterBase implementation is
 * responsible for taking that per-cell integral's derivative (with
 * respect to local solution unknowns) and accumulating it into the
 * appropriate entries of the response's global gradient (`dgdx`) or,
 * if enabled, Hessian-vector product (`d2gdx2`).
 */
class FunctionalScatterBase {
public:
  virtual ~FunctionalScatterBase() {}

  /**
   * \brief Accumulates dg/dx contributions from one workset's cell integrals into the global gradient vector(s).
   * \param cellIntegral per-cell functional integral for this workset, carrying Jacobian (first-derivative) sensitivities.
   * \param workset the workset being processed.
   * \param wda accessor used to pull cell/block details out of \p workset.
   * \param dgdx one global gradient vector per global indexer (block), to be accumulated into.
   */
  virtual void scatterDerivative(const PHX::MDField<const panzer::Traits::Jacobian::ScalarT,panzer::Cell> & cellIntegral,
                                 panzer::Traits::EvalData workset,
                                 WorksetDetailsAccessor& wda,
                                 const std::vector<Teuchos::ArrayRCP<double> > & dgdx) const = 0;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  /**
   * \brief Accumulates d2g/dx2 contributions from one workset's cell integrals into the global Hessian-vector product(s). Only defined when built with Hessian support.
   * \param cellIntegral per-cell functional integral for this workset, carrying Hessian (second-derivative) sensitivities.
   * \param workset the workset being processed.
   * \param wda accessor used to pull cell/block details out of \p workset.
   * \param d2gdx2 one global Hessian-vector product vector per global indexer (block), to be accumulated into.
   */
  virtual void scatterHessian(const PHX::MDField<const panzer::Traits::Hessian::ScalarT,panzer::Cell> & cellIntegral,
                              panzer::Traits::EvalData workset,
                              WorksetDetailsAccessor& wda,
                              const std::vector<Teuchos::ArrayRCP<double> > & d2gdx2) const = 0;
#endif
};

/**
 * \brief Concrete FunctionalScatterBase that scatters using one or more panzer::GlobalIndexer objects.
 * \tparam LO local ordinal type.
 * \tparam GO global ordinal type.
 */
template <typename LO,typename GO>
class FunctionalScatter : public FunctionalScatterBase {
public:
   /// \brief Construct for a single global indexer (single-block problem); no-op if \p globalIndexer is null.
   FunctionalScatter(const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer)
   {
     if(globalIndexer!=Teuchos::null)
       ugis_.push_back(globalIndexer);
   }

   /// \brief Construct for a blocked problem with one global indexer per block.
   FunctionalScatter(const std::vector<Teuchos::RCP<const panzer::GlobalIndexer> > & ugis)
     : ugis_(ugis) {}

   /// \copydoc FunctionalScatterBase::scatterDerivative
   void scatterDerivative(const PHX::MDField<const panzer::Traits::Jacobian::ScalarT,panzer::Cell> & cellIntegral,
                         panzer::Traits::EvalData workset,
                         WorksetDetailsAccessor& wda,
                         const std::vector<Teuchos::ArrayRCP<double> > & dgdx) const;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
   /// \copydoc FunctionalScatterBase::scatterHessian
   void scatterHessian(const PHX::MDField<const panzer::Traits::Hessian::ScalarT,panzer::Cell> & cellIntegral,
                       panzer::Traits::EvalData workset,
                       WorksetDetailsAccessor& wda,
                       const std::vector<Teuchos::ArrayRCP<double> > & d2gdx2) const;
#endif

private:

   /// One global indexer per block; a single entry for non-blocked problems.
   std::vector<Teuchos::RCP<const panzer::GlobalIndexer> > ugis_;
};

/** This class handles responses with values aggregated
  * on each finite element cell.
  *
  * It depends on a "Cell Integral" field (an integrand already
  * integrated over each cell) having been computed upstream in the
  * evaluation DAG under the integrand name; evaluateFields() sums that
  * field into the associated Response_Functional, and for the
  * Jacobian/Hessian evaluation types, hands its derivatives off to a
  * FunctionalScatterBase to accumulate into the response's global
  * gradient / Hessian-vector product.
  */
template<typename EvalT, typename Traits>
class ResponseScatterEvaluator_Functional : public panzer::EvaluatorWithBaseImpl<Traits>,
                                            public PHX::EvaluatorDerived<EvalT, Traits>  {
public:

  //! A constructor with concrete arguments instead of a parameter list.
  /** \brief Construct where the response name and the upstream integrand field name are the same.
    *
    * \param[in] name name of both the upstream "Cell Integral" field to sum and the response it feeds.
    * \param[in] cd cell data (dimension, workset size) for the field data layout.
    * \param[in] functionalScatter used to scatter derivative/Hessian information into the response's global vectors.
    */
  ResponseScatterEvaluator_Functional(const std::string & name,const CellData & cd,
                                      const Teuchos::RCP<FunctionalScatterBase> & functionalScatter);
  /** \brief Construct with distinct names for the upstream integrand field and the response it feeds.
    *
    * \param[in] integrandName name of the upstream "Cell Integral" field to sum.
    * \param[in] responseName name of the response this evaluator contributes to.
    * \param[in] cd cell data (dimension, workset size) for the field data layout.
    * \param[in] functionalScatter used to scatter derivative/Hessian information into the response's global vectors.
    */
  ResponseScatterEvaluator_Functional(const std::string & integrandName,const std::string & responseName,const CellData & cd,
                                      const Teuchos::RCP<FunctionalScatterBase> & functionalScatter);

  /// \brief Sums this workset's per-cell integrand values (and scatters derivatives, via #scatterObj_) into the response.
  void evaluateFields(typename Traits::EvalData d);

  /// \brief Looks up the Response_Functional to accumulate into for the upcoming fill.
  void preEvaluate(typename Traits::PreEvalData d);

private:
  typedef typename EvalT::ScalarT ScalarT;

  /// Name of the response this evaluator contributes to.
  std::string responseName_;
  /// The response object accumulated into by evaluateFields().
  Teuchos::RCP<Response_Functional<EvalT> > responseObj_;

  Teuchos::RCP<PHX::FieldTag> scatterHolder_; // dummy target
  PHX::MDField<const ScalarT,panzer::Cell> cellIntegral_; // holds cell integrals
  /// Used to scatter derivative/Hessian information into the response's global vectors. \sa FunctionalScatterBase
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
