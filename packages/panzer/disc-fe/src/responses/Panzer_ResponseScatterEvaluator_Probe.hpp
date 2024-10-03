// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ResponseScatterEvaluator_Probe_hpp__
#define __Panzer_ResponseScatterEvaluator_Probe_hpp__

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Response_Probe.hpp"
#include "Panzer_GlobalIndexer.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

class ProbeScatterBase {
public:
  virtual ~ProbeScatterBase() {}

  virtual void scatterDerivative(
    const panzer::Traits::Jacobian::ScalarT& probeValue,
    const size_t cell_index,
    const bool has_probe,
    panzer::Traits::EvalData workset,
    WorksetDetailsAccessor& wda,
    Teuchos::ArrayRCP<double> & dgdx) const = 0;
};

template <typename LO,typename GO>
class ProbeScatter : public ProbeScatterBase {
public:
   ProbeScatter(const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer)
     : globalIndexer_(globalIndexer) { }

   void scatterDerivative(
     const panzer::Traits::Jacobian::ScalarT& probeValue,
     const size_t cell_index,
     const bool has_probe,
     panzer::Traits::EvalData workset,
     WorksetDetailsAccessor& wda,
     Teuchos::ArrayRCP<double> & dgdx) const;

private:

   Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
};

/** This class handles calculation of a DOF at a single point in space
 */
template<typename EvalT, typename Traits, typename LO, typename GO>
class ResponseScatterEvaluator_ProbeBase :
    public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>  {
public:

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator_ProbeBase(
    const std::string & responseName,
    const std::string & fieldName,
    const int fieldComponent,
    const Teuchos::Array<double>& point,
    const IntegrationRule & ir,
    const Teuchos::RCP<const PureBasis>& basis,
    const Teuchos::RCP<const panzer::GlobalIndexer>& indexer,
    const Teuchos::RCP<ProbeScatterBase> & probeScatter);

  void evaluateFields(typename Traits::EvalData d);

  void postRegistrationSetup(typename Traits::SetupData,
                             PHX::FieldManager<Traits>&);

  void preEvaluate(typename Traits::PreEvalData d);

  // Should be protected, but is public for cuda lambda support
  bool findCellAndComputeBasisValues(typename Traits::EvalData d);

protected:
  typedef typename EvalT::ScalarT ScalarT;

  std::string responseName_;
  std::string fieldName_;
  int fieldComponent_;
  Teuchos::Array<double> point_;
  Teuchos::RCP<const panzer::PureBasis> basis_;
  Teuchos::RCP<Response_Probe<EvalT> > responseObj_;
  Teuchos::RCP<const shards::CellTopology> topology_;
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;

  Teuchos::RCP<PHX::FieldTag> scatterHolder_; // dummy target
  PHX::MDField<const ScalarT,Cell,BASIS> field_; // holds field values
  Teuchos::RCP<ProbeScatterBase> scatterObj_;

  bool haveProbe_;
  int cellIndex_;
  size_t workset_id_;
  size_t num_basis, num_dim;
  Kokkos::DynRankView<double,PHX::Device> basis_values_;
};

/** This class handles calculation of a DOF at a single point in space
 */
template<typename EvalT, typename Traits, typename LO, typename GO>
class ResponseScatterEvaluator_Probe :
    public ResponseScatterEvaluator_ProbeBase<EvalT,Traits,LO,GO>  {
public:

  typedef ResponseScatterEvaluator_ProbeBase<EvalT,Traits,LO,GO> Base;

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator_Probe(
    const std::string & responseName,
    const std::string & fieldName,
    const int fieldComponent,
    const Teuchos::Array<double>& point,
    const IntegrationRule & ir,
    const Teuchos::RCP<const PureBasis>& basis,
    const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
    const Teuchos::RCP<ProbeScatterBase> & probeScatter) :
    Base(responseName, fieldName, fieldComponent, point,
         ir, basis, indexer, probeScatter) {}
};

/** This class handles calculation of a DOF at a single point in space
 */
template<typename LO, typename GO>
class ResponseScatterEvaluator_Probe<panzer::Traits::Jacobian,panzer::Traits,LO,GO> :
    public ResponseScatterEvaluator_ProbeBase<panzer::Traits::Jacobian,panzer::Traits,LO,GO>  {
public:

  typedef ResponseScatterEvaluator_ProbeBase<panzer::Traits::Jacobian,panzer::Traits,LO,GO> Base;

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator_Probe(
    const std::string & responseName,
    const std::string & fieldName,
    const int fieldComponent,
    const Teuchos::Array<double>& point,
    const IntegrationRule & ir,
    const Teuchos::RCP<const PureBasis>& basis,
    const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
    const Teuchos::RCP<ProbeScatterBase> & probeScatter) :
    Base(responseName, fieldName, fieldComponent, point,
         ir, basis, indexer, probeScatter) {}

  void evaluateFields(typename panzer::Traits::EvalData d);
};

template <typename LO,typename GO>
void ProbeScatter<LO,GO>::scatterDerivative(
  const panzer::Traits::Jacobian::ScalarT& probeValue,
  const size_t cell_index,
  const bool has_probe,
  panzer::Traits::EvalData workset,
  WorksetDetailsAccessor& wda,
  Teuchos::ArrayRCP<double> & dgdx) const
{

  if (has_probe) {
    PHX::View<const LO*> LIDs = globalIndexer_->getElementLIDs(cell_index);

    // loop over basis functions
    for(std::size_t i=0; i<LIDs.size(); ++i) {
      dgdx[LIDs[i]] += probeValue.dx(i);
    }
  }
}

}

#endif
