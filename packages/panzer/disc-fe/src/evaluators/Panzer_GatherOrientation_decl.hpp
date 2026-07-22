// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_GATHER_ORIENTATION_DECL_HPP
#define PANZER_EVALUATOR_GATHER_ORIENTATION_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

class GlobalIndexer; //forward declaration

/** \brief Gathers orientations per field from the global indexer and
    stores them in the field manager.
*/
template<typename EvalT, typename TRAITS,typename LO,typename GO> 
class GatherOrientation
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<EvalT, TRAITS>,
    public panzer::CloneableEvaluator  {
   
public:

  GatherOrientation(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer)
  { indexers_.push_back(indexer); }

  GatherOrientation(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                        const Teuchos::ParameterList& p);

  GatherOrientation(const std::vector<Teuchos::RCP<const GlobalIndexer>> & indexers)
     : indexers_(indexers) {}

  GatherOrientation(const std::vector<Teuchos::RCP<const GlobalIndexer>> & indexers,
                    const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);
  
  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherOrientation<EvalT,TRAITS,LO,GO>(indexers_,pl)); }
  
private:

  typedef typename EvalT::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;

  std::vector<Teuchos::RCP<const GlobalIndexer>> indexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFieldOrientations_;

  Teuchos::RCP<std::vector<std::string> > indexerNames_;

  GatherOrientation();
};

}

// **************************************************************
#endif
