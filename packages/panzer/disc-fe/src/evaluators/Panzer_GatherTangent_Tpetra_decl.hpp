// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_GATHER_TANGENT_TPETRA_DECL_HPP
#define PANZER_EVALUATOR_GATHER_TANGENT_TPETRA_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_KokkosViewOfViews.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"

#include"Panzer_NodeType.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

class GlobalIndexer; //forward declaration

/** \brief Gathers tangent vectors dx/dp for computing df/dx*dx/dp + df/dp into
    the nodal fields of the field manager.

    This evaluator is very similar to GatherSolution, however it always gathers
    into fields of type double, and it is a no-op if the global evaluation data
    container does not exist (which is an error for GatherSolution).

    Currently makes an assumption that the stride is constant for dofs
    and that the nmber of dofs is equal to the size of the solution
    names vector.
*/
template<typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT=panzer::TpetraNodeType>
class GatherTangent_Tpetra
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<EvalT, TRAITS>,
    public panzer::CloneableEvaluator  {
public:

  GatherTangent_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer) :
     globalIndexer_(indexer) {}

  GatherTangent_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                        const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherTangent_Tpetra<EvalT,TRAITS,LO,GO,NodeT>(globalIndexer_,pl)); }

  // for testing purposes
  const PHX::FieldTag & getFieldTag(int i) const
  { TEUCHOS_ASSERT(i < Teuchos::as<int>(gatherFields_.size())); return gatherFields_[i].fieldTag(); }

private:

  // We always use RealType for gathering as we never compute derivatives for this evaluator
  typedef typename panzer::Traits::RealType ScalarT;
  //typedef typename EvalT::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;
  PHX::ViewOfViews<1,PHX::View<ScalarT**>> gatherFieldsVoV_;

  Teuchos::RCP<std::vector<std::string> > indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const TpetraLinearObjContainer<double,LO,GO,NodeT> > tpetraContainer_;

  GatherTangent_Tpetra();
};

}

// **************************************************************
#endif
