// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ReorderADValues_Evaluator_decl_hpp__
#define __Panzer_ReorderADValues_Evaluator_decl_hpp__

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

class GlobalIndexer;

/** \brief Reorders the ad values of a specified field to match a different
           unique global indexer.

    This is neccessary primarily when gathering with one unique global indexer    
    and scattering with a second unique global indexer.
*/
template<typename EvalT, typename TRAITS> 
class ReorderADValues_Evaluator  
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<EvalT, TRAITS> {
public:
  typedef typename EvalT::ScalarT ScalarT;

  ReorderADValues_Evaluator(const std::string & outPrefix,
                            const std::vector<std::string> & inFieldNames,
                            const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                            const std::string & elementBlock,
                            const GlobalIndexer & indexerSrc,
                            const GlobalIndexer & indexerDest);

  ReorderADValues_Evaluator(const std::string & outPrefix,
                            const std::vector<std::string> & inFieldNames,
                            const std::vector<std::string> & inDOFs,
                            const std::vector<std::string> & outDOFs,
                            const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                            const std::string & elementBlock,
                            const GlobalIndexer & indexerSrc,
                            const GlobalIndexer & indexerDest);

  void evaluateFields(typename TRAITS::EvalData d);

private:
  // fields to be modified
  std::vector< PHX::MDField<const ScalarT> > inFields_;

  // fields that need to be modified
  std::vector< PHX::MDField<ScalarT> > outFields_;
};

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Jacobian 
// **************************************************************
template<typename TRAITS>
class ReorderADValues_Evaluator<typename TRAITS::Jacobian,TRAITS>  
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<typename TRAITS::Jacobian, TRAITS> {
  
public:

  ReorderADValues_Evaluator(const std::string & outPrefix,
                            const std::vector<std::string> & inFieldNames,
                            const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                            const std::string & elementBlock,
                            const GlobalIndexer & indexerSrc,
                            const GlobalIndexer & indexerDest);

  ReorderADValues_Evaluator(const std::string & outPrefix,
                            const std::vector<std::string> & inFieldNames,
                            const std::vector<std::string> & inDOFs,
                            const std::vector<std::string> & outDOFs,
                            const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                            const std::string & elementBlock,
                            const GlobalIndexer & indexerSrc,
                            const GlobalIndexer & indexerDest);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
private:
  typedef typename TRAITS::Jacobian::ScalarT ScalarT;

  void buildSrcToDestMap(const std::string & elementBlock,
                         const GlobalIndexer & indexerSrc,
                         const GlobalIndexer & indexerDest);

  // Build a source to destination map using all the pairs
  // of field numers in the <code>fieldNumberMaps</code>
  void buildSrcToDestMap(const std::string & elementBlock,
                         const std::map<int,int> & fieldNumberMaps,
                         const GlobalIndexer & indexerSrc,
                         const GlobalIndexer & indexerDest);

  // fields to be modified
  std::vector< PHX::MDField<const ScalarT> > inFields_;

  // fields that need to be modified
  std::vector< PHX::MDField<ScalarT> > outFields_;

  // This allows indexing into a destination sized vector and 
  // maps to a source vector. If a value is less then 0
  // then that implies that value is not mapped. That is a strange
  // case but this structure supports it
  Kokkos::View<int*> dstFromSrcMapView_;

  ReorderADValues_Evaluator() {}
  ReorderADValues_Evaluator(const ReorderADValues_Evaluator &) {}
};

}

// **************************************************************
#endif
