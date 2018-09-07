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

class UniqueGlobalIndexerBase;

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
                            const UniqueGlobalIndexerBase & indexerSrc,
                            const UniqueGlobalIndexerBase & indexerDest);

  ReorderADValues_Evaluator(const std::string & outPrefix,
                            const std::vector<std::string> & inFieldNames,
                            const std::vector<std::string> & inDOFs,
                            const std::vector<std::string> & outDOFs,
                            const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                            const std::string & elementBlock,
                            const UniqueGlobalIndexerBase & indexerSrc,
                            const UniqueGlobalIndexerBase & indexerDest);

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
                            const UniqueGlobalIndexerBase & indexerSrc,
                            const UniqueGlobalIndexerBase & indexerDest);

  ReorderADValues_Evaluator(const std::string & outPrefix,
                            const std::vector<std::string> & inFieldNames,
                            const std::vector<std::string> & inDOFs,
                            const std::vector<std::string> & outDOFs,
                            const std::vector<Teuchos::RCP<PHX::DataLayout> > & fieldLayouts,
                            const std::string & elementBlock,
                            const UniqueGlobalIndexerBase & indexerSrc,
                            const UniqueGlobalIndexerBase & indexerDest);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
private:
  typedef typename TRAITS::Jacobian::ScalarT ScalarT;

  void buildSrcToDestMap(const std::string & elementBlock,
                         const UniqueGlobalIndexerBase & indexerSrc,
                         const UniqueGlobalIndexerBase & indexerDest);

  // Build a source to destination map using all the pairs
  // of field numers in the <code>fieldNumberMaps</code>
  void buildSrcToDestMap(const std::string & elementBlock,
                         const std::map<int,int> & fieldNumberMaps,
                         const UniqueGlobalIndexerBase & indexerSrc,
                         const UniqueGlobalIndexerBase & indexerDest);

  // fields to be modified
  std::vector< PHX::MDField<const ScalarT> > inFields_;

  // fields that need to be modified
  std::vector< PHX::MDField<ScalarT> > outFields_;

  // This allows indexing into a destination sized vector and 
  // maps to a source vector. If a value is less then 0
  // then that implies that value is not mapped. That is a strange
  // case but this structure supports it
  std::vector<int> dstFromSrcMap_;

  ReorderADValues_Evaluator() {}
  ReorderADValues_Evaluator(const ReorderADValues_Evaluator &) {}
};

}

// **************************************************************
#endif
