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

#ifdef HAVE_STOKHOS

#ifndef PANZER_EVALUATOR_GATHER_SOLUTION_EPETRA_SG_DECL_HPP
#define PANZER_EVALUATOR_GATHER_SOLUTION_EPETRA_SG_DECL_HPP


//
// Note: This file is included in Panzer_GatherSolution_Epetra.hpp
//       so many of the required includes and data types are defined
//       there
//

namespace panzer {


// **************************************************************
// SGResidual 
// **************************************************************
template<typename Traits,typename LO,typename GO>
class GatherSolution_Epetra<panzer::Traits::SGResidual,Traits,LO,GO>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::SGResidual, Traits>,
    public panzer::CloneableEvaluator  {
   
  
public:
  
  GatherSolution_Epetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer) :
     globalIndexer_(indexer) {}

  GatherSolution_Epetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
                        const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_Epetra<panzer::Traits::SGResidual,Traits,LO,GO>(globalIndexer_,pl)); }
  
private:

  typedef typename panzer::Traits::SGResidual::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  Teuchos::RCP<std::vector<std::string> > indexerNames_;
  bool useTimeDerivativeSolutionVector_;

  GatherSolution_Epetra();
};

// **************************************************************
// SGJacobian
// **************************************************************
template<typename Traits,typename LO,typename GO>
class GatherSolution_Epetra<panzer::Traits::SGJacobian,Traits,LO,GO>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::SGJacobian, Traits>,
    public panzer::CloneableEvaluator  {
  
public:
  GatherSolution_Epetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer) :
     globalIndexer_(indexer) {}
  
  GatherSolution_Epetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
                        const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_Epetra<panzer::Traits::SGJacobian,Traits,LO,GO>(globalIndexer_,pl)); }
  
private:

  typedef typename panzer::Traits::SGJacobian::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  Teuchos::RCP<std::vector<std::string> > indexerNames_;
  bool useTimeDerivativeSolutionVector_;

  GatherSolution_Epetra();
};

}

// **************************************************************
#endif
#endif
