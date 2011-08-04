#ifdef HAVE_STOKHOS

#ifndef PANZER_EVALUATOR_GATHER_SOLUTION_EPETRA_SG_HPP
#define PANZER_EVALUATOR_GATHER_SOLUTION_EPETRA_SG_HPP


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

#include "Panzer_GatherSolution_EpetraSGT.hpp"

// **************************************************************
#endif
#endif
