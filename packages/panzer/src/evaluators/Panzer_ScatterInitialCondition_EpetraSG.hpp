#ifdef HAVE_STOKHOS
#ifndef PANZER_SCATTER_INITIAL_CONDITION_EPETRA_SG_HPP
#define PANZER_SCATTER_INITIAL_CONDITION_EPETRA_SG_HPP

namespace panzer {

// **************************************************************
// SGResidual 
// **************************************************************
template<typename Traits,typename LO,typename GO>
class ScatterInitialCondition_Epetra<panzer::Traits::SGResidual,Traits,LO,GO>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::SGResidual, Traits>,
    public panzer::CloneableEvaluator {
  
public:
  ScatterInitialCondition_Epetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer) 
     : globalIndexer_(indexer) {}
  
  ScatterInitialCondition_Epetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
                         const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterInitialCondition_Epetra<panzer::Traits::SGResidual,Traits,LO,GO>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::SGResidual::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

};

// **************************************************************
// SGJacobian
// **************************************************************
template<typename Traits,typename LO,typename GO>
class ScatterInitialCondition_Epetra<panzer::Traits::SGJacobian,Traits,LO,GO>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::SGJacobian, Traits>, 
    public panzer::CloneableEvaluator {
  
public:
  
  ScatterInitialCondition_Epetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer) 
     : globalIndexer_(indexer) {}

  ScatterInitialCondition_Epetra(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
                         const Teuchos::ParameterList& pl);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterInitialCondition_Epetra<panzer::Traits::SGJacobian,Traits,LO,GO>(globalIndexer_,pl)); }

private:

  typedef typename panzer::Traits::SGJacobian::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  ScatterInitialCondition_Epetra();
};

}

// **************************************************************

#endif
#endif
