// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ScatterResidual_Epetra_Hessian_hpp__
#define __Panzer_ScatterResidual_Epetra_Hessian_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// Epetra scatter residual file

namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class ScatterResidual_Epetra<panzer::Traits::Hessian,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Hessian, TRAITS>,
    public panzer::CloneableEvaluator {
  
public:
  ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                         const Teuchos::RCP<const panzer::GlobalIndexer> & cIndexer=Teuchos::null,
                         bool useDiscreteAdjoint=false)
     : globalIndexer_(indexer), colGlobalIndexer_(cIndexer), useDiscreteAdjoint_(useDiscreteAdjoint)  {}
  
  ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                         const Teuchos::RCP<const panzer::GlobalIndexer> & cIndexer,
                         const Teuchos::ParameterList& p,bool=false);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm) ;

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_Epetra<panzer::Traits::Hessian,TRAITS,LO,GO>(globalIndexer_,colGlobalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Hessian::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_, colGlobalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const EpetraLinearObjContainer> epetraContainer_;

  ScatterResidual_Epetra();

  bool useDiscreteAdjoint_;
};

}

// **************************************************************
#endif

#endif
