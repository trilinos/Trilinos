// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ScatterDirichletResidual_BlockedEpetra_Hessian_hpp__
#define __Panzer_ScatterDirichletResidual_BlockedEpetra_Hessian_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// blocked Epetra scatter dirichlet residual file


namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Hessian, TRAITS>,
    public panzer::CloneableEvaluator  {
  
public:
  ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > & rIndexers,
                                         const std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > & cIndexers)
     : rowIndexers_(rIndexers), colIndexers_(cIndexers) {}
  
  ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > & rIndexers,
                                         const std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > & cIndexers,
                                         const Teuchos::ParameterList& p,
                                         bool useDiscreteAdjoint=false);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>(rowIndexers_,colIndexers_,pl)); }

private:
  typedef typename panzer::Traits::Hessian::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > rowIndexers_;
  std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > colIndexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::size_t num_nodes;
  std::size_t num_eq;

  std::size_t side_subcell_dim_;
  std::size_t local_side_id_;

  Teuchos::RCP<Thyra::ProductVectorBase<double> > dirichletCounter_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > Jac_;

  //! If set to true, allows runtime disabling of dirichlet BCs on node-by-node basis
  bool checkApplyBC_;

  // Allows runtime disabling of dirichlet BCs on node-by-node basis
  std::vector< PHX::MDField<const bool,Cell,NODE> > applyBC_;

  ScatterDirichletResidual_BlockedEpetra();
};

}

// **************************************************************
#endif

#endif
