// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ScatterDirichletResidual_BlockedTpetra_Hessian_hpp__
#define __Panzer_ScatterDirichletResidual_BlockedTpetra_Hessian_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// blocked Tpetra scatter dirichlet file

namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Hessian, TRAITS>,
    public panzer::CloneableEvaluator  {
  
public:
  ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer)
     : globalIndexer_(indexer) { }
  
  ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
                                  const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
 
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Hessian::ScalarT ScalarT;
  typedef typename TRAITS::RealType RealType;

  typedef BlockedTpetraLinearObjContainer<RealType,LO,GO,NodeT> ContainerType;
  typedef Tpetra::Vector<RealType,LO,GO,NodeT> VectorType;
  typedef Tpetra::CrsMatrix<RealType,LO,GO,NodeT> CrsMatrixType;
  typedef Tpetra::CrsGraph<LO,GO,NodeT> CrsGraphType;
  typedef Tpetra::Map<LO,GO,NodeT> MapType;
  typedef Tpetra::Import<LO,GO,NodeT> ImportType;
  typedef Tpetra::Export<LO,GO,NodeT> ExportType;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::BlockedDOFManager> globalIndexer_;

  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::size_t num_nodes;

  std::size_t side_subcell_dim_;
  std::size_t local_side_id_;

  Teuchos::RCP<Thyra::ProductVectorBase<double> > dirichletCounter_;
  std::string globalDataKey_; // what global data does this fill?
  Teuchos::RCP<const BlockedTpetraLinearObjContainer<RealType,LO,GO,NodeT> > blockedContainer_;

  //! If set to true, allows runtime disabling of dirichlet BCs on node-by-node basis
  bool checkApplyBC_;

  // If set to true, scattering an initial condition
  bool scatterIC_;

  // Allows runtime disabling of dirichlet BCs on node-by-node basis
  std::vector< PHX::MDField<const bool,Cell,NODE> > applyBC_;

  ScatterDirichletResidual_BlockedTpetra() {}
};

}

// **************************************************************
#endif

#endif
