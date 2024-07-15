// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_GatherSolution_BlockedTpetra_Hessian_hpp__
#define __Panzer_GatherSolution_BlockedTpetra_Hessian_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// blocked Tpetra gather solution file

namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
class GatherSolution_BlockedTpetra<panzer::Traits::Hessian,TRAITS,S,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Hessian, TRAITS>,
    public panzer::CloneableEvaluator  {


public:

   GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer)
     : gidIndexer_(indexer) {}

   GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & /* indexer */,
                                const Teuchos::ParameterList& /* p */) {}

  void postRegistrationSetup(typename TRAITS::SetupData /* d */,
                             PHX::FieldManager<TRAITS>& /* vm */) {}

  void preEvaluate(typename TRAITS::PreEvalData /* d */) {}

  void evaluateFields(typename TRAITS::EvalData /* d */) {}

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_BlockedTpetra<panzer::Traits::Hessian,TRAITS,S,LO,GO>(gidIndexer_,pl)); }


private:
  typedef typename panzer::Traits::Hessian::ScalarT ScalarT;

  typedef BlockedTpetraLinearObjContainer<S,LO,GO,NodeT> ContainerType;
  typedef Tpetra::Vector<S,LO,GO,NodeT> VectorType;
  typedef Tpetra::CrsMatrix<S,LO,GO,NodeT> CrsMatrixType;
  typedef Tpetra::CrsGraph<LO,GO,NodeT> CrsGraphType;
  typedef Tpetra::Map<LO,GO,NodeT> MapType;
  typedef Tpetra::Import<LO,GO,NodeT> ImportType;
  typedef Tpetra::Export<LO,GO,NodeT> ExportType;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const BlockedDOFManager> gidIndexer_;

  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  Teuchos::RCP<std::vector<std::string> > indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const BlockedTpetraLinearObjContainer<S,LO,GO,NodeT> > blockedContainer_;

  // Fields for storing tangent components dx/dp of solution vector x
  // These are not actually used by the residual specialization of this evaluator,
  // even if they are supplied, but it is useful to declare them as dependencies anyway
  // when saving the tangent components to the output file
  bool has_tangent_fields_;
  std::vector< std::vector< PHX::MDField<ScalarT,Cell,NODE> > > tangentFields_;

  GatherSolution_BlockedTpetra();
};

}

// **************************************************************
#endif

#endif
