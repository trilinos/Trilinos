// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ScatterDirichletResidual_Tpetra_Hessian_hpp__
#define __Panzer_ScatterDirichletResidual_Tpetra_Hessian_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// Tpetra scatter dirichlet residual file

namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************

/**
 * \brief Hessian specialization of ScatterDirichletResidual_Tpetra.
 *
 * Enforces Dirichlet boundary values at the affected degrees of
 * freedom (overwriting, rather than accumulating into, the
 * corresponding residual/matrix entries) while additionally
 * propagating second-derivative (Hessian) information, so that
 * Dirichlet BCs are correctly reflected in Hessian-vector products
 * used by e.g. optimization or uncertainty quantification drivers.
 */
template<typename TRAITS,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidual_Tpetra<panzer::Traits::Hessian,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Hessian, TRAITS>,
    public panzer::CloneableEvaluator  {

public:
  /// \brief Construct with a global indexer only; scattered field names must be set later (e.g. via clone()).
  ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer)
     : globalIndexer_(indexer) {}

  /// \brief Construct from a global indexer and a ParameterList describing the fields to scatter and the Dirichlet BC side set.
  ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
                                  const Teuchos::ParameterList& p);

  /// \brief Looks up and caches the Tpetra linear object container and field offsets needed by evaluateFields(). Called once before the first evaluation.
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  /// \brief Fetches the Tpetra containers (residual/matrix/Dirichlet counter) to scatter into for the upcoming fill.
  void preEvaluate(typename TRAITS::PreEvalData d);

  /// \brief Overwrites the Dirichlet-affected residual, matrix, and Hessian-vector-product entries for this workset with the scattered field values and their derivatives.
  void evaluateFields(typename TRAITS::EvalData workset);

  /// \brief Creates a copy of this evaluator configured from a new ParameterList, sharing the same global indexer.
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_Tpetra<panzer::Traits::Hessian,TRAITS,LO,GO>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Hessian::ScalarT ScalarT;
  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::size_t num_nodes;

  std::size_t side_subcell_dim_;
  std::size_t local_side_id_;

  ScatterDirichletResidual_Tpetra() {}

  Teuchos::RCP<typename LOC::VectorType> dirichletCounter_;

  std::string globalDataKey_; // what global data does this fill?
  Teuchos::RCP<const LOC> tpetraContainer_;

  //! If set to true, allows runtime disabling of dirichlet BCs on node-by-node basis
  bool checkApplyBC_;

  // If set to true, scattering an initial condition
  bool scatterIC_;

  // Allows runtime disabling of dirichlet BCs on node-by-node basis
  std::vector< PHX::MDField<const bool,Cell,NODE> > applyBC_;
};

}

// **************************************************************
#endif // end hessian support

#endif
