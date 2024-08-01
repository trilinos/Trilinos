// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_SCATTER_DIRICHLET_RESIDUAL_TPETRA_DECL_HPP
#define PANZER_EVALUATOR_SCATTER_DIRICHLET_RESIDUAL_TPETRA_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"

#include "Panzer_NodeType.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

class GlobalIndexer;

/** \brief Pushes residual values into the residual vector for a 
           Newton-based solve

    Currently makes an assumption that the stride is constant for dofs
    and that the number of dofs is equal to the size of the solution
    names vector.

*/
template<typename EvalT, typename Traits,typename LO,typename GO,typename NodeT=panzer::TpetraNodeType>
class ScatterDirichletResidual_Tpetra;

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual 
// **************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidual_Tpetra<panzer::Traits::Residual,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {
  
public:
  ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer)
     : globalIndexer_(indexer) {}
  
  ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
                                  const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_Tpetra<panzer::Traits::Residual,TRAITS,LO,GO>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;
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

// **************************************************************
// Tangent 
// **************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Tangent, TRAITS>,
    public panzer::CloneableEvaluator  {
  
public:
  ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer)
     : globalIndexer_(indexer) {}
  
  ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
                                  const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_Tpetra<panzer::Traits::Tangent,TRAITS,LO,GO>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Tangent::ScalarT ScalarT;
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
  std::vector< Teuchos::ArrayRCP<double> > dfdp_vectors_;

  //! If set to true, allows runtime disabling of dirichlet BCs on node-by-node basis
  bool checkApplyBC_;

  // If set to true, scattering an initial condition
  bool scatterIC_;

  // Allows runtime disabling of dirichlet BCs on node-by-node basis
  std::vector< PHX::MDField<const bool,Cell,NODE> > applyBC_;
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>,
    public panzer::CloneableEvaluator  {
  
public:
  ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer)
     : globalIndexer_(indexer) {}
  
  ScatterDirichletResidual_Tpetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
                                  const Teuchos::ParameterList& p);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);
  
  void evaluateFields(typename TRAITS::EvalData workset);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_Tpetra<panzer::Traits::Jacobian,TRAITS,LO,GO>(globalIndexer_,pl)); }
  
private:

  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;
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
  std::size_t num_eq;

  std::size_t side_subcell_dim_;
  std::size_t local_side_id_;

  Teuchos::RCP<typename LOC::VectorType> dirichletCounter_;

  ScatterDirichletResidual_Tpetra();

  std::string globalDataKey_; // what global data does this fill?
  Teuchos::RCP<const TpetraLinearObjContainer<double,LO,GO,NodeT> > tpetraContainer_;

  //! If set to true, allows runtime disabling of dirichlet BCs on node-by-node basis
  bool checkApplyBC_;

  // Allows runtime disabling of dirichlet BCs on node-by-node basis
  std::vector< PHX::MDField<const bool,Cell,NODE> > applyBC_;
};

}

// optionally include hessian support
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_ScatterDirichletResidual_Tpetra_Hessian.hpp"
#endif

// **************************************************************
#endif
