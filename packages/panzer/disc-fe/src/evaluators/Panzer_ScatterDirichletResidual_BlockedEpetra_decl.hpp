// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_SCATTER_DIRICHLET_RESIDUAL_BLOCKEDEPETRA_DECL_HPP
#define PANZER_EVALUATOR_SCATTER_DIRICHLET_RESIDUAL_BLOCKEDEPETRA_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace Thyra {
   template <typename ScalarT> class ProductVectorBase;
   template <typename ScalarT> class BlockedLinearOpBase;
}

namespace panzer {

class GlobalIndexer;

/** \brief Pushes residual values into the residual vector for a
           Newton-based solve

    Currently makes an assumption that the stride is constant for dofs
    and that the number of dofs is equal to the size of the solution
    names vector.

*/
template<typename EvalT, typename TRAITS,typename LO,typename GO> class ScatterDirichletResidual_BlockedEpetra
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {
public:
  typedef typename EvalT::ScalarT ScalarT;

  ScatterDirichletResidual_BlockedEpetra(const Teuchos::ParameterList& p) {}

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_BlockedEpetra<EvalT,TRAITS,LO,GO>(pl)); }

  void postRegistrationSetup(typename TRAITS::SetupData d, PHX::FieldManager<TRAITS>& vm)
   { }
  void evaluateFields(typename TRAITS::EvalData d)
   { std::cout << "unspecialized version of \"ScatterDirichletResidual_BlockedEpetra::evaluateFields\" on \""+PHX::print<EvalT>()+"\" should not be used!" << std::endl;
    TEUCHOS_ASSERT(false); }
};

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {

public:

  ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                                         const std::vector<Teuchos::RCP<const GlobalIndexer> > & /* cIndexers */)
     : rowIndexers_(rIndexers) {}

  ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                                         const std::vector<Teuchos::RCP<const GlobalIndexer> > & cIndexers,
                                         const Teuchos::ParameterList& p,
                                         bool useDiscreteAdjoint=false);

  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData workset);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual,TRAITS,LO,GO>(rowIndexers_,colIndexers_,pl)); }

private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;
  using BCFieldType = PHX::MDField<const bool,Cell,NODE>;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  std::vector<Teuchos::RCP<const GlobalIndexer> > rowIndexers_;
  std::vector<Teuchos::RCP<const GlobalIndexer> > colIndexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::size_t num_nodes;

  std::size_t side_subcell_dim_;
  std::size_t local_side_id_;

  Teuchos::RCP<Thyra::ProductVectorBase<double> > dirichletCounter_;
  Teuchos::RCP<Thyra::ProductVectorBase<double> > r_;
  std::string globalDataKey_; // what global data does this fill?

  //! If set to true, allows runtime disabling of dirichlet BCs on node-by-node basis
  bool checkApplyBC_;

  // If set to true, scattering an initial condition
  bool scatterIC_;

  // Allows runtime disabling of dirichlet BCs on node-by-node basis
  std::vector< BCFieldType > applyBC_;

  ScatterDirichletResidual_BlockedEpetra() {}
};

// **************************************************************
// Tangent
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Tangent,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Tangent, TRAITS>,
    public panzer::CloneableEvaluator  {

public:

  ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                                         const std::vector<Teuchos::RCP<const GlobalIndexer> > & /* cIndexers */)
     : rowIndexers_(rIndexers) {}

  ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                                         const std::vector<Teuchos::RCP<const GlobalIndexer> > & cIndexers,
                                         const Teuchos::ParameterList& p,
                                         bool useDiscreteAdjoint=false);

  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData workset);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Tangent,TRAITS,LO,GO>(rowIndexers_,colIndexers_,pl)); }

private:
  typedef typename panzer::Traits::Tangent::ScalarT ScalarT;
  using BCFieldType = PHX::MDField<const bool,Cell,NODE>;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  std::vector<Teuchos::RCP<const GlobalIndexer> > rowIndexers_;
  std::vector<Teuchos::RCP<const GlobalIndexer> > colIndexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::size_t num_nodes;

  std::size_t side_subcell_dim_;
  std::size_t local_side_id_;

  Teuchos::RCP<Thyra::ProductVectorBase<double> > dirichletCounter_;
  Teuchos::RCP<Thyra::ProductVectorBase<double> > r_;
  std::string globalDataKey_; // what global data does this fill?

  //! If set to true, allows runtime disabling of dirichlet BCs on node-by-node basis
  bool checkApplyBC_;

  // If set to true, scattering an initial condition
  bool scatterIC_;

  // Allows runtime disabling of dirichlet BCs on node-by-node basis
  std::vector< BCFieldType > applyBC_;

  ScatterDirichletResidual_BlockedEpetra() {}
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>,
    public panzer::CloneableEvaluator  {

public:

  ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                                         const std::vector<Teuchos::RCP<const GlobalIndexer> > & cIndexers)
     : rowIndexers_(rIndexers), colIndexers_(cIndexers) {}

  ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                                         const std::vector<Teuchos::RCP<const GlobalIndexer> > & cIndexers,
                                         const Teuchos::ParameterList& p,
                                         bool useDiscreteAdjoint=false);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void evaluateFields(typename TRAITS::EvalData workset);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian,TRAITS,LO,GO>(rowIndexers_,colIndexers_,pl)); }

private:

  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;
  using BCFieldType = PHX::MDField<const bool,Cell,NODE>;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  std::vector<Teuchos::RCP<const GlobalIndexer> > rowIndexers_;
  std::vector<Teuchos::RCP<const GlobalIndexer> > colIndexers_;

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

  Teuchos::RCP<Thyra::ProductVectorBase<double> > r_;
  Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > Jac_;

  //! If set to true, allows runtime disabling of dirichlet BCs on node-by-node basis
  bool checkApplyBC_;

  // Allows runtime disabling of dirichlet BCs on node-by-node basis
  std::vector< BCFieldType > applyBC_;

  ScatterDirichletResidual_BlockedEpetra();
};

}

// optionally include hessian support
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_ScatterDirichletResidual_BlockedEpetra_Hessian.hpp"
#endif

// **************************************************************
#endif
