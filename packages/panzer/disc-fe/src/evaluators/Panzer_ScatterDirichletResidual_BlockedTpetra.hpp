// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_SCATTER_DIRICHLET_RESIDUAL_BLOCKEDTPETRA_HPP
#define PANZER_EVALUATOR_SCATTER_DIRICHLET_RESIDUAL_BLOCKEDTPETRA_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace Thyra {
   template <typename ScalarT> class ProductVectorBase;
}

namespace panzer {

class BlockedDOFManager;
class GlobalIndexer;

/** \brief Pushes residual values into the residual vector for a
           Newton-based solve

    Currently makes an assumption that the stride is constant for dofs
    and that the number of dofs is equal to the size of the solution
    names vector.

*/
template <typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT=panzer::TpetraNodeType>
class ScatterDirichletResidual_BlockedTpetra
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {
public:
   typedef typename EvalT::ScalarT ScalarT;
   ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & /* indexer */)
   { }

   ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & gidProviders,
                                          const Teuchos::ParameterList& p);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_BlockedTpetra<EvalT,TRAITS,LO,GO,NodeT>(Teuchos::null,pl)); }

  void postRegistrationSetup(typename TRAITS::SetupData /* d */, PHX::FieldManager<TRAITS>& /* vm */)
   { }
  void evaluateFields(typename TRAITS::EvalData /* d */)
   { std::cout << "unspecialized version of \"ScatterDirichletResidual_BlockedTpetra::evaluateFields\" on \""+PHX::print<EvalT>()+"\" should not be used!" << std::endl;
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
template <typename TRAITS,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {

public:
  ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer)
     : globalIndexer_(indexer) {}

  ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
                                  const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData workset);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual,TRAITS,LO,GO,NodeT>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;
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

  //! Vector of global indexers, one for each scattered field
  //! respectively. This is the global indexer for the Thyra
  //! ProductVector sub-block.
  std::vector<Teuchos::RCP<const panzer::GlobalIndexer>> fieldGlobalIndexers_;

  //! Field IDs in the local product vector block (not global field id)
  std::vector<int> fieldIds_;

  //! Returns the index into the Thyra ProductVector sub-block. Size
  //! of number of fields to scatter.
  std::vector<int> productVectorBlockIndex_;

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  //! Local indices for unknowns
  PHX::View<LO**> worksetLIDs_;

  //! Offset into the cell lids for each field
  std::vector<PHX::View<int*>> fieldOffsets_;

  //! The local basis index corresponding to the fieldOffset_. Used to
  //! index into the basis index of MDFields. This is only required
  //! for tangent/normal BCs.
  std::vector<PHX::View<int*>> basisIndexForMDFieldOffsets_;

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

// **************************************************************
// Jacobian
// **************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
class ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>,
    public panzer::CloneableEvaluator  {

public:
  ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer)
     : globalIndexer_(indexer) {}

  ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
                                  const Teuchos::ParameterList& p);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void evaluateFields(typename TRAITS::EvalData workset);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;
  typedef typename TRAITS::RealType RealType;

  typedef BlockedTpetraLinearObjContainer<RealType,LO,GO,NodeT> ContainerType;
  typedef Tpetra::Operator<RealType,LO,GO,NodeT> OperatorType;
  typedef Tpetra::CrsMatrix<RealType,LO,GO,NodeT> CrsMatrixType;
  typedef Tpetra::Map<LO,GO,NodeT> MapType;

  typedef Thyra::TpetraLinearOp<RealType,LO,GO,NodeT> ThyraLinearOp;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::BlockedDOFManager> globalIndexer_;

  std::vector<int> fieldIds_; // field IDs needing mapping

  //! Returns the index into the Thyra ProductVector sub-block. Size
  //! of number of fields to scatter
  std::vector<int> productVectorBlockIndex_;

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::size_t side_subcell_dim_;
  std::size_t local_side_id_;

  Teuchos::RCP<Thyra::ProductVectorBase<double> > dirichletCounter_;
  std::string globalDataKey_; // what global data does this fill?
  Teuchos::RCP<const BlockedTpetraLinearObjContainer<RealType,LO,GO,NodeT> > blockedContainer_;

  //! Local indices for unknowns
  PHX::View<LO**> worksetLIDs_;

  //! Offset into the cell lids for each field. Size of number of fields to scatter.
  std::vector<PHX::View<int*>> fieldOffsets_;

  //! The local basis index corresponding to the fieldOffset_. Used to
  //! index into the basis index of MDFields. This is only required
  //! for tangent/normal BCs.
  std::vector<PHX::View<int*>> basisIndexForMDFieldOffsets_;

  //! The offset values of the blocked DOFs per element. Size of number of blocks in the product vector + 1. The plus one is a sentinel.
  PHX::View<LO*> blockOffsets_;

  //! If set to true, allows runtime disabling of dirichlet BCs on node-by-node basis
  bool checkApplyBC_;

  // Allows runtime disabling of dirichlet BCs on node-by-node basis
  std::vector< PHX::MDField<const bool,Cell,NODE> > applyBC_;

  //! Used to allocate temporary space on device
  static constexpr int maxDerivativeArraySize_ = 256;

  ScatterDirichletResidual_BlockedTpetra();
};

}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_ScatterDirichletResidual_BlockedTpetra_Hessian.hpp"
#endif

// **************************************************************
#endif
