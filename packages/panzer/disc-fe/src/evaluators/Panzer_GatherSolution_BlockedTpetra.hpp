// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_GATHER_SOLUTION_BLOCKEDTPETRA_DECL_HPP
#define PANZER_EVALUATOR_GATHER_SOLUTION_BLOCKEDTPETRA_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
class BlockedTpetraLinearObjContainer;

class GlobalIndexer;
class BlockedDOFManager;

/** \brief Gathers solution values from the Newton solution vector into
    the nodal fields of the field manager

    Currently makes an assumption that the stride is constant for dofs
    and that the nmber of dofs is equal to the size of the solution
    names vector.
*/
template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT=panzer::TpetraNodeType>
class GatherSolution_BlockedTpetra
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {
public:
   typedef typename EvalT::ScalarT ScalarT;

   GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer)
   { }

   GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & gidProviders,
                                const Teuchos::ParameterList& p);

   virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp(new GatherSolution_BlockedTpetra<EvalT,TRAITS,S,LO,GO>(Teuchos::null,pl)); }

   void postRegistrationSetup(typename TRAITS::SetupData d, PHX::FieldManager<TRAITS>& vm)
   { }

   void evaluateFields(typename TRAITS::EvalData d)
   { std::cout << "unspecialized version of \"GatherSolution_BlockedTpetra::evaluateFields\" on \""+PHX::print<EvalT>()+"\" should not be used!" << std::endl;
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
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
class GatherSolution_BlockedTpetra<panzer::Traits::Residual,TRAITS,S,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {


public:

   GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer)
     : globalIndexer_(indexer) {}

   GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
                                const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_BlockedTpetra<panzer::Traits::Residual,TRAITS,S,LO,GO>(globalIndexer_,pl)); }


private:
  typedef typename panzer::Traits::Residual EvalT;
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;

  typedef BlockedTpetraLinearObjContainer<S,LO,GO,NodeT> ContainerType;
  typedef Tpetra::Vector<S,LO,GO,NodeT> VectorType;
  typedef Tpetra::CrsMatrix<S,LO,GO,NodeT> CrsMatrixType;
  typedef Tpetra::CrsGraph<LO,GO,NodeT> CrsGraphType;
  typedef Tpetra::Map<LO,GO,NodeT> MapType;
  typedef Tpetra::Import<LO,GO,NodeT> ImportType;
  typedef Tpetra::Export<LO,GO,NodeT> ExportType;

  //! Maps the local (field,element,basis) triplet to a global ID for
  // scattering
  Teuchos::RCP<const BlockedDOFManager> globalIndexer_;

  //! Field IDs in the local product vector block (not global field id)
  std::vector<int> fieldIds_;

  //! Vector of global indexers, one for each field to gather, respectively
  std::vector<Teuchos::RCP<const panzer::GlobalIndexer>> fieldGlobalIndexers_;

  //! Returns the index to the Thyra ProductVector sub-block. Size
  //! of number of fields to gather
  std::vector<int> productVectorBlockIndex_;

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  std::vector<std::string> indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const BlockedTpetraLinearObjContainer<S,LO,GO,NodeT> > blockedContainer_;

  // Fields for storing tangent components dx/dp of solution vector x
  // These are not actually used by the residual specialization of this evaluator,
  // even if they are supplied, but it is useful to declare them as dependencies anyway
  // when saving the tangent components to the output file
  bool has_tangent_fields_;
  std::vector< std::vector< PHX::MDField<const ScalarT,Cell,NODE> > > tangentFields_;

  //! Local indices for unknowns
  PHX::View<LO**> worksetLIDs_;

  //! Offset into the cell lids for each field
  std::vector<PHX::View<int*>> fieldOffsets_;

  GatherSolution_BlockedTpetra();
};

// **************************************************************
// Tangent
// **************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
class GatherSolution_BlockedTpetra<panzer::Traits::Tangent,TRAITS,S,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Tangent, TRAITS>,
    public panzer::CloneableEvaluator  {


public:

   GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer)
     : gidIndexer_(indexer) {}

   GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
                                const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_BlockedTpetra<panzer::Traits::Tangent,TRAITS,S,LO,GO>(gidIndexer_,pl)); }


private:
  typedef typename panzer::Traits::Tangent EvalT;
  typedef typename panzer::Traits::Tangent::ScalarT ScalarT;
  //typedef typename panzer::Traits::RealType RealT;

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

  std::vector<std::string> indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const BlockedTpetraLinearObjContainer<S,LO,GO,NodeT> > blockedContainer_;

  // Fields for storing tangent components dx/dp of solution vector x
  bool has_tangent_fields_;
  std::vector< std::vector< PHX::MDField<const ScalarT,Cell,NODE> > > tangentFields_;

  GatherSolution_BlockedTpetra();
};

// **************************************************************
// Jacobian
// **************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
class GatherSolution_BlockedTpetra<panzer::Traits::Jacobian,TRAITS,S,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>,
    public panzer::CloneableEvaluator  {

public:
  GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer)
     : globalIndexer_(indexer) {}

  GatherSolution_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
                               const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_BlockedTpetra<panzer::Traits::Jacobian,TRAITS,S,LO,GO>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Jacobian EvalT;
  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;
  typedef typename TRAITS::RealType RealType;

  typedef BlockedTpetraLinearObjContainer<S,LO,GO,NodeT> ContainerType;
  typedef Tpetra::Vector<S,LO,GO,NodeT> VectorType;
  typedef Tpetra::CrsMatrix<S,LO,GO,NodeT> CrsMatrixType;
  typedef Tpetra::CrsGraph<LO,GO,NodeT> CrsGraphType;
  typedef Tpetra::Map<LO,GO,NodeT> MapType;
  typedef Tpetra::Import<LO,GO,NodeT> ImportType;
  typedef Tpetra::Export<LO,GO,NodeT> ExportType;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const BlockedDOFManager> globalIndexer_;

  std::vector<int> fieldIds_; // field IDs needing mapping

  //! Returns the index into the Thyra ProductVector sub-block. Size
  //! of number of fields to scatter
  std::vector<int> productVectorBlockIndex_;

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  std::vector<std::string> indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  bool disableSensitivities_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const BlockedTpetraLinearObjContainer<S,LO,GO,NodeT> > blockedContainer_;

  //! Local indices for unknowns
  PHX::View<LO**> worksetLIDs_;

  //! Offset into the cell lids for each field. Size of number of fields to scatter.
  std::vector<PHX::View<int*>> fieldOffsets_;

  //! The offset values of the blocked DOFs per element. Size of number of blocks in the product vector + 1. The plus one is a sentinel.
  PHX::View<LO*> blockOffsets_;

  GatherSolution_BlockedTpetra();
};

}

#include "Panzer_GatherSolution_BlockedTpetra_impl.hpp"

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_GatherSolution_BlockedTpetra_Hessian.hpp"
#endif

// **************************************************************
#endif
