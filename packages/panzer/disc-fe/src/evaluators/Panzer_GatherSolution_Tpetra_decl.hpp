// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_GATHER_SOLUTION_TPETRA_DECL_HPP
#define PANZER_EVALUATOR_GATHER_SOLUTION_TPETRA_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_KokkosViewOfViews.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"

#include"Panzer_NodeType.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

class GlobalIndexer; //forward declaration

/** \brief Gathers solution values from the Newton solution vector into
    the nodal fields of the field manager

    Currently makes an assumption that the stride is constant for dofs
    and that the nmber of dofs is equal to the size of the solution
    names vector.
*/
template<typename EvalT, typename Traits,typename LO,typename GO,typename NodeT=panzer::TpetraNodeType>
class GatherSolution_Tpetra;

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual
// **************************************************************

/** \brief Residual specialization of GatherSolution_Tpetra.
 *
 * Reads cell-local values directly out of the Tpetra solution vector
 * (no derivative information), using #globalIndexer_ to map each
 * (field, element, basis) triplet to its global unknown.
 */
template<typename TRAITS,typename LO,typename GO,typename NodeT>
class GatherSolution_Tpetra<panzer::Traits::Residual,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {


public:

  /// \brief Construct with a global indexer only; solution/parameter names must be set later (e.g. via clone()).
  GatherSolution_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer) :
     globalIndexer_(indexer) {}

  /// \brief Construct from a global indexer and a ParameterList giving the DOF names to gather (see the class-level Panzer_GatherSolution_Input options).
  GatherSolution_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                        const Teuchos::ParameterList& p);

  /// \brief Looks up and caches the Tpetra linear object container and field offsets needed by evaluateFields(). Called once before the first evaluation.
  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  /// \brief Fetches the Tpetra solution vector to gather from for the upcoming fill.
  void preEvaluate(typename TRAITS::PreEvalData d);

  /// \brief Copies solution values for this workset from the Tpetra vector into #gatherFields_.
  void evaluateFields(typename TRAITS::EvalData d);

  /// \brief Creates a copy of this evaluator configured from a new ParameterList, sharing the same global indexer.
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_Tpetra<panzer::Traits::Residual,TRAITS,LO,GO,NodeT>(globalIndexer_,pl)); }

  // for testing purposes
  /// \brief Returns the FieldTag of the i-th gathered field. For testing purposes.
  const PHX::FieldTag & getFieldTag(int i) const
  { TEUCHOS_ASSERT(i < Teuchos::as<int>(gatherFields_.size())); return gatherFields_[i].fieldTag(); }

private:

  typedef typename panzer::Traits::Residual EvalT;
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  std::vector<std::string> indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const TpetraLinearObjContainer<double,LO,GO,NodeT> > tpetraContainer_;

  // Fields for storing tangent components dx/dp of solution vector x
  // These are not actually used by the residual specialization of this evaluator,
  // even if they are supplied, but it is useful to declare them as dependencies anyway
  // when saving the tangent components to the output file
  bool has_tangent_fields_;
  std::vector< std::vector< PHX::MDField<const ScalarT,Cell,NODE> > > tangentFields_;

  PHX::View<int**> scratch_lids_;
  std::vector<PHX::View<int*> > scratch_offsets_;

  GatherSolution_Tpetra();
};

// **************************************************************
// Tangent
// **************************************************************

/** \brief Tangent specialization of GatherSolution_Tpetra.
 *
 * Gathers solution values as in the Residual specialization, and
 * additionally exposes the tangent fields dx/dp (derivatives of the
 * solution with respect to parameters) so they can be carried through
 * a forward sensitivity (tangent) evaluation.
 */
template<typename TRAITS,typename LO,typename GO,typename NodeT>
class GatherSolution_Tpetra<panzer::Traits::Tangent,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Tangent, TRAITS>,
    public panzer::CloneableEvaluator  {


public:

  /// \brief Construct with a global indexer only; solution/parameter names must be set later (e.g. via clone()).
  GatherSolution_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer) :
     globalIndexer_(indexer) {}

  /// \brief Construct from a global indexer and a ParameterList giving the DOF names to gather (see the class-level Panzer_GatherSolution_Input options).
  GatherSolution_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                        const Teuchos::ParameterList& p);

  /// \brief Looks up and caches the Tpetra linear object container and field offsets needed by evaluateFields(). Called once before the first evaluation.
  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  /// \brief Fetches the Tpetra solution vector (and tangent vectors, if present) to gather from for the upcoming fill.
  void preEvaluate(typename TRAITS::PreEvalData d);

  /// \brief Copies solution and tangent values for this workset from the Tpetra vector(s) into #gatherFields_ / #tangentFields_.
  void evaluateFields(typename TRAITS::EvalData d);

  /// \brief Creates a copy of this evaluator configured from a new ParameterList, sharing the same global indexer.
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_Tpetra<panzer::Traits::Tangent,TRAITS,LO,GO,NodeT>(globalIndexer_,pl)); }

private:

  typedef typename panzer::Traits::Tangent EvalT;
  typedef typename panzer::Traits::Tangent::ScalarT ScalarT;
  typedef typename panzer::Traits::RealType RealT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;
  PHX::ViewOfViews<1,PHX::View<ScalarT**>> gatherFieldsVoV_;

  std::vector<std::string> indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const TpetraLinearObjContainer<double,LO,GO,NodeT> > tpetraContainer_;

  // Fields for storing tangent components dx/dp of solution vector x
  bool has_tangent_fields_;
  std::vector< std::vector< PHX::MDField<const RealT,Cell,NODE> > > tangentFields_;
  PHX::ViewOfViews<2,PHX::View<const RealT**>> tangentFieldsVoV_;
  PHX::View<size_t*> tangentInnerVectorSizes_;

  GatherSolution_Tpetra();
};

// **************************************************************
// Jacobian
// **************************************************************

/** \brief Jacobian specialization of GatherSolution_Tpetra.
 *
 * Gathers solution values from the Tpetra vector and seeds each one
 * as an independent AD variable (Sacado FAD), so that derivatives
 * with respect to the local solution unknowns propagate through the
 * rest of the field evaluation DAG to produce the Jacobian.
 */
template<typename TRAITS,typename LO,typename GO,typename NodeT>
class GatherSolution_Tpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>,
    public panzer::CloneableEvaluator  {

public:
  /// \brief Construct with a global indexer only; solution/parameter names must be set later (e.g. via clone()).
  GatherSolution_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer) :
     globalIndexer_(indexer) {}

  /// \brief Construct from a global indexer and a ParameterList giving the DOF names to gather (see the class-level Panzer_GatherSolution_Input options).
  GatherSolution_Tpetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                        const Teuchos::ParameterList& p);

  /// \brief Looks up and caches the Tpetra linear object container and field offsets needed by evaluateFields(). Called once before the first evaluation.
  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  /// \brief Fetches the Tpetra solution vector to gather from and determines whether sensitivities should be seeded for the upcoming fill.
  void preEvaluate(typename TRAITS::PreEvalData d);

  /// \brief Launches the device functor (operator()) over all cells in the workset to gather and, if enabled, seed solution values with derivatives.
  void evaluateFields(typename TRAITS::EvalData d);

  /// \brief Creates a copy of this evaluator configured from a new ParameterList, sharing the same global indexer.
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_Tpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>(globalIndexer_,pl)); }

  /// \brief Device functor: gathers one cell's solution values into #functor_data.field, seeding derivative components per #functor_data.seed_value.
  KOKKOS_INLINE_FUNCTION
  void operator()(const int cell) const;


  /// \brief Tag type selecting the operator() overload that gathers values without seeding any AD derivatives.
  struct NoSeed {};
  /// \brief Device functor: gathers one cell's solution values into #functor_data.field with no derivative seeding. \sa NoSeed
  KOKKOS_INLINE_FUNCTION
  void operator()(const NoSeed,const int cell) const;

private:

  typedef typename panzer::Traits::Jacobian EvalT;
  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  std::vector<std::string> indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  bool disableSensitivities_;     // This disables sensitivities absolutely
  std::string sensitivitiesName_; // This sets which gather operations have sensitivities
  bool applySensitivities_;       // This is a local variable that is used by evaluateFields
                                  // to turn on/off a certain set of sensitivities
  std::string globalDataKey_; // what global data does this fill?
  int gatherSeedIndex_; // what gather seed in the workset to use
                        // if less than zero then use alpha or beta
                        // as appropriate

  Teuchos::RCP<const TpetraLinearObjContainer<double,LO,GO,NodeT> > tpetraContainer_;
  Teuchos::RCP<typename TpetraLinearObjContainer<double,LO,GO,NodeT>::VectorType> x_vector;

  GatherSolution_Tpetra();

  PHX::View<int**> scratch_lids_;
  std::vector<PHX::View<int*> > scratch_offsets_;

  // functor data
  struct {
    // input values
    PHX::View<const LO**> lids;    // local indices for unknowns
    PHX::View<const int*> offsets; // how to get a particular field
    Kokkos::View<const double**, Kokkos::LayoutLeft,PHX::Device> x_data;
    double seed_value;                            // AD seed information
    int dos;	                                  // Offset for special interface bc

    // output fields
    PHX::MDField<ScalarT,Cell,NODE> field;
  } functor_data;
};

}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_GatherSolution_Tpetra_Hessian.hpp"
#endif

// **************************************************************
#endif
