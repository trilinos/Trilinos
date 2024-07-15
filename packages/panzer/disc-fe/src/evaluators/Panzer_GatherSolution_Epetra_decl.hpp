// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_GatherSolution_Epetra_decl_hpp__
#define   __Panzer_GatherSolution_Epetra_decl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_Traits.hpp"

// Phalanx
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Forward Declarations
//
///////////////////////////////////////////////////////////////////////////////

namespace panzer
{
  class GlobalIndexer;
}

namespace panzer
{
  /**
   *  \brief Gathers solution values from the Newton solution vector into the
   *         nodal fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the number of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename EvalT, typename TRAITS, typename LO, typename GO>
  class GatherSolution_Epetra;

  /**
   *  \brief GatherSolution_Epetra (Residual Specialization).
   *
   *  Gathers solution values from the Newton solution vector into the nodal
   *  fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the number of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename TRAITS, typename LO, typename GO>
  class GatherSolution_Epetra<panzer::Traits::Residual, TRAITS, LO, GO>
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator
  {
    public:

      /**
       *  \brief Constructor.
       *
       *  Simply saves the input `indexer` as this object's `globalIndexer_`.
       *
       *  \param[in] indexer The `GlobalIndexer` that handles the global
       *                     unknown numbering.
       */
      GatherSolution_Epetra(
        const Teuchos::RCP<const panzer::GlobalIndexer>& indexer)
        :
        globalIndexer_(indexer)
      {
      } // end of Constructor

      /**
       *  \brief Initializing Constructor.
       *
       *  Saves the input `indexer` as this object's `globalIndexer_`,
       *  allocates fields, sets up dependent tangent fields (if requested),
       *  and determines the first active name.
       *
       *  \param[in] indexer The `GlobalIndexer` that handles the global
       *                     unknown numbering.
       *  \param[in] p       A `ParameterList` used as input for
       *                     `GatherSolution_Input`.
       */
      GatherSolution_Epetra(
        const Teuchos::RCP<const panzer::GlobalIndexer>& indexer,
        const Teuchos::ParameterList&                                  p);

      /**
       *  \brief Post-Registration Setup.
       *
       *  Loops over the `gatherFields_` and sets the `fieldIds_`.
       *
       *  \param[in] d  Unused.
       *  \param[in] fm Unused.
       */
      void
      postRegistrationSetup(
        typename TRAITS::SetupData d,
        PHX::FieldManager<TRAITS>& fm);

      /**
       *  \brief Pre-Evaluate:  Sets the solution vector.
       *
       *  If using an `EpetraVector_ReadOnly_GlobalEvaluationData`, this sets
       *  the `GlobalEvaluationData` containing both the owned and ghosted
       *  solution vectors.  If using the older `EpetraLinearObjContainer`,
       *  this sets the solution vector.
       *
       *  \param[in] d The `PreEvalData` containing the
       *               `GlobalEvaluationDataContainer` and the `first_` and
       *               `second_sensitivities_name`s.
       */
      void
      preEvaluate(
        typename TRAITS::PreEvalData d);

      /**
       *  \brief Evaluate Fields:  Gather operation.
       *
       *  Loops over the cells in the workset, the fields to be gathered, and
       *  the basis functions, and fills in the fields.
       *
       *  \param[in] d The `Workset` on which we're going to do all the work.
       */
      void
      evaluateFields(
        typename TRAITS::EvalData d);

      /**
       *  \brief Create a copy.
       *
       *  Creates a `GatherSolution_Epetra` using the Initializing Constructor
       *  and the current object's `globalIndexer_`.
       *
       *  \param[in] pl A `ParameterList` used as input for
       *                `GatherSolution_Input`.
       *
       *  \returns A `GatherSolution_Epetra` constructed with this object's
       *           `globalIndexer_` and the input `ParameterList`.
       */
      virtual Teuchos::RCP<CloneableEvaluator>
      clone(
        const Teuchos::ParameterList& pl) const
      {
        using panzer::Traits;
        using Teuchos::rcp;
        return rcp(new GatherSolution_Epetra<Traits::Residual, TRAITS, LO, GO>
          (globalIndexer_, pl));
      } // end of clone()

    private:

      /**
       *  \brief The evaluation type.
       */
      typedef typename panzer::Traits::Residual EvalT;

      /**
       *  \brief The scalar type.
       */
      typedef typename panzer::Traits::Residual::ScalarT ScalarT;

      /**
       *  \brief Maps the local (field, element, basis) triplet to a global ID
       *         for scattering.
       */
      Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;

      /**
       *  \brief Field IDs, which need to be mapped.
       */
      std::vector<int> fieldIds_;

      /**
       *  \brief The fields to be gathered.
       */
      std::vector<PHX::MDField<ScalarT, Cell, NODE>> gatherFields_;

      /**
       *  \brief A list of the names of the fields to be gathered.
       */
      std::vector<std::string> indexerNames_;

      /**
       *  \brief A flag indicating whether we're to be working with \f$ x \f$
       *         or \f$ \dot{x} \f$.
       */
      bool useTimeDerivativeSolutionVector_;

      /**
       *  \brief The key identifying the `GlobalEvaluationData`.
       */
      std::string globalDataKey_;

      /**
       *  \brief The solution vector.
       */
      Teuchos::RCP<Epetra_Vector> x_;

      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted vectors.
       */
      Teuchos::RCP<panzer::EpetraVector_ReadOnly_GlobalEvaluationData>
      xEvRoGed_;

      /**
       *  \brief A flag indicating whether or not we have tangent fields.
       */
      bool hasTangentFields_;

      /**
       *  \brief Fields for storing the tangent components
       *         \f$\left(\frac{dx}{dp}\right)\f$ of the solution vector
       *         \f$(x)\f$.
       *
       *  These are not actually used by the residual specialization of this
       *  evaluator, even if they are supplied, but it is useful to declare
       *  them as dependencies anyway when saving the tangent components to the
       *  output file.
       */
      std::vector<std::vector<PHX::MDField<const ScalarT, Cell, NODE>>>
        tangentFields_;

      /**
       *  \brief Default Constructor (disabled).
       */
      GatherSolution_Epetra();

  }; // end of class GatherSolution_Epetra (Residual specialization)

  /**
   *  \brief GatherSolution_Epetra (Tangent Specialization).
   *
   *  Gathers solution values from the Newton solution vector into the nodal
   *  fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the number of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename TRAITS, typename LO, typename GO>
  class GatherSolution_Epetra<panzer::Traits::Tangent, TRAITS, LO, GO>
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Tangent, TRAITS>,
    public panzer::CloneableEvaluator
  {
    public:

      /**
       *  \brief Constructor.
       *
       *  Simply saves the input `indexer` as this object's `globalIndexer_`.
       *
       *  \param[in] indexer The `GlobalIndexer` that handles the global
       *                     unknown numbering.
       */
      GatherSolution_Epetra(
        const Teuchos::RCP<const panzer::GlobalIndexer>& indexer)
        :
        globalIndexer_(indexer)
      {
      } // end of Constructor

      /**
       *  \brief Initializing Constructor.
       *
       *  Saves the input `indexer` as this object's `globalIndexer_`,
       *  allocates fields, sets up dependent tangent fields (if requested),
       *  and determines the first active name.
       *
       *  \param[in] indexer The `GlobalIndexer` that handles the global
       *                     unknown numbering.
       *  \param[in] p       A `ParameterList` used as input for
       *                     `GatherSolution_Input`.
       */
      GatherSolution_Epetra(
        const Teuchos::RCP<const panzer::GlobalIndexer>& indexer,
        const Teuchos::ParameterList&                                  p);

      /**
       *  \brief Post-Registration Setup.
       *
       *  Loops over the `gatherFields_` and sets the `fieldIds_`.
       *
       *  \param[in] d  Unused.
       *  \param[in] fm Unused.
       */
      void
      postRegistrationSetup(
        typename TRAITS::SetupData d,
        PHX::FieldManager<TRAITS>& fm);

      /**
       *  \brief Pre-Evaluate:  Sets the solution vector.
       *
       *  If using an `EpetraVector_ReadOnly_GlobalEvaluationData`, this sets
       *  the `GlobalEvaluationData` containing both the owned and ghosted
       *  solution vectors.  If using the older `EpetraLinearObjContainer`,
       *  this sets the solution vector.
       *
       *  \param[in] d The `PreEvalData` containing the
       *               `GlobalEvaluationDataContainer` and the `first_` and
       *               `second_sensitivities_name`s.
       */
      void
      preEvaluate(
        typename TRAITS::PreEvalData d);

      /**
       *  \brief Evaluate Fields:  Gather operation.
       *
       *  Loops over the cells in the workset, the fields to be gathered, and
       *  the basis functions, and fills in the fields.  Also sets derivative
       *  information if tangent fields are enabled.
       *
       *  \param[in] d The `Workset` on which we're going to do all the work.
       */
      void
      evaluateFields(
        typename TRAITS::EvalData d);

      /**
       *  \brief Create a copy.
       *
       *  Creates a `GatherSolution_Epetra` using the Initializing Constructor
       *  and the current object's `globalIndexer_`.
       *
       *  \param[in] pl A `ParameterList` used as input for
       *                `GatherSolution_Input`.
       *
       *  \returns A `GatherSolution_Epetra` constructed with this object's
       *           `globalIndexer_` and the input `ParameterList`.
       */
      virtual Teuchos::RCP<CloneableEvaluator>
      clone(
        const Teuchos::ParameterList& pl) const
      {
        using panzer::Traits;
        using Teuchos::rcp;
        return rcp(new GatherSolution_Epetra<Traits::Tangent, TRAITS, LO, GO>
          (globalIndexer_, pl));
      } // end of clone()

    private:

      /**
       *  \brief The evaluation type.
       */
      typedef typename panzer::Traits::Tangent EvalT;

      /**
       *  \brief The scalar type.
       */
      typedef typename panzer::Traits::Tangent::ScalarT ScalarT;

      /**
       *  \brief Maps the local (field, element, basis) triplet to a global ID
       *         for scattering.
       */
      Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;

      /**
       *  \brief Field IDs, which need to be mapped.
       */
      std::vector<int> fieldIds_;

      /**
       *  \brief The fields to be gathered.
       */
      std::vector<PHX::MDField<ScalarT, Cell, NODE>> gatherFields_;

      /**
       *  \brief A list of the names of the fields to be gathered.
       */
      std::vector<std::string> indexerNames_;

      /**
       *  \brief A flag indicating whether we should be working with \f$ x \f$
       *         or \f$ \dot{x} \f$.
       */
      bool useTimeDerivativeSolutionVector_;

      /**
       *  \brief The key identifying the `GlobalEvaluationData`.
       */
      std::string globalDataKey_;

      /**
       *  \brief The solution vector.
       */
      Teuchos::RCP<Epetra_Vector> x_;

      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted vectors.
       */
      Teuchos::RCP<panzer::EpetraVector_ReadOnly_GlobalEvaluationData>
      xEvRoGed_;

      /**
       *  \brief A flag indicating whether or not we have tangent fields.
       */
      bool hasTangentFields_;

      /**
       *  \brief Fields for storing the tangent components
       *         \f$\left(\frac{dx}{dp}\right)\f$ of the solution vector
       *         \f$(x)\f$.
       */
      std::vector<std::vector<PHX::MDField<const ScalarT, Cell, NODE>>>
        tangentFields_;

      /**
       *  \brief Default Constructor (disabled).
       */
      GatherSolution_Epetra();

  }; // end of class GatherSolution_Epetra (Tangent specialization)

  /**
   *  \brief GatherSolution_Epetra (Jacobian Specialization).
   *
   *  Gathers solution values from the Newton solution vector into the nodal
   *  fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the number of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename TRAITS, typename LO, typename GO>
  class GatherSolution_Epetra<panzer::Traits::Jacobian, TRAITS, LO, GO>
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>,
    public panzer::CloneableEvaluator
  {
    public:

      /**
       *  \brief Constructor.
       *
       *  Simply saves the input `indexer` as this object's `globalIndexer_`.
       *
       *  \param[in] indexer The `GlobalIndexer` that handles the global
       *                     unknown numbering.
       */
      GatherSolution_Epetra(
        const Teuchos::RCP<const panzer::GlobalIndexer>& indexer)
        :
        globalIndexer_(indexer)
      {
      } // end of Constructor

      /**
       *  \brief Initializing Constructor.
       *
       *  Saves the input `indexer` as this object's `globalIndexer_`,
       *  allocates fields, and determines the first active name.
       *
       *  \param[in] indexer The `GlobalIndexer` that handles the global
       *                     unknown numbering.
       *  \param[in] p       A `ParameterList` used as input for
       *                     `GatherSolution_Input`.
       */
      GatherSolution_Epetra(
        const Teuchos::RCP<const panzer::GlobalIndexer>& indexer,
        const Teuchos::ParameterList&                                  p);

      /**
       *  \brief Post-Registration Setup.
       *
       *  Loops over the `gatherFields_` and sets the `fieldIds_`.
       *
       *  \param[in] d  Unused.
       *  \param[in] fm Unused.
       */
      void
      postRegistrationSetup(
        typename TRAITS::SetupData d,
        PHX::FieldManager<TRAITS>& fm);

      /**
       *  \brief Pre-Evaluate:  Sets the solution vector.
       *
       *  If using an `EpetraVector_ReadOnly_GlobalEvaluationData`, this sets
       *  the `GlobalEvaluationData` containing both the owned and ghosted
       *  solution vectors.  If using the older `EpetraLinearObjContainer`,
       *  this sets the solution vector.  Also determines whether or not to
       *  apply sensitivities.
       *
       *  \param[in] d The `PreEvalData` containing the
       *               `GlobalEvaluationDataContainer` and the `first_` and
       *               `second_sensitivities_name`s.
       */
      void
      preEvaluate(
        typename TRAITS::PreEvalData d);

      /**
       *  \brief Evaluate Fields:  Gather operation.
       *
       *  Loops over the cells in the workset, the fields to be gathered, and
       *  the basis functions, and fills in the fields.  If sensitivities are
       *  to be applied, this also seeds the derivatives.
       *
       *  \param[in] d The `Workset` on which we're going to do all the work.
       */
      void
      evaluateFields(
        typename TRAITS::EvalData d);

      /**
       *  \brief Create a copy.
       *
       *  Creates a `GatherSolution_Epetra` using the Initializing Constructor
       *  and the current object's `globalIndexer_`.
       *
       *  \param[in] pl A `ParameterList` used as input for
       *                `GatherSolution_Input`.
       *
       *  \returns A `GatherSolution_Epetra` constructed with this object's
       *           `globalIndexer_` and the input `ParameterList`.
       */
      virtual Teuchos::RCP<CloneableEvaluator>
      clone(
        const Teuchos::ParameterList& pl) const
      {
        using panzer::Traits;
        using Teuchos::rcp;
        return rcp(new GatherSolution_Epetra<Traits::Jacobian, TRAITS, LO, GO>
          (globalIndexer_, pl));
      } // end of clone()

    private:

      /**
       *  \brief The evaluation type.
       */
      typedef typename panzer::Traits::Jacobian EvalT;

      /**
       *  \brief The scalar type.
       */
      typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;

      /**
       *  \brief Maps the local (field, element, basis) triplet to a global ID
       *         for scattering.
       */
      Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;

      /**
       *  \brief Field IDs, which need to be mapped.
       */
      std::vector<int> fieldIds_;

      /**
       *  \brief The fields to be gathered.
       */
      std::vector<PHX::MDField<ScalarT, Cell, NODE>> gatherFields_;

      /**
       *  \brief A list of the names of the fields to be gathered.
       */
      std::vector<std::string> indexerNames_;

      /**
       *  \brief A flag indicating whether we're to be working with \f$ x \f$
       *         or \f$ \dot{x} \f$.
       */
      bool useTimeDerivativeSolutionVector_;

      /**
       *  \brief Flag to disable sensitivities absolutely.
       */
      bool disableSensitivities_;

      /**
       *  \brief Sets which gather operations have sensitivities.
       */
      std::string sensitivitiesName_;

      /**
       *  \brief Used by `evaluateFields()` to turn on/off a certain set of
       *         sensitivities.
       */
      bool applySensitivities_;

      /**
       *  \brief The key identifying the `GlobalEvaluationData`.
       */
      std::string globalDataKey_;

      /**
       *  \brief Which gather seed in the workset to use.
       *
       *  If it's less than zero, then use alpha or beta as appropriate.
       */
      int gatherSeedIndex_;

      /**
       *  \brief The solution vector.
       */
      Teuchos::RCP<Epetra_Vector> x_;

      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted vectors.
       */
      Teuchos::RCP<panzer::EpetraVector_ReadOnly_GlobalEvaluationData>
      xEvRoGed_;

      /**
       *  \brief Default Constructor (disabled).
       */
      GatherSolution_Epetra();

  }; // end of class GatherSolution_Epetra (Jacobian specialization)

} // end of namespace panzer

#ifdef    Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_GatherSolution_Epetra_Hessian.hpp"
#endif // Panzer_BUILD_HESSIAN_SUPPORT

#endif // __Panzer_GatherSolution_Epetra_decl_hpp__
