// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_GatherSolution_Epetra_Hessian_hpp__
#define   __Panzer_GatherSolution_Epetra_Hessian_hpp__

// Only do this if required by the user.
#ifdef    Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main
// Epetra gather solution file

namespace panzer
{
  /**
   *  \brief GatherSolution_Epetra (Hessian Specialization)
   *
   *  Gathers solution values from the Newton solution vector into the nodal
   *  fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the number of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename TRAITS, typename LO, typename GO>
  class GatherSolution_Epetra<panzer::Traits::Hessian, TRAITS, LO, GO>
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Hessian, TRAITS>,
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
        const Teuchos::ParameterList& p);

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
        PHX::FieldManager<TRAITS>& vm);

      /**
       *  \brief Pre-Evaluate:  Sets the solution vector.
       *
       *  If using an `EpetraVector_ReadOnly_GlobalEvaluationData`, this sets
       *  the `GlobalEvaluationData`(s) containing both the owned and ghosted
       *  solution (and, if applicaple, derivative) vectors.  If using the
       *  older `EpeteraLinearObjContainer`, this sets the solution vector. 
       *  Also determines whether or not to apply sensitivities.
       *
       *  \param[in] d The `PreEvalData` containing the
       *               `GlobalEvaluationDataContainer` and the `first_` and
       *               `second_sensitivities_name`s.
       *
       *  \throws std::logic_error If it's unable to find the solution or (if
       *                           applicable) derivative vectors.
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
        return rcp(new GatherSolution_Epetra<Traits::Hessian, TRAITS, LO, GO>
          (globalIndexer_, pl));
      } // end of clone()

    private:

      /**
       *  \brief The evaluation type.
       */
      typedef typename panzer::Traits::Hessian EvalT;

      /**
       *  \brief The scalar type.
       */
      typedef typename panzer::Traits::Hessian::ScalarT ScalarT;

      /**
       *  \brief Maps the local (field, element, basis) triplet to a global ID
       *         for scattering.
       */
      Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;

      /**
       *  \brief A list of the names of the fields to be gathered.
       */
      std::vector<std::string> indexerNames_;

      /**
       *  \brief Field IDs, which need to be mapped.
       */
      std::vector<int> fieldIds_;

      /**
       *  \brief The fields to be gathered.
       */
      std::vector< PHX::MDField<ScalarT, Cell, NODE>> gatherFields_;

      /**
       *  \brief The sensitivity fields.
       */
      std::vector< PHX::MDField<ScalarT, Cell, NODE>> sensFields_;

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
       *  \brief Sets which gather operations have sensitivities.
       */
      std::string sensitivitiesName_;

      /**
       *  \brief Which gather seed in the workset to use.
       *
       *  If it's less than zero, then use alpha or beta as appropriate.
       */
      int gatherSeedIndex_;

      /**
       *  \brief A flag indicating whether or not we're to be working with the
       *         first derivative sensitivities.
       */
      bool firstSensitivitiesAvailable_;

      /**
       *  \brief Used by `evaluateFields()` to turn on/off the first derivative
       *         sensitivities.
       */
      bool firstApplySensitivities_;

      /**
       *  \brief The prefix for the field containing the second sensitivities.
       */
      std::string sensitivities2ndPrefix_;

      /**
       *  \brief A flag indicating whether or not we're to be working with the
       *         second derivative sensitivities.
       */
      bool secondSensitivitiesAvailable_;

      /**
       *  \brief Used by `evaluateFields()` to turn on/off the second
       *         derivative sensitivities.
       */
      bool secondApplySensitivities_;

      /**
       *  \brief The solution vector.
       */
      Teuchos::RCP<Epetra_Vector> x_;

      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted solution vectors.
       */
      Teuchos::RCP<panzer::EpetraVector_ReadOnly_GlobalEvaluationData>
      xEvRoGed_;

      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted derivative vectors.
       */
      Teuchos::RCP<panzer::EpetraVector_ReadOnly_GlobalEvaluationData>
      dxEvRoGed_;

      /**
       *  \brief Default Constructor (disabled).
       */
      GatherSolution_Epetra();

  }; // end of class GatherSolution_Epetra (Hessian Specialization)

} // end of namespace panzer

#endif // Panzer_BUILD_HESSIAN_SUPPORT

#endif // __Panzer_GatherSolution_Epetra_Hessian_hpp__
