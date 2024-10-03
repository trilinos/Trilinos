// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_GatherSolution_BlockedEpetra_Hessian_hpp__
#define   __Panzer_GatherSolution_BlockedEpetra_Hessian_hpp__

// Only do this if required by the user.
#ifdef    Panzer_BUILD_HESSIAN_SUPPORT

namespace panzer
{
  /**
   *  \brief `GatherSolution_BlockedEpetra` (Hessian Specialization).
   *
   *  Gathers solution values from the Newton solution vector into the nodal
   *  fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the nmber of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename TRAITS, typename LO, typename GO>
  class GatherSolution_BlockedEpetra<panzer::Traits::Hessian, TRAITS, LO, GO>
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Hessian, TRAITS>,
    public panzer::CloneableEvaluator
  {
    public:
    
      /**
       *  \brief Constructor.
       *
       *  Simply saves the input `indexers` as this object's `indexers_`.
       *
       *  \param[in] indexers The `vector` of `GlobalIndexer`s that
       *                      handle the global unknown numbering.
       */
      GatherSolution_BlockedEpetra(
        const std::vector<Teuchos::RCP<const GlobalIndexer<LO, int>>>&
          indexers)
        :
        indexers_(indexers)
      {
      } // end of Constructor
    
      /**
       *  \brief Initializing Constructor.
       *
       *  Saves the input `indexers` as this object's `indexers_`, allocates
       *  fields, sets up dependent tangent fields (if requested), and
       *  determines the first active name.
       *
       *  \param[in] indexers The `vector` of `GlobalIndexer`s that
       *                      handle the global unknown numbering.
       *  \param[in] p        A `ParameterList` used as input for
       *                      `GatherSolution_Input`.
       */
      GatherSolution_BlockedEpetra(
        const std::vector<Teuchos::RCP<const GlobalIndexer<LO, int>>>&
          indexers,
        const Teuchos::ParameterList& p);
    
      /**
       *  \brief Post-Registration Setup.
       *
       *  Loops over the `gatherFields_` and sets the `indexerIds_` and
       *  `subFieldIds_`.
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
       *  If using a `BlockedVector_ReadOnly_GlobalEvaluationData`, this saves
       *  it for use later in `evaluateFields()`.  If using the older
       *  `BlockedEpetraLinearObjContainer`, this sets the solution vector.
       *
       *  \param[in] d The `PreEvalData` containing the
       *               `GlobalEvaluationDataContainer`.
       */
      void
      preEvaluate(
        typename TRAITS::PreEvalData d);
    
      /**
       *  \brief Evaluate Fields:  Gather operation.
       *
       *  Loops over the fields to be gathered, the cells in the workset, and
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
       *  Creates a `GatherSolution_BlockedEpetra` using the Initializing
       *  Constructor and the current object's `indexers_`.
       *
       *  \param[in] pl A `ParameterList` used as input for
       *                `GatherSolution_Input`.
       *
       *  \returns A `GatherSolution_BlockedEpetra` constructed with this
       *           object's `indexers_` and the input `ParameterList`.
       */
      virtual Teuchos::RCP<CloneableEvaluator>
      clone(
        const Teuchos::ParameterList& pl) const
      {
        using panzer::GatherSolution_BlockedEpetra;
        using panzer::Traits;
        using Teuchos::rcp;
        return rcp(new
          GatherSolution_BlockedEpetra<Traits::Hessian, TRAITS, LO, GO>
          (indexers_, pl));
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
       *  \brief A list of the names of the fields to be gathered.
       */
      std::vector<std::string> indexerNames_;
    
      /**
       *  \brief The key identifying the `GlobalEvaluationData`.
       */
      std::string globalDataKey_;
    
      /**
       *  \brief These map the local (field, element, basis) triplet to a
       *         global ID for scattering.
       */
      std::vector<Teuchos::RCP<const GlobalIndexer<LO, int>>> indexers_;

      /**
       *  \brief The block index into `indexers_`.
       */
      std::vector<int> indexerIds_;

      /**
       *  \brief Sub-field IDs, which need to be mapped.
       */
      std::vector<int> subFieldIds_;
    
      /**
       *  \brief The fields to be gathered.
       */
      std::vector< PHX::MDField<ScalarT, Cell, NODE>> gatherFields_;
    
      /**
       *  \brief A flag indicating whether we should be working with \f$ x \f$
       *         or \f$ \dot{x} \f$.
       */
      bool useTimeDerivativeSolutionVector_;
    
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
       *  \brief A flag indicating whether or not first sensitivities
       *         information is available.
       */
      bool firstSensitivitiesAvailable_;

      /**
       *  \brief Used by `evaluateFields()` to turn on/off the first
       *         sensitivities.
       */
      bool firstApplySensitivities_;
                            
      /**
       *  \brief The prefix for the field containing the second sensitivities.
       */
      std::string sensitivities2ndPrefix_;

      /**
       *  \brief A flag indicating whether or not second sensitivities
       *         information is available.
       */
      bool secondSensitivitiesAvailable_;

      /**
       *  \brief Used by `evaluateFields()` to turn on/off the second
       *         sensitivities.
       */
      bool secondApplySensitivities_;
    
      /**
       *  \brief The solution vector.
       */
      Teuchos::RCP<Thyra::ProductVectorBase<double>> x_;

      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted solution vectors.
       */
      Teuchos::RCP<panzer::BlockedVector_ReadOnly_GlobalEvaluationData>
      xBvRoGed_;

      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted derivative vectors.
       */
      Teuchos::RCP<panzer::BlockedVector_ReadOnly_GlobalEvaluationData>
      dxBvRoGed_;
    
      /**
       *  \brief Default Constructor (disabled)
       */
      GatherSolution_BlockedEpetra();

  }; // end of class GatherSolution_BlockedEpetra (Hessian Specialization)

} // end of namespace panzer

#endif // Panzer_BUILD_HESSIAN_SUPPORT

#endif // __Panzer_GatherSolution_BlockedEpetra_Hessian_hpp__
