// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_GatherSolution_BlockedEpetra_decl_hpp__
#define   __Panzer_GatherSolution_BlockedEpetra_decl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_Dimension.hpp"
#include "PanzerDiscFE_config.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_Traits.hpp"

// Phalanx
#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

// Teuchos
#include "Teuchos_ParameterList.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Forward Declarations
//
///////////////////////////////////////////////////////////////////////////////

namespace panzer
{
  class GlobalIndexer;
}
namespace Thyra
{
  template<typename> class ProductVectorBase;
}

namespace panzer
{
  /**
   *  \brief Gathers solution values from the Newton solution vector into the
   *         nodal fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the nmber of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename EvalT, typename TRAITS, typename LO, typename GO>
  class GatherSolution_BlockedEpetra
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator
  {
    public:

      /**
       *  \brief The scalar type.
       */
      typedef typename EvalT::ScalarT ScalarT;
  
      /**
       *  \brief Constructor.
       *
       *  For this unspecialized class, this constructor is empty; that is, it
       *  doesn't do anything with the input `ParameterList`.
       *
       *  \param[in] p A `ParameterList` that isn't used.
       */
      GatherSolution_BlockedEpetra(
        const Teuchos::ParameterList& p)
      {
      } // end of Constructor
  
      /**
       *  \brief Create a copy.
       *
       *  For this unspecialized class, this actually just calls the Default
       *  Constructor and doesn't copy anything.
       *
       *  \param[in] pl A `ParameterList`, which is passed to the Default
       *                Constructor.
       *
       *  \returns A `GatherSolution_BlockedEpetra` created by the Default
       *           Constructor.
       */
      virtual Teuchos::RCP<CloneableEvaluator>
      clone(
        const Teuchos::ParameterList& pl) const
      {
        return Teuchos::rcp(new
          GatherSolution_BlockedEpetra<EvalT, TRAITS, LO, GO>(pl));
      } // end of clone()
  
      /**
       *  \brief Post-Registration Setup.
       *
       *  For this unspecialized class, this routine does nothing.
       *
       *  \param[in] d  Unused.
       *  \param[in] fm Unused.
       */
      void
      postRegistrationSetup(
        typename TRAITS::SetupData d,
        PHX::FieldManager<TRAITS>& fm)
      {
      } // end of postRegistrationSetup()
  
      /**
       *  \brief Evaluate Fields.
       *
       *  For this unspecialized class, this routine does nothing, other than
       *  tell you that you can't use it.
       *
       *  \param[in] d Unused.
       */
      void
      evaluateFields(
        typename TRAITS::EvalData d)
      {
        using PHX::print;
        using std::cout;
        using std::endl;
        cout << "Unspecialized version of \"GatherSolution_BlockedEpetra::"   \
          "evaluateFields\" on " + print<EvalT>() + "\" should not "   \
          "be used!" << endl;
        TEUCHOS_ASSERT(false);
      } // end of evaluateFields()

  }; // end of class GatherSolution_BlockedEpetra

  /**
   *  \brief GatherSolution_BlockedEpetra (Residual Specialization).
   *
   *  Gathers solution values from the Newton solution vector into the nodal
   *  fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the number of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename TRAITS, typename LO, typename GO>
  class GatherSolution_BlockedEpetra<panzer::Traits::Residual, TRAITS, LO, GO>
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
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
        const std::vector<Teuchos::RCP<const GlobalIndexer>>&
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
        const std::vector<Teuchos::RCP<const GlobalIndexer>>&
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
        PHX::FieldManager<TRAITS>& fm);
       
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
        using panzer::Traits;
        using Teuchos::rcp;
        return rcp(new
          GatherSolution_BlockedEpetra<Traits::Residual, TRAITS, LO, GO>
          (indexers_, pl));
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
       *  \brief These map the local (field, element, basis) triplet to a
       *         global ID for scattering.
       */
      std::vector<Teuchos::RCP<const GlobalIndexer>> indexers_;
       
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
      Teuchos::RCP<Thyra::ProductVectorBase<double>> x_;
       
      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted solution vectors.
       */
      Teuchos::RCP<panzer::BlockedVector_ReadOnly_GlobalEvaluationData>
      xBvRoGed_;
  
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
       *  \brief Default Constructor (disabled)
       */
      GatherSolution_BlockedEpetra();

  }; // end of class GatherSolution_BlockedEpetra (Residual Specialization)

  /**
   *  \brief GatherSolution_BlockedEpetra (Tangent Specialization).
   *
   *  Gathers solution values from the Newton solution vector into the nodal
   *  fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the nmber of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename TRAITS, typename LO, typename GO>
  class GatherSolution_BlockedEpetra<panzer::Traits::Tangent, TRAITS, LO, GO>
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Tangent, TRAITS>,
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
        const std::vector<Teuchos::RCP<const GlobalIndexer>>&
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
        const std::vector<Teuchos::RCP<const GlobalIndexer>>&
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
        PHX::FieldManager<TRAITS>& fm);
  
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
        using panzer::Traits;
        using Teuchos::rcp;
        return rcp(new
          GatherSolution_BlockedEpetra<Traits::Tangent, TRAITS, LO, GO>
          (indexers_, pl));
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
       *  \brief These map the local (field, element, basis) triplet to a
       *         global ID for scattering.
       */
      std::vector<Teuchos::RCP<const GlobalIndexer>> indexers_;
  
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
      Teuchos::RCP<Thyra::ProductVectorBase<double>> x_;
  
      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted solution vectors.
       */
      Teuchos::RCP<panzer::BlockedVector_ReadOnly_GlobalEvaluationData>
      xBvRoGed_;

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
       *  \brief Default Constructor (disabled)
       */
      GatherSolution_BlockedEpetra();

  }; // end of class GatherSolution_BlockedEpetra (Tangent Specialization)

  /**
   *  \brief GatherSolution_BlockedEpetra (Jacobian Specialization).
   *
   *  Gathers solution values from the Newton solution vector into the nodal
   *  fields of the field manager.
   *
   *  Currently makes an assumption that the stride is constant for degrees of
   *  freedom (DOFs) and that the nmber of DOFs is equal to the size of the
   *  solution names vector.
   */
  template<typename TRAITS, typename LO, typename GO>
  class GatherSolution_BlockedEpetra<panzer::Traits::Jacobian, TRAITS, LO, GO>
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>,
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
        const std::vector<Teuchos::RCP<const GlobalIndexer>>&
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
        const std::vector<Teuchos::RCP<const GlobalIndexer>>&
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
        PHX::FieldManager<TRAITS>& fm);

      /**
       *  \brief Pre-Evaluate:  Sets the solution vector.
       *
       *  If using a `BlockedVector_ReadOnly_GlobalEvaluationData`, this saves
       *  it for use later in `evaluateFields()`.  If using the older
       *  `BlockedEpetraLinearObjContainer`, this sets the solution vector.
       *  Also determines whether or not to apply sensitivities.
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
        using panzer::Traits;
        using Teuchos::rcp;
        return rcp(new
          GatherSolution_BlockedEpetra<Traits::Jacobian, TRAITS, LO, GO>
          (indexers_, pl));
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
       *  \brief These map the local (field, element, basis) triplet to a
       *         global ID for scattering.
       */
      std::vector<Teuchos::RCP<const GlobalIndexer>> indexers_;

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
      Teuchos::RCP<Thyra::ProductVectorBase<double>> x_;

      /**
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted solution vectors.
       */
      Teuchos::RCP<panzer::BlockedVector_ReadOnly_GlobalEvaluationData>
      xBvRoGed_;

      /**
       *  \brief Default Constructor (disabled)
       */
      GatherSolution_BlockedEpetra();

  }; // end of class GatherSolution_BlockedEpetra (Jacobian Specialization)

} // end of namespace panzer

#ifdef    Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_GatherSolution_BlockedEpetra_Hessian.hpp"
#endif // Panzer_BUILD_HESSIAN_SUPPORT

#endif // __Panzer_GatherSolution_BlockedEpetra_decl_hpp__
