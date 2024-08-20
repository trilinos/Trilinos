// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_GatherTangent_BlockedEpetra_decl_hpp__
#define   __Panzer_GatherTangent_BlockedEpetra_decl_hpp__

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

namespace panzer
{
  /**
   *  \brief GatherTangent_BlockedEpetra.
   *
   *  Gathers tangent vectors \f$\left(\frac{dx}{dp}\right)\f$ for computing
   *  \f$\frac{df}{dx}\frac{dx}{dp} + \frac{df}{dp}\f$ into the nodal fields of
   *  the field manager.
   *
   *  This evaluator is very similar to `GatherSolution`, however it always
   *  gathers into fields of type `double`, and it is a no-op if the global
   *  evaluation data container does not exist (which is an error for
   *  `GatherSolution`).
   *
   *  Currently makes an assumption that the stride is constant for DOFs and
   *  that the number of DOFs is equal to the size of the solution names
   *  vector.
   */
  template<typename EvalT, typename TRAITS, typename LO, typename GO>
  class GatherTangent_BlockedEpetra
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<EvalT, TRAITS>,
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
      GatherTangent_BlockedEpetra(
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
       *  \param[in] p        The input parameters.
       */
      GatherTangent_BlockedEpetra(
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
        PHX::FieldManager<TRAITS>& vm);

      /**
       *  \brief Pre-Evaluate:  Sets the tangent vector.
       *
       *  Sets the owned and ghosted tangent vectors.
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
       *  Creates a `GatherTangent_BlockedEpetra` using the Initializing
       *  Constructor and the current object's `indexers_`.
       *
       *  \param[in] pl The input parameters.
       *
       *  \returns A `GatherTangent_BlockedEpetra` constructed with this
       *           object's `indexers_` and the input `ParameterList`.
       */
      virtual Teuchos::RCP<CloneableEvaluator>
      clone(
        const Teuchos::ParameterList& pl) const
      {
        using Teuchos::rcp;
        return rcp(new
          GatherTangent_BlockedEpetra<EvalT, TRAITS, LO, GO>(indexers_, pl));
      } // end of clone()

    private:

      /**
       *  \brief The scalar type.
       */
      typedef typename EvalT::ScalarT ScalarT;

      /**
       *  \brief These map the local (field, element, basis) triplet to a
       *         global ID for scattering.
       */
      std::vector<Teuchos::RCP<const GlobalIndexer>> indexers_;

      /**
       *  \brief The block index into `indexers_`.
       */
      std::vector<int> indexerIds_;   // block index

      /**
       *  \brief Sub-field IDs, which need to be mapped.
       */
      std::vector<int> subFieldIds_; // sub field numbers

      /**
       *  \brief The fields to be gathered.
       */
      std::vector< PHX::MDField<ScalarT, Cell, NODE>> gatherFields_;

      /**
       *  \brief A list of the names of the fields to be gathered.
       */
      Teuchos::RCP<std::vector<std::string>> indexerNames_;

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
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted tangent vectors.
       */
      Teuchos::RCP<panzer::BlockedVector_ReadOnly_GlobalEvaluationData>
      xBvRoGed_;

      /**
       *  \brief Default Constructor (disabled)
       */
      GatherTangent_BlockedEpetra();

  }; // end of class GatherTangent_BlockedEpetra

} // end of namespace panzer

#endif // __Panzer_GatherTangent_BlockedEpetra_decl_hpp__
