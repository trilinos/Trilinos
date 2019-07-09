// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef   __Panzer_GatherTangent_Epetra_decl_hpp__
#define   __Panzer_GatherTangent_Epetra_decl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"

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
   *  \brief GatherTangent_Epetra.
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
  class GatherTangent_Epetra
    :
    public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<EvalT, TRAITS>,
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
      GatherTangent_Epetra(
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
       *  \param[in] p       The input parameters.
       */
      GatherTangent_Epetra(
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
       *  \brief Pre-Evaluate:  Sets the tangent vector.
       *
       *  Sets the `GlobalEvaluationData` containing the owned and ghosted
       *  tangent vectors.
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
       *  Creates a `GatherTangent_Epetra` using the Initializing Constructor
       *  and the current object's `globalIndexer_`.
       *
       *  \param[in] pl The input parameters.
       *
       *  \returns A `GatherTangent_Epetra` constructed with this object's
       *           `globalIndexer_` and the input `ParameterList`.
       */
      virtual Teuchos::RCP<CloneableEvaluator>
      clone(
        const Teuchos::ParameterList& pl) const
      {
        using Teuchos::rcp;
        return rcp(new
          GatherTangent_Epetra<EvalT, TRAITS, LO, GO>(globalIndexer_, pl));
      } // end of clone()

    private:

      /**
       *  \brief The scalar type.
       */
      typedef typename EvalT::ScalarT ScalarT;

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
      Teuchos::RCP<std::vector<std::string>> indexerNames_;

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
       *  \brief The `GlobalEvaluationData` containing both the owned and
       *         ghosted tangent vectors.
       */
      Teuchos::RCP<panzer::EpetraVector_ReadOnly_GlobalEvaluationData>
      dxdpEvRoGed_;

      /**
       *  \brief Default Constructor (disabled).
       */
      GatherTangent_Epetra();

  }; // end of class GatherTangent_Epetra

} // end of namespace panzer

#endif // __Panzer_GatherTangent_Epetra_decl_hpp__
