// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   PANZER_ZEROCONTRIBUTEDFIELD_HPP
#define   PANZER_ZEROCONTRIBUTEDFIELD_HPP

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_Evaluator_WithBaseImpl.hpp"

// Phalanx
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer
{
  /**
   *  \brief Build a field using a specified data layout, and set it to zero.
   *
   *  Used to initialize a field to zero before other `Evaluator`s contribute
   *  to it.
   */
  template<typename EvalT, typename Traits>
  class ZeroContributedField
    :
    public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
  {
    public:

      /**
       *  \brief Constructor.
       *
       *  Given the field name and layout, create the field to be initialized
       *  to zero.
       *
       *  \param[in] fieldName The name of the field to be initialized to zero.
       *  \param[in] layout    The data layout to use when creating the field.
       */
      ZeroContributedField(
        const std::string& fieldName,
        PHX::DataLayout&   layout);

      /**
       *  \brief Evaluate the field.
       *
       *  Set the field to zero.
       */
      void
      evaluateFields(
        typename Traits::EvalData d) override;

    private:

      /**
       *  \brief The scalar data type.
       */
      using ScalarT = typename EvalT::ScalarT;

      /**
       *  \brief The field being initialized to zero.
       */
      PHX::MDField<ScalarT> field_;

  }; // end of class ZeroContributedField

} // end of namespace panzer

#endif // PANZER_ZEROCONTRIBUTEDFIELD_HPP
