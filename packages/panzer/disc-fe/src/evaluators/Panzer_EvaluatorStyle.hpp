// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_EvaluatorStyle_hpp__
#define   __Panzer_EvaluatorStyle_hpp__

namespace panzer
{
  /**
   *  \brief An indication of how an `Evaluator` will behave.
   *
   *  An `Evaluator` will compute the result of its evaluation and then
   *  behave according to the table below.
   */
  enum class EvaluatorStyle
  {
    CONTRIBUTES, /*!< Contribute the result to a specified residual, not saving
                      anything. */
    EVALUATES,   /*!< Save the result under a specified name for future use. */
  }; // end of enum class EnumeratorStyle

} // end of namespace panzer

#endif // __Panzer_EvaluatorStyle_hpp__
