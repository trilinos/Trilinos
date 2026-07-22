// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __mySourceTerm_hpp__
#define   __mySourceTerm_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <string>

// Panzer
#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

// Phalanx
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_FieldManager.hpp"

/**
 *  \brief Our source term.
 *
 *  This class represents our source term, \f$ \sin(2 \pi x) \sin(2 \pi y) \f$.
 */
template<typename EvalT, typename Traits>
class MySourceTerm
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    /**
     * 	\brief Default Constructor.
     *
     * 	This routine creates the field that will store the values of our source
     * 	term at all cells and integration points and tells Panzer that his
     * 	`Evaluator` will be evaluating that field.
     *
     * 	\param[in] name The name for our field.
     * 	\param[in] ir   The integration rule for our field.
     */
    MySourceTerm(
      const std::string&             name,
      const panzer::IntegrationRule& ir);

    /**
     * 	\brief Post-Registration Setup.
     *
     * 	This gets the index pointing to the integration rule degree in the
     * 	workset, which will be used later in grabbing the x- and y-locations
     * 	that correspond to integration points.
     *
     * 	\param[in] d  Essentially a list of `Workset`s, which are groups of
     * 	              cells (elements) that all reside on the same processor.
     * 	\param[in] fm This is unused in this routine, though it is a part of
     * 	              the `EvaluatorWithBaseImpl` interface.
     */
    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    /**
     * 	\brief Evaluate Fields.
     *
     * 	Loop over the cells in the workset, then over the integration points in
     * 	each cell, get the x- and y-locations corresponding to the integration
     * 	point, and compute our source term,
     * 	\f$ \sin(2 \pi x) \sin(2 \pi y) \f$.
     *
     * 	\param[in] d The `Workset` (collection of elements on a single
     * 	             processor) on which we're going to do some work.
     */
    void
    evaluateFields(
      typename Traits::EvalData d);

  private:
    typedef typename EvalT::ScalarT ScalarT;

    /**
     * 	\brief The field storing the values of our source term.
     */
    PHX::MDField<ScalarT, panzer::Cell, panzer::Point> result_;

    /**
     * 	\brief The degree of the integration rule.
     */
    int irDegree_;

    /**
     * 	\brief The index pointing to the integration rule degree in the
     * 	       workset.
     *
     * 	This is used in grabbing the x- and y-locations corresponding to the
     * 	integration points.
     */
    int irIndex_;
}; // end of class MySourceTerm

#endif // __mySourceTerm_hpp__
