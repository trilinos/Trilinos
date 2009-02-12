/*
// @HEADER
// ***********************************************************************
// 
//    GlobiPack: Collection of Scalar 1D globalizaton utilities
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef GLOBIPACK_BRENTS_LINE_SEARCH_DECL_HPP
#define GLOBIPACK_BRENTS_LINE_SEARCH_DECL_HPP


#include "GlobiPack_LineSearchBase.hpp"
#include "GlobiPack_GoldenQuadInterpBracket.hpp"
#include "GlobiPack_Brents1DMinimization.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace GlobiPack {


/** \brief Linesearch subclass implementing a function-value-only approximate
 * minimization algorithm using bracketing and then Brent's 1D minimization
 * method.
 *
 * This lineserach class is designed for more accurate linesearches and it
 * will march forward (as well as backward) from the given initial guess for
 * the step length in order to find it.  This linesearch is therefore more
 * appropriate for optimization algorithms like steppest decent and nonlinear
 * CG what require more accurate linesearches and where the scaling of the
 * step is not well know.  Also, this linesearch likely satisifies the Strong
 * Wolf Conditions but there are not checks for this at all so in the end it
 * may not.
 *
 * ToDo: Finish Documentation!
 */
template<typename Scalar>
class BrentsLineSearch
  : public LineSearchBase<Scalar>,
    protected Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \name Constructor/Initializers/Accessors */
  //@{

  /** \brief Construct with default parameters. */
  BrentsLineSearch();

  /** \brief For unit testing only . */
  const GoldenQuadInterpBracket<Scalar>& bracket() const;

  /** \brief For unit testing only . */
  const Brents1DMinimization<Scalar>& brentsMin() const;

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions) */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Overrridden from LineSearchBase. */
  //@{

  /** \brief Returns true. */
  virtual bool requiresBaseDeriv() const;

  /** \brief Returns false. */
  virtual bool requiresDerivEvals() const;

  /** \brief . */
  virtual bool doLineSearch(
    const MeritFunc1DBase<Scalar> &phi,
    const PointEval1D<Scalar> &point_k,
    const Ptr<PointEval1D<Scalar> > &point_kp1,
    const Ptr<int> &numIters
    ) const;

  //@}

private:

  // //////////////////////
  // Private data members

  GoldenQuadInterpBracket<Scalar> bracket_;
  Brents1DMinimization<Scalar> brentsMin_;

};


/** \brief Nonmember constructor.
 *
 * \relates BrentsLineSearch
 */
template<typename Scalar>
const RCP<BrentsLineSearch<Scalar> > brentsLineSearch()
{
  return Teuchos::rcp(new BrentsLineSearch<Scalar>());
}


// Default values are exposed here for unit testing purposes


namespace BrentsLineSearchUtils {


const std::string bracket_name = "Bracket";

const std::string minimize_name = "Minimize";


} // namespace BrentsLineSearchUtils



} // namespace GlobiPack


#endif // GLOBIPACK_BRENTS_LINE_SEARCH_DECL_HPP
