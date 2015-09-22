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
