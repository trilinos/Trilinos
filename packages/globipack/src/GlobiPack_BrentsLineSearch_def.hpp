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

#ifndef GLOBIPACK_BRENTS_LINE_SEARCH_DEF_HPP
#define GLOBIPACK_BRENTS_LINE_SEARCH_DEF_HPP


#include "GlobiPack_BrentsLineSearch_decl.hpp"
#include "Teuchos_TabularOutputter.hpp"


namespace GlobiPack {


// Constructor/Initializers/Accessors


template<typename Scalar>
BrentsLineSearch<Scalar>::BrentsLineSearch()
{}


template<typename Scalar>
const GoldenQuadInterpBracket<Scalar>&
BrentsLineSearch<Scalar>::bracket() const
{
  return bracket_;
}


template<typename Scalar>
const Brents1DMinimization<Scalar>&
BrentsLineSearch<Scalar>::brentsMin() const
{
  return brentsMin_;
}


// Overridden from ParameterListAcceptor


template<class Scalar>
void BrentsLineSearch<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  //typedef ScalarTraits<Scalar> ST; // unused
  namespace BLSU = BrentsLineSearchUtils;
  using Teuchos::sublist;
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  bracket_.setParameterList(sublist(paramList, BLSU::bracket_name, true));
  brentsMin_.setParameterList(sublist(paramList, BLSU::minimize_name, true));
  setMyParamList(paramList);
}


template<class Scalar>
RCP<const ParameterList>
BrentsLineSearch<Scalar>::getValidParameters() const
{
  namespace BLSU = BrentsLineSearchUtils;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList());
    pl->sublist(BLSU::bracket_name).setParameters(
      *bracket_.getValidParameters()
      ).disableRecursiveValidation();
    pl->sublist(BLSU::minimize_name).setParameters(
      *brentsMin_.getValidParameters()
      ).disableRecursiveValidation();
    validPL = pl;
  }
  return validPL;
}


// Overrridden from LineSearchBase


template<typename Scalar>
bool BrentsLineSearch<Scalar>::requiresBaseDeriv() const
{
  return false;
}


template<typename Scalar>
bool BrentsLineSearch<Scalar>::requiresDerivEvals() const
{
  return false;
}


template<typename Scalar>
bool BrentsLineSearch<Scalar>::doLineSearch(
  const MeritFunc1DBase<Scalar> &phi,
  const PointEval1D<Scalar> &point_k,
  const Ptr<PointEval1D<Scalar> > &point_kp1,
  const Ptr<int> &numIters
  ) const
{

  using Teuchos::as;
  using Teuchos::OSTab;
  using Teuchos::outArg;
  using Teuchos::inOutArg;

#ifdef TEUCHOS_DEBUG
  typedef ScalarTraits<Scalar> ST;
  typedef PointEval1D<Scalar> PE1D;

  TEUCHOS_ASSERT_EQUALITY(point_k.alpha, ST::zero());
  TEUCHOS_ASSERT_INEQUALITY(point_k.phi, !=, PE1D::valNotGiven());
  TEUCHOS_ASSERT_EQUALITY(point_k.Dphi, PE1D::valNotGiven());
  TEUCHOS_ASSERT(!is_null(point_kp1));
  TEUCHOS_ASSERT_INEQUALITY(point_kp1->alpha, >, ST::zero());
  TEUCHOS_ASSERT_INEQUALITY(point_kp1->phi, !=, PE1D::valNotGiven());
  TEUCHOS_ASSERT_EQUALITY(point_kp1->Dphi, PE1D::valNotGiven());
#endif

  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  bracket_.setOStream(out);
  brentsMin_.setOStream(out);

  *out << "\nStarting bracketing and brents 1D minimization linesearch ...\n";

  OSTab tab(out);

  int totalNumIters = 0;

  PointEval1D<Scalar> p_l = point_k;
  PointEval1D<Scalar> &p_m = *point_kp1; // Update in place!
  PointEval1D<Scalar> p_u;

  bool success = true;

  // A) Bracket the minimum

  int numBracketIters = -1;

  const bool bracketSuccess = bracket_.bracketMinimum(
    phi, inOutArg(p_l), inOutArg(p_m), outArg(p_u), outArg(numBracketIters) );

  if (!bracketSuccess) success = false;

  totalNumIters += numBracketIters;

  // B) Do approximate mimimization in the bracket

  if (bracketSuccess) {

    int numBrentsIters = -1;

    const bool brentsSuccess = brentsMin_.approxMinimize(
      phi, p_l, inOutArg(p_m), p_u, outArg(numBrentsIters) );

    if (!brentsSuccess) success = false;

    totalNumIters += numBrentsIters;

  }

  // C) Overall success?

  if (!is_null(numIters))
    *numIters = totalNumIters;

  return success;

}


} // namespace GlobiPack


#endif // GLOBIPACK_BRENTS_LINE_SEARCH_DEF_HPP
