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

#ifndef GLOBIPACK_GOLDEN_BRACKET_QUAD_INTERP_DEF_HPP
#define GLOBIPACK_GOLDEN_BRACKET_QUAD_INTERP_DEF_HPP


#include "GlobiPack_GoldenQuadInterpBracket_decl.hpp"
#include "Teuchos_TabularOutputter.hpp"


namespace GlobiPack {


// Constructor/Initializers/Accessors


template<typename Scalar>
GoldenQuadInterpBracket<Scalar>::GoldenQuadInterpBracket()
{}


// Overridden from ParameterListAcceptor (simple forwarding functions)


template<typename Scalar>
void GoldenQuadInterpBracket<Scalar>::setParameterList(RCP<ParameterList> const& paramList)
{
  //typedef ScalarTraits<Scalar> ST; // unused
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  // ToDo: Add parameters!
  setMyParamList(paramList);
}


template<typename Scalar>
RCP<const ParameterList> GoldenQuadInterpBracket<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList());
    // ToDo: Add parameters!
    validPL = pl;
  }
  return validPL;
}


// Bracket


template<typename Scalar>
bool GoldenQuadInterpBracket<Scalar>::bracketMinimum(
  const MeritFunc1DBase<Scalar> &phi,
  const Ptr<PointEval1D<Scalar> > &pointLower,
  const Ptr<PointEval1D<Scalar> > &pointMiddle,
  const Ptr<PointEval1D<Scalar> > &pointUpper,
  const Ptr<int> &numIters
  ) const
{

  using Teuchos::as;
  using Teuchos::TabularOutputter;
  typedef Teuchos::TabularOutputter TO;
  typedef ScalarTraits<Scalar> ST;
  using Teuchos::OSTab;
#ifdef TEUCHOS_DEBUG
  typedef PointEval1D<Scalar> PE1D;

  TEUCHOS_TEST_FOR_EXCEPT(is_null(pointLower));
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pointUpper));
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pointMiddle));
  TEUCHOS_ASSERT_INEQUALITY(pointLower->alpha, <, pointMiddle->alpha);
  TEUCHOS_ASSERT_INEQUALITY(pointLower->phi, !=, PE1D::valNotGiven());
  TEUCHOS_ASSERT_INEQUALITY(pointMiddle->phi, !=, PE1D::valNotGiven());
#endif

  const RCP<Teuchos::FancyOStream> out = this->getOStream();

  // ToDo: Make these variable!
  const Scalar GOLDEN_RATIO = 1.618033988749895;
  const Scalar SMALL_DIV = 1e-20;
  const Scalar MAX_EXTRAP_FACTOR = 100.0;
  const int MAX_TOTAL_ITERS = 30;

  *out << "\nStarting golden quadratic interpolating bracketing of the minimum ...\n\n";

  // Repeatedly evaluate the function along the search direction until
  // we know we've bracketed a minimum.

  Scalar &alpha_l = pointLower->alpha, &phi_l = pointLower->phi;
  Scalar &alpha_m = pointMiddle->alpha, &phi_m = pointMiddle->phi;
  Scalar &alpha_u = pointUpper->alpha = ST::nan(), &phi_u = pointUpper->phi = ST::nan();

  Scalar tmp = ST::nan(), q = ST::nan(), r = ST::nan();

  const Scalar zero = ST::zero();

  // This does a simple backtracking
  alpha_u = zero;
  const Scalar goldinv = 1.0/(1.0+GOLDEN_RATIO);

  TabularOutputter tblout(out);

  tblout.pushFieldSpec("itr", TO::INT);
  tblout.pushFieldSpec("alpha_l", TO::DOUBLE);
  tblout.pushFieldSpec("alpha_m", TO::DOUBLE);
  tblout.pushFieldSpec("alpha_u", TO::DOUBLE);
  tblout.pushFieldSpec("phi_l", TO::DOUBLE);
  tblout.pushFieldSpec("phi_m", TO::DOUBLE);
  tblout.pushFieldSpec("phi_u", TO::DOUBLE);
  tblout.pushFieldSpec("step type             ", TO::STRING);

  tblout.outputHeader();

  int icount = 0;

  std::string stepType = "";

  //
  // A) Find phi_l > phi_m first
  //

  tblout.outputField("-");
  tblout.outputField(alpha_l);
  tblout.outputField(alpha_m);
  tblout.outputField("-");
  tblout.outputField(phi_l);
  tblout.outputField(phi_m);
  tblout.outputField("-");
  tblout.outputField("start");
  tblout.nextRow();

  for (; icount < MAX_TOTAL_ITERS; ++icount) {

    // ToDo: Put in a check for NAN and backtrack if you find it!

    if (phi_l > phi_m) {
      break;
    }

    stepType = "golden back";
    alpha_u = alpha_m;
    phi_u = phi_m;
    alpha_m = goldinv * (alpha_u + GOLDEN_RATIO*alpha_l);
    phi_m = computeValue<Scalar>(phi, alpha_m);

    tblout.outputField(icount);
    tblout.outputField(alpha_l);
    tblout.outputField(alpha_m);
    tblout.outputField(alpha_u);
    tblout.outputField(phi_l);
    tblout.outputField(phi_m);
    tblout.outputField(phi_u);
    tblout.outputField(stepType);
    tblout.nextRow();

  }

  if (alpha_u == zero) {
    // The following factor of gold was reduced to (GOLDEN_RATIO-1) to save
    // one function evaluation near convergence.
    alpha_u = alpha_m + (GOLDEN_RATIO-1.0) * (alpha_m-alpha_l);
    phi_u = computeValue<Scalar>(phi, alpha_u);
  }

  //
  // B) Quadratic interpolation iterations
  //

  bool bracketedMin = false;

  for (; icount < MAX_TOTAL_ITERS; ++icount) {

    if (phi_m < phi_u) {
      bracketedMin = true;
      break;
    }

    // find the extremum alpha_quad of a quadratic model interpolating there
    // points
    q = (phi_m-phi_l)*(alpha_m-alpha_u);
    r = (phi_m-phi_u)*(alpha_m-alpha_l);

    // avoid division by small (q-r) by bounding with signed minimum
    tmp = ST::magnitude(q-r);
    tmp = (tmp > SMALL_DIV ? tmp : SMALL_DIV);
    tmp = (q-r >= 0  ? tmp : -tmp);

    Scalar alpha_quad =
      alpha_m - (q*(alpha_m-alpha_u) - r*(alpha_m-alpha_l))/(2.0*tmp);

    // maximum point for which we trust the interpolation
    const Scalar alpha_lim = alpha_m + MAX_EXTRAP_FACTOR * (alpha_u-alpha_m);

    // now detect which interval alpha_quad is in and act accordingly
    bool skipToNextIter = false;
    Scalar phi_quad = ST::nan();
    if ( (alpha_m-alpha_quad)*(alpha_quad-alpha_u) > zero ) {  // [alpha_m, alpha_u]
      phi_quad = computeValue<Scalar>(phi, alpha_quad);
      if (phi_quad < phi_u) {                        // use points [b, alpha_quad, c]
        alpha_l = alpha_m;
        phi_l = phi_m;
        alpha_m = alpha_quad;
        phi_m = phi_quad;
        skipToNextIter = true;
        stepType = "alpha_quad middle";
      }
      else if (phi_quad > phi_m) {                   // use points [a, b, alpha_quad]
        alpha_u = alpha_quad;
        phi_u = phi_quad;
        skipToNextIter = true;
        stepType = "alpha_quad upper";
      }
      else {
        alpha_quad = alpha_u + GOLDEN_RATIO*(alpha_u-alpha_m);
        phi_quad = computeValue<Scalar>(phi, alpha_quad);
      }
    }

    if (!skipToNextIter) {

      if ((alpha_u-alpha_quad)*(alpha_quad-alpha_lim) > zero) {  // [alpha_u, alpha_lim]
        phi_quad = computeValue<Scalar>(phi, alpha_quad);
        stepType = "[alpha_u, alpha_lim]";
        if (phi_quad < phi_u) {
          alpha_m = alpha_u;
          alpha_u = alpha_quad;
          alpha_quad = alpha_u + GOLDEN_RATIO*(alpha_u-alpha_m);
          phi_m = phi_u;
          phi_u = phi_quad;
          phi_quad = computeValue<Scalar>(phi, alpha_quad);
          stepType = "phi_quad < phi_u";
        }
      }
      else if ((alpha_quad-alpha_lim)*(alpha_lim-alpha_u) >= zero ) { // [alpha_lim, inf]
        alpha_quad = alpha_lim;
        phi_quad = computeValue<Scalar>(phi, alpha_quad);
        stepType = "[alpha_lim, inf]";
      }
      else {                                  // [0,alpha_m]
        alpha_quad = alpha_u + GOLDEN_RATIO*(alpha_u-alpha_m);
        phi_quad = computeValue<Scalar>(phi, alpha_quad);
        stepType = "[0, alpha_m]";
      }

      // shift to newest 3 points before loop
      alpha_l = alpha_m;
      phi_l = phi_m;
      alpha_m = alpha_u;
      phi_m = phi_u;
      alpha_u = alpha_quad;
      phi_u = phi_quad;

    }

    tblout.outputField(icount);
    tblout.outputField(alpha_l);
    tblout.outputField(alpha_m);
    tblout.outputField(alpha_u);
    tblout.outputField(phi_l);
    tblout.outputField(phi_m);
    tblout.outputField(phi_u);
    tblout.outputField(stepType);
    tblout.nextRow();

  }  // end for loop

  if (icount >= MAX_TOTAL_ITERS) {
    *out <<"\nExceeded maximum number of iterations!.\n";
  }

  if (!is_null(numIters)) {
    *numIters = icount;
  }

  *out << "\n";

  return bracketedMin;

}


} // namespace GlobiPack


#endif // GLOBIPACK_GOLDEN_BRACKET_QUAD_INTERP_DEF_HPP
