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

#ifndef GLOBIPACK_BRENTS_1D_MINIMIZATION_DEF_HPP
#define GLOBIPACK_BRENTS_1D_MINIMIZATION_DEF_HPP


#include "GlobiPack_Brents1DMinimization_decl.hpp"
#include "Teuchos_TabularOutputter.hpp"


namespace GlobiPack {


// Constructor/Initializers/Accessors


template<typename Scalar>
Brents1DMinimization<Scalar>::Brents1DMinimization()
  :rel_tol_(Brents1DMinimizationUtils::rel_tol_default),
   bracket_tol_(Brents1DMinimizationUtils::bracket_tol_default),
   max_iters_(Brents1DMinimizationUtils::max_iters_default)
{}


// Overridden from ParameterListAcceptor (simple forwarding functions)


template<typename Scalar>
void Brents1DMinimization<Scalar>::setParameterList(RCP<ParameterList> const& paramList)
{
  typedef ScalarTraits<Scalar> ST;
  namespace BMU = Brents1DMinimizationUtils;
  using Teuchos::getParameter;
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  rel_tol_ = getParameter<double>(*paramList, BMU::rel_tol_name);
  bracket_tol_ = getParameter<double>(*paramList, BMU::bracket_tol_name);
  max_iters_ = getParameter<int>(*paramList, BMU::max_iters_name);
  TEUCHOS_ASSERT_INEQUALITY( rel_tol_, >, ST::zero() );
  TEUCHOS_ASSERT_INEQUALITY( bracket_tol_, >, ST::zero() );
  TEUCHOS_ASSERT_INEQUALITY( max_iters_, >=, 0 );
  setMyParamList(paramList);
}


template<typename Scalar>
RCP<const ParameterList> Brents1DMinimization<Scalar>::getValidParameters() const
{
  namespace BMU = Brents1DMinimizationUtils;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList());
    pl->set( BMU::rel_tol_name, BMU::rel_tol_default );
    pl->set( BMU::bracket_tol_name, BMU::bracket_tol_default );
    pl->set( BMU::max_iters_name, BMU::max_iters_default );
    validPL = pl;
  }
  return validPL;
}


// Bracket


template<typename Scalar>
bool Brents1DMinimization<Scalar>::approxMinimize(
  const MeritFunc1DBase<Scalar> &phi,
  const PointEval1D<Scalar> &pointLower,
  const Ptr<PointEval1D<Scalar> > &pointMiddle,
  const PointEval1D<Scalar> &pointUpper,
  const Ptr<int> &numIters
  ) const
{

  using Teuchos::as;
  using Teuchos::TabularOutputter;
  typedef Teuchos::TabularOutputter TO;
  typedef ScalarTraits<Scalar> ST;
  using Teuchos::OSTab;
  typedef PointEval1D<Scalar> PE1D;
  using std::min;
  using std::max;
  
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(pointMiddle));
  TEUCHOS_ASSERT_INEQUALITY(pointLower.alpha, <, pointMiddle->alpha);
  TEUCHOS_ASSERT_INEQUALITY(pointMiddle->alpha, <, pointUpper.alpha);
  TEUCHOS_ASSERT_INEQUALITY(pointLower.phi, !=, PE1D::valNotGiven());
  TEUCHOS_ASSERT_INEQUALITY(pointMiddle->phi, !=, PE1D::valNotGiven());
  TEUCHOS_ASSERT_INEQUALITY(pointUpper.phi, !=, PE1D::valNotGiven());
#endif

  const RCP<Teuchos::FancyOStream> out = this->getOStream();

  *out << "\nStarting Brent's 1D minimization algorithm ...\n\n";
  
  TabularOutputter tblout(out);
  
  tblout.pushFieldSpec("itr", TO::INT);
  tblout.pushFieldSpec("alpha_a", TO::DOUBLE);
  tblout.pushFieldSpec("alpha_min", TO::DOUBLE);
  tblout.pushFieldSpec("alpha_b", TO::DOUBLE);
  tblout.pushFieldSpec("phi(alpha_min)", TO::DOUBLE);
  tblout.pushFieldSpec("alpha_b - alpha_a", TO::DOUBLE);
  tblout.pushFieldSpec("alpha_min - alpha_avg", TO::DOUBLE);
  tblout.pushFieldSpec("tol", TO::DOUBLE);

  tblout.outputHeader();
  
  const Scalar INV_GOLD2=0.3819660112501051518; // (1/golden-ratio)^2
  const Scalar TINY = ST::squareroot(ST::eps());
  
  const Scalar alpha_l = pointLower.alpha, phi_l = pointLower.phi;
  Scalar &alpha_m = pointMiddle->alpha, &phi_m = pointMiddle->phi;
  const Scalar alpha_u = pointUpper.alpha, phi_u = pointUpper.phi;

  Scalar d = ST::nan();
  Scalar e = ST::nan();
  Scalar u = ST::nan();
  
  Scalar phi_w = min(phi_l, phi_u);

  Scalar alpha_v = ST::nan();
  Scalar alpha_w = ST::nan();
  Scalar phi_v = ST::nan();

  if (phi_w == phi_l){  
    alpha_w  = alpha_l;
    alpha_v  = alpha_u;
    phi_v = phi_u;
  }
  else {
    alpha_w  = alpha_u;
    alpha_v  = alpha_l;
    phi_v = phi_l;
  }

  Scalar alpha_min = alpha_m;
  Scalar phi_min = phi_m;
  Scalar alpha_a = alpha_l;
  Scalar alpha_b = alpha_u;
  
  bool foundMin = false;

  int iteration = 0;

  for ( ; iteration <= max_iters_; ++iteration) {

    if (iteration < 2)
      e = 2.0 * (alpha_b - alpha_a);

    const Scalar alpha_avg = 0.5 *(alpha_a + alpha_b);
    const Scalar tol1 = rel_tol_ * ST::magnitude(alpha_min) + TINY;
    const Scalar tol2 = 2.0 * tol1;

    const Scalar step_diff = alpha_min - alpha_avg;
    const Scalar step_diff_tol = (tol2 + bracket_tol_ * (alpha_b - alpha_a));

    // 2009/02/11: rabartl: Above, I changed from (tol2-0.5*(alpha_b-alpha_a)) which is
    // actually in Brent's netlib code!  This gives a negative tolerence when
    // the solution alpha_min is near a minimum so you will max out the iters because
    // a possitive number can never be smaller than a negative number.  The
    // above convergence criteria makes sense to me.

    tblout.outputField(iteration);
    tblout.outputField(alpha_a);
    tblout.outputField(alpha_min);
    tblout.outputField(alpha_b);
    tblout.outputField(phi_min);
    tblout.outputField(alpha_b - alpha_a);
    tblout.outputField(step_diff);
    tblout.outputField(step_diff_tol);
    tblout.nextRow();

    // If the difference between current point and the middle of the shrinking
    // interval [alpha_a, alpha_b] is relatively small, then terminate the
    // algorithm.  Also, terminate the algorithm if this difference is small
    // relative to the size of alpha.  Does this make sense?  However, don't
    // terminate on the very first iteration because we have to take at least
    // one step.
    if (
      ST::magnitude(step_diff) <= step_diff_tol
      && iteration > 0
      )
    {
      foundMin = true;
      break;
    }
    // 2009/02/11: rabartl: Above, I added the iteration > 0 condition because
    // the original version that I was given would terminate on the first
    // first iteration if the initial guess for alpha happened to be too close
    // to the midpoint of the bracketing interval!  Is that crazy or what!

    if (ST::magnitude(e) > tol1 || iteration < 2) {

      const Scalar r = (alpha_min - alpha_w) * (phi_min - phi_v);
      Scalar q = (alpha_min - alpha_v) * (phi_min - phi_w);
      Scalar p = (alpha_min - alpha_v) * q - (alpha_min - alpha_w) * r;
      q = 2.0 * (q - r);
      if (q > ST::zero())
        p = -p;
      q = ST::magnitude(q);
      const Scalar etemp = e;
      e = d;

      if ( ST::magnitude(p) >= ST::magnitude(0.5 * q * etemp)
        || p <= q * (alpha_a - alpha_min)
        || p >= q * (alpha_b - alpha_min)
        )
      {
        e = (alpha_min >= alpha_avg ? alpha_a - alpha_min : alpha_b - alpha_min);
        d = INV_GOLD2 * e;
      }
      else {
        d = p/q;
        u = alpha_min + d;
        if (u - alpha_a < tol2 || alpha_b - u < tol2) 
           // sign(tol1,alpha_avg-alpha_min)
          d = ( alpha_avg - alpha_min > ST::zero()
            ? ST::magnitude(tol1)
            : -ST::magnitude(tol1) );
      }

    }
    else {

      e = (alpha_min >= alpha_avg ? alpha_a - alpha_min : alpha_b - alpha_min);
      d = INV_GOLD2 * e;

    }
    
    u = ( ST::magnitude(d) >= tol1
      ? alpha_min + d
      : alpha_min + (d >= 0 ? ST::magnitude(tol1) : -ST::magnitude(tol1))
      ); 
    
    const Scalar phi_eval_u = computeValue<Scalar>(phi, u);

    if (phi_eval_u <= phi_min) {

      if (u >= alpha_min)
        alpha_a = alpha_min;
      else
        alpha_b = alpha_min;

      alpha_v = alpha_w;
      phi_v = phi_w;
      alpha_w = alpha_min;
      phi_w = phi_min;
      alpha_min = u;
      phi_min = phi_eval_u;

    }
    else {

      if (u < alpha_min)
        alpha_a = u;
      else
        alpha_b = u;

      if (phi_eval_u <= phi_w || alpha_w == alpha_min) {
        alpha_v = alpha_w;
        phi_v = phi_w;
        alpha_w = u;
        phi_w = phi_eval_u;
      }
      else if (phi_eval_u <= phi_v || alpha_v == alpha_min || alpha_v == alpha_w) {
        alpha_v = u;
        phi_v = phi_eval_u;
      }

    }
  }

  alpha_m = alpha_min;
  phi_m = phi_min;
  if (!is_null(numIters))
    *numIters = iteration;
  
  if (foundMin) {
    *out <<"\nFound the minimum alpha="<<alpha_m<<", phi(alpha)="<<phi_m<<"\n";
  }
  else {
    *out <<"\nExceeded maximum number of iterations!\n";
  }

  *out << "\n";

  return foundMin;

}


} // namespace GlobiPack


#endif // GLOBIPACK_BRENTS_1D_MINIMIZATION_DEF_HPP
