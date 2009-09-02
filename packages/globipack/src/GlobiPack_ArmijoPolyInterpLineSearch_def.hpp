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

#ifndef GLOBIPACK_POLY_INTERP_LINE_SEARCH_DEF_HPP
#define GLOBIPACK_POLY_INTERP_LINE_SEARCH_DEF_HPP


#include "GlobiPack_ArmijoPolyInterpLineSearch_decl.hpp"
#include "Teuchos_TabularOutputter.hpp"


namespace GlobiPack {


// Constructor/Initializers/Accessors


template<typename Scalar>
ArmijoPolyInterpLineSearch<Scalar>::ArmijoPolyInterpLineSearch()
  : eta_(ArmijoPolyInterpLineSearchUtils::eta_default),
    minFrac_(ArmijoPolyInterpLineSearchUtils::minFrac_default),
    maxFrac_(ArmijoPolyInterpLineSearchUtils::maxFrac_default),
    minIters_(ArmijoPolyInterpLineSearchUtils::minIters_default),
    maxIters_(ArmijoPolyInterpLineSearchUtils::maxIters_default),
    doMaxIters_(ArmijoPolyInterpLineSearchUtils::doMaxIters_default)
{}


template<typename Scalar>
Scalar ArmijoPolyInterpLineSearch<Scalar>::eta() const
{
  return eta_;
}


template<typename Scalar>
Scalar ArmijoPolyInterpLineSearch<Scalar>::minFrac() const
{
  return minFrac_;
}


template<typename Scalar>
Scalar ArmijoPolyInterpLineSearch<Scalar>::maxFrac() const
{
  return maxFrac_;
}


template<typename Scalar>
int	ArmijoPolyInterpLineSearch<Scalar>::minIters() const
{
  return minIters_;
}


template<typename Scalar>
int	ArmijoPolyInterpLineSearch<Scalar>::maxIters() const
{
  return maxIters_;
}


template<typename Scalar>
bool ArmijoPolyInterpLineSearch<Scalar>::doMaxIters() const
{
  return doMaxIters_;
}


// Overridden from ParameterListAcceptor


template<class Scalar>
void ArmijoPolyInterpLineSearch<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  typedef ScalarTraits<Scalar> ST;
  namespace AQLSU = ArmijoPolyInterpLineSearchUtils;
  using Teuchos::getParameter;
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  eta_ = getParameter<double>(*paramList, AQLSU::eta_name);
  minFrac_ = getParameter<double>(*paramList, AQLSU::minFrac_name);
  maxFrac_ = getParameter<double>(*paramList, AQLSU::maxFrac_name);
  minIters_ = getParameter<int>(*paramList, AQLSU::minIters_name);
  maxIters_ = getParameter<int>(*paramList, AQLSU::maxIters_name);
  doMaxIters_ = getParameter<bool>(*paramList, AQLSU::doMaxIters_name);
  TEUCHOS_ASSERT_INEQUALITY( eta_, >=, ST::zero() );
  TEUCHOS_ASSERT_INEQUALITY( eta_, <, ST::one() );
  TEUCHOS_ASSERT_INEQUALITY( minFrac_, >=, ST::zero() );
  TEUCHOS_ASSERT_INEQUALITY( minFrac_, <, maxFrac_ );
  TEUCHOS_ASSERT_INEQUALITY( minIters_, >=, 0 );
  TEUCHOS_ASSERT_INEQUALITY( minIters_, <=, maxIters_ );
  setMyParamList(paramList);
}


template<class Scalar>
RCP<const ParameterList>
ArmijoPolyInterpLineSearch<Scalar>::getValidParameters() const
{
  namespace AQLSU = ArmijoPolyInterpLineSearchUtils;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList());
    pl->set( AQLSU::eta_name, AQLSU::eta_default );
    pl->set( AQLSU::minFrac_name, AQLSU::minFrac_default );
    pl->set( AQLSU::maxFrac_name, AQLSU::maxFrac_default );
    pl->set( AQLSU::minIters_name, AQLSU::minIters_default );
    pl->set( AQLSU::maxIters_name, AQLSU::maxIters_default );
    pl->set( AQLSU::doMaxIters_name, AQLSU::doMaxIters_default );
    validPL = pl;
  }
  return validPL;
}


// Overrridden from LineSearchBase


template<typename Scalar>
bool ArmijoPolyInterpLineSearch<Scalar>::requiresBaseDeriv() const
{
  return true;
}


template<typename Scalar>
bool ArmijoPolyInterpLineSearch<Scalar>::requiresDerivEvals() const
{
  return false;
}


template<typename Scalar>
bool ArmijoPolyInterpLineSearch<Scalar>::doLineSearch(
  const MeritFunc1DBase<Scalar> &phi,
  const PointEval1D<Scalar> &point_k,
  const Ptr<PointEval1D<Scalar> > &point_kp1,
  const Ptr<int> &numIters_out
  ) const
{

  using Teuchos::as;
  using Teuchos::ptrFromRef;
  using Teuchos::TabularOutputter;
  typedef Teuchos::TabularOutputter TO;
  typedef ScalarTraits<Scalar> ST;
  typedef PointEval1D<Scalar> PE1D;
  using std::min;
  using std::max;

  const RCP<Teuchos::FancyOStream> out = this->getOStream();

  *out << "\nStarting Armijo Quadratic interpolation linesearch ...\n";

#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(point_k.alpha, ST::zero());
  TEUCHOS_ASSERT_INEQUALITY(point_k.phi, !=, PE1D::valNotGiven());
  TEUCHOS_ASSERT_INEQUALITY(point_k.Dphi, !=, PE1D::valNotGiven());
  TEUCHOS_ASSERT(!is_null(point_kp1));
  TEUCHOS_ASSERT_INEQUALITY(point_kp1->alpha, >, ST::zero());
  TEUCHOS_ASSERT_INEQUALITY(point_kp1->phi, !=, PE1D::valNotGiven());
  TEUCHOS_ASSERT_EQUALITY(point_kp1->Dphi, PE1D::valNotGiven());
#endif

  const Scalar phi_k = point_k.phi;
  Scalar &alpha_k = point_kp1->alpha;
  Scalar &phi_kp1 = point_kp1->phi;

  // Loop initialization (technically the first iteration)

  const Scalar Dphi_k = point_k.Dphi;

  // output header

  *out
    << "\nDphi_k = " << Dphi_k
    << "\nphi_k = " << phi_k
    << "\n";
  if (minIters_ > 0)
    *out << "\nminIters == " << minIters_ << "\n";
  if (doMaxIters_)
    *out << "\ndoMaxIters == true, maxing out maxIters = "
         <<maxIters_<<" iterations!\n"; 
  *out << "\n";
  
  TabularOutputter tblout(out);
  
  tblout.pushFieldSpec("itr", TO::INT);
  tblout.pushFieldSpec("alpha_k", TO::DOUBLE);
  tblout.pushFieldSpec("phi_kp1", TO::DOUBLE);
  tblout.pushFieldSpec("phi_kp1-frac_phi", TO::DOUBLE);
  tblout.pushFieldSpec("alpha_interp", TO::DOUBLE);
  tblout.pushFieldSpec("alpha_min", TO::DOUBLE);
  tblout.pushFieldSpec("alpha_max", TO::DOUBLE);

  tblout.outputHeader();

  // Check that this is a decent direction
  TEST_FOR_EXCEPTION( !(Dphi_k < ST::zero()), Exceptions::NotDescentDirection,
    "ArmijoPolyInterpLineSearch::doLineSearch(): "
    "The given descent direction for the given "
    "phi Dphi_k="<<Dphi_k<<" >= 0!" );
  
  // keep memory of the best value
  Scalar best_alpha = alpha_k;
  Scalar best_phi = phi_kp1;

  // Perform linesearch.
  bool success = false;
  int numIters = 0;
  for ( ; numIters < maxIters_; ++numIters ) {

    // Print out this iteration.

    Scalar frac_phi = phi_k + eta_ * alpha_k * Dphi_k;
    tblout.outputField(numIters);
    tblout.outputField(alpha_k);
    tblout.outputField(phi_kp1);
    tblout.outputField(((phi_kp1)-frac_phi));
    
    if (ST::isnaninf(phi_kp1)) {

      // Cut back the step to minFrac * alpha_k
      alpha_k = minFrac_ * alpha_k;
      best_alpha = ST::zero();
      best_phi = phi_k;
    }
    else {		

      // Armijo condition
      if (phi_kp1 < frac_phi) {
        // We have found an acceptable point
        success = true;
        if (numIters < minIters_) {
          // Keep going!
        }
        else if ( !doMaxIters_ || ( doMaxIters_ && numIters == maxIters_ - 1 ) ) {
          tblout.nextRow(true);
          break;	// get out of the loop, we are done!
        }
      }
      else {
        success = false;
      }

      // Select a new alpha to try:
      //   alpha_k = ( minFracalpha_k <= quadratic interpolation <= maxFracalpha_k )

      // Quadratic interpolation of alpha_k that minimizes phi.
      // We know the values of phi at the initail point and alpha_k and
      // the derivate of phi w.r.t alpha at the initial point and
      // that's enough information for a quandratic interpolation.
      
      Scalar alpha_quad =
        ( -as<Scalar>(0.5) * Dphi_k * alpha_k * alpha_k )
        /
        ( (phi_kp1) - phi_k - alpha_k * Dphi_k );

      tblout.outputField(alpha_quad);

      // 2009/01/29: rabartl: ToDo: Call the interpolation and then add
      // support for other types of interpolations based on various points.

      const Scalar alpha_min = minFrac_ * alpha_k;
      const Scalar alpha_max = maxFrac_ * alpha_k;

      tblout.outputField(alpha_min);
      tblout.outputField(alpha_max);

      alpha_k =
        min(
          max(alpha_min, alpha_quad),
          alpha_max
          );

    }

    tblout.nextRow(true);

    
    // Evaluate the point

    phi_kp1 = computeValue<Scalar>(phi, alpha_k);

    // Save the best point found
    if (phi_kp1 < best_phi) {
      best_phi = phi_kp1;
      best_alpha = alpha_k;
    }

  }

  if (!is_null(numIters_out))
    *numIters_out = numIters;

  if( success ) {
    *out << "\nLine search success!\n";
    return true;
  }

  // Line search failure.  Return the best point found and let the client
  // decide what to do.
  alpha_k = best_alpha;
  phi_kp1 = computeValue<Scalar>(phi, best_alpha);
  *out << "\nLine search failure!\n";
  return false; 
  
}


} // namespace GlobiPack


#endif // GLOBIPACK_POLY_INTERP_LINE_SEARCH_DEF_HPP
