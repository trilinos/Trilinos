// ***********************************************************************
// 
//                           Stokhos Package
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
// Questions? Contact Christopher W. Miller (cmiller@math.umd.edu).
// 
// ***********************************************************************
// @HEADER

template <typename ordinal_type, typename value_type>
Stokhos::DiscretizedStieltjesBasis<ordinal_type,value_type>::
DiscretizedStieltjesBasis(const std::string& label,
			  const ordinal_type& p, 
			  value_type (*weightFn)(const value_type&),
			  const value_type& leftEndPt,
			  const value_type& rightEndPt, 
			  bool normalize, 
			  Stokhos::GrowthPolicy growth) :
  RecurrenceBasis<ordinal_type,value_type>(std::string("DiscretizedStieltjes -- ") + label, p, normalize, growth),
  scaleFactor(1),
  leftEndPt_(leftEndPt),
  rightEndPt_(rightEndPt),
  weightFn_(weightFn)
{   
  // Set up quadrature points for discretized stieltjes procedure
  Teuchos::RCP<const Stokhos::LegendreBasis<ordinal_type,value_type> > quadBasis = 
    Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(this->p));
  quadBasis->getQuadPoints(200*this->p, quad_points, quad_weights, quad_values);

  // Setup rest of recurrence basis
  this->setup();
}

template <typename ordinal_type, typename value_type>
Stokhos::DiscretizedStieltjesBasis<ordinal_type,value_type>::
DiscretizedStieltjesBasis(const ordinal_type& p, 
			  const DiscretizedStieltjesBasis& basis) :
  RecurrenceBasis<ordinal_type,value_type>(p, basis),
  scaleFactor(basis.scaleFactor),
  leftEndPt_(basis.leftEndPt_),
  rightEndPt_(basis.rightEndPt_),
  weightFn_(basis.weightFn_)
{   
  // Set up quadrature points for discretized stieltjes procedure
  Teuchos::RCP<const Stokhos::LegendreBasis<ordinal_type,value_type> > quadBasis = 
    Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(this->p));
  quadBasis->getQuadPoints(200*this->p, quad_points, quad_weights, quad_values);

  // Compute coefficients in 3-term recurrsion
  computeRecurrenceCoefficients(p+1, this->alpha, this->beta, this->delta,
				this->gamma);

  // Setup rest of recurrence basis
  this->setup();
}

template <typename ordinal_type, typename value_type>
Stokhos::DiscretizedStieltjesBasis<ordinal_type,value_type>::
~DiscretizedStieltjesBasis()
{
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::DiscretizedStieltjesBasis<ordinal_type,value_type>::
computeRecurrenceCoefficients(ordinal_type k,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta,
			      Teuchos::Array<value_type>& gamma) const
{
  //The Discretized Stieltjes polynomials are defined by a recurrance relation,
  //P_n+1 = \gamma_n+1[(x-\alpha_n) P_n - \beta_n P_n-1].
  //The alpha and beta coefficients are generated first using the
  //discritized stilges procidure described in "On the Calculation of DiscretizedStieltjes Polynomials and Quadratures",
  //Robin P. Sagar, Vedene H. Smith.  The gamma coefficients are then optionally set so that each
  //polynomial has norm 1.  If normalization is not enabled then the gammas are set to 1.

  scaleFactor = 1;
   
  //First renormalize the weight function so that it has measure 1.
  value_type oneNorm = expectedValue_J_nsquared(0, alpha, beta);
  //future evaluations of the weight function will scale it by this factor.
  scaleFactor = 1/oneNorm;
  
  
  value_type integral2;
  //NOTE!! This evaluation of 'expectedValue_J_nsquared(0)' is different
  //from the one above since we rescaled the weight.  Don't combine
  //the two!!!

  value_type past_integral = expectedValue_J_nsquared(0, alpha, beta);
  alpha[0] = expectedValue_tJ_nsquared(0, alpha, beta)/past_integral;
  //beta[0] := \int_-c^c w(x) dx.
  beta[0] = 1;
  delta[0] = 1;
  gamma[0] = 1;
  //These formulas are from the above reference.
  for (ordinal_type n = 1; n<k; n++){
    integral2 = expectedValue_J_nsquared(n, alpha, beta);
    alpha[n] = expectedValue_tJ_nsquared(n, alpha, beta)/integral2;
    beta[n] = integral2/past_integral;
    past_integral = integral2;
    delta[n] = 1.0;
    gamma[n] = 1.0;
  }

  return false;
}

template <typename ordinal_type, typename value_type>
value_type 
Stokhos::DiscretizedStieltjesBasis<ordinal_type,value_type>::
evaluateWeight(const value_type& x) const
{
  return (x < leftEndPt_ || x > rightEndPt_) ? 0: scaleFactor*weightFn_(x);
} 

template <typename ordinal_type, typename value_type>
value_type
Stokhos::DiscretizedStieltjesBasis<ordinal_type,value_type>::
expectedValue_tJ_nsquared(const ordinal_type& order, 
			  const Teuchos::Array<value_type>& alpha,
			  const Teuchos::Array<value_type>& beta) const
{
  //Impliments a gaussian quadrature routine to evaluate the integral,
  // \int_-c^c J_n(x)^2w(x)dx.  This is needed to compute the recurrance coefficients.
  value_type integral = 0;
  for(ordinal_type quadIdx = 0; 
      quadIdx < static_cast<ordinal_type>(quad_points.size()); quadIdx++) {
    value_type x = (rightEndPt_ - leftEndPt_)*.5*quad_points[quadIdx] + 
      (rightEndPt_ + leftEndPt_)*.5;
    value_type val = evaluateRecurrence(x,order,alpha,beta);
    integral += x*val*val*evaluateWeight(x)*quad_weights[quadIdx];
  }
  
  return integral*(rightEndPt_ - leftEndPt_);
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::DiscretizedStieltjesBasis<ordinal_type,value_type>::
expectedValue_J_nsquared(const ordinal_type& order, 
			 const Teuchos::Array<value_type>& alpha,
			 const Teuchos::Array<value_type>& beta) const
{
  //Impliments a gaussian quadrature routineroutine to evaluate the integral,
  // \int_-c^c J_n(x)^2w(x)dx.  This is needed to compute the recurrance coefficients.
  value_type integral = 0;
  for(ordinal_type quadIdx = 0; 
      quadIdx < static_cast<ordinal_type>(quad_points.size()); quadIdx++){
    value_type x = (rightEndPt_ - leftEndPt_)*.5*quad_points[quadIdx] + 
      (rightEndPt_ + leftEndPt_)*.5;
    value_type val = evaluateRecurrence(x,order,alpha,beta);
    integral += val*val*evaluateWeight(x)*quad_weights[quadIdx];
  }
  
  return integral*(rightEndPt_ - leftEndPt_);
}

template <typename ordinal_type, typename value_type>
value_type 
Stokhos::DiscretizedStieltjesBasis<ordinal_type, value_type>::
eval_inner_product(const ordinal_type& order1, const ordinal_type& order2) const
{
   //Impliments a gaussian quadrature routine to evaluate the integral,
  // \int_-c^c J_n(x)J_m w(x)dx.  This method is intended to allow the user to
  // test for orthogonality and proper normalization.
  value_type integral = 0;
  for(ordinal_type quadIdx = 0; 
      quadIdx < static_cast<ordinal_type>(quad_points.size()); quadIdx++){
    value_type x = (rightEndPt_ - leftEndPt_)*.5*quad_points[quadIdx] + 
      (rightEndPt_ + leftEndPt_)*.5;
    integral += this->evaluate(x,order1)*this->evaluate(x,order2)*evaluateWeight(x)*quad_weights[quadIdx];
  }
  
  return integral*(rightEndPt_ - leftEndPt_);
}

template <typename ordinal_type, typename value_type>
value_type 
Stokhos::DiscretizedStieltjesBasis<ordinal_type, value_type>::
evaluateRecurrence(const value_type& x, 
		   ordinal_type k, 
		   const Teuchos::Array<value_type>& alpha,
		   const Teuchos::Array<value_type>& beta) const
{
  if (k == 0)
    return value_type(1.0);
  else if (k == 1)
    return x-alpha[0];

  value_type v0 = value_type(1.0);
  value_type v1 = x-alpha[0]*v0;
  value_type v2 = value_type(0.0);
  for (ordinal_type i=2; i<=k; i++) {
    v2 = (x-alpha[i-1])*v1 - beta[i-1]*v0;
    v0 = v1;
    v1 = v2;
  }

  return v2;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::DiscretizedStieltjesBasis<ordinal_type, value_type>::
cloneWithOrder(ordinal_type p) const
{
   return Teuchos::rcp(new Stokhos::DiscretizedStieltjesBasis<ordinal_type,value_type>(p,*this));
}
