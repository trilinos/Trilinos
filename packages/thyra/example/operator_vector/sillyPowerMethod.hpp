// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SILLY_POWER_METHOD_HPP
#define THYRA_SILLY_POWER_METHOD_HPP

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp" // DEBUG STUFF

/** \brief Silly little example power method abstract numerical algorithm (ANA).
 *
 * This little function is just a silly little ANA that implements the
 * power method for estimating the dominate eigenvalue for a matrix
 * using the \ref Thyra_Op_Vec_foundational_interfaces_sec "foundational Thyra operator/vector interfaces".
 *
 * This function is small and is meant to be looked at so study its
 * implementation by clicking on the below link to its definition.
 *
 * \ingroup Thyra_Op_Vec_examples_power_method_grp
 */
template<class Scalar>
bool sillyPowerMethod(
  const Thyra::LinearOpBase<Scalar>                              &A
  ,const int                                                     maxNumIters
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType   tolerance
  ,Scalar                                                        *lambda
  ,std::ostream                                                  *out          = NULL
  )
{
  if(out) *out << std::setprecision(16); // DEBUG STUFF
  // Create some typedefs and some other stuff to make the code cleaner
  typedef Teuchos::ScalarTraits<Scalar> ST; typedef typename ST::magnitudeType ScalarMag;
  const Scalar one = ST::one(); using Thyra::NOTRANS;
  typedef Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > VectorSpacePtr;
  typedef Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > VectorPtr;
  // Initialize
  if(out) *out << "\nStarting power method ...\n\n";
  VectorPtr q = createMember(A.domain()), z = createMember(A.range()), r = createMember(A.range());
  Thyra::seed_randomize<Scalar>(0);
  Thyra::randomize( Scalar(-one), Scalar(+one), &*z );
  if(out) *out << "z=\n" << *z; // DEBUG STUFF
  // Perform iterations
  for( int iter = 0; iter < maxNumIters; ++iter ) {
    if(out) *out << "itrer="<<iter<<"\n"; // DEBUG STUFF
    const ScalarMag z_nrm = norm(*z);       // Compute natural norm of z
    if(out) *out << "||z||="<<z_nrm<<"\n"; // DEBUG STUFF
    V_StV( &*q, Scalar(one/z_nrm), *z );    // q = (1/||z}*z 
    if(out) *out << "q=\n" << *q; // DEBUG STUFF
    apply( A, NOTRANS , *q, &*z );          // z = A*q
    if(out) *out << "z=\n" << *z; // DEBUG STUFF
    *lambda = scalarProd(*q,*z);            // lambda = <q,z>    : Approximate maximum absolute eigenvalue
    if(out) *out << "lambda="<<*lambda<<"\n"; // DEBUG STUFF
    if( iter%(maxNumIters/10) == 0 || iter+1 == maxNumIters ) {
      V_StVpV(&*r,Scalar(-*lambda),*q,*z);  // r = -lambda*q + z : Compute residual of eigenvalue equation
      if(out) *out << "r=\n" << *r; // DEBUG STUFF
      const ScalarMag r_nrm = norm(*r);     // Compute natural norm of r
      if(out) *out << "Iter = " << iter << ", lambda = " << (*lambda) << ", ||A*q-lambda*q|| = " << r_nrm << std::endl;
      if( r_nrm < tolerance ) return true;  // Success!
    }
  }
  return false; // Failure
} // end sillyPowerMethod

#endif // THYRA_SILLY_POWER_METHOD_HPP
