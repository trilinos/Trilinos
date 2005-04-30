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

/** \brief Silly little example power method abstract numerical algorithm (ANA).
 *
 * This little function is just a silly little ANA that implements the
 * power method for estimating the dominate eigienvalue for a matrix
 * using the \ref Thyra_Op_Vec_foundational_interfaces_sec "foundational Thyra operator/vector interfaces".
 *
 * This function is small and is ment to be looked at so study its
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
  typedef Teuchos::ScalarTraits<Scalar>   ST;         // We need to use ScalarTraits to support arbitrary types!         
  typedef typename ST::magnitudeType      ScalarMag;  // This is the type for a norm
  using Teuchos::arrayArg;
  if(out) *out << "\nStarting power method ...\n\n";
  // Create workspace vectors
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    q = createMember(A.domain()),
    z = createMember(A.range()),
    r = createMember(A.range());
  // Randomize initial z
  Thyra::seed_randomize<Scalar>(0); // Make repeated runs with same data unique
  Thyra::randomize( Scalar(-ST::one()), Scalar(+ST::one()), &*z );
  // Perform iterations
  for( int iter = 0; iter < maxNumIters; ++iter ) {
    const ScalarMag z_nrm = Thyra::norm(*z);          // Compute natural norm of z
    Thyra::V_StV( &*q, Scalar(ST::one()/z_nrm), *z ); // q = (1/||z}*z 
    A.apply( Thyra::NOTRANS , *q, &*z );              // z = A*q
    *lambda = A.range()->scalarProd(*q,*z);           // lambda = <q,z> : Approximate maximum absolute eigenvalue
    if( iter%(maxNumIters/10) == 0 || iter+1 == maxNumIters ) {
      Thyra::linear_combination(                      // r = z - lambda*q : Compute residual of eigenvalue equation
        2,arrayArg<Scalar>(ST::one(),-*lambda)()
        ,arrayArg<const Thyra::VectorBase<Scalar>*>(&*z,&*q)()
        ,ST::zero(),&*r
        );
      const ScalarMag r_nrm = Thyra::norm(*r);        // Compute natural norm of r
      if(out) *out << "Iter = " << iter << ", lambda = " << (*lambda) << ", ||A*q-lambda*q|| = " << r_nrm << std::endl;
      if( r_nrm < tolerance ) return true;            // Success!
    }
  }
  return false; // Failure
} // end sillyPowerMethod

#endif // THYRA_SILLY_POWER_METHOD_HPP
