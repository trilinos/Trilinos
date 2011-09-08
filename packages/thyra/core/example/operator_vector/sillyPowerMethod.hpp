// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
  const Thyra::LinearOpBase<Scalar> &A,
  const int maxNumIters,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance,
  const Teuchos::Ptr<Scalar> &lambda,
  std::ostream &out
  )
{

  // Create some typedefs and some other stuff to make the code cleaner
  typedef Teuchos::ScalarTraits<Scalar> ST; typedef typename ST::magnitudeType ScalarMag;
  using Thyra::apply;
  const Scalar one = ST::one(); using Thyra::NOTRANS;
  typedef Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > VectorSpacePtr;
  typedef Teuchos::RCP<Thyra::VectorBase<Scalar> > VectorPtr;

  // Initialize
  out << "\nStarting power method (target tolerance = "<<tolerance<<") ...\n\n";
  VectorPtr q = createMember(A.domain()), z = createMember(A.range()), r = createMember(A.range());
  Thyra::seed_randomize<Scalar>(0);
  Thyra::randomize( Scalar(-one), Scalar(+one), z.ptr() );

  // Perform iterations
  for( int iter = 0; iter < maxNumIters; ++iter ) {
    const ScalarMag z_nrm = norm(*z);           // Compute natural norm of z
    V_StV( q.ptr(), Scalar(one/z_nrm), *z );    // q = (1/||z||)*z 
    apply<Scalar>( A, NOTRANS , *q, z.ptr() );  // z = A*q
    *lambda = scalarProd(*q,*z);                // lambda = <q,z>
    if( iter%(maxNumIters/10) == 0 || iter+1 == maxNumIters ) {
      V_StVpV(r.ptr(),Scalar(-*lambda),*q,*z);  // r = -lambda*q + z
      const ScalarMag r_nrm = norm(*r);         // Compute natural norm of r
      out << "Iter = " << iter << ", lambda = " << (*lambda)
          << ", ||A*q-lambda*q|| = " << r_nrm << std::endl;
      if( r_nrm < tolerance )
        return true;  // Success!
    }
  }

  out << "\nMaximum number of iterations exceeded with ||-lambda*q + z||"
    " > tolerence = " << tolerance << "\n";
  return false; // Failure

} // end sillyPowerMethod


#endif // THYRA_SILLY_POWER_METHOD_HPP
