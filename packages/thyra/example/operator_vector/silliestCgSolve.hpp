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

#ifndef THYRA_SILLIEST_CG_SOLVE_HPP
#define THYRA_SILLIEST_CG_SOLVE_HPP

#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_VectorImpl.hpp"
#include "Thyra_AssertOp.hpp"

/** \brief Silly little example unpreconditioned CG solver that uses handles.
 *
 * This little function is just a silly little ANA that implements the
 * CG (conjugate gradient) method for solving symmetric positive definite
 * systems using the \ref Thyra_Op_Vec_foundational_interfaces_sec "foundational Thyra operator/vector interfaces".
 *
 * This function is small and is meant to be looked at so study its
 * implementation by clicking on the below link to its definition.
 *
 * \ingroup Thyra_Op_Vec_examples_cg_grp
 */
template<class Scalar>
bool silliestCgSolve(
  Thyra::ConstLinearOperator<Scalar>                           const& A
  ,Thyra::ConstVector<Scalar>                                  const& b
  ,const int                                                          maxNumIters
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType        tolerance
  ,Thyra::Vector<Scalar>                                            & x
  ,std::ostream                                                     * out = NULL
  )
{
  // Create some typedefs and inject some names into local namespace
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar one = ST::one(), zero = ST::zero();
  using Thyra::VectorSpace;
  using Thyra::Vector;
  // Initialization of the algorithm
  const VectorSpace<Scalar>     space = A.domain();
  Vector<Scalar>                r = b - A*x;
  ScalarMag                     r0_nrm = norm(r);
  if(r0_nrm==zero) return true;
  Vector<Scalar>                p(space), q(space);
  Scalar                        rho_old = -one;
  // Perform the iterations
  for( int iter = 0; iter <= maxNumIters; ++iter ) {
    // Check convergence and output iteration
    const ScalarMag    r_nrm = norm(r);
    const bool         isConverged = ( (r_nrm/r0_nrm) <= tolerance );
    if( ( iter%(maxNumIters/10+1) == 0 || iter == maxNumIters || isConverged ) && out )
      *out<<"Iter = "<<iter<<", ||b-A*x||/||b-A*x0|| = "<<(r_nrm/r0_nrm)<<std::endl;
    if( r_nrm/r0_nrm < tolerance ) return true; // Success!
    // Compute the iteration
    const Scalar      rho = inner(r,r);
    if(iter==0)       copyInto(r,p);
    else              p = Scalar(rho/rho_old)*p + r;
    q = A*p;
    const Scalar      alpha = rho/inner(p,q);
    x += Scalar(+alpha)*p;
    r += Scalar(-alpha)*q;
    rho_old = rho;
  }
  return false; // Failure
} // end silliestCgSolve

#endif // THYRA_SILLIEST_CG_SOLVE_HPP
