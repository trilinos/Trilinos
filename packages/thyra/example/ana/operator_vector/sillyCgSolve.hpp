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

#ifndef THYRA_SILLY_CG_SOLVE_HPP
#define THYRA_SILLY_CG_SOLVE_HPP

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"

/** \brief Silly little example unpreconditioned CG solver
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
bool sillyCgSolve(
  const Thyra::LinearOpBase<Scalar>                              &A
  ,const Thyra::VectorBase<Scalar>                               &b
  ,const int                                                     maxNumIters
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType   tolerance
  ,Thyra::VectorBase<Scalar>                                     *x
  ,std::ostream                                                  *out          = NULL
  )
{
  using Teuchos::RefCountPtr;
  typedef Teuchos::ScalarTraits<Scalar>   ST;         // We need to use ScalarTraits to support arbitrary types.         
  typedef typename ST::magnitudeType      ScalarMag;  // This is the type returned from a vector norm.
  // Validate input
  TEST_FOR_EXCEPT(x==NULL);
  THYRA_ASSERT_LINEAR_OP_VEC_APPLY_SPACES("sillyCgSolve()",A,Thyra::NOTRANS,*x,&b); // A*x - b agree?
  Teuchos::EVerbosityLevel vl = Teuchos::VERB_MEDIUM; // Set the verbosity level
  if(out) *out << "\nStarting CG solver ...\n" << std::scientific << "\ndescribe A:\n"<<describe(A,vl)
               << "\ndescribe b:\n"<<describe(b,vl)<<"\ndescribe x:\n"<<describe(*x,vl)<<"\n";
  // Get the vector space (domain and range spaces should be the same)
  RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > space = A.domain();
  // Compute initial residual : r = b - A*x
  RefCountPtr<Thyra::VectorBase<Scalar> > r = createMember(space);                     
  Thyra::assign(&*r,b);                               // r = b
  A.apply(Thyra::NOTRANS,*x,&*r,-ST::one(),ST::one());// r = -A*x + r
  const ScalarMag r0_nrm = Thyra::norm(*r);           // Compute ||r0|| = sqrt(<r0,r0>) for convergence test
  if(r0_nrm == ST::zero()) return true;               // Trivial RHS and initial LHS guess?
  // Create workspace vectors and scalars
  RefCountPtr<Thyra::VectorBase<Scalar> > p = createMember(space), q = createMember(space);
  Scalar rho_old;
  // Perform the iterations
  for( int iter = 0; iter <= maxNumIters; ++iter ) {
    // Check convergence and output iteration
    const ScalarMag r_nrm = Thyra::norm(*r);          // Compute ||r|| = sqrt(<r,r>)
    const bool isConverged = r_nrm/r0_nrm <= tolerance;
    if( iter%(maxNumIters/10+1) == 0 || iter == maxNumIters || isConverged ) {
      if(out) *out << "Iter = " << iter << ", ||b-A*x||/||b-A*x0|| = " << (r_nrm/r0_nrm) << std::endl;
      if( r_nrm/r0_nrm < tolerance ) return true;     // Converged to tolerance, Success!
    }
    // Compute iteration
    const Scalar rho = space->scalarProd(*r,*r);      // <r,r>              -> rho
    if(iter==0) Thyra::assign(&*p,*r);                // r                  -> p   (iter == 0)
    else Thyra::Vp_V( &*p, *r, Scalar(rho/rho_old) ); // r+(rho/rho_old)*p  -> p   (iter  > 0)
    A.apply(Thyra::NOTRANS,*p,&*q);                   // A*p                -> q
    const Scalar alpha = rho/space->scalarProd(*p,*q);// rho/<p,q>          -> alpha
    Thyra::Vp_StV( x,   Scalar(+alpha), *p );         // +alpha*p + x       -> x
    Thyra::Vp_StV( &*r, Scalar(-alpha), *q );         // -alpha*q + r       -> r
    rho_old = rho;                                    // rho                -> rho_old (remember rho for next iter)
  }
  return false; // Failure
} // end sillyCgSolve

#endif // THYRA_SILLY_CG_SOLVE_HPP
