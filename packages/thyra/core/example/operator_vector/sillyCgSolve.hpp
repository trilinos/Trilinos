// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SILLY_CG_SOLVE_HPP
#define THYRA_SILLY_CG_SOLVE_HPP

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"


/** \brief Silly little example unpreconditioned CG solver.
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
  const Thyra::LinearOpBase<Scalar> &A,
  const Thyra::VectorBase<Scalar> &b,
  const int maxNumIters,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x,
  std::ostream &out
  )
{

  // Create some typedefs and some other stuff to make the code cleaner
  typedef Teuchos::ScalarTraits<Scalar> ST; typedef typename ST::magnitudeType ScalarMag;
  const Scalar one = ST::one(), zero = ST::zero();  using Teuchos::as;
  using Teuchos::RCP; using Thyra::VectorSpaceBase; using Thyra::VectorBase;
  using Thyra::NOTRANS; using Thyra::V_V; using Thyra::apply;


  // Validate input
  THYRA_ASSERT_LINEAR_OP_VEC_APPLY_SPACES("sillyCgSolve()", A, Thyra::NOTRANS, *x, &b);
  Teuchos::EVerbosityLevel vl = Teuchos::VERB_MEDIUM;

  std::ios::fmtflags fmt(out.flags());

  out << "\nStarting CG solver ...\n" << std::scientific << "\ndescribe A:\n"<<describe(A, vl)
      << "\ndescribe b:\n"<<describe(b, vl)<<"\ndescribe x:\n"<<describe(*x, vl)<<"\n";

  // Initialization
  const RCP<const VectorSpaceBase<Scalar> > space = A.domain();
  const RCP<VectorBase<Scalar> > r = createMember(space);
  // r = -A*x + b
  V_V(r.ptr(), b); apply<Scalar>(A, NOTRANS, *x, r.ptr(), -one, one);
  const ScalarMag r0_nrm = norm(*r);
  if (r0_nrm==zero) return true;
  const RCP<VectorBase<Scalar> > p = createMember(space), q = createMember(space);
  Scalar rho_old = -one;

  // Perform the iterations
  for( int iter = 0; iter <= maxNumIters; ++iter ) {

    // Check convergence and output iteration
    const ScalarMag r_nrm = norm(*r);
    const bool isConverged = r_nrm/r0_nrm <= tolerance;
    if( iter%(maxNumIters/10+1) == 0 || iter == maxNumIters || isConverged ) {
      out << "Iter = " << iter << ", ||b-A*x||/||b-A*x0|| = " << (r_nrm/r0_nrm) << std::endl;
      if( r_nrm/r0_nrm < tolerance ) {
        out.flags(fmt);
        return true; // Success!
      }
    }

    // Compute iteration
    const Scalar rho = inner(*r, *r);        // <r,r>              -> rho
    if (iter==0) V_V(p.ptr(), *r);           // r                  -> p   (iter == 0)
    else Vp_V( p.ptr(), *r, rho/rho_old );   // r+(rho/rho_old)*p  -> p   (iter  > 0)
    apply<Scalar>(A, NOTRANS, *p, q.ptr());  // A*p                -> q
    const Scalar alpha = rho/inner(*p, *q);  // rho/<p,q>          -> alpha
    Vp_StV( x, +alpha, *p );                 // +alpha*p + x       -> x
    Vp_StV( r.ptr(), -alpha, *q );           // -alpha*q + r       -> r
    rho_old = rho;                           // rho                -> rho_old (for next iter)

  }

  out.flags(fmt);
  return false; // Failure
} // end sillyCgSolve

#endif // THYRA_SILLY_CG_SOLVE_HPP
