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

#ifndef THYRA_SILLIER_CG_SOLVE_HPP
#define THYRA_SILLIER_CG_SOLVE_HPP

#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_VectorImpl.hpp"
#include "Thyra_AssertOp.hpp"
#include "silliestCgSolve.hpp"


/** \brief Silly little example unpreconditioned CG solver (calls templated code).
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
bool sillierCgSolve(
  const Thyra::LinearOpBase<Scalar> &A_in,
  const Thyra::VectorBase<Scalar> &b_in,
  const int maxNumIters,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x_inout,
  std::ostream &out
  )
{

  // Validate the input
  THYRA_ASSERT_LINEAR_OP_VEC_APPLY_SPACES("sillyCgSolve()", A_in,
    Thyra::NOTRANS, *x_inout, &b_in);

  // Create handle wrappers to facilitate the use of operator overloading
  const Thyra::ConstLinearOperator<Scalar> A(Teuchos::rcpFromRef(A_in));
  const Thyra::ConstVector<Scalar> b(Teuchos::rcpFromRef(b_in));
  Thyra::Vector<Scalar> x(Teuchos::rcpFromPtr(x_inout));

  // Describe the arguments
  Teuchos::EVerbosityLevel vl = Teuchos::VERB_MEDIUM;
  out << "\nStarting CG solver ...\n" << std::scientific << "\ndescribe A:\n"<<describe(A,vl)
      << "\ndescribe b:\n"<<describe(b,vl)<<"\ndescribe x:\n"<<describe(x,vl)<<"\n";
  
  return silliestCgSolve(A, b, maxNumIters, tolerance, x, out);

} // end sillierCgSolve


#endif // THYRA_SILLIER_CG_SOLVE_HPP
