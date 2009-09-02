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

#ifndef THYRA_SILLY_MODIFIED_GRAM_SHMIDT_HPP
#define THYRA_SILLY_MODIFIED_GRAM_SHMIDT_HPP

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

/** \brief Silly little implementation of the modified Gram-Schmidt algorithm
 * to compute a QR factorization V=Q*R of a multi-vector V.
 *
 * \param   V   [in/out] On input, contains the columns to compute the factorization
 *              for.  On output, contains the columns of Q.
 * \param   R   [out] On output, contains the upper triangular matrix R.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_examples_cg_grp
 */
template<class Scalar>
void sillyModifiedGramSchmidt(
  Thyra::MultiVectorBase<Scalar> *V_inout
  ,Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > *R_out
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT(V_inout==NULL);
  Thyra::MultiVectorBase<Scalar> &V = *V_inout;
  const int n = V.domain()->dim();
  *R_out = Thyra::createMembers(V.domain(),n);
  //Thyra::assign(&*(*R_out),ST::zero());
  Thyra::DetachedMultiVectorView<Scalar> R(*(*R_out));
  for( int k = 0; k < n; ++k ) {
    R(k,k) = Thyra::norm(*V.col(k));
    Thyra::scale(Scalar(ST::one()/R(k,k)),&*V.col(k));
    for( int j = k+1; j < n; ++j ) {
      R(k,j) = Thyra::scalarProd(*V.col(k),*V.col(j));
      Thyra::update( Scalar(-R(k,j)), *V.col(k), &*V.col(j) );
    }
  }

} // end sillyModifiedGramSchmidt

#endif // THYRA_SILLY_MODIFIED_GRAM_SHMIDT_HPP
