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

#ifndef THYRA_LINEAROPERATOR_DECL_HPP
#define THYRA_LINEAROPERATOR_DECL_HPP

#include "Teuchos_Handle.hpp"
#include "Thyra_ConfigDefs.hpp"
#include "Thyra_SingleScalarLinearOpBaseDecl.hpp"

namespace Thyra
{
  /** 
   *
   */
  template <class Scalar>
  class ConstLinearOperator 
    : public virtual Teuchos::ConstHandle<SingleScalarLinearOpBase<Scalar> >
  {
  public:
    TEUCHOS_CONST_HANDLE_CTORS(ConstLinearOperator<Scalar>, 
                               SingleScalarLinearOpBase<Scalar>);

    /** Return the domain */
    const VectorSpace<Scalar> domain() const ;

    /** Return the range */
    const VectorSpace<Scalar> range() const ;


    /** 
     * Compute
     * \code
     * out = beta*out + alpha*op*in;
     * \endcode
     **/
    void apply(const ConstVector<Scalar>& in,
	       Vector<Scalar>& out,
	       const Scalar& alpha = 1.0,
	       const Scalar& beta = 0.0) const ;

    /**  
     * Compute
     * \code
     * out = beta*out + alpha*op^T*in;
     * \endcode
     **/
    void applyTranspose(const ConstVector<Scalar>& in,
			Vector<Scalar>& out,
			const Scalar& alpha = 1.0,
			const Scalar& beta = 0.0) const ;

  };


  /** 
   *
   */
  template <class Scalar>
  class LinearOperator 
    : public Teuchos::Handle<SingleScalarLinearOpBase<Scalar> >,
      public ConstLinearOperator<Scalar>
  {
  public:
    TEUCHOS_HANDLE_CTORS(LinearOperator<Scalar>, 
                         SingleScalarLinearOpBase<Scalar>);
  };

  
}

#endif
