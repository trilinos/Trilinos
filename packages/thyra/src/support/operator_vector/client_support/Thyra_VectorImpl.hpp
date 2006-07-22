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

#ifndef THYRA_VECTOR_IMPL_HPP
#define THYRA_VECTOR_IMPL_HPP

#include "Thyra_VectorDecl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearCombinationImpl.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_DetachedVectorView.hpp"


namespace Thyra
{
  /** */
  template <class Scalar> inline
  Scalar ConstVector<Scalar>::operator[](Index globalIndex) const 
  {
    ConstDetachedVectorView<Scalar> view(this->ptr(), Range1D(0, dim(*this)-1));
    return view[globalIndex];
  }

  template <class Scalar> inline
  bool ConstVector<Scalar>::containsVector(const Thyra::VectorBase<Scalar>* vec) const
  {
    return this->ptr().get()==vec;
  }
  

  template <class Scalar> inline
  void ConstVector<Scalar>::evalInto(Vector<Scalar>& acceptor) const
  {
    acceptor.acceptCopyOf(*this);
  }

  template <class Scalar> inline
  void ConstVector<Scalar>::addInto(Vector<Scalar>& acceptor, LCSign sign) const
  {
    axpy(sign, *this, acceptor);
  }


  /* */
  template <class Scalar> inline
  std::ostream& operator<<(std::ostream& os, const ConstVector<Scalar>& v)
  {
    os << v.description() ;
    return os;
  }

  template <class Scalar> inline
  Vector<Scalar>& Vector<Scalar>::acceptCopyOf(const ConstVector<Scalar>& other)
  {
    Thyra::VectorBase<Scalar>* p = this->ptr().get();
    const Thyra::VectorBase<Scalar>* px = other.ptr().get();
    
    if (p==0) 
      {
        Vector<Scalar> me = space(other).createMember();
        this->ptr() = me.ptr();
      }
    Thyra::assign(p, *px);
    return *this;
  }


  /** 
   * \relates ConstVector
   * Return the dimension of the vector 
   */
  template <class Scalar> inline
  Index dim(const ConstVector<Scalar>& x) 
  {
    return x.ptr()->space()->dim();
  }

  
  /** 
   * \relates ConstVector
   * Return the vector space in which the vector lives */
  template <class Scalar> inline
  VectorSpace<Scalar> space(const ConstVector<Scalar>& x) 
  {
    return x.ptr()->space();
  }

  /* copy */
  THYRA_UNARY_VECTOR_OP(copy, copyInto, assign, "copy")


  //===========================================================================
  template <class Scalar> inline
  ConstVector<Scalar> ConstVector<Scalar>::getBlock(int i) const
  {
    const Thyra::ProductVectorBase<Scalar>* pv = 
      dynamic_cast <const Thyra::ProductVectorBase<Scalar>* >(this->ptr().get());
    if (pv==0) 
      {
        TEST_FOR_EXCEPTION(i != 0, runtime_error,
                           "Nonzero block index " << i << " into a vector that is not "
                           "a product vector");
        return *this;
      }
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > b = pv->getVectorBlock(i);
    return b;
  }

  //===========================================================================
  template <class Scalar> inline
  Vector<Scalar> Vector<Scalar>::getBlock(int i) 
  {
    Thyra::ProductVectorBase<Scalar>* pv = 
      dynamic_cast <Thyra::ProductVectorBase<Scalar>* >(this->ptr().get());
    if (pv==0) 
      {
        TEST_FOR_EXCEPTION(i != 0, runtime_error,
                           "Nonzero block index " << i << " into a vector that is not "
                           "a product vector");
        return *this;
      }
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > b = pv->getNonconstVectorBlock(i);
    return b;
  }

  /** \relates Vector */
  template <class Scalar> inline  
  void setElement(const Vector<Scalar>& x, Index i, const Scalar& x_i);

  /** \relates ConstVector */
  template <class Scalar> inline  
  const Scalar& getElement(const ConstVector<Scalar>& x, Index i);

  
  


}

#endif
