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

#ifndef THYRA_LINEAROPERATOR_IMPL_HPP
#define THYRA_LINEAROPERATOR_IMPL_HPP

#include "Thyra_LinearOperatorDecl.hpp"
#include "Thyra_ConfigDefs.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"

namespace Thyra
{

  template<class Scalar> inline 
  const VectorSpace<Scalar> ConstLinearOperator<Scalar>
  ::domain() const 
  {return this->constPtr()->domain();}
    
  template<class Scalar> inline 
  const VectorSpace<Scalar> ConstLinearOperator<Scalar>
  ::range() const 
  {return this->constPtr()->range();}

  template<class Scalar> inline 
  void ConstLinearOperator<Scalar>
  ::apply(const ConstVector<Scalar>& in,
          Vector<Scalar>& out,
          const Scalar& alpha,
          const Scalar& beta) const
  {
    /* the result std::vector might not be initialized. If it's null,
     * create a new std::vector in the range space */
    if (out.ptr().get()==0)
      {
        out = this->range().createMember();
      }
    this->constPtr()->apply(NOTRANS, *(in.constPtr().get()),
                       out.ptr().ptr(), alpha, beta);
  }

  template<class Scalar> inline 
  void ConstLinearOperator<Scalar>
  ::applyTranspose(const ConstVector<Scalar>& in,
                   Vector<Scalar>& out,
                   const Scalar& alpha,
                   const Scalar& beta) const
  {
    /* the result std::vector might not be initialized. If it's null,
     * create a new std::vector in the domain space (i.e., the range space
     * of the transpose operator */
    if (out.ptr().get()==0)
      {
        out = this->domain().createMember();
      }
    this->constPtr()->applyTranspose(NONCONJ_ELE, *(in.constPtr().get()),
                                out.ptr().get(), 
                                alpha, beta);
  }
  
  template<class Scalar> inline 
  int ConstLinearOperator<Scalar>::numBlockRows() const
  {
    return range().numBlocks();
  }
  
  template<class Scalar> inline 
  int ConstLinearOperator<Scalar>::numBlockCols() const
  {
    return domain().numBlocks();
  }

  template<class Scalar> inline 
  ConstLinearOperator<Scalar> 
  ConstLinearOperator<Scalar>::getBlock(int blockRow, 
                                                           int blockCol) const
  {
    const Thyra::BlockedLinearOpBase<Scalar>* p = 
      dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar>* >(this->constPtr().get());
    TEST_FOR_EXCEPTION(p==0 && blockRow != 0, std::runtime_error,
                       "request for block row=" << blockRow << " in a non-block operator");
    
    TEST_FOR_EXCEPTION(p==0 && blockCol != 0, std::runtime_error,
                       "request for block col=" << blockCol << " in a non-block operator");
    
    if (p != 0)
      {
        return p->getBlock(blockRow, blockCol);
      }
    return *this;
  }

  template<class Scalar> inline 
  LinearOperator<Scalar> 
  LinearOperator<Scalar>::getBlock(int blockRow, 
                                                      int blockCol) 
  {
    Thyra::BlockedLinearOpBase<Scalar>* p = 
      dynamic_cast<Thyra::BlockedLinearOpBase<Scalar>* >(this->ptr().get());
    TEST_FOR_EXCEPTION(p==0 && blockRow != 0, std::runtime_error,
                       "request for block row=" << blockRow << " in a non-block operator");
    
    TEST_FOR_EXCEPTION(p==0 && blockCol != 0, std::runtime_error,
                       "request for block col=" << blockCol << " in a non-block operator");
    
    if (p != 0)
      {
        return p->getNonconstBlock(blockRow, blockCol);
      }
    return *this;
  }

  //
  // Nonmember functions
  //

  template <class Scalar> inline 
  ConstLinearOperator<Scalar>
  operator*(const Scalar& a, 
            const ConstLinearOperator<Scalar>& A)
  {
    return Thyra::scale<Scalar>(a,A.constPtr());
  }

  template <class Scalar> inline 
  LinearOperator<Scalar>
  operator*(const Scalar& a, 
            const LinearOperator<Scalar>& A)
  {
    return Thyra::nonconstScale<Scalar>(a,A.ptr());
  }
  
  template <class Scalar> inline 
  ConstLinearOperator<Scalar>
  operator*(const ConstLinearOperator<Scalar>& A,
            const ConstLinearOperator<Scalar>& B)
  {
    return Thyra::multiply<Scalar>(A.constPtr(),B.constPtr());
  }
  
  template <class Scalar> inline 
  LinearOperator<Scalar>
  operator*(const LinearOperator<Scalar>& A,
            const LinearOperator<Scalar>& B)
  {
    return Thyra::nonconstMultiply<Scalar>(A.ptr(),B.ptr());
  }
 
  template <class Scalar> inline 
  ConstLinearOperator<Scalar>
  operator+(const ConstLinearOperator<Scalar>& A,
            const ConstLinearOperator<Scalar>& B)
  {
    return Thyra::add<Scalar>(A.constPtr(),B.constPtr());
  }
  
  template <class Scalar> inline 
  LinearOperator<Scalar>
  operator+(const LinearOperator<Scalar>& A,
            const LinearOperator<Scalar>& B)
  {
    return Thyra::nonconstAdd<Scalar>(A.ptr(),B.ptr());
  }
  
  template <class Scalar> inline 
  ConstLinearOperator<Scalar>
  block2x2(const ConstLinearOperator<Scalar>& A00,
           const ConstLinearOperator<Scalar>& A01,
           const ConstLinearOperator<Scalar>& A10,
           const ConstLinearOperator<Scalar>& A11)
  {
    return block2x2(A00.constPtr(),A01.constPtr(),A10.constPtr(),A11.constPtr());
  }

  template <class Scalar> inline 
  ConstLinearOperator<Scalar>
  block2x1(const ConstLinearOperator<Scalar>& A00,
           const ConstLinearOperator<Scalar>& A10)
  {
    return block2x1(A00.constPtr(), A10.constPtr());
  }

  template <class Scalar> inline 
  ConstLinearOperator<Scalar>
  block1x2(const ConstLinearOperator<Scalar>& A00,
           const ConstLinearOperator<Scalar>& A01)
  {
    return block2x1(A00.constPtr(), A01.constPtr());
  }

  template <class Scalar> inline 
  LinearOperator<Scalar>
  block2x2(const LinearOperator<Scalar>& A00,
           const LinearOperator<Scalar>& A01,
           const LinearOperator<Scalar>& A10,
           const LinearOperator<Scalar>& A11)
  {
    return nonconstBlock2x2(A00.ptr(),A01.ptr(),A10.ptr(),A11.ptr());
  }

  template <class Scalar> inline 
  LinearOperator<Scalar>
  block2x1(const LinearOperator<Scalar>& A00,
           const LinearOperator<Scalar>& A10)
  {
    return nonconstBlock2x1(A00.ptr(),A10.ptr());
  }

  template <class Scalar> inline 
  LinearOperator<Scalar>
  block1x2(const LinearOperator<Scalar>& A00,
           const LinearOperator<Scalar>& A01)
  {
    return nonconstBlock2x1(A00.ptr(),A01.ptr());
  }
  
} // namespace Thyra


template<class Scalar>
Thyra::ConstLinearOperator<Scalar>
Thyra::identity( const VectorSpace<Scalar> &space )
{
  return identity(space.constPtr());
}

template<class Scalar>
Thyra::ConstLinearOperator<Scalar>
Thyra::zero( const VectorSpace<Scalar> &range, const VectorSpace<Scalar> &domain )
{
  return zero(range.constPtr(),domain.constPtr());
}

#endif
