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

#include "Teuchos_Handle.hpp"
#include "Thyra_ConfigDefs.hpp"
#include "Thyra_LinearOperatorDecl.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

namespace Thyra
{
  /* Return the domain */
  template <class RangeScalar, class DomainScalar> inline 
  const VectorSpace<DomainScalar> ConstLinearOperator<RangeScalar, DomainScalar>
  ::domain() const 
  {return this->constPtr()->domain();}
    
  
  /* Return the range */
  template <class RangeScalar, class DomainScalar> inline 
  const VectorSpace<RangeScalar> ConstLinearOperator<RangeScalar, DomainScalar>
  ::range() const 
  {return this->constPtr()->range();}

  template <class RangeScalar, class DomainScalar> inline 
  void ConstLinearOperator<RangeScalar, DomainScalar>
  ::apply(const ConstVector<DomainScalar>& in,
          Vector<RangeScalar>& out,
          const RangeScalar& alpha,
          const RangeScalar& beta) const
  {
    /* the result vector might not be initialized. If it's null,
     * create a new vector in the range space */
    if (out.ptr().get()==0)
      {
        out = this->range().createMember();
      }
    this->constPtr()->apply(NONCONJ_ELE, *(in.constPtr().get()),
                       out.ptr().get(), alpha, beta);
  }


  template <class RangeScalar, class DomainScalar> inline 
  void ConstLinearOperator<RangeScalar, DomainScalar>
  ::applyTranspose(const ConstVector<RangeScalar>& in,
                   Vector<DomainScalar>& out,
                   const DomainScalar& alpha,
                   const DomainScalar& beta) const
  {
    /* the result vector might not be initialized. If it's null,
     * create a new vector in the domain space (i.e., the range space
     * of the transpose operator */
    if (out.ptr().get()==0)
      {
        out = this->domain().createMember();
      }
    this->constPtr()->applyTranspose(NONCONJ_ELE, *(in.constPtr().get()),
                                out.ptr().get(), 
                                alpha, beta);
  }
  
  
  template <class RangeScalar, class DomainScalar> inline 
  int ConstLinearOperator<RangeScalar, DomainScalar>::numBlockRows() const
  {
    return range().numBlocks();
  }
  
  
  template <class RangeScalar, class DomainScalar> inline 
  int ConstLinearOperator<RangeScalar, DomainScalar>::numBlockCols() const
  {
    return domain().numBlocks();
  }

  template <class RangeScalar, class DomainScalar> inline 
  ConstLinearOperator<RangeScalar, DomainScalar> 
  ConstLinearOperator<RangeScalar, DomainScalar>::getBlock(int blockRow, 
                                                           int blockCol) const
  {
    const Thyra::BlockedLinearOpBase<RangeScalar, DomainScalar>* p = 
      dynamic_cast<const Thyra::BlockedLinearOpBase<RangeScalar, DomainScalar>* >(this->constPtr().get());
    TEST_FOR_EXCEPTION(p==0 && blockRow != 0, runtime_error,
                       "request for block row=" << blockRow << " in a non-block operator");
    
    TEST_FOR_EXCEPTION(p==0 && blockCol != 0, runtime_error,
                       "request for block col=" << blockCol << " in a non-block operator");
    
    if (p != 0)
      {
        return p->getBlock(blockRow, blockCol);
      }
    return *this;
  }
  


  template <class RangeScalar, class DomainScalar> inline 
  LinearOperator<RangeScalar, DomainScalar> 
  LinearOperator<RangeScalar, DomainScalar>::getBlock(int blockRow, 
                                                      int blockCol) 
  {
    Thyra::BlockedLinearOpBase<RangeScalar, DomainScalar>* p = 
      dynamic_cast<Thyra::BlockedLinearOpBase<RangeScalar, DomainScalar>* >(this->ptr().get());
    TEST_FOR_EXCEPTION(p==0 && blockRow != 0, runtime_error,
                       "request for block row=" << blockRow << " in a non-block operator");
    
    TEST_FOR_EXCEPTION(p==0 && blockCol != 0, runtime_error,
                       "request for block col=" << blockCol << " in a non-block operator");
    
    if (p != 0)
      {
        return p->getNonconstBlock(blockRow, blockCol);
      }
    return *this;
  }
  

  template <class Scalar> inline 
  ConstLinearOperator<Scalar>
  block2x2(const ConstLinearOperator<Scalar>& A00,
           const ConstLinearOperator<Scalar>& A01,
           const ConstLinearOperator<Scalar>& A10,
           const ConstLinearOperator<Scalar>& A11)
  {
    return block2x2(A00.constPtr(), A01.constPtr(), A10.constPtr(), A11.constPtr());
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
    return block2x2(A00.ptr(), 
                    A01.ptr(), 
                    A10.ptr(), 
                    A11.ptr());
  }
  

  template <class Scalar> inline 
  LinearOperator<Scalar>
  block2x1(const LinearOperator<Scalar>& A00,
           const LinearOperator<Scalar>& A10)
  {
    return block2x1(A00.ptr(), 
                    A10.ptr());
  }
  

  template <class Scalar> inline 
  LinearOperator<Scalar>
  block1x2(const LinearOperator<Scalar>& A00,
           const LinearOperator<Scalar>& A01)
  {
    return block2x1(A00.ptr(), 
                    A01.ptr());
  }
  

  


  
  
}

#endif
