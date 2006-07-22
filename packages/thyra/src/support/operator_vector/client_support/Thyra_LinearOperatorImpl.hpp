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

namespace Thyra
{
  /* Return the domain */
  template <class RangeScalar, class DomainScalar> inline 
  const VectorSpace<DomainScalar> ConstLinearOperator<RangeScalar, DomainScalar>
  ::domain() const 
  {return this->ptr()->domain();}
    
  
  /* Return the range */
  template <class RangeScalar, class DomainScalar> inline 
  const VectorSpace<RangeScalar> ConstLinearOperator<RangeScalar, DomainScalar>
  ::range() const 
  {return this->ptr()->domain();}

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
    this->ptr()->apply(NONCONJ_ELE, *(in.ptr().get()),
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
    this->ptr()->applyTranspose(NONCONJ_ELE, *(in.ptr().get()),
                                out.ptr().get(), alpha, beta);
  }
  
 
  
}

#endif
