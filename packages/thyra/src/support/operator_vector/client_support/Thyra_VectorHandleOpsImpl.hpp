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

#ifndef THYRA_VECTOR_HANDLE_OPS_IMPL_HPP
#define THYRA_VECTOR_HANDLE_OPS_IMPL_HPP

#include "Thyra_ConfigDefs.hpp"
#include "Thyra_VectorHandleOpsDecl.hpp"

namespace Thyra
{
  
  /* */
  template <class Scalar> inline  
  void setToConstant(Vector<Scalar>& x, const Scalar& a)
  {
    put_scalar(a, x.ptr().get());
  }

  /* */
  template <class Scalar> inline  
  void zeroOut(Vector<Scalar>& x)
  {
    setToConstant(x, 0.0);
  }

  

  /* abs */
  THYRA_UNARY_VECTOR_OP(abs, absInto, abs, "absolute value")
    
  /* reciprocal */
    THYRA_UNARY_VECTOR_OP(reciprocal, reciprocalInto, reciprocal, "reciprocal")
    
  /* */
    THYRA_UNARY_VECTOR_ROP(norm1, norm_1, "1-norm")

  /* */
    THYRA_UNARY_VECTOR_ROP(norm2, norm_2, "2-norm")

  /* */
    THYRA_UNARY_VECTOR_ROP(normInf, norm_inf, "inf-norm")

  /* */ 
    THYRA_UNARY_VECTOR_ROP(sum, sum, "sum of the elements")

  /* */ 
    THYRA_UNARY_VECTOR_ROP(max, max, "max element")

  /* */ 
    THYRA_UNARY_VECTOR_ROP(min, min, "min element")

  /**
   * \relates ConstVector
   * 
   * Return the max of a vector and its location
   */
    template <class Scalar> inline 
  Scalar max(const Converter<Scalar, ConstVector<Scalar> >& x, int& index) 
  {
    Scalar maxEl;
    Scalar* maxElP = &maxEl;
    int* indexP = &index;
    Thyra::max(*(toVector(x).ptr()), maxElP, indexP); 
    return maxEl;
  }

  /*
   * Return the min of a vector and its location
   */
  template <class Scalar> inline 
  Scalar min(const Converter<Scalar, ConstVector<Scalar> >& x, int& index) 
  {
    Scalar minEl;
    Scalar* minElP = &minEl;
    int* indexP = &index;
    Thyra::min(*(toVector(x).ptr()), minElP, indexP); 
    return minEl;
  }

  /*
   *
   */
  template <class Scalar> inline  
  Scalar min(const Converter<Scalar, ConstVector<Scalar> >& x, 
             const Scalar& bound, int& index)
  {
    Scalar minEl;
    Scalar* minElP = &minEl;
    int* indexP = &index;
    Thyra::minGreaterThanBound(*(toVector(x).ptr()), bound, minElP, indexP); 
    return minEl;
  }


  /* 
   * 
   */
  template <class Scalar> inline  
  Scalar max(const Converter<Scalar, ConstVector<Scalar> >& x, 
             const Scalar& bound, int& index)
  {
    Scalar maxEl;
    Scalar* maxElP = &maxEl;
    int* indexP = &index;
    Thyra::maxLessThanBound(*(toVector(x).ptr()), bound, maxElP, indexP); 
    return maxEl;
  }

  /*
   * 
   */
  template <class Scalar> inline  
  void dotStarInto(const Converter<Scalar, ConstVector<Scalar> >& x, 
                   const Converter<Scalar, ConstVector<Scalar> >& y,
                   Vector<Scalar>& result)
    {
      zeroOut(result);
      Thyra::ele_wise_prod(1.0, *(toVector(x).ptr()), *(toVector(y).ptr()), 
                           result.ptr().get());
    }

  /*
   * 
   */
  template <class Scalar> inline  
  Vector<Scalar> dotStar(const Converter<Scalar, ConstVector<Scalar> >& xIn, 
                         const Converter<Scalar, ConstVector<Scalar> >& y)
    {
      ConstVector<Scalar> x = toVector(xIn);
      Vector<Scalar> result = space(x).createMember();
      dotStarInto(x, toVector(y), result);
      return result;
    }



  /*
   * 
   */
  template <class Scalar> inline  
  void dotSlashInto(const Converter<Scalar, ConstVector<Scalar> >& x, 
                   const Converter<Scalar, ConstVector<Scalar> >& y,
                   Vector<Scalar>& result)
    {
      Thyra::ele_wise_divide(1.0, *(toVector(x).ptr()), *(toVector(y).ptr()), 
                             result.ptr().get());
    }

  /*
   * 
   */
  template <class Scalar> inline  
  Vector<Scalar> dotSlash(const Converter<Scalar, ConstVector<Scalar> >& xIn, 
                          const Converter<Scalar, ConstVector<Scalar> >& y)
    {
      ConstVector<Scalar> x = toVector(xIn);
      Vector<Scalar> result = space(x).createMember();
      dotSlashInto(x, toVector(y), result);
      return result;
    }

  
  /*
   * 
   */
  template <class Scalar> inline  
  void axpy(const double& a, const Converter<Scalar, ConstVector<Scalar> >& xIn, 
            Vector<Scalar>& y)
  {
    ConstVector<Scalar> x = toVector(xIn);
    VectorBase<Scalar>* p = y.ptr().get();
    const VectorBase<Scalar>* px = x.ptr().get();
    Vp_StV(p, a, *px);
  }

  /*
   * 
   */
  template <class Scalar> inline  
  void scale(Vector<Scalar>& x, const double& a)
  {
    VectorBase<Scalar>* p = x.ptr().get();
    Thyra::scale(a, p);
  }

  /*
   * 
   */
  template <class Scalar> inline     
  void scaleInto(const Converter<Scalar, ConstVector<Scalar> >& x,
                 const double& alpha, Vector<Scalar>& result)
  {
    x.evalInto(result);
    result.scale(alpha);
  }


  
 
}

#endif
