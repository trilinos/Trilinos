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

#ifndef THYRA_VECTOR_HANDLE_OPS_DECL_HPP
#define THYRA_VECTOR_HANDLE_OPS_DECL_HPP

#include "Thyra_ConfigDefs.hpp"
#include "Thyra_VectorDecl.hpp"
#include "Thyra_VecOpMacros.hpp"

namespace Thyra 
{
  /** \relates Vector */
  template <class Scalar> inline  
  void setToConstant(const Vector<Scalar>& x, const Scalar& a);

  /** \relates Vector */
  template <class Scalar> inline  
  void zeroOut(const Vector<Scalar>& x);

  /* abs */
  THYRA_UNARY_VECTOR_OP_DECL(abs, absInto, abs, "absolute value");

  /* reciprocal */
  THYRA_UNARY_VECTOR_OP_DECL(reciprocal, reciprocalInto, reciprocal, "reciprocal");

  /* */
  THYRA_UNARY_VECTOR_ROP_DECL(norm1, norm_1, "1-norm");

  /* */
  THYRA_UNARY_VECTOR_ROP_DECL(norm2, norm_2, "2-norm");

  /* */
  THYRA_UNARY_VECTOR_ROP_DECL(normInf, norm_inf, "inf-norm");
  //template <class Scalar> 
  //  Scalar normInf(const Thyra::ConvertibleToVector<Scalar>& x);

  /* */ 
  THYRA_UNARY_VECTOR_ROP_DECL(sum, sum, "sum of the elements");

  /* */ 
  THYRA_UNARY_VECTOR_ROP_DECL(max, max, "max element");

  /* */ 
  THYRA_UNARY_VECTOR_ROP_DECL(min, min, "min element");

  /** 
   * \relates Vector
   */
  template <class Scalar>   
  Vector<Scalar> dotStar(const Converter<Scalar, ConstVector<Scalar> >& x, 
                         const Converter<Scalar, ConstVector<Scalar> >& y);


  /** 
   * \relates Vector
   */
  template <class Scalar>   
  Vector<Scalar> dotSlash(const Converter<Scalar, ConstVector<Scalar> >& x, 
                          const Converter<Scalar, ConstVector<Scalar> >& y);

  /**
   * \relates Vector
   * 
   * Return the max of a vector and its location
   */
  template <class Scalar>  
  Scalar maxloc(const Converter<Scalar, ConstVector<Scalar> >& x, int& index) ;

  /**
   * \relates Vector
   * \brief.
   * Return the min of a vector and its location
   */
  template <class Scalar>  
  Scalar minloc(const Converter<Scalar, ConstVector<Scalar> >& x, int& index) ;

  /** 
   * \relates Vector<Scalar>
   * \brief.
   *
   */
  template <class Scalar>   
  Scalar minloc(const Converter<Scalar, ConstVector<Scalar> >& x, 
             const Scalar& bound, int& index) ;


  /** 
   * \relates Vector
   * \brief.
   * 
   */
  template <class Scalar>   
  Scalar maxloc(const Converter<Scalar, ConstVector<Scalar> >& x, 
                const Scalar& bound, int& index);


  /** */
  template <class Scalar>   
  void scaleInto(const Converter<Scalar, ConstVector<Scalar> >& x,
                 const double& alpha, Vector<Scalar>& result);

  /** */
  template <class Scalar>   
  void scale(Vector<Scalar>& x, const double& alpha);

  /** */
  template <class Scalar>   
  void axpy(const double& a, const Converter<Scalar, ConstVector<Scalar> >& x, 
            Vector<Scalar>& y);
  
}

#endif
