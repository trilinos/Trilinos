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

#ifndef THYRA_VECOPMACROS_HPP
#define THYRA_VECOPMACROS_HPP

#include "Thyra_ConfigDefs.hpp"


/** \brief Boilerplate code for defining unary vector transformation
 * operations.
 */ 
#define THYRA_UNARY_VECTOR_OP(opName, intoOpName, thyraOpName, descr)   \
  template <class Scalar> inline                                        \
  Thyra::Vector<Scalar> opName(const Converter<Scalar, ConstVector<Scalar> >& x) \
  {                                                                     \
    Thyra::ConstVector<Scalar> y=Thyra::toVector(x);                    \
    Thyra::Vector<Scalar> rtn = space(y).createMember();                \
    intoOpName(y, rtn);                                                 \
    return rtn;                                                         \
  }                                                                     \
  template <class Scalar> inline                                        \
  void intoOpName(const Converter<Scalar, ConstVector<Scalar> >& donorSource,   \
                  Vector<Scalar>& acceptor)                             \
  {                                                                     \
    Thyra::ConstVector<Scalar> donor=Thyra::toVector(donorSource);      \
    if (acceptor.ptr().get()==0)                                        \
      {                                                                 \
        acceptor = space(donor).createMember();                         \
      }                                                                 \
    Thyra::VectorBase<Scalar>* p = acceptor.ptr().get();                \
    const Thyra::VectorBase<Scalar>* px = donor.constPtr().get();            \
    Thyra::thyraOpName(p, *px);                                         \
  }      

      
/** \brief Boilerplate code for declaring unary vector transformation
 * operations.
 */                                                         
#define THYRA_UNARY_VECTOR_OP_DECL(opName, intoOpName, thyraOpName, descr) \
  template <class Scalar>                                               \
  Vector<Scalar> opName(const Converter<Scalar, ConstVector<Scalar> >& x);      \
  template <class Scalar>                                               \
  void intoOpName(const Converter<Scalar, ConstVector<Scalar> >& donor,         \
                  Vector<Scalar>& acceptor)

/** \brief Boilerplate code for defining unary vector reduction operations.
 */ 
#define THYRA_UNARY_VECTOR_ROP(funcName, ROpName, descr)    \
  template <class Scalar> inline                            \
  Scalar funcName(const Converter<Scalar, ConstVector<Scalar> >& x) \
  {                                                         \
    Thyra::ConstVector<Scalar> y = Thyra::toVector(x);      \
    return Thyra::ROpName(*(y.constPtr()));                      \
  }                                                                    

/** \brief Boilerplate code for declaring unary vector reduction operations.
 */ 
#define THYRA_UNARY_VECTOR_ROP_DECL(funcName, ROpName, descr) \
  template <class Scalar>                                     \
  Scalar funcName(const Converter<Scalar, ConstVector<Scalar> >& x)           

#endif
