// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.  
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses, 
//         templates : new C++ techniques 
//            for scientific computing 
// 
//********************************************************
//
//  A short implementation ( not all operators and 
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#include "Sacado_ConfigDefs.h"
#include "Sacado_mpl_range_c.hpp"
#include "Sacado_mpl_for_each.hpp"

template <typename T, int Num> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
Expr(const int sz, const T & x) : val_(x)
{
#ifdef SACADO_DEBUG
  if (sz != Num)
    throw "SELRFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
#endif

  ss_array<T>::zero(dx_, Num); 
}

template <typename T, int Num> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
Expr(const int sz, const int i, const T & x) : val_(x) 
{
#ifdef SACADO_DEBUG
  if (sz != Num)
    throw "SELRFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
  if (i >= Num)
    throw "SELRFad::SFad() Error:  Invalid derivative index.";
#endif

  ss_array<T>::zero(dx_, Num);
  dx_[i]=1.; 
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
Expr(const Expr<S>& x) : val_(x.val())
{
#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SELRFad::SFad() Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  // Compute partials
  LocalAccumOp< Expr<S> > op(x);
  
  // Compute each tangent direction
  for(int i=0; i<Num; ++i) {
    op.t = T(0.);
    op.i = i;

    // Automatically unrolled loop that computes
    // for (int j=0; j<N; j++)
    //   op.t += op.partials[j] * x.getTangent<j>(i);
    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);

    dx_[i] = op.t;

  }
}


template <typename T, int Num> 
inline void 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
diff(const int ith, const int n) 
{ 
#ifdef SACADO_DEBUG
  if (n != Num)
    throw "SELRFad::diff() Error:  Supplied derivative dimension does not match compile time length.";
#endif

  ss_array<T>::zero(dx_, Num);
  dx_[ith] = T(1.);
}

template <typename T, int Num> 
inline void 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
resize(int sz)
{
#ifdef SACADO_DEBUG
  if (sz != Num)
    throw "SELRFad::resize() Error:  Cannot resize fixed derivative array dimension";
#endif
}

template <typename T, int Num> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator=(const T& v) 
{
  val_ = v;
  ss_array<T>::zero(dx_, Num);

  return *this;
}

template <typename T, int Num> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator=(const Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& x) 
{
  // Copy value
  val_ = x.val_;

  // Copy dx_
  //ss_array<T>::copy(x.dx_, dx_, Num);
  for (int i=0; i<Num; i++)
    dx_[i] = x.dx_[i];
  
  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator=(const Expr<S>& x) 
{
#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SELRFad::operator=() Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  // Compute partials
  LocalAccumOp< Expr<S> > op(x);
  
  // Compute each tangent direction
  for(int i=0; i<Num; ++i) {
    op.t = T(0.);
    op.i = i;

    // Automatically unrolled loop that computes
    // for (int j=0; j<N; j++)
    //   op.t += op.partials[j] * x.getTangent<j>(i);
    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
    
    dx_[i] = op.t;
  }
  
  // Compute value
  val_ = x.val();
  
  return *this;
}

template <typename T, int Num> 
inline  Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator += (const T& v)
{
  val_ += v;

  return *this;
}

template <typename T, int Num> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator -= (const T& v)
{
  val_ -= v;

  return *this;
}

template <typename T, int Num> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator *= (const T& v)
{
  val_ *= v;

  for (int i=0; i<Num; ++i)
    dx_[i] *= v;

  return *this;
}

template <typename T, int Num> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator /= (const T& v)
{
  val_ /= v;

  for (int i=0; i<Num; ++i)
    dx_[i] /= v;

  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator += (const Sacado::ELRFad::Expr<S>& x)
{
#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SELRFad::operator+=() Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  // Compute partials
  LocalAccumOp< Expr<S> > op(x);

  // Compute each tangent direction
  for(int i=0; i<Num; ++i) {
    op.t = T(0.);
    op.i = i;
	
    // Automatically unrolled loop that computes
    // for (int j=0; j<N; j++)
    //   op.t += op.partials[j] * x.getTangent<j>(i);
    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
    
    dx_[i] += op.t;
  }
 
  // Compute value
  val_ += x.val();

  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator -= (const Sacado::ELRFad::Expr<S>& x)
{
#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SELRFad::operator-=() Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  // Compute partials
  LocalAccumOp< Expr<S> > op(x);

  // Compute each tangent direction
  for(int i=0; i<Num; ++i) {
    op.t = T(0.);
    op.i = i;
	
    // Automatically unrolled loop that computes
    // for (int j=0; j<N; j++)
    //   op.t += op.partials[j] * x.getTangent<j>(i);
    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	
    dx_[i] -= op.t;
  }

  // Compute value
  val_ -= x.val();

  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator *= (const Sacado::ELRFad::Expr<S>& x)
{
  T xval = x.val();

#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SELRFad::operator*=() Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  // Compute partials
  LocalAccumOp< Expr<S> > op(x);

  // Compute each tangent direction
  for(int i=0; i<Num; ++i) {
    op.t = T(0.);
    op.i = i;
	
    // Automatically unrolled loop that computes
    // for (int j=0; j<N; j++)
    //   op.t += op.partials[j] * x.getTangent<j>(i);
    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
    
    dx_[i] = val_ * op.t + dx_[i] * xval;
  }
 
  // Compute value
  val_ *= xval;

  return *this;
}

template <typename T, int Num>
template <typename S> 
inline Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >& 
Sacado::ELRFad::Expr< Sacado::ELRFad::SFadExprTag<T,Num> >::
operator /= (const Sacado::ELRFad::Expr<S>& x)
{
  T xval = x.val();

#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SELRFad::operator/=() Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  // Compute partials
  LocalAccumOp< Expr<S> > op(x);

  T xval2 = xval*xval;
    
  // Compute each tangent direction
  for(int i=0; i<Num; ++i) {
    op.t = T(0.);
    op.i = i;
    
    // Automatically unrolled loop that computes
    // for (int j=0; j<N; j++)
    //   op.t += op.partials[j] * x.getTangent<j>(i);
    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
    
    dx_[i] = (dx_[i] * xval - val_ * op.t) / xval2;
  }

  // Compute value
  val_ /= xval;

  return *this;
}

