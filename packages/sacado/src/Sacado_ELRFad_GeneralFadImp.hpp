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

// template <typename T, typename Storage> 
// template <typename S> 
// inline Sacado::ELRFad::GeneralFad<T,Storage>::GeneralFad(const Expr<S>& x) :
//   s_(T(0.))
// {
//   int sz = x.size();

//   if (sz != s_.size()) 
//     s_.resize(sz);

//   if (sz) {

//     // Number of arguments
//     const int N = Expr<S>::num_args;

//     // Array to store partials
//     T partials[N];

//     // Array to store tangents
//     T dots[N];

//     // Compute partials
//     x.computePartials(T(1.), partials);
  
//     // Compute each tangent direction
//     for(int i=0; i<sz; ++i) {
//       T t = T(0.);
//       x.getTangents(i,dots);
//       for (int j=0; j<N; j++) {
// 	t += partials[j]*dots[j];
//       }
//       s_.dx_[i] = t;
//     }

//   }
  
//   // Compute value
//   s_.val_ = x.val();
// }

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>::GeneralFad(const Expr<S>& x) :
  s_(T(0.))
{
  int sz = x.size();

  if (sz != s_.size()) 
    s_.resize(sz);

  if (sz) {

    // Number of arguments
    const int N = Expr<S>::num_args;

    // Compute partials
    LocalAccumOp< Expr<S> > op(x);
  
    // Compute each tangent direction
    for(int i=0; i<sz; ++i) {
      op.t = T(0.);
      op.i = i;

      // Automatically unrolled loop that computes
      // for (int j=0; j<N; j++)
      //   op.t += op.partials[j] * x.getTangent<j>(i);
      Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);

      s_.dx_[i] = op.t;
    }

  }
  
  // Compute value
  s_.val_ = x.val();
}

template <typename T, typename Storage> 
inline void 
Sacado::ELRFad::GeneralFad<T,Storage>::diff(const int ith, const int n) 
{ 
  if (s_.size() == 0) 
    s_.resize(n);

  s_.zero();
  s_.dx_[ith] = T(1.);

}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator=(const T& v) 
{
  s_.val_ = v;

  if (s_.size()) {
    s_.zero();    // We need to zero out the array for future resizes
    s_.resize(0);
  }

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator=(
			       const Sacado::ELRFad::GeneralFad<T,Storage>& x) 
{
  // Copy value & dx_
  s_.operator=(x.s_);
  
  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator=(const Expr<S>& x) 
{
  int sz = x.size();

  if (sz != s_.size()) 
    s_.resize(sz);

  if (sz) {

    // Number of arguments
    const int N = Expr<S>::num_args;

    // Compute partials
    LocalAccumOp< Expr<S> > op(x);
  
    // Compute each tangent direction
    for(int i=0; i<sz; ++i) {
      op.t = T(0.);
      op.i = i;

      // Automatically unrolled loop that computes
      // for (int j=0; j<N; j++)
      //   op.t += op.partials[j] * x.getTangent<j>(i);
      Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);

      s_.dx_[i] = op.t;
    }

  }
  
  // Compute value
  s_.val_ = x.val();
  
  return *this;
}

template <typename T, typename Storage> 
inline  Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator += (const T& v)
{
  s_.val_ += v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator -= (const T& v)
{
  s_.val_ -= v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator *= (const T& v)
{
  int sz = s_.size();

  s_.val_ *= v;
  for (int i=0; i<sz; ++i)
    s_.dx_[i] *= v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator /= (const T& v)
{
  int sz = s_.size();

  s_.val_ /= v;
  for (int i=0; i<sz; ++i)
    s_.dx_[i] /= v;

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator += (
					     const Sacado::ELRFad::Expr<S>& x)
{
  int xsz = x.size(), sz = s_.size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  if (xsz) {

    // Compute partials
    LocalAccumOp< Expr<S> > op(x);

    if (sz) {

      // Compute each tangent direction
      for(int i=0; i<xsz; ++i) {
	op.t = T(0.);
	op.i = i;
	
	// Automatically unrolled loop that computes
	// for (int j=0; j<N; j++)
	//   op.t += op.partials[j] * x.getTangent<j>(i);
	Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	
	s_.dx_[i] += op.t;
      }

    }

    else {

      s_.resize(xsz);

      // Compute each tangent direction
      for(int i=0; i<xsz; ++i) {
	op.t = T(0.);
	op.i = i;
	
	// Automatically unrolled loop that computes
	// for (int j=0; j<N; j++)
	//   op.t += op.partials[j] * x.getTangent<j>(i);
	Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	
	s_.dx_[i] = op.t;
      }

    }

  }
  
  // Compute value
  s_.val_ += x.val();

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator -= (
					      const Sacado::ELRFad::Expr<S>& x)
{
  int xsz = x.size(), sz = s_.size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  if (xsz) {

    // Compute partials
    LocalAccumOp< Expr<S> > op(x);

    if (sz) {

      // Compute each tangent direction
      for(int i=0; i<xsz; ++i) {
	op.t = T(0.);
	op.i = i;
	
	// Automatically unrolled loop that computes
	// for (int j=0; j<N; j++)
	//   op.t += op.partials[j] * x.getTangent<j>(i);
	Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	
	s_.dx_[i] -= op.t;
      }

    }

    else {

      s_.resize(xsz);

      // Compute each tangent direction
      for(int i=0; i<xsz; ++i) {
	op.t = T(0.);
	op.i = i;
	
	// Automatically unrolled loop that computes
	// for (int j=0; j<N; j++)
	//   op.t += op.partials[j] * x.getTangent<j>(i);
	Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	
	s_.dx_[i] = -op.t;
      }

    }

  }

  s_.val_ -= x.val();

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator *= (
					     const Sacado::ELRFad::Expr<S>& x)
{
  int xsz = x.size(), sz = s_.size();
  T xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  if (xsz) {

    // Compute partials
    LocalAccumOp< Expr<S> > op(x);

    if (sz) {

      // Compute each tangent direction
      for(int i=0; i<xsz; ++i) {
	op.t = T(0.);
	op.i = i;
	
	// Automatically unrolled loop that computes
	// for (int j=0; j<N; j++)
	//   op.t += op.partials[j] * x.getTangent<j>(i);
	Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	
	s_.dx_[i] = s_.val_ * op.t + s_.dx_[i] * xval;
      }

    }

    else {

      s_.resize(xsz);

      // Compute each tangent direction
      for(int i=0; i<xsz; ++i) {
	op.t = T(0.);
	op.i = i;
	
	// Automatically unrolled loop that computes
	// for (int j=0; j<N; j++)
	//   op.t += op.partials[j] * x.getTangent<j>(i);
	Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	
	s_.dx_[i] = s_.val_ * op.t;
      }

    }

  }

  else {

    if (sz) {
      for (int i=0; i<sz; ++i)
	s_.dx_[i] *= xval;
    }

  }

  s_.val_ *= xval;

  return *this;
}

template <typename T, typename Storage>
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator /= (
					     const Sacado::ELRFad::Expr<S>& x)
{
  int xsz = x.size(), sz = s_.size();
  T xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  // Number of arguments
  const int N = Expr<S>::num_args;

  if (xsz) {

    // Compute partials
    LocalAccumOp< Expr<S> > op(x);

    T xval2 = xval*xval;

    if (sz) {

      // Compute each tangent direction
      for(int i=0; i<xsz; ++i) {
	op.t = T(0.);
	op.i = i;
	
	// Automatically unrolled loop that computes
	// for (int j=0; j<N; j++)
	//   op.t += op.partials[j] * x.getTangent<j>(i);
	Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	
	s_.dx_[i] = (s_.dx_[i] * xval - s_.val_ * op.t) / xval2;
      }

    }

    else {

      s_.resize(xsz);

      // Compute each tangent direction
      for(int i=0; i<xsz; ++i) {
	op.t = T(0.);
	op.i = i;
	
	// Automatically unrolled loop that computes
	// for (int j=0; j<N; j++)
	//   op.t += op.partials[j] * x.getTangent<j>(i);
	Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	
	s_.dx_[i] = -s_.val_ * op.t / xval2;
      }

    }

  }

  else {

    if (sz) {
      for (int i=0; i<sz; ++i)
	s_.dx_[i] /= xval;
    }

  }

  s_.val_ /= xval;

  return *this;
}

