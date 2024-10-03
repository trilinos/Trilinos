// Emacs will be in -*- Mode: c++ -*-
// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

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
//  functions are overloaded ) of 1st-order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
#ifndef _tfad_h_
#define _tfad_h_

// ANSI C++ algorithm include
//#include <algobase.h>
#include <algorithm>

// C include
#include <cmath>
#include <cstring>

// type promotion include
#include "utils/promote.h"

namespace FAD {

template <class L, class R> class NumericalTraits;

template <class T> class TFadExpr;
template <class T> class TFadCst;

template <class T> class TFadUnaryPlus;
template <class T> class TFadUnaryMin;

template <class L, class R> class TFadBinaryAdd;
template <class L, class R> class TFadBinaryMinus;
template <class L, class R> class TFadBinaryMul;
template <class L, class R> class TFadBinaryDiv;

template <int Num, class T=float> class TFad {
public:
  typedef T value_type;

  void copy(const TFad<Num,T>& rhs);
protected:
  T val_;
  T dx_[Num];

public:

  TFad() : val_( T(0.f)) {memset(dx_,0,Num*sizeof(T));}
  TFad(const T & x) : val_(x) {memset(dx_,0,Num*sizeof(T));}
  TFad(const T & x, const int i) : val_(x)
    {memset(dx_,0,Num*sizeof(T));dx_[i]=1.0;}
  TFad(const int sz, const T & x) : val_(x)
    {memset(dx_,0,Num*sizeof(T));}
  TFad(const int sz, const int i, const T & x) : val_(x)
    {memset(dx_,0,Num*sizeof(T));dx_[i]=1.0;}
  TFad(const TFad<Num,T> & x);
  template <class ExprT> TFad(const TFadExpr<ExprT>& fadexpr);
  ~TFad() {}

  void diff(const int ith, const int n = Num);

  const T& val() const { return val_; }
  T& val() { return val_; }

  bool hasFastAccess() const { return true; }

  T& fastAccessDx(int i) { return dx_[i]; }
  const T& fastAccessDx(int i) const { return dx_[i]; }
  const T& d(int i) const { return dx_[i]; }
  const T& dx(int i) const { return dx_[i]; }

  int size() const { return Num;}

  TFad<Num,T> & operator=(const T& val);
  TFad<Num,T> & operator=(const TFad<Num,T>& rhs);
  template <class ExprT> TFad<Num,T> & operator=(const TFadExpr<ExprT>& fadexpr);

  TFadExpr< TFadUnaryPlus< TFad<Num,T> > > operator+ () const;
  TFadExpr< TFadUnaryMin< TFad<Num,T> > > operator- () const;

  TFad<Num,T>& operator+= (const T& x);
  TFad<Num,T>& operator-= (const T& x);
  TFad<Num,T>& operator*= (const T& x);
  TFad<Num,T>& operator/= (const T& x);

  TFad<Num,T>& operator+= (const TFad<Num,T>& x);
  TFad<Num,T>& operator-= (const TFad<Num,T>& x);
  TFad<Num,T>& operator*= (const TFad<Num,T>& x);
  TFad<Num,T>& operator/= (const TFad<Num,T>& x);

  template <class ExprT> TFad<Num,T>& operator*= (const TFadExpr<ExprT>& fadexpr);
  template <class ExprT> TFad<Num,T>& operator/= (const TFadExpr<ExprT>& fadexpr);
  template <class ExprT> TFad<Num,T>& operator+= (const TFadExpr<ExprT>& fadexpr);
  template <class ExprT> TFad<Num,T>& operator-= (const TFadExpr<ExprT>& fadexpr);

};

template <int Num,class T> inline  void TFad<Num,T>::copy(const TFad<Num,T>& rhs)
{
  //dx_ = rhs.dx_;
  for (int i=0; i<Num; ++i)
    dx_[i] = rhs.dx_[i];
  val_ = rhs.val_;
}

template <int Num,class T> inline  TFad<Num,T>::TFad(const TFad<Num,T> & rhs)
{
  copy(rhs);
}

template <int Num,class T> template <class ExprT> inline
TFad<Num,T>::TFad(const TFadExpr<ExprT>& fadexpr)
: val_(fadexpr.val())
{

  for(int i=0; i<Num; ++i)
    dx_[i] = fadexpr.dx(i);

}

template <int Num,class T> inline  void TFad<Num,T>::diff(const int ith, const int n)
{

  for(int i=0; i<Num; ++i)
    dx_[i] = T(0);
  dx_[ith] = T(1.);

}

template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator=(const T& v)
{
  val_ = v;

  for(int i=0; i<Num; ++i)
    dx_[i] = T(0);

  return *this;
}

template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator=(const TFad<Num,T>& rhs)
{
  if ( this != &rhs ) copy(rhs);

  return *this;
}

template <int Num,class T> template <class ExprT> inline TFad<Num,T> & TFad<Num,T>::operator=(const TFadExpr<ExprT>& fadexpr)
{
  for(int i=0; i<Num; ++i)
    dx_[i] = fadexpr.dx(i);

  val_ = fadexpr.val();

  return *this;
}

template <int Num,class T> inline TFadExpr< TFadUnaryPlus< TFad<Num,T> > >
TFad<Num,T>::operator+ () const
{
  return TFadExpr< TFadUnaryPlus< TFad<Num,T> > >(*this);
}

template <int Num,class T> inline TFadExpr< TFadUnaryMin< TFad<Num,T> > >
TFad<Num,T>::operator- () const
{
  return TFadExpr< TFadUnaryMin< TFad<Num,T> > >(*this);
}


template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator+= (const T& v)
{
  val_ += v;

  return *this;
}

template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator-= (const T& v)
{
  val_ -= v;

  return *this;
}

template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator*= (const T& v)
{
  val_ *= v;

  for (int i=0; i<Num;++i)
      dx_[i] *= v;

  return *this;
}

template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator/= (const T& v)
{
  val_ /= v;

  for (int i=0; i<Num;++i)
      dx_[i] /= v;

  return *this;
}


template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator+= (const TFad<Num,T>& x)
{
  for (int i=0; i<Num; ++i)
    dx_[i] += x.dx_[i];

  val_ += x.val_;

  return *this;
}

template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator-= (const TFad<Num,T>& x)
{
  for (int i=0; i<Num; ++i)
    dx_[i] -= x.dx_[i];

  val_ -= x.val_;

  return *this;
}

template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator*= (const TFad<Num,T>& x)
{
  for (int i=0; i<Num; ++i)
    dx_[i] = val_ * x.dx_[i] + dx_[i] * x.val_;

  val_ *= x.val_;

  return *this;
}

template <int Num,class T> inline  TFad<Num,T> & TFad<Num,T>::operator/= (const TFad<Num,T>& x)
{
  for (int i=0; i<Num; ++i)
    dx_[i] = (dx_[i]*x.val_ - val_*x.dx_[i])/ (x.val_*x.val_);

  val_ /= x.val_;

  return *this;
}

template <int Num,class T> template <class ExprT> inline  TFad<Num,T> & TFad<Num,T>::operator+= (const TFadExpr<ExprT>& x)
{
  for (int i=0; i<Num; ++i)
    dx_[i] += x.dx(i);

  val_ += x.val();

  return *this;
}

template <int Num,class T> template <class ExprT> inline  TFad<Num,T> & TFad<Num,T>::operator-= (const TFadExpr<ExprT>& x)
{
  for (int i=0; i<Num; ++i)
    dx_[i] -= x.dx(i);
  val_ -= x.val();


  return *this;
}

template <int Num,class T> template <class ExprT> inline  TFad<Num,T> & TFad<Num,T>::operator*= (const TFadExpr<ExprT>& x)
{
  T xval = x.val();
  for (int i=0; i<Num; ++i)
    dx_[i] = val_ * x.dx(i) + dx_[i] * xval;

  val_ *= xval;

  return *this;
}

template <int Num,class T> template <class ExprT> inline  TFad<Num,T> & TFad<Num,T>::operator/= (const TFadExpr<ExprT>& x)
{
  T xval = x.val();

  for (int i=0; i<Num; ++i)
    dx_[i] = ( dx_[i]*xval - val_*x.dx(i) )/ (xval*xval);

  val_ /= xval;

  return *this;
}

//------------------------------- TFad ostream operator ------------------------------------------
// template <int Num,class T> inline ostream& operator << (ostream& os, const TFad<Num,T>& a)
// {
//   os.setf(ios::fixed,ios::floatfield);
//   os.width(12);
//   os << a.val() << "  [";


//   for (int i=0; i< Num; i++) {
//      os.width(12);
//      os << a.dx(i);
//   }

//   os << "]\n";
//   return os;
// }

template <int Num,class T> inline std::ostream& operator << (std::ostream& os, const TFad<Num,T>& a)
{
  os << a.val();
  return os;
}

//------------------------------- TFad expression ------------------------------------------
template < class T > class TFadExpr {
public:
  typedef typename T::value_type value_type;

protected:
  TFadExpr() {}

  T fadexpr_;

public:
  explicit TFadExpr(const T& fadexpr) : fadexpr_(fadexpr) {;}

  value_type val()     const { return fadexpr_.val();}
  value_type dx(int i) const { return fadexpr_.dx(i);}
  int size() const {return fadexpr_.size();}

  bool hasFastAccess() const { return fadexpr_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return fadexpr_.fastAccessDx(i);}
};

//------------------------------- TFad constant ------------------------------------------
template < class T > class TFadCst {
public:
  typedef T value_type;

protected:
  TFadCst() {}

  const T constant_;

public:
  explicit TFadCst(const T& value) : constant_(value) {;}

  const value_type& val()     const { return constant_;}
  const value_type dx(int i) const { return value_type(0);}
  int size() const {return 0;}

  bool hasFastAccess() const { return 1;}
  value_type& fastAccessDx(int i) const { return value_type(0);}
};

//------------------------------- TFad unary + ------------------------------------------
template < class T > class TFadUnaryPlus {
public:
  typedef typename T::value_type value_type;

protected:

  const T& expr_;

public:
  TFadUnaryPlus(const T& value) : expr_(value) {;}

  const value_type val()     const { return expr_.val();}
  const value_type dx(int i) const { return expr_.dx(i);}
  int size() const {return expr_.size();}

  bool hasFastAccess() const { return expr_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return expr_.fastAccessDx(i);}
};

//------------------------------- TFad unary - ------------------------------------------
template < class T > class TFadUnaryMin {
public:
  typedef typename T::value_type value_type;

protected:

  const T& expr_;

public:
  TFadUnaryMin(const T& value) : expr_(value) {;}

  const value_type val()     const { return - expr_.val();}
  const value_type dx(int i) const { return - expr_.dx(i);}
  int size() const {return expr_.size();}

  bool hasFastAccess() const { return expr_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return - expr_.fastAccessDx(i);}
};

template <class T> inline
TFadExpr< TFadUnaryPlus< TFadExpr<T> > >
operator + (const TFadExpr<T>& expr)
{
  typedef TFadUnaryPlus< TFadExpr<T> > expr_t;

  return TFadExpr< expr_t >( expr_t(expr) );
}

template <class T> inline
TFadExpr< TFadUnaryMin< TFadExpr<T> > >
operator - (const TFadExpr<T>& expr)
{
  typedef TFadUnaryMin< TFadExpr<T> > expr_t;

  return TFadExpr< expr_t >( expr_t(expr) );
}

template <class T> inline std::ostream& operator << (std::ostream& os, const TFadExpr<T>& a)
{
  os << a.val();
  return os;
}

#include "TinyFadET/tfadlog.h"
#include "TinyFadET/tfadop.h"
#include "TinyFadET/tfadfunc.h"

} // namespace FAD

#endif
