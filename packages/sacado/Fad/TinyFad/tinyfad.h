// Emacs will be in -*- Mode: c++ -*-
// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Tools for Automatic Differentiaton (order 1)
// RELEASE: 0.1
// USAGE  : You may copy freely these files and use it for
//          teaching or research. These or part of these may
//          not be sold or used for a commercial purpose with-
//          out our consent : fax (33)1 44 27 72 00
//
// AUTHOR : Nicolas Di cesare
// ORG    :
// E-MAIL : Nicolas.Dicesare@ann.jussieu.fr
//
// ORIG-DATE: September 97
// LAST-MOD : 28/07/98
// ********************************************************
// FILE   : tinyfad.h
// ********************************************************
#ifndef _tinyfad_h_
#define _tinyfad_h_

// C include
#include <cmath>

// ANSI C++ include
#include <iostream>

// include to manage error
#include "utils/error.h"

namespace FAD {

template <class A, class B> class NumericalTraits;

class No_Initialization {
public:
  No_Initialization() {;}
  ~No_Initialization() {;}
};



template <int Num, class T = float> class TinyFad {
protected:

  int n;        // number of the control or zero for independant variables
  T val_;       // value

  T df_[Num];   // vector of partial derivatives


public:
  typedef T value_type;

  // controls constructor
  TinyFad(const T& ind, const int ini);

  // expressions constructors
  TinyFad();
  TinyFad(const No_Initialization &): n(0) {;}
  TinyFad(const T& in);
  TinyFad(const TinyFad& in);

  // destructor
  ~TinyFad() ;

  void diff(const int ith, const int n = Num);

  // acces functions
  int    N()       const {return n-1;}

  const T& val()     const { return val_;}
  T& val()                 { return val_;}


  const T& d(int i) const { return df_[i];}

  const T& dx(int i) const;
  T& dx(int i);


  // operators
  TinyFad<Num,T> & operator = (const TinyFad<Num,T> & in);
  TinyFad<Num,T> & operator = (const T & in);

  TinyFad<Num,T> & operator += (const TinyFad<Num,T> & in);
  TinyFad<Num,T> & operator -= (const TinyFad<Num,T> & in);
  TinyFad<Num,T> & operator *= (const TinyFad<Num,T> & in);
  TinyFad<Num,T> & operator /= (const TinyFad<Num,T> & in);

  TinyFad<Num,T> & operator += (const T & in);
  TinyFad<Num,T> & operator -= (const T & in);
  TinyFad<Num,T> & operator *= (const T & in);
  TinyFad<Num,T> & operator /= (const T & in);

  TinyFad<Num,T> operator++(int);
  TinyFad<Num,T> operator--(int);
  TinyFad<Num,T> & operator++();
  TinyFad<Num,T> & operator--();

};


template <int Num, class T> inline
TinyFad<Num,T>::TinyFad(const T& ind, const int ini) : n(ini+1), val_(ind)
{
#ifdef DEBUG_TINYFAD
  cout << "TinyFad::TinyFad(const T& ind, const int ini), ind = " << ind << ", ini = " << ini << endl;
#endif
  if ( ini >= 2 ) error("control number ini, out of bound ");

  for (int i=0; i<Num; ++i)
    df_[i] = 0.;

  df_[ini] = 1.;
}

template <int Num, class T> inline
TinyFad<Num,T>::TinyFad() : n(0), val_(0.)
{
#ifdef DEBUG_TINYFAD
  cout << "TinyFad::TinyFad()" << endl;
#endif
  for (int i=0; i<Num; ++i)
    df_[i] = 0.;
}

template <int Num, class T> inline
TinyFad<Num,T>::TinyFad(const T& in) : n(0), val_(in)
{
#ifdef DEBUG_TINYFAD
  cout << "TinyFad::TinyFad(const T& in)"  << endl;
#endif
  for (int i=0; i<Num; ++i)
    df_[i] = 0.;
}

template <int Num, class T> inline
TinyFad<Num,T>::TinyFad(const TinyFad<Num,T>& in) : n(0), val_(in.val_)
{
#ifdef DEBUG_TINYFAD
  cout << "TinyFad::TinyFad(const TinyFad& in)"  << endl;
#endif

  for (int i=0; i<Num; ++i)
    df_[i] = in.df_[i];
}

template <int Num, class T> inline
TinyFad<Num,T>::~TinyFad()
{
#ifdef DEBUG_TINYFAD
  cout << "TinyFad::~TinyFad()"  << endl;
#endif
}

template <int Num, class T> inline
void TinyFad<Num,T>::diff(const int ith, const int n_)
{

  df_ = T(0.);
  df_[ith] = T(1.);

}

template <int Num, class T>
const T& TinyFad<Num,T>::dx(int i) const {
  if ( (i<0) || (i>=Num ) ) {
    std::cerr << "i = " << i << std::endl;
    error("df_, partial derivative undefined");
  }

  return df_[i];
}

template <int Num, class T>
T& TinyFad<Num,T>::dx(int i) {
  if ( (i<0) || (i>=Num ) ) {
    std::cerr << "i = " << i << std::endl;
    error("df_, partial derivative undefined");
  }

  return df_[i];
}



template <int Num, class T> inline
TinyFad<Num,T> & TinyFad<Num,T>::operator = (const TinyFad<Num,T> & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad & TinyFad::operator = (const TinyFad & in), you do not change the value of control");
#endif

  val_ = in.val_;

  for (int i=0; i<Num; ++i)
    df_[i] = in.df_[i];

  return *this;
}



template <int Num, class T> TinyFad<Num,T> & TinyFad<Num,T>::operator = (const T & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad & TinyFad::operator = (const T & in), you do not change the value of control");
#endif
  val_ = in;

  for (int i=0; i<Num; i++)
      df_[i] = 0.;

  return *this;
}


template <int Num, class T> TinyFad<Num,T> & TinyFad<Num,T>::operator += (const TinyFad<Num,T> & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad & TinyFad::operator += (const TinyFad & in), you do not change the value of control");
#endif

  for (int i=0; i< Num; i++)
      df_[i] += in.df_[i] ;

  val_ += in.val_;

  return *this;
}

template <int Num, class T> TinyFad<Num,T> & TinyFad<Num,T>::operator -= (const TinyFad<Num,T> & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad & TinyFad::operator -= (const TinyFad & in), you do not change the value of control");
#endif

  for (int i=0; i< Num; i++)
      df_[i] -= in.df_[i] ;

  val_ -= in.val_;

  return *this;
}

template <int Num, class T> TinyFad<Num,T> & TinyFad<Num,T>::operator *= (const TinyFad<Num,T> & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad & TinyFad::operator *= (const TinyFad & in), you do not change the value of control");
#endif

  for (int i=0; i< Num; i++)
	df_[i] = df_[i] * in.val_ + val_ * in.df_[i];

  val_ *= in.val_;

  return *this;
}

template <int Num, class T> TinyFad<Num,T> & TinyFad<Num,T>::operator /= (const TinyFad<Num,T> & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad & TinyFad::operator /= (const TinyFad & in), you do not change the value of control");
#endif

  if (in.val_ == 0.) error("TinyFad & TinyFad::operator /= (const TinyFad & in), dividing by 0");

  for (int i=0; i< Num; i++)
      df_[i] = ( df_[i]* in.val_ - val_ * in.df_[i] ) / in.val_ / in.val_ ;

  val_ /= in.val_;

  return *this;
}

template <int Num, class T> TinyFad<Num,T> & TinyFad<Num,T>::operator += (const T & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad<Num,T> & TinyFad<Num,T>::operator += (const T & in), you do not change the value of control");
#endif

  val_ += in;

  return *this;
}

template <int Num, class T> TinyFad<Num,T> & TinyFad<Num,T>::operator -= (const T & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad & TinyFad::operator -= (const T & in), you do not change the value of control");
#endif

  val_ -= in;

  return *this;
}

template <int Num, class T> TinyFad<Num,T> & TinyFad<Num,T>::operator *= (const T & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad & TinyFad::operator *= (const T & in), you do not change the value of control");
#endif

  val_ *= in;

  for (int i=0; i< Num; i++)
         df_[i] *= in;

  return *this;
}

template <int Num, class T> TinyFad<Num,T> & TinyFad<Num,T>::operator /= (const T & in)
{
#ifdef CHECK_VAR
  if (n) error("TinyFad & TinyFad::operator /= (const T & in), you do not change the value of control");
#endif

  if (in == 0.) error("TinyFad & TinyFad::operator /= (const T & in), dividing by 0");

  val_ /= in;

  for (int i=0; i< Num; i++)
      df_[i] /= in;


  return *this;
}


template <int Num, class T> inline
TinyFad<Num,T> TinyFad<Num,T>::operator++(int)
{
  TinyFad<Num,T> tmp(*this);
  tmp.val_++;
  return tmp;
}

template <int Num, class T> inline
TinyFad<Num,T> TinyFad<Num,T>::operator--(int)
{
  TinyFad<Num,T> tmp(*this);
  tmp.val_--;
  return tmp;
}

template <int Num, class T> inline
TinyFad<Num,T> & TinyFad<Num,T>::operator++()
{
  ++val_;
  return *this;
}

template <int Num, class T> inline
TinyFad<Num,T> & TinyFad<Num,T>::operator--()
{
  --val_;
  return *this;
}


//----------------------------------------------------------------------------------------------//
// unary operators
template <int Num, class T> inline TinyFad<Num,T> operator + (const TinyFad<Num,T>& in)
{
  return TinyFad<Num,T>(in);
}

template <int Num, class T> inline TinyFad<Num,T> operator - (const TinyFad<Num,T>& in)
{
  TinyFad<Num,T> tmp;
  tmp -= in;

  return tmp;
}


template <int Num, class T> std::ostream& operator << (std::ostream& os, const TinyFad<Num,T>& a)
{
  os.setf(std::ios::fixed,std::ios::floatfield);
  os.width(12);
  os << a.val() << "  [";

  for (int i=0; i< Num; i++) {
     os.width(12);
     os << a.dx(i);
  }

  os << "]\n";
  return os;
}



// logical operators
#include "TinyFad/tinyfadlog.h"
// binary operators
#include "TinyFad/tinyfadbin.h"
// math functions
#include "TinyFad/tinyfadfunc.h"

// specializations
#include "TinyFad/Specializations/tinyfadone.h"
#include "TinyFad/Specializations/tinyfadtwo.h"
#include "TinyFad/Specializations/tinyfadthree.h"
#include "TinyFad/Specializations/tinyfadfour.h"
#include "TinyFad/Specializations/tinyfadfive.h"
#include "TinyFad/Specializations/tinyfadsix.h"
#include "TinyFad/Specializations/tinyfadseven.h"
#include "TinyFad/Specializations/tinyfadeight.h"
#include "TinyFad/Specializations/tinyfadnine.h"
#include "TinyFad/Specializations/tinyfadten.h"
#include "TinyFad/Specializations/tinyfadeleven.h"
#include "TinyFad/Specializations/tinyfadtwelve.h"
#include "TinyFad/Specializations/tinyfadthirteen.h"
#include "TinyFad/Specializations/tinyfadfourteen.h"
#include "TinyFad/Specializations/tinyfadfifteen.h"
#include "TinyFad/Specializations/tinyfadsixteen.h"
#include "TinyFad/Specializations/tinyfadseventeen.h"
#include "TinyFad/Specializations/tinyfadeighteen.h"
#include "TinyFad/Specializations/tinyfadnineteen.h"
#include "TinyFad/Specializations/tinyfadtwenty.h"

} // namespace FAD

#endif
