// Emacs will be in -*- Mode: c++ -*-
// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// ***************** DO NOT REMOVE THIS BANNER *****************
//
// SUMMARY: Templatized Oriented Object Finte Element Method
//          TOOFEM
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
// ORIG-DATE: September 98
// LAST-MOD : 15/09/98
// ************************************************************
#ifndef _tinyvec_h
#define _tinyvec_h

namespace FAD {
  template < class T,  int Num > class TinyVector;
}

#include "utils/tinyveccpy.h"

//------------------------------------------------------------------------------------------------//

namespace FAD {

template < class T,  int Num > class TinyVector {

public:
  typedef T value_type;
// Constructors
  TinyVector() { copy( T() );}
  TinyVector(const T& val) { copy( val );}
  TinyVector(const TinyVector< T, Num >& a) { copy(a); }

// destructor
  ~TinyVector() {;}

// Operators
  inline       T&      operator [] (int i)
    { CheckSize(i); return ptr_to_data[i]; }

  inline const T&      operator [] (int i) const
    { CheckSize(i); return ptr_to_data[i]; }

  inline       T&      operator () ( int i)
    { CheckSize(i); return ptr_to_data[i]; }

  inline const T&      operator () ( int i) const
    { CheckSize(i); return ptr_to_data[i]; }


  TinyVector< T, Num >& operator=(const T & val) { copy(val); return *this;}
  TinyVector< T, Num >& operator=(const TinyVector< T, Num >& a) { copy(a); return *this;}

// Member functions
  inline       int     size()                      const { return capacity();}
  inline       int     capacity()                  const { return Num;}
  inline       T*      begin()                     const { return ptr_to_data;}
  inline       int     no(const T * ptr)           const { return (ptr - begin());}

  inline       TinyVector< T, Num >& operator += (const TinyVector< T, Num >& v);
  inline       TinyVector< T, Num >& operator -= (const TinyVector< T, Num >& v);

private:
  void CheckSize(int i) const
    {
#ifdef CHECK_SIZE
      if ( !( (i >= 0) && (i < Num) ) ) error("TinyVector<>::CheckSize(int i), index out of bound");
#endif
    }


  void    copy( const TinyVector< T, Num >& a )
    { Copy(*this,a); }

  void    copy( const T & val)
    { Copy(*this,val); }

private:
  T     ptr_to_data[Num];
};




template < class T,  int Num > inline
TinyVector< T, Num > operator * (const T& val, const TinyVector< T, Num >& v)
{
  TinyVector< T, Num > tmp(v);

  for (int i=0; i<Num; ++i)
    tmp[i] *= val;

  return tmp;
}


template < class T,  int Num > inline
TinyVector< T, Num >& TinyVector< T, Num >::operator += (const TinyVector< T, Num >& v)
{

  for (int i=0; i<Num; ++i)
    ptr_to_data[i] += v.ptr_to_data[i];

  return *this;
}


template < class T,  int Num > inline
TinyVector< T, Num >& TinyVector< T, Num >::operator -= (const TinyVector< T, Num >& v)
{

  for (int i=0; i<Num; ++i)
    ptr_to_data[i] -= v.ptr_to_data[i];

  return *this;
}

template < class T,  int Num > inline
std::ostream& operator << (std::ostream& os, const TinyVector< T, Num >& v)
{
  os.setf(std::ios::fixed,std::ios::floatfield);

  for (int i=0; i<Num; ++i)
    os << std::setw(12) << v[i];
  os << std::endl;

  return os;
}

} // namespace FAD

#endif
