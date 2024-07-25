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
#ifndef _tinyveccpy_h
#define _tinyveccpy_h

namespace FAD {

// Recursive definition
template < int N> struct TinyCopy {
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const TinyVector< T, Num >& x)
    {
      y[N] = x[N];
      TinyCopy<N-1>::eval(y,x);
    }
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const T & val)
    {
      y[N] = val;
      TinyCopy<N-1>::eval(y,val);
    }
};

// Specialization
template <> struct TinyCopy<3> {
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const TinyVector< T, Num >& x)
    {
      y[0] = x[0]; y[1] = x[1]; y[2] = x[2]; y[3] = x[3];
    }
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const T & val)
    {
      y[0] = y[1] = y[2] = y[3] = val;
    }
};

template <> struct TinyCopy<2> {
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const TinyVector< T, Num >& x)
    {
      y[0] = x[0]; y[1] = x[1]; y[2] = x[2]; 
    }
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const T & val)
    {
      y[0] = y[1] = y[2] = val;
    }
};

template <> struct TinyCopy<1> {
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const TinyVector< T, Num >& x)
    {
      y[0] = x[0]; y[1] = x[1]; 
    }
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const T & val)
    {
      y[0] = y[1] = val;
    }
};

template <> struct TinyCopy<0> {
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const TinyVector< T, Num >& x)
    {
      y[0] = x[0];
    }
  template <int Num, class T> static inline void eval( TinyVector< T, Num >& y, const T & val)
    {
      y[0] = val;
    }
};



// Recursive definition
template < int Num, class T> inline
void Copy( TinyVector< T, Num >& y, const TinyVector< T, Num >& x)
{
   for (int i=0; i<Num; ++i)
     y[i] = x[i];
  //  TinyCopy<Num-1>::eval(y,x);
}

template < int Num, class T> inline
void Copy( TinyVector< T, Num >& y, const T & val)
{
   for (int i=0; i<Num; ++i)
     y[i] = val;
   //  TinyCopy<Num-1>::eval(y,val);
}


//Specializaions
template < class T> inline
void Copy( TinyVector< T, 6 >& y, const TinyVector< T, 6 >& x)
{
  y[0] = x[0]; y[1] = x[1]; y[2] = x[2]; 
  y[3] = x[3]; y[4] = x[4]; y[5] = x[5];
}

template < class T> inline
void Copy( TinyVector< T, 6>& y, const T & val)
{
  y[0] = y[1] = y[2] = y[3] = y[4] = y[5] = val;
}


template < class T> inline
void Copy( TinyVector< T, 5 >& y, const TinyVector< T, 5 >& x)
{
  y[0] = x[0]; y[1] = x[1]; y[2] = x[2]; y[3] = x[3]; y[4] = x[4];
}

template < class T> inline
void Copy( TinyVector< T, 5>& y, const T & val)
{
  y[0] = y[1] = y[2] = y[3] = y[4] = val;
}


template < class T> inline
void Copy( TinyVector< T, 4 >& y, const TinyVector< T, 4 >& x)
{
  y[0] = x[0]; y[1] = x[1]; y[2] = x[2]; y[3] = x[3];
}

template < class T> inline
void Copy( TinyVector< T, 4>& y, const T & val)
{
  y[0] = y[1] = y[2] = y[3] = val;
}


template < class T> inline
void Copy( TinyVector< T, 3 >& y, const TinyVector< T, 3 >& x)
{
  y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
}

template < class T> inline
void Copy( TinyVector< T, 3>& y, const T & val)
{
  y[0] = y[1] = y[2] = val;
}


template < class T> inline
void Copy( TinyVector< T, 2 >& y, const TinyVector< T, 2 >& x)
{
  y[0] = x[0]; y[1] = x[1];
}

template < class T> inline
void Copy( TinyVector< T, 2>& y, const T & val)
{
  y[0] = y[1] = val;
}


template < class T> inline
void Copy( TinyVector< T, 1 >& y, const TinyVector< T, 1 >& x)
{
  y[0] = x[0];
}

template < class T> inline
void Copy( TinyVector< T, 1>& y, const T & val)
{
  y[0] = val;
}

} // namespace FAD

#endif
