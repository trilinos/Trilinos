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
#ifndef _vectors_h
#define _vectors_h

#include <cstdio>
#include <cstddef>

#include <memory>
#include <iostream>
#include <iomanip>


#include "utils/error.h"
#include "utils/tinyvec.h"

#define RESTRICT

namespace FAD {

// Copy ala Blitz
template <class T> class MEM_CPY {
 public:
  static void copy(T* RESTRICT dest, const T* src, const int N)
    {
      // Unwind the inner loop, four elements at a time
      int Nmod4 = N & 0x03;

      int i=0;
      for (; i< Nmod4; ++i)
  dest[i] = src[i];

      for (; i<N; i+=4)
        {
    // Common subexpression elimination: avoid doing i+1, i+2, i+3
    // multiple times (compilers *won't* do this cse automatically)
    int i1 = i+1;
    int i2 = i+2;
    int i3 = i+3;

    // Reading all the results first avoids aliasing
    // ambiguities, and provides more flexibility in
    // register allocation & instruction scheduling
    T tmp0,tmp1, tmp2, tmp3;
    tmp0 = src[i];
    tmp1 = src[i1];
    tmp2 = src[i2];
    tmp3 = src[i3];

    dest[i]  = tmp0;
    dest[i1] = tmp1;
    dest[i2] = tmp2;
    dest[i3] = tmp3;
        }
    }
};


//-------------------------------------------------------------------
//                         Tableau 1D
//-------------------------------------------------------------------
//interface

template < class T > class Vector {

public:
  typedef T value_type;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef value_type* iterator;
  typedef const value_type* const_iterator;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef size_t size_type;
  typedef std::ptrdiff_t difference_type;

// Constructors
  inline               Vector();
  inline               Vector(int csize);
                       Vector(int csize, const T& val);
  inline               Vector(const Vector< T >& a);

// destructor
                      ~Vector();

// Operators
  inline       T&      operator [] (int );
  inline const T&      operator [] (int )         const;
  inline       T&      operator () ( int );
  inline const T&      operator () ( int )         const;
               Vector< T >& operator=(const Vector< T >& a);
               Vector< T >& operator=(const T & val);

  Vector< T >& operator+=(const Vector< T >& a)
    {
      int sz = size();
      if ( sz != a.size() ) error("operator-=(const Vector< T >& a), size error");
      for (int i=0; i<sz; ++i) (*this)[i] += a[i];
      return *this;
    }

  Vector< T >& operator-=(const Vector< T >& a)
    {
      int sz = size();
      if ( sz != a.size() ) error("operator-=(const Vector< T >& a), size error");
      for (int i=0; i<sz; ++i) (*this)[i] -= a[i];
      return *this;
    }


// Member functions
  inline       int  empty()  const { return (!size());}
  inline       void reserve(int );
  inline       void resize(int );
  inline       int  size() const { return capacity();}
  inline       int  length() const { return capacity();}
  inline       int  capacity() const;
  inline       T*   begin() const;
  inline       T*   end() const { return ptr_to_data + num_elts;}
  inline       void destroy();
  inline       int  no(const T * ptr) const { return (ptr - begin());}

private:

               void copy(const Vector< T >& a);

private:

  int   num_elts;
  T* RESTRICT  ptr_to_data;
};

template < class T> Vector<T> operator+(const Vector<T> & x, const Vector<T> & y)
{
  int sz = x.size();
  if ( sz!=y.size() ) error("operator+(const Vector<T> & x, const Vector<T> & y), size error");

  Vector<T> z(sz);

  for (int i=0; i<sz; ++i)
    z[i] = x[i] + y[i];

  return z;
}

template < class T> Vector<T> operator-(const Vector<T> & x, const Vector<T> & y)
{
  int sz = x.size();
  if ( sz!=y.size() ) error("operator+(const Vector<T> & x, const Vector<T> & y), size error");

  Vector<T> z(sz);

  for (int i=0; i<sz; ++i)
    z[i] = x[i] - y[i];

  return z;
}

template < class T> Vector<T> operator*(const T& a, const Vector<T> & x)
{
  int sz = x.size();
  Vector<T> y(sz);

  for (int i=0; i<sz; ++i)
    y[i] = a*x[i];

  return y;
}

template <class T> inline std::ostream& operator << (std::ostream& os, const Vector<T>& x)
{
  os.setf(std::ios::fixed,std::ios::floatfield);
  os.width(12);
  int sz = x.size();
  os << sz << std::endl;


  for (int i=0; i< sz; ++i) {
     os.width(12);
     os << x[i] << std::endl;
  }

  os << "\n";

  return os;
}





template< class T > inline Vector< T >::Vector() : num_elts(0), ptr_to_data((T*)0)
// default constructor
{
#ifdef FAD_DEBUG
    cerr << "Vector<>::Vector()               |default constructor|";
#endif

#ifdef FAD_DEBUG
    fprintf(stderr," adr =  %x \n", ptr_to_data);
#endif
}


template< class T > inline Vector< T >::Vector(int csize) : num_elts(0), ptr_to_data((T*)0)
// size constructor
{
#ifdef FAD_DEBUG
    cerr << "Vector<>::Vector(int )          |size constructor   |";
#endif

    if (csize > 0 ){
  num_elts = csize;
  ptr_to_data = new T[num_elts];

  if ( ptr_to_data == 0 ) error("Vector<>::Vector(), Memoire insuffisante");
    }

#ifdef FAD_DEBUG
    fprintf(stderr," adr =  %x \n", ptr_to_data);
#endif
}


template< class T > inline Vector< T >::Vector(const Vector< T >& a) : num_elts(a.num_elts), ptr_to_data((T*)0)
// copy constructor
{
#ifdef FAD_DEBUG
    cerr << "Vector<>::Vector(const Vector< T >& ) |copy constructor   |";
#endif

    if ( num_elts != 0 ) {
      ptr_to_data = new T[num_elts];
      copy(a);
    }

#ifdef FAD_DEBUG
    fprintf(stderr," adr =  %x \n", ptr_to_data);
#endif
}



template< class T > inline void Vector< T >::destroy()
{
#ifdef FAD_DEBUG
    fprintf(stderr,"Vector<>::destroy()         |                   | adr =  %x \n", ptr_to_data);
#endif
    if (ptr_to_data != 0)
  delete [] ptr_to_data;
    else {
    }

    num_elts = 0; ptr_to_data = (T*)0;
}



template< class T >  inline T&   Vector< T >::operator [] (int i)
{
#ifdef FAD_DEBUG
  //    cerr << "T&   Vector< T >::operator [] (int i)" << endl;
#endif

#ifdef CHECK_SIZE
    if ( ptr_to_data == 0 ) error("Vector<>::operator[], empty array");

    if ( !( (i >= 0) && (i < num_elts) ) ) error("Vector<>::operator[], index out of bound");
#endif

    return *(ptr_to_data + i);
}


template< class T >  inline const T&   Vector< T >::operator [] (int i)  const
{
#ifdef FAD_DEBUG
//    cerr << "const T&   Vector< T >::operator [] (int i)  const" << endl;
#endif

#ifdef CHECK_SIZE
    if ( ptr_to_data == 0 ) error("Vector<>::operator[], empty array");

    if ( !( (i >= 0) && (i < num_elts) ) ) error("Vector<>::operator[], index out of bound");
#endif

    return *(ptr_to_data + i);
}

template< class T >  inline T&   Vector< T >::operator () (int i)
{
  return Vector<T>::operator[](i);
}

template< class T >  inline const T&   Vector< T >::operator () (int i)  const
{
  return Vector<T>::operator[](i);
}



template< class T > inline void Vector< T >::reserve(int ssize)
{
#ifdef CHECK_SIZE
    if ( ptr_to_data != 0 ) error("Vector<>::reserve(), array already allocted");

    if ( ssize <0 ) error("Vector<>::reserve(), negative size");
#endif

    if (ssize != 0){
  num_elts = ssize;
  ptr_to_data = new T[num_elts];

  if ( ptr_to_data == 0 ) error("Vector<>::reserve(), not enough memory");
    }
#ifdef FAD_DEBUG
    fprintf(stderr,"Vector<>::reserve()         |                   | adr =  %x \n", ptr_to_data);
#endif
}


template< class T > inline void Vector< T >::resize(int ssize)
{
#ifdef CHECK_SIZE
    if ( ssize <0 ) error("Vector<>::reserve(), negative size");
#endif

    if ( ssize != 0){
      if ( num_elts != 0) destroy();

      num_elts = ssize;
      ptr_to_data = new T[num_elts];

      if ( ptr_to_data == 0 ) error("Vector<>::reserve(), not enough memory");

      for (int i=0; i<num_elts; i++) ptr_to_data[i] = T(0.0);

    }
    else
      if ( num_elts != 0) destroy();

#ifdef FAD_DEBUG
    fprintf(stderr,"Vector<>::reserve()         |                   | adr =  %x \n", ptr_to_data);
#endif
}


template< class T > inline  int Vector< T >::capacity() const
{
    return num_elts;
}



template< class T > inline  T* Vector< T >::begin() const
{
    return ptr_to_data;
}


template< class T > Vector< T >::Vector(int csize, const T& val) : num_elts(0), ptr_to_data((T*)0)
{
#ifdef FAD_DEBUG
    cerr << "Vector<>::Vector(int, const T& )|size constructor   |";
#endif

    if (csize > 0 ){
  num_elts = csize;
  ptr_to_data = new T[num_elts];

  if ( ptr_to_data == 0 ) error("Vector<>::Vector(), Memoire insuffisante");
    }

  T* p = ptr_to_data + num_elts;
  while ( p > ptr_to_data )
    *--p = val;//pb si = n'est pas surcharge pour le type T

#ifdef FAD_DEBUG
    fprintf(stderr," adr =  %x \n", ptr_to_data);
#endif
}


template< class T > void Vector< T >::copy(const Vector< T >& a)
{
  MEM_CPY<T>::copy(ptr_to_data,a.ptr_to_data, num_elts);
}


template< class T > Vector< T >::~Vector()// destructeur
{
#ifdef FAD_DEBUG
    cerr << "Vector<>::~Vector()              |destructor         |" << endl;;
#endif

    destroy();

}


//le return dans le if empeche de mettre la fonction en ligne
template< class T > Vector< T >& Vector< T >::operator = (const Vector< T >& a)
{
  if ( this != &a) {
    // Cas ou le pointeur est non alloue
    if ( num_elts == 0 ) {
      if ( a.num_elts == 0 ) return (*this);
      else reserve(a.num_elts);
    }
    // Cas ou le pointeur est deja alloue
    if ( num_elts != a.num_elts ) {
      //cerr << "size = " << num_elts << ", a.size = " << a.num_elts << endl;
      //cerr << "ptr_to_data = " << (void*)ptr_to_data << ", a.ptr_to_data = " << (void*)a.ptr_to_data << endl;
      //error("Vector::operator=, tailles incompatibles");
      destroy();//on desalloue
      reserve(a.num_elts);
    }

    copy(a);
  }

  return *this;
}


template< class T > Vector< T >& Vector< T >::operator = (const T& val)
{
  if (num_elts != 0)
    for (int i=0; i<num_elts; ++i)
      ptr_to_data[i] = val;
  else error("Vector< T >& Vector< T >::operator = (const T& val), unallocated vector");

  return *this;
}


} // Namespace FAD

#endif
