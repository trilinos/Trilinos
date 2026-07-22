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
// FILE   : tinyfadlog.h
// ********************************************************
#ifndef _tinyfadlog_h_
#define _tinyfadlog_h_



template <int Num, class T> inline
bool operator != (const TinyFad<Num,T>& x, const TinyFad<Num,T>& y) {
  return x.val() != y.val();
}

template <int Num, class T> inline
bool operator != (const T& x, const TinyFad<Num,T>& y) {
  return x != y.val();
}

template <int Num, class T> inline
bool operator != (const TinyFad<Num,T>& x, const T& y) {
  return x.val() != y;
}


template <int Num, class T> inline
bool operator == (const TinyFad<Num,T>& x, const TinyFad<Num,T>& y) {
  return x.val() == y.val();
}

template <int Num, class T> inline
bool operator == (const T& x, const TinyFad<Num,T>& y) {
  return x == y.val();
}

template <int Num, class T> inline
bool operator == (const TinyFad<Num,T>& x, const T& y) {
  return x.val() == y;
}

template <int Num, class T> inline
bool operator > (const TinyFad<Num,T>& x, const TinyFad<Num,T>& y) {
  return x.val() > y.val();
}

template <int Num, class T> inline
bool operator > (const T& x, const TinyFad<Num,T>& y) {
  return x > y.val();
}

template <int Num, class T> inline
bool operator > (const TinyFad<Num,T>& x, const T& y) {
  return x.val() > y;
}


template <int Num, class T> inline
bool operator < (const TinyFad<Num,T>& x, const TinyFad<Num,T>& y) {
  return x.val() < y.val();
}

template <int Num, class T> inline
bool operator < (const T& x, const TinyFad<Num,T>& y) {
  return x < y.val();
}

template <int Num, class T> inline
bool operator < (const TinyFad<Num,T>& x, const T& y) {
  return x.val() < y;
}

template <int Num, class T> inline
bool operator >= (const TinyFad<Num,T>& x, const TinyFad<Num,T>& y) {
  return x.val() >= y.val();
}

template <int Num, class T> inline
bool operator >= (const T& x, const TinyFad<Num,T>& y) {
  return x >= y.val();
}

template <int Num, class T> inline
bool operator >= (const TinyFad<Num,T>& x, const T& y) {
  return x.val() >= y;
}


template <int Num, class T> inline
bool operator <= (const TinyFad<Num,T>& x, const TinyFad<Num,T>& y) {
  return x.val() <= y.val();
}

template <int Num, class T> inline
bool operator <= (const T& x, const TinyFad<Num,T>& y) {
  return x <= y.val();
}

template <int Num, class T> inline
bool operator <= (const TinyFad<Num,T>& x, const T& y) {
  return x.val() <= y;
}


#endif
