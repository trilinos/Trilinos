// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef FAD_EXPR_FUNCS_HPP
#define FAD_EXPR_FUNCS_HPP

#include "Sacado.hpp"
#include "Sacado_Fad_SimpleFad.hpp"

// ADOL-C includes
#ifdef HAVE_ADOLC
#ifdef PACKAGE
#undef PACKAGE
#endif
#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif
#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif
#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif
#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif
#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif
#ifdef VERSION
#undef VERSION
#endif
//#define ADOLC_TAPELESS
#define NUMBER_DIRECTIONS 100
#include "adolc/adouble.h"
#include "adolc/drivers/drivers.h"
#include "adolc/interfaces.h"
#include "adolc/taping.h"
#endif

struct ExprFuncs {
  static const int nfunc = 8;
  static const char* mult_names[nfunc];
  static const char* nest_names[nfunc];
  static const char* add_names[nfunc];
  static const int nx_max = 21;
  
  template <typename T, int N> struct mult {};
  template <typename T, int N> struct mult_base { static const int n = N+1; };

  template <typename T, int N> struct add {};
  template <typename T, int N> struct add_base { static const int n = N+1; };

  template <typename T, int N> struct nest {};
  template <typename T, int N> struct nest_base { static const int n = 1; };
};

template <typename T> struct ExprFuncs::mult<T,1> : public mult_base<T,1> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::mult<T,2> : public mult_base<T,2> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::mult<T,3> : public mult_base<T,3> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::mult<T,4> : public mult_base<T,4> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::mult<T,5> : public mult_base<T,5> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::mult<T,10> : public mult_base<T,10> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::mult<T,15> : public mult_base<T,15> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::mult<T,20> : public mult_base<T,20> { 
  void operator()(const T x[], T& y) const; };

template <typename T> struct ExprFuncs::add<T,1> : public add_base<T,1> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::add<T,2> : public add_base<T,2> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::add<T,3> : public add_base<T,3> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::add<T,4> : public add_base<T,4> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::add<T,5> : public add_base<T,5> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::add<T,10> : public add_base<T,10> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::add<T,15> : public add_base<T,15> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::add<T,20> : public add_base<T,20> { 
  void operator()(const T x[], T& y) const; };

  
template <typename T> struct ExprFuncs::nest<T,1> : public nest_base<T,1> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::nest<T,2> : public nest_base<T,2> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::nest<T,3> : public nest_base<T,3> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::nest<T,4> : public nest_base<T,4> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::nest<T,5> : public nest_base<T,5> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::nest<T,10> : public nest_base<T,10> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::nest<T,15> : public nest_base<T,15> { 
  void operator()(const T x[], T& y) const; };
template <typename T> struct ExprFuncs::nest<T,20> : public nest_base<T,20> { 
  void operator()(const T x[], T& y) const; };

#endif
