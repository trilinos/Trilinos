// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "fad_expr_funcs.hpp"

using std::sin;

template <typename T>
void ExprFuncs::mult<T,1>::operator()(const T x[], T& y) const {
  y = 
    x[0]*x[1];
}

template <typename T>
void ExprFuncs::mult<T,2>::operator()(const T x[], T& y) const {
  y = 
    x[0]*x[1]*x[2];
}

template <typename T>
void ExprFuncs::mult<T,3>::operator()(const T x[], T& y) const {
  y = 
    x[0]*x[1]*x[2]*x[3];
}

template <typename T>
void ExprFuncs::mult<T,4>::operator()(const T x[], T& y) const {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4];
}

template <typename T>
void ExprFuncs::mult<T,5>::operator()(const T x[], T& y) const {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4]*x[5];
}

template <typename T>
void ExprFuncs::mult<T,10>::operator()(const T x[], T& y) const {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*
    x[8]*x[9]*x[10];
}

template <typename T>
void ExprFuncs::mult<T,15>::operator()(const T x[], T& y) const {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*
    x[8]*x[9]*x[10]*x[11]*x[12]*x[13]*x[14]*x[15];
}

template <typename T>
void ExprFuncs::mult<T,20>::operator()(const T x[], T& y) const {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*
    x[8]*x[9]*x[10]*x[11]*x[12]*x[13]*x[14]*x[15]*
    x[16]*x[17]*x[18]*x[19]* x[20];
}

template <typename T>
void ExprFuncs::nest<T,1>::operator()(const T x[], T& y) const {
  y = sin(x[0]);
}

template <typename T>
void ExprFuncs::nest<T,2>::operator()(const T x[], T& y) const {
  y = sin(sin(x[0]));
}

template <typename T>
void ExprFuncs::nest<T,3>::operator()(const T x[], T& y) const {
  y = sin(sin(sin(x[0])));
}

template <typename T>
void ExprFuncs::nest<T,4>::operator()(const T x[], T& y) const {
  y = sin(sin(sin(sin(x[0]))));
}

template <typename T>
void ExprFuncs::nest<T,5>::operator()(const T x[], T& y) const {
  y = sin(sin(sin(sin(sin(x[0])))));
}

template <typename T>
void ExprFuncs::nest<T,10>::operator()(const T x[], T& y) const {
  y = sin(sin(sin(sin(sin(
    sin(sin(sin(sin(sin(x[0]))))))))));
}

template <typename T>
void ExprFuncs::nest<T,15>::operator()(const T x[], T& y) const {
  y = sin(sin(sin(sin(sin(
    sin(sin(sin(sin(sin(
    sin(sin(sin(sin(sin(x[0])))))))))))))));
}

template <typename T>
void ExprFuncs::nest<T,20>::operator()(const T x[], T& y) const {
  y = sin(sin(sin(sin(sin(
    sin(sin(sin(sin(sin(
    sin(sin(sin(sin(sin(
    sin(sin(sin(sin(sin(x[0]))))))))))))))))))));
}

template <typename T>
void ExprFuncs::add<T,1>::operator()(const T x[], T& y) const {
  y = 
    x[0]+x[1];
}

template <typename T>
void ExprFuncs::add<T,2>::operator()(const T x[], T& y) const {
  y = 
    x[0]+x[1]+x[2];
}

template <typename T>
void ExprFuncs::add<T,3>::operator()(const T x[], T& y) const {
  y = 
    x[0]+x[1]+x[2]+x[3];
}

template <typename T>
void ExprFuncs::add<T,4>::operator()(const T x[], T& y) const {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4];
}

template <typename T>
void ExprFuncs::add<T,5>::operator()(const T x[], T& y) const {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4]+x[5];
}

template <typename T>
void ExprFuncs::add<T,10>::operator()(const T x[], T& y) const {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+
    x[8]+x[9]+x[10];
}

template <typename T>
void ExprFuncs::add<T,15>::operator()(const T x[], T& y) const {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+
    x[8]+x[9]+x[10]+x[11]+x[12]+x[13]+x[14]+x[15];
}

template <typename T>
void ExprFuncs::add<T,20>::operator()(const T x[], T& y) const {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+
    x[8]+x[9]+x[10]+x[11]+x[12]+x[13]+x[14]+x[15]+
    x[16]+x[17]+x[18]+x[19]+x[20];
}

const char* ExprFuncs::mult_names[ExprFuncs::nfunc] = { 
  "mult1",
  "mult2",
  "mult3",
  "mult4",
  "mult5",
  "mult10",
  "mult15",
  "mult20"
};

const char* ExprFuncs::nest_names[ExprFuncs::nfunc] = { 
  "nest1",
  "nest2",
  "nest3",
  "nest4",
  "nest5",
  "nest10",
  "nest15",
  "nest20"
};

const char* ExprFuncs::add_names[ExprFuncs::nfunc] = { 
  "add1",
  "add2",
  "add3",
  "add4",
  "add5",
  "add10",
  "add15",
  "add20"
};

#define INSTANTIATE_FUNCS(TYPE) \
  template struct ExprFuncs::mult< TYPE,1>;	\
  template struct ExprFuncs::mult< TYPE,2>;	\
  template struct ExprFuncs::mult< TYPE,3>;	\
  template struct ExprFuncs::mult< TYPE,4>;	\
  template struct ExprFuncs::mult< TYPE,5>;	\
  template struct ExprFuncs::mult< TYPE,10>;	\
  template struct ExprFuncs::mult< TYPE,15>;	\
  template struct ExprFuncs::mult< TYPE,20>;	\
						\
  template struct ExprFuncs::nest< TYPE,1>;	\
  template struct ExprFuncs::nest< TYPE,2>;	\
  template struct ExprFuncs::nest< TYPE,3>;	\
  template struct ExprFuncs::nest< TYPE,4>;	\
  template struct ExprFuncs::nest< TYPE,5>;	\
  template struct ExprFuncs::nest< TYPE,10>;	\
  template struct ExprFuncs::nest< TYPE,15>;	\
  template struct ExprFuncs::nest< TYPE,20>;	\
						\
  template struct ExprFuncs::add< TYPE,1>;	\
  template struct ExprFuncs::add< TYPE,2>;	\
  template struct ExprFuncs::add< TYPE,3>;	\
  template struct ExprFuncs::add< TYPE,4>;	\
  template struct ExprFuncs::add< TYPE,5>;	\
  template struct ExprFuncs::add< TYPE,10>;	\
  template struct ExprFuncs::add< TYPE,15>;	\
  template struct ExprFuncs::add< TYPE,20>;

INSTANTIATE_FUNCS(double)
INSTANTIATE_FUNCS(Sacado::Fad::DFad<double>)
INSTANTIATE_FUNCS(Sacado::ELRFad::DFad<double>)
INSTANTIATE_FUNCS(Sacado::CacheFad::DFad<double>)
INSTANTIATE_FUNCS(Sacado::ELRCacheFad::DFad<double>)
INSTANTIATE_FUNCS(Sacado::Fad::SimpleFad<double>)
#ifdef HAVE_ADOLC
INSTANTIATE_FUNCS(adouble)
#endif

#undef INSTANTIATE_FUNCS
