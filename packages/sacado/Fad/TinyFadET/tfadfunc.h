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
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
#ifndef _tfadfunc_h_
#define _tfadfunc_h_

#include <cmath>

#define FAD_FUNC_MACRO(NAME,FCT,FCTSTD,GRD,FASTGRD)                               \
template <class Expr> class NAME                                          \
{                                                                         \
public:                                                                   \
  typedef typename Expr::value_type value_type;                           \
protected:                                                                \
  NAME () {}                                                              \
                                                                          \
  Expr expr_;                                                             \
public:                                                                   \
    NAME (const Expr & expr) : expr_(expr) {;}                            \
    inline value_type val() const  { return FCTSTD(expr_.val());}         \
    inline value_type dx(int i) const {return GRD ; }                     \
    inline int size() const { return expr_.size();}                       \
                                                                          \
    bool hasFastAccess() const { return expr_.hasFastAccess();}           \
    value_type fastAccessDx(int i) const { return FASTGRD;}               \
};                                                                        \
                                                                          \
template <class Expr> inline TFadExpr< NAME< TFadExpr<Expr> > >             \
FCT (const TFadExpr<Expr>& expr)                                           \
{                                                                         \
    typedef NAME< TFadExpr<Expr> > expr_t;                                 \
    return TFadExpr< expr_t >(  expr_t(expr) );                            \
}                                                                         \
                                                                          \
template <class T,int Num> inline TFadExpr< NAME< TFad<Num,T> > >                       \
FCT (const TFad<Num,T>& x)                                                     \
{                                                                         \
    typedef NAME< TFad<Num,T> > expr_t;                                        \
    return TFadExpr< expr_t >(  expr_t(x) );                               \
}

FAD_FUNC_MACRO(TFadFuncCos,
	       cos,
	       std::cos,
	       -expr_.dx(i)*std::sin( expr_.val() ),
	       -expr_.fastAccessDx(i)*std::sin( expr_.val() ) )
FAD_FUNC_MACRO(TFadFuncSin,
	       sin,
	       std::sin,
	       expr_.dx(i)*std::cos(expr_.val()),
	       expr_.fastAccessDx(i)*std::cos(expr_.val()) )
FAD_FUNC_MACRO(TFadFuncTan,
	       tan,
	       std::tan,
	       expr_.dx(i)*(1.+std::tan(expr_.val())*std::tan(expr_.val())),
	       expr_.fastAccessDx(i)*(1.+std::tan(expr_.val())*std::tan(expr_.val())))
FAD_FUNC_MACRO(TFadFuncAcos,
	       acos,
	       std::acos,
	       -expr_.dx(i)/std::sqrt(1.-expr_.val()*expr_.val()),
	       -expr_.fastAccessDx(i)/std::sqrt(1.-expr_.val()*expr_.val()))
FAD_FUNC_MACRO(TFadFuncAsin,
	       asin,
	       std::asin,
	       expr_.dx(i)/std::sqrt(1.-expr_.val()*expr_.val()),
	       expr_.fastAccessDx(i)/std::sqrt(1.-expr_.val()*expr_.val()))
FAD_FUNC_MACRO(TFadFuncAtan,
	       atan,
	       std::atan,
	       expr_.dx(i)/(1.+expr_.val()*expr_.val()),
	       expr_.fastAccessDx(i)/(1.+expr_.val()*expr_.val()))
FAD_FUNC_MACRO(TFadFuncCosh,
	       cosh,
	       std::cosh,
	       expr_.dx(i)*std::sinh(expr_.val()),
	       expr_.fastAccessDx(i)*std::sinh(expr_.val()))
FAD_FUNC_MACRO(TFadFuncSinh,
	       sinh,
	       std::sinh,
	       expr_.dx(i)*std::cosh(expr_.val()),
	       expr_.fastAccessDx(i)*std::cosh(expr_.val()))
FAD_FUNC_MACRO(TFadFuncTanh,
	       tanh,
	       std::tanh,
	       expr_.dx(i)/(std::cosh(expr_.val())*std::cosh(expr_.val())),
	       expr_.fastAccessDx(i)/(std::cosh(expr_.val())*std::cosh(expr_.val())))
FAD_FUNC_MACRO(TFadFuncAcosh,
	       acish,
	       acosh,
	       expr_.dx(i)/std::sqrt((expr_.val()-1.)/(expr_.val()+1.)),
	       expr_.fastAccessDx(i)/std::sqrt((expr_.val()-1.)/(expr_.val()+1.)))
FAD_FUNC_MACRO(TFadFuncAsinh,
	       asinh,
	       asinh,
	       expr_.dx(i)/std::sqrt(1.+expr_.val()*expr_.val()),
	       expr_.fastAccessDx(i)/std::sqrt(1.+expr_.val()*expr_.val()))
FAD_FUNC_MACRO(TFadFuncAtanh,
	       atanh,
	       atanh,
	       expr_.dx(i)/(1.-expr_.val()*expr_.val()),
	       expr_.fastAccessDx(i)/(1.-expr_.val()*expr_.val()))
FAD_FUNC_MACRO(TFadFuncSqrt,
	       sqrt,
	       std::sqrt,
	       expr_.dx(i)/(2.*std::sqrt(expr_.val())),
	       expr_.fastAccessDx(i)/(2.*std::sqrt(expr_.val())))
// 	       sqrt,
// 	       expr_.dx(i)/(2.*sqrt(expr_.val())),
// 	       expr_.fastAccessDx(i)/(2.*sqrt(expr_.val())))
FAD_FUNC_MACRO(TFadFuncExp,
	       exp,
	       std::exp,
	       expr_.dx(i)*std::exp(expr_.val()),
	       expr_.fastAccessDx(i)*std::exp(expr_.val()))
FAD_FUNC_MACRO(TFadFuncLog,
	       log,
	       std::log,
	       expr_.dx(i)/expr_.val(),
	       expr_.fastAccessDx(i)/expr_.val())
FAD_FUNC_MACRO(TFadFuncLog10,
	       log10,
	       std::log10,
	       expr_.dx(i)/(std::log(value_type(10))*expr_.val()),
	       expr_.fastAccessDx(i)/(std::log(value_type(10))*expr_.val()))
FAD_FUNC_MACRO(TFadFuncAbs,
	       abs,
	       std::abs,
	       expr_.val() >= 0 ? expr_.dx(i) : -expr_.dx(i),
	       expr_.val() >= 0 ? expr_.fastAccessDx(i) : -expr_.fastAccessDx(i))
FAD_FUNC_MACRO(TFadFuncFabs,
	       fabs,
	       std::fabs,
	       expr_.val() >= 0 ? expr_.dx(i) : -expr_.dx(i),
	       expr_.val() >= 0 ? expr_.fastAccessDx(i) : -expr_.fastAccessDx(i))


#undef FAD_FUNC_MACRO



#endif
