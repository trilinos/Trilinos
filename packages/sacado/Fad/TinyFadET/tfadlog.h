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
#ifndef _tfadlog_h_
#define _tfadlog_h_

#define FAD_LOG_MACRO(OP)                                    \
template <int Num,class T> inline bool                       \
operator OP(const TFad<Num,T> &a, const TFad<Num,T> &b)      \
{                                                            \
    return (a.val() OP b.val());                             \
}                                                            \
                                                             \
template <int Num,class T> inline bool                       \
operator OP(const TFad<Num,T> &a, const T &b)                \
{                                                            \
    return (a.val() OP b);                                   \
}                                                            \
                                                             \
template <int Num,class T> inline bool                       \
operator OP(const T &a, const TFad<Num,T> &b)                \
{                                                            \
    return (a OP b.val());                                   \
}                                                            \
                                                             \
template <class T> inline bool                               \
operator OP(const TFadExpr<T> &a, const TFadExpr<T> &b)        \
{                                                            \
    return (a.val() OP b.val());                             \
}                                                            \
                                                             \
template <class T> inline bool                               \
operator OP(const T &a, const TFadExpr<T> &b)                 \
{                                                            \
    return (a OP b.val());                                   \
}                                                            \
                                                             \
template <class T> inline bool                               \
operator OP(const TFadExpr<T> &a, const T &b)                 \
{                                                            \
    return (a.val() OP b);                                   \
}
                                                 

FAD_LOG_MACRO(==)
FAD_LOG_MACRO(!=)
FAD_LOG_MACRO(<)
FAD_LOG_MACRO(>)
FAD_LOG_MACRO(<=)
FAD_LOG_MACRO(>=)
FAD_LOG_MACRO(<<=)
FAD_LOG_MACRO(>>=)
FAD_LOG_MACRO(&)

#undef FAD_LOG_MACRO


template <int Num,class T> inline bool operator !(const TFad<Num,T> &a) {
    return ( !a.val() );
}

template <class T> inline bool operator !(const TFadExpr<T> &a) {
    return ( !a.val() );
}

#endif
