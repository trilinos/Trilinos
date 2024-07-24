// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_PROMOTION_TRAITS_HPP_
#define _TEUCHOS_PROMOTION_TRAITS_HPP_

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

template <class A, class B>
class PromotionTraits
{
public:
};

//Specialization
template <class T> class PromotionTraits<T,T> {
public:
  typedef T promote;
};

#define PT_SPEC(type1,type2,type3) \
template <> class PromotionTraits< type1 , type2 > { \
public: \
    typedef type3 promote; \
}; \
template <> class PromotionTraits< type2 , type1 > { \
public: \
    typedef type3 promote; \
};

#ifdef HAVE_TEUCHOS_COMPLEX
PT_SPEC(double,std::complex<float>,std::complex<double>)
PT_SPEC(float,std::complex<double>,std::complex<double>)
PT_SPEC(float,std::complex<float>,std::complex<float>)
PT_SPEC(double,std::complex<double>,std::complex<double>)
#endif // HAVE_TEUCHOS_COMPLEX
PT_SPEC(double,float,double)
PT_SPEC(double,long,double)
PT_SPEC(double,int,double)
PT_SPEC(float,long,float)
PT_SPEC(float,int,float)

// ToDo: Add specializations for extended precision types!

#undef PT_SPEC

} // Teuchos namespace

#endif // _TEUCHOS_PROMOTION_TRAITS_HPP_
