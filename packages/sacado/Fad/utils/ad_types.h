// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef AD_TYPES_H
#define AD_TYPES_H

#include "Fad/fad.h"
#include "TinyFadET/tfad.h"
#include "TinyFad/tinyfad.h"

namespace FAD {

  template <typename T>
  class ADTypes {
  public:
  };

  template <>
  class ADTypes<double> {
  public:
    typedef double ConstADType;
    typedef double IndepADType;
  };

  template <typename T>
  class ADTypes< Fad<T> > {
  public:
    typedef Fad<T> ConstADType;
    typedef Fad<T> IndepADType;
  };

  template <int Num, typename T>
  class ADTypes< TFad<Num,T> > {
  public:
    typedef TFad<Num,T> ConstADType;
    typedef TFad<Num,T> IndepADType;
  };

  template <int Num, typename T>
  class ADTypes< TinyFad<Num,T> > {
  public:
    typedef TinyFad<Num,T> ConstADType;
    typedef TinyFad<Num,T> IndepADType;
  };

}

#endif
