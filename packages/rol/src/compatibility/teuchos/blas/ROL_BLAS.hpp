// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BLAS_H
#define ROL_BLAS_H

#include "Teuchos_BLAS.hpp"
/** \class ROL::BLAS
  \brief Provides interface to BLAS
  */


namespace ROL { 

template<typename OrdinalType, typename Real>
using BLAS = Teuchos::BLAS<OrdinalType, Real>;

}

#endif
