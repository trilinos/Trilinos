// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LAPACK_H
#define ROL_LAPACK_H

#include "Teuchos_LAPACK.hpp"
/** \class ROL::LAPACK
  \brief Provides interface to Lapack
  */


namespace ROL { 

template<typename OrdinalType, typename Real>
using LAPACK = Teuchos::LAPACK<OrdinalType, Real>;

}

#endif
