// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_LINEARALGEBRA_HPP
#define ROL_LINEARALGEBRA_HPP


#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

/** \file  ROL_LinearAlgebra.hpp
    \brief Provides basic capabilities in solving dense
           linear systems and eigenvalue problems using
           Eigen to provide the implementation */


namespace ROL {

namespace LA {

template<typename Real>
using Vector = Teuchos::SerialDenseVector<int, Real>;

template<typename Real>
using Matrix = Teuchos::SerialDenseMatrix<int, Real>;

using Teuchos::ETransp;
using Teuchos::NO_TRANS;
using Teuchos::CONJ_TRANS;
using Teuchos::DataAccess;
using Teuchos::Copy;
using Teuchos::View;

} // namespace LA

} // namespace ROL

#endif // ROL_LINEARAPGEBRA_HPP
