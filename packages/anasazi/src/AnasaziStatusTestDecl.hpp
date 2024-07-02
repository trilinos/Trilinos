// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//                 Belos: Block Linear Solvers Package
//

#ifndef ANASAZI_STATUS_TEST_DECL_HPP
#define ANASAZI_STATUS_TEST_DECL_HPP

/*!
  \file AnasaziStatusTestDecl.hpp
  \brief Forward declaration of pure virtual base class Anasazi::StatusTest.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

namespace Anasazi {

  /// \class StatusTest
  /// \brief Common interface of stopping criteria for Anasazi's solvers.
  ///
  /// StatusTest is an interface that can be implemented to create
  /// convergence tests for all Anasazi solvers. Almost any kind of test
  /// can be expressed using this mechanism, including composite tests
  /// (see StatusTestCombo).
  template <class ScalarType, class MV, class OP>
  class StatusTest;
}


#endif
