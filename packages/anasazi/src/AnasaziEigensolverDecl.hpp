// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_EIGENSOLVER_DECL_HPP
#define ANASAZI_EIGENSOLVER_DECL_HPP

/*! \file AnasaziEigensolverDecl.hpp
    \brief Forward declaration of the virtual base class Anasazi::Eigensolver.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

namespace Anasazi {

  /*! \class Eigensolver
    \brief The Eigensolver is a templated virtual base class that defines the
     basic interface that any eigensolver will support.
  
     This interface is mainly concerned with providing a set of eigensolver status method that
     can be requested from any eigensolver by an StatusTest object.
  */
  template<class ScalarType, class MV, class OP>
  class Eigensolver;
}

#endif
