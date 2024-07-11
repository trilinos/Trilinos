// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_LOCA_OPERATORTRAITS_H
#define ANASAZI_LOCA_OPERATORTRAITS_H

// Anasazi includes
#include "AnasaziOperatorTraits.hpp"

// NOX/LOCA includes
#include "NOX_Abstract_MultiVector.H"
#include "LOCA_AnasaziOperator_AbstractStrategy.H"

namespace Anasazi {

  /*!
   * \brief Implementation of the Anasazi::OperatorTraits for LOCA operators.
   */
  /*!
   * This class provides an implementation of the Anasazi::OperatorTraits class
   * using LOCA operators.  The operator is represented by the
   * LOCA::AnasaziOperator::AbstractStrategy class which serves as a base class
   * for many different types of operators, including shift-invert and
   * Generalized Cayley transformations.
   */
  template <>
  class OperatorTraits<double,
               NOX::Abstract::MultiVector,
               LOCA::AnasaziOperator::AbstractStrategy>
  {
  public:

    //! Synony for double
    typedef double ScalarType;

    //! Synonym for NOX::Abstract::MultiVector
    typedef NOX::Abstract::MultiVector MV;

    //! Synonym for LOCA::AnasaziOperator::AbstractStrategy
    typedef LOCA::AnasaziOperator::AbstractStrategy OP;

    //@{ Matrix/Operator application method.

    //! Application method which performs operation <b>y = Op*x</b>
    /*!
     *\return Status of the operation.
     */
    static void Apply ( const OP& Op, const MV& x, MV& y )
    {
      Op.apply(x,y);
    };

    //@}

  }; // class MultiVecTraits

} // namespace Anasazi

#endif
