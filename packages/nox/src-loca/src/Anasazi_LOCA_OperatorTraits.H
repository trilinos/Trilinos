// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
