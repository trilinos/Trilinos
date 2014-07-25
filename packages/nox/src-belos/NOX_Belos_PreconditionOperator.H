// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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

#ifndef NOX_BELOS_PRECONDITIONOPERATOR_H
#define NOX_BELOS_PRECONDITIONOPERATOR_H

#include "BelosOperator.hpp"            // base class
#include "NOX_Abstract_Group.H"         // data member

namespace NOX {

  namespace Belos {

    /*!
     * \brief Implementation of Belos::Operator for a preconditioner
     * defined by a NOX::Abstract::Group
     */
    /*!
     * This class implements the Belos::Operator interface representing a
     * preconditioner operation using the NOX::Abstract::Group.  It
     * implements the operator Apply() method via
     * applyRightPreconditioningMultiVector.
     */
    class PreconditionOperator : public ::Belos::Operator<double> {

      public:

      //! Constructor.
      PreconditionOperator(NOX::Abstract::Group& g,
               NOX::Parameter::List& preconditionerParameters);

      //! Destructor.
      virtual ~PreconditionOperator();

      //! Apply operator
      /*!
       * Calls underlying group's applyRightPreconditioningMultiVector
       */
      virtual ::Belos::ReturnType
      Apply(const ::Belos::MultiVec<double>& x,
        ::Belos::MultiVec<double>& y,
        ::Belos::ETrans trans=::Belos::NOTRANS) const;

      //! ApplyInverse operator
      /*!
       * Returns ::Belos::Undefined
       */
      virtual ::Belos::ReturnType
      ApplyInverse(const ::Belos::MultiVec<double>& x,
           ::Belos::MultiVec<double>& y,
           ::Belos::ETrans trans=::Belos::NOTRANS) const;

    protected:

      /*!
       * \brief Converts NOX::Abstract::Group::ReturnType to the corresponding
       * Belos type.
       */
      ::Belos::ReturnType
      noxReturnTypeToBelos(NOX::Abstract::Group::ReturnType noxStatus) const;

    protected:

      //! Underlying group
      NOX::Abstract::Group& grp;

      //! Parameter list for applyRightPreconditioning
      NOX::Parameter::List& precondParams;

    };
  } // namespace Belos
} // namespace LOCA

#endif
