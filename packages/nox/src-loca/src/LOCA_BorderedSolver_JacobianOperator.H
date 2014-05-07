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

#ifndef LOCA_BORDEREDSOLVER_JACOBIANOPERATOR_H
#define LOCA_BORDEREDSOLVER_JACOBIANOPERATOR_H

#include "Teuchos_RCP.hpp"

#include "LOCA_BorderedSolver_AbstractOperator.H" // base class

namespace LOCA {

  namespace BorderedSolver {

    /*!
     * \brief Bordered solver operator representing the Jacobian as
     * implemented in the NOX::Abstract::Group.
     */
    class JacobianOperator : public LOCA::BorderedSolver::AbstractOperator {

    public:

      //! Constructor
      JacobianOperator(
         const Teuchos::RCP<const NOX::Abstract::Group>& grp);

      //! Destructor
      virtual ~JacobianOperator();

      //! Get group pointer
      virtual Teuchos::RCP<const NOX::Abstract::Group>
      getGroup() const;

      //! Apply the operator
      virtual NOX::Abstract::Group::ReturnType
      apply(const NOX::Abstract::MultiVector& X,
        NOX::Abstract::MultiVector& Y) const;

      //! Apply transpose of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyTranspose(const NOX::Abstract::MultiVector& X,
             NOX::Abstract::MultiVector& Y) const;

      //! Apply inverse of the operator
      virtual NOX::Abstract::Group::ReturnType
      applyInverse(Teuchos::ParameterList& params,
           const NOX::Abstract::MultiVector& B,
           NOX::Abstract::MultiVector& X) const;

      //! Apply inverse transpose of the operator
      /*!
       * Group must be of type LOCA::Abstract::TransposeSolveGroup
       * for this method to be defined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyInverseTranspose(Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& B,
                NOX::Abstract::MultiVector& X) const;

    private:

      //! Private to prohibit copying
      JacobianOperator(const JacobianOperator&);

      //! Private to prohibit copying
      JacobianOperator& operator = (const JacobianOperator&);

    protected:

      //! Pointer to group
      Teuchos::RCP<const NOX::Abstract::Group> grpPtr;

    }; // Class JacobianOperator

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_JACOBIANOPERATOR_H
