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

#ifndef LOCA_EPETRA_TRANSPOSELINEARSYSTEM_LEFTPRECONDITIONING_H
#define LOCA_EPETRA_TRANSPOSELINEARSYSTEM_LEFTPRECONDITIONING_H

#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H" // base class

// forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  namespace Epetra {
    class LinearSystem;
  }
}
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

  namespace Epetra {

    namespace TransposeLinearSystem {

      /*!
       * \brief Method for solving the transpose of a linear system by
       * transposing the preconditioner and switching to left
       * preconditioning.
       */
      /*
       * Note:  This class does not support scaling.
       */
      class LeftPreconditioning : public virtual AbstractStrategy {

      public:

    //! Constructor
    LeftPreconditioning(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
        const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys);

    //! Destructor
    virtual ~LeftPreconditioning();

    /*!
     * \brief Applies the inverse of the Jacobian matrix transpose to the
     * given input vector and puts the answer in result.
     *
     * Computes
     * \f[ v = J^{-T} u, \f]
     * where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector,
     * and \f$v\f$ is the result vector.
     *
     * The parameter list contains the linear solver options.
     */
    virtual bool
    applyJacobianTransposeInverse(Teuchos::ParameterList &params,
                      const NOX::Epetra::Vector &input,
                      NOX::Epetra::Vector &result);



    //! Evaluates the Jacobian-transpose based on the solution vector x.
    /*!
     * Note:  For flexibility, this method does not compute the original
     * Jacobian matrix.  It uses whatever is currently stored in the
     * linear system.
     */
    virtual bool
    createJacobianTranspose();

    /*!
     * \brief Explicitly constructs a preconditioner based on the solution
     * vector x and the parameter list p.
     */
    /*!
     * Note:  \c x is only needed for user-supplied preconditioners.
     * When using a built-in preconditioner (e.g., Ifpack), \c x will
     * note be used.
     */
    virtual bool
    createTransposePreconditioner(const NOX::Epetra::Vector& x,
                      Teuchos::ParameterList& p);

    //! Get Jacobian-transpose operator
    virtual Teuchos::RCP<Epetra_Operator>
    getJacobianTransposeOperator();

    //! Get transpose-preconditioner
    virtual Teuchos::RCP<Epetra_Operator>
    getTransposePreconditioner();

    //! Set Jacobian-transpose operator
    virtual void setJacobianTransposeOperator(
           const Teuchos::RCP<Epetra_Operator>& new_jac_trans);

    //! Set transpose-preconditioner
    virtual void setTransposePreconditioner(
          const Teuchos::RCP<Epetra_Operator>& new_prec_trans);

      protected:

    //! Global data object
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! Linear system solver
    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;

    //! Jacobian operator
    Teuchos::RCP<Epetra_Operator> jac;

    //! Preconditioner
    Teuchos::RCP<Epetra_Operator> prec;

      }; // class TransposePreconditioner

    } // namespace TransposeLinearSystem

  } // namespace Epetra

} // namespace LOCA


#endif
