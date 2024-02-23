// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_EPETRAOPERATOR_HPP
#define MUELU_EPETRAOPERATOR_HPP

//! @file

#include <Epetra_Operator.h>
#include "MueLu_Hierarchy.hpp"
// TODO: Kokkos headers

#if defined(HAVE_MUELU_SERIAL) and defined(HAVE_MUELU_EPETRA)

namespace MueLu {

/*! @class EpetraOperator
    @brief Turns a MueLu::Hierarchy into a Epetra_Operator.
    It allows MueLu to be used as a preconditioner for AztecOO (for instance).
*/
class EpetraOperator : public Epetra_Operator {
  typedef double SC;
  typedef int LO;
  typedef int GO;
  typedef Xpetra::EpetraNode NO;

  typedef Xpetra::Matrix<SC, LO, GO, NO> Matrix;
  typedef MueLu::Hierarchy<SC, LO, GO, NO> Hierarchy;
  typedef MueLu::Utilities<SC, LO, GO, NO> Utils;

 public:
  //! @name Constructor/Destructor
  //@{

  //! Constructor
  EpetraOperator(const RCP<Hierarchy>& H)
    : Hierarchy_(H) {}

  //! Destructor.
  virtual ~EpetraOperator() {}

  //@}

  int SetUseTranspose(bool /* UseTransposeBool */) { return -1; }

  //! @name Mathematical functions
  //@{

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*!
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const { return -1; }

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*!
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must
    support the case where X and Y are the same object.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

     \warning This method must not be called unless HasNormInf() returns true.
  */
  double NormInf() const { return 0; }
  //@}

  //! @name Attribute access functions
  //@{

  //! Returns a character string describing the operator
  const char* Label() const { return "MueLu::Hierarchy"; }

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const { return false; }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const { return 0; }

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm& Comm() const;

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map& OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map& OperatorRangeMap() const;

  //@}

  //! @name MueLu specific
  //@{

  //! Direct access to the underlying MueLu::Hierarchy.
  RCP<Hierarchy> GetHierarchy() const { return Hierarchy_; }

  //@}

 private:
  RCP<Hierarchy> Hierarchy_;
};

}  // namespace MueLu

#endif  // HAVE_MUELU_EPETRA and HAVE_MUELU_SERIAL

#endif  // MUELU_EPETRAOPERATOR_HPP
