// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_ADAPTERS_AZTECOO_MUELU_AZTECEPETRAOPERATOR_HPP_
#define PACKAGES_MUELU_ADAPTERS_AZTECOO_MUELU_AZTECEPETRAOPERATOR_HPP_

#include <Epetra_Operator.h>

#include "Xpetra_Operator.hpp"

#if defined(HAVE_MUELU_SERIAL) and defined(HAVE_MUELU_EPETRA)

namespace MueLu {

/*! @class AztecEpetraOperator
    @brief Turns a Xpetra::Operator into a Epetra_Operator.
    It allows an Xpetra::Operator to be used as a preconditioner for AztecOO (for instance).

    Currently only used for RefMaxwell.
*/
class AztecEpetraOperator : public Epetra_Operator {
  typedef double SC;
  typedef int LO;
  typedef int GO;
  typedef Xpetra::EpetraNode NO;

  typedef Xpetra::Map<LO, GO, NO> Map;
  typedef Xpetra::EpetraMapT<GO, NO> EpetraMap;
  typedef Xpetra::Operator<SC, LO, GO, NO> Operator;

 public:
  //! @name Constructor/Destructor
  //@{

  //! Constructor
  AztecEpetraOperator(const Teuchos::RCP<Operator>& Op)
    : xOp_(Op) {}

  //! Destructor.
  virtual ~AztecEpetraOperator() {}

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
  const char* Label() const { return "MueLu::AztecEpetraOperator"; }

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

  //! Direct access to the underlying Xpetra::Operator.
  Teuchos::RCP<Operator> GetOperator() const { return xOp_; }

  //@}

 private:
  Teuchos::RCP<Operator> xOp_;
};

}  // namespace MueLu

#endif

#endif /* PACKAGES_MUELU_ADAPTERS_AZTECOO_MUELU_AZTECEPETRAOPERATOR_HPP_ */
