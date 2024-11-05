#ifndef ML_INVERSEOPERATOR_H
#define ML_INVERSEOPERATOR_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_InverseOperator.h

\brief Base class for smoothers and coarse solvers.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"
#include "ml_MultiLevelPreconditioner.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Operator.h"
#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos {
  class List;
}
class Ifpack_Preconditioner;

namespace MLAPI {

class MultiLevel;

/*!
 * \class InverseOperator
 *
 * \brief InverseOperator: basic class to define smoother and coarse solvers.
 *
 *
 * \author Marzio Sala, D-INFK/ETHZ.
 *
 * \date Last updated on Mar-06.
 */

class InverseOperator : public BaseOperator, public CompObject, public TimeObject {

public:
  // @{ \name Constructors and destructors.

  //! Empty constructor.
  InverseOperator() {}

  //! Constructor for a given Operator and type, and default parameters.
  InverseOperator(const Operator& Op, const std::string Type);

  //! Constructor for a given Operator, type and parameters.
  InverseOperator(const Operator& Op, const std::string Type,
                  Teuchos::ParameterList& List);

  //! Copy constructor.
  InverseOperator(const InverseOperator& RHS);

  //! Destructor.
  ~InverseOperator()
  {}

  // @}
  // @{ \name Overloaded operators

  //! Operator =.
  InverseOperator& operator=(const InverseOperator& RHS);

  // @}
  // @{ Reshaping methods

  //! Resets \c this object.
  void Reshape();

  //! Reshapes the object with default values.
  void Reshape(const Operator& Op, const std::string Type);

  //! Reshapes the object by setting the Operator and the specified type.
  void Reshape(const Operator& Op, const std::string Type,
               Teuchos::ParameterList& List,
               Teuchos::ParameterList* pushlist = NULL);

  //! Reshape with preconstructed smoother as Ifpack_Preconditioner
  void Reshape(Ifpack_Preconditioner* prec, const Operator& Op,
               const bool ownership);

  // @}
  // @{ Get and Set methods.

  //! Returns a reference to the range space of \c this object.
  const Space GetOperatorRangeSpace() const;

  //! Returns a reference to the domain space of \c this object.
  const Space GetOperatorDomainSpace() const;

  //! Returns a reference to the range space of \c this object.
  const Space GetRangeSpace() const;

  //! Returns a reference to the domain space of \c this object.
  const Space GetDomainSpace() const;

  //! Returns pointer of the internally stored ML_Epetra::RowMatrix object.
  const Teuchos::RefCountPtr<Epetra_RowMatrix> RCPRowMatrix() const;

  //! Returns pointer of the internally stored ML_Epetra::RowMatrix object.
  Epetra_RowMatrix* RowMatrix() const;

  //! Returns a reference to the Operator of which \c this object defines the inverse.
  const Operator& GetOperator() const;

  //! Returns a pointer to the internally stored IFPACK preconditioner.
  Teuchos::RefCountPtr<Ifpack_Preconditioner>& GetRCPData();

  //! Returns a pointer to the internally stored IFPACK preconditioner.
  Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner>& GetRCPMLPrec();

  //! Returns a pointer to the internally stored IFPACK preconditioner.
  const Teuchos::RefCountPtr<Ifpack_Preconditioner>& GetRCPData() const;

  //! Returns a pointer to the internally stored ML preconditioner.
  const Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner>& GetRCPMLPrec() const;

  // @}
  // @{ Mathematical methods

  //! Applies \c this object to vector \c lhs, returns values in \c rhs.
  int Apply(const MultiVector& x, MultiVector& y) const;

  //! Applies the operator to LHS, returns the results.
  MultiVector operator()(const MultiVector& LHS);

  //! Applies the operator to LHS using RHS as initial solution, returns the results.
  MultiVector operator()(const MultiVector& LHS,
                         const MultiVector& RHS);

  // @}
  // @{ \name Miscellaneous methods

  //! Prints out basic information about \c this object.
  std::ostream& Print(std::ostream& os, const bool verbose = true) const;

private:

  // @}
  // @{ \name Private data and methods

  void Destroy();

  //! Operator of which \c this object define the inverse.
  Operator Op_;
  //! Wrapper for IFPACK
  Teuchos::RefCountPtr<Epetra_RowMatrix> RCPRowMatrix_;
  //! IFPACK preconditioner.
  Teuchos::RefCountPtr<Ifpack_Preconditioner> RCPData_;
  //! ML preconditioner
  Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner>  RCPMLPrec_;
  // @}

}; // InverseOperator

} // namespace MLAPI

#endif // ML_INVERSEOPERATOR_H
