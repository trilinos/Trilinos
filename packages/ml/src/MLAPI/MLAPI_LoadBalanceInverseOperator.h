#ifndef ML_LOADBALANCEINVERSEOPERATOR_H
#define ML_LOADBALANCEINVERSEOPERATOR_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_LoadBalanceInverseOperator.h

\brief wraps an MLAPI inverseoperator with zero rows on some processors.

\author Michael Gee, TU Munich.

\date Last updated on July-10.
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
#include "MLAPI_InverseOperator.h"
#include "MLAPI_LoadBalanceOperator.h"
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

class LoadBalanceInverseOperator : public InverseOperator {

public:
  // @{ \name Constructors and destructors.

  //! Empty constructor.
  LoadBalanceInverseOperator() : InverseOperator() {}

  //! Copy constructor.
  LoadBalanceInverseOperator(const LoadBalanceInverseOperator& RHS);

  //! Destructor.
  virtual ~LoadBalanceInverseOperator()
  {Destroy();}

  // @}
  // @{ \name Overloaded operators

  //! Operator =.
  LoadBalanceInverseOperator& operator=(const LoadBalanceInverseOperator& RHS);

  // @}
  // @{ Reshaping methods

  //! Resets \c this object.
  void Reshape();

  //! Reshape with preconstructed smoother as Ifpack_Preconditioner
  void Reshape(Ifpack_Preconditioner* prec, const LoadBalanceOperator& Op,
               const bool ownership);

  // @}
  // @{ Get and Set methods.

  //! Returns a bool indicating whether this proc participates in the operator application
  virtual inline bool GetParticipation() const {
    return(GetOperator().GetParticipation());
  }

  //! Returns a reference to the range space of \c this object.
  const Space GetOperatorRangeSpace() const;

  //! Returns a reference to the domain space of \c this object.
  const Space GetOperatorDomainSpace() const;

  //! Returns a reference to the range space of \c this object.
  const Space GetRangeSpace() const;

  //! Returns a reference to the domain space of \c this object.
  const Space GetDomainSpace() const;

  //! Returns pointer of the internally stored ML_Epetra::RowMatrix object.
  const Teuchos::RCP<Epetra_RowMatrix> RCPRowMatrix() const;

  //! Returns pointer of the internally stored ML_Epetra::RowMatrix object.
  Epetra_RowMatrix* RowMatrix() const;

  //! Returns a reference to the Operator of which \c this object defines the inverse.
  const LoadBalanceOperator& GetOperator() const;

  //! Returns a pointer to the internally stored IFPACK preconditioner.
  Teuchos::RCP<Ifpack_Preconditioner>& GetRCPData();

  //! Returns a pointer to the internally stored IFPACK preconditioner.
  const Teuchos::RCP<Ifpack_Preconditioner>& GetRCPData() const;

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
  LoadBalanceOperator Op_;
  //! Wrapper for IFPACK
  Teuchos::RCP<Epetra_RowMatrix> RCPRowMatrix_;
  //! IFPACK preconditioner.
  Teuchos::RCP<Ifpack_Preconditioner> RCPData_;
  // @}

}; // InverseOperator

} // namespace MLAPI

#endif // ML_INVERSEOPERATOR_H
