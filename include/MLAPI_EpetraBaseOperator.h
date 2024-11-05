#ifndef MLAPI_EPETRAPRECONDITIONER_H
#define MLAPI_EPETRAPRECONDITIONER_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_EpetraBaseOperator.h

\brief Basic class to wrap MLAPI::InverseOperator into Epetra_Operator.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"

#include "Epetra_Operator.h"
#include "MLAPI_Error.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Workspace.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "ml_epetra.h"

namespace MLAPI {

/*!
\class EpetraBaseOperator

\brief Basic class to wrap MLAPI::InverseOperator into Epetra_Operator.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.

*/

class EpetraBaseOperator : public Epetra_Operator {

public:

  //! Constructor.
  EpetraBaseOperator(const Epetra_Map& inMap,
                     const BaseOperator& Op) :
    Map_(inMap),
    Op_(Op)
  {}

  //! Destructor.
  virtual ~EpetraBaseOperator() {}

  //! Applies the operator to \c X, returns the results in \c Y.
  /*! \note Apply() and ApplyInverse() are the SAME function!
   */
  int ApplyInverse(const Epetra_MultiVector& X_Epetra,
                   Epetra_MultiVector& Y_Epetra) const
  {
    return(Apply(X_Epetra, Y_Epetra));
  }

  //! Sets the use of tranpose (NOT IMPLEMENTED).
  virtual int SetUseTranspose(bool /* UseTransposeFlag */)
  {
    ML_CHK_ERR(-1);
  }

  //! Applies the operator to \c X, returns the results in \c Y.
  virtual int Apply(const Epetra_MultiVector& X_Epetra,
                    Epetra_MultiVector& Y_Epetra) const
  {
    // NOTE: no checks on maps. These checks can be
    // expensive, and I prefer to skip them.

    if (X_Epetra.NumVectors() != Y_Epetra.NumVectors())
      ML_THROW("X.NumVectors() != Y.NumVectors(), " +
               GetString(X_Epetra.NumVectors()) + " vs. " +
               GetString(Y_Epetra.NumVectors()), -1);

    // FIXME: this is not the most efficient way (though for
    // ML should be the same, as there is not native support
    // for multivectors I am using)

    for (int v = 0 ; v < X_Epetra.NumVectors() ; ++v) {

      MultiVector X_ML(Op_.GetOperatorDomainSpace(),(double**)&(X_Epetra[v]), 1);

      // need additional vector for AztecOO
      MultiVector Y_ML(Op_.GetOperatorRangeSpace(), 1);

      ML_CHK_ERR(Op_.Apply(X_ML,Y_ML));

      int n = Y_Epetra.MyLength();
      int incr = 1;
      DCOPY_F77(&n, Y_ML.GetValues(0), &incr, &(Y_Epetra[v][0]), &incr);
    }

    return(0);
  }

  //! NOT IMPLEMENTED.
  virtual double NormInf() const
  {
    return(-1.0);
  }

  //! Returns the label of \c this object.
  virtual const char* Label() const
  {
    return(Op_.GetLabel().c_str());
  }

  //! Returns \c false.
  virtual bool UseTranspose() const
  {
    return(false);
  }

  //! NOT IMPLEMENTED.
  virtual bool HasNormInf() const
  {
    return(false);
  }

  //! Returns a reference to the communicator object.
  virtual const Epetra_Comm& Comm() const
  {
    return(GetEpetra_Comm());
  }

  //! Returns a reference to the OperatorDomainMap.
  virtual const Epetra_Map& OperatorDomainMap() const
  {
    return(Map_);
  }

  //! Returns a reference to the OperatorRangeMap.
  virtual const Epetra_Map& OperatorRangeMap() const
  {
    return(Map_);
  }

  //! Returns a reference to the Map of \c this object.
  virtual const Epetra_Map& Map() const
  {
    return(Map_);
  }

  const BaseOperator& GetBaseOperator() const
  {
    return(Op_);
  }

private:

  //! Copy constructor (should not be used).
  EpetraBaseOperator(const EpetraBaseOperator& rhs) :
    Map_(rhs.Map()),
    Op_(rhs.GetBaseOperator())
  { }

  //! operator= (should not be used).
  EpetraBaseOperator& operator=(const EpetraBaseOperator& /* rhs */)
  {
    return(*this);
  }

  //! Reference to the Map.
  const Epetra_Map& Map_;
  //! Reference to the MLAPI BaseOperator.
  const BaseOperator& Op_;

}; // Epetra_MultiVector

} // namespace MLAPI

#endif
