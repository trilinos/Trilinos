#ifndef MLAPI_EPETRAPRECONDITIONER_H
#define MLAPI_EPETRAPRECONDITIONER_H

#include "Epetra_Operator.h"
#include "MLAPI_Workspace.h"

class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Comm;

namespace MLAPI {

class Preconditioner;

/*!
\class EpetraPreconditioner

\brief Basic class to wrap MLAPI::Preconditioner into Epetra_Operator.

\author Marzio Sala, SNL 9214.

\date Last updated on 07-Jan-05.

*/

class EpetraPreconditioner : public Epetra_Operator {

public:

  //! Constructor.
  EpetraPreconditioner(const Epetra_Map& Map,
                       const Preconditioner& Prec) :
    Map_(Map),
    Prec_(Prec)
  {}

  //! Destructor.
  virtual ~EpetraPreconditioner() {}

  //! Applies the preconditioner to LHS, returns the result in RHS.
  /*! \note this is AztecOO compliant (that is, LHS and RHS can point
   *        to the same memory location.
   */
  int ApplyInverse(const Epetra_MultiVector& LHS,
                   Epetra_MultiVector& RHS) const
  {

    assert (LHS.NumVectors() == 1);

    DoubleVector LHS2(Prec_.DomainSpace(),&(LHS[0][0]));

    // need additional vector for AztecOO
    DoubleVector RHS2(Prec_.RangeSpace());

    ML_CHK_ERR(Prec_.Solve(LHS2,RHS2));
    
    int n = RHS.MyLength();
    int incr = 1;
    DCOPY_F77(&n, RHS2.Values(), &incr, RHS[0], &incr);

    return(0);
  }
  
  //! Sets the use of tranpose (NOT IMPLEMENTED).
  virtual int SetUseTranspose(bool UseTranspose)
  {
    ML_CHK_ERR(-1);
  }

  //! NOT IMPLEMENTED.
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  {
    ML_CHK_ERR(-1);
  }

  //! NOT IMPLEMENTED.
  virtual double NormInf() const
  {
    return(-1.0);
  }

  //! Returns the label of \c this object.
  virtual const char* Label() const
  {
    return("MLAPI::EpetraPreconditioner");
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
private:
  //! Reference to the Map.
  const Epetra_Map& Map_;
  //! Reference to the MLAPI Preconditioner.
  const Preconditioner& Prec_;

}; // Epetra_MultiVector

} // namespace MLAPI

#endif

