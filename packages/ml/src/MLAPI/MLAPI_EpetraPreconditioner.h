#ifndef MLAPI_EPETRAPRECONDITIONER_H
#define MLAPI_EPETRAPRECONDITIONER_H

#include "Epetra_Operator.h"

class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Comm;

namespace MLAPI {

class Preconditioner;

class EpetraPreconditioner : public Epetra_Operator {

public:

  EpetraPreconditioner(const Epetra_Comm& Comm,
                       const Epetra_Map& Map,
                       const Preconditioner& Prec) :
    Comm_(Comm),
    Map_(Map),
    Prec_(Prec)
  {}

  virtual ~EpetraPreconditioner() {}

  // this is AztecOO compliant
  int ApplyInverse(const Epetra_MultiVector& LHS,
                   Epetra_MultiVector& RHS) const
  {

    assert (LHS.NumVectors() == 1);

    DoubleVector LHS2(Prec_.DomainSpace(),&(LHS[0][0]));
    // need additional vector for AztecOO
    DoubleVector RHS2(Prec_.RangeSpace());

    ML_CHK_ERR(Prec_.Solve(LHS2,RHS2));
    
    // copy back the result
    for (int i = 0 ; i < Map_.NumMyElements() ; ++i) {
      RHS[0][i] = RHS2(i);
    }

    return(0);
  }
  
  virtual int SetUseTranspose(bool UseTranspose)
  {
    ML_CHK_ERR(-1);
  }

  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  {
    ML_CHK_ERR(-1);
  }

  virtual double NormInf() const
  {
    return(-1.0);
  }

  virtual const char* Label() const
  {
    return("MLAPI::EpetraPreconditioner");
  }

  virtual bool UseTranspose() const
  {
    return(false);
  }

  virtual bool HasNormInf() const 
  {
    return(false);
  }

  virtual const Epetra_Comm & Comm() const
  {
    return(Comm_);
  }

  virtual const Epetra_Map & OperatorDomainMap() const
  {
    return(Map_);
  }

  virtual const Epetra_Map & OperatorRangeMap() const
  {
    return(Map_);
  }

  virtual const Epetra_Map& Map() const
  {
    return(Map_);
  }
private:
  const Epetra_Comm& Comm_;
  const Epetra_Map& Map_;
  const Preconditioner& Prec_;

}; // Epetra_MultiVector

} // namespace MLAPI

#endif

