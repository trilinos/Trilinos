#ifndef IFPACK_DIAG_PRECONDITIONER_H
#define IFPACK_DIAG_PRECONDITIONER_H

#include "Ifpack_ConfigDefs.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
class Epetra_BlockMap;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Comm;

using namespace std;

//! Ifpack_DiagPreconditioner: a class for diagonal preconditioning.
/*
Ifpack_DiagPreconditioner: a class to wrap a vector as diagonal preconditioner. The preconditioner is simply defined by
\f[
z_i = D_i r_i,
\f]
where \f$r,z\f$ are the vector to be preconditioned and the preconditioned vector, and \f$D_i\f$ is the i-th element of the scaling vector.

\author Marzio Sala, ETHZ/D-INFK

\date Last updated on 17-Apr-06

 */
class Ifpack_DiagPreconditioner : public Epetra_Operator
{
  public:

    //! ctor
    Ifpack_DiagPreconditioner(const Epetra_Map& DomainMap,
                              const Epetra_Map& RangeMap,
                              const Epetra_Vector& diag);

    //! dtor
    ~Ifpack_DiagPreconditioner();

    int SetUseTranspose(bool UseTranspose)
    {
      UseTranspose_ = UseTranspose;
      return(0);
    }

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    double NormInf() const
    {
      return(-1.0);
    }

    const char* Label() const
    {
      return("Ifpack_DiagPreconditioner");
    }

    bool UseTranspose() const
    {
      return(UseTranspose_);
    }

    bool HasNormInf() const
    {
      return(false);
    }

    const Epetra_Comm& Comm() const
    {
      return(diag_.Comm());
    }

    const Epetra_Map& OperatorDomainMap() const
    {
      return(RangeMap_);
    }

    const Epetra_Map& OperatorRangeMap() const
    {
      return(DomainMap_);
    }

    const Epetra_BlockMap& Map() const
    {
      return(diag_.Map());
    }

  private:
    bool UseTranspose_;
    const Epetra_Map& DomainMap_;
    const Epetra_Map& RangeMap_;
    const Epetra_Vector& diag_;
};

#endif
