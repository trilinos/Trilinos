#ifndef ML_MATRIX_FREE_PRECONDITIONER
#define ML_MATRIX_FREE_PRECONDITIONER

#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT)
#include <string>
#include "ml_epetra.h"
#include "Epetra_Operator.h"
#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"

class Epetra_Time;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_CrsGraph;
class Epetra_MultiVector;
class Epetra_CrsMatrix;
class Epetra_FECrsMatrix;

namespace ML_Epetra {

class MultiLevelPreconditioner;

class MatrixFreePreconditioner : public Epetra_Operator 
{
  public:
    //! Constructor
    MatrixFreePreconditioner(const Epetra_Operator& Operator,
                             const Epetra_CrsGraph& Graph,
                             Teuchos::ParameterList& List,
                             Epetra_MultiVector& NullSpace);

    //! destructor
    ~MatrixFreePreconditioner();

    int SetUseTranspose(bool UseTranspose);

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    double NormInf() const
    {
      return(-1.0);
    }

    const char * Label() const
    {
      return(Label_.c_str());
    }

    bool UseTranspose() const
    {
      return(false);
    }

    bool HasNormInf() const
    {
      return(false);
    }

    //! Returns a reference to the communicator object.
    const Epetra_Comm & Comm() const
    {
      return(Comm_);
    }

    //! Returns the domain map of the operator.
    const Epetra_Map & OperatorDomainMap() const
    {
      return(Operator_.OperatorDomainMap());
    }

    //! Returns the range map of the operator.
    const Epetra_Map & OperatorRangeMap() const
    {
      return(Operator_.OperatorRangeMap());
    }

    const Epetra_CrsMatrix& C() const
    {
      return(*C_);
    }

    const Epetra_FECrsMatrix& R() const
    {
      return(*R_);
    }

    int Coarsen(ML_Operator* A, ML_Aggregate** aggr, ML_Operator** P, 
                ML_Operator** R, ML_Operator** C, int = 1, int = 1, double* = NULL);

    ML_Comm* Comm_ML()
    {
      return(Comm_ML_);
    }

    inline int MyPID() const
    {
      return(Comm().MyPID());
    }

    inline int NumProc() const
    {
      return(Comm().NumProc());
    }

    bool IsComputed() const
    {
      return(IsComputed_);
    }

  private:
    Epetra_Time& Time()
    {
      return(*Time_);
    }

    ML_Comm* Comm_ML_;

    //! Computes the preconditioner.
    int Compute(Epetra_MultiVector& NullSpace);

    bool IsComputed_;
    //! Label of this object
    std::string Label_; 
    //! Communicator object
    const Epetra_Comm& Comm_;
    //! Fine-level operator
    const Epetra_Operator& Operator_;
    //! Fine-level graph
    const Epetra_CrsGraph& Graph_;
    //! List containing all the parameters
    Teuchos::ParameterList List_;
    //! Prolongator from coarse to fine
    Epetra_CrsMatrix* P_;
    Epetra_FECrsMatrix* R_;
    Epetra_CrsMatrix* C_;

    Epetra_Time* Time_;

    MultiLevelPreconditioner* MLP_;

}; // class MatrixFreePreconditioner

} // namespace ML_Epetra

#endif // HAVE_ML_EPETRA
#endif // ML_MATRIX_FREE_PRECONDITIONER
