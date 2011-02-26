
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Amesos_BaseSolver.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

class HyperLU_Probing_Operator : public virtual Epetra_Operator
{

    public:
    HyperLU_Probing_Operator(Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *LocalDRowMap);

    int SetUseTranspose(bool useTranspose);

    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    double NormInf() const;

    const char *Label() const;

    bool UseTranspose() const;

    bool HasNormInf() const;

    const Epetra_Comm& Comm() const;

    const Epetra_Map& OperatorDomainMap() const;

    const Epetra_Map& OperatorRangeMap() const;

    Epetra_CrsMatrix *G_;
    Epetra_CrsMatrix *R_;
    Epetra_LinearProblem *LP_;
    Amesos_BaseSolver *solver_;
    Epetra_CrsMatrix *C_;
    Epetra_Map *localDRowMap_;

};
