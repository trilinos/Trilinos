/** \file Ifpack_Amesos_Schur.cpp

    \brief Use Amesos within IFPACK.

    \author Siva Rajamanickam

*/

#include "Ifpack_Amesos_Schur.h"

Ifpack_Amesos_Schur::Ifpack_Amesos_Schur(Epetra_CrsMatrix* A):
    A_(A),
    IsParallel_(true),
    IsInitialized_(false),
    IsComputed_(false),
    Label_(),
    NumApplyInverse_(0),
    Time_(A_->Comm())
{
}

void Ifpack_Amesos_Schur::Destroy()
{
    if (IsInitialized_)
    {
    }
    if (IsComputed_)
    {
    }
}

int Ifpack_Amesos_Schur::Initialize()
{
    if(Comm().NumProc() != 1) 
        IsParallel_ = true;
    else 
        IsParallel_ = false;

    LP_ = Teuchos::RCP<Epetra_LinearProblem> (new Epetra_LinearProblem());
    Amesos Factory;
    char* SolverType = "Amesos_Klu";
    bool IsAvailable = Factory.Query(SolverType);
    assert(IsAvailable == true);
    Solver_ = Teuchos::RCP<Amesos_BaseSolver> (Factory.Create(SolverType,
                                 *LP_));

    Teuchos::ParameterList aList;
    aList.set("Reindex", true);

    LP_->SetOperator(A_);
    Solver_->SetParameters(aList);

    Solver_->SymbolicFactorization();

    IsInitialized_ = true;
    return 0;
}

int Ifpack_Amesos_Schur::SetParameters(Teuchos::ParameterList& parameterlist)
{
    List_ = parameterlist;
    return 0;
}

int Ifpack_Amesos_Schur::Compute()
{
    IsComputed_ = true;
    Solver_->NumericFactorization();
    return 0;
}

int Ifpack_Amesos_Schur::ApplyInverse(const Epetra_MultiVector& X,
    Epetra_MultiVector& Y) const
{
    NumApplyInverse_++;
    LP_->SetRHS(const_cast<Epetra_MultiVector *>(&X));
    LP_->SetLHS(&Y);
    Solver_->Solve();
    return 0;
}

//! Computes the estimated condition number and returns the value.
double Ifpack_Amesos_Schur::Condest(const Ifpack_CondestType CT,
     const int MaxIters, const double Tol, Epetra_RowMatrix* Matrix_in)
{
    return -1.0;
}

//! Prints on stream basic information about \c this object.
ostream& Ifpack_Amesos_Schur::Print(ostream& os) const
{
    os << " !!!!!!!!! " << endl;
    return os;
}
