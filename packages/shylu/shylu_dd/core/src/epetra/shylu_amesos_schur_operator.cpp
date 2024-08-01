// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file AmesosSchurOperator.cpp

    \brief Use Amesos within IFPACK.

    \author Siva Rajamanickam

*/

#include "shylu_amesos_schur_operator.h"

AmesosSchurOperator::AmesosSchurOperator(Epetra_CrsMatrix* A):
    A_(A),
    IsParallel_(true),
    IsInitialized_(false),
    IsComputed_(false),
    Label_(),
    NumApplyInverse_(0),
    Time_(A_->Comm())
{
}

void AmesosSchurOperator::Destroy()
{
    if (IsInitialized_)
    {
    }
    if (IsComputed_)
    {
    }
}

int AmesosSchurOperator::Initialize()
{
    if(Comm().NumProc() != 1)
        IsParallel_ = true;
    else
        IsParallel_ = false;

    LP_ = Teuchos::RCP<Epetra_LinearProblem> (new Epetra_LinearProblem());
    Amesos Factory;
    const char* SolverType = "Amesos_Klu";
    // mfh 25 May 2015: Remember that in a release build (NDEBUG not
    // defined), assert() gets defined to nothing.  This results in an
    // unused variable warning for IsAvailable.  I've rewritten the
    // code so that in a release build, the query result is ignored.
    // This is still a bad idea -- inexpensive error checks in
    // non-performance-critical code should stay in a release build --
    // but it's not my place to rewrite this code that I didn't write.
#ifdef NDEBUG
    (void) Factory.Query(SolverType);
#else
    bool IsAvailable = Factory.Query(SolverType);
    assert(IsAvailable == true);
#endif // NDEBUG
    Solver_ = Teuchos::RCP<Amesos_BaseSolver> (Factory.Create(SolverType,
                                 *LP_));

    Teuchos::ParameterList aList;
    aList.set("Reindex", true);

    LP_->SetOperator(A_);
    Solver_->SetParameters(aList);
    LP_->SetLHS(0); LP_->SetRHS(0);

    Solver_->SymbolicFactorization();

    IsInitialized_ = true;
    return 0;
}

int AmesosSchurOperator::SetParameters(Teuchos::ParameterList& parameterlist)
{
    List_ = parameterlist;
    return 0;
}

int AmesosSchurOperator::Compute()
{
    IsComputed_ = true;
    Solver_->NumericFactorization();
    return 0;
}

int AmesosSchurOperator::ApplyInverse(const Epetra_MultiVector& X,
    Epetra_MultiVector& Y) const
{
    NumApplyInverse_++;
    LP_->SetRHS(const_cast<Epetra_MultiVector *>(&X));
    LP_->SetLHS(&Y);
    Solver_->Solve();
    return 0;
}

//! Prints on stream basic information about \c this object.
std::ostream& AmesosSchurOperator::Print(std::ostream& os) const
{
    os << " !!!!!!!!! " << std::endl;
    return os;
}
