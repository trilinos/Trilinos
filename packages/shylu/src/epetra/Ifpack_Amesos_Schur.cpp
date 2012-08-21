
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

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
    LP_->SetLHS(0); LP_->SetRHS(0);

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
