
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


#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Amesos_BaseSolver.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_Time.hpp>
#include "shylu_probing_operator.h"

// TODO: 1. ltemp is not needed in the all local case.

ShyLU_Probing_Operator::ShyLU_Probing_Operator(
    shylu_symbolic *ssym,   // symbolic structure
    Epetra_CrsMatrix *G,
    Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap, int nvectors)
{
    ssym_ = ssym;
    G_ = G;
    R_ = R;
    LP_ = LP;
    solver_ = solver;
    C_ = C;
    localDRowMap_ = localDRowMap;
    ResetTempVectors(nvectors);
    nvectors_ = nvectors;

    cntApply = 0;

#ifdef TIMING_OUTPUT
    using Teuchos::RCP;
    using Teuchos::Time;
    matvec_time_ = RCP<Time>(new Time("First 2 Matvec time in prober: "));
    localize_time_ = RCP<Time>(new Time("Localize time in prober: "));
    trisolve_time_ = RCP<Time>(new Time("Triangular solve time in prober: "));
    dist_time_ = RCP<Time>(new Time("Distribute time in prober: "));
    matvec2_time_ = RCP<Time>(new Time("Second 1 matvec time in prober: "));
    apply_time_ = RCP<Time>(new Time("Apply time in prober: "));
    update_time_ = RCP<Time>(new Time("Update time in prober: "));
#endif
}

int ShyLU_Probing_Operator::SetUseTranspose(bool useTranspose)
{
    return 0;
}

int ShyLU_Probing_Operator::Apply(const Epetra_MultiVector &X,
            Epetra_MultiVector &Y) const
{
#ifdef TIMING_OUTPUT
    apply_time_->start();
#endif

    int nvectors = X.NumVectors();
    bool local = (C_->Comm().NumProc() == 1);
    int err;
    //cout << "No of colors after probing" << nvectors << endl;

#ifdef TIMING_OUTPUT
    matvec_time_->start();
#endif

    err = G_->Multiply(false, X, *temp2);
    assert(err == 0);
    if (!local)
        err = C_->Multiply(false, X, *temp);
    else
    {
        // localize X
        double *values;
        int mylda;
        X.ExtractView(&values, &mylda);

       Epetra_SerialComm LComm;        // Use Serial Comm for the local blocks.
       Epetra_Map SerialMap(X.Map().NumMyElements(), X.Map().NumMyElements(),
                   X.Map().MyGlobalElements(), 0, LComm);
       Epetra_MultiVector Xl(View, SerialMap, values, mylda, X.NumVectors());
       err = C_->Multiply(false, Xl, *temp);
    }
    assert(err == 0);

#ifdef TIMING_OUTPUT
    matvec_time_->stop();
#endif

    int nrows = C_->RowMap().NumMyElements();

#ifdef DEBUG
    cout << "DEBUG MODE" << endl;
    assert(nrows == localDRowMap_->NumGlobalElements());

    int gids[nrows], gids1[nrows];
    C_->RowMap().MyGlobalElements(gids);
    localDRowMap_->MyGlobalElements(gids1);

    for (int i = 0; i < nrows; i++)
    {
       assert(gids[i] == gids1[i]);
    }
#endif

#ifdef TIMING_OUTPUT
    localize_time_->start();
#endif

    //int err;
    int lda;
    double *values;
    if (!local)
    {
        err = temp->ExtractView(&values, &lda);
        assert (err == 0);

        // copy to local vector //TODO: OMP parallel
        assert(lda == nrows);

    //#pragma omp parallel for shared(nvectors, nrows, values)
        for (int v = 0; v < nvectors; v++)
        {
           for (int i = 0; i < nrows; i++)
           {
               err = ltemp->ReplaceMyValue(i, v, values[i+v*lda]);
               assert (err == 0);
           }
        }
    }

#ifdef TIMING_OUTPUT
    localize_time_->stop();
    trisolve_time_->start();
#endif

    if (!local)
    {
        LP_->SetRHS(ltemp.getRawPtr());
    }
    else
    {
        //LP_->SetRHS(temp.getRawPtr());
    }
    //LP_->SetLHS(localX.getRawPtr());

    //TODO: Why not just in Reset(). Check the distr path.
    ssym_->OrigLP->SetLHS(localX.getRawPtr());
    ssym_->OrigLP->SetRHS(temp.getRawPtr());
    ssym_->ReIdx_LP->fwd();
    solver_->Solve();

#ifdef TIMING_OUTPUT
    trisolve_time_->stop();
    dist_time_->start();
#endif

    if (!local)
    {
        err = localX->ExtractView(&values, &lda);
        assert (err == 0);

        //Copy back to dist vector //TODO: OMP parallel
    //#pragma omp parallel for
        for (int v = 0; v < nvectors; v++)
        {
           for (int i = 0; i < nrows; i++)
           {
               err = temp->ReplaceMyValue(i, v, values[i+v*lda]);
               assert (err == 0);
           }
        }
    }

#ifdef TIMING_OUTPUT
    dist_time_->stop();
    matvec2_time_->start();
#endif

    if (!local)
    {
        R_->Multiply(false, *temp, Y);
    }
    else
    {
        // Should Y be localY in Multiply and then exported to Y ?? TODO:
        // Use view mode ?
        double *values;
        int mylda;
        Y.ExtractView(&values, &mylda);

       Epetra_SerialComm LComm;        // Use Serial Comm for the local blocks.
       Epetra_Map SerialMap(Y.Map().NumMyElements(), Y.Map().NumMyElements(),
                   Y.Map().MyGlobalElements(), 0, LComm);
       Epetra_MultiVector Yl(View, SerialMap, values, mylda, Y.NumVectors());
        R_->Multiply(false, *localX, Yl);
    }

#ifdef TIMING_OUTPUT
    matvec2_time_->stop();
    update_time_->start();
#endif

    err = Y.Update(1.0, *temp2, -1.0);
    //cout << Y.MyLength() << " " << temp2.MyLength() << endl;
    assert(err == 0);

#ifdef TIMING_OUTPUT
    update_time_->stop();
    apply_time_->stop();
#endif
    cntApply++;
    return 0;
}

void ShyLU_Probing_Operator::PrintTimingInfo()
{
#ifdef TIMING_OUTPUT
    cout << matvec_time_->name()<< matvec_time_->totalElapsedTime() << endl;
    cout << localize_time_->name()<< localize_time_->totalElapsedTime() << endl;
    cout << trisolve_time_->name()<< trisolve_time_->totalElapsedTime() << endl;
    cout << dist_time_->name() <<dist_time_->totalElapsedTime() << endl;
    cout << matvec2_time_->name() <<matvec2_time_->totalElapsedTime() << endl;
    cout << update_time_->name() <<update_time_->totalElapsedTime() << endl;
    cout << apply_time_->name() <<apply_time_->totalElapsedTime() << endl;
#endif
}

void ShyLU_Probing_Operator::ResetTempVectors(int nvectors)
{
    using Teuchos::RCP;
    nvectors_ = nvectors;
    if (nvectors <= ssym_->Drhs->NumVectors())
    {
        // If vectors were created already, they will be freed.
        ltemp = Teuchos::RCP<Epetra_MultiVector>
            (new Epetra_MultiVector(View, *(ssym_->Drhs), 0,  nvectors));
        localX = Teuchos::RCP<Epetra_MultiVector>
            (new Epetra_MultiVector(View, *(ssym_->Dlhs), 0,  nvectors));
        temp2 = Teuchos::RCP<Epetra_MultiVector>
            (new Epetra_MultiVector(View, *(ssym_->Gvec), 0,  nvectors));
    }
    else
    {
        ltemp = Teuchos::RCP<Epetra_MultiVector>
            (new Epetra_MultiVector((*(ssym_->Drhs)).Map(), nvectors));
        localX = Teuchos::RCP<Epetra_MultiVector>
            (new Epetra_MultiVector((*(ssym_->Dlhs)).Map(), nvectors));
        temp2 = Teuchos::RCP<Epetra_MultiVector>
            (new Epetra_MultiVector((*(ssym_->Gvec)).Map(), nvectors));
    }
    temp = RCP<Epetra_MultiVector>(new Epetra_MultiVector(C_->RowMap(),
                                     nvectors));
    ssym_->OrigLP->SetLHS(localX.getRawPtr());
    ssym_->OrigLP->SetRHS(temp.getRawPtr());
    ssym_->ReIdx_LP->fwd();
}

int ShyLU_Probing_Operator::ApplyInverse(const Epetra_MultiVector &X,
             Epetra_MultiVector &Y) const
{
    return 0;
}

double ShyLU_Probing_Operator::NormInf() const
{
   return -1;
}

const char* ShyLU_Probing_Operator::Label() const
{
    return "Shylu probing";
}

bool ShyLU_Probing_Operator::UseTranspose() const
{
    return false;
}

bool ShyLU_Probing_Operator::HasNormInf() const
{
    return false;
}

const Epetra_Comm& ShyLU_Probing_Operator::Comm() const
{
    return G_->Comm();
}

const Epetra_Map& ShyLU_Probing_Operator::OperatorDomainMap() const
{
    //return C_->ColMap();
    return G_->DomainMap();
}

const Epetra_Map& ShyLU_Probing_Operator::OperatorRangeMap() const
{
    //return R_->RowMap();
    return G_->RangeMap();
}
