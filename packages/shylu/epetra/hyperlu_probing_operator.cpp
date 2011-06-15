
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Amesos_BaseSolver.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_Time.hpp>
#include "hyperlu_probing_operator.h"

// TODO: 1. ltemp is not needed in the all local case.

HyperLU_Probing_Operator::HyperLU_Probing_Operator(Epetra_CrsMatrix *G, 
    Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap, int nvectors)
{
    G_ = G;
    R_ = R;
    LP_ = LP;
    solver_ = solver;
    C_ = C;
    localDRowMap_ = localDRowMap;
    ResetTempVectors(nvectors);
    nvectors_ = nvectors;

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

int HyperLU_Probing_Operator::SetUseTranspose(bool useTranspose)
{
    return 0;
}

int HyperLU_Probing_Operator::Apply(const Epetra_MultiVector &X,
            Epetra_MultiVector &Y) const
{
#ifdef TIMING_OUTPUT
    apply_time_->start();
#endif

    int nvectors = X.NumVectors();
    bool local = (G_->Comm().NumProc() == 1);
    //cout << "No of colors after probing" << nvectors << endl;

#ifdef TIMING_OUTPUT
    matvec_time_->start();
#endif

    G_->Multiply(false, X, *temp2);
    C_->Multiply(false, X, *temp);

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

    int err;
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
        LP_->SetRHS(temp.getRawPtr());
    }
    LP_->SetLHS(localX.getRawPtr());
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
        R_->Multiply(false, *localX, Y);
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
    return 0;
}

void HyperLU_Probing_Operator::PrintTimingInfo()
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

void HyperLU_Probing_Operator::ResetTempVectors(int nvectors)
{
    using Teuchos::RCP;
    nvectors_ = nvectors;
    // If vectors were created already, they will be freed.
    temp = RCP<Epetra_MultiVector>(new Epetra_MultiVector(C_->RowMap(),
                                     nvectors));
    temp2 = RCP<Epetra_MultiVector>(new Epetra_MultiVector(G_->RowMap(),
                                     nvectors));
    ltemp = RCP<Epetra_MultiVector>(new Epetra_MultiVector(*localDRowMap_,
                                     nvectors));
    localX = RCP<Epetra_MultiVector>(new Epetra_MultiVector(*localDRowMap_,
                                     nvectors));
}

int HyperLU_Probing_Operator::ApplyInverse(const Epetra_MultiVector &X,
             Epetra_MultiVector &Y) const
{
    return 0;
}

double HyperLU_Probing_Operator::NormInf() const
{
   return -1;
}

const char* HyperLU_Probing_Operator::Label() const
{
    return "Hyperlu probing";
}

bool HyperLU_Probing_Operator::UseTranspose() const
{
    return false;
}

bool HyperLU_Probing_Operator::HasNormInf() const
{
    return false;
}

const Epetra_Comm& HyperLU_Probing_Operator::Comm() const
{
    return C_->Comm();
}

const Epetra_Map& HyperLU_Probing_Operator::OperatorDomainMap() const
{
    //return C_->ColMap();
    return C_->DomainMap();
}

const Epetra_Map& HyperLU_Probing_Operator::OperatorRangeMap() const
{
    //return R_->RowMap();
    return R_->RangeMap();
}
