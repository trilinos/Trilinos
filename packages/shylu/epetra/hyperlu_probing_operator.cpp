
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Amesos_BaseSolver.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_Time.hpp>
#include "hyperlu_probing_operator.h"

HyperLU_Probing_Operator::HyperLU_Probing_Operator(Epetra_CrsMatrix *G, 
    Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap)
{
    G_ = G;
    R_ = R;
    LP_ = LP;
    solver_ = solver;
    C_ = C;
    localDRowMap_ = localDRowMap;
}

int HyperLU_Probing_Operator::SetUseTranspose(bool useTranspose)
{
    return 0;
}

int HyperLU_Probing_Operator::Apply(const Epetra_MultiVector &X,
            Epetra_MultiVector &Y) const
{
    int nvectors = X.NumVectors();
    //cout << "No of colors after probing" << nvectors << endl;
    Epetra_MultiVector temp(C_->RowMap(), nvectors);
    Epetra_MultiVector temp2(G_->RowMap(), nvectors);

#ifdef TIMING_OUTPUT_2
    Teuchos::Time ftime("setup time");
    ftime.start();
#endif
    G_->Multiply(false, X, temp2);
    C_->Multiply(false, X, temp);
#ifdef TIMING_OUTPUT_2
    ftime.stop();
    cout << "Time to Compute 2 matvecs" << ftime.totalElapsedTime() << endl;
    ftime.reset();
#endif

    Epetra_MultiVector ltemp(*localDRowMap_, nvectors);
    Epetra_MultiVector localX(*localDRowMap_, nvectors);

    int nrows = C_->RowMap().NumMyElements();
#ifdef DEBUG
    assert(nrows == localDRowMap_->NumGlobalElements());

    int gids[nrows], gids1[nrows];
    C_->RowMap().MyGlobalElements(gids);
    localDRowMap_->MyGlobalElements(gids1);

    for (int i = 0; i < nrows; i++)
    {
       assert(gids[i] == gids1[i]);
    }
#endif
    //cout << "Map check done" << endl;
    // ]

#ifdef TIMING_OUTPUT_2
    ftime.start();
#endif
    int lda;
    double *values;
    int err = temp.ExtractView(&values, &lda);
    assert (err == 0);

    // copy to local vector //TODO: OMP parallel
    assert(lda == nrows);

//#pragma omp parallel for shared(nvectors, nrows, values)
    for (int v = 0; v < nvectors; v++)
    {
       for (int i = 0; i < nrows; i++)
       {
           err = ltemp.ReplaceMyValue(i, v, values[i+v*lda]);
           assert (err == 0);
       }
    }

#ifdef TIMING_OUTPUT_2
    ftime.stop();
    cout << "Time to localize vector" << ftime.totalElapsedTime() << endl;
    ftime.reset();
#endif
    LP_->SetRHS(&ltemp);
    LP_->SetLHS(&localX);
#ifdef TIMING_OUTPUT_2
    ftime.start();
#endif
    solver_->Solve();
#ifdef TIMING_OUTPUT_2
    ftime.stop();
    cout << "Time to do triangular solve" << ftime.totalElapsedTime() << endl;
    ftime.reset();
    ftime.start();
#endif
    err = localX.ExtractView(&values, &lda);
    assert (err == 0);

    //Copy back to dist vector //TODO: OMP parallel
//#pragma omp parallel for
    for (int v = 0; v < nvectors; v++)
    {
       for (int i = 0; i < nrows; i++)
       {
           err = temp.ReplaceMyValue(i, v, values[i+v*lda]);
           assert (err == 0);
       }
    }
#ifdef TIMING_OUTPUT_2
    ftime.stop();
    cout << "Time to distribute vector" << ftime.totalElapsedTime() << endl;
    ftime.reset();
    ftime.start();
#endif

    R_->Multiply(false, temp, Y);
#ifdef TIMING_OUTPUT_2
    ftime.stop();
    cout << "Time to do 1 matvec" << ftime.totalElapsedTime() << endl;
    ftime.reset();
#endif
    err = Y.Update(1.0, temp2, -1.0);
    //cout << Y.MyLength() << " " << temp2.MyLength() << endl;
    assert(err == 0);
    return 0;
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
    return C_->ColMap();
}

const Epetra_Map& HyperLU_Probing_Operator::OperatorRangeMap() const
{
    return R_->RowMap();
}
