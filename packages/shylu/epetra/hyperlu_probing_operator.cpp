
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Amesos_BaseSolver.h>
#include <Epetra_MultiVector.h>
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

    G_->Multiply(false, X, temp2);
    C_->Multiply(false, X, temp);

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

    int lda;
    double *values;
    int err = temp.ExtractView(&values, &lda);
    assert (err == 0);

    // copy to local vector //TODO: OMP parallel
    assert(lda == nrows);
    for (int v = 0; v < nvectors; v++)
    {
       for (int i = 0; i < nrows; i++)
       {
           err = ltemp.ReplaceMyValue(i, v, values[i+v*lda]);
           assert (err == 0);
       }
    }

    LP_->SetRHS(&ltemp);
    LP_->SetLHS(&localX);
    solver_->Solve();
    err = localX.ExtractView(&values, &lda);
    assert (err == 0);

    //Copy back to dist vector //TODO: OMP parallel
    for (int v = 0; v < nvectors; v++)
    {
       for (int i = 0; i < nrows; i++)
       {
           err = temp.ReplaceMyValue(i, v, values[i+v*lda]);
           assert (err == 0);
       }
    }

    R_->Multiply(false, temp, Y);
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
