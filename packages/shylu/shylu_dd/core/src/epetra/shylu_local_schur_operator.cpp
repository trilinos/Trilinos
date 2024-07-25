// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Amesos_BaseSolver.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_Time.hpp>
#include "shylu_local_schur_operator.h"
#include "shylu_util.h"


ShyLU_Local_Schur_Operator::ShyLU_Local_Schur_Operator(
    shylu_config *config,
    shylu_symbolic *ssym,   // symbolic structure
    Epetra_CrsMatrix *G,
    Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver,
    Ifpack_Preconditioner *ifSolver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap, int nvectors)
{
    ssym_ = ssym;
    config_ = config;
    G_ = G;
    R_ = R;
    LP_ = LP;
    solver_ = solver;
    ifSolver_ = ifSolver;
    C_ = C;
    localDRowMap_ = localDRowMap;
    orig_lhs_ = ssym_->OrigLP->GetLHS();
    orig_rhs_ = ssym_->OrigLP->GetRHS();
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

int ShyLU_Local_Schur_Operator::SetUseTranspose(bool useTranspose)
{
    return 0;
}

int ShyLU_Local_Schur_Operator::Apply(const Epetra_MultiVector &X,
            Epetra_MultiVector &Y) const
{
#ifdef TIMING_OUTPUT
    apply_time_->start();
#endif

    //int nvectors = X.NumVectors();
    // TODO: For local_scur_operation this will always be local !!!
    // Remove these cases or merge with schur operator.
    //bool local = (G_->Comm().NumProc() == 1);
    int err;

#ifdef TIMING_OUTPUT
    matvec_time_->start();
#endif

    err = G_->Multiply(false, X, *temp2);
    assert((err == 0));
    err = C_->Multiply(false, X, *temp);
    assert((err == 0));

#ifdef TIMING_OUTPUT
    matvec_time_->stop();
#endif


#ifdef DEBUG
    int nrows = C_->RowMap().NumMyElements();
    cout << "DEBUG MODE" << endl;
    SHYLU_CORE_ASSERT((nrows == localDRowMap_->NumGlobalElements()));

    int gids[nrows], gids1[nrows];
    C_->RowMap().MyGlobalElements(gids);
    localDRowMap_->MyGlobalElements(gids1);

    for (int i = 0; i < nrows; i++)
    {
       SHYLU_CORE_ASSERT((gids[i] == gids1[i]));
    }
#endif

#ifdef TIMING_OUTPUT
    localize_time_->start();
#endif

#ifdef TIMING_OUTPUT
    localize_time_->stop();
    trisolve_time_->start();
#endif

    // The lhs and rhs is set in ResetTempVectors()
    if (config_->amesosForDiagonal)
        solver_->Solve();
    else
        ifSolver_->ApplyInverse(*temp, *localX);

#ifdef TIMING_OUTPUT
    trisolve_time_->stop();
    dist_time_->start();
#endif

#ifdef TIMING_OUTPUT
    dist_time_->stop();
    matvec2_time_->start();
#endif

    R_->Multiply(false, *localX, Y);

#ifdef TIMING_OUTPUT
    matvec2_time_->stop();
    update_time_->start();
#endif

    err = Y.Update(1.0, *temp2, -1.0);
    //cout << Y.MyLength() << " " << temp2.MyLength() << endl;
    assert((err == 0));

#ifdef TIMING_OUTPUT
    update_time_->stop();
    apply_time_->stop();
#endif
    //cout << "Out of local schur's Apply" << endl;
    cntApply++;
    return err;
}

void ShyLU_Local_Schur_Operator::PrintTimingInfo()
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

void ShyLU_Local_Schur_Operator::ResetTempVectors(int nvectors)
{
    using Teuchos::RCP;
    nvectors_ = nvectors;
    // If vectors were created already, they will be freed.
    temp = Teuchos::RCP<Epetra_MultiVector>
            (new Epetra_MultiVector(View, *(ssym_->Drhs), 0,  nvectors));
    localX = Teuchos::RCP<Epetra_MultiVector>
            (new Epetra_MultiVector(View, *(ssym_->Dlhs), 0,  nvectors));
    temp2 = Teuchos::RCP<Epetra_MultiVector>
            (new Epetra_MultiVector(View, *(ssym_->Gvec), 0,  nvectors));
    ssym_->OrigLP->SetLHS(localX.getRawPtr());
    ssym_->OrigLP->SetRHS(temp.getRawPtr());
    ssym_->ReIdx_LP->fwd();

}

int ShyLU_Local_Schur_Operator::ApplyInverse(const Epetra_MultiVector &X,
             Epetra_MultiVector &Y) const
{
    return 0;
}

double ShyLU_Local_Schur_Operator::NormInf() const
{
   return -1;
}

const char* ShyLU_Local_Schur_Operator::Label() const
{
    return "Shylu Local Schur Operator (Wide Separator)";
}

bool ShyLU_Local_Schur_Operator::UseTranspose() const
{
    return false;
}

bool ShyLU_Local_Schur_Operator::HasNormInf() const
{
    return false;
}

const Epetra_Comm& ShyLU_Local_Schur_Operator::Comm() const
{
    return G_->Comm();
}

const Epetra_Map& ShyLU_Local_Schur_Operator::OperatorDomainMap() const
{
    //return C_->ColMap();
    return G_->DomainMap();
}

const Epetra_Map& ShyLU_Local_Schur_Operator::OperatorRangeMap() const
{
    //return R_->RowMap();
    return G_->RangeMap();
}
