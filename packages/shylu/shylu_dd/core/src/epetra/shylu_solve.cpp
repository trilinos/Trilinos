// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu_solve.cpp

    \brief Solve the system for a given r.h.s (approximately)

    \author Siva Rajamanickam

*/

#include "shylu_util.h"
#include "shylu.h"

#include <Ifpack.h>

static int shylu_dist_solve(
    shylu_symbolic *ssym,
    shylu_data *data,
    shylu_config *config,
    const Epetra_MultiVector& X,
    Epetra_MultiVector& Y
)
{
    int err = 0;
    AztecOO *solver = 0;
    assert(X.Map().SameAs(Y.Map()));
    //assert(X.Map().SameAs(A_->RowMap()));
    const Epetra_MultiVector *newX;
    newX = &X;
    //rd_->redistribute(X, newX);

    int nvectors = newX->NumVectors();

    // May have to use importer/exporter
    Epetra_Map BsMap(-1, data->Snr, data->SRowElems, 0, X.Comm());
    Epetra_Map BdMap(-1, data->Dnr, data->DRowElems, 0, X.Comm());

    Epetra_MultiVector Bs(BsMap, nvectors);
    Epetra_Import BsImporter(BsMap, newX->Map());

    assert(BsImporter.SourceMap().SameAs(newX->Map()));
    assert((newX->Map()).SameAs(BsImporter.SourceMap()));

    Bs.Import(*newX, BsImporter, Insert);
    Epetra_MultiVector Xs(BsMap, nvectors);

    Epetra_SerialComm LComm;        // Use Serial Comm for the local vectors.
    Epetra_Map LocalBdMap(-1, data->Dnr, data->DRowElems, 0, LComm);
    Epetra_MultiVector localrhs(LocalBdMap, nvectors);
    Epetra_MultiVector locallhs(LocalBdMap, nvectors);

    Epetra_MultiVector Z(BdMap, nvectors);

    Epetra_MultiVector Bd(BdMap, nvectors);
    Epetra_Import BdImporter(BdMap, newX->Map());
    assert(BdImporter.SourceMap().SameAs(newX->Map()));
    assert((newX->Map()).SameAs(BdImporter.SourceMap()));
    Bd.Import(*newX, BdImporter, Insert);

    int lda;
    double *values;
    err = Bd.ExtractView(&values, &lda);
    assert (err == 0);
    int nrows = ssym->C->RowMap().NumMyElements();

    // copy to local vector //TODO: OMP ?
    assert(lda == nrows);
    for (int v = 0; v < nvectors; v++)
    {
       for (int i = 0; i < nrows; i++)
       {
           err = localrhs.ReplaceMyValue(i, v, values[i+v*lda]);
           assert (err == 0);
       }
    }

    // TODO : Do we need to reset the lhs and rhs here ?
    if (config->amesosForDiagonal)
    {
        ssym->LP->SetRHS(&localrhs);
        ssym->LP->SetLHS(&locallhs);
        ssym->Solver->Solve();
    }
    else
    {
        ssym->ifSolver->ApplyInverse(localrhs, locallhs);
    }

    err = locallhs.ExtractView(&values, &lda);
    assert (err == 0);

    // copy to distributed vector //TODO: OMP ?
    assert(lda == nrows);
    for (int v = 0; v < nvectors; v++)
    {
       for (int i = 0; i < nrows; i++)
       {
           err = Z.ReplaceMyValue(i, v, values[i+v*lda]);
           assert (err == 0);
       }
    }

    Epetra_MultiVector temp1(BsMap, nvectors);
    ssym->R->Multiply(false, Z, temp1);
    Bs.Update(-1.0, temp1, 1.0);

    Xs.PutScalar(0.0);

    Epetra_LinearProblem Problem(data->Sbar.get(), &Xs, &Bs);
    if (config->schurSolver == "Amesos")
    {
        Amesos_BaseSolver *solver2 = data->dsolver;
        data->LP2->SetLHS(&Xs);
        data->LP2->SetRHS(&Bs);
        //cout << "Calling solve *****************************" << endl;
        solver2->Solve();
        //cout << "Out of solve *****************************" << endl;
    }
    else
    {
        if (config->libName == "Belos")
        {
            solver = data->innersolver;
            solver->SetLHS(&Xs);
            solver->SetRHS(&Bs);
        }
        else
        {
            // See the comment above on why we are not able to reuse the solver
            // when outer solve is AztecOO as well.
            solver = new AztecOO();
            //solver.SetPrecOperator(precop_);
            solver->SetAztecOption(AZ_solver, AZ_gmres);
            // Do not use AZ_none
            solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
            //solver->SetAztecOption(AZ_precond, AZ_none);
            //solver->SetAztecOption(AZ_precond, AZ_Jacobi);
            ////solver->SetAztecOption(AZ_precond, AZ_Neumann);
            //solver->SetAztecOption(AZ_overlap, 3);
            //solver->SetAztecOption(AZ_subdomain_solve, AZ_ilu);
            //solver->SetAztecOption(AZ_output, AZ_all);
            //solver->SetAztecOption(AZ_diagnostics, AZ_all);
            solver->SetProblem(Problem);
        }

        // What should be a good inner_tolerance :-) ?
        solver->Iterate(config->inner_maxiters, config->inner_tolerance);
    }

    Epetra_MultiVector temp(BdMap, nvectors);
    ssym->C->Multiply(false, Xs, temp);
    temp.Update(1.0, Bd, -1.0);

    //Epetra_SerialComm LComm;        // Use Serial Comm for the local vectors.
    //Epetra_Map LocalBdMap(-1, data->Dnr, data->DRowElems, 0, LComm);
    //Epetra_MultiVector localrhs(LocalBdMap, nvectors);
    //Epetra_MultiVector locallhs(LocalBdMap, nvectors);

    //int lda;
    //double *values;
    err = temp.ExtractView(&values, &lda);
    assert (err == 0);
    //int nrows = data->Cptr->RowMap().NumMyElements();

    // copy to local vector //TODO: OMP ?
    assert(lda == nrows);
    for (int v = 0; v < nvectors; v++)
    {
       for (int i = 0; i < nrows; i++)
       {
           err = localrhs.ReplaceMyValue(i, v, values[i+v*lda]);
           assert (err == 0);
       }
    }

    if (config->amesosForDiagonal)
    {
        ssym->LP->SetRHS(&localrhs);
        ssym->LP->SetLHS(&locallhs);
        ssym->Solver->Solve();
    }
    else
    {
        ssym->ifSolver->ApplyInverse(localrhs, locallhs);
    }

    err = locallhs.ExtractView(&values, &lda);
    assert (err == 0);

    // copy to distributed vector //TODO: OMP ?
    assert(lda == nrows);
    for (int v = 0; v < nvectors; v++)
    {
       for (int i = 0; i < nrows; i++)
       {
           err = temp.ReplaceMyValue(i, v, values[i+v*lda]);
           assert (err == 0);
       }
    }

    // For checking faults
    //if (NumApplyInverse_ == 5)  temp.ReplaceMyValue(0, 0, 0.0);

    Epetra_Export XdExporter(BdMap, Y.Map());
    Y.Export(temp, XdExporter, Insert);

    Epetra_Export XsExporter(BsMap, Y.Map());
    Y.Export(Xs, XsExporter, Insert);

    if (config->libName == "Belos" || config->schurSolver == "Amesos")
    {
        // clean up
    }
    else
    {
        delete solver;
    }
    return err;
}

static int shylu_local_solve(
    shylu_symbolic *ssym,
    shylu_data *data,
    shylu_config *config,
    const Epetra_MultiVector& X,
    Epetra_MultiVector& Y
)
{
    int err = 0;
#ifndef NDEBUG
    int nvectors = X.NumVectors();
    assert (nvectors == data->localrhs->NumVectors());
#endif // NDEBUG

    // Initialize the X vector for iterative solver
    data->Xs->PutScalar(0.0);

    // Get local portion of X
    data->localrhs->Import(X, *(data->BdImporter), Insert);

    // locallhs is z in paper
    if (config->amesosForDiagonal) {
        ssym->OrigLP->SetRHS((data->localrhs).getRawPtr());
        ssym->OrigLP->SetLHS((data->locallhs).getRawPtr());
        ssym->ReIdx_LP->fwd();
        ssym->Solver->Solve();
    }
    else {
        ssym->ifSolver->ApplyInverse(*(data->localrhs), *(data->locallhs));
    }

    err = ssym->R->Multiply(false, *(data->locallhs), *(data->temp1));
    assert (err == 0);

    // Export temp1 to a dist vector - temp2
    data->temp2->Import(*(data->temp1), *(data->DistImporter), Insert);

    //Epetra_MultiVector Bs(SMap, nvectors); // b_2 - R * z in ShyLU paper
    data->Bs->Import(X, *(data->BsImporter), Insert);
    data->Bs->Update(-1.0, *(data->temp2), 1.0);

    AztecOO *solver = 0;
    Epetra_LinearProblem Problem(data->Sbar.get(),
                             (data->Xs).getRawPtr(), (data->Bs).getRawPtr());
    if ((config->schurSolver == "G") || (config->schurSolver == "IQR"))
    {
        IFPACK_CHK_ERR(data->iqrSolver->Solve(*(data->schur_op),
                                 *(data->Bs), *(data->Xs)));
    }
    else if (config->schurSolver == "Amesos")
    {
        Amesos_BaseSolver *solver2 = data->dsolver;
        data->OrigLP2->SetLHS((data->Xs).getRawPtr());
        data->OrigLP2->SetRHS((data->Bs).getRawPtr());
        data->ReIdx_LP2->fwd();
        //cout << "Calling solve *****************************" << endl;
        solver2->Solve();
        //cout << "Out of solve *****************************" << endl;
    }
    else
    {
        if (config->libName == "Belos")
        {
            solver = data->innersolver;
            solver->SetLHS((data->Xs).getRawPtr());
            solver->SetRHS((data->Bs).getRawPtr());
        }
        else
        {
            // See the comment above on why we are not able to reuse the solver
            // when outer solve is AztecOO as well.
            solver = new AztecOO();
            //solver.SetPrecOperator(precop_);
            solver->SetAztecOption(AZ_solver, AZ_gmres);
            // Do not use AZ_none
            solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
            //solver->SetAztecOption(AZ_precond, AZ_none);
            //solver->SetAztecOption(AZ_precond, AZ_Jacobi);
            ////solver->SetAztecOption(AZ_precond, AZ_Neumann);
            //solver->SetAztecOption(AZ_overlap, 3);
            //solver->SetAztecOption(AZ_subdomain_solve, AZ_ilu);
            //solver->SetAztecOption(AZ_output, AZ_all);
            //solver->SetAztecOption(AZ_diagnostics, AZ_all);
            solver->SetProblem(Problem);
        }

        // What should be a good inner_tolerance :-) ?
        solver->Iterate(config->inner_maxiters, config->inner_tolerance);
    }

    // Import Xs locally
    data->LocalXs->Import(*(data->Xs), *(data->XsImporter), Insert);

    err = ssym->C->Multiply(false, *(data->LocalXs), *(data->temp3));
    assert (err == 0);
    data->temp3->Update(1.0, *(data->localrhs), -1.0);

    if (config->amesosForDiagonal) {
        ssym->OrigLP->SetRHS((data->temp3).getRawPtr());
        ssym->OrigLP->SetLHS((data->locallhs).getRawPtr());
        ssym->ReIdx_LP->fwd();
        ssym->Solver->Solve();
    }
    else {
        ssym->ifSolver->ApplyInverse(*(data->temp3), *(data->locallhs));
    }

    Y.Export(*(data->locallhs), *(data->XdExporter), Insert);
    Y.Export(*(data->LocalXs), *(data->XsExporter), Insert);

    if (config->libName == "Belos" || config->schurSolver == "Amesos")
    {
        // clean up
    }
    else
    {
        delete solver;
    }
    return err;
}

int shylu_solve(
    shylu_symbolic *ssym,
    shylu_data *data,
    shylu_config *config,
    const Epetra_MultiVector& X,
    Epetra_MultiVector& Y
)
{
    if (config->sep_type != 1)
        shylu_dist_solve(ssym, data, config, X, Y);
    else
        shylu_local_solve(ssym, data, config, X, Y);

    return 0;
}
