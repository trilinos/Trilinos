/** \file shylu_solve.cpp

    \brief Solve the system for a given r.h.s (approximately)

    \author Siva Rajamanickam

*/

#include "shylu_util.h"
#include "shylu.h"

int shylu_solve(shylu_data *data, shylu_config *config,
    const Epetra_MultiVector& X, Epetra_MultiVector& Y)
{
    int err;
    AztecOO *solver;
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
    int nrows = data->Cptr->RowMap().NumMyElements();

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

    data->LP->SetRHS(&localrhs);
    data->LP->SetLHS(&locallhs);
    data->Solver->Solve();

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
    data->Rptr->Multiply(false, Z, temp1);
    Bs.Update(-1.0, temp1, 1.0);

    Xs.PutScalar(0.0);

    Epetra_LinearProblem Problem(data->Sbar.get(), &Xs, &Bs);
    if (config->schurSolver == "Amesos")
    {
        Amesos_BaseSolver *solver = data->dsolver;
        data->LP2->SetLHS(&Xs);
        data->LP2->SetRHS(&Bs);
        //cout << "Calling solve *****************************" << endl;
        solver->Solve();
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
    data->Cptr->Multiply(false, Xs, temp);
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

    data->LP->SetRHS(&localrhs);
    data->LP->SetLHS(&locallhs);
    data->Solver->Solve();

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
    return 0;
}
