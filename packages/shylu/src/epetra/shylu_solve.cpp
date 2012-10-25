
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

/** \file shylu_solve.cpp

    \brief Solve the system for a given r.h.s (approximately)

    \author Siva Rajamanickam

*/

#include "shylu_util.h"
#include "shylu.h"

static int shylu_dist_solve(
    shylu_symbolic *ssym,
    shylu_data *data,
    shylu_config *config,
    const Epetra_MultiVector& X,
    Epetra_MultiVector& Y
)
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
    ssym->LP->SetRHS(&localrhs);
    ssym->LP->SetLHS(&locallhs);
    ssym->Solver->Solve();

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

    ssym->LP->SetRHS(&localrhs);
    ssym->LP->SetLHS(&locallhs);
    ssym->Solver->Solve();

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

static int shylu_local_solve(
    shylu_symbolic *ssym,
    shylu_data *data,
    shylu_config *config,
    const Epetra_MultiVector& X,
    Epetra_MultiVector& Y
)
{
    int err;
    int nvectors = X.NumVectors();
    Epetra_SerialComm LComm;        // Use Serial Comm for the local blocks.
    //cout <<" In local solve " << endl;

    Epetra_Map LocalDMap(-1, data->Dnr, data->DRowElems, 0, LComm);
    Epetra_Map LocalSMap(-1, data->Snr, data->SRowElems, 0, LComm);
    Epetra_Map SMap(-1, data->Snr, data->SRowElems, 0, X.Comm());

    Epetra_Import BdImporter(LocalDMap, X.Map()); // TODO: Construct only once

    // Get local portion of X
    Epetra_MultiVector localrhs(LocalDMap, nvectors);
    localrhs.Import(X, BdImporter, Insert);

    Epetra_MultiVector locallhs (View, *(ssym->Dlhs), 0,  nvectors); // z in
                                                                    // paper
    Epetra_MultiVector temp3 (View, *(ssym->Drhs), 0,  nvectors);

    ssym->OrigLP->SetRHS(&localrhs);
    ssym->OrigLP->SetLHS(&locallhs);
    ssym->ReIdx_LP->fwd();
    ssym->Solver->Solve();

    Epetra_MultiVector temp1(LocalSMap, nvectors);
    err = ssym->R->Multiply(false, locallhs, temp1);
    assert (err == 0);

    // Export temp1 to a dist vector - temp2
    Epetra_MultiVector temp2(SMap, nvectors);
    Epetra_Import DistImporter(SMap, LocalSMap);
    //temp2.Import(X, DistImporter, Insert);
    temp2.Import(temp1, DistImporter, Insert);

    Epetra_MultiVector Bs(SMap, nvectors); // b_2 - R * z in ShyLU paper
    Epetra_MultiVector Xs(SMap, nvectors);
    Epetra_Import BsImporter(SMap, X.Map());
    Bs.Import(X, BsImporter, Insert);

    Bs.Update(-1.0, temp2, 1.0);

    AztecOO *solver;
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

    // Import Xs locally
    Epetra_MultiVector LocalXs(LocalSMap, nvectors);
    Epetra_Import XsImporter(LocalSMap, SMap);
    LocalXs.Import(Xs, XsImporter, Insert);

    err = ssym->C->Multiply(false, LocalXs, temp3);
    assert (err == 0);
    temp3.Update(1.0, localrhs, -1.0);

    ssym->OrigLP->SetRHS(&temp3);
    ssym->OrigLP->SetLHS(&locallhs);
    ssym->ReIdx_LP->fwd();
    ssym->Solver->Solve();

    Epetra_Export XdExporter(LocalDMap, Y.Map());
    Y.Export(locallhs, XdExporter, Insert);

    Epetra_Export XsExporter2(LocalSMap, Y.Map());
    Y.Export(LocalXs, XsExporter2, Insert);

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
}
