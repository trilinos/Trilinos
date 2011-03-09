/** \file Ifpack_HyperLU.cpp

    \brief Use HyperLU as a preconditioner within IFPACK.

    \author Siva Rajamanickam

*/

#include "hyperlu.h"
#include "Ifpack_HyperLU.h"
#include "Teuchos_Time.hpp"

//#define DUMP_MATRICES

Ifpack_HyperLU::Ifpack_HyperLU(Epetra_CrsMatrix* A):
    A_(A),
    IsParallel_(true),
    IsInitialized_(false),
    IsComputed_(false),
    Label_(),
    NumApplyInverse_(0),
    Time_(A_->Comm())
{
}

void Ifpack_HyperLU::Destroy()
{
    if (IsInitialized_)
    {
        //delete A_;
        //delete partitioner_;
        //delete rd_;
    }
    if (IsComputed_)
    {
        delete LP_;
        delete Solver_;
        delete[] DRowElems_;
        delete[] SRowElems_;
    }
}

int Ifpack_HyperLU::Initialize()
{
    if(Comm().NumProc() != 1) 
        IsParallel_ = true;
    else 
        IsParallel_ = false;

    // TODO:
    // Need to enable partitioning here in Initialize once moving to Belos

    // Cannot call this method , need the partitioner around TODO : Can we 
    // avoid this
    //A_ = balanceAndRedistribute(A_, List_);

    // ==================== Symbolic factorization =========================
    // 1. Partition and redistribute [
    //partitioner_ = new Isorropia::Epetra::Partitioner(A_, List_, false);
    //partitioner_->partition();

    //rd_ = new Isorropia::Epetra::Redistributor(partitioner_);
    //Epetra_CrsMatrix *newA;
    //rd_->redistribute(*A_, newA);
    //A_ = newA;
    // ]

    double Sdiagfactor = 0.05; // hard code the diagonals
    HyperLU_factor(A_, 1, LP_, Solver_, C_, Dnr_, DRowElems_, Snr_, SRowElems_,
                    Sbar_, Sdiagfactor);

    IsInitialized_ = true;
    return 0;
}

int Ifpack_HyperLU::SetParameters(Teuchos::ParameterList& parameterlist)
{
    List_ = parameterlist;
}

int Ifpack_HyperLU::Compute()
{
    Teuchos::Time ftime("setup time");
    ftime.start();

    solver_  = new AztecOO() ;
    int err = solver_->SetUserMatrix(Sbar_.get());
    assert (err == 0);
    //err = solver_->SetPrecMatrix(Sbar_.get());
    //assert (err == 0);
    solver_->SetAztecOption(AZ_solver, AZ_gmres);
    solver_->SetAztecOption(AZ_precond, AZ_dom_decomp);

    solver_->SetAztecOption(AZ_keep_info, 1);
    double condest;
    err = solver_->ConstructPreconditioner(condest);
    assert (err == 0);
    cout << "Condition number of inner Sbar" << condest << endl;

    // Do a dummy iterate to construct the preconditioner
    Epetra_Map BsMap(-1, Snr_, SRowElems_, 0, A_->Comm());
    Xs_ = new Epetra_MultiVector(BsMap, 1);
    Bs_ = new Epetra_MultiVector(BsMap, 1);
    Bs_->PutScalar(1.0);
    solver_->SetLHS(Xs_);
    solver_->SetRHS(Bs_);
    solver_->Iterate(30, 1e-10);

    /* solver_->SetAztecOption(AZ_pre_calc, AZ_reuse);*/

    ftime.stop();
    cout << "Time to ConstructPreconditioner" << ftime.totalElapsedTime() 
            << endl;
    IsComputed_ = true;
    cout << solver_ << endl;
    cout << " Done with the compute" << endl ;
    return 0;
}

int Ifpack_HyperLU::JustTryIt()
{
    cout << "Entering JustTryIt" << endl;
    cout << solver_ << endl;
    solver_->Iterate(30, 1e-10);
}

int Ifpack_HyperLU::ApplyInverse(const Epetra_MultiVector& X, 
    Epetra_MultiVector& Y) const
{
#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0)
    {
        EpetraExt::MultiVectorToMatlabFile("X.mat", X);
        cout << X;
    }
#endif
    cout << "Entering ApplyInvers" << endl;
    int err;
    //if(NumApplyInverse_ == 0) cout << X;
    assert(X.Map().SameAs(Y.Map()));
    assert(X.Map().SameAs(A_->RowMap()));
    const Epetra_MultiVector *newX; 
    newX = &X;
    //rd_->redistribute(X, newX);

    int nvectors = newX->NumVectors();

    // May have to use importer/exporter
    Epetra_Map BsMap(-1, Snr_, SRowElems_, 0, X.Comm());
    Epetra_Map BdMap(-1, Dnr_, DRowElems_, 0, X.Comm());

    Epetra_MultiVector Bs(BsMap, nvectors);
    Epetra_Import BsImporter(BsMap, newX->Map());

    assert(BsImporter.SourceMap().SameAs(newX->Map()));
    assert((newX->Map()).SameAs(BsImporter.SourceMap()));

    Bs.Import(*newX, BsImporter, Insert);
    Epetra_MultiVector Xs(BsMap, nvectors);
    Xs.PutScalar(0.0);
    //cout << " Done with first import " << endl;
    //if(NumApplyInverse_ == 0) cout << Bs;

#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0 && A_->Comm().MyPID() == 1)
    {
        EpetraExt::MultiVectorToMatlabFile("localXs.mat", localXs);
    }
#endif

    Epetra_LinearProblem Problem(Sbar_.get(), &Xs, &Bs);
    //assert(solver_ != NULL);


#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0 )
    {
        //cout << *(Sbar_.get()) << endl;
        EpetraExt::RowMatrixToMatlabFile("Sbar.mat", *(Sbar_.get()));
    }
#endif

    // TODO: Can I create the solver in setup ?
    AztecOO solver;
    //solver.SetPrecOperator(precop_);
    solver.SetAztecOption(AZ_solver, AZ_gmres);
    // Do not use AZ_none
    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    //solver.SetAztecOption(AZ_precond, AZ_none);
    //solver.SetAztecOption(AZ_precond, AZ_Jacobi);
    ////solver.SetAztecOption(AZ_precond, AZ_Neumann);
    //solver.SetAztecOption(AZ_overlap, 3);
    //solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
    //solver.SetAztecOption(AZ_output, AZ_all);
    //solver.SetAztecOption(AZ_diagnostics, AZ_all);
    solver.SetProblem(Problem);

    //solver_->SetLHS(&Xs);
    //solver_->SetRHS(&Bs);
    //err = solver_->CheckInput();
    //assert (err == 0);

    /*const int *az_options = solver_->GetAllAztecOptions();
    for (int az_i = 0; az_i < AZ_OPTIONS_SIZE; az_i++)
        cout << az_options[az_i] << " ";
    cout << endl;*/

    solver.Iterate(30, 1e-10);

    //if(NumApplyInverse_ == 0) cout << "X vector from inner iteration" << Xs;

    Epetra_MultiVector Bd(BdMap, nvectors);
    Epetra_Import BdImporter(BdMap, newX->Map());
    assert(BdImporter.SourceMap().SameAs(newX->Map()));
    assert((newX->Map()).SameAs(BdImporter.SourceMap()));
    Bd.Import(*newX, BdImporter, Insert);

    Epetra_MultiVector temp(BdMap, nvectors);
    C_->Multiply(false, Xs, temp);
    temp.Update(1.0, Bd, -1.0);

    Epetra_SerialComm LComm;        // Use Serial Comm for the local vectors.
    Epetra_Map LocalBdMap(-1, Dnr_, DRowElems_, 0, LComm);
    Epetra_MultiVector localrhs(LocalBdMap, nvectors);
    Epetra_MultiVector locallhs(LocalBdMap, nvectors);

    int lda;
    double *values;
    err = temp.ExtractView(&values, &lda);
    assert (err == 0);
    int nrows = C_->RowMap().NumMyElements();

    // copy to local vector
    assert(lda == nrows);
    for (int v = 0; v < nvectors; v++)
    {
       for (int i = 0; i < nrows; i++)
       {
           err = localrhs.ReplaceMyValue(i, v, values[i+v*lda]);
           assert (err == 0);
       }
    }

    LP_->SetRHS(&localrhs);
    LP_->SetLHS(&locallhs);
    Solver_->Solve();

    err = locallhs.ExtractView(&values, &lda);
    assert (err == 0);

    // copy to distributed vector
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


#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0)
    {
        EpetraExt::MultiVectorToMatlabFile("Y.mat", Y);
        cout << Y;
    }
#endif
    NumApplyInverse_++;
    //delete newX;
    cout << "Leaving ApplyInvers" << endl;
    return 0;
}

//! Computes the estimated condition number and returns the value.
double Ifpack_HyperLU::Condest(const Ifpack_CondestType CT, 
     const int MaxIters, const double Tol, Epetra_RowMatrix* Matrix_in)
{
    return -1.0;
}

//! Prints on stream basic information about \c this object.
ostream& Ifpack_HyperLU::Print(ostream& os) const
{
    os << " !!!!!!!!! " << endl;
    return os;
}
