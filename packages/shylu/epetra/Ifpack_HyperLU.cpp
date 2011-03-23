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
        delete hlu_data_.LP;
        delete hlu_data_.Solver;
        delete hlu_data_.innersolver;
        delete[] hlu_data_.DRowElems;
        delete[] hlu_data_.SRowElems;
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

    /*double Sdiagfactor = 0.05; // hard code the diagonals
    HyperLU_factor(A_, 1, LP_, Solver_, C_, Dnr_, DRowElems_, Snr_, SRowElems_,
                    Sbar_, Sdiagfactor);*/
    hlu_config_.sym =  Teuchos::getParameter<int>(List_,
                                                "Symmetry");
    hlu_config_.libName = Teuchos::getParameter<string>(List_,
                                                "Outer Solver Library");
    string schurApproxMethod = Teuchos::getParameter<string>(List_,
                                                "Schur Approximation Method");
    hlu_config_.relative_threshold =  0.0;
    hlu_config_.Sdiagfactor =  0.05;
    if (schurApproxMethod == "A22AndBlockDiagonals")
    {
        hlu_config_.schurApproxMethod = 1;
        hlu_config_.Sdiagfactor =  Teuchos::getParameter<double>(List_,
                                                    "Diagonal Factor");
    }
    else if (schurApproxMethod == "Threshold")
    {
        hlu_config_.schurApproxMethod = 2;
        hlu_config_.relative_threshold =  Teuchos::getParameter<double>(List_,
                                                    "Relative Threshold");
    }

    hlu_config_.inner_tolerance =  Teuchos::getParameter<double>(List_,
                                                "Inner Solver Tolerance");
    hlu_config_.inner_maxiters =  Teuchos::getParameter<int>(List_,
                                                "Inner Solver MaxIters");
    HyperLU_factor(A_, &hlu_data_, &hlu_config_);

    IsInitialized_ = true;
    return 0;
}

int Ifpack_HyperLU::SetParameters(Teuchos::ParameterList& parameterlist)
{
    List_ = parameterlist;
}

int Ifpack_HyperLU::Compute()
{
    AztecOO *solver;
    Teuchos::Time ftime("setup time");
    ftime.start();

    hlu_config_.libName = Teuchos::getParameter<string>(List_,
                                                "Outer Solver Library");
    if (hlu_config_.libName == "Belos")
    {
        solver  = new AztecOO() ;
        int err = solver->SetUserMatrix(hlu_data_.Sbar.get());
        assert (err == 0);
        solver->SetAztecOption(AZ_solver, AZ_gmres);
        solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
        solver->SetAztecOption(AZ_keep_info, 1);
        solver->SetMatrixName(999);

        double condest;
        err = solver->ConstructPreconditioner(condest);
        assert (err == 0);
        //cout << "Condition number of inner Sbar" << condest << endl;
    }
    else
    {
        // I suspect there is a bug in AztecOO. Doing what we do in the if case
        // here will cause an error when we use the solver in ApplyInverse
        // The error will not happen when we call the dummy JustTryIt() below
        solver = NULL;
    }

    ftime.stop();
    //cout << "Time to ConstructPreconditioner" << ftime.totalElapsedTime() 
            //<< endl;
    hlu_data_.innersolver = solver;
    IsComputed_ = true;
    //cout << " Done with the compute" << endl ;
    return 0;
}

int Ifpack_HyperLU::JustTryIt()
{
    // Dummy function, To show the error in AztecOO, This works
    //cout << "Entering JustTryIt" << endl;
    AztecOO *solver;
    solver = hlu_data_.innersolver;
    //cout << solver_ << endl;
    Epetra_Map BsMap(-1, hlu_data_.Snr, hlu_data_.SRowElems, 0, A_->Comm());
    Epetra_MultiVector Xs(BsMap, 1);
    Epetra_MultiVector Bs(BsMap, 1);
    Xs.PutScalar(0.0);
    solver->SetLHS(&Xs);
    solver->SetRHS(&Bs);
    solver->Iterate(30, 1e-10);
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
    //cout << "Entering ApplyInvers" << endl;

    hyperlu_solve(&hlu_data_, &hlu_config_, X, Y);
#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0)
    {
        EpetraExt::MultiVectorToMatlabFile("Y.mat", Y);
        cout << Y;
    }
#endif
    NumApplyInverse_++;
    //cout << "Leaving ApplyInvers" << endl;
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
