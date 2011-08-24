/** \file Ifpack_ShyLU.cpp

    \brief Use ShyLU as a preconditioner within IFPACK.

    \author Siva Rajamanickam

*/

#include "Ifpack_ShyLU.h"
#include "Teuchos_Time.hpp"

//#define DUMP_MATRICES

Ifpack_ShyLU::Ifpack_ShyLU(Epetra_CrsMatrix* A):
    A_(A),
    IsParallel_(true),
    IsInitialized_(false),
    IsComputed_(false),
    Label_(),
    NumApplyInverse_(0),
    Time_(A_->Comm())
{
}

void Ifpack_ShyLU::Destroy()
{
    if (IsInitialized_)
    {
        //delete A_;
        //delete partitioner_;
        //delete rd_;

        // I would rather explicitly delete
        /*delete slu_sym_.LP;
        delete slu_sym_.Solver;
        delete slu_sym_.C;
        delete slu_sym_.R;
        delete slu_sym_.D;
        delete slu_sym_.G;
        delete slu_sym_.Sg;
        delete slu_sym_.prober;*/

        if (slu_config_.schurSolver == "Amesos")
        {
            delete slu_data_.LP2;
            delete slu_data_.dsolver;
        }
        else if (slu_config_.libName == "Belos")
            delete slu_data_.innersolver;

        delete[] slu_data_.DRowElems;
        delete[] slu_data_.SRowElems;
        delete[] slu_data_.DColElems;
        delete[] slu_data_.gvals;
    }
    if (IsComputed_)
    {
        // I would rather explicitly delete
        // delete slu_data_.Sbar;
    }
}

int Ifpack_ShyLU::Initialize()
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
    shylu_factor(A_, 1, LP_, Solver_, C_, Dnr_, DRowElems_, Snr_, SRowElems_,
                    Sbar_, Sdiagfactor);*/
    slu_config_.sym =  Teuchos::getParameter<int>(List_,
                                                "Symmetry");
    slu_config_.libName = Teuchos::getParameter<string>(List_,
                                                "Outer Solver Library");
    string schurApproxMethod = Teuchos::getParameter<string>(List_,
                                                "Schur Approximation Method");
    slu_config_.schurSolver = Teuchos::getParameter<string>(List_,
                                                "Schur Complement Solver");
    slu_config_.relative_threshold =  0.0;
    slu_config_.Sdiagfactor =  0.05;
    if (schurApproxMethod == "A22AndBlockDiagonals")
    {
        slu_config_.schurApproxMethod = 1;
        slu_config_.Sdiagfactor =  Teuchos::getParameter<double>(List_,
                                                    "Diagonal Factor");
    }
    else if (schurApproxMethod == "Threshold")
    {
        slu_config_.schurApproxMethod = 2;
        slu_config_.relative_threshold =  Teuchos::getParameter<double>(List_,
                                                    "Relative Threshold");
    }

    slu_config_.inner_tolerance =  Teuchos::getParameter<double>(List_,
                                                "Inner Solver Tolerance");
    slu_config_.inner_maxiters =  Teuchos::getParameter<int>(List_,
                                                "Inner Solver MaxIters");
    string sep_type = Teuchos::getParameter<string>(List_,
                                                    "Separator Type");
    int dl =  List_.get<int>("Debug Level", 0);
    //slu_config_.dm.setDebugLevel(dl); 

    if (sep_type == "Wide")
        slu_config_.sep_type = 1;
    else
        slu_config_.sep_type = 2;

    shylu_symbolic_factor(A_, &slu_sym_, &slu_data_, &slu_config_);

    if (slu_config_.schurSolver == "Amesos")
    {
        slu_data_.LP2 = new Epetra_LinearProblem();
        Amesos Factory;
        std::string SolverType = "Amesos_Klu";
        bool IsAvailable = Factory.Query(SolverType);
        assert(IsAvailable == true);
        slu_data_.dsolver = Factory.Create(SolverType, *(slu_data_.LP2));
    }
    else
    {
        if (slu_config_.libName == "Belos")
            slu_data_.innersolver = new AztecOO() ;
        else
            slu_data_.innersolver = NULL;
    }

    IsInitialized_ = true;
    return 0;
}

int Ifpack_ShyLU::SetParameters(Teuchos::ParameterList& parameterlist)
{
    List_ = parameterlist;
    return 0;
}

int Ifpack_ShyLU::Compute()
{
    Teuchos::Time ftime("setup time");
    ftime.start();

    shylu_factor(A_, &slu_sym_, &slu_data_, &slu_config_);

    slu_config_.libName = Teuchos::getParameter<string>(List_,
                                                "Outer Solver Library");
    if (slu_config_.schurSolver == "Amesos")
    {
        Teuchos::ParameterList aList;
        aList.set("Reindex", true);

        slu_data_.LP2->SetOperator(slu_data_.Sbar.get());
        slu_data_.LP2->SetLHS(0); slu_data_.LP2->SetRHS(0);
        slu_data_.dsolver->SetParameters(aList);
        //cout << "Created the direct Schur  Solver" << endl;

        slu_data_.dsolver->SymbolicFactorization();
        //cout << "Symbolic Factorization done for schur complement" << endl;

        //cout << "In Numeric Factorization of Schur complement" << endl;
        slu_data_.dsolver->NumericFactorization();
        //cout << "Numeric Factorization done for schur complement" << endl;

    }
    else
    {
        if (slu_config_.libName == "Belos")
        {
            int err = slu_data_.innersolver->SetUserMatrix
                        (slu_data_.Sbar.get());
            assert (err == 0);
            slu_data_.innersolver->SetAztecOption(AZ_solver, AZ_gmres);
            slu_data_.innersolver->SetAztecOption(AZ_precond, AZ_dom_decomp);
            slu_data_.innersolver->SetAztecOption(AZ_keep_info, 1);
            //solver->SetAztecOption(AZ_overlap, 3);
            //solver->SetAztecOption(AZ_subdomain_solve, AZ_ilu);
            slu_data_.innersolver->SetMatrixName(999);

            double condest;
            err = slu_data_.innersolver->ConstructPreconditioner(condest);
            assert (err == 0);
            //cout << "Condition number of inner Sbar" << condest << endl;
        }
        else
        {
            // I suspect there is a bug in AztecOO. Doing what we do in the if
            // case
            // here will cause an error when we use the solver in ApplyInverse
            // The error will not happen when we call the dummy JustTryIt()
            // below
            slu_data_.innersolver = NULL;
        }
    }

    ftime.stop();
    //cout << "Time to ConstructPreconditioner" << ftime.totalElapsedTime() 
            //<< endl;
    IsComputed_ = true;
    //cout << " Done with the compute" << endl ;
    return 0;
}

int Ifpack_ShyLU::JustTryIt()
{
    // Dummy function, To show the error in AztecOO, This works
    //cout << "Entering JustTryIt" << endl;
    AztecOO *solver;
    solver = slu_data_.innersolver;
    //cout << solver_ << endl;
    Epetra_Map BsMap(-1, slu_data_.Snr, slu_data_.SRowElems, 0, A_->Comm());
    Epetra_MultiVector Xs(BsMap, 1);
    Epetra_MultiVector Bs(BsMap, 1);
    Xs.PutScalar(0.0);
    solver->SetLHS(&Xs);
    solver->SetRHS(&Bs);
    solver->Iterate(30, 1e-10);
    return 0;
}

int Ifpack_ShyLU::ApplyInverse(const Epetra_MultiVector& X,
    Epetra_MultiVector& Y) const
{
#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0)
    {
        EpetraExt::MultiVectorToMatlabFile("X.mat", X);
    }
#endif
    //cout << "Entering ApplyInvers" << endl;

    shylu_solve(&slu_sym_, &slu_data_, &slu_config_, X, Y);
#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0)
    {
        EpetraExt::MultiVectorToMatlabFile("Y.mat", Y);
    }
#endif
    NumApplyInverse_++;
    //cout << "Leaving ApplyInvers" << endl;
    return 0;
}

//! Computes the estimated condition number and returns the value.
double Ifpack_ShyLU::Condest(const Ifpack_CondestType CT,
     const int MaxIters, const double Tol, Epetra_RowMatrix* Matrix_in)
{
    return -1.0;
}

//! Prints on stream basic information about \c this object.
ostream& Ifpack_ShyLU::Print(ostream& os) const
{
    os << " !!!!!!!!! " << endl;
    return os;
}
