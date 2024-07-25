// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file Ifpack_ShyLU.cpp

    \brief Use ShyLU as a preconditioner within IFPACK.

    \author Siva Rajamanickam

*/

#include "Ifpack_ShyLU.h"
#include "Teuchos_Time.hpp"
#include "Ifpack_config.h"

#include <IQRSolver.h>

//#define DUMP_MATRICES

Ifpack_ShyLU::Ifpack_ShyLU(Epetra_CrsMatrix* A):
    A_(A),
    IsParallel_(true),
    IsInitialized_(false),
    IsComputed_(false),
    Label_("Ifpack_ShyLU"),
    NumCompute_(0),
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
        slu_data_.localSbargraph = Teuchos::null;
        slu_data_.guided_prober = Teuchos::null;
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
    slu_config_.libName = Teuchos::getParameter<std::string>(List_,
                                                "Outer Solver Library");
    std::string schurApproxMethod = Teuchos::getParameter<std::string>(List_,
                                                "Schur Approximation Method");
    slu_config_.schurSolver = Teuchos::getParameter<std::string>(List_,
                                                "Schur Complement Solver");
    slu_config_.schurAmesosSolver = List_.get<std::string>("Schur Amesos Solver",
                                                "Amesos_Klu");
    //slu_config_.diagonalBlockSolver =
            //Teuchos::getParameter<std::string>(List_, "Diagonal Block Solver");
    slu_config_.diagonalBlockSolver = List_.get<std::string>("Diagonal Block Solver",
                                                    "Amesos_Klu");


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
    else if (schurApproxMethod == "Guided Probing")
    {
        slu_config_.schurApproxMethod = 3;
        slu_config_.reset_iter =  List_.get<int>("Schur Recompute Iteration", 10);
        slu_config_.relative_threshold =  Teuchos::getParameter<double>(List_,
                                                    "Relative Threshold");
    }
    if (schurApproxMethod == "IQR")
    {
        slu_config_.schurSolver = "IQR";
        slu_config_.schurApproxMethod = 4;

        List_.set<bool>("Use full IQR", true, "");
        slu_data_.iqrSolver.reset(new IQR::IQRSolver(List_));
    }
    if (schurApproxMethod == "G")
    {
        slu_config_.schurSolver = "G";
        slu_config_.schurApproxMethod = 5;

        List_.set<bool>("Use full IQR", false, "");
        slu_data_.iqrSolver.reset(new IQR::IQRSolver(List_));
    }

    slu_config_.schurPreconditioner = List_.get<std::string>("Schur Preconditioner",
                                                        "ILU stand-alone");
    slu_config_.silent_subiter = List_.get<bool>("Silent subiterations",
                                                 true);

    slu_config_.inner_tolerance =  Teuchos::getParameter<double>(List_,
                                                "Inner Solver Tolerance");
    slu_config_.inner_maxiters =  Teuchos::getParameter<int>(List_,
                                                "Inner Solver MaxIters");
    slu_config_.overlap =  List_.get<int>("Schur Preconditioner Overlap", 0);
    std::string sep_type = Teuchos::getParameter<std::string>(List_,
                                                    "Separator Type");
    // mfh 25 May 2015: dl is unused (don't declare it, or you'll get
    // build warnings).  However, two-argument get() has side effects,
    // so I don't want to comment it out.
    //
    //int dl =  List_.get<int>("Debug Level", 0);
    (void) List_.get<int>("Debug Level", 0);
    //slu_config_.dm.setDebugLevel(dl);

    if (sep_type == "Wide")
        slu_config_.sep_type = 1;
    else
        slu_config_.sep_type = 2;

    shylu_symbolic_factor(A_, &slu_sym_, &slu_data_, &slu_config_);

    if (slu_config_.schurSolver == "Amesos")
    {
        slu_data_.LP2 = Teuchos::rcp(new Epetra_LinearProblem());
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
#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
    List_ = parameterlist.sublist("ShyLU list", true);
#else
    List_ = parameterlist;
#endif
    return 0;
}

int Ifpack_ShyLU::Compute()
{
    Teuchos::Time ftime("setup time");
    ftime.start();

    slu_data_.num_compute = NumCompute_;

    shylu_factor(A_, &slu_sym_, &slu_data_, &slu_config_);

    ftime.stop();
    IsComputed_ = true;
    NumCompute_++;
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
std::ostream& Ifpack_ShyLU::Print(std::ostream& os) const
{
    os << " !!!!!!!!! " << std::endl;
    return os;
}
