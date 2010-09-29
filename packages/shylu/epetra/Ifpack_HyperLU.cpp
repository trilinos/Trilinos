/** \file Ifpack_HyperLU.cpp

    \brief Use HyperLU as a preconditioner within IFPACK.

    \author Siva Rajamanickam

*/

#include "hyperlu.h"
#include "Ifpack_HyperLU.h"

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
        delete partitioner_;
        delete rd_;
    }
    if (IsComputed_)
    {
        delete localS_;
        delete LP_;
        delete Solver_;
        delete[] DRowElems_;
        delete[] SRowElems_;
        delete[] piv_;
        delete CMV_;
    }
}

int Ifpack_HyperLU::Initialize()
{
    if(Comm().NumProc() != 1) 
        IsParallel_ = true;
    else 
        IsParallel_ = false;

    // Cannot call this method , need the partitioner around TODO : Can we 
    // avoid this
    //A_ = balanceAndRedistribute(A_, List_);

    // ==================== Symbolic factorization =========================
    // 1. Partition and redistribute [
    partitioner_ = new Isorropia::Epetra::Partitioner(A_, List_, false);
    partitioner_->partition();

    rd_ = new Isorropia::Epetra::Redistributor(partitioner_);
    //Epetra_CrsMatrix *newA;
    //rd_->redistribute(*A_, newA);
    //A_ = newA;
    // ]

    IsInitialized_ = true;
    return 0;
}

int Ifpack_HyperLU::SetParameters(Teuchos::ParameterList& parameterlist)
{
    List_ = parameterlist;
}

int Ifpack_HyperLU::Compute()
{
    HyperLU_factor(A_, 1, localS_, LP_, Solver_, Dnr_, DRowElems_, Snr_,
               SRowElems_, piv_, CMV_);
    IsComputed_ = true;
    return 0;
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
    assert(X.Map().SameAs(Y.Map()));
    assert(X.Map().SameAs(A_->RowMap()));
    //cout << "Entering ApplyInvers" << endl;
    const Epetra_MultiVector *newX; 
    newX = &X;
    //rd_->redistribute(X, newX);

    int nvectors = newX->NumVectors();

    // May have to use importer/exporter
    Epetra_SerialComm LComm;        // Use Serial Comm for the local vectors.
    Epetra_Map localXsMap(-1, Snr_, SRowElems_, 0, LComm);
    Epetra_Map localXdMap(-1, Dnr_, DRowElems_, 0, LComm);

    Epetra_MultiVector localXs(localXsMap, nvectors);
    Epetra_Import XsImporter(localXsMap, newX->Map());
    assert(XsImporter.SourceMap().SameAs(newX->Map()));
    assert((newX->Map()).SameAs(XsImporter.SourceMap()));
    localXs.Import(*newX, XsImporter, Insert);
    //cout << " Done with first import " << endl;

    Teuchos::LAPACK<int, double> lapack;
    int info;
    // TODO : Need to check Snr for non blk diagonal Si.
    lapack.GETRS('N', Snr_, nvectors, localS_->Values(), Snr_, piv_,
                localXs.Values(), Snr_, &info);

#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0 && A_->Comm().MyPID() == 1)
    {
        EpetraExt::MultiVectorToMatlabFile("localXs.mat", localXs);
    }
#endif

    Epetra_MultiVector temp(localXdMap, nvectors);
    temp.Multiply('N', 'N', 1.0, *CMV_, localXs, 0.0);

    Epetra_MultiVector localXd(localXdMap, nvectors);
    Epetra_Import XdImporter(localXdMap, newX->Map());
    localXd.Import(*newX, XdImporter, Insert);
    assert(XdImporter.SourceMap().SameAs(newX->Map()));
    assert((newX->Map()).SameAs(XdImporter.SourceMap()));

    temp.Update(1.0, localXd, -1.0);

    LP_->SetRHS(&temp);
    LP_->SetLHS(&localXd);
    Solver_->Solve();

#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0 && A_->Comm().MyPID() == 1)
    {
        EpetraExt::MultiVectorToMatlabFile("localXd.mat", localXd);
    }
#endif

    Epetra_Export XdExporter(localXdMap, Y.Map());
    Y.Export(localXd, XdExporter, Insert);

    Epetra_Export XsExporter(localXsMap, Y.Map());
    Y.Export(localXs, XsExporter, Insert);

#ifdef DUMP_MATRICES
    if (NumApplyInverse_ == 0)
    {
        EpetraExt::MultiVectorToMatlabFile("Y.mat", Y);
        cout << Y;
    }
#endif
    NumApplyInverse_++;
    //delete newX;
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
