#include <cstdlib>
#include <iostream>
#include <vector>
#include <functional>

#include <mpi.h>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>

#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>
#include <Galeri_Utils.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>

//#include <Ifpack.h>

#include <gmres_tools.h>
#include <gmres.h>
#include <ProjectionPreconditioner.hpp>

namespace {
typedef double Scalar;
typedef std::vector<Scalar> LocalVector;
typedef std::vector<LocalVector> LocalMatrix;
typedef IQR::GMRESManager<Epetra_BlockMap, Epetra_MultiVector, LocalMatrix, LocalVector> MyGMRESManager;
typedef IQR::ProjectionPreconditioner<Epetra_BlockMap,
									  Epetra_Operator,
									  Epetra_MultiVector>
		MyProjectionPreconditioner;

struct IdPreconditioner
{
    void ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)
    {
        Y = X;
    }
};

std::function<int(const Epetra_Operator&, Epetra_MultiVector&, const Epetra_MultiVector&,
                  IdPreconditioner*, IdPreconditioner*, MyGMRESManager&, int &, Scalar &)>
        gmresSolve1 = &IQR::GMRES<Epetra_Operator, Epetra_MultiVector,
                                         IdPreconditioner, IdPreconditioner, MyGMRESManager,
                                         LocalVector, Scalar>;
std::function<int(const Epetra_Operator&, Epetra_MultiVector&, const Epetra_MultiVector&,
                  IdPreconditioner*, MyGMRESManager*, MyGMRESManager&, int &, Scalar &)>
        gmresSolve2 = &IQR::GMRES<Epetra_Operator, Epetra_MultiVector,
                                         IdPreconditioner, MyGMRESManager, MyGMRESManager,
                                         LocalVector, Scalar>;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    bool leader = ! Comm.MyPID();

    Teuchos::RCP<Epetra_Map> Map;
    Teuchos::RCP<Epetra_CrsMatrix> A;

    Teuchos::ParameterList GaleriList;
    int n = std::atoi(argv[1]);
    GaleriList.set("nx", n);
    GaleriList.set("ny", n);
    GaleriList.set("nz", n);

    Map = Teuchos::rcp(Galeri::CreateMap("Cartesian3D", Comm, GaleriList));
    if (leader) std::cout << "Created Epetra map." << std::endl;

    // Creates a diagonal matrix with 1's on the diagonal
    A   = Teuchos::rcp(Galeri::CreateCrsMatrix("Laplace3D", &*Map, GaleriList));

    Teuchos::RCP<Epetra_MultiVector> LHS(new Epetra_MultiVector(*Map, 1, true));
    Teuchos::RCP<Epetra_MultiVector> LHS2(new Epetra_MultiVector(*Map, 1, true));

    Teuchos::RCP<Epetra_MultiVector> RHS(new Epetra_MultiVector(*Map, 1, false));
    RHS->PutScalar(1.0);
    if (leader) std::cout << "Created linear system." << std::endl;

    EpetraExt::RowMatrixToMatlabFile("matrix.m", *A);
    EpetraExt::MultiVectorToMatlabFile("rhs.m", *RHS);
    if (leader) std::cout << "Linear system exported to Matlab format." << std::endl;

    // Preconditioner for A;
//     Ifpack precFactory;
//     const char precType[] = "ILU";
//     Teuchos::RCP<Ifpack_Preconditioner> prec(precFactory.Create(precType, &*A, 1));
//     prec->Initialize();
//     prec->Compute();
//     if (leader) std::cout << "Computed preconditioner: " << precType << "." << std::endl;

    // Let's try to solve it with Simone's GMRES

    int krylovDim = std::atoi(argv[2]);
    int maxIter1 = std::atoi(argv[3]);
    int maxIter2 = std::atoi(argv[4]);
    int inexactIter = std::atoi(argv[5]);
    Scalar tol = 1e-8;
    Scalar tolAchieved = tol;

    if (leader) {
        std::cout << "Attempting GMRES solve with parameters:"<< std::endl
                  << "    max iter (stage 1)     = " << maxIter1 << std::endl
                  << "    max iter (stage 2)     = " << maxIter2 << std::endl
                  << "    Krylov space dimension = " << krylovDim << std::endl
                  << "    tolerance    = " << tol << std::endl;
    }


    RHS->SetSeed(10);
    RHS->Random();

    // Computing the larger mode;
    Scalar rescaledLHS;
    for (int i(0); i < 0; ++i)
    {
        A->Apply(*RHS,*LHS); // x = A*b
        LHS->Norm2(&rescaledLHS);
        LHS->Scale(1./ rescaledLHS);
        *RHS = *LHS ; //  b = x/ ||x|| = Ab/ ||Ab||
    }

    // Testing ProjectionPreconditioner
    Teuchos::ParameterList projPL;
    projPL.set("relative subspace dimension", 0.01, "");
    MyProjectionPreconditioner projectionPreconditioner(A, A, projPL);
    projectionPreconditioner.Setup();
    projectionPreconditioner.ApplyInverse(*RHS, *LHS);

    // Stage one
    MyGMRESManager gmresManager1(A->Map(), maxIter1, false);
    IdPreconditioner L;
    IdPreconditioner M;

    tolAchieved = tol/10.;

    Teuchos::Time t1("time to compute prec");
    t1.start();
    int err1 = gmresSolve1(*A, *LHS, *RHS, &L, &M, gmresManager1, maxIter1, tolAchieved);
    t1.stop();

    double solutionNorm1;
    LHS->Norm2(&solutionNorm1);
    if (leader) {
        std::cout << "Stage 1" << std::endl;
        std::cout << "Number of GMRES iterations: " << maxIter1 << std::endl;
        std::cout << "Norm of solution: " << solutionNorm1 << std::endl;
        std::cout << "Time: " << t1.totalElapsedTime() << std::endl;
    }


    // Stage two
    tolAchieved = tol;
    RHS->SetSeed(10);
    RHS->Random();

    LHS->PutScalar(0.0);
    MyGMRESManager gmresManager2(A->Map(), krylovDim, false);
    t1.reset();t1.start();
    int err2 = gmresSolve1(*A, *LHS, *RHS, &L, &M, gmresManager2, maxIter2, tolAchieved);
    t1.stop();

    double solutionNorm2;
    LHS->Norm2(&solutionNorm2);
    if (leader) {
        std::cout << "Stage 2" << std::endl;
        std::cout << "Number of GMRES iterations: " << maxIter2 << std::endl;
        std::cout << "Norm of solution: " << solutionNorm2 << std::endl;
        std::cout << "Time: " << t1.totalElapsedTime() << std::endl;
    }


    // Stage three
    tolAchieved = tol;
    RHS->SetSeed(10);
    RHS->Random();

    LHS->PutScalar(0.0);

    gmresManager2.start();
    t1.reset(); t1.start();
    err2 = gmresSolve2(*A, *LHS, *RHS, &L, &gmresManager1, gmresManager2, maxIter2, tolAchieved);
    t1.stop();

    LHS->Norm2(&solutionNorm2);
    if (leader) {
        std::cout << "Stage 3" << std::endl;
        std::cout << "Number of GMRES iterations: " << maxIter2 << std::endl;
        std::cout << "Norm of solution: " << solutionNorm2 << std::endl;
        std::cout << "Time: " << t1.totalElapsedTime() << std::endl;
    }

    // Stage four, inexact solve
    tolAchieved = tol;
    RHS->SetSeed(10);
    RHS->Random();
    LHS2->PutScalar(0.0);

    gmresManager2.start();
    t1.reset(); t1.start();
//    gmresManager1.ApplyInverse(*RHS,*LHS2);
    err2 = gmresSolve2(*A, *LHS2, *RHS, &L, &gmresManager1, gmresManager2, inexactIter, tolAchieved);
    t1.stop();

    LHS2->Update(-1.0, *LHS, 1.0);
    LHS2->Norm2(&solutionNorm1);
    if (leader) {
        std::cout << "Stage 4, inexact GMRES solve, iterations : " << inexactIter <<  std::endl;
        std::cout << "Relative error: " << solutionNorm1/solutionNorm2 << std::endl;
        std::cout << "Time: " << t1.totalElapsedTime() << std::endl;
    }

    // Stage five, inexact solve
    tolAchieved = tol;
    RHS->SetSeed(10);
    RHS->Random();
    LHS2->PutScalar(0.0);

    gmresManager2.start();
    t1.reset(); t1.start();
    gmresManager1.ApplyInverse(*RHS,*LHS2);
    t1.stop();

    LHS2->Update(-1.0, *LHS, 1.0);
    LHS2->Norm2(&solutionNorm1);
    if (leader) {
        std::cout << "Stage 5, inexact solve" << std::endl;
        std::cout << "Relative error: " << solutionNorm1/solutionNorm2 << std::endl;
        std::cout << "Time: " << t1.totalElapsedTime() << std::endl;
    }


    MPI_Finalize();

    return 0;
}
