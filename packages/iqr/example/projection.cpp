#include <cstdlib>
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

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
                  IdPreconditioner*, MyProjectionPreconditioner*, MyGMRESManager&, int &, Scalar &)>
        gmresSolveProj = &IQR::GMRES<Epetra_Operator, Epetra_MultiVector,
                                         IdPreconditioner, MyProjectionPreconditioner, MyGMRESManager,
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
    Teuchos::RCP<Epetra_MultiVector> RHS(new Epetra_MultiVector(*Map, 1, false));
    RHS->SetSeed(10);
    RHS->Random();

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
    double subspaceDim = std::atof(argv[2]);
    int maxIter = std::atoi(argv[3]);

    Scalar tol = 1e-7;
    Scalar tolAchieved = tol;

    if (leader) {
        std::cout << "Attempting GMRES solve with parameters:"<< std::endl
                  << "  Max iter: " << maxIter << std::endl
                  << "  Projection subspace dimension: "
                  << static_cast<int>(std::floor(subspaceDim * n * n * n))
                  << std::endl
                  << "  Tolerance: " << tol << std::endl;
    }


    // Testing ProjectionPreconditioner
    IdPreconditioner L;
    Teuchos::ParameterList projPL;
    projPL.set("relative subspace dimension", subspaceDim, "");
    MyProjectionPreconditioner projectionPreconditioner(A, A, projPL);
    projectionPreconditioner.Setup();
    MyGMRESManager gmresManagerP(A->Map(), maxIter, false);
    int errP = gmresSolveProj(*A, *LHS, *RHS, &L, &projectionPreconditioner, gmresManagerP, maxIter, tolAchieved);

    if (leader) {
        std::cout << "Achieved tolerance: " << tolAchieved << std::endl
                  << "Number of iterations: " << maxIter << std::endl;
    }

    MPI_Finalize();

    return 0;
}
