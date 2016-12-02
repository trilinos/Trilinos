//GPU_VERS should be 1 if you are *not* using a GPU 
//and 0 if you are. 
//REORDER controls whether the reordering code will be 
//compiled or not. This currently does not work on 
//GPU platforms.
#define __GPU_VERS__ 1
#define __REORDER__ 0

#include <cstdlib>
#include <iostream>
#include <random>
#include <algorithm>
#include <utility>
#include <vector>
#include <string>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <Tpetra_MatrixIO.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Zoltan2_OrderingProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Ifpack2_ReorderFilter.hpp>
#include <Ifpack2_RILUK.hpp>


#include "fast_ilu.hpp"
#include "Ifpack2_filu.hpp"
#include "fast_ic.hpp"
#include "Ifpack2_fic.hpp"
#include "fast_ldl.hpp"
#include "Ifpack2_fldl.hpp"

//typedef Kokkos::OpenMP MyExecSpace;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ScalarTraits;
using Tpetra::MultiVector;

#if __GPU_VERS__
typedef Kokkos::OpenMP MyExecSpace;
#else 
typedef Kokkos::Cuda MyExecSpace;
#endif

typedef double Scalar;
typedef int Ordinal;

typedef Kokkos::View<Ordinal *, MyExecSpace> OrdinalArray;
typedef Kokkos::View<const Ordinal *, MyExecSpace> ConstOrdinalArray;
typedef Kokkos::View<Scalar *, MyExecSpace> ScalarArray;
typedef Kokkos::CrsMatrix<Scalar, Ordinal, MyExecSpace> KCrsMatrix;

typedef Teuchos::ScalarTraits<Scalar>::magnitudeType MagnitudeType;

typedef Tpetra::CrsMatrix<Scalar, Ordinal, Ordinal> TCrsMatrix;
typedef Tpetra::Map<>::node_type Node;
typedef Tpetra::Operator<Scalar, Ordinal> OperatorType;
typedef Tpetra::MultiVector<Scalar, Ordinal> TMultiVector;

typedef Belos::MultiVecTraits<Scalar, MultiVector<Scalar, Ordinal> > BMultiVecTraits;



#if __GPU_VERS__
typedef Zoltan2::XpetraCrsMatrixAdapter<TCrsMatrix> ZMatrixAdapter;
#endif

using std::pair;
using std::vector;

void print_usage(char **argv)
{
    std::cout << argv[0] << " <Matrix Market filename> <factor sweeps> <trisol sweeps> <ic flag> <level> <ex-flag> <solver-flag> <omega> <shift> <guess-flag> <blkSz>" << std::endl;
}

#if __GPU_VERS__
bool compPair(pair<Ordinal, Scalar> p1, pair<Ordinal, Scalar> p2)
{
    return (p1.first < p2.first);
}

void reorderAndCopy(KCrsMatrix::row_map_type &aRowMap, OrdinalArray &aColIdx, ScalarArray &aVal, Ordinal *perm, 
        Ordinal *permInv, Ordinal nRows, OrdinalArray &aRowMapO, OrdinalArray &aColIdxO, ScalarArray &aValO)
{
    vector<pair<Ordinal, Scalar>> newRowCopy(nRows);
    aRowMapO[0] = 0;

    for(Ordinal row = 0; row < nRows; row++) 
    {
        Ordinal newRow = perm[row];
        aRowMapO[row + 1] = aRowMapO[row] + (aRowMap[newRow + 1] - aRowMap[newRow]);
        //Now copy newRow to position of row and re-index the columns
        Ordinal ptr = 0;
        for (Ordinal k = aRowMap[newRow]; k < aRowMap[newRow + 1]; k++) 
        {
            newRowCopy[ptr].first = permInv[aColIdx[k]];
            newRowCopy[ptr].second = aVal[k];
            ptr++;
        }
        std::sort(newRowCopy.begin(), newRowCopy.begin() + (aRowMap[newRow+1] - aRowMap[newRow]) , compPair);
        ptr = 0;
        for (Ordinal k = aRowMapO[row]; k < aRowMapO[row + 1]; k++) 
        {
            aColIdxO[k] = newRowCopy[ptr].first;
            aValO[k] = newRowCopy[ptr].second;
            ptr++;
        }
    }
    assert(aRowMapO[nRows] == aRowMap[nRows]);
}

void createReorderedTpetraMatrix(RCP<TCrsMatrix> &aIn, RCP<TCrsMatrix> &aOut, Ordinal *perm, Ordinal *permInv, 
        RCP<Node> oNode)
{

    typedef Tpetra::Map<Ordinal, Ordinal, Node> OMapType;
    Teuchos::ArrayRCP<size_t> rowPtr(aIn->getNodeNumRows() + 1, 0);
    Teuchos::ArrayRCP<Ordinal> colIdx(aIn->getNodeNumEntries(), 0);
    Teuchos::ArrayRCP<Scalar> values(aIn->getNodeNumEntries(), 0.0);
    RCP<const OMapType> outRowMap = aIn->getRowMap()->clone<Node>(oNode);
    RCP<const OMapType> outColMap = aIn->getColMap()->clone<Node>(oNode);
    KCrsMatrix aLocal = aIn->getLocalMatrix();
    Ordinal nRows = aLocal.numRows();
    Ordinal nNZ = aLocal.nnz();
    ScalarArray aVal = aLocal.values;
    OrdinalArray aColIdx = aLocal.graph.entries;
    KCrsMatrix::row_map_type aRowMap = aLocal.graph.row_map;
    auto aRowMapO = OrdinalArray("aRowMapO", nRows + 1);;
    auto aColIdxO = OrdinalArray("aColIdxO", nNZ);
    auto aValO = ScalarArray("aValO", nNZ);
    RCP<const OMapType> clonedRangeMap = aIn->getRangeMap()->clone<Node> (oNode);
    RCP<const OMapType> clonedDomainMap = aIn->getDomainMap()->clone<Node> (oNode);

    reorderAndCopy(aRowMap, aColIdx, aVal, perm, permInv, nRows, aRowMapO, aColIdxO, aValO);

    for (int i = 0; i <= nRows; i++)
    {
        rowPtr[i] = aRowMapO[i];
    }
    for(int i = 0; i < nRows; i++) 
    {
        for (int k = rowPtr[i]; k < rowPtr[i+1]; k++) {
            colIdx[k] = aColIdxO[k];
            values[k] = aValO[k];
    //        std::cout << "Debug: (" << i << "," << colIdx[k] << ") = " << values[k] << std::endl;
        }
    }

    aOut = rcp(new TCrsMatrix(outRowMap, outColMap, rowPtr, colIdx, values, Teuchos::null));  
    aOut->fillComplete(clonedDomainMap, clonedRangeMap);

}
#endif

int main(int argc, char **argv)
{
    //Initialize MPI 
    Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);
    Kokkos::initialize(argc, argv);
    MyExecSpace::print_configuration(std::cout);


    if (argc != 12) {
        print_usage(argv);
        return -1;
    }

    auto filename = std::string(argv[1]);
    int nFact = atoi(argv[2]);
    int nTrisol = atoi(argv[3]);
    int icFlag = atoi(argv[4]);
    int level = atoi(argv[5]);
    int exFlag = atoi(argv[6]);
    int solverFlag = atoi(argv[7]);
    double omega = atof(argv[8]);
    double shift = atof(argv[9]);
    int guessFlag = atoi(argv[10]);
    int blkSz = atoi(argv[11]);

    //Extract the matrix name from the file path
    //Note this is not portable only works for Unix paths
    unsigned int pos1 = filename.rfind('.');   
    unsigned int pos2 = filename.rfind('/');   
    std::string matrixName = filename.substr(pos2 + 1, pos1 - pos2 - 1);
    
    RCP<Node> node = Tpetra::DefaultPlatform::getDefaultPlatform().getNode();
    
    //Get a communicator object (required only to read the matrix, for now)
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

    //Highly obfuscated line to read the matrix
    RCP<TCrsMatrix> a(Tpetra::MatrixMarket::Reader<TCrsMatrix>::readSparseFile(filename, comm, node, 
                true, true, true));
    RCP<TCrsMatrix> aOrd;

    //Zoltan2 for ordering
#if __GPU_VERS__
#if __REORDER__
    Teuchos::ParameterList zoltanParams;
    zoltanParams.set("order_method", "rcm");
    ZMatrixAdapter aAdapted(a); 
    Zoltan2::OrderingProblem<ZMatrixAdapter> orderingProblem(&aAdapted, &zoltanParams);
    orderingProblem.solve();
    std::cout << "Zoltan2 ordering solution computed" << std::endl;


    Zoltan2::OrderingSolution<Ordinal, Ordinal> *ordSoln = orderingProblem.getSolution();
    Ordinal* perm = ordSoln->getPermutation();
    ordSoln->computeInverse();
    Ordinal* invPerm = ordSoln->getPermutation(true);

    //This code prints out the wrong thing.
    createReorderedTpetraMatrix(a, aOrd, perm, invPerm, node);
#else 
    aOrd = a;
#endif 
#else
    aOrd = a;
#endif

   // Fic<Scalar, Ordinal, Ordinal, Node, MyExecSpace> iterIC(aOrd, nFact, nTrisol, level);
    Fldl<Scalar, Ordinal, Ordinal, Node, MyExecSpace> iterLDL(aOrd, nFact, nTrisol, level, omega, 
            shift, guessFlag, blkSz);
    Filu<Scalar, Ordinal, Ordinal, Node, MyExecSpace> iterILU(aOrd, nFact, nTrisol, level, omega, 
            shift, guessFlag, blkSz);
    Fic<Scalar, Ordinal, Ordinal, Node, MyExecSpace> iterIC(aOrd, nFact, nTrisol, level, omega, 
            shift, guessFlag, blkSz);

    if (icFlag == 1)
    {
        if (exFlag == 0)
        {
            std::cout << "Inexact: LDL("<< level <<")" << std::endl;
            iterLDL.initialize();
            iterLDL.compute();
            iterLDL.printStatus();
        }
    }
    else if (icFlag == 0) 
    {
        if (exFlag == 0)
        {
            std::cout<< "Inexact: ILU(" << level << ")" << std::endl;
            iterILU.initialize();
            iterILU.compute();
            iterILU.printStatus();
        }
    }
    else if (icFlag == 2)
    {
        if (exFlag == 0)
        {
            std::cout << "Inexact: IC("<< level <<")" << std::endl;
            iterIC.initialize();
            iterIC.compute();
            iterIC.printStatus();
        }
    }

   

    //Exact ILU(k)
#if __GPU_VERS__
    RCP<const TCrsMatrix> constA = aOrd;
    Teuchos::ParameterList rilukList;
    rilukList.set("fact: iluk level-of-fill", level);
    Ifpack2::RILUK<TCrsMatrix> exILU(constA);
    if (exFlag == 1)
    {
        std::cout << "Exact: RILU(k)" << std::endl;
        exILU.setParameters(rilukList);
        exILU.initialize();
        exILU.compute();
    }
#endif

    //Initialize a Belos Linear Problem
    bool verbose = false;
    int frequency = 500;
    int numRhs = 1; 
    int blockSize = 1;
    int numBlocks = 50;
    int maxIters = 2500;
    int maxRestarts = 600;
    MagnitudeType tol = 1e-6;
    RCP<const Tpetra::Map<Ordinal>> map = aOrd->getDomainMap();
    RCP<const Tpetra::Map<Ordinal>> colMap = aOrd->getRangeMap();
    RCP<MultiVector<Scalar, Ordinal>> B, X;
    X = rcp(new MultiVector<Scalar, int>(map, numRhs));
    B = rcp(new MultiVector<Scalar, int>(map, numRhs));
    //BMultiVecTraits::MvRandom(*B);
    BMultiVecTraits::MvInit(*B, 1.0);
    BMultiVecTraits::MvInit(*X, 0.0);
    Teuchos::ParameterList belosList;
    belosList.set( "Block Size", blockSize );              // Blocksize to be used by iterative solver
    belosList.set( "Num Blocks", numBlocks );
    belosList.set( "Maximum Iterations", maxIters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Output Frequency", frequency);
    belosList.set( "Maximum Restarts", maxRestarts);
    int verbLevel = Belos::Errors + Belos::Warnings;
    verbLevel += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
    belosList.set( "Verbosity", verbLevel);

    Belos::LinearProblem<Scalar, TMultiVector, OperatorType> problem(aOrd, X, B);

    bool set = problem.setProblem();
    //Set right preconditioner for the problem
    if (exFlag == 0)
    {
        if (icFlag == 0)
        {
            problem.setRightPrec(rcp(&iterILU, false));
        }
        else if (icFlag == 1)
        {
            problem.setRightPrec(rcp(&iterLDL, false));
        }
        else if (icFlag == 2)
        {
            problem.setRightPrec(rcp(&iterIC, false));
        }
       
    }
    else 
    {
#if __GPU_VERS__
        problem.setRightPrec(rcp(&exILU, false));
#endif
    }
    Belos::ReturnType ret;
    int solverIters = 0;
    if (solverFlag == 0)
    {
        Belos::BlockCGSolMgr<Scalar, TMultiVector, OperatorType> solver(rcp(&problem, false), rcp(&belosList, false));
        ret = solver.solve();
        solverIters = solver.getNumIters();
    }
    else if (solverFlag == 1)
    {
        Belos::BlockGmresSolMgr<Scalar, TMultiVector, OperatorType> solver(rcp(&problem, false), 
                rcp(&belosList, false));
        ret = solver.solve();
        solverIters = solver.getNumIters();
    }

    if (ret == Belos::Converged)
    {
        std::cout << "(" << matrixName << "," << nFact << "," << nTrisol << "," << level << "," << omega << "," << shift << ")" << " :" << " iters= " << solverIters << std::endl;
    }
    else 
    {
        std::cout << "(" << matrixName << "," << nFact << "," << nTrisol << "," << level << "," << omega << "," << shift << ")" << " :" << " iters= *" << std::endl;

    }

    //End
    Kokkos::finalize();
    return 0;
}
