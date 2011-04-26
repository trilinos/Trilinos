#include <Teuchos_UnitTestHarness.hpp>

// include all the Tpetra headers twice to check if headers are protected by #ifndef FILE_HPP #define ... #endif
#include <Tpetra_BlockMap.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrixSolveOp.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_BlockCrsGraph.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Directory.hpp>
#include <Tpetra_CrsMatrixMultiplyOp.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_VbrMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <TpetraExt_BlockExtraction.hpp>
#include <Tpetra_MMHelpers.hpp>
#include <Tpetra_MatrixMatrix.hpp>
#include <Tpetra_MatrixIO.hpp>

#include <Tpetra_BlockMap.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrixSolveOp.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_BlockCrsGraph.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Directory.hpp>
#include <Tpetra_CrsMatrixMultiplyOp.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_VbrMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <TpetraExt_BlockExtraction.hpp>
#include <Tpetra_MMHelpers.hpp>
#include <Tpetra_MatrixMatrix.hpp>
#include <Tpetra_MatrixIO.hpp>


// this test do not compile if Tpetra objects are declared outside of the namespace Tpetra.
typedef int BinaryFunctorAdapter;
typedef int BinaryFunctorAdapterWithAlphaBeta;
typedef int BinaryOp;
typedef int BlockCrsGraph;
typedef int BlockMap;
typedef int BlockMultiVector;
typedef int CrsGraph;
typedef int CrsMatrix;
typedef int CrsMatrixMultiplyOp;
typedef int CrsMatrixSolveOp;
typedef int DefaultPlatform;
typedef int Directory;
typedef int DistObject;
typedef int Distributor;
typedef int Export;
typedef int FDStencil;
typedef int HybridPlatform;
typedef int Import;
typedef int KernelOp;
typedef int Map;
typedef int MpiPlatform;
typedef int mprec_mult;
typedef int MultiVector;
typedef int OneOp;
typedef int Operator;
typedef int ReductionGlob;
typedef int RowGraph;
typedef int RowMatrix;
typedef int RowMatrixTransposer;
typedef int RTIPreTransformReductionAdapter;
typedef int RTIReductionAdapter;
typedef int ScaleKernel;
typedef int SerialPlatform;
typedef int StdOpKernel;
typedef int TeuchosValueTypeReductionOpAdapter;
typedef int TransformReductionGlob;
typedef int UnaryFunctorAdapter;
typedef int VbrMatrix;
typedef int Vector;
typedef int ZeroOp;

TEUCHOS_UNIT_TEST( Compilation, Bug_ClassDeclarationOutsideOfNamespace )
{
   
}

