#include <examples/Kokkos_DummySparseKernelClass.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_CrsMatrixSolveOp.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

// definitions
#include "Tpetra_MultiVector_def.hpp"
#include "Tpetra_Vector_def.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_CrsGraph_def.hpp"
#include "Tpetra_CrsMatrixMultiplyOp_def.hpp"
#include "Tpetra_CrsMatrixSolveOp_def.hpp"

typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType DefNode;
typedef KokkosExamples::DummySparseKernel<DefNode>             SparseOp;

namespace Tpetra {

  // need float multivector
#if !defined(HAVE_TPETRA_INST_FLOAT)
  TPETRA_MULTIVECTOR_INSTANT(float,int,int,DefNode)
#endif

  // float matrix support for dummy sparse kernel class
  template class CrsMatrix < float , int , int , DefNode , SparseOp >;
}

#endif
