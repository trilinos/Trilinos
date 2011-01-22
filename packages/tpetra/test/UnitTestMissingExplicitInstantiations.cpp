#include "Tpetra_Map.hpp"
#include "Tpetra_BlockMap.hpp"
#include "Tpetra_Directory.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_CrsMatrixSolveOp.hpp"
#include "TpetraExt_BlockExtraction.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

// definitions
#include "Tpetra_Map_def.hpp"
#include "Tpetra_BlockMap_def.hpp"
#include "Tpetra_Directory_def.hpp"
#include "Tpetra_MultiVector_def.hpp"
#include "Tpetra_Vector_def.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_CrsGraph_def.hpp"
#include "Tpetra_CrsMatrixMultiplyOp_def.hpp"
#include "Tpetra_CrsMatrixSolveOp_def.hpp"
#include "TpetraExt_BlockExtraction_def.hpp"
// nodes
#include <Kokkos_SerialNode.hpp>
#if defined(HAVE_KOKKOS_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
#  include <Kokkos_TPINode.hpp>
#endif
#if defined(HAVE_KOKKOS_THRUST)
#  include <Kokkos_ThrustGPUNode.hpp>
#endif

/* the unit tests require some explicit instantiations that is not enabled in the build of the library
   specifically,

   usually weird stuff, like double wrappers around int matrices
 */

namespace Tpetra {

  // int matrix support for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_DOUBLE) || defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
  TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOS_TBB)
    TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::TPINode)
# endif
# if defined(HAVE_KOKKOS_THRUST)
    TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::ThrustGPUNode)
# endif
#endif

  // mixed: double wrapper for int matrix for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_DOUBLE) 
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(double,int,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOS_TBB)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(double,int,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(double,int,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,int,int,int,Kokkos::TPINode)
# endif
# if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_DOUBLE)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(double,int,int,int,Kokkos::ThrustGPUNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,int,int,int,Kokkos::ThrustGPUNode)
# endif
#endif

  // mixed: float wrapper for int matrix for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_FLOAT) 
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOS_TBB)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::TPINode)
# endif
# if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_FLOAT)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::ThrustGPUNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::ThrustGPUNode)
# endif
#endif

  // mixed: complex<float> wrapper for int matrix for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<float>,int,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<float>,int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOS_TBB)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<float>,int,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<float>,int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<float>,int,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<float>,int,int,int,Kokkos::TPINode)
# endif
# if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_COMPLEX_FLOAT)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<float>,int,int,int,Kokkos::ThrustGPUNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<float>,int,int,int,Kokkos::ThrustGPUNode)
# endif
#endif

  // mixed: complex<double> wrapper for int matrix for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<double>,int,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<double>,int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOS_TBB)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<double>,int,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<double>,int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<double>,int,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<double>,int,int,int,Kokkos::TPINode)
# endif
# if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_COMPLEX_DOUBLE)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<double>,int,int,int,Kokkos::ThrustGPUNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<double>,int,int,int,Kokkos::ThrustGPUNode)
# endif
#endif

typedef Kokkos::DefaultNode::DefaultNodeType Node;

template Teuchos::RCP< const Map<int,long,Kokkos::DefaultNode::DefaultNodeType> >
createContigMap<int,long>(size_t numElements, size_t numLocalElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

TPETRA_CRSMATRIX_INSTANT(double,int,long,Node)
TPETRA_MAP_INSTANT(int,long,Node)
TPETRA_BLOCKMAP_INSTANT(int,long,Node)
TPETRA_DIRECTORY_INSTANT(int,long,Node)

namespace Ext {
TPETRAEXT_BLOCKEXTRACTION_INSTANT(double,int,long,Node)
}

}

#endif
