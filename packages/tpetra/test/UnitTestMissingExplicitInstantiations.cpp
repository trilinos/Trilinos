#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_CrsMatrixSolveOp.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

// definitions
#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_CrsMatrixMultiplyOp_def.hpp"
#include "Tpetra_CrsMatrixSolveOp_def.hpp"
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


/* the unit tests require some explicit instantiations that may not be enabled in the build of the library
   specifically,

    CrsMatrix<float,int,int,N>

    CrsMatrixMultiplyOp<float,float,int,int,N> 
    CrsMatrixSolveOp   <float,float,int,int,N>

    CrsMatrixMultiplyOp<float,int  ,int,int,N> 
    CrsMatrixSolveOp   <float,int  ,int,int,N> 

 */

namespace Tpetra {

#if !defined(HAVE_TPETRA_INST_FLOAT)
  TPETRA_CRSMATRIX_INSTANT(float,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,float,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,float,int,int,Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
    TPETRA_CRSMATRIX_INSTANT(float,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,float,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,float,int,int,Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_CRSMATRIX_INSTANT(float,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,float,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,float,int,int,Kokkos::TPINode)
#endif
#if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_FLOAT)
    TPETRA_CRSMATRIX_INSTANT(float,int,int,Kokkos::ThrustGPUNode)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,float,int,int,Kokkos::ThrustGPUNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,float,int,int,Kokkos::ThrustGPUNode)
#endif
#endif

  // int matrix
  TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
    TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::TPINode)
#endif
#if defined(HAVE_KOKKOS_THRUST)
    TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::ThrustGPUNode)
#endif

  // mixed: float wrapper for int matrix

  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::TPINode)
#endif
#if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_FLOAT)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::ThrustGPUNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::ThrustGPUNode)
#endif


}

#endif
