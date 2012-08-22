/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER
*/

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
#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
#  include <Kokkos_OpenMPNode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
#  include <Kokkos_TPINode.hpp>
#endif

/* the unit tests require some explicit instantiations that is not enabled in the build of the library
   specifically,

   usually weird stuff, like double wrappers around int matrices
  
   we currently have no mixed scalar support on CUDA GPUs (neither Cusp nor CUSPARSE provide it)
 */

namespace Tpetra {

  TPETRA_CRSMATRIX_CONVERT_INSTANT(double,int,int,int,Kokkos::SerialNode)
  // int matrix support for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_DOUBLE) || defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
  TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOSCLASSIC_TBB)
    TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_OPENMP)
    TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::OpenMPNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
    TPETRA_CRSMATRIX_INSTANT(int,int,int,Kokkos::TPINode)
# endif
#endif

  // mixed: double wrapper for int matrix for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_DOUBLE) 
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(double,int,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOSCLASSIC_OPENMP)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(double,int,int,int,Kokkos::OpenMPNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,int,int,int,Kokkos::OpenMPNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_TBB)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(double,int,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(double,int,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(double,int,int,int,Kokkos::TPINode)
# endif
#endif

  // mixed: float wrapper for int matrix for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_FLOAT) 
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOSCLASSIC_OPENMP)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::OpenMPNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::OpenMPNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_TBB)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(float,int,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(float,int,int,int,Kokkos::TPINode)
# endif
#endif

  // mixed: complex<float> wrapper for int matrix for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<float>,int,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<float>,int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOSCLASSIC_OPENMP)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<float>,int,int,int,Kokkos::OpenMPNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<float>,int,int,int,Kokkos::OpenMPNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_TBB)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<float>,int,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<float>,int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<float>,int,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<float>,int,int,int,Kokkos::TPINode)
# endif
#endif

  // mixed: complex<double> wrapper for int matrix for CrsMatrix unit test
#if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<double>,int,int,int,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<double>,int,int,int,Kokkos::SerialNode)
# if defined(HAVE_KOKKOSCLASSIC_OPENMP)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<double>,int,int,int,Kokkos::OpenMPNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<double>,int,int,int,Kokkos::OpenMPNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_TBB)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<double>,int,int,int,Kokkos::TBBNode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<double>,int,int,int,Kokkos::TBBNode)
# endif
# if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
    TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(std::complex<double>,int,int,int,Kokkos::TPINode)
    TPETRA_CRSMATRIX_SOLVEOP_INSTANT(std::complex<double>,int,int,int,Kokkos::TPINode)
# endif
#endif

typedef Kokkos::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_TPETRA_INST_INT_LONG
template Teuchos::RCP< const Map<int,long,Kokkos::DefaultNode::DefaultNodeType> >
createContigMap<int,long>(size_t numElements, size_t numLocalElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

TPETRA_CRSMATRIX_INSTANT(double,int,long,Node)
TPETRA_MAP_INSTANT(int,long,Node)
TPETRA_BLOCKMAP_INSTANT(int,long,Node)
TPETRA_DIRECTORY_INSTANT(int,long,Node)

namespace Ext {
TPETRAEXT_BLOCKEXTRACTION_INSTANT(double,int,long,Node)
}
#endif

}

#endif
