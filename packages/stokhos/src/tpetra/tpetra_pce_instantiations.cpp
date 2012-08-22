// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_Sacado.hpp"

typedef Stokhos::StandardStorage<int,double> Storage;
typedef Sacado::PCE::OrthogPoly<double,Storage> pce_type;
typedef pce_type Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"

#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include <Kokkos_SerialNode.hpp>
#if defined(HAVE_KOKKOSCLASSIC_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
#  include <Kokkos_TPINode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
#  include <Kokkos_OpenMPNode.hpp>
#endif
// #if defined(HAVE_KOKKOSCLASSIC_THRUST)
// #  include <Kokkos_ThrustGPUNode.hpp>
// #endif

#include "Tpetra_MultiVector_def.hpp"
#include "Tpetra_Vector_def.hpp"
#include "Tpetra_CrsGraph_def.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_CrsMatrixMultiplyOp_def.hpp"

#include "TpetraExt_MatrixMatrix_def.hpp"
#include "Tpetra_RowMatrixTransposer_def.hpp"

namespace Tpetra {

  TPETRA_MULTIVECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::SerialNode)
  TPETRA_VECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::SerialNode)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(Scalar,Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::SerialNode)
  TPETRA_MATRIXMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::SerialNode)
  TPETRA_ROWMATRIXTRANSPOSE_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::SerialNode)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_MULTIVECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TBBNode)
  TPETRA_VECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TBBNode)
  TPETRA_CRSMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TBBNode)
TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(Scalar,Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TBBNode)
  TPETRA_MATRIXMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TBBNode)
  TPETRA_ROWMATRIXTRANSPOSE_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TBBNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_MULTIVECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TPINode)
  TPETRA_VECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TPINode)
  TPETRA_CRSMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TPINode)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(Scalar,Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TPINode)
  TPETRA_MATRIXMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TPINode)
  TPETRA_ROWMATRIXTRANSPOSE_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_MULTIVECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::OpenMPNode)
  TPETRA_VECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::OpenMPNode)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(Scalar,Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::OpenMPNode)
  TPETRA_MATRIXMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::OpenMPNode)
  TPETRA_ROWMATRIXTRANSPOSE_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::OpenMPNode)
#endif

  // Not sure if we support thrust yet
// #if defined(HAVE_KOKKOSCLASSIC_THRUST) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
//   TPETRA_MULTIVECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode)
//   TPETRA_VECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode)
//   TPETRA_CRSMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode)
  // TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(Scalar,Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode)
  // TPETRA_MATRIXMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode)
  // TPETRA_ROWMATRIXTRANSPOSE_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::ThrustGPUNode)
// #endif
  
}

#endif
