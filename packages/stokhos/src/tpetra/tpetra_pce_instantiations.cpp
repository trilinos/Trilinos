// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Stokhos_Sacado.hpp"
#include "kokkos_pce_specializations.hpp"

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

#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_MultiVector_def.hpp"
#include "Tpetra_Vector_def.hpp"
#include "Tpetra_CrsGraph_def.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_CrsMatrixMultiplyOp_def.hpp"

#include "TpetraExt_MatrixMatrix_def.hpp"
#include "Tpetra_RowMatrixTransposer_def.hpp"

namespace Tpetra {

#if defined(HAVE_KOKKOSCLASSIC_SERIAL)
  TPETRA_MULTIVECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::DoNotUse::SerialNode)
  TPETRA_VECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::DoNotUse::SerialNode)
  TPETRA_CRSMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::DoNotUse::SerialNode)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(Scalar,Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::DoNotUse::SerialNode)
  TPETRA_MATRIXMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::DoNotUse::SerialNode)
  TPETRA_ROWMATRIXTRANSPOSE_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::DoNotUse::SerialNode)
#endif // defined(HAVE_KOKKOSCLASSIC_SERIAL)

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  TPETRA_MULTIVECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TBBNode)
  TPETRA_VECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TBBNode)
  TPETRA_CRSMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TBBNode)
TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(Scalar,Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TBBNode)
  TPETRA_MATRIXMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TBBNode)
  TPETRA_ROWMATRIXTRANSPOSE_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TBBNode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  TPETRA_MULTIVECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode)
  TPETRA_VECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode)
  TPETRA_CRSMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(Scalar,Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode)
  TPETRA_MATRIXMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode)
  TPETRA_ROWMATRIXTRANSPOSE_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
  TPETRA_MULTIVECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::OpenMPNode)
  TPETRA_VECTOR_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::OpenMPNode)
  TPETRA_CRSMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::OpenMPNode)
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(Scalar,Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::OpenMPNode)
  TPETRA_MATRIXMATRIX_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::OpenMPNode)
  TPETRA_ROWMATRIXTRANSPOSE_INSTANT(Scalar,LocalOrdinal,GlobalOrdinal,KokkosClassic::OpenMPNode)
#endif

}

#endif
