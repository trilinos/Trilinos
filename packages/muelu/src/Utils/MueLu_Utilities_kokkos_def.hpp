// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_UTILITIES_KOKKOS_DEF_HPP
#define MUELU_UTILITIES_KOKKOS_DEF_HPP

#include <algorithm>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_EPETRA
# ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
# endif
#endif

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_getDiagCopy.hpp>

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_MultiVectorIn.h>
#include <EpetraExt_BlockMapIn.h>
#include <Xpetra_EpetraUtils.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <EpetraExt_BlockMapOut.h>
#endif

#ifdef HAVE_MUELU_TPETRA
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#endif

#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraMap.hpp>
#endif

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <MueLu_Utilities_kokkos_decl.hpp>
#include <MueLu_Utilities.hpp>

#include <KokkosKernels_Handle.hpp>
#include <KokkosGraph_RCM.hpp>


namespace MueLu {


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  GetMatrixDiagonalInverse(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
                           typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol, const bool doLumped) {
    Teuchos::TimeMonitor MM = *Teuchos::TimeMonitor::getNewTimer("Utilities_kokkos::GetMatrixDiagonalInverse");
    // Some useful type definitions
    using Matrix            = Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using Map               = Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>;
    using Vector            = Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using VectorFactory     = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using local_matrix_type = typename Matrix::local_matrix_type;
    // using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
    using value_type        = typename local_matrix_type::value_type;
    using ordinal_type      = typename local_matrix_type::ordinal_type;
    using execution_space   = typename local_matrix_type::execution_space;
    // using memory_space      = typename local_matrix_type::memory_space;
    // Be careful with this one, if using Kokkos::ArithTraits<Scalar>
    // you are likely to run into errors when handling std::complex<>
    // a good way to work around that is to use the following:
    // using KAT = Kokkos::ArithTraits<Kokkos::ArithTraits<Scalar>::val_type> >
    // here we have: value_type = Kokkos::ArithTraits<Scalar>::val_type
    using KAT               = Kokkos::ArithTraits<value_type>;

    // Get/Create distributed objects
    RCP<const Map> rowMap = A.getRowMap();
    RCP<Vector> diag      = VectorFactory::Build(rowMap,false);

    // Now generate local objects
    local_matrix_type localMatrix = A.getLocalMatrixDevice();
    auto diagVals = diag->getDeviceLocalView(Xpetra::Access::ReadWrite);

    ordinal_type numRows = localMatrix.graph.numRows();

    // Note: 2019-11-21, LBV
    // This could be implemented with a TeamPolicy over the rows
    // and a TeamVectorRange over the entries in a row if performance
    // becomes more important here.
    if (!doLumped)
      Kokkos::parallel_for("Utilities_kokkos::GetMatrixDiagonalInverse",
                           Kokkos::RangePolicy<ordinal_type, execution_space>(0, numRows),
                           KOKKOS_LAMBDA(const ordinal_type rowIdx) {
                             bool foundDiagEntry = false;
                             auto myRow = localMatrix.rowConst(rowIdx);
                             for(ordinal_type entryIdx = 0; entryIdx < myRow.length; ++entryIdx) {
                               if(myRow.colidx(entryIdx) == rowIdx) {
                                 foundDiagEntry = true;
                                 if(KAT::magnitude(myRow.value(entryIdx)) > KAT::magnitude(tol)) {
                                   diagVals(rowIdx, 0) = KAT::one() / myRow.value(entryIdx);
                                 } else {
                                   diagVals(rowIdx, 0) = KAT::zero();
                                 }
                                 break;
                               }
                             }

                             if(!foundDiagEntry) {diagVals(rowIdx, 0) = KAT::zero();}
                           });
    else
      Kokkos::parallel_for("Utilities_kokkos::GetMatrixDiagonalInverse",
                           Kokkos::RangePolicy<ordinal_type, execution_space>(0, numRows),
                           KOKKOS_LAMBDA(const ordinal_type rowIdx) {
                             auto myRow = localMatrix.rowConst(rowIdx);
                             for(ordinal_type entryIdx = 0; entryIdx < myRow.length; ++entryIdx) {
                               diagVals(rowIdx, 0) += KAT::magnitude(myRow.value(entryIdx));
                             }
                             if(KAT::magnitude(diagVals(rowIdx, 0)) > KAT::magnitude(tol))
                               diagVals(rowIdx, 0) = KAT::one() / diagVals(rowIdx, 0);
                             else
                               diagVals(rowIdx, 0) = KAT::zero();

                           });

    return diag;
  } //GetMatrixDiagonalInverse

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal,Node> >
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal,Node>::
  GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol, const bool doLumped)
  {
    return MueLu::GetMatrixDiagonalInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A,tol,doLumped);
  }

  template <class Node>
  Teuchos::RCP<Xpetra::Vector<double,int,int,Node> >
  Utilities_kokkos<double,int,int,Node>::
  GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol, const bool doLumped)
  {
    return MueLu::GetMatrixDiagonalInverse<double, int, int, Node>(A,tol,doLumped);
  }


} //namespace MueLu

#define MUELU_UTILITIES_KOKKOS_SHORT
#endif // MUELU_UTILITIES_KOKKOS_DEF_HPP
