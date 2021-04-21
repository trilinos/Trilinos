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

#include <KokkosKernels_Handle.hpp>
#include <KokkosGraph_RCM.hpp>


namespace MueLu {


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar> Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GetMatrixDiagonal(const Matrix& A) {
    // FIXME_KOKKOS

    size_t numRows = A.getRowMap()->getNodeNumElements();
    Teuchos::ArrayRCP<SC> diag(numRows);

    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i = 0; i < numRows; ++i) {
      A.getLocalRowView(i, cols, vals);

      LO j = 0;
      for (; j < cols.size(); ++j) {
        if (Teuchos::as<size_t>(cols[j]) == i) {
          diag[i] = vals[j];
          break;
        }
      }
      if (j == cols.size()) {
        // Diagonal entry is absent
        diag[i] = Teuchos::ScalarTraits<SC>::zero();
      }
    }

    return diag;
  } //GetMatrixDiagonal

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol, const bool doLumped) {
    Teuchos::TimeMonitor MM = *Teuchos::TimeMonitor::getNewTimer("Utilities_kokkos::GetMatrixDiagonalInverse");
    // Some useful type definitions
    using local_matrix_type = typename Matrix::local_matrix_type;
    using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
    using value_type        = typename local_matrix_type::value_type;
    using ordinal_type      = typename local_matrix_type::ordinal_type;
    using execution_space   = typename local_matrix_type::execution_space;
    using memory_space      = typename local_matrix_type::memory_space;
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
    local_matrix_type localMatrix = A.getLocalMatrix();
    auto diagVals = diag->getDeviceLocalView();

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
  RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GetMatrixOverlappedDiagonal(const Matrix& A) {
    // FIXME_KOKKOS
    RCP<const Map> rowMap = A.getRowMap(), colMap = A.getColMap();
    RCP<Vector>    localDiag     = VectorFactory::Build(rowMap);

    try {
       const CrsMatrixWrap* crsOp = dynamic_cast<const CrsMatrixWrap*>(&A);
       if (crsOp == NULL) {
         throw Exceptions::RuntimeError("cast to CrsMatrixWrap failed");
       }
       Teuchos::ArrayRCP<size_t> offsets;
       crsOp->getLocalDiagOffsets(offsets);
       crsOp->getLocalDiagCopy(*localDiag,offsets());
    }
    catch (...) {
      ArrayRCP<SC>   localDiagVals = localDiag->getDataNonConst(0);
      Teuchos::ArrayRCP<SC> diagVals = GetMatrixDiagonal(A);
      for (LO i = 0; i < localDiagVals.size(); i++)
        localDiagVals[i] = diagVals[i];
      localDiagVals = diagVals = null;
    }

    RCP<Vector> diagonal = VectorFactory::Build(colMap);
    RCP< const Import> importer;
    importer = A.getCrsGraph()->getImporter();
    if (importer == Teuchos::null) {
      importer = ImportFactory::Build(rowMap, colMap);
    }
    diagonal->doImport(*localDiag, *(importer), Xpetra::INSERT);

    return diagonal;
  } //GetMatrixOverlappedDiagonal

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const SC>& scalingVector,
                   bool doInverse, bool doFillComplete, bool doOptimizeStorage)
  {
    SC one = Teuchos::ScalarTraits<SC>::one();
    Teuchos::ArrayRCP<SC> sv(scalingVector.size());
    if (doInverse) {
      for (int i = 0; i < scalingVector.size(); ++i)
        sv[i] = one / scalingVector[i];
    } else {
      for (int i = 0; i < scalingVector.size(); ++i)
        sv[i] = scalingVector[i];
    }

    switch (Op.getRowMap()->lib()) {
      case Xpetra::UseTpetra:
        MyOldScaleMatrix_Tpetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;

      case Xpetra::UseEpetra:
        // FIXME_KOKKOS
        // Utils2_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
        throw std::runtime_error("FIXME");
#ifndef __NVCC__ //prevent nvcc warning
        break;
#endif

      default:
        throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled.");
#ifndef __NVCC__ //prevent nvcc warning
        break;
#endif
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Epetra(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& /* Op */, const Teuchos::ArrayRCP<Scalar>& /* scalingVector */, bool /* doFillComplete */, bool /* doOptimizeStorage */) {
    throw Exceptions::RuntimeError("MyOldScaleMatrix_Epetra: Epetra needs SC=double and LO=GO=int.");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
#ifdef HAVE_MUELU_TPETRA
    try {
      Tpetra::CrsMatrix<SC,LO,GO,NO>& tpOp = Op2NonConstTpetraCrs(Op);

      const RCP<const Tpetra::Map<LO,GO,NO> > rowMap    = tpOp.getRowMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > domainMap = tpOp.getDomainMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > rangeMap  = tpOp.getRangeMap();

      size_t maxRowSize = tpOp.getNodeMaxNumRowEntries();
      if (maxRowSize == Teuchos::as<size_t>(-1)) // hasn't been determined yet
        maxRowSize = 20;

      std::vector<SC> scaledVals(maxRowSize);
      if (tpOp.isFillComplete())
        tpOp.resumeFill();

      if (Op.isLocallyIndexed() == true) {
        Teuchos::ArrayView<const LO> cols;
        Teuchos::ArrayView<const SC> vals;

        for (size_t i = 0; i < rowMap->getNodeNumElements(); ++i) {
          tpOp.getLocalRowView(i, cols, vals);
          size_t nnz = tpOp.getNumEntriesInLocalRow(i);
          if (nnz > maxRowSize) {
            maxRowSize = nnz;
            scaledVals.resize(maxRowSize);
          }
          for (size_t j = 0; j < nnz; ++j)
            scaledVals[j] = vals[j]*scalingVector[i];

          if (nnz > 0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0], nnz);
            tpOp.replaceLocalValues(i, cols, valview);
          }
        } //for (size_t i=0; ...

      } else {
        Teuchos::ArrayView<const GO> cols;
        Teuchos::ArrayView<const SC> vals;

        for (size_t i = 0; i < rowMap->getNodeNumElements(); ++i) {
          GO gid = rowMap->getGlobalElement(i);
          tpOp.getGlobalRowView(gid, cols, vals);
          size_t nnz = tpOp.getNumEntriesInGlobalRow(gid);
          if (nnz > maxRowSize) {
            maxRowSize = nnz;
            scaledVals.resize(maxRowSize);
          }
          // FIXME FIXME FIXME FIXME FIXME FIXME
          for (size_t j = 0; j < nnz; ++j)
            scaledVals[j] = vals[j]*scalingVector[i]; //FIXME i or gid?

          if (nnz > 0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0], nnz);
            tpOp.replaceGlobalValues(gid, cols, valview);
          }
        } //for (size_t i=0; ...
      }

      if (doFillComplete) {
        if (domainMap == Teuchos::null || rangeMap == Teuchos::null)
          throw Exceptions::RuntimeError("In Utils_kokkos::Scaling: cannot fillComplete because the domain and/or range map hasn't been defined");

        RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
        params->set("Optimize Storage",    doOptimizeStorage);
        params->set("No Nonlocal Changes", true);
        Op.fillComplete(Op.getDomainMap(), Op.getRangeMap(), params);
      }
    } catch(...) {
      throw Exceptions::RuntimeError("Only Tpetra::CrsMatrix types can be scaled (Err.1)");
    }
#else
    throw Exceptions::RuntimeError("Matrix scaling is not possible because Tpetra has not been enabled.");
#endif
  } //MyOldScaleMatrix_Tpetra()


  template <class SC, class LO, class GO, class NO>
  Kokkos::View<bool*, typename NO::device_type>
  DetectDirichletRows(const Xpetra::Matrix<SC,LO,GO,NO>& A,
                      const typename Teuchos::ScalarTraits<SC>::magnitudeType& tol,
                      const bool count_twos_as_dirichlet) {
    using impl_scalar_type = typename Kokkos::ArithTraits<SC>::val_type;
    using ATS        = Kokkos::ArithTraits<impl_scalar_type>;
    using range_type = Kokkos::RangePolicy<LO, typename NO::execution_space>;

    auto localMatrix = A.getLocalMatrix();
    LO   numRows     = A.getNodeNumRows();

    Kokkos::View<bool*, typename NO::device_type> boundaryNodes(Kokkos::ViewAllocateWithoutInitializing("boundaryNodes"), numRows);
    if (count_twos_as_dirichlet)
      Kokkos::parallel_for("MueLu:Utils::DetectDirichletRows_Twos_As_Dirichlet", range_type(0,numRows),
                           KOKKOS_LAMBDA(const LO row) {
                             auto rowView = localMatrix.row(row);
                             auto length  = rowView.length;

                             boundaryNodes(row) = true;
                             if (length > 2) {
                               decltype(length) colID;
                               for (colID = 0; colID < length; colID++)
                                 if ((rowView.colidx(colID) != row) &&
                                     (ATS::magnitude(rowView.value(colID)) > tol)) {
                                   if (!boundaryNodes(row))
                                     break;
                                   boundaryNodes(row) = false;
                                 }
                               if (colID == length)
                                 boundaryNodes(row) = true;
                             }
                           });
    else
      Kokkos::parallel_for("MueLu:Utils::DetectDirichletRows", range_type(0,numRows),
                           KOKKOS_LAMBDA(const LO row) {
                             auto rowView = localMatrix.row(row);
                             auto length  = rowView.length;

                             boundaryNodes(row) = true;
                             for (decltype(length) colID = 0; colID < length; colID++)
                               if ((rowView.colidx(colID) != row) &&
                                   (ATS::magnitude(rowView.value(colID)) > tol)) {
                                 boundaryNodes(row) = false;
                                 break;
                               }
                           });

    return boundaryNodes;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Kokkos::View<bool*, typename Node::device_type>
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  DetectDirichletRows(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& tol, const bool count_twos_as_dirichlet) {
    return MueLu::DetectDirichletRows<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A, tol, count_twos_as_dirichlet);
  }

  template <class Node>
  Kokkos::View<bool*, typename Node::device_type>
  Utilities_kokkos<double,int,int,Node>::
  DetectDirichletRows(const Xpetra::Matrix<double,int,int,Node>& A, const typename Teuchos::ScalarTraits<double>::magnitudeType& tol, const bool count_twos_as_dirichlet) {
    return MueLu::DetectDirichletRows<double,int,int,Node>(A, tol,count_twos_as_dirichlet);
  }


  template <class SC, class LO, class GO, class NO>
  Kokkos::View<bool*, typename NO::device_type>
  DetectDirichletCols(const Xpetra::Matrix<SC,LO,GO,NO>& A,
                      const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows) {
    using ATS        = Kokkos::ArithTraits<SC>;
    using impl_ATS = Kokkos::ArithTraits<typename ATS::val_type>;
    using range_type = Kokkos::RangePolicy<LO, typename NO::execution_space>;

    SC zero = ATS::zero();
    SC one = ATS::one();

    auto localMatrix = A.getLocalMatrix();
    LO   numRows     = A.getNodeNumRows();

    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > domMap = A.getDomainMap();
    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > colMap = A.getColMap();
    Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > myColsToZero = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(colMap,1);
    myColsToZero->putScalar(zero);
    auto myColsToZeroView = myColsToZero->getDeviceLocalView();
    // Find all local column indices that are in Dirichlet rows, record in myColsToZero as 1.0
    Kokkos::parallel_for("MueLu:Utils::DetectDirichletCols1", range_type(0,numRows),
                         KOKKOS_LAMBDA(const LO row) {
                           if (dirichletRows(row)) {
                             auto rowView = localMatrix.row(row);
                             auto length  = rowView.length;

                             for (decltype(length) colID = 0; colID < length; colID++)
                               myColsToZeroView(rowView.colidx(colID),0) = one;
                           }
                         });

    Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > globalColsToZero = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(domMap,1);
    globalColsToZero->putScalar(zero);
    Teuchos::RCP<Xpetra::Export<LO,GO,NO> > exporter = Xpetra::ExportFactory<LO,GO,NO>::Build(colMap,domMap);
    // export to domain map
    globalColsToZero->doExport(*myColsToZero,*exporter,Xpetra::ADD);
    // import to column map
    myColsToZero->doImport(*globalColsToZero,*exporter,Xpetra::INSERT);

    auto myCols = myColsToZero->getDeviceLocalView();
    size_t numColEntries = colMap->getNodeNumElements();
    Kokkos::View<bool*, typename NO::device_type> dirichletCols(Kokkos::ViewAllocateWithoutInitializing("dirichletCols"), numColEntries);
    const typename ATS::magnitudeType eps = 2.0*ATS::eps();

    Kokkos::parallel_for("MueLu:Utils::DetectDirichletCols2", range_type(0,numColEntries),
                         KOKKOS_LAMBDA (const size_t i) {
                           dirichletCols(i) = impl_ATS::magnitude(myCols(i,0))>eps;
                         });
    return dirichletCols;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Kokkos::View<bool*, typename Node::device_type>
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  DetectDirichletCols(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
                      const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows) {
    return MueLu::DetectDirichletCols<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A, dirichletRows);
  }

  template <class Node>
  Kokkos::View<bool*, typename Node::device_type>
  Utilities_kokkos<double,int,int,Node>::
  DetectDirichletCols(const Xpetra::Matrix<double,int,int,Node>& A,
                      const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows) {
    return MueLu::DetectDirichletCols<double,int,int,Node>(A, dirichletRows);
  }


  // Zeros out rows
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  ZeroDirichletRows(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
                    const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows,
                    Scalar replaceWith) {
    using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

    auto localMatrix = A->getLocalMatrix();
    LocalOrdinal numRows = A->getNodeNumRows();

    Kokkos::parallel_for("MueLu:Utils::ZeroDirichletRows", range_type(0,numRows),
                         KOKKOS_LAMBDA(const LocalOrdinal row) {
                           if (dirichletRows(row)) {
                             auto rowView = localMatrix.row(row);
                             auto length  = rowView.length;
                             for (decltype(length) colID = 0; colID < length; colID++)
                               rowView.value(colID) = replaceWith;
                           }
                         });
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ZeroDirichletRows(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                    const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows,
                    Scalar replaceWith) {
    MueLu::ZeroDirichletRows<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A, dirichletRows, replaceWith);
  }

  template <class Node>
  void
  Utilities_kokkos<double,int,int,Node>::
  ZeroDirichletRows(RCP<Xpetra::Matrix<double, int, int, Node> >& A,
                    const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows,
                    double replaceWith) {
    return MueLu::ZeroDirichletRows<double,int,int,Node>(A, dirichletRows, replaceWith);
  }


  // Zeros out rows
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  ZeroDirichletRows(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& X,
                    const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows,
                    Scalar replaceWith) {
    using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;
    auto myCols = X->getDeviceLocalView();
    size_t numVecs = X->getNumVectors();
    Kokkos::parallel_for("MueLu:Utils::ZeroDirichletRows_MV", range_type(0,dirichletRows.size()),
                         KOKKOS_LAMBDA(const size_t i) {
                           if (dirichletRows(i)) {
                             for(size_t j=0; j<numVecs; j++)
                               myCols(i,j) = replaceWith;
                           }
                         });
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ZeroDirichletRows(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& X,
                    const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows,
                    Scalar replaceWith) {
    MueLu::ZeroDirichletRows<Scalar,LocalOrdinal,GlobalOrdinal,Node>(X, dirichletRows, replaceWith);
  }

  template <class Node>
  void
  Utilities_kokkos<double,int,int,Node>::
  ZeroDirichletRows(RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, int, int, Node> >& X,
                    const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows,
                    double replaceWith) {
    return MueLu::ZeroDirichletRows<double,int,int,Node>(X, dirichletRows, replaceWith);
  }


  // Zeros out columns
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  ZeroDirichletCols(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                    const Kokkos::View<const bool*, typename Node::device_type>& dirichletCols,
                    Scalar replaceWith) {
    using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

    auto localMatrix = A->getLocalMatrix();
    LocalOrdinal numRows = A->getNodeNumRows();

    Kokkos::parallel_for("MueLu:Utils::ZeroDirichletCols", range_type(0,numRows),
                         KOKKOS_LAMBDA(const LocalOrdinal row) {
                           auto rowView = localMatrix.row(row);
                           auto length  = rowView.length;
                           for (decltype(length) colID = 0; colID < length; colID++)
                             if (dirichletCols(rowView.colidx(colID))) {
                               rowView.value(colID) = replaceWith;
                             }
                         });
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ZeroDirichletCols(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                    const Kokkos::View<const bool*, typename NO::device_type>& dirichletCols,
                    Scalar replaceWith) {
    MueLu::ZeroDirichletCols<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A, dirichletCols, replaceWith);
  }

  template <class Node>
  void
  Utilities_kokkos<double,int,int,Node>::
  ZeroDirichletCols(RCP<Xpetra::Matrix<double,int,int,Node> >& A,
                    const Kokkos::View<const bool*, typename Node::device_type>& dirichletCols,
                    double replaceWith) {
    return MueLu::ZeroDirichletCols<double,int,int,Node>(A, dirichletCols, replaceWith);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  RealValuedToScalarMultiVector(RCP<RealValuedMultiVector > X) {
    RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Xscalar;
#if defined(HAVE_XPETRA_TPETRA) && (defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) || defined(HAVE_TPETRA_INST_COMPLEX_FLOAT))
    using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

    // Need to cast the real-valued multivector to Scalar=complex
    if ((typeid(Scalar).name() == typeid(std::complex<double>).name()) ||
        (typeid(Scalar).name() == typeid(std::complex<float>).name())) {
      size_t numVecs = X->getNumVectors();
      Xscalar = Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(X->getMap(),numVecs);
      auto XVec = X->getDeviceLocalView();
      auto XVecScalar = Xscalar->getDeviceLocalView();

      Kokkos::parallel_for("MueLu:Utils::RealValuedToScalarMultiVector", range_type(0,X->getLocalLength()),
                           KOKKOS_LAMBDA(const size_t i) {
                             for (size_t j=0; j<numVecs; j++)
                               XVecScalar(i,j) = XVec(i,j);
                           });
    } else
#endif
      Xscalar = rcp_dynamic_cast<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(X);
    return Xscalar;
  }

  template <class Node>
  RCP<Xpetra::MultiVector<double,int,int,Node> >
  Utilities_kokkos<double,int,int,Node>::
  RealValuedToScalarMultiVector(RCP<Xpetra::MultiVector<Magnitude,int,int,Node> > X) {
    return X;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > ReverseCuthillMcKee(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Op) {
    using local_matrix_type = typename Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type;
    using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
    using lno_nnz_view_t    = typename local_graph_type::entries_type::non_const_type;
    using device            = typename local_graph_type::device_type;
    using execution_space   = typename local_matrix_type::execution_space;
    using ordinal_type      = typename local_matrix_type::ordinal_type;

    local_graph_type localGraph = Op.getLocalMatrix().graph;

    lno_nnz_view_t rcmOrder = KokkosGraph::Experimental::graph_rcm
      <device, typename local_graph_type::row_map_type, typename local_graph_type::entries_type, lno_nnz_view_t>
      (localGraph.row_map, localGraph.entries);

    RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > retval = 
      Xpetra::VectorFactory<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>::Build(Op.getRowMap());

    // Copy out and reorder data
    auto view1D = Kokkos::subview(retval->getDeviceLocalView(),Kokkos::ALL (), 0);
    Kokkos::parallel_for("Utilities_kokkos::ReverseCuthillMcKee",
                         Kokkos::RangePolicy<ordinal_type, execution_space>(0, localGraph.numRows()),
                         KOKKOS_LAMBDA(const ordinal_type rowIdx) {
                           view1D(rcmOrder(rowIdx)) = rowIdx;
                         });
    return retval;
  }
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > CuthillMcKee(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Op) {
    using local_matrix_type = typename Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type;
    using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
    using lno_nnz_view_t    = typename local_graph_type::entries_type::non_const_type;
    using device            = typename local_graph_type::device_type;
    using execution_space   = typename local_matrix_type::execution_space;
    using ordinal_type      = typename local_matrix_type::ordinal_type;

    local_graph_type localGraph = Op.getLocalMatrix().graph;
    LocalOrdinal numRows = localGraph.numRows();

    lno_nnz_view_t rcmOrder = KokkosGraph::Experimental::graph_rcm
      <device, typename local_graph_type::row_map_type, typename local_graph_type::entries_type, lno_nnz_view_t>
      (localGraph.row_map, localGraph.entries);

    RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > retval = 
      Xpetra::VectorFactory<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>::Build(Op.getRowMap());

    // Copy out data
    auto view1D = Kokkos::subview(retval->getDeviceLocalView(),Kokkos::ALL (), 0);
    // Since KokkosKernels produced RCM, also reverse the order of the view to get CM
    Kokkos::parallel_for("Utilities_kokkos::ReverseCuthillMcKee",
                         Kokkos::RangePolicy<ordinal_type, execution_space>(0, numRows),
                         KOKKOS_LAMBDA(const ordinal_type rowIdx) {
                           view1D(rcmOrder(numRows - 1 - rowIdx)) = rowIdx;
                         });
    return retval;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> >
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReverseCuthillMcKee(const Matrix &Op) {
    return MueLu::ReverseCuthillMcKee<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Op);
  }

  template <class Node>
  Teuchos::RCP<Xpetra::Vector<int,int,int,Node> >  
  Utilities_kokkos<double,int,int,Node>::ReverseCuthillMcKee(const Matrix &Op) {
    return MueLu::ReverseCuthillMcKee<double,int,int,Node>(Op);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> >
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CuthillMcKee(const Matrix &Op) {
    return MueLu::CuthillMcKee<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Op);
  }

  template <class Node>
  Teuchos::RCP<Xpetra::Vector<int,int,int,Node> >  
  Utilities_kokkos<double,int,int,Node>::CuthillMcKee(const Matrix &Op) {
    return MueLu::CuthillMcKee<double,int,int,Node>(Op);
  }

  // Applies Ones-and-Zeros to matrix rows
  // Takes a Boolean array.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                       const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows) {
    TEUCHOS_ASSERT(A->isFillComplete());
    using ATS        = Kokkos::ArithTraits<Scalar>;
    using impl_ATS = Kokkos::ArithTraits<typename ATS::val_type>;
    using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

    RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > domMap = A->getDomainMap();
    RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > ranMap = A->getRangeMap();
    RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > Rmap = A->getRowMap();
    RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > Cmap = A->getColMap();

    TEUCHOS_ASSERT(static_cast<size_t>(dirichletRows.size()) == Rmap->getNodeNumElements());

    const Scalar one  = impl_ATS::one();
    const Scalar zero = impl_ATS::zero();

    auto localMatrix = A->getLocalMatrix();
    auto localRmap = Rmap->getLocalMap();
    auto localCmap = Cmap->getLocalMap();

    Kokkos::parallel_for("MueLu::Utils::ApplyOAZ",range_type(0,dirichletRows.extent(0)),
                         KOKKOS_LAMBDA(const LocalOrdinal row) {
                           if (dirichletRows(row)){
                             auto rowView = localMatrix.row(row);
                             auto length = rowView.length;
                             auto row_gid = localRmap.getGlobalElement(row);
                             auto row_lid = localCmap.getLocalElement(row_gid);

                             for (decltype(length) colID = 0; colID < length; colID++)
                               if (rowView.colidx(colID) == row_lid)
                                 rowView.value(colID) = one;
                               else
                                 rowView.value(colID) = zero;
                           }
                         });
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Utilities_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                       const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows) {
    MueLu::ApplyOAZToMatrixRows<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A, dirichletRows);
  }

  template <class Node>
  void
  Utilities_kokkos<double,int,int,Node>::
  ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<double,int,int,Node> >& A,
                       const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows) {
    MueLu::ApplyOAZToMatrixRows<double,int,int,Node>(A, dirichletRows);
  }

} //namespace MueLu

#define MUELU_UTILITIES_KOKKOS_SHORT
#endif // MUELU_UTILITIES_KOKKOS_DEF_HPP
