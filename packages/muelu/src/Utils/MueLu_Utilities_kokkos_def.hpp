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

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_EPETRA
# ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
# endif
#endif

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
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <MueLu_Utilities_kokkos_decl.hpp>

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
  using Xpetra::EpetraCrsMatrix;   // TODO: mv in Xpetra_UseShortNamesScalar
  using Xpetra::EpetraMultiVector;
#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar> Utils_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixDiagonal(const Matrix& A) {
    // FIXME Kokkos

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
  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Utils_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol) {
    // FIXME Kokkos
    RCP<const Map> rowMap = A.getRowMap();
    RCP<Vector> diag      = VectorFactory::Build(rowMap);
    ArrayRCP<SC> diagVals = diag->getDataNonConst(0);

    size_t numRows = rowMap->getNodeNumElements();

    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i = 0; i < numRows; ++i) {
      A.getLocalRowView(i, cols, vals);

      LO j = 0;
      for (; j < cols.size(); ++j) {
        if (Teuchos::as<size_t>(cols[j]) == i) {
          if(Teuchos::ScalarTraits<SC>::magnitude(vals[j]) > tol)
            diagVals[i] = Teuchos::ScalarTraits<SC>::one() / vals[j];
          else
            diagVals[i]=Teuchos::ScalarTraits<SC>::zero();
          break;
        }
      }
      if (j == cols.size()) {
        // Diagonal entry is absent
        diagVals[i]=Teuchos::ScalarTraits<SC>::zero();
      }
    }
    diagVals=null;

    return diag;
  } //GetMatrixDiagonalInverse

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar> Utils_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetLumpedMatrixDiagonal(const Matrix &A) {
    // FIXME: Kokkos
    size_t numRows = A.getRowMap()->getNodeNumElements();
    Teuchos::ArrayRCP<SC> diag(numRows);

    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i = 0; i < numRows; ++i) {
      A.getLocalRowView(i, cols, vals);

      diag[i] = Teuchos::ScalarTraits<Scalar>::zero();
      for (LO j = 0; j < cols.size(); ++j) {
        diag[i] += Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]);
      }
    }

    return diag;
  } //GetMatrixDiagonal

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Utils_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixOverlappedDiagonal(const Matrix& A) {
    // FIXME: Kokkos
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
  void Utils_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doInverse) {
    // FIXME: Kokkos
#ifdef HAVE_MUELU_TPETRA
    try {
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpOp = Op2NonConstTpetraCrs(Op);

      Tpetra::Vector<SC,LO,GO,NO> x(tpOp.getRowMap(), scalingVector());
      if (doInverse){
        Tpetra::Vector<SC,LO,GO,NO> xi(tpOp.getRowMap());
        xi.reciprocal(x);
        tpOp.leftScale(xi);

      } else {
        tpOp.leftScale(x);
      }
    } catch(...) {
      throw Exceptions::RuntimeError("Matrix scaling has not been implemented Epetra");
    }
#else
    throw Exceptions::RuntimeError("Matrix scaling has not been implemented Epetra");
#endif // HAVE_MUELU_TPETRA
  } //ScaleMatrix()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const SC>& scalingVector, bool doInverse,
                               bool doFillComplete,
                               bool doOptimizeStorage)
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
        // FIXME?
        // Utils2_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
        throw std::runtime_error("FIXME");
        break;

      default:
        throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled.");
        break;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
#ifdef HAVE_MUELU_TPETRA
    try {
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpOp = Op2NonConstTpetraCrs(Op);

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
  ArrayRCP<const bool> Utils_kokkos<SC, LO, GO, NO>::DetectDirichletRows(const Matrix& A, const typename Teuchos::ScalarTraits<SC>::magnitudeType& tol) {
    LO numRows = A.getNodeNumRows();

    typedef Teuchos::ScalarTraits<SC> STS;

    ArrayRCP<bool> boundaryNodes(numRows, true);
    for (LO row = 0; row < numRows; row++) {
      ArrayView<const LO> indices;
      ArrayView<const SC> vals;
      A.getLocalRowView(row, indices, vals);

      size_t nnz = A.getNumEntriesInLocalRow(row);
      if (nnz > 1)
        for (size_t col = 0; col < nnz; col++)
          if ( (indices[col] != row) && STS::magnitude(vals[col]) > tol) {
            boundaryNodes[row] = false;
            break;
          }
    }

    return boundaryNodes;
  }

  template <class SC, class LO, class GO, class NO>
  void Utils_kokkos<SC, LO, GO, NO>::findDirichletRows(Teuchos::RCP<Matrix> A,
                                                     std::vector<LO>& dirichletRows) {
    dirichletRows.resize(0);
    for(size_t i=0; i<A->getNodeNumRows(); i++) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> values;
      A->getLocalRowView(i,indices,values);
      int nnz=0;
      for (int j=0; j<indices.size(); j++) {
        if (abs(values[j]) > 1.0e-16) {
          nnz++;
        }
      }
      if (nnz == 1 || nnz == 2) {
        dirichletRows.push_back(i);
      }
    }
  }

  template<class SC, class LO, class GO, class NO>
  void Utils_kokkos<SC, LO, GO, NO>::findDirichletCols(Teuchos::RCP<Matrix> A,
                                                     std::vector<LO>& dirichletRows,
                                                     std::vector<LO>& dirichletCols) {
    Teuchos::RCP<const Map> domMap = A->getDomainMap();
    Teuchos::RCP<const Map> colMap = A->getColMap();
    Teuchos::RCP< Xpetra::Export<LO,GO,NO> > exporter
      = Xpetra::ExportFactory<LO,GO,NO>::Build(colMap,domMap);
    Teuchos::RCP<MultiVector> myColsToZero = MultiVectorFactory::Build(colMap,1);
    Teuchos::RCP<MultiVector> globalColsToZero = MultiVectorFactory::Build(domMap,1);
    myColsToZero->putScalar((SC)0.0);
    globalColsToZero->putScalar((SC)0.0);
    for(size_t i=0; i<dirichletRows.size(); i++) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> values;
      A->getLocalRowView(dirichletRows[i],indices,values);
      for(int j=0; j<indices.size(); j++)
        myColsToZero->replaceLocalValue(indices[j],0,(SC)1.0);
    }
    globalColsToZero->doExport(*myColsToZero,*exporter,Xpetra::ADD);
    myColsToZero->doImport(*globalColsToZero,*exporter,Xpetra::INSERT);
    Teuchos::ArrayRCP<const SC> myCols = myColsToZero->getData(0);
    dirichletCols.resize(colMap->getNodeNumElements());
    for(size_t i=0; i<colMap->getNodeNumElements(); i++) {
      if(abs(myCols[i])>0.0)
        dirichletCols[i]=1;
      else
        dirichletCols[i]=0;
    }
  }

  template<class SC, class LO, class GO, class NO>
  void Utils_kokkos<SC, LO, GO, NO>::Apply_BCsToMatrixRows(Teuchos::RCP<Matrix>& A,
                                                         std::vector<LO>& dirichletRows) {
    for(size_t i=0; i<dirichletRows.size(); i++) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> values;
      A->getLocalRowView(dirichletRows[i],indices,values);
      std::vector<SC> vec;
      vec.resize(indices.size());
      Teuchos::ArrayView<SC> zerovalues(vec);
      for(int j=0; j<indices.size(); j++)
        zerovalues[j]=(SC)1.0e-32;
      A->replaceLocalValues(dirichletRows[i],indices,zerovalues);
    }
  }

  template<class SC, class LO, class GO, class NO>
  void Utils_kokkos<SC, LO, GO, NO>::Apply_BCsToMatrixCols(Teuchos::RCP<Matrix>& A,
                                                         std::vector<LO>& dirichletCols) {
    for(size_t i=0; i<A->getNodeNumRows(); i++) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> values;
      A->getLocalRowView(i,indices,values);
      std::vector<SC> vec;
      vec.resize(indices.size());
      Teuchos::ArrayView<SC> zerovalues(vec);
      for(int j=0; j<indices.size(); j++) {
        if(dirichletCols[indices[j]]==1)
          zerovalues[j]=(SC)1.0e-32;
        else
          zerovalues[j]=values[j];
      }
      A->replaceLocalValues(i,indices,zerovalues);
    }
  }

  template<class SC, class LO, class GO, class NO>
  void Utils_kokkos<SC, LO, GO, NO>::Remove_Zeroed_Rows(Teuchos::RCP<Matrix>& A, double tol) {
    Teuchos::RCP<const Map> rowMap = A->getRowMap();
    RCP<Matrix> DiagMatrix = MatrixFactory::Build(rowMap,1);
    RCP<Matrix> NewMatrix = MatrixFactory::Build(rowMap,1);
    for(size_t i=0; i<A->getNodeNumRows(); i++) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> values;
      A->getLocalRowView(i,indices,values);
      int nnz=0;
      for (int j=0; j<indices.size(); j++) {
        if (abs(values[j]) > tol) {
          nnz++;
        }
      }
      SC one = (SC)1.0;
      SC zero = (SC)0.0;
      GO row = rowMap->getGlobalElement(i);
      if (nnz == 0) {
        DiagMatrix->insertGlobalValues(row,
                                       Teuchos::ArrayView<GO>(&row,1),
                                       Teuchos::ArrayView<SC>(&one,1));
      }
      else {
        DiagMatrix->insertGlobalValues(row,
                                       Teuchos::ArrayView<GO>(&row,1),
                                       Teuchos::ArrayView<SC>(&zero,1));
      }
    }
    DiagMatrix->fillComplete();
    A->fillComplete();
    // add matrices together
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Utils2_kokkos<SC,LO,GO,NO>::TwoMatrixAdd(*DiagMatrix,false,(SC)1.0,*A,false,(SC)1.0,NewMatrix,*out);
    NewMatrix->fillComplete();
    A=NewMatrix;
  }

} //namespace MueLu

#define MUELU_UTILITIES_KOKKOS_SHORT
#endif // MUELU_UTILITIES_KOKKOS_DEF_HPP
