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
#ifndef MUELU_UTILITIES_DEF_HPP
#define MUELU_UTILITIES_DEF_HPP

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
//#include <Xpetra_DefaultPlatform.hpp>
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

#include <Xpetra_MatrixMatrix.hpp>

#include <MueLu_Utilities_decl.hpp>
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_ML)
#include <ml_operator.h>
#include <ml_epetra_utils.h>
#endif

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
  using Xpetra::EpetraCrsMatrix;   // TODO: mv in Xpetra_UseShortNamesScalar
  using Xpetra::EpetraMultiVector;
#endif

#ifdef HAVE_MUELU_EPETRA
  //defined after Utils class
  template<typename SC,typename LO,typename GO,typename NO>
  RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB);
#endif

#ifdef HAVE_MUELU_EPETRA
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Epetra_MultiVector> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2EpetraMV(const RCP<MultiVector> Vec) {
    RCP<const EpetraMultiVector > tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
    return tmpVec->getEpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Epetra_MultiVector> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstEpetraMV(RCP<MultiVector> Vec) {
    RCP<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
    return tmpVec->getEpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Epetra_MultiVector& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstEpetraMV(MultiVector &Vec) {
    const EpetraMultiVector& tmpVec = dynamic_cast<const EpetraMultiVector&>(Vec);
    return *(tmpVec.getEpetra_MultiVector());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Epetra_MultiVector& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2EpetraMV(const MultiVector& Vec) {
    const EpetraMultiVector& tmpVec = dynamic_cast<const EpetraMultiVector&>(Vec);
    return *(tmpVec.getEpetra_MultiVector());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2EpetraCrs(RCP<const Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const EpetraCrsMatrix>& tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrix();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(RCP<Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Epetra_CrsMatrix& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2EpetraCrs(const Matrix& Op) {
    try {
      const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
      try {
        const EpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<const EpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getEpetra_CrsMatrix();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Epetra_CrsMatrix& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(Matrix& Op) {
    try {
      CrsMatrixWrap& crsOp = dynamic_cast<CrsMatrixWrap&>(Op);
      try {
        EpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<EpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getEpetra_CrsMatrixNonConst();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Epetra_Map& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map2EpetraMap(const Map& map) {
    RCP<const Xpetra::EpetraMap> xeMap = rcp_dynamic_cast<const Xpetra::EpetraMap>(rcpFromRef(map));
    if (xeMap == Teuchos::null)
      throw Exceptions::BadCast("Utils::Map2EpetraMap : Cast from Xpetra::Map to Xpetra::EpetraMap failed");
    return xeMap->getEpetra_Map();
  }
#endif

#ifdef HAVE_MUELU_TPETRA
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2TpetraMV(RCP<MultiVector> const Vec) {
    RCP<const TpetraMultiVector > tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
    return tmpVec->getTpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstTpetraMV(RCP<MultiVector> Vec) {
    RCP<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
    return tmpVec->getTpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstTpetraMV(MultiVector& Vec) {
    const TpetraMultiVector& tmpVec = dynamic_cast<const TpetraMultiVector&>(Vec);
    return *(tmpVec.getTpetra_MultiVector());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstTpetraMV2(MultiVector &Vec) {
    const TpetraMultiVector& tmpVec = dynamic_cast<const TpetraMultiVector&>(Vec);
    return tmpVec.getTpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2TpetraMV(const MultiVector& Vec) {
    const TpetraMultiVector& tmpVec = dynamic_cast<const TpetraMultiVector&>(Vec);
    return *(tmpVec.getTpetra_MultiVector());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op) {
    // Get the underlying Tpetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_CrsMatrix();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(RCP<Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(const Matrix& Op) {
    try {
      const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
      try {
        const TpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<const TpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_CrsMatrix();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(Matrix& Op) {
    try {
      CrsMatrixWrap& crsOp = dynamic_cast<CrsMatrixWrap&>(Op);
      try {
        TpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<TpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_CrsMatrixNonConst();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraRow(RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> crsMat = crsOp->getCrsMatrix();
    const RCP<const TpetraCrsMatrix> tmp_Crs = rcp_dynamic_cast<const TpetraCrsMatrix>(crsMat);
    RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmp_BlockCrs;
    if(!tmp_Crs.is_null()) {
      return tmp_Crs->getTpetra_CrsMatrixNonConst();
    }
    else {
      tmp_BlockCrs= rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(crsMat);
      if (tmp_BlockCrs.is_null())
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix and Xpetra::TpetraBlockCrsMatrix failed");
      return tmp_BlockCrs->getTpetra_BlockCrsMatrixNonConst();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraRow(RCP<Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> crsMat = crsOp->getCrsMatrix();
    const RCP<const TpetraCrsMatrix> tmp_Crs = rcp_dynamic_cast<const TpetraCrsMatrix>(crsMat);
    RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmp_BlockCrs;
    if(!tmp_Crs.is_null()) {
      return tmp_Crs->getTpetra_CrsMatrixNonConst();
    }
    else {
      tmp_BlockCrs= rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(crsMat);
      if (tmp_BlockCrs.is_null())
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix and Xpetra::TpetraBlockCrsMatrix failed");
      return tmp_BlockCrs->getTpetra_BlockCrsMatrixNonConst();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal,Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map2TpetraMap(const Map& map) {
    const RCP<const TpetraMap>& tmp_TMap = rcp_dynamic_cast<const TpetraMap>(rcpFromRef(map));
    if (tmp_TMap == Teuchos::null)
      throw Exceptions::BadCast("Utils::Map2TpetraMap : Cast from Xpetra::Map to Xpetra::TpetraMap failed");
    return tmp_TMap->getTpetra_Map();
  }
#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op) {
    if (Op.is_null())
      return Teuchos::null;

    return rcp(new CrsMatrixWrap(Op));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixDiagonal(const Matrix& A) {

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
  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixDiagonalInverse(const Matrix& A,Magnitude tol) {
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
  Teuchos::ArrayRCP<Scalar> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetLumpedMatrixDiagonal(const Matrix &A) {
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
  RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixOverlappedDiagonal(const Matrix& A) {
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
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doInverse) {
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
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ResidualNorm(const Operator& Op, const MultiVector& X, const MultiVector& RHS) {
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
    const size_t numVecs = X.getNumVectors();

    RCP<MultiVector> RES = Residual(Op, X, RHS);
    Teuchos::Array<Magnitude> norms(numVecs);
    RES->norm2(norms);

    return norms;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Residual(const Operator& Op, const MultiVector& X, const MultiVector& RHS) {
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
    const size_t numVecs = X.getNumVectors();

    SC one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();

    RCP<MultiVector> RES = MultiVectorFactory::Build(Op.getRangeMap(), numVecs, false); // no need to initialize to zero
    Op.apply(X, *RES, Teuchos::NO_TRANS, one, zero);
    RES->update(one, RHS, negone);

    return RES;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(const std::string& fileName, const Matrix& Op) {
    Write("rowmap_"    + fileName, *(Op.getRowMap()));
    Write("colmap_"    + fileName, *(Op.getColMap()));
    Write("domainmap_" + fileName, *(Op.getDomainMap()));
    Write("rangemap_"  + fileName, *(Op.getRangeMap()));

    const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp.getCrsMatrix();
#if defined(HAVE_MUELU_EPETRA)
    const RCP<const EpetraCrsMatrix>& tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx != Teuchos::null) {
#if defined(HAVE_MUELU_EPETRAEXT)
      RCP<const Epetra_CrsMatrix> A = tmp_ECrsMtx->getEpetra_CrsMatrix();
      int rv = EpetraExt::RowMatrixToMatrixMarketFile(fileName.c_str(), *A);
      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::RowMatrixToMatrixMarketFile return value of " + toString(rv));
#else
      throw Exceptions::RuntimeError("Compiled without EpetraExt");
#endif
      return;
    }
#endif

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraCrsMatrix>& tmp_TCrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_TCrsMtx != Teuchos::null) {
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A = tmp_TCrsMtx->getTpetra_CrsMatrix();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeSparseFile(fileName, A);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw Exceptions::BadCast("Could not cast to EpetraCrsMatrix or TpetraCrsMatrix in matrix writing");
  } //Write

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(const std::string& fileName, const RCP<const Map>& map){
    Xpetra::UnderlyingLib lib = map->lib();

    if (lib == Xpetra::UseEpetra) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, ::Xpetra::Exceptions::BadCast, "Epetra can only be used with Scalar=double and Ordinal=int");

    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>                          reader_type;
      typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>                            map_type;
      typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>            multivector_type;

      RCP<const map_type>   temp = toTpetra(map);
      RCP<multivector_type> TMV  = reader_type::readDenseFile(fileName,map->getComm(),map->getNode(),temp);
      RCP<MultiVector>      rmv  = Xpetra::toXpetra(TMV);
      return rmv;
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }

    return Teuchos::null;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >
  Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMap(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm) {
    if (lib == Xpetra::UseEpetra) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, ::Xpetra::Exceptions::BadCast, "Epetra can only be used with Scalar=double and Ordinal=int");
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>                          reader_type;

      RCP<Node> node = rcp(new Node());

      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > tMap = reader_type::readMapFile(fileName, comm, node);
      if (tMap.is_null())
        throw Exceptions::RuntimeError("The Tpetra::Map returned from readSparseFile() is null.");

      return Xpetra::toXpetra(tMap);
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }

    return Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm, bool binary) {
    if (binary == false) {
      // Matrix Market file format (ASCII)
      if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
        Epetra_CrsMatrix *eA;
        const RCP<const Epetra_Comm> epcomm = Xpetra::toEpetra(comm);
        int rv = EpetraExt::MatrixMarketFileToCrsMatrix(fileName.c_str(), *epcomm, eA);
        if (rv != 0)
          throw Exceptions::RuntimeError("EpetraExt::MatrixMarketFileToCrsMatrix return value of " + toString(rv));

        RCP<Epetra_CrsMatrix> tmpA = rcp(eA);
        RCP<Matrix>           A    = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA);

        return A;
#else
        throw Exceptions::RuntimeError("MueLu has not been compiled with Epetra and EpetraExt support.");
#endif
      } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
        typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;

        typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;

        //RCP<Node> node = Xpetra::DefaultPlatform::getDefaultPlatform().getNode();
        Teuchos::ParameterList pl = Teuchos::ParameterList();
        RCP<Node> node = rcp(new Node(pl));
        bool callFillComplete = true;

        RCP<sparse_matrix_type> tA = reader_type::readSparseFile(fileName, comm, node, callFillComplete);

        if (tA.is_null())
          throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from readSparseFile() is null.");

        RCP<TpetraCrsMatrix> tmpA1 = rcp(new TpetraCrsMatrix(tA));
        RCP<CrsMatrix>       tmpA2 = rcp_implicit_cast<CrsMatrix>(tmpA1);
        RCP<Matrix>          A     = rcp(new CrsMatrixWrap(tmpA2));

        return A;
#else
        throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
      } else {
        throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
      }
    } else {
      // Custom file format (binary)
      std::ifstream ifs(fileName.c_str(), std::ios::binary);
      TEUCHOS_TEST_FOR_EXCEPTION(!ifs.good(), Exceptions::RuntimeError, "Can not read \"" << fileName << "\"");
      int m, n, nnz;
      ifs.read(reinterpret_cast<char*>(&m),   sizeof(m));
      ifs.read(reinterpret_cast<char*>(&n),   sizeof(n));
      ifs.read(reinterpret_cast<char*>(&nnz), sizeof(nnz));

      int myRank = comm->getRank();

      GO indexBase = 0;
      RCP<Map>    rowMap = MapFactory::Build(lib, m, (myRank == 0 ? m : 0), indexBase, comm), rangeMap  = rowMap;
      RCP<Map>    colMap = MapFactory::Build(lib, n, (myRank == 0 ? n : 0), indexBase, comm), domainMap = colMap;
      RCP<Matrix> A   = MatrixFactory::Build(rowMap, colMap, 1);

      TEUCHOS_TEST_FOR_EXCEPTION(sizeof(int) != sizeof(GO), Exceptions::RuntimeError, "Incompatible sizes");

      if (myRank == 0) {
        Teuchos::Array<GO> inds;
        Teuchos::Array<SC> vals;
        for (int i = 0; i < m; i++) {
          int row, rownnz;
          ifs.read(reinterpret_cast<char*>(&row),    sizeof(row));
          ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));
          inds.resize(rownnz);
          vals.resize(rownnz);
          for (int j = 0; j < rownnz; j++) {
            int index;
            ifs.read(reinterpret_cast<char*>(&index), sizeof(index));
            inds[j] = Teuchos::as<GO>(index);
          }
          for (int j = 0; j < rownnz; j++) {
            double value;
            ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
            vals[j] = Teuchos::as<SC>(value);
          }
          A->insertGlobalValues(row, inds, vals);
        }
      }

      A->fillComplete(domainMap, rangeMap);

      return A;
    }

    return Teuchos::null;

  } //Read()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(const std::string&   fileName,
                                                                      const RCP<const Map> rowMap,
                                                                            RCP<const Map> colMap,
                                                                      const RCP<const Map> domainMap,
                                                                      const RCP<const Map> rangeMap,
                                                                      const bool           callFillComplete,
                                                                      const bool           binary,
                                                                      const bool           tolerant,
                                                                      const bool           debug
                                                                     ) {
    TEUCHOS_TEST_FOR_EXCEPTION(rowMap.is_null(), Exceptions::RuntimeError, "Utils::Read() : rowMap cannot be null");

    RCP<const Map> domain = (domainMap.is_null() ? rowMap : domainMap);
    RCP<const Map> range  = (rangeMap .is_null() ? rowMap : rangeMap);

    const Xpetra::UnderlyingLib lib = rowMap->lib();
    if (binary == false) {
      if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
        Epetra_CrsMatrix *eA;
        const RCP<const Epetra_Comm> epcomm = Xpetra::toEpetra(rowMap->getComm());
        const Epetra_Map& epetraRowMap    = Map2EpetraMap(*rowMap);
        const Epetra_Map& epetraDomainMap = (domainMap.is_null() ? epetraRowMap : Map2EpetraMap(*domainMap));
        const Epetra_Map& epetraRangeMap  = (rangeMap .is_null() ? epetraRowMap : Map2EpetraMap(*rangeMap));
        int rv;
        if (colMap.is_null()) {
          rv = EpetraExt::MatrixMarketFileToCrsMatrix(fileName.c_str(), epetraRowMap, epetraRangeMap, epetraDomainMap, eA);

        } else {
          const Epetra_Map& epetraColMap  = Map2EpetraMap(*colMap);
          rv = EpetraExt::MatrixMarketFileToCrsMatrix(fileName.c_str(), epetraRowMap, epetraColMap, epetraRangeMap, epetraDomainMap, eA);
        }

        if (rv != 0)
          throw Exceptions::RuntimeError("EpetraExt::MatrixMarketFileToCrsMatrix return value of " + toString(rv));

        RCP<Epetra_CrsMatrix> tmpA = rcp(eA);
        RCP<Matrix>           A    = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA);

        return A;
#else
        throw Exceptions::RuntimeError("MueLu has not been compiled with Epetra and EpetraExt support.");
#endif
      } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
        typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
        typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>             reader_type;
        typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>               map_type;

        const RCP<const map_type> tpetraRowMap    = Map2TpetraMap(*rowMap);
        RCP<const map_type>       tpetraColMap    = (colMap.is_null()    ? Teuchos::null : Map2TpetraMap(*colMap));
        const RCP<const map_type> tpetraRangeMap  = (rangeMap.is_null()  ? tpetraRowMap  : Map2TpetraMap(*rangeMap));
        const RCP<const map_type> tpetraDomainMap = (domainMap.is_null() ? tpetraRowMap  : Map2TpetraMap(*domainMap));

        RCP<sparse_matrix_type> tA = reader_type::readSparseFile(fileName, tpetraRowMap, tpetraColMap, tpetraDomainMap, tpetraRangeMap,
                                                                 callFillComplete, tolerant, debug);
        if (tA.is_null())
          throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from readSparseFile() is null.");

        RCP<TpetraCrsMatrix> tmpA1 = rcp(new TpetraCrsMatrix(tA));
        RCP<CrsMatrix>       tmpA2 = rcp_implicit_cast<CrsMatrix>(tmpA1);
        RCP<Matrix>          A     = rcp(new CrsMatrixWrap(tmpA2));

        return A;
#else
        throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
      } else {
        throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
      }
    } else {
      // Custom file format (binary)
      std::ifstream ifs(fileName.c_str(), std::ios::binary);
      TEUCHOS_TEST_FOR_EXCEPTION(!ifs.good(), Exceptions::RuntimeError, "Can not read \"" << fileName << "\"");
      int m, n, nnz;
      ifs.read(reinterpret_cast<char*>(&m),   sizeof(m));
      ifs.read(reinterpret_cast<char*>(&n),   sizeof(n));
      ifs.read(reinterpret_cast<char*>(&nnz), sizeof(nnz));

      RCP<Matrix> A = MatrixFactory::Build(rowMap, colMap, 1);

      TEUCHOS_TEST_FOR_EXCEPTION(sizeof(int) != sizeof(GO), Exceptions::RuntimeError, "Incompatible sizes");

      Teuchos::ArrayView<const GO> rowElements = rowMap->getNodeElementList();
      Teuchos::ArrayView<const GO> colElements = colMap->getNodeElementList();

      Teuchos::Array<GO> inds;
      Teuchos::Array<SC> vals;
      for (int i = 0; i < m; i++) {
        int row, rownnz;
        ifs.read(reinterpret_cast<char*>(&row),    sizeof(row));
        ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));
        inds.resize(rownnz);
        vals.resize(rownnz);
        for (int j = 0; j < rownnz; j++) {
          int index;
          ifs.read(reinterpret_cast<char*>(&index), sizeof(index));
          inds[j] = colElements[Teuchos::as<LO>(index)];
        }
        for (int j = 0; j < rownnz; j++) {
          double value;
          ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
          vals[j] = Teuchos::as<SC>(value);
        }
        A->insertGlobalValues(rowElements[row], inds, vals);
      }
      A->fillComplete(domainMap, rangeMap);
      return A;
    }

    return Teuchos::null;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(const std::string& fileName, const MultiVector& x) {

    std::string mapfile = "map_" + fileName;
    Write(mapfile, *(x.getMap()));

    RCP<const MultiVector> tmp_Vec = rcpFromRef(x);
#ifdef HAVE_MUELU_EPETRA
    const RCP<const EpetraMultiVector>& tmp_EVec = rcp_dynamic_cast<const EpetraMultiVector>(tmp_Vec);
    if (tmp_EVec != Teuchos::null) {
#ifdef HAVE_MUELU_EPETRAEXT
      int rv = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(), *(tmp_EVec->getEpetra_MultiVector()));
      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::RowMatrixToMatrixMarketFile return value of " + toString(rv));
#else
      throw Exceptions::RuntimeError("Compiled without EpetraExt");
#endif
      return;
    }
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraMultiVector> &tmp_TVec = rcp_dynamic_cast<const TpetraMultiVector>(tmp_Vec);
    if (tmp_TVec != Teuchos::null) {
      RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > TVec = tmp_TVec->getTpetra_MultiVector();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeDenseFile(fileName, TVec);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in multivector writing");

  } //Write

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(const std::string& fileName, const Map& M) {
    RCP<const Map> tmp_Map = rcpFromRef(M);
#ifdef HAVE_MUELU_EPETRAEXT
    const RCP<const Xpetra::EpetraMap>& tmp_EMap = rcp_dynamic_cast<const Xpetra::EpetraMap>(tmp_Map);
    if (tmp_EMap != Teuchos::null) {
#ifdef HAVE_MUELU_EPETRAEXT
      int rv = EpetraExt::BlockMapToMatrixMarketFile(fileName.c_str(), tmp_EMap->getEpetra_Map());
      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::BlockMapToMatrixMarketFile() return value of " + toString(rv));
#else
      throw(Exceptions::RuntimeError("Compiled without EpetraExt"));
#endif
      return;
    }
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraMap> &tmp_TMap = rcp_dynamic_cast<const TpetraMap>(tmp_Map);
    if (tmp_TMap != Teuchos::null) {
      RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > TMap = tmp_TMap->getTpetra_Map();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeMapFile(fileName, *TMap);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw Exceptions::BadCast("Could not cast to EpetraMap or TpetraMap in map writing");

  } //Write

#ifndef _WIN32
#include <unistd.h>

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PauseForDebugger() {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    int myPID = comm->getRank();
    int pid   = getpid();

    char hostname[80];
    for (int i = 0; i <comm->getSize(); i++) {
      if (i == myPID) {
        gethostname(hostname, sizeof(hostname));
        std::cout << "Host: " << hostname << "\tMPI rank: " << myPID << ",\tPID: " << pid << "\n\tattach " << pid << std::endl;
        sleep(1);
      }
    }

    if (myPID == 0) {
      std::cout << "** Enter a character to continue > " << std::endl;
      char go = ' ';
      int r = scanf("%c", &go);
      (void)r;
      assert(r > 0);
    }
    comm->barrier();
  } //PauseForDebugger
#else
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PauseForDebugger() {
      throw(Exceptions::RuntimeError("MueLu Utils: PauseForDebugger not implemented on Windows."));
  }

#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  PowerMethod(const Matrix& A, bool scaleByDiag, LO niters, Magnitude tolerance, bool verbose, unsigned int seed) {
    TEUCHOS_TEST_FOR_EXCEPTION(!(A.getRangeMap()->isSameAs(*(A.getDomainMap()))), Exceptions::Incompatible,
        "Utils::PowerMethod: operator must have domain and range maps that are equivalent.");

    // Create three vectors, fill z with random numbers
    RCP<Vector> q = VectorFactory::Build(A.getDomainMap());
    RCP<Vector> r = VectorFactory::Build(A.getRangeMap());
    RCP<Vector> z = VectorFactory::Build(A.getRangeMap());

    z->setSeed(seed);  // seed random number generator
    z->randomize(true);// use Xpetra implementation: -> same results for Epetra and Tpetra

    Teuchos::Array<Magnitude> norms(1);

    typedef Teuchos::ScalarTraits<SC> STS;

    const SC zero = STS::zero(), one = STS::one();

    SC lambda = zero;
    Magnitude residual = STS::magnitude(zero);

    // power iteration
    RCP<Vector> diagInvVec;
    if (scaleByDiag) {
      RCP<Vector> diagVec = VectorFactory::Build(A.getRowMap());
      A.getLocalDiagCopy(*diagVec);
      diagInvVec = VectorFactory::Build(A.getRowMap());
      diagInvVec->reciprocal(*diagVec);
    }

    for (int iter = 0; iter < niters; ++iter) {
      z->norm2(norms);                                  // Compute 2-norm of z
      q->update(one/norms[0], *z, zero);                // Set q = z / normz
      A.apply(*q, *z);                                  // Compute z = A*q
      if (scaleByDiag)
        z->elementWiseMultiply(one, *diagInvVec, *z, zero);
      lambda = q->dot(*z);                              // Approximate maximum eigenvalue: lamba = dot(q,z)

      if (iter % 100 == 0 || iter + 1 == niters) {
        r->update(1.0, *z, -lambda, *q, zero);          // Compute A*q - lambda*q
        r->norm2(norms);
        residual = STS::magnitude(norms[0] / lambda);
        if (verbose) {
          std::cout << "Iter = " << iter
                    << "  Lambda = " << lambda
                    << "  Residual of A*q - lambda*q = " << residual
                    << std::endl;
        }
      }
      if (residual < tolerance)
        break;
    }

    return lambda;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const SC>& scalingVector, bool doInverse,
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
        Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;

      default:
        throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled.");
        break;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
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
          throw Exceptions::RuntimeError("In Utils::Scaling: cannot fillComplete because the domain and/or range map hasn't been defined");

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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Teuchos::FancyOStream> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MakeFancy(std::ostream& os) {
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(os));
    return fancy;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Distance2(const MultiVector& v, LocalOrdinal i0, LocalOrdinal i1) {
    size_t numVectors = v.getNumVectors();

    Scalar d = Teuchos::ScalarTraits<Scalar>::zero();
    for (size_t j = 0; j < numVectors; j++) {
      Teuchos::ArrayRCP<const Scalar> vv = v.getData(j);
      d += (vv[i0] - vv[i1])*(vv[i0] - vv[i1]);
    }

    return Teuchos::ScalarTraits<SC>::magnitude(d);
  }

  template <class SC, class LO, class GO, class NO>
  ArrayRCP<const bool> Utils<SC, LO, GO, NO>::DetectDirichletRows(const Matrix& A, const typename Teuchos::ScalarTraits<SC>::magnitudeType& tol) {
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

  //pulled directly from ml_utils.cpp
  template <class SC, class LO, class GO, class NO>
  void Utils<SC, LO, GO, NO>::SetRandomSeed(const Teuchos::Comm<int> &comm) {
    // Distribute the seeds evenly in [1,maxint-1].  This guarantees nothing
    // about where in random number stream we are, but avoids overflow situations
    // in parallel when multiplying by a PID.  It would be better to use
    // a good parallel random number generator.

    double one = 1.0;
    int maxint = INT_MAX; //= 2^31-1 = 2147483647 for 32-bit integers
    int mySeed = Teuchos::as<int>((maxint-1) * (one -(comm.getRank()+1)/(comm.getSize()+one)) );
    if (mySeed < 1 || mySeed == maxint) {
      std::ostringstream errStr;
      errStr << "Error detected with random seed = " << mySeed << ". It should be in the interval [1,2^31-2].";
      throw Exceptions::RuntimeError(errStr.str());
    }

    std::srand(mySeed);

    // For Tpetra, we could use Kokkos' random number generator here.
    Teuchos::ScalarTraits<SC>::seedrandom(mySeed);

    // Epetra
    //   MultiVector::Random() -> Epetra_Util::RandomDouble() -> Epetra_Utils::RandomInt()
    // Its own random number generator, based on Seed_. Seed_ is initialized in Epetra_Util constructor with std::rand()
    // So our setting std::srand() affects that too
  }

  template <class SC, class LO, class GO, class NO>
  void Utils<SC, LO, GO, NO>::findDirichletRows(Teuchos::RCP<Matrix> A,
                                                     std::vector<LO>& dirichletRows) {
    dirichletRows.resize(0);
    for(size_t i=0; i<A->getNodeNumRows(); i++) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> values;
      A->getLocalRowView(i,indices,values);
      int nnz=0;
      for (int j=0; j<indices.size(); j++) {
        // FIXME (mfh 12 Sep 2015) I just replaced abs with the
        // appropriate ScalarTraits call.  However, this is NOT
        // correct for arbitrary scalar types!!!  I'm guessing you
        // should use the equivalent of LAPACK's SFMIN or machine
        // epsilon here.
        if (Teuchos::ScalarTraits<SC>::magnitude(values[j]) > 1.0e-16) {
          nnz++;
        }
      }
      if (nnz == 1 || nnz == 2) {
        dirichletRows.push_back(i);
      }
    }
  }

  template<class SC, class LO, class GO, class NO>
  void Utils<SC, LO, GO, NO>::findDirichletCols(Teuchos::RCP<Matrix> A,
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
      if(Teuchos::ScalarTraits<SC>::magnitude(myCols[i])>0.0)
        dirichletCols[i]=1;
      else
        dirichletCols[i]=0;
    }
  }

  template<class SC, class LO, class GO, class NO>
  void Utils<SC, LO, GO, NO>::Apply_BCsToMatrixRows(Teuchos::RCP<Matrix>& A,
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
  void Utils<SC, LO, GO, NO>::Apply_BCsToMatrixCols(Teuchos::RCP<Matrix>& A,
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
  void Utils<SC, LO, GO, NO>::Remove_Zeroed_Rows(Teuchos::RCP<Matrix>& A,
                                                      double tol) {
    Teuchos::RCP<const Map> rowMap = A->getRowMap();
    RCP<Matrix> DiagMatrix = MatrixFactory::Build(rowMap,1);
    RCP<Matrix> NewMatrix = MatrixFactory::Build(rowMap,1);
    for(size_t i=0; i<A->getNodeNumRows(); i++) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> values;
      A->getLocalRowView(i,indices,values);
      int nnz=0;
      for (int j=0; j<indices.size(); j++) {
        if (Teuchos::ScalarTraits<SC>::magnitude(values[j]) > tol) {
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
    Xpetra::MatrixMatrix<SC,LO,GO,NO>::TwoMatrixAdd(*DiagMatrix,false,(SC)1.0,*A,false,(SC)1.0,NewMatrix,*out);
    NewMatrix->fillComplete();
    A=NewMatrix;

  }



  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Transpose (Matrix& Op, bool optimizeTranspose,const std::string & label) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    std::string TorE = "epetra";
#else
    std::string TorE = "tpetra";
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    try {
      const Epetra_CrsMatrix& epetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(Op);
      (void) epetraOp; // silence "unused variable" compiler warning
    } catch (...) {
      TorE = "tpetra";
    }
#endif

#ifdef HAVE_MUELU_TPETRA
    if (TorE == "tpetra") {
      try {
        const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(Op);

        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
        Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp),label); //more than meets the eye
        A = transposer.createTranspose();

        RCP<TpetraCrsMatrix> AA = rcp(new TpetraCrsMatrix(A) );
        RCP<CrsMatrix> AAA = rcp_implicit_cast<CrsMatrix>(AA);
        RCP<Matrix> AAAA = rcp( new CrsMatrixWrap(AAA) );
        if (!AAAA->isFillComplete())
          AAAA->fillComplete(Op.getRangeMap(),Op.getDomainMap());

        return AAAA;

      } catch (std::exception& e) {
        std::cout << "threw exception '" << e.what() << "'" << std::endl;
        throw Exceptions::RuntimeError("Utils::Transpose failed, perhaps because matrix is not a Crs matrix");
      }
    } //if
#endif

    // Epetra case
    std::cout << "Utilities::Transpose() not implemented for Epetra" << std::endl;
    return Teuchos::null;

  } // Transpose

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Epetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MyOldScalematrix and Epetra cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
  }
} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DEF_HPP

//  LocalWords:  LocalOrdinal
