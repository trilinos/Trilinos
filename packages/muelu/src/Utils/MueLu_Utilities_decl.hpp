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
#ifndef MUELU_UTILITIES_DECL_HPP
#define MUELU_UTILITIES_DECL_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Xpetra_TpetraBlockCrsMatrix_fwd.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include <Xpetra_MatrixMatrix.hpp>

#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>

// needed because of inlined function
// TODO: remove inline function?
#include <Xpetra_EpetraCrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>

#endif

#include "MueLu_Exceptions.hpp"

#ifdef HAVE_MUELU_EPETRAEXT
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Vector;
#include "EpetraExt_Transpose_RowMatrix.h"
#endif

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#include <Tpetra_FECrsMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_FEMultiVector.hpp>
#include <Xpetra_TpetraCrsMatrix_fwd.hpp>
#include <Xpetra_TpetraMultiVector_fwd.hpp>

#include <MueLu_UtilitiesBase.hpp>

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
// defined after Utilities class
template <typename SC, typename LO, typename GO, typename NO>
RCP<Xpetra::CrsMatrixWrap<SC, LO, GO, NO> >
Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix>& epAB);

template <typename SC, typename LO, typename GO, typename NO>
RCP<Xpetra::Matrix<SC, LO, GO, NO> >
EpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Epetra_CrsMatrix>& A);

template <typename SC, typename LO, typename GO, typename NO>
RCP<Xpetra::MultiVector<SC, LO, GO, NO> >
EpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Epetra_MultiVector>& V);
#endif

template <typename SC, typename LO, typename GO, typename NO>
RCP<Xpetra::Matrix<SC, LO, GO, NO> >
TpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> >& Atpetra);

template <typename SC, typename LO, typename GO, typename NO>
RCP<Xpetra::Matrix<SC, LO, GO, NO> >
TpetraFECrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::FECrsMatrix<SC, LO, GO, NO> >& Atpetra);

template <typename SC, typename LO, typename GO, typename NO>
RCP<Xpetra::MultiVector<SC, LO, GO, NO> >
TpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> >& Vtpetra);

template <typename SC, typename LO, typename GO, typename NO>
RCP<Xpetra::MultiVector<SC, LO, GO, NO> >
TpetraFEMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Tpetra::FEMultiVector<SC, LO, GO, NO> >& Vtpetra);

template <typename SC, typename LO, typename GO, typename NO>
void leftRghtDofScalingWithinNode(const Xpetra::Matrix<SC, LO, GO, NO>& Atpetra, size_t blkSize, size_t nSweeps, Teuchos::ArrayRCP<SC>& rowScaling, Teuchos::ArrayRCP<SC>& colScaling);

/*!
  @class Utilities
  @brief MueLu utility class.

  This class provides a number of static helper methods. Some are temporary and will eventually
  go away, while others should be moved to Xpetra.
  */
template <class Scalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class Utilities : public UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_UTILITIES_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;

#ifdef HAVE_MUELU_EPETRA
  //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
  // @{
  static RCP<const Epetra_MultiVector> MV2EpetraMV(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > const vec);
  static RCP<Epetra_MultiVector> MV2NonConstEpetraMV(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > vec);

  static const Epetra_MultiVector& MV2EpetraMV(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec);
  static Epetra_MultiVector& MV2NonConstEpetraMV(Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec);

  static RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op);
  static RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op);

  static const Epetra_CrsMatrix& Op2EpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op);
  static Epetra_CrsMatrix& Op2NonConstEpetraCrs(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op);

  static const Epetra_Map& Map2EpetraMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map);
  // @}
#endif

  //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
  static RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > MV2TpetraMV(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > const vec);
  static RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > MV2NonConstTpetraMV(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > vec);
  static RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > MV2NonConstTpetraMV2(Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec);

  static const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& MV2TpetraMV(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec);
  static Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& MV2NonConstTpetraMV(Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec);

  static RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2TpetraCrs(RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op);
  static RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2NonConstTpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op);

  static const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op2TpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op);
  static Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op2NonConstTpetraCrs(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op);

  static RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2TpetraBlockCrs(RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op);
  static RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2NonConstTpetraBlockCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op);

  static const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op2TpetraBlockCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op);
  static Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op2NonConstTpetraBlockCrs(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op);

  static RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2TpetraRow(RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op);
  static RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2NonConstTpetraRow(RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op);

  static const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > Map2TpetraMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map);

  static void MyOldScaleMatrix(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const Teuchos::ArrayRCP<const Scalar>& scalingVector, bool doInverse = true,
                               bool doFillComplete = true, bool doOptimizeStorage = true);

  static void MyOldScaleMatrix_Epetra(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
                                      bool doFillComplete, bool doOptimizeStorage);
  static void MyOldScaleMatrix_Tpetra(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
                                      bool doFillComplete, bool doOptimizeStorage);

  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Transpose(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, bool optimizeTranspose = false, const std::string& label = std::string(), const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > RealValuedToScalarMultiVector(RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > X);

  static RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> > ExtractCoordinatesFromParameterList(ParameterList& paramList);

};  // class Utilities

///////////////////////////////////////////

#ifdef HAVE_MUELU_EPETRA
/*!
  @class Utilities
  @brief MueLu utility class (specialization SC=double and LO=GO=int).

  This class provides a number of static helper methods. Some are temporary and will eventually
  go away, while others should be moved to Xpetra.

Note: this is the implementation for Epetra. Tpetra throws if TPETRA_INST_INT_INT is disabled!
*/
template <>
class Utilities<double, int, int, Xpetra::EpetraNode> : public UtilitiesBase<double, int, int, Xpetra::EpetraNode> {
 public:
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef Xpetra::EpetraNode Node;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
#undef MUELU_UTILITIES_SHORT
#include "MueLu_UseShortNames.hpp"

 private:
  using EpetraMap         = Xpetra::EpetraMapT<GlobalOrdinal, Node>;
  using EpetraMultiVector = Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>;
  // using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node>;
 public:
  //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
  // @{
  static RCP<const Epetra_MultiVector> MV2EpetraMV(RCP<MultiVector> const vec) {
    RCP<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
    return tmpVec->getEpetra_MultiVector();
  }
  static RCP<Epetra_MultiVector> MV2NonConstEpetraMV(RCP<MultiVector> vec) {
    RCP<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
    return tmpVec->getEpetra_MultiVector();
  }

  static const Epetra_MultiVector& MV2EpetraMV(const MultiVector& vec) {
    const EpetraMultiVector& tmpVec = dynamic_cast<const EpetraMultiVector&>(vec);
    return *(tmpVec.getEpetra_MultiVector());
  }
  static Epetra_MultiVector& MV2NonConstEpetraMV(MultiVector& vec) {
    const EpetraMultiVector& tmpVec = dynamic_cast<const EpetraMultiVector&>(vec);
    return *(tmpVec.getEpetra_MultiVector());
  }

  static RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<const Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const EpetraCrsMatrix>& tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrix();
  }
  static RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const EpetraCrsMatrix>& tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
  }

  static const Epetra_CrsMatrix& Op2EpetraCrs(const Matrix& Op) {
    try {
      const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
      try {
        const EpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<const EpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getEpetra_CrsMatrix();
      } catch (std::bad_cast&) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      }
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }
  static Epetra_CrsMatrix& Op2NonConstEpetraCrs(Matrix& Op) {
    try {
      CrsMatrixWrap& crsOp = dynamic_cast<CrsMatrixWrap&>(Op);
      try {
        EpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<EpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getEpetra_CrsMatrixNonConst();
      } catch (std::bad_cast&) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      }
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  static const Epetra_Map& Map2EpetraMap(const Map& map) {
    RCP<const EpetraMap> xeMap = rcp_dynamic_cast<const EpetraMap>(rcpFromRef(map));
    if (xeMap == Teuchos::null)
      throw Exceptions::BadCast("Utilities::Map2EpetraMap : Cast from Xpetra::Map to Xpetra::EpetraMap failed");
    return xeMap->getEpetra_Map();
  }
  // @}

  //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
  static RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > MV2TpetraMV(RCP<MultiVector> const vec) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("MV2TpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#else
    RCP<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmpVec = rcp_dynamic_cast<Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
    return tmpVec->getTpetra_MultiVector();
#endif
  }
  static RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > MV2NonConstTpetraMV(RCP<MultiVector> vec) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("MV2NonConstTpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#else
    RCP<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmpVec = rcp_dynamic_cast<Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
    return tmpVec->getTpetra_MultiVector();
#endif
  }
  static RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > MV2NonConstTpetraMV2(MultiVector& vec) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("MV2NonConstTpetraMV2: Tpetra has not been compiled with support for LO=GO=int.");
#else
    const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(vec);
    return tmpVec.getTpetra_MultiVector();
#endif
  }

  static const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& MV2TpetraMV(const MultiVector& vec) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("MV2TpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#else
    const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(vec);
    return *(tmpVec.getTpetra_MultiVector());
#endif
  }
  static Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& MV2NonConstTpetraMV(MultiVector& vec) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("MV2NonConstTpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#else
    const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(vec);
    return *(tmpVec.getTpetra_MultiVector());
#endif
  }

  static RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2TpetraCrs(RCP<const Matrix> Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2TpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
    // Get the underlying Tpetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_CrsMatrix();
#endif
  }
  static RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2NonConstTpetraCrs(RCP<Matrix> Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2NonConstTpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
#endif
  };

  static const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op2TpetraCrs(const Matrix& Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2TpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
    try {
      const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
      try {
        const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_CrsMatrix();
      } catch (std::bad_cast&) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      }
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
#endif
  }
  static Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op2NonConstTpetraCrs(Matrix& Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2NonConstTpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
    try {
      CrsMatrixWrap& crsOp = dynamic_cast<CrsMatrixWrap&>(Op);
      try {
        Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_CrsMatrixNonConst();
      } catch (std::bad_cast&) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      }
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
#endif
  }

  static RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2TpetraBlockCrs(RCP<const Matrix> Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2TpetraBlockCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
    // Get the underlying Tpetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_BlockCrsMatrix();
#endif
  }

  static RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2NonConstTpetraBlockCrs(RCP<Matrix> Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2NonConstTpetraBlockCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_BlockCrsMatrixNonConst();
#endif
  };

  static const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op2TpetraBlockCrs(const Matrix& Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2TpetraBlockCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
    try {
      const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
      try {
        const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_BlockCrsMatrix();
      } catch (std::bad_cast&) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
      }
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
#endif
  }
  static Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op2NonConstTpetraBlockCrs(Matrix& Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2NonConstTpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
    try {
      CrsMatrixWrap& crsOp = dynamic_cast<CrsMatrixWrap&>(Op);
      try {
        Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_BlockCrsMatrixNonConst();
      } catch (std::bad_cast&) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
      }
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
#endif
  }

  static RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2TpetraRow(RCP<const Operator> Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2TpetraRow: Tpetra has not been compiled with support for LO=GO=int.");
#else
    RCP<const Matrix> mat = rcp_dynamic_cast<const Matrix>(Op);
    if (!mat.is_null()) {
      RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(mat);
      if (crsOp == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

      RCP<const CrsMatrix> crsMat                                                                  = crsOp->getCrsMatrix();
      const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_Crs = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(crsMat);
      RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_BlockCrs;
      if (!tmp_Crs.is_null()) {
        return tmp_Crs->getTpetra_CrsMatrixNonConst();
      } else {
        tmp_BlockCrs = rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(crsMat);
        if (tmp_BlockCrs.is_null())
          throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix and Xpetra::TpetraBlockCrsMatrix failed");
        return tmp_BlockCrs->getTpetra_BlockCrsMatrixNonConst();
      }
    } else {
      RCP<const TpetraOperator> tpOp                                                = rcp_dynamic_cast<const TpetraOperator>(Op, true);
      RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tOp   = tpOp->getOperatorConst();
      RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tRow = rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tOp, true);
      return tRow;
    }
#endif
  }

  static RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op2NonConstTpetraRow(RCP<Operator> Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Op2NonConstTpetraRow: Tpetra has not been compiled with support for LO=GO=int.");
#else
    RCP<Matrix> mat = rcp_dynamic_cast<Matrix>(Op);
    if (!mat.is_null()) {
      RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(mat);
      if (crsOp == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

      RCP<const CrsMatrix> crsMat                                                                  = crsOp->getCrsMatrix();
      const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_Crs = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(crsMat);
      RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tmp_BlockCrs;
      if (!tmp_Crs.is_null()) {
        return tmp_Crs->getTpetra_CrsMatrixNonConst();
      } else {
        tmp_BlockCrs = rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(crsMat);
        if (tmp_BlockCrs.is_null())
          throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix and Xpetra::TpetraBlockCrsMatrix failed");
        return tmp_BlockCrs->getTpetra_BlockCrsMatrixNonConst();
      }
    } else {
      RCP<TpetraOperator> tpOp                                                = rcp_dynamic_cast<TpetraOperator>(Op, true);
      RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tOp   = tpOp->getOperator();
      RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tRow = rcp_dynamic_cast<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tOp, true);
      return tRow;
    }
#endif
  };

  static const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > Map2TpetraMap(const Map& map) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Map2TpetraMap: Tpetra has not been compiled with support for LO=GO=int.");
#else
    const RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >& tmp_TMap = rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >(rcpFromRef(map));
    if (tmp_TMap == Teuchos::null)
      throw Exceptions::BadCast("Utilities::Map2TpetraMap : Cast from Xpetra::Map to Xpetra::TpetraMap failed");
    return tmp_TMap->getTpetra_Map();
#endif
  };

  static void MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const Scalar>& scalingVector, bool doInverse = true,
                               bool doFillComplete = true, bool doOptimizeStorage = true) {
    Scalar one = Teuchos::ScalarTraits<Scalar>::one();
    Teuchos::ArrayRCP<Scalar> sv(scalingVector.size());
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
        MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;

      default:
        throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled.");
    }
  }

  // TODO This is the <double,int,int> specialization
  static void MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
                                      bool doFillComplete, bool doOptimizeStorage) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw Exceptions::RuntimeError("Matrix scaling is not possible because Tpetra has not been compiled with support for LO=GO=int.");
#else
    try {
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpOp = Op2NonConstTpetraCrs(Op);

      const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap    = tpOp.getRowMap();
      const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > domainMap = tpOp.getDomainMap();
      const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rangeMap  = tpOp.getRangeMap();

      size_t maxRowSize = tpOp.getLocalMaxNumRowEntries();
      if (maxRowSize == Teuchos::as<size_t>(-1))  // hasn't been determined yet
        maxRowSize = 20;

      std::vector<Scalar> scaledVals(maxRowSize);
      if (tpOp.isFillComplete())
        tpOp.resumeFill();

      if (Op.isLocallyIndexed() == true) {
        typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type cols;
        typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type vals;
        for (size_t i = 0; i < rowMap->getLocalNumElements(); ++i) {
          tpOp.getLocalRowView(i, cols, vals);
          size_t nnz = tpOp.getNumEntriesInLocalRow(i);
          if (nnz > maxRowSize) {
            maxRowSize = nnz;
            scaledVals.resize(maxRowSize);
          }
          for (size_t j = 0; j < nnz; ++j)
            scaledVals[j] = vals[j] * scalingVector[i];

          if (nnz > 0) {
            Teuchos::ArrayView<const LocalOrdinal> cols_view(cols.data(), nnz);
            Teuchos::ArrayView<const Scalar> valview(&scaledVals[0], nnz);
            tpOp.replaceLocalValues(i, cols_view, valview);
          }
        }  // for (size_t i=0; ...

      } else {
        typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type cols;
        typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type vals;

        for (size_t i = 0; i < rowMap->getLocalNumElements(); ++i) {
          GlobalOrdinal gid = rowMap->getGlobalElement(i);
          tpOp.getGlobalRowView(gid, cols, vals);
          size_t nnz = tpOp.getNumEntriesInGlobalRow(gid);
          if (nnz > maxRowSize) {
            maxRowSize = nnz;
            scaledVals.resize(maxRowSize);
          }
          // FIXME FIXME FIXME FIXME FIXME FIXME
          for (size_t j = 0; j < nnz; ++j)
            scaledVals[j] = vals[j] * scalingVector[i];  // FIXME i or gid?

          if (nnz > 0) {
            Teuchos::ArrayView<const LocalOrdinal> cols_view(cols.data(), nnz);
            Teuchos::ArrayView<const Scalar> valview(&scaledVals[0], nnz);
            tpOp.replaceGlobalValues(gid, cols_view, valview);
          }
        }  // for (size_t i=0; ...
      }

      if (doFillComplete) {
        if (domainMap == Teuchos::null || rangeMap == Teuchos::null)
          throw Exceptions::RuntimeError("In Utilities::Scaling: cannot fillComplete because the domain and/or range map hasn't been defined");

        RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
        params->set("Optimize Storage", doOptimizeStorage);
        params->set("No Nonlocal Changes", true);
        Op.fillComplete(Op.getDomainMap(), Op.getRangeMap(), params);
      }
    } catch (...) {
      throw Exceptions::RuntimeError("Only Tpetra::CrsMatrix types can be scaled (Err.1)");
    }
#endif
  }

  static void MyOldScaleMatrix_Epetra(Matrix& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector, bool /* doFillComplete */, bool /* doOptimizeStorage */) {
#ifdef HAVE_MUELU_EPETRA
    try {
      // const Epetra_CrsMatrix& epOp = Utilities<double,int,int>::Op2NonConstEpetraCrs(Op);
      const Epetra_CrsMatrix& epOp = Op2NonConstEpetraCrs(Op);

      Epetra_Map const& rowMap = epOp.RowMap();
      int nnz;
      double* vals;
      int* cols;

      for (int i = 0; i < rowMap.NumMyElements(); ++i) {
        epOp.ExtractMyRowView(i, nnz, vals, cols);
        for (int j = 0; j < nnz; ++j)
          vals[j] *= scalingVector[i];
      }

    } catch (...) {
      throw Exceptions::RuntimeError("Only Epetra_CrsMatrix types can be scaled");
    }
#else
    throw Exceptions::RuntimeError("Matrix scaling is not possible because Epetra has not been enabled.");
#endif  // HAVE_MUELU_EPETRA
  }

  /*! @brief Transpose a Xpetra::Matrix

      Note: Currently, an error is thrown if the matrix isn't a Tpetra::CrsMatrix or Epetra_CrsMatrix.
      In principle, however, we could allow any Epetra_RowMatrix because the Epetra transposer does.
  */
  static RCP<Matrix> Transpose(Matrix& Op, bool /* optimizeTranspose */ = false, const std::string& label = std::string(), const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
    switch (Op.getRowMap()->lib()) {
      case Xpetra::UseTpetra: {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
        throw Exceptions::RuntimeError("Utilities::Transpose: Tpetra is not compiled with LO=GO=int. Add TPETRA_INST_INT_INT:BOOL=ON to your configuration!");
#else
        using Helpers = Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
        /***************************************************************/
        if (Helpers::isTpetraCrs(Op)) {
          const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpetraOp = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(Op);

          // Compute the transpose A of the Tpetra matrix tpetraOp.
          RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
          Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp), label);

          {
            using Teuchos::ParameterList;
            using Teuchos::rcp;
            RCP<ParameterList> transposeParams = params.is_null() ? rcp(new ParameterList) : rcp(new ParameterList(*params));
            transposeParams->set("sort", false);
            A = transposer.createTranspose(transposeParams);
          }

          RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AA = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A));
          RCP<CrsMatrix> AAA                                                          = rcp_implicit_cast<CrsMatrix>(AA);
          RCP<Matrix> AAAA                                                            = rcp(new CrsMatrixWrap(AAA));

          if (Op.IsView("stridedMaps"))
            AAAA->CreateView("stridedMaps", Teuchos::rcpFromRef(Op), true /*doTranspose*/);

          return AAAA;
        }
        /***************************************************************/
        else if (Helpers::isTpetraBlockCrs(Op)) {
          using BCRS = Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
          // using CRS  = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
          const BCRS& tpetraOp = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraBlockCrs(Op);
          RCP<BCRS> At;
          {
            Tpetra::BlockCrsMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp), label);

            using Teuchos::ParameterList;
            using Teuchos::rcp;
            RCP<ParameterList> transposeParams = params.is_null() ? rcp(new ParameterList) : rcp(new ParameterList(*params));
            transposeParams->set("sort", false);
            At = transposer.createTranspose(transposeParams);
          }

          RCP<Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AA = rcp(new Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(At));
          RCP<CrsMatrix> AAA                                                               = rcp_implicit_cast<CrsMatrix>(AA);
          RCP<Matrix> AAAA                                                                 = rcp(new CrsMatrixWrap(AAA));

          if (Op.IsView("stridedMaps"))
            AAAA->CreateView("stridedMaps", Teuchos::rcpFromRef(Op), true /*doTranspose*/);

          return AAAA;

        }
        /***************************************************************/
        else {
          throw Exceptions::RuntimeError("Utilities::Transpose failed, perhaps because matrix is not a Crs or BlockCrs matrix");
        }
#endif
      }
      case Xpetra::UseEpetra: {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
        Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("ZZ Entire Transpose"));
        // Epetra case
        Epetra_CrsMatrix& epetraOp = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(Op);
        EpetraExt::RowMatrix_Transpose transposer;
        Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(&transposer(epetraOp));
        transposer.ReleaseTranspose();  // So we can keep A in Muelu...

        RCP<Epetra_CrsMatrix> rcpA(A);
        RCP<EpetraCrsMatrix> AA = rcp(new EpetraCrsMatrix(rcpA));
        RCP<CrsMatrix> AAA      = rcp_implicit_cast<CrsMatrix>(AA);
        RCP<Matrix> AAAA        = rcp(new CrsMatrixWrap(AAA));

        if (Op.IsView("stridedMaps"))
          AAAA->CreateView("stridedMaps", Teuchos::rcpFromRef(Op), true /*doTranspose*/);

        return AAAA;
#else
        throw Exceptions::RuntimeError("Epetra (Err. 2)");
#endif
      }
      default:
        throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be transposed.");
    }

    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

  static RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> >
  RealValuedToScalarMultiVector(RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > X) {
    RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Xscalar = rcp_dynamic_cast<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(X, true);
    return Xscalar;
  }

  /*! @brief Extract coordinates from parameter list and return them in a Xpetra::MultiVector
   */
  static RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> > ExtractCoordinatesFromParameterList(ParameterList& paramList) {
    RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> > coordinates = Teuchos::null;

    // check whether coordinates are contained in parameter list
    if (paramList.isParameter("Coordinates") == false)
      return coordinates;

#if (defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT)) || \
    (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))

      // define Tpetra::MultiVector type with Scalar=float only if
      // * ETI is turned off, since then the compiler will instantiate it automatically OR
      // * Tpetra is instantiated on Scalar=float
#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_TPETRA_INST_FLOAT)
    typedef Tpetra::MultiVector<float, LocalOrdinal, GlobalOrdinal, Node> tfMV;
    RCP<tfMV> floatCoords = Teuchos::null;
#endif

    // define Tpetra::MultiVector type with Scalar=double only if
    // * ETI is turned off, since then the compiler will instantiate it automatically OR
    // * Tpetra is instantiated on Scalar=double
    typedef Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> tdMV;
    RCP<tdMV> doubleCoords = Teuchos::null;
    if (paramList.isType<RCP<tdMV> >("Coordinates")) {
      // Coordinates are stored as a double vector
      doubleCoords = paramList.get<RCP<tdMV> >("Coordinates");
      paramList.remove("Coordinates");
    }
#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_TPETRA_INST_FLOAT)
    else if (paramList.isType<RCP<tfMV> >("Coordinates")) {
      // check if coordinates are stored as a float vector
      floatCoords = paramList.get<RCP<tfMV> >("Coordinates");
      paramList.remove("Coordinates");
      doubleCoords = rcp(new tdMV(floatCoords->getMap(), floatCoords->getNumVectors()));
      deep_copy(*doubleCoords, *floatCoords);
    }
#endif
    // We have the coordinates in a Tpetra double vector
    if (doubleCoords != Teuchos::null) {
      coordinates = Teuchos::rcp(new Xpetra::TpetraMultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>(doubleCoords));
      TEUCHOS_TEST_FOR_EXCEPT(doubleCoords->getNumVectors() != coordinates->getNumVectors());
    }
#endif  // Tpetra instantiated on GO=int and EpetraNode

#if defined(HAVE_MUELU_EPETRA)
    RCP<Epetra_MultiVector> doubleEpCoords;
    if (paramList.isType<RCP<Epetra_MultiVector> >("Coordinates")) {
      doubleEpCoords = paramList.get<RCP<Epetra_MultiVector> >("Coordinates");
      paramList.remove("Coordinates");
      RCP<Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node> > epCoordinates = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>(doubleEpCoords));
      coordinates                                                         = rcp_dynamic_cast<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> >(epCoordinates);
      TEUCHOS_TEST_FOR_EXCEPT(doubleEpCoords->NumVectors() != Teuchos::as<int>(coordinates->getNumVectors()));
    }
#endif

    // check for Xpetra coordinates vector
    if (paramList.isType<decltype(coordinates)>("Coordinates")) {
      coordinates = paramList.get<decltype(coordinates)>("Coordinates");
    }

    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(coordinates));
    return coordinates;
  }

};  // class Utilities (specialization SC=double LO=GO=int)

#endif  // HAVE_MUELU_EPETRA

/*!
\brief Extract non-serializable data from level-specific sublists and move it to a separate parameter list

Look through the level-specific sublists form \c inList, extract non-serializable data and move it to \c nonSerialList.
Everything else is copied to the \c serialList.

\note Data is considered "non-serializable" if it is not the same on every rank/processor.

Non-serializable data to be moved:
- Operator "A"
- Prolongator "P"
- Restrictor "R"
- "M"
- "Mdiag"
- "K"
- Nullspace information "Nullspace"
- Coordinate information "Coordinates"
- "Node Comm"
- Primal-to-dual node mapping "DualNodeID2PrimalNodeID"
- "Primal interface DOF map"
- "pcoarsen: element to node map

@param[in] inList List with all input parameters/data as provided by the user
@param[out] serialList All serializable data from the input list
@param[out] nonSerialList All non-serializable, i.e. rank-specific data from the input list

@return This function returns the level number of the highest level for which non-serializable data was provided.

*/
long ExtractNonSerializableData(const Teuchos::ParameterList& inList, Teuchos::ParameterList& serialList, Teuchos::ParameterList& nonSerialList);

/*! Tokenizes a (comma)-separated string, removing all leading and trailing whitespace
WARNING: This routine is not threadsafe on most architectures
*/
void TokenizeStringAndStripWhiteSpace(const std::string& stream, std::vector<std::string>& tokenList, const char* token = ",");

/*! Returns true if a parameter name is a valid Muemex custom level variable, e.g. "MultiVector myArray"
 */
bool IsParamMuemexVariable(const std::string& name);

/*! Returns true if a parameter name is a valid user custom level variable, e.g. "MultiVector myArray"
 */
bool IsParamValidVariable(const std::string& name);

#ifdef HAVE_MUELU_EPETRA
/*! \fn EpetraCrs_To_XpetraMatrix
    @brief Helper function to convert a Epetra::CrsMatrix to an Xpetra::Matrix
    TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
EpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Epetra_CrsMatrix>& A) {
  typedef Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node> XECrsMatrix;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XCrsMatrix;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> XCrsMatrixWrap;

  RCP<XCrsMatrix> Atmp = rcp(new XECrsMatrix(A));
  return rcp(new XCrsMatrixWrap(Atmp));
}

/*! \fn EpetraMultiVector_To_XpetraMultiVector
  @brief Helper function to convert a Epetra::MultiVector to an Xpetra::MultiVector
  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
EpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Epetra_MultiVector>& V) {
  return rcp(new Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>(V));
}
#endif

/*! \fn TpetraCrs_To_XpetraMatrix
  @brief Helper function to convert a Tpetra::CrsMatrix to an Xpetra::Matrix
  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Atpetra) {
  typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XTCrsMatrix;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XCrsMatrix;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> XCrsMatrixWrap;

  RCP<XCrsMatrix> Atmp = rcp(new XTCrsMatrix(Atpetra));
  return rcp(new XCrsMatrixWrap(Atmp));
}

/*! \fn leftRghtDofScalingWithinNode
  @brief Helper function computes 2k left/right matrix scaling coefficients for PDE system with k x k blocks

  Heuristic algorithm computes rowScaling and colScaling so that one can effectively derive matrices
  rowScalingMatrix and colScalingMatrix such that the abs(rowsums) and abs(colsums) of

            rowScalingMatrix * Amat * colScalingMatrix

  are roughly constant. If D = diag(rowScalingMatrix), then

     D(i:blkSize:end) = rowScaling(i)   for i=1,..,blkSize .

  diag(colScalingMatrix) is defined analogously. This function only computes rowScaling/colScaling.
  You will need to copy them into a tpetra vector to use tpetra functions such as leftScale() and rightScale()
  via some kind of loop such as

  rghtScaleVec = Teuchos::rcp(new Tpetra::Vector<SC,LO,GO,NO>(tpetraMat->getColMap()));
  rghtScaleData  = rghtScaleVec->getDataNonConst(0);
  size_t itemp = 0;
  for (size_t i = 0; i < tpetraMat->getColMap()->getLocalNumElements(); i++) {
    rghtScaleData[i] = rghtDofPerNodeScale[itemp++];
    if (itemp == blkSize) itemp = 0;
  }
  followed by tpetraMat->rightScale(*rghtScaleVec);

  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void leftRghtDofScalingWithinNode(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Amat, size_t blkSize, size_t nSweeps, Teuchos::ArrayRCP<Scalar>& rowScaling, Teuchos::ArrayRCP<Scalar>& colScaling) {
  LocalOrdinal nBlks = (Amat.getRowMap()->getLocalNumElements()) / blkSize;

  Teuchos::ArrayRCP<Scalar> rowScaleUpdate(blkSize);
  Teuchos::ArrayRCP<Scalar> colScaleUpdate(blkSize);

  for (size_t i = 0; i < blkSize; i++) rowScaling[i] = 1.0;
  for (size_t i = 0; i < blkSize; i++) colScaling[i] = 1.0;

  for (size_t k = 0; k < nSweeps; k++) {
    LocalOrdinal row = 0;
    for (size_t i = 0; i < blkSize; i++) rowScaleUpdate[i] = 0.0;

    for (LocalOrdinal i = 0; i < nBlks; i++) {
      for (size_t j = 0; j < blkSize; j++) {
        Teuchos::ArrayView<const LocalOrdinal> cols;
        Teuchos::ArrayView<const Scalar> vals;
        Amat.getLocalRowView(row, cols, vals);

        for (size_t kk = 0; kk < Teuchos::as<size_t>(vals.size()); kk++) {
          size_t modGuy = (cols[kk] + 1) % blkSize;
          if (modGuy == 0) modGuy = blkSize;
          modGuy--;
          rowScaleUpdate[j] += rowScaling[j] * (Teuchos::ScalarTraits<Scalar>::magnitude(vals[kk])) * colScaling[modGuy];
        }
        row++;
      }
    }
    // combine information across processors
    Teuchos::ArrayRCP<Scalar> tempUpdate(blkSize);
    Teuchos::reduceAll(*(Amat.getRowMap()->getComm()), Teuchos::REDUCE_SUM, (LocalOrdinal)blkSize, rowScaleUpdate.getRawPtr(), tempUpdate.getRawPtr());
    for (size_t i = 0; i < blkSize; i++) rowScaleUpdate[i] = tempUpdate[i];

    /* We want to scale by sqrt(1/rowScaleUpdate), but we'll         */
    /* normalize things by the minimum rowScaleUpdate. That is, the  */
    /* largest scaling is always one (as normalization is arbitrary).*/

    Scalar minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude((rowScaleUpdate[0] / rowScaling[0]) / rowScaling[0]);

    for (size_t i = 1; i < blkSize; i++) {
      Scalar temp = (rowScaleUpdate[i] / rowScaling[i]) / rowScaling[i];
      if (Teuchos::ScalarTraits<Scalar>::magnitude(temp) < Teuchos::ScalarTraits<Scalar>::magnitude(minUpdate))
        minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude(temp);
    }
    for (size_t i = 0; i < blkSize; i++) rowScaling[i] *= sqrt(minUpdate / rowScaleUpdate[i]);

    row = 0;
    for (size_t i = 0; i < blkSize; i++) colScaleUpdate[i] = 0.0;

    for (LocalOrdinal i = 0; i < nBlks; i++) {
      for (size_t j = 0; j < blkSize; j++) {
        Teuchos::ArrayView<const LocalOrdinal> cols;
        Teuchos::ArrayView<const Scalar> vals;
        Amat.getLocalRowView(row, cols, vals);
        for (size_t kk = 0; kk < Teuchos::as<size_t>(vals.size()); kk++) {
          size_t modGuy = (cols[kk] + 1) % blkSize;
          if (modGuy == 0) modGuy = blkSize;
          modGuy--;
          colScaleUpdate[modGuy] += colScaling[modGuy] * (Teuchos::ScalarTraits<Scalar>::magnitude(vals[kk])) * rowScaling[j];
        }
        row++;
      }
    }
    Teuchos::reduceAll(*(Amat.getRowMap()->getComm()), Teuchos::REDUCE_SUM, (LocalOrdinal)blkSize, colScaleUpdate.getRawPtr(), tempUpdate.getRawPtr());
    for (size_t i = 0; i < blkSize; i++) colScaleUpdate[i] = tempUpdate[i];

    /* We want to scale by sqrt(1/colScaleUpdate), but we'll         */
    /* normalize things by the minimum colScaleUpdate. That is, the  */
    /* largest scaling is always one (as normalization is arbitrary).*/

    minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude((colScaleUpdate[0] / colScaling[0]) / colScaling[0]);

    for (size_t i = 1; i < blkSize; i++) {
      Scalar temp = (colScaleUpdate[i] / colScaling[i]) / colScaling[i];
      if (Teuchos::ScalarTraits<Scalar>::magnitude(temp) < Teuchos::ScalarTraits<Scalar>::magnitude(minUpdate))
        minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude(temp);
    }
    for (size_t i = 0; i < blkSize; i++) colScaling[i] *= sqrt(minUpdate / colScaleUpdate[i]);
  }
}

/*! \fn TpetraCrs_To_XpetraMatrix
  @brief Helper function to convert a Tpetra::FECrsMatrix to an Xpetra::Matrix
  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraFECrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Atpetra) {
  typedef typename Tpetra::FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::crs_matrix_type tpetra_crs_matrix_type;
  typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XTCrsMatrix;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XCrsMatrix;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> XCrsMatrixWrap;

  RCP<XCrsMatrix> Atmp = rcp(new XTCrsMatrix(rcp_dynamic_cast<tpetra_crs_matrix_type>(Atpetra)));
  return rcp(new XCrsMatrixWrap(Atmp));
}

/*! \fn TpetraMultiVector_To_XpetraMultiVector
  @brief Helper function to convert a Tpetra::MultiVector to an Xpetra::MultiVector
  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Vtpetra) {
  return rcp(new Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Vtpetra));
}

/*! \fn TpetraFEMultiVector_To_XpetraMultiVector
@brief Helper function to convert a Tpetra::FEMultiVector to an Xpetra::MultiVector
  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraFEMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Tpetra::FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Vtpetra) {
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
  RCP<const MV> Vmv = Teuchos::rcp_dynamic_cast<const MV>(Vtpetra);
  return rcp(new Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Vmv));
}

//! Little helper function to convert non-string types to strings
template <class T>
std::string toString(const T& what) {
  std::ostringstream buf;
  buf << what;
  return buf.str();
}

#ifdef HAVE_MUELU_EPETRA
/*! \fn EpetraCrs_To_XpetraMatrix
  @brief Helper function to convert a Epetra::CrsMatrix to an Xpetra::Matrix
  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
EpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Epetra_CrsMatrix>& A);

/*! \fn EpetraMultiVector_To_XpetraMultiVector
  @brief Helper function to convert a Epetra::MultiVector to an Xpetra::MultiVector
  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
EpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Epetra_MultiVector>& V);
#endif

/*! \fn TpetraCrs_To_XpetraMatrix
  @brief Helper function to convert a Tpetra::CrsMatrix to an Xpetra::Matrix
  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Atpetra);

/*! \fn TpetraMultiVector_To_XpetraMultiVector
  @brief Helper function to convert a Tpetra::MultiVector to an Xpetra::MultiVector
  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Vtpetra);

// Generates a communicator whose only members are other ranks of the baseComm on my node
Teuchos::RCP<const Teuchos::Comm<int> > GenerateNodeComm(RCP<const Teuchos::Comm<int> >& baseComm, int& NodeId, const int reductionFactor);

// Lower case string
std::string lowerCase(const std::string& s);

}  // namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif  // MUELU_UTILITIES_DECL_HPP
