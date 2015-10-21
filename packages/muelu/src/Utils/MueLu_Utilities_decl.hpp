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

#include <unistd.h> //necessary for "sleep" function in debugging methods
#include <string>

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MUELU_TPETRA
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#endif
#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_ExportFactory.hpp>

#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraCrsMatrix_fwd.hpp>

// needed because of inlined function
//TODO: remove inline function?
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#endif

#include "MueLu_Exceptions.hpp"

#ifdef HAVE_MUELU_EPETRAEXT
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Vector;
#include "EpetraExt_Transpose_RowMatrix.h"
#endif

#ifdef HAVE_MUELU_TPETRA
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Xpetra_TpetraCrsMatrix_fwd.hpp>
#include <Xpetra_TpetraMultiVector_fwd.hpp>
#endif

#include <MueLu_UtilitiesBase.hpp>


namespace MueLu {
// MPI helpers
#define MueLu_sumAll(rcpComm, in, out)                                        \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out))
#define MueLu_minAll(rcpComm, in, out)                                        \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out))
#define MueLu_maxAll(rcpComm, in, out)                                        \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out))

#ifdef HAVE_MUELU_EPETRA
//defined after Utilities class
template<typename SC,typename LO,typename GO,typename NO>
RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO> >
Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB);

template<typename SC,typename LO,typename GO,typename NO>
RCP<Xpetra::Matrix<SC, LO, GO, NO> >
EpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Epetra_CrsMatrix>& A);

template<typename SC,typename LO,typename GO,typename NO>
RCP<Xpetra::MultiVector<SC, LO, GO, NO> >
EpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Epetra_MultiVector>& V);
#endif

#ifdef HAVE_MUELU_TPETRA
template<typename SC,typename LO,typename GO,typename NO>
RCP<Xpetra::Matrix<SC, LO, GO, NO> >
TpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> >& Atpetra);


template<typename SC,typename LO,typename GO,typename NO>
RCP<Xpetra::MultiVector<SC, LO, GO, NO> >
TpetraCrs_To_XpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> >& Vtpetra);
#endif

/*!
    @class Utilities
    @brief MueLu utility class.

    This class provides a number of static helper methods. Some are temporary and will eventually
    go away, while others should be moved to Xpetra.
 */
template <class Scalar,
class LocalOrdinal  = int,
class GlobalOrdinal = LocalOrdinal,
class Node          = KokkosClassic::DefaultNode::DefaultNodeType>
class Utilities : public UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_UTILITIES_SHORT
  //#include "MueLu_UseShortNames.hpp"

public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;

#ifdef HAVE_MUELU_EPETRA
  //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
  // @{
  static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > const Vec);
  static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Vec);

  static const Epetra_MultiVector&                        MV2EpetraMV(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Vec);
  static       Epetra_MultiVector&                        MV2NonConstEpetraMV(Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Vec);

  static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);
  static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);

  static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);
  static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);

  static const Epetra_Map&                                Map2EpetraMap(const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& map);
  // @}
#endif

#ifdef HAVE_MUELU_TPETRA
  //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
  // @{
  static RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2TpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > const Vec);
  static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2NonConstTpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Vec);
  static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2NonConstTpetraMV2(Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Vec);

  static const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&      MV2TpetraMV(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Vec);
  static       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&      MV2NonConstTpetraMV(Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Vec);

  static RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraCrs(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);
  static RCP<      Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraCrs(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);

  static const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2TpetraCrs(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);
  static       Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2NonConstTpetraCrs(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);

  static RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraRow(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);
  static RCP<      Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraRow(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);


  static const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >        Map2TpetraMap(const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& map);
#endif

  static RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >      Crs2Op(RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op) { return UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(Op); }
  static Teuchos::ArrayRCP<Scalar>                                         GetMatrixDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonal(A); }
  static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >      GetMatrixDiagonalInverse(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonalInverse(A,tol); }
  static Teuchos::ArrayRCP<Scalar>                                         GetLumpedMatrixDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetLumpedMatrixDiagonal(A); }
  static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >      GetMatrixOverlappedDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixOverlappedDiagonal(A); }
  static Teuchos::Array<Magnitude>                                         ResidualNorm(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ResidualNorm(Op,X,RHS); }
  static RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Residual(Op,X,RHS); }
  static void                                                              PauseForDebugger() { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PauseForDebugger(); }
  static RCP<Teuchos::FancyOStream>                                        MakeFancy(std::ostream& os) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MakeFancy(os); }
  static typename Teuchos::ScalarTraits<Scalar>::magnitudeType             Distance2(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& v, LocalOrdinal i0, LocalOrdinal i1) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Distance2(v,i0,i1); }
  static Teuchos::ArrayRCP<const bool>                                     DetectDirichletRows(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletRows(A,tol); }
  static void                                                              SetRandomSeed(const Teuchos::Comm<int> &comm) { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetRandomSeed(comm); }

  static Scalar PowerMethod(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, bool scaleByDiag = true,
      LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
    return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PowerMethod(A,scaleByDiag,niters,tolerance,verbose,seed);
  }

  static void MyOldScaleMatrix(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Teuchos::ArrayRCP<const Scalar>& scalingVector, bool doInverse = true,
      bool doFillComplete = true, bool doOptimizeStorage = true);

  static void MyOldScaleMatrix_Epetra(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
      bool doFillComplete, bool doOptimizeStorage);
  static void MyOldScaleMatrix_Tpetra(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
      bool doFillComplete, bool doOptimizeStorage);

  static RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Transpose(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, bool optimizeTranspose = false,const std::string & label = std::string());
}; // class Utilities

///////////////////////////////////////////

// This probably should be the double,int,int case
/*!
      @class Utilities
      @brief MueLu utility class (specialization SC=double and LO=GO=int).

      This class provides a number of static helper methods. Some are temporary and will eventually
      go away, while others should be moved to Xpetra.

      Note: this is the implementation for Epetra. Tpetra throws if TPETRA_INST_INT_INT is disabled!
 */
template <class Node>
class Utilities<double,int,int,Node> : public UtilitiesBase<double,int,int,Node> {
public:
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;

private:
  typedef Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsMatrixWrap;
  typedef Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsMatrix;
  typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Matrix;
  typedef Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Vector;
  typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MultiVector;
  typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
public:

#ifdef HAVE_MUELU_EPETRA
  //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
  // @{
  static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<MultiVector> const Vec) {
    RCP<const Xpetra::EpetraMultiVector > tmpVec = rcp_dynamic_cast<Xpetra::EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
    return tmpVec->getEpetra_MultiVector();
  }
  static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<MultiVector> Vec) {
    RCP<const Xpetra::EpetraMultiVector> tmpVec = rcp_dynamic_cast<Xpetra::EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
    return tmpVec->getEpetra_MultiVector();
  }

  static const Epetra_MultiVector&                        MV2EpetraMV(const MultiVector& Vec) {
    const Xpetra::EpetraMultiVector& tmpVec = dynamic_cast<const Xpetra::EpetraMultiVector&>(Vec);
    return *(tmpVec.getEpetra_MultiVector());
  }
  static       Epetra_MultiVector&                        MV2NonConstEpetraMV(MultiVector& Vec) {
    const Xpetra::EpetraMultiVector& tmpVec = dynamic_cast<const Xpetra::EpetraMultiVector&>(Vec);
    return *(tmpVec.getEpetra_MultiVector());
  }

  static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::EpetraCrsMatrix>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrix();
  }
  static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
  }

  static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Matrix& Op) {
    try {
      const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
      try {
        const Xpetra::EpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<const Xpetra::EpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getEpetra_CrsMatrix();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }
  static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Matrix& Op) {
    try {
      CrsMatrixWrap& crsOp = dynamic_cast<CrsMatrixWrap&>(Op);
      try {
        Xpetra::EpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<Xpetra::EpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getEpetra_CrsMatrixNonConst();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  static const Epetra_Map&                                Map2EpetraMap(const Map& map) {
    RCP<const Xpetra::EpetraMap> xeMap = rcp_dynamic_cast<const Xpetra::EpetraMap>(rcpFromRef(map));
    if (xeMap == Teuchos::null)
      throw Exceptions::BadCast("Utilities::Map2EpetraMap : Cast from Xpetra::Map to Xpetra::EpetraMap failed");
    return xeMap->getEpetra_Map();
  }
  // @}
#endif

#ifdef HAVE_MUELU_TPETRA
  //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
  // @{
  static RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2TpetraMV(RCP<MultiVector> const Vec)   {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    RCP<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmpVec = rcp_dynamic_cast<Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
    return tmpVec->getTpetra_MultiVector();
#else
    throw Exceptions::RuntimeError("MV2TpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  }
  static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2NonConstTpetraMV(RCP<MultiVector> Vec) {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    RCP<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmpVec = rcp_dynamic_cast<Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
    return tmpVec->getTpetra_MultiVector();
#else
    throw Exceptions::RuntimeError("MV2NonConstTpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#endif

  }
  static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2NonConstTpetraMV2(MultiVector& Vec)    {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(Vec);
    return tmpVec.getTpetra_MultiVector();
#else
    throw Exceptions::RuntimeError("MV2NonConstTpetraMV2: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  }

  static const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&      MV2TpetraMV(const MultiVector& Vec)   {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(Vec);
    return *(tmpVec.getTpetra_MultiVector());
#else
    throw Exceptions::RuntimeError("MV2TpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  }
  static       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&      MV2NonConstTpetraMV(MultiVector& Vec) {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(Vec);
    return *(tmpVec.getTpetra_MultiVector());
#else
    throw Exceptions::RuntimeError("MV2NonConstTpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  }

  static RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraCrs(RCP<const Matrix> Op)  {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    // Get the underlying Tpetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_CrsMatrix();
#else
    throw Exceptions::RuntimeError("Op2TpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  }
  static RCP<      Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraCrs(RCP<Matrix> Op){
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
#else
    throw Exceptions::RuntimeError("Op2NonConstTpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  };

  static const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2TpetraCrs(const Matrix& Op)   {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    try {
      const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
      try {
        const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& tmp_ECrsMtx = dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_CrsMatrix();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
#else
    throw Exceptions::RuntimeError("Op2TpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  }
  static       Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2NonConstTpetraCrs(Matrix& Op) {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    try {
      CrsMatrixWrap& crsOp = dynamic_cast<CrsMatrixWrap&>(Op);
      try {
        Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& tmp_ECrsMtx = dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_CrsMatrixNonConst();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
#else
    throw Exceptions::RuntimeError("Op2NonConstTpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  }

  static RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraRow(RCP<const Matrix> Op)   {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> crsMat = crsOp->getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmp_Crs = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(crsMat);
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
#else
    throw Exceptions::RuntimeError("Op2TpetraRow: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  }
  static RCP<      Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraRow(RCP<Matrix> Op) {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> crsMat = crsOp->getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmp_Crs = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(crsMat);
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
#else
    throw Exceptions::RuntimeError("Op2NonConstTpetraRow: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  };


  static const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >          Map2TpetraMap(const Map& map) {
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    const RCP<const Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>>& tmp_TMap = rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> >(rcpFromRef(map));
    if (tmp_TMap == Teuchos::null)
      throw Exceptions::BadCast("Utilities::Map2TpetraMap : Cast from Xpetra::Map to Xpetra::TpetraMap failed");
    return tmp_TMap->getTpetra_Map();
#else
    throw Exceptions::RuntimeError("Map2TpetraMap: Tpetra has not been compiled with support for LO=GO=int.");
#endif
  };
#endif

  static RCP<Matrix>                                                       Crs2Op(RCP<CrsMatrix> Op) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(Op); }
  static Teuchos::ArrayRCP<Scalar>                                         GetMatrixDiagonal(const Matrix& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonal(A); }
  static RCP<Vector>                                                       GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonalInverse(A,tol); }
  static Teuchos::ArrayRCP<Scalar>                                         GetLumpedMatrixDiagonal(const Matrix& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetLumpedMatrixDiagonal(A); }
  static RCP<Vector>                                                       GetMatrixOverlappedDiagonal(const Matrix& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixOverlappedDiagonal(A); }
  static Teuchos::Array<Magnitude>                                         ResidualNorm(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ResidualNorm(Op,X,RHS); }
  static RCP<MultiVector>                                                  Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Residual(Op,X,RHS); }
  static void                                                              PauseForDebugger() { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PauseForDebugger(); }
  static RCP<Teuchos::FancyOStream>                                        MakeFancy(std::ostream& os) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MakeFancy(os); }
  static typename Teuchos::ScalarTraits<Scalar>::magnitudeType             Distance2(const MultiVector& v, LocalOrdinal i0, LocalOrdinal i1) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Distance2(v,i0,i1); }
  static Teuchos::ArrayRCP<const bool>                                     DetectDirichletRows(const Matrix& A, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletRows(A,tol); }
  static void                                                              SetRandomSeed(const Teuchos::Comm<int> &comm) { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetRandomSeed(comm); }

  static Scalar PowerMethod(const Matrix& A, bool scaleByDiag = true,
      LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
    return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PowerMethod(A,scaleByDiag,niters,tolerance,verbose,seed);
  }

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
      break;
    }
  }

  // TODO This is the <double,int,int> specialization
  static void MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
      bool doFillComplete, bool doOptimizeStorage) {
#ifdef HAVE_MUELU_TPETRA
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
    try {
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpOp = Op2NonConstTpetraCrs(Op);

      const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap    = tpOp.getRowMap();
      const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > domainMap = tpOp.getDomainMap();
      const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rangeMap  = tpOp.getRangeMap();

      size_t maxRowSize = tpOp.getNodeMaxNumRowEntries();
      if (maxRowSize == Teuchos::as<size_t>(-1)) // hasn't been determined yet
        maxRowSize = 20;

      std::vector<Scalar> scaledVals(maxRowSize);
      if (tpOp.isFillComplete())
        tpOp.resumeFill();

      if (Op.isLocallyIndexed() == true) {
        Teuchos::ArrayView<const LocalOrdinal> cols;
        Teuchos::ArrayView<const Scalar> vals;

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
            Teuchos::ArrayView<const Scalar> valview(&scaledVals[0], nnz);
            tpOp.replaceLocalValues(i, cols, valview);
          }
        } //for (size_t i=0; ...

      } else {
        Teuchos::ArrayView<const GlobalOrdinal> cols;
        Teuchos::ArrayView<const Scalar> vals;

        for (size_t i = 0; i < rowMap->getNodeNumElements(); ++i) {
          GlobalOrdinal gid = rowMap->getGlobalElement(i);
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
            Teuchos::ArrayView<const Scalar> valview(&scaledVals[0], nnz);
            tpOp.replaceGlobalValues(gid, cols, valview);
          }
        } //for (size_t i=0; ...
      }

      if (doFillComplete) {
        if (domainMap == Teuchos::null || rangeMap == Teuchos::null)
          throw Exceptions::RuntimeError("In Utilities::Scaling: cannot fillComplete because the domain and/or range map hasn't been defined");

        RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
        params->set("Optimize Storage",    doOptimizeStorage);
        params->set("No Nonlocal Changes", true);
        Op.fillComplete(Op.getDomainMap(), Op.getRangeMap(), params);
      }
    } catch(...) {
      throw Exceptions::RuntimeError("Only Tpetra::CrsMatrix types can be scaled (Err.1)");
    }
#else
    throw Exceptions::RuntimeError("Matrix scaling is not possible because Tpetra has not been compiled with support for LO=GO=int.");
#endif
#else
    throw Exceptions::RuntimeError("Matrix scaling is not possible because Tpetra has not been enabled.");
#endif
  }

  static void MyOldScaleMatrix_Epetra (Matrix& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector, bool doFillComplete, bool doOptimizeStorage) {
#ifdef HAVE_MUELU_EPETRA
    try {
      //const Epetra_CrsMatrix& epOp = Utilities<double,int,int>::Op2NonConstEpetraCrs(Op);
      const Epetra_CrsMatrix& epOp = Op2NonConstEpetraCrs(Op);

      Epetra_Map const &rowMap = epOp.RowMap();
      int nnz;
      double *vals;
      int *cols;

      for (int i = 0; i < rowMap.NumMyElements(); ++i) {
        epOp.ExtractMyRowView(i, nnz, vals, cols);
        for (int j = 0; j < nnz; ++j)
          vals[j] *= scalingVector[i];
      }

    } catch (...){
      throw Exceptions::RuntimeError("Only Epetra_CrsMatrix types can be scaled");
    }
#else
    throw Exceptions::RuntimeError("Matrix scaling is not possible because Epetra has not been enabled.");
#endif // HAVE_MUELU_EPETRA
  }

  /*! @brief Transpose a Xpetra::Matrix

    Note: Currently, an error is thrown if the matrix isn't a Tpetra::CrsMatrix or Epetra_CrsMatrix.
    In principle, however, we could allow any Epetra_RowMatrix because the Epetra transposer does.
   */
  static RCP<Matrix> Transpose(Matrix& Op, bool optimizeTranspose = false,const std::string & label = std::string()) {
    switch (Op.getRowMap()->lib()) {
    case Xpetra::UseTpetra:
    {
#ifdef HAVE_MUELU_TPETRA
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
      try {
        const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpetraOp = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(Op);

        // Compute the transpose A of the Tpetra matrix tpetraOp.
        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
        Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp),label);
        A = transposer.createTranspose();
        RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AA   = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A));
        RCP<CrsMatrix>                                                           AAA  = rcp_dynamic_cast<CrsMatrix>(AA);
        RCP<CrsMatrixWrap>                                                       AAAA = rcp( new CrsMatrixWrap(AAA));

        return AAAA;
      }
      catch (std::exception& e) {
        std::cout << "threw exception '" << e.what() << "'" << std::endl;
        throw Exceptions::RuntimeError("Utilities::Transpose failed, perhaps because matrix is not a Crs matrix");
      }
#else
      throw Exceptions::RuntimeError("Utilities::Transpose: Tpetra is not compiled with LO=GO=int. Add TPETRA_INST_INT_INT:BOOL=ON to your configuration!");
#endif
#else
      throw Exceptions::RuntimeError("Utilities::Transpose: Tpetra is not compiled!");
#endif
      break;
    }
    case Xpetra::UseEpetra:
    {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("ZZ Entire Transpose"));
      // Epetra case
      Epetra_CrsMatrix& epetraOp = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(Op);
      EpetraExt::RowMatrix_Transpose transposer;
      Epetra_CrsMatrix * A = dynamic_cast<Epetra_CrsMatrix*>(&transposer(epetraOp));
      transposer.ReleaseTranspose(); // So we can keep A in Muelu...

      RCP<Epetra_CrsMatrix> rcpA(A);
      RCP<Xpetra::EpetraCrsMatrix>            AA   = rcp(new Xpetra::EpetraCrsMatrix(rcpA));
      RCP<CrsMatrix>                          AAA  = rcp_dynamic_cast<CrsMatrix>(AA);
      RCP<CrsMatrixWrap>                      AAAA = rcp( new CrsMatrixWrap(AAA));
      AAAA->fillComplete(Op.getRangeMap(), Op.getDomainMap());

      return AAAA;
#else
      throw Exceptions::RuntimeError("Epetra (Err. 2)");
#endif
      break;
    }
    default:
      throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be transposed.");
      break;
    }

    return Teuchos::null;
  }


}; // class Utilities (specialization SC=double LO=GO=int)

///////////////////////////////////////////



/*! Removes the following non-serializable data (A,P,R,Nullspace,Coordinates) from level-specific sublists from inList
    and moves it to nonSerialList.  Everything else is copied to serialList.  This function returns the level number of the highest level
    for which non-serializable data was provided.
 */
long ExtractNonSerializableData(const Teuchos::ParameterList& inList, Teuchos::ParameterList& serialList, Teuchos::ParameterList& nonSerialList);


/*! Tokenizes a (comma)-separated string, removing all leading and trailing whitespace
    WARNING: This routine is not threadsafe on most architectures
 */
void TokenizeStringAndStripWhiteSpace(const std::string & stream, std::vector<std::string> & tokenList, const char* token = ",");

/*! Returns true if a parameter name is a valid Muemex custom level variable, e.g. "MultiVector myArray"
 */
bool IsParamMuemexVariable(const std::string& name);

#ifdef HAVE_MUELU_EPETRA
//This non-member templated function exists so that the matrix-matrix multiply will compile if Epetra, Tpetra, and ML are enabled.
template<class SC,class LO,class GO,class NO>
RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO> >
Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap (RCP<Epetra_CrsMatrix> &epAB)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
  return Teuchos::null;
}

typedef KokkosClassic::DefaultNode::DefaultNodeType KDNT;

//specialization for the case of ScalarType=double and LocalOrdinal=GlobalOrdinal=int
template<>
inline RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<double,int,int,KDNT > (RCP<Epetra_CrsMatrix> &epAB) {
  RCP<Xpetra::EpetraCrsMatrix> tmpC1 = rcp(new Xpetra::EpetraCrsMatrix(epAB));
  RCP<Xpetra::CrsMatrix<double,int,int,KDNT> > tmpC2 = rcp_dynamic_cast<Xpetra::CrsMatrix<double,int,int,KDNT> >(tmpC1);
  RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT> > tmpC3 = rcp(new Xpetra::CrsMatrixWrap<double,int,int,KDNT>(tmpC2));
  return tmpC3;
}

/*! \fn EpetraCrs_To_XpetraMatrix
    @brief Helper function to convert a Epetra::CrsMatrix to an Xpetra::Matrix
    TODO move this function to an Xpetra utility file
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
EpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Epetra_CrsMatrix>& A) {
  typedef Xpetra::EpetraCrsMatrix                                            XECrsMatrix;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>       XCrsMatrix;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>   XCrsMatrixWrap;

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
  return rcp(new Xpetra::EpetraMultiVector(V));
}
#endif

#ifdef HAVE_MUELU_TPETRA
/*! \fn TpetraCrs_To_XpetraMatrix
    @brief Helper function to convert a Tpetra::CrsMatrix to an Xpetra::Matrix
    TODO move this function to an Xpetra utility file
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Atpetra) {
  typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XTCrsMatrix;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>       XCrsMatrix;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>   XCrsMatrixWrap;

  RCP<XCrsMatrix> Atmp = rcp(new XTCrsMatrix(Atpetra));
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
#endif

//! Little helper function to convert non-string types to strings
template<class T>
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

#ifdef HAVE_MUELU_TPETRA
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
#endif

} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DECL_HPP
