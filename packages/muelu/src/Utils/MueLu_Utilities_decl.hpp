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
#endif

#ifdef HAVE_MUELU_TPETRA
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Xpetra_TpetraCrsMatrix_fwd.hpp>
#include <Xpetra_TpetraMultiVector_fwd.hpp>
#endif


namespace MueLu {
// MPI helpers
#define MueLu_sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out))
#define MueLu_minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out))
#define MueLu_maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out))

#ifdef HAVE_MUELU_EPETRA
  //defined after Utils class
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
    @class Utils
    @brief MueLu utility class.

    This class provides a number of static helper methods. Some are temporary and will eventually
    go away, while others should be moved to Xpetra.
  */
  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = KokkosClassic::DefaultNode::DefaultNodeType>
  class Utils {
#undef MUELU_UTILITIES_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

#ifdef HAVE_MUELU_EPETRA
    //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
    // @{
    static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<MultiVector> const Vec);
    static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<MultiVector> Vec);

    static const Epetra_MultiVector&                        MV2EpetraMV(const MultiVector& Vec);
    static       Epetra_MultiVector&                        MV2NonConstEpetraMV(MultiVector& Vec);

    static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Matrix> Op);
    static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Matrix> Op);

    static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Matrix& Op);
    static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Matrix& Op);

    static const Epetra_Map&                                Map2EpetraMap(const Map& map);
    // @}
#endif

#ifdef HAVE_MUELU_TPETRA
    //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
    // @{
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> >     MV2TpetraMV(RCP<MultiVector> const Vec);
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV(RCP<MultiVector> Vec);
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV2(MultiVector& Vec);

    static const Tpetra::MultiVector<SC,LO,GO,NO>&          MV2TpetraMV(const MultiVector& Vec);
    static       Tpetra::MultiVector<SC,LO,GO,NO>&          MV2NonConstTpetraMV(MultiVector& Vec);

    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2TpetraCrs(RCP<const Matrix> Op);
    static RCP<      Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraCrs(RCP<Matrix> Op);

    static const Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2TpetraCrs(const Matrix& Op);
    static       Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2NonConstTpetraCrs(Matrix& Op);

    static RCP<const Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2TpetraRow(RCP<const Matrix> Op);
    static RCP<      Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraRow(RCP<Matrix> Op);


    static const RCP<const Tpetra::Map<LO, GO, NO> >        Map2TpetraMap(const Map& map);
#endif

    static RCP<Xpetra::Matrix<SC,LO,GO,NO> >                Crs2Op(RCP<CrsMatrix> Op);




      /*! @brief Extract Matrix Diagonal

    Returns Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetMatrixDiagonal(const Matrix& A);

    /*! @brief Extract Matrix Diagonal

    Returns inverse of the Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = Teuchos::ScalarTraits<SC>::eps()*100);



    /*! @brief Extract Matrix Diagonal of lumped matrix

    Returns Matrix diagonal of lumped matrix in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetLumpedMatrixDiagonal(const Matrix& A);

    /*! @brief Extract Overlapped Matrix Diagonal

    Returns overlapped Matrix diagonal in ArrayRCP.

    The local overlapped diagonal has an entry for each index in A's column map.
    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixOverlappedDiagonal(const Matrix& A);

    /*! @brief Left scale matrix by an arbitrary vector.

    Algorithmically, this left scales a matrix by a diagonal matrix.
    The inverse of a diagonal matrix can also be applied.

    @param Op matrix to be scaled
    @param scalingVector vector that represents diagonal matrix
    @doInverse Indicates whether the inverse of the diagonal matrix should be applied.  (Default is to use inverse.)
    */
    static void ScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doInverse = true);


    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
    static Teuchos::Array<Magnitude> ResidualNorm(const Operator& Op, const MultiVector& X, const MultiVector& RHS);

    static RCP<MultiVector> Residual(const Operator& Op, const MultiVector& X, const MultiVector& RHS);

    static void PauseForDebugger();

    /*! @brief Simple transpose for Tpetra::CrsMatrix types

        Note:  This is very inefficient, as it inserts one entry at a time.
    */

    /*! @brief Power method.

    @param A matrix
    @param scaleByDiag if true, estimate the largest eigenvalue of \f$ D^; A \f$.
    @param niters maximum number of iterations
    @param tolerance stopping tolerance
    @verbose if true, print iteration information

    (Shamelessly grabbed from tpetra/examples.)
    */
    static Scalar PowerMethod(const Matrix& A, bool scaleByDiag = true,
                              LO niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123);

    static void MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const SC>& scalingVector, bool doInverse = true,
                                 bool doFillComplete = true, bool doOptimizeStorage = true);

    static void MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                                        bool doFillComplete, bool doOptimizeStorage);

    static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream& os);

    /*! @brief Squared distance between two rows in a multivector

       Used for coordinate vectors.
    */
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType Distance2(const MultiVector& v, LocalOrdinal i0, LocalOrdinal i1);

    /*! @brief Detect Dirichlet rows

        @param[in] A matrix
        @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.

        @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
    */
    static Teuchos::ArrayRCP<const bool> DetectDirichletRows(const Matrix& A, const Magnitude& tol = Teuchos::ScalarTraits<SC>::zero());

    /*! @brief Set seed for random number generator.

      Distribute the seeds evenly in [1,INT_MAX-1].  This guarantees nothing
      about where random number streams on difference processes will intersect.
      This does avoid overflow situations in parallel when multiplying by a PID.
      It also avoids the pathological case of having the *same* random number stream
      on each process.
    */

    static void SetRandomSeed(const Teuchos::Comm<int> &comm);

  }; // class Utils

  ///////////////////////////////////////////

#ifndef HAVE_MUELU_TPETRA_INST_INT_INT
  /*!
    @class Utils
    @brief MueLu utility class (specialization LO=GO=int).

    This class provides a number of static helper methods. Some are temporary and will eventually
    go away, while others should be moved to Xpetra.
  */
  //template <class Scalar, class Node>
  //class Utils<Scalar, int, int, Node> {
  template <class Node>
  class Utils<double,int,int,Node> {
  public:
    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;

#ifdef HAVE_MUELU_EPETRA
    //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
    // @{
    static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > const Vec) {
      RCP<const Xpetra::EpetraMultiVector > tmpVec = rcp_dynamic_cast<Xpetra::EpetraMultiVector>(Vec);
        if (tmpVec == Teuchos::null)
          throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
        return tmpVec->getEpetra_MultiVector();
    }
    static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Vec) {
      RCP<const Xpetra::EpetraMultiVector> tmpVec = rcp_dynamic_cast<Xpetra::EpetraMultiVector>(Vec);
      if (tmpVec == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
      return tmpVec->getEpetra_MultiVector();
    }

    static const Epetra_MultiVector&                        MV2EpetraMV(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Vec) {
      const Xpetra::EpetraMultiVector& tmpVec = dynamic_cast<const Xpetra::EpetraMultiVector&>(Vec);
      return *(tmpVec.getEpetra_MultiVector());
    }
    static       Epetra_MultiVector&                        MV2NonConstEpetraMV(Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Vec) {
      const Xpetra::EpetraMultiVector& tmpVec = dynamic_cast<const Xpetra::EpetraMultiVector&>(Vec);
      return *(tmpVec.getEpetra_MultiVector());
    }

    static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op) {
      RCP<const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(Op);
      if (crsOp == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
      const RCP<const Xpetra::EpetraCrsMatrix>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(crsOp->getCrsMatrix());
      if (tmp_ECrsMtx == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      return tmp_ECrsMtx->getEpetra_CrsMatrix();
    }
    static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op) {
      RCP<const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(Op);
      if (crsOp == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
      const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(crsOp->getCrsMatrix());
      if (tmp_ECrsMtx == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
    }

    static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Op) {
      try {
        const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>& crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(Op);
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
    static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Op) {
      try {
        Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(Op);
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

    static const Epetra_Map&                                Map2EpetraMap(const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & map) {
      RCP<const Xpetra::EpetraMap> xeMap = rcp_dynamic_cast<const Xpetra::EpetraMap>(rcpFromRef(map));
      if (xeMap == Teuchos::null)
        throw Exceptions::BadCast("Utils::Map2EpetraMap : Cast from Xpetra::Map to Xpetra::EpetraMap failed");
      return xeMap->getEpetra_Map();
    }
    // @}
#endif

#ifdef HAVE_MUELU_TPETRA
    //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
    // @{
    static RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     MV2TpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > const Vec)   { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::MultiVector<int,int>", "Xpetra::MultiVector<int,int>", "int"); return Teuchos::null; };
    static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     MV2NonConstTpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Vec) { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::MultiVector<int,int>", "Xpetra::MultiVector<int,int>", "int"); return Teuchos::null; };
    static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     MV2NonConstTpetraMV2(Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Vec)    { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::MultiVector<int,int>", "Xpetra::MultiVector<int,int>", "int"); return Teuchos::null; };

    static const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&          MV2TpetraMV(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Vec)   { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::MultiVector<int,int>", "Xpetra::MultiVector<int,int>", "int");  };
    static       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&          MV2NonConstTpetraMV(Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Vec) { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::MultiVector<int,int>", "Xpetra::MultiVector<int,int>", "int");  };

    static RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >       Op2TpetraCrs(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op)  { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::Matrix<int,int>", "Xpetra::Matrix<int,int>", "int"); return Teuchos::null; };
    static RCP<      Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >       Op2NonConstTpetraCrs(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op){ MUELU_TPETRA_ETI_EXCEPTION("Xpetra::Matrix<int,int>", "Xpetra::Matrix<int,int>", "int"); return Teuchos::null; };

    static const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&            Op2TpetraCrs(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Op)   { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::Matrix<int,int>", "Xpetra::Matrix<int,int>", "int");  };
    static       Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&            Op2NonConstTpetraCrs(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Op) { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::Matrix<int,int>", "Xpetra::Matrix<int,int>", "int");  };

    static RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >       Op2TpetraRow(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op)   { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::Matrix<int,int>", "Xpetra::Matrix<int,int>", "int"); return Teuchos::null; };
    static RCP<      Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >       Op2NonConstTpetraRow(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op) { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::Matrix<int,int>", "Xpetra::Matrix<int,int>", "int"); return Teuchos::null; };


    static const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >              Map2TpetraMap(const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & map) { MUELU_TPETRA_ETI_EXCEPTION("Xpetra::Map<int,int>", "Xpetra::Map<int,int>", "int"); return Teuchos::null; };
#endif

    static RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >                Crs2Op(RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op) {
      if (Op.is_null())
        return Teuchos::null;
      return rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Op));
    }

    /*! @brief Extract Matrix Diagonal

    Returns Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<Scalar> GetMatrixDiagonal(const Xpetra::Matrix<Scalar,int,int,Node>& A) {
      size_t numRows = A.getRowMap()->getNodeNumElements();
      Teuchos::ArrayRCP<Scalar> diag(numRows);
      Teuchos::ArrayView<const LocalOrdinal> cols;
      Teuchos::ArrayView<const Scalar> vals;
      for (size_t i = 0; i < numRows; ++i) {
        A.getLocalRowView(i, cols, vals);
        LocalOrdinal j = 0;
        for (; j < cols.size(); ++j) {
          if (Teuchos::as<size_t>(cols[j]) == i) {
            diag[i] = vals[j];
            break;
          }
        }
        if (j == cols.size()) {
          // Diagonal entry is absent
          diag[i] = Teuchos::ScalarTraits<Scalar>::zero();
        }
      }
      return diag;
    }

    /*! @brief Extract Matrix Diagonal

    Returns inverse of the Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetMatrixDiagonalInverse(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100) {
      RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap = A.getRowMap();
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > diag = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap);
      ArrayRCP<Scalar> diagVals = diag->getDataNonConst(0);
      size_t numRows = rowMap->getNodeNumElements();
      Teuchos::ArrayView<const LocalOrdinal> cols;
      Teuchos::ArrayView<const Scalar> vals;
      for (size_t i = 0; i < numRows; ++i) {
        A.getLocalRowView(i, cols, vals);
        LocalOrdinal j = 0;
        for (; j < cols.size(); ++j) {
          if (Teuchos::as<size_t>(cols[j]) == i) {
            if(Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]) > tol)
              diagVals[i] = Teuchos::ScalarTraits<Scalar>::one() / vals[j];
            else
              diagVals[i]=Teuchos::ScalarTraits<Scalar>::zero();
            break;
          }
        }
        if (j == cols.size()) {
          // Diagonal entry is absent
          diagVals[i]=Teuchos::ScalarTraits<Scalar>::zero();
        }
      }
      diagVals=null;
      return diag;
    }



    /*! @brief Extract Matrix Diagonal of lumped matrix

    Returns Matrix diagonal of lumped matrix in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<Scalar> GetLumpedMatrixDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) {
      size_t numRows = A.getRowMap()->getNodeNumElements();
      Teuchos::ArrayRCP<Scalar> diag(numRows);
      Teuchos::ArrayView<const LocalOrdinal> cols;
      Teuchos::ArrayView<const Scalar> vals;
      for (size_t i = 0; i < numRows; ++i) {
        A.getLocalRowView(i, cols, vals);
        diag[i] = Teuchos::ScalarTraits<Scalar>::zero();
        for (LocalOrdinal j = 0; j < cols.size(); ++j) {
          diag[i] += Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]);
        }
      }
      return diag;
    }

    /*! @brief Extract Overlapped Matrix Diagonal

    Returns overlapped Matrix diagonal in ArrayRCP.

    The local overlapped diagonal has an entry for each index in A's column map.
    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetMatrixOverlappedDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) {
      RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap = A.getRowMap(), colMap = A.getColMap();
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > localDiag = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap);

      try {
         const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>* crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>*>(&A);
         if (crsOp == NULL) {
           throw Exceptions::RuntimeError("cast to CrsMatrixWrap failed");
         }
         Teuchos::ArrayRCP<size_t> offsets;
         crsOp->getLocalDiagOffsets(offsets);
         crsOp->getLocalDiagCopy(*localDiag,offsets());
      }
      catch (...) {
        ArrayRCP<Scalar>   localDiagVals = localDiag->getDataNonConst(0);
        Teuchos::ArrayRCP<Scalar> diagVals = GetMatrixDiagonal(A);
        for (LocalOrdinal i = 0; i < localDiagVals.size(); i++)
          localDiagVals[i] = diagVals[i];
        localDiagVals = diagVals = null;
      }

      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > diagonal = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(colMap);
      RCP< const Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer;
      importer = A.getCrsGraph()->getImporter();
      if (importer == Teuchos::null) {
        importer = Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap, colMap);
      }
      diagonal->doImport(*localDiag, *(importer), Xpetra::INSERT);
      return diagonal;
    }

    /*! @brief Left scale matrix by an arbitrary vector.

    Algorithmically, this left scales a matrix by a diagonal matrix.
    The inverse of a diagonal matrix can also be applied.

    @param Op matrix to be scaled
    @param scalingVector vector that represents diagonal matrix
    @doInverse Indicates whether the inverse of the diagonal matrix should be applied.  (Default is to use inverse.)
    */
    static void ScaleMatrix(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector, bool doInverse = true) {
      // TODO this is incomplete.
    #ifdef HAVE_MUELU_TPETRA
          throw Exceptions::RuntimeError("This is the stub implementation of ScaleMatirx for Teptra with GO=int disabled. There is no Epetra implementation, too.");
    #else
        // TODO why is this not implemented in Epetra??
        throw Exceptions::RuntimeError("Matrix scaling has not been implemented Epetra");
    #endif // HAVE_MUELU_TPETRA
    }


    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
    static Teuchos::Array<Magnitude> ResidualNorm(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) {
      TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
       const size_t numVecs = X.getNumVectors();
       RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > RES = Residual(Op, X, RHS);
       Teuchos::Array<Magnitude> norms(numVecs);
       RES->norm2(norms);
       return norms;
    }

    static RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) {
      TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
        const size_t numVecs = X.getNumVectors();
        Scalar one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();
        RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > RES = Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Op.getRangeMap(), numVecs, false); // no need to initialize to zero
        Op.apply(X, *RES, Teuchos::NO_TRANS, one, zero);
        RES->update(one, RHS, negone);
        return RES;
    }

#ifndef _WIN32
#include <unistd.h>
    static void PauseForDebugger() {
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
    }
#else
    static void PauseForDebugger() {
         throw(Exceptions::RuntimeError("MueLu Utils: PauseForDebugger not implemented on Windows."));
     }
#endif

    /*! @brief Simple transpose for Tpetra::CrsMatrix types

        Note:  This is very inefficient, as it inserts one entry at a time.
    */

    /*! @brief Power method.

    @param A matrix
    @param scaleByDiag if true, estimate the largest eigenvalue of \f$ D^; A \f$.
    @param niters maximum number of iterations
    @param tolerance stopping tolerance
    @verbose if true, print iteration information

    (Shamelessly grabbed from tpetra/examples.)
    */
    static Scalar PowerMethod(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, bool scaleByDiag = true,
                              LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      TEUCHOS_TEST_FOR_EXCEPTION(!(A.getRangeMap()->isSameAs(*(A.getDomainMap()))), Exceptions::Incompatible,
          "Utils::PowerMethod: operator must have domain and range maps that are equivalent.");

      // Create three vectors, fill z with random numbers
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > q = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getDomainMap());
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > r = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getRangeMap());
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > z = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getRangeMap());

      z->setSeed(seed);  // seed random number generator
      z->randomize(true);// use Xpetra implementation: -> same results for Epetra and Tpetra

      Teuchos::Array<Magnitude> norms(1);

      typedef Teuchos::ScalarTraits<Scalar> STS;

      const Scalar zero = STS::zero(), one = STS::one();

      Scalar lambda = zero;
      Magnitude residual = STS::magnitude(zero);

      // power iteration
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > diagInvVec;
      if (scaleByDiag) {
        RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > diagVec = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getRowMap());
        A.getLocalDiagCopy(*diagVec);
        diagInvVec = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getRowMap());
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

    static void MyOldScaleMatrix(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Teuchos::ArrayRCP<const Scalar>& scalingVector, bool doInverse = true,
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

      /*switch (Op.getRowMap()->lib()) {
      // TODO fix me: put code directly in here!
        case Xpetra::UseEpetra:
          Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
          break;

        default:
          throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled.");
          break;
      }*/

      // no Tpetra version, since no Tpetra with GO=int enabled...
      // TODO what about scalars??
#ifdef HAVE_MUELU_EPETRA
    try {
      const Epetra_CrsMatrix& epOp = MueLu::Utils<double,int,int>::Op2NonConstEpetraCrs(Op);

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

    // This is the <int,int> specialization with TPETRA_INST_INT_INT disabled!
    static void MyOldScaleMatrix_Tpetra(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
                                        bool doFillComplete, bool doOptimizeStorage) {}

    static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream& os) {
      RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(os));
      return fancy;
    }

    /*! @brief Squared distance between two rows in a multivector

       Used for coordinate vectors.
    */
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType Distance2(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& v, LocalOrdinal i0, LocalOrdinal i1) {
      size_t numVectors = v.getNumVectors();

      Scalar d = Teuchos::ScalarTraits<Scalar>::zero();
      for (size_t j = 0; j < numVectors; j++) {
        Teuchos::ArrayRCP<const Scalar> vv = v.getData(j);
        d += (vv[i0] - vv[i1])*(vv[i0] - vv[i1]);
      }
      return Teuchos::ScalarTraits<Scalar>::magnitude(d);
    }

    /*! @brief Detect Dirichlet rows

        @param[in] A matrix
        @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.

        @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
    */
    static Teuchos::ArrayRCP<const bool> DetectDirichletRows(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::zero()) {
      LocalOrdinal numRows = A.getNodeNumRows();
      typedef Teuchos::ScalarTraits<Scalar> STS;
      ArrayRCP<bool> boundaryNodes(numRows, true);
      for (LocalOrdinal row = 0; row < numRows; row++) {
        ArrayView<const LocalOrdinal> indices;
        ArrayView<const Scalar> vals;
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

    /*! @brief Set seed for random number generator.

      Distribute the seeds evenly in [1,INT_MAX-1].  This guarantees nothing
      about where random number streams on difference processes will intersect.
      This does avoid overflow situations in parallel when multiplying by a PID.
      It also avoids the pathological case of having the *same* random number stream
      on each process.
    */

    static void SetRandomSeed(const Teuchos::Comm<int> &comm) {
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
      Teuchos::ScalarTraits<Scalar>::seedrandom(mySeed);
      // Epetra
      //   MultiVector::Random() -> Epetra_Util::RandomDouble() -> Epetra_Utils::RandomInt()
      // Its own random number generator, based on Seed_. Seed_ is initialized in Epetra_Util constructor with std::rand()
      // So our setting std::srand() affects that too
    }


  }; // class Utils (specialization LO=GO=int)
#endif



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
    RCP<Xpetra::CrsMatrix<double,int,int,KDNT> > tmpC2 = rcp_implicit_cast<Xpetra::CrsMatrix<double,int,int,KDNT> >(tmpC1);
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

  /*!
    @class Utils2
    @brief MueLu utility class.

    Separate class for utilities that need a specialization for Epetra.
  */
  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = KokkosClassic::DefaultNode::DefaultNodeType>
  class Utils2 {

#include "MueLu_UseShortNames.hpp"

  public:

    /*! @brief Transpose a Xpetra::Matrix

    Note: Currently, an error is thrown if the matrix isn't a Tpetra::CrsMatrix or Epetra_CrsMatrix.
    In principle, however, we could allow any Epetra_RowMatrix because the Epetra transposer does.
    */
    static RCP<Matrix> Transpose(Matrix& Op, bool optimizeTranspose = false,const std::string & label = std::string());

    //! Scale an Epetra matrix.
    static void MyOldScaleMatrix_Epetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doFillComplete, bool doOptimizeStorage);
  }; // class Utils2

  // specialization Utils2 for SC=double, LO=GO=int
  template<>
  class Utils2<double,int,int> {
    typedef double                                                     SC;
    typedef int                                                        LO;
    typedef int                                                        GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType                NO;
    typedef Xpetra::Map<int,int,NO>                                    Map;
    typedef Xpetra::Matrix<double,int,int,NO>                          Matrix;
    typedef Xpetra::MultiVector<double,int,int,NO>                     MultiVector;

  public:

    static RCP<Matrix>      Transpose               (Matrix& Op, bool optimizeTranspose = false,const std::string & label = std::string());
    static void             MyOldScaleMatrix_Epetra (Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doFillComplete, bool doOptimizeStorage);
  };


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
//#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
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
//#endif
#endif

} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DECL_HPP
