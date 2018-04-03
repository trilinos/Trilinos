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
    static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > const vec);
    static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > vec);

    static const Epetra_MultiVector&                        MV2EpetraMV(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& vec);
    static       Epetra_MultiVector&                        MV2NonConstEpetraMV(Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& vec);

    static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);
    static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);

    static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);
    static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);

    static const Epetra_Map&                                Map2EpetraMap(const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& map);
    // @}
#endif

#ifdef HAVE_MUELU_TPETRA
    //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
    static RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2TpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > const vec);
    static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2NonConstTpetraMV(RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > vec);
    static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2NonConstTpetraMV2(Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& vec);

    static const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&      MV2TpetraMV(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& vec);
    static       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&      MV2NonConstTpetraMV(Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& vec);

    static RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraCrs(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);
    static RCP<      Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraCrs(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);

    static const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2TpetraCrs(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);
    static       Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2NonConstTpetraCrs(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);

    static RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraRow(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);
    static RCP<      Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraRow(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);


    static const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >        Map2TpetraMap(const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& map);
#endif

    static RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >          Crs2Op(RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op) { return UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(Op); }
    static Teuchos::ArrayRCP<Scalar>                                             GetMatrixDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonal(A); }
    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >          GetMatrixDiagonalInverse(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonalInverse(A,tol); }
    static Teuchos::ArrayRCP<Scalar>                                             GetLumpedMatrixDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetLumpedMatrixDiagonal(A); }
    static Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetLumpedMatrixDiagonal(Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetLumpedMatrixDiagonal(A); }
    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >          GetMatrixOverlappedDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixOverlappedDiagonal(A); }
    static Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetInverse(Teuchos::RCP<const Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > v, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100, Scalar tolReplacement = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetInverse(v,tol,tolReplacement); }
    static Teuchos::Array<Magnitude>                                             ResidualNorm(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ResidualNorm(Op,X,RHS); }
    static RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Residual(Op,X,RHS); }
    static void Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op,  const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,  const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Resid) { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Residual(Op,X,RHS,Resid);}
    static void                                                                  PauseForDebugger() { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PauseForDebugger(); }
    static RCP<Teuchos::FancyOStream>                                            MakeFancy(std::ostream& os) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MakeFancy(os); }
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType                 Distance2(const Teuchos::Array<Teuchos::ArrayRCP<const Scalar>>& v, LocalOrdinal i0, LocalOrdinal i1) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Distance2(v,i0,i1); }
    static Teuchos::ArrayRCP<const bool>                                         DetectDirichletRows(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::magnitude(0.), const bool count_twos_as_dirichlet=false) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletRows(A,tol,count_twos_as_dirichlet); }
    static Teuchos::ArrayRCP<const bool>                                         DetectDirichletRowsExt(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, bool & bHasZeroDiagonal, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletRowsExt(A,bHasZeroDiagonal,tol); }

    static void                                                                  SetRandomSeed(const Teuchos::Comm<int> &comm) { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetRandomSeed(comm); }

    static Scalar PowerMethod(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, bool scaleByDiag = true,
                              LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PowerMethod(A,scaleByDiag,niters,tolerance,verbose,seed);
    }

    static Scalar Frobenius(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B) {
      return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Frobenius(A, B);
    }

    static void MyOldScaleMatrix(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Teuchos::ArrayRCP<const Scalar>& scalingVector, bool doInverse = true,
                                 bool doFillComplete = true, bool doOptimizeStorage = true);

    static void MyOldScaleMatrix_Epetra(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
                                        bool doFillComplete, bool doOptimizeStorage);
    static void MyOldScaleMatrix_Tpetra(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
                                        bool doFillComplete, bool doOptimizeStorage);

    static RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Transpose(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, bool optimizeTranspose = false,const std::string & label = std::string(),const Teuchos::RCP<Teuchos::ParameterList> &params=Teuchos::null);

    static RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > RealValuedToScalarMultiVector(RCP<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> > X);

    static RCP<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> > ExtractCoordinatesFromParameterList(ParameterList& paramList);

    static void FindDirichletRows(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & A, std::vector<LocalOrdinal>& dirichletRows,bool count_twos_as_dirichlet=false) {
      MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::FindDirichletRows(A,dirichletRows,count_twos_as_dirichlet);
    }


    static void ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,const std::vector<LocalOrdinal>& dirichletRows) {
      MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyOAZToMatrixRows(A,dirichletRows);
    }

    static void ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,const Teuchos::ArrayRCP<const bool>& dirichletRows) {
      MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyOAZToMatrixRows(A,dirichletRows);
    }

    static void ZeroDirichletRows(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,const std::vector<LocalOrdinal>& dirichletRows, Scalar replaceWith=Teuchos::ScalarTraits<Scalar>::zero()) {
      MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ZeroDirichletRows(A,dirichletRows,replaceWith);
    }

    static void ZeroDirichletRows(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,const Teuchos::ArrayRCP<const bool>& dirichletRows, Scalar replaceWith=Teuchos::ScalarTraits<Scalar>::zero()) {
      MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ZeroDirichletRows(A,dirichletRows,replaceWith);
    }

    static void ZeroDirichletRows(Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& X,const Teuchos::ArrayRCP<const bool>& dirichletRows,Scalar replaceWith=Teuchos::ScalarTraits<Scalar>::zero()) {
      MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ZeroDirichletRows(X,dirichletRows,replaceWith);
    }

    static void ZeroDirichletCols(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,const Teuchos::ArrayRCP<const bool>& dirichletCols, Scalar replaceWith=Teuchos::ScalarTraits<Scalar>::zero()) {
      MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ZeroDirichletCols(A,dirichletCols,replaceWith);
    }

  }; // class Utilities

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
  class Utilities<double,int,int,Xpetra::EpetraNode> : public UtilitiesBase<double,int,int,Xpetra::EpetraNode> {
  public:
    typedef double              Scalar;
    typedef int                 LocalOrdinal;
    typedef int                 GlobalOrdinal;
    typedef Xpetra::EpetraNode  Node;
    typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;

  private:
    typedef Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsMatrixWrap;
    typedef Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CrsMatrix;
    typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Matrix;
    typedef Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Vector;
    typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MultiVector;
    typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
#ifdef HAVE_MUELU_EPETRA
    typedef Xpetra::EpetraMapT<GlobalOrdinal,Node> EpetraMap;
    typedef Xpetra::EpetraMultiVectorT<GlobalOrdinal,Node> EpetraMultiVector;
    typedef Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node> EpetraCrsMatrix;
#endif
  public:

#ifdef HAVE_MUELU_EPETRA
    //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
    // @{
    static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<MultiVector> const vec) {
      RCP<const EpetraMultiVector > tmpVec = rcp_dynamic_cast<EpetraMultiVector>(vec);
      if (tmpVec == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
      return tmpVec->getEpetra_MultiVector();
    }
    static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<MultiVector> vec) {
      RCP<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(vec);
      if (tmpVec == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
      return tmpVec->getEpetra_MultiVector();
    }

    static const Epetra_MultiVector&                        MV2EpetraMV(const MultiVector& vec) {
      const EpetraMultiVector& tmpVec = dynamic_cast<const EpetraMultiVector&>(vec);
      return *(tmpVec.getEpetra_MultiVector());
    }
    static       Epetra_MultiVector&                        MV2NonConstEpetraMV(MultiVector& vec) {
      const EpetraMultiVector& tmpVec = dynamic_cast<const EpetraMultiVector&>(vec);
      return *(tmpVec.getEpetra_MultiVector());
    }

    static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Matrix> Op) {
      RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
      if (crsOp == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
      const RCP<const EpetraCrsMatrix>& tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(crsOp->getCrsMatrix());
      if (tmp_ECrsMtx == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      return tmp_ECrsMtx->getEpetra_CrsMatrix();
    }
    static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Matrix> Op) {
      RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
      if (crsOp == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
      const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(crsOp->getCrsMatrix());
      if (tmp_ECrsMtx == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
    }

    static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Matrix& Op) {
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
    static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Matrix& Op) {
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

    static const Epetra_Map&                                Map2EpetraMap(const Map& map) {
      RCP<const EpetraMap> xeMap = rcp_dynamic_cast<const EpetraMap>(rcpFromRef(map));
      if (xeMap == Teuchos::null)
        throw Exceptions::BadCast("Utilities::Map2EpetraMap : Cast from Xpetra::Map to Xpetra::EpetraMap failed");
      return xeMap->getEpetra_Map();
    }
    // @}
#endif

#ifdef HAVE_MUELU_TPETRA
    //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
    static RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2TpetraMV(RCP<MultiVector> const vec)   {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("MV2TpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#else
      RCP<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmpVec = rcp_dynamic_cast<Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(vec);
      if (tmpVec == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
      return tmpVec->getTpetra_MultiVector();
#endif
    }
    static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2NonConstTpetraMV(RCP<MultiVector> vec) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("MV2NonConstTpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#else
      RCP<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmpVec = rcp_dynamic_cast<Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(vec);
      if (tmpVec == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
      return tmpVec->getTpetra_MultiVector();
#endif

    }
    static RCP<      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MV2NonConstTpetraMV2(MultiVector& vec)    {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("MV2NonConstTpetraMV2: Tpetra has not been compiled with support for LO=GO=int.");
#else
      const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(vec);
      return tmpVec.getTpetra_MultiVector();
#endif
    }

    static const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&      MV2TpetraMV(const MultiVector& vec)   {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("MV2TpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#else
      const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(vec);
      return *(tmpVec.getTpetra_MultiVector());
#endif
    }
    static       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&      MV2NonConstTpetraMV(MultiVector& vec) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("MV2NonConstTpetraMV: Tpetra has not been compiled with support for LO=GO=int.");
#else
      const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(vec);
      return *(tmpVec.getTpetra_MultiVector());
#endif
    }

    static RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraCrs(RCP<const Matrix> Op)  {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Op2TpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
      // Get the underlying Tpetra Mtx
      RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
      if (crsOp == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
      const RCP<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(crsOp->getCrsMatrix());
      if (tmp_ECrsMtx == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      return tmp_ECrsMtx->getTpetra_CrsMatrix();
#endif
    }
    static RCP<      Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraCrs(RCP<Matrix> Op){
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Op2NonConstTpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
      RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
      if (crsOp == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
      const RCP<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(crsOp->getCrsMatrix());
      if (tmp_ECrsMtx == Teuchos::null)
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      return tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
#endif
    };

    static const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2TpetraCrs(const Matrix& Op)   {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Op2TpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
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
#endif
    }
    static       Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2NonConstTpetraCrs(Matrix& Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Op2NonConstTpetraCrs: Tpetra has not been compiled with support for LO=GO=int.");
#else
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
#endif
    }

    static RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraRow(RCP<const Matrix> Op)   {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Op2TpetraRow: Tpetra has not been compiled with support for LO=GO=int.");
#else
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
#endif
    }
    static RCP<      Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraRow(RCP<Matrix> Op) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Op2NonConstTpetraRow: Tpetra has not been compiled with support for LO=GO=int.");
#else
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
#endif
    };


    static const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >          Map2TpetraMap(const Map& map) {
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Map2TpetraMap: Tpetra has not been compiled with support for LO=GO=int.");
#else
      const RCP<const Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>>& tmp_TMap = rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> >(rcpFromRef(map));
      if (tmp_TMap == Teuchos::null)
        throw Exceptions::BadCast("Utilities::Map2TpetraMap : Cast from Xpetra::Map to Xpetra::TpetraMap failed");
      return tmp_TMap->getTpetra_Map();
#endif
    };
#endif

    static RCP<Matrix>                                                           Crs2Op(RCP<CrsMatrix> Op) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(Op); }
    static Teuchos::ArrayRCP<Scalar>                                             GetMatrixDiagonal(const Matrix& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonal(A); }
    static RCP<Vector>                                                           GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonalInverse(A,tol); }
    static Teuchos::ArrayRCP<Scalar>                                             GetLumpedMatrixDiagonal(const Matrix& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetLumpedMatrixDiagonal(A); }
    static Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetLumpedMatrixDiagonal(Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetLumpedMatrixDiagonal(A); }
    static RCP<Vector>                                                           GetMatrixOverlappedDiagonal(const Matrix& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixOverlappedDiagonal(A); }
    static RCP<Vector>                                                           GetInverse(Teuchos::RCP<const Vector> v, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100, Scalar tolReplacement = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetInverse(v,tol,tolReplacement); }
    static Teuchos::Array<Magnitude>                                             ResidualNorm(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ResidualNorm(Op,X,RHS); }
    static RCP<MultiVector>                                                      Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Residual(Op,X,RHS); }
    static void Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op,  const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,  const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Resid) { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Residual(Op,X,RHS,Resid);}
    static void                                                                  PauseForDebugger() { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PauseForDebugger(); }
    static RCP<Teuchos::FancyOStream>                                            MakeFancy(std::ostream& os) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MakeFancy(os); }
    static Teuchos::ScalarTraits<Scalar>::magnitudeType                 Distance2(const Teuchos::Array<Teuchos::ArrayRCP<const Scalar>>& v, LocalOrdinal i0, LocalOrdinal i1) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Distance2(v,i0,i1); }
    static Teuchos::ArrayRCP<const bool>                                         DetectDirichletRows(const Matrix& A, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::zero(), const bool count_twos_as_dirichlet=false) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletRows(A,tol,count_twos_as_dirichlet); }
    static Teuchos::ArrayRCP<const bool>                                         DetectDirichletRowsExt(const Matrix& A, bool & bHasZeroDiagonal, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletRowsExt(A,bHasZeroDiagonal,tol); }
    static void                                                                  SetRandomSeed(const Teuchos::Comm<int> &comm) { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetRandomSeed(comm); }

    static Scalar PowerMethod(const Matrix& A, bool scaleByDiag = true,
                              LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PowerMethod(A,scaleByDiag,niters,tolerance,verbose,seed);
    }

    static Scalar Frobenius(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B) {
      return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Frobenius(A, B);
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
      }
    }

    // TODO This is the <double,int,int> specialization
    static void MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
                                        bool doFillComplete, bool doOptimizeStorage) {
#ifdef HAVE_MUELU_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Matrix scaling is not possible because Tpetra has not been compiled with support for LO=GO=int.");
#else
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
    static RCP<Matrix> Transpose(Matrix& Op, bool optimizeTranspose = false,const std::string & label = std::string(),const Teuchos::RCP<Teuchos::ParameterList> &params=Teuchos::null) {
      switch (Op.getRowMap()->lib()) {
        case Xpetra::UseTpetra: {
#ifdef HAVE_MUELU_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
            throw Exceptions::RuntimeError("Utilities::Transpose: Tpetra is not compiled with LO=GO=int. Add TPETRA_INST_INT_INT:BOOL=ON to your configuration!");
#else
            try {
              const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpetraOp = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(Op);

              // Compute the transpose A of the Tpetra matrix tpetraOp.
              RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
              Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp),label);
              A = transposer.createTranspose(params);
              RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AA   = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A));
              RCP<CrsMatrix>                                                           AAA  = rcp_implicit_cast<CrsMatrix>(AA);
              RCP<Matrix>                                                              AAAA = rcp( new CrsMatrixWrap(AAA));

              if (Op.IsView("stridedMaps"))
                AAAA->CreateView("stridedMaps", Teuchos::rcpFromRef(Op), true/*doTranspose*/);

              return AAAA;
            }
            catch (std::exception& e) {
              std::cout << "threw exception '" << e.what() << "'" << std::endl;
              throw Exceptions::RuntimeError("Utilities::Transpose failed, perhaps because matrix is not a Crs matrix");
            }
#endif
#else
            throw Exceptions::RuntimeError("Utilities::Transpose: Tpetra is not compiled!");
#endif
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
            RCP<EpetraCrsMatrix> AA   = rcp(new EpetraCrsMatrix(rcpA));
            RCP<CrsMatrix>       AAA  = rcp_implicit_cast<CrsMatrix>(AA);
            RCP<Matrix>          AAAA = rcp( new CrsMatrixWrap(AAA));

            if (Op.IsView("stridedMaps"))
              AAAA->CreateView("stridedMaps", Teuchos::rcpFromRef(Op), true/*doTranspose*/);

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

    static RCP<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> >
    RealValuedToScalarMultiVector(RCP<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> > X) {
      RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Xscalar = rcp_dynamic_cast<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(X);
      return Xscalar;
    }

    /*! @brief Extract coordinates from parameter list and return them in a Xpetra::MultiVector
    */
    static RCP<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> > ExtractCoordinatesFromParameterList(ParameterList& paramList) {
      RCP<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> > coordinates = Teuchos::null;

      // check whether coordinates are contained in parameter list
      if(paramList.isParameter ("Coordinates") == false)
        return coordinates;

  #if defined(HAVE_MUELU_TPETRA)
  #if ( defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT)) || \
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
      typedef Tpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> tdMV;
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
      if(doubleCoords != Teuchos::null) {
        coordinates = Teuchos::rcp(new Xpetra::TpetraMultiVector<double, LocalOrdinal, GlobalOrdinal, Node>(doubleCoords));
        TEUCHOS_TEST_FOR_EXCEPT(doubleCoords->getNumVectors() != coordinates->getNumVectors());
      }
  #endif // Tpetra instantiated on GO=int and EpetraNode
  #endif // endif HAVE_TPETRA

  #if defined(HAVE_MUELU_EPETRA)
      RCP<Epetra_MultiVector> doubleEpCoords;
      if (paramList.isType<RCP<Epetra_MultiVector> >("Coordinates")) {
        doubleEpCoords = paramList.get<RCP<Epetra_MultiVector> >("Coordinates");
        paramList.remove("Coordinates");
        RCP<Xpetra::EpetraMultiVectorT<GlobalOrdinal,Node> > epCoordinates = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GlobalOrdinal,Node>(doubleEpCoords));
        coordinates = rcp_dynamic_cast<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> >(epCoordinates);
        TEUCHOS_TEST_FOR_EXCEPT(doubleEpCoords->NumVectors() != Teuchos::as<int>(coordinates->getNumVectors()));
      }
  #endif

      // check for Xpetra coordinates vector
      if(paramList.isType<decltype(coordinates)>("Coordinates")) {
        coordinates = paramList.get<decltype(coordinates)>("Coordinates");
      }

      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(coordinates));
      return coordinates;
    }

  }; // class Utilities (specialization SC=double LO=GO=int)

#endif // HAVE_MUELU_EPETRA



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
    typedef Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>                      XECrsMatrix;
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
    return rcp(new Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>(V));
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
