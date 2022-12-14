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

#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_CrsGraphFactory_fwd.hpp>
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

#include "MueLu_Exceptions.hpp"

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

  template<typename SC,typename LO,typename GO,typename NO>
  RCP<Xpetra::Matrix<SC, LO, GO, NO> >
  TpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> >& Atpetra);

  template<typename SC,typename LO,typename GO,typename NO>
  RCP<Xpetra::Matrix<SC, LO, GO, NO> >
  TpetraFECrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::FECrsMatrix<SC, LO, GO, NO> >& Atpetra);

  template<typename SC,typename LO,typename GO,typename NO>
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> >
  TpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> >& Vtpetra);

  template<typename SC,typename LO,typename GO,typename NO>
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> >
  TpetraFEMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Tpetra::FEMultiVector<SC, LO, GO, NO> >& Vtpetra);

  template<typename SC,typename LO,typename GO,typename NO>
  void leftRghtDofScalingWithinNode(const Xpetra::Matrix<SC,LO,GO,NO> & Atpetra, size_t blkSize, size_t nSweeps, Teuchos::ArrayRCP<SC> & rowScaling, Teuchos::ArrayRCP<SC> & colScaling);

  /*!
    @class Utilities
    @brief MueLu utility class.

    This class provides a number of static helper methods. Some are temporary and will eventually
    go away, while others should be moved to Xpetra.
    */
  template <class Scalar,
            class LocalOrdinal = DefaultLocalOrdinal,
            class GlobalOrdinal = DefaultGlobalOrdinal,
            class Node = DefaultNode>
  class Utilities : public UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_UTILITIES_SHORT
    //#include "MueLu_UseShortNames.hpp"

  public:
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;

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

    static RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraBlockCrs(RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);
   static RCP<      Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraBlockCrs(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);

   static const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2TpetraBlockCrs(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);
    static       Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&        Op2NonConstTpetraBlockCrs(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op);



    static RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2TpetraRow(RCP<const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);
    static RCP<      Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   Op2NonConstTpetraRow(RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op);


    static const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >        Map2TpetraMap(const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& map);

    static RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >          Crs2Op(RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op) { return UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(Op); }
    static RCP<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> >   GetThresholdedMatrix(const RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A, const Scalar threshold, const bool keepDiagonal, const GlobalOrdinal expectedNNZperRow) {return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetThresholdedMatrix(A, threshold, keepDiagonal, expectedNNZperRow); }
    static RCP<Xpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >               GetThresholdedGraph(const RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A, const Magnitude threshold, const GlobalOrdinal expectedNNZperRow) {return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetThresholdedGraph(A, threshold, expectedNNZperRow); }
    static Teuchos::ArrayRCP<Scalar>                                             GetMatrixDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonal(A); }
    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >          GetMatrixDiagonalInverse(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100, Scalar valReplacement = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixDiagonalInverse(A,tol,valReplacement); }
    static Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetLumpedMatrixDiagonal(Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> const &A, const bool doReciprocal=false, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100, Scalar tolReplacement = Teuchos::ScalarTraits<Scalar>::zero(), const bool replaceSingleEntryRowWithZero = false, const bool useAverageAbsDiagVal = false)
    { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetLumpedMatrixDiagonal(A, doReciprocal, tol, tolReplacement, replaceSingleEntryRowWithZero, useAverageAbsDiagVal); }
    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >          GetMatrixOverlappedDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixOverlappedDiagonal(A); }
    static Teuchos::ArrayRCP<Magnitude>                                          GetMatrixMaxMinusOffDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixMaxMinusOffDiagonal(A); }
    static Teuchos::ArrayRCP<Magnitude>                                          GetMatrixMaxMinusOffDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> &BlockNumber) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetMatrixMaxMinusOffDiagonal(A,BlockNumber); }

    static Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetInverse(Teuchos::RCP<const Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > v, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100, Scalar tolReplacement = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetInverse(v,tol,tolReplacement); }
    static Teuchos::Array<Magnitude>                                             ResidualNorm(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ResidualNorm(Op,X,RHS); }
    static Teuchos::Array<Magnitude>                                             ResidualNorm(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Resid) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ResidualNorm(Op,X,RHS,Resid); }
    static RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Residual(Op,X,RHS); }
    static void Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op,  const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,  const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS, Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Resid) { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Residual(Op,X,RHS,Resid);}
    static RCP<Teuchos::FancyOStream>                                            MakeFancy(std::ostream& os) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MakeFancy(os); }
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType                 Distance2(const Teuchos::Array<Teuchos::ArrayRCP<const Scalar>>& v, LocalOrdinal i0, LocalOrdinal i1) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Distance2(v,i0,i1); }
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType                 Distance2(const Teuchos::ArrayView<double> & weight,const Teuchos::Array<Teuchos::ArrayRCP<const Scalar>>& v, LocalOrdinal i0, LocalOrdinal i1) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Distance2(weight,v,i0,i1); }
    static Teuchos::ArrayRCP<const bool>                                         DetectDirichletRows(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::magnitude(0.), const bool count_twos_as_dirichlet=false) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletRows(A,tol,count_twos_as_dirichlet); }
    static Teuchos::ArrayRCP<const bool>                                         DetectDirichletRowsExt(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, bool & bHasZeroDiagonal, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletRowsExt(A,bHasZeroDiagonal,tol); }
    static void                                                                  FindNonZeros(const Teuchos::ArrayRCP<const Scalar> vals,Teuchos::ArrayRCP<bool> nonzeros) {return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::FindNonZeros(vals,nonzeros);}
    static void                                                                  DetectDirichletColsAndDomains(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Teuchos::ArrayRCP<bool>& dirichletRows,  Teuchos::ArrayRCP<bool> dirichletCols, Teuchos::ArrayRCP<bool> dirichletDomains) {MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletColsAndDomains(A,dirichletRows,dirichletCols,dirichletDomains); }

    static void                                                                  ApplyRowSumCriterion(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Magnitude rowSumTol, Teuchos::ArrayRCP<bool>& dirichletRows) {return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyRowSumCriterion(A,rowSumTol,dirichletRows); }
    static void                                                                  ApplyRowSumCriterion(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> &BlockNumber,const Magnitude rowSumTol, Teuchos::ArrayRCP<bool>& dirichletRows) {return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyRowSumCriterion(A,BlockNumber,rowSumTol,dirichletRows); }

    static void                                                                  SetRandomSeed(const Teuchos::Comm<int> &comm) { MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetRandomSeed(comm); }

    static Scalar PowerMethod(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, bool scaleByDiag = true,
                              LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PowerMethod(A,scaleByDiag,niters,tolerance,verbose,seed);
    }

    static Scalar PowerMethod(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> &invDiag,
                              LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PowerMethod(A,invDiag,niters,tolerance,verbose,seed);
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

    static RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > RealValuedToScalarMultiVector(RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node> > X);

    static RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LocalOrdinal,GlobalOrdinal,Node> > ExtractCoordinatesFromParameterList(ParameterList& paramList);

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
  void TokenizeStringAndStripWhiteSpace(const std::string & stream, std::vector<std::string> & tokenList, const char* token = ",");

  /*! Returns true if a parameter name is a valid Muemex custom level variable, e.g. "MultiVector myArray"
  */
  bool IsParamMuemexVariable(const std::string& name);

  /*! Returns true if a parameter name is a valid user custom level variable, e.g. "MultiVector myArray"
  */
  bool IsParamValidVariable(const std::string& name);

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
  void leftRghtDofScalingWithinNode(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Amat, size_t blkSize, size_t nSweeps, Teuchos::ArrayRCP<Scalar> & rowScaling, Teuchos::ArrayRCP<Scalar> & colScaling) {

     LocalOrdinal     nBlks = (Amat.getRowMap()->getLocalNumElements())/blkSize;

     Teuchos::ArrayRCP<Scalar>   rowScaleUpdate(blkSize);
     Teuchos::ArrayRCP<Scalar>   colScaleUpdate(blkSize);


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
             size_t modGuy = (cols[kk]+1)%blkSize;
             if (modGuy == 0) modGuy = blkSize;
             modGuy--;
             rowScaleUpdate[j] += rowScaling[j]*(Teuchos::ScalarTraits<Scalar>::magnitude(vals[kk]))*colScaling[modGuy];
           }
           row++;
         }
       }
       // combine information across processors
       Teuchos::ArrayRCP<Scalar>   tempUpdate(blkSize);
       Teuchos::reduceAll(*(Amat.getRowMap()->getComm()), Teuchos::REDUCE_SUM, (LocalOrdinal) blkSize, rowScaleUpdate.getRawPtr(), tempUpdate.getRawPtr());
       for (size_t i = 0; i < blkSize; i++) rowScaleUpdate[i] = tempUpdate[i];

       /* We want to scale by sqrt(1/rowScaleUpdate), but we'll         */
       /* normalize things by the minimum rowScaleUpdate. That is, the  */
       /* largest scaling is always one (as normalization is arbitrary).*/

       Scalar minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude((rowScaleUpdate[0]/rowScaling[0])/rowScaling[0]);

       for (size_t i = 1; i < blkSize; i++) {
          Scalar  temp = (rowScaleUpdate[i]/rowScaling[i])/rowScaling[i];
          if ( Teuchos::ScalarTraits<Scalar>::magnitude(temp) < Teuchos::ScalarTraits<Scalar>::magnitude(minUpdate))
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
             size_t modGuy = (cols[kk]+1)%blkSize;
             if (modGuy == 0) modGuy = blkSize;
             modGuy--;
             colScaleUpdate[modGuy] += colScaling[modGuy]* (Teuchos::ScalarTraits<Scalar>::magnitude(vals[kk])) *rowScaling[j];
           }
           row++;
         }
       }
       Teuchos::reduceAll(*(Amat.getRowMap()->getComm()), Teuchos::REDUCE_SUM, (LocalOrdinal) blkSize, colScaleUpdate.getRawPtr(), tempUpdate.getRawPtr());
       for (size_t i = 0; i < blkSize; i++) colScaleUpdate[i] = tempUpdate[i];

       /* We want to scale by sqrt(1/colScaleUpdate), but we'll         */
       /* normalize things by the minimum colScaleUpdate. That is, the  */
       /* largest scaling is always one (as normalization is arbitrary).*/


       minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude((colScaleUpdate[0]/colScaling[0])/colScaling[0]);

       for (size_t i = 1; i < blkSize; i++) {
          Scalar  temp = (colScaleUpdate[i]/colScaling[i])/colScaling[i];
          if ( Teuchos::ScalarTraits<Scalar>::magnitude(temp) < Teuchos::ScalarTraits<Scalar>::magnitude(minUpdate))
            minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude(temp);
       }
       for (size_t i = 0; i < blkSize; i++) colScaling[i] *= sqrt(minUpdate/colScaleUpdate[i]);
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
    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>       XCrsMatrix;
    typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>   XCrsMatrixWrap;

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
  template<class T>
  std::string toString(const T& what) {
    std::ostringstream buf;
    buf << what;
    return buf.str();
  }

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
  Teuchos::RCP<const Teuchos::Comm<int> > GenerateNodeComm(RCP<const Teuchos::Comm<int> > & baseComm, int &NodeId, const int reductionFactor);

  // Lower case string
  std::string lowerCase (const std::string& s);

} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DECL_HPP
