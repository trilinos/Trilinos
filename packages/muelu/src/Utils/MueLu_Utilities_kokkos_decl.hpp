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
#ifndef MUELU_UTILITIES_KOKKOS_DECL_HPP
#define MUELU_UTILITIES_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <string>

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Xpetra_BlockedCrsMatrix_fwd.hpp"
#include "Xpetra_CrsMatrix_fwd.hpp"
#include "Xpetra_CrsMatrixWrap_fwd.hpp"
#include "Xpetra_ExportFactory.hpp"
#include "Xpetra_ImportFactory_fwd.hpp"
#include "Xpetra_MapFactory_fwd.hpp"
#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_MatrixFactory_fwd.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_MultiVector_fwd.hpp"
#include "Xpetra_Operator_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"
#include "Xpetra_Vector_fwd.hpp"

#include "Xpetra_IO.hpp"

#include "Kokkos_ArithTraits.hpp"

#ifdef HAVE_MUELU_EPETRA
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Xpetra_EpetraCrsMatrix_fwd.hpp"
#include "Xpetra_EpetraMultiVector_fwd.hpp"
#endif

#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_UtilitiesBase.hpp"

#ifdef HAVE_MUELU_TPETRA
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Xpetra_TpetraCrsMatrix_fwd.hpp"
#include "Xpetra_TpetraMultiVector_fwd.hpp"
#endif


namespace MueLu {

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
  class Utilities_kokkos : public MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
#undef MUELU_UTILITIES_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    using TST                   = Teuchos::ScalarTraits<SC>;
    using Magnitude             = typename TST::magnitudeType;
    using CoordinateType        = typename TST::coordinateType;
    using RealValuedMultiVector = Xpetra::MultiVector<CoordinateType,LO,GO,NO>;

#ifdef HAVE_MUELU_EPETRA
    //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
    // @{
    static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<MultiVector> const vec)     { return Utilities::MV2EpetraMV(vec); }
    static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<MultiVector> vec)   { return Utilities::MV2NonConstEpetraMV(vec); }

    static const Epetra_MultiVector&                        MV2EpetraMV(const MultiVector& vec)         { return Utilities::MV2EpetraMV(vec); }
    static       Epetra_MultiVector&                        MV2NonConstEpetraMV(MultiVector& vec)       { return Utilities::MV2NonConstEpetraMV(vec); }

    static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Matrix> Op)          { return Utilities::Op2EpetraCrs(Op); }
    static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Matrix> Op)        { return Utilities::Op2NonConstEpetraCrs(Op); }

    static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Matrix& Op)              { return Utilities::Op2EpetraCrs(Op); }
    static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Matrix& Op)            { return Utilities::Op2NonConstEpetraCrs(Op); }

    static const Epetra_Map&                                Map2EpetraMap(const Map& map)               { return Utilities::Map2EpetraMap(map); }
    // @}
#endif

#ifdef HAVE_MUELU_TPETRA
    //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
    // @{
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> >     MV2TpetraMV(RCP<MultiVector> const vec)     { return Utilities::MV2TpetraMV(vec); }
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV(RCP<MultiVector> vec)   { return Utilities::MV2NonConstTpetraMV(vec); }
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV2(MultiVector& vec)      { return Utilities::MV2NonConstTpetraMV2(vec); }

    static const Tpetra::MultiVector<SC,LO,GO,NO>&          MV2TpetraMV(const MultiVector& vec)         { return Utilities::MV2TpetraMV(vec); }
    static       Tpetra::MultiVector<SC,LO,GO,NO>&          MV2NonConstTpetraMV(MultiVector& vec)       { return Utilities::MV2NonConstTpetraMV(vec); }

    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2TpetraCrs(RCP<const Matrix> Op)          { return Utilities::Op2TpetraCrs(Op); }
    static RCP<      Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraCrs(RCP<Matrix> Op)        { return Utilities::Op2NonConstTpetraCrs(Op); }

    static const Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2TpetraCrs(const Matrix& Op)              { return Utilities::Op2TpetraCrs(Op); }
    static       Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2NonConstTpetraCrs(Matrix& Op)            { return Utilities::Op2NonConstTpetraCrs(Op); }

    static RCP<const Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2TpetraRow(RCP<const Matrix> Op)          { return Utilities::Op2TpetraRow(Op); }
    static RCP<      Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraRow(RCP<Matrix> Op)        { return Utilities::Op2NonConstTpetraRow(Op); }

    static const RCP<const Tpetra::Map<LO, GO, NO> >        Map2TpetraMap(const Map& map)               { return Utilities::Map2TpetraMap(map); }
#endif

    static RCP<Xpetra::Matrix<SC,LO,GO,NO> >                Crs2Op(RCP<CrsMatrix> Op)                   { return Utilities::Crs2Op(Op); }

     /*! @brief Extract Matrix Diagonal

    Returns Matrix diagonal in RCP<Vector>.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixDiagonal(const Matrix& A); // FIXME

    /*! @brief Extract Matrix Diagonal

    Returns inverse of the Matrix diagonal in RCP<Vector>.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = TST::eps()*100, const bool doLumped = false); // FIXME



    /*! @brief Extract Overlapped Matrix Diagonal

    Returns overlapped Matrix diagonal in ArrayRCP.

    The local overlapped diagonal has an entry for each index in A's column map.
    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixOverlappedDiagonal(const Matrix& A); // FIXME

    /*! @brief Return vector containing inverse of input vector
     *
     * @param[in] v: input vector
     * @param[in] tol: tolerance. If entries of input vector are smaller than tolerance they are replaced by tolReplacement (see below). The default value for tol is 100*eps (machine precision)
     * @param[in] tolReplacement: Value put in for undefined entries in output vector (default: 0.0)
     * @ret: vector containing inverse values of input vector v
    */
    static Teuchos::RCP<Vector> GetInverse(Teuchos::RCP<const Vector> v, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100, Scalar tolReplacement = Teuchos::ScalarTraits<Scalar>::zero()) { return MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetInverse(v,tol,tolReplacement); }

    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
    static Teuchos::Array<Magnitude> ResidualNorm(const Operator& Op, const MultiVector& X, const MultiVector& RHS) {
      return Utilities::ResidualNorm(Op, X, RHS);
    }

    static RCP<MultiVector> Residual(const Operator& Op, const MultiVector& X, const MultiVector& RHS) {
      return Utilities::Residual(Op, X, RHS);
    }

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
    static SC PowerMethod(const Matrix& A, bool scaleByDiag = true,
                          LO niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      return Utilities::PowerMethod(A, scaleByDiag, niters, tolerance, verbose, seed);
    }

    static SC PowerMethod(const Matrix& A, const Teuchos::RCP<Vector> &invDiag,
                          LO niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      return Utilities::PowerMethod(A, invDiag, niters, tolerance, verbose, seed);
    }

    static void MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const SC>& scalingVector, bool doInverse = true,
                                 bool doFillComplete = true, bool doOptimizeStorage = true); // FIXME

    static void MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                                        bool doFillComplete, bool doOptimizeStorage); // FIXME

    static void MyOldScaleMatrix_Epetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                                            bool doFillComplete, bool doOptimizeStorage); // FIXME

    static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream& os)       { return Utilities::MakeFancy(os); }

    /*! @brief Detect Dirichlet rows

        @param[in] A matrix
        @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.

        @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
    */
    static Kokkos::View<bool*, typename NO::device_type> DetectDirichletRows(const Matrix& A, const Magnitude& tol = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<SC>::magnitudeType>::zero(), const bool count_twos_as_dirichlet=false);


    /*! @brief Find non-zero values in an ArrayRCP
      Compares the value to 2 * machine epsilon

      @param[in]  vals - ArrayRCP<const Scalar> of values to be tested
      @param[out] nonzeros - ArrayRCP<bool> of true/false values for whether each entry in vals is nonzero
    */
    static void FindNonZeros(const typename Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dual_view_type::t_dev_const_um vals,
                             Kokkos::View<bool*, typename Node::device_type> nonzeros);


    /*! @brief Detect Dirichlet columns based on Dirichlet rows

        The routine finds all column indices that are in Dirichlet rows, where Dirichlet rows are described by dirichletRows,
        as returned by DetectDirichletRows.

        @param[in] A matrix
        @param[in] dirichletRows array of Dirichlet rows.

        @return boolean array.  The ith entry is true iff row i is a Dirichlet column.
    */
    static Kokkos::View<bool*, typename NO::device_type> DetectDirichletCols(const Matrix& A, const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows);


    /*! @brief Detects Dirichlet columns & domains from a list of Dirichlet rows

      @param[in] A - Matrix on which to apply Dirichlet column detection
      @param[in] dirichletRows - View<bool> of indicators as to which rows are Dirichlet
      @param[out] dirichletCols - View<bool> of indicators as to which cols are Dirichlet
      @param[out] dirichletDomain - View<bool> of indicators as to which domains are Dirichlet
    */
    static void DetectDirichletColsAndDomains(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
                                              const Kokkos::View<bool*, typename Node::device_type> & dirichletRows,
                                              Kokkos::View<bool*, typename Node::device_type> dirichletCols,
                                              Kokkos::View<bool*, typename Node::device_type> dirichletDomain);


    static void ZeroDirichletRows(RCP<Matrix>& A, const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows, SC replaceWith=Teuchos::ScalarTraits<SC>::zero());

    static void ZeroDirichletRows(RCP<MultiVector>& X, const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows, SC replaceWith=Teuchos::ScalarTraits<SC>::zero());

    static void ZeroDirichletCols(RCP<Matrix>& A, const Kokkos::View<const bool*, typename NO::device_type>& dirichletCols, SC replaceWith=Teuchos::ScalarTraits<SC>::zero());

    static void ApplyRowSumCriterion(const Matrix& A,
                                     const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol,
                                     Kokkos::View<bool*, typename NO::device_type> & dirichletRows);

    static RCP<MultiVector> RealValuedToScalarMultiVector(RCP<RealValuedMultiVector> X);

    /*! @brief Set seed for random number generator.

      Distribute the seeds evenly in [1,INT_MAX-1].  This guarantees nothing
      about where random number streams on difference processes will intersect.
      This does avoid overflow situations in parallel when multiplying by a PID.
      It also avoids the pathological case of having the *same* random number stream
      on each process.
    */

    static void SetRandomSeed(const Teuchos::Comm<int> &comm)                       { return Utilities::SetRandomSeed(comm); }

    /*! @brief Transpose a Xpetra::Matrix

    Note: Currently, an error is thrown if the matrix isn't a Tpetra::CrsMatrix or Epetra_CrsMatrix.
    In principle, however, we could allow any Epetra_RowMatrix because the Epetra transposer does.
    */
    static RCP<Matrix> Transpose(Matrix& Op, bool optimizeTranspose = false, const std::string & label = std::string()) {
      return Utilities::Transpose(Op, optimizeTranspose, label);
    }

    static RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LocalOrdinal,GlobalOrdinal,Node> > ExtractCoordinatesFromParameterList(ParameterList& paramList) {
      return Utilities::ExtractCoordinatesFromParameterList(paramList);
    }

    /*! Perform a Cuthill-McKee (CM) or Reverse Cuthill-McKee (RCM) ordering of the local component of the matrix.
      Kokkos-Kernels has an RCM implementation, so we reverse that here if we call CM.
    */
    static RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > CuthillMcKee(const Matrix &Op);

    /*! Perform a Reverse Cuthill-McKee (RCM) ordering of the local component of the matrix.
     */
    static RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > ReverseCuthillMcKee(const Matrix &Op);

    static void ApplyOAZToMatrixRows(RCP<Matrix>& A, const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows);

  }; // class Utils


  /*!
        @class Utilities
        @brief MueLu utility class (specialization SC=double and LO=GO=int).

        This class provides a number of static helper methods. Some are temporary and will eventually
        go away, while others should be moved to Xpetra.

        Note: this is the implementation for Epetra. Tpetra throws if TPETRA_INST_INT_INT is disabled!
   */
  template <class Node>
  class Utilities_kokkos<double,int,int,Node> : public UtilitiesBase<double,int,int,Node> {
  public:
    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
    using CoordinateType        = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
    using RealValuedMultiVector = Xpetra::MultiVector<CoordinateType,LocalOrdinal,GlobalOrdinal,Node>;

  private:
#undef MUELU_UTILITIES_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

#ifdef HAVE_MUELU_EPETRA
    //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
    // @{
    static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<MultiVector> const vec)     { return Utilities::MV2EpetraMV(vec); }
    static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<MultiVector> vec)   { return Utilities::MV2NonConstEpetraMV(vec); }

    static const Epetra_MultiVector&                        MV2EpetraMV(const MultiVector& vec)         { return Utilities::MV2EpetraMV(vec); }
    static       Epetra_MultiVector&                        MV2NonConstEpetraMV(MultiVector& vec)       { return Utilities::MV2NonConstEpetraMV(vec); }

    static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Matrix> Op)          { return Utilities::Op2EpetraCrs(Op); }
    static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Matrix> Op)        { return Utilities::Op2NonConstEpetraCrs(Op); }

    static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Matrix& Op)              { return Utilities::Op2EpetraCrs(Op); }
    static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Matrix& Op)            { return Utilities::Op2NonConstEpetraCrs(Op); }

    static const Epetra_Map&                                Map2EpetraMap(const Map& map)               { return Utilities::Map2EpetraMap(map); }
    // @}
#endif

#ifdef HAVE_MUELU_TPETRA
    //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
    // @{
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> >     MV2TpetraMV(RCP<MultiVector> const vec)     { return Utilities::MV2TpetraMV(vec); }
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV(RCP<MultiVector> vec)   { return Utilities::MV2NonConstTpetraMV(vec); }
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV2(MultiVector& vec)      { return Utilities::MV2NonConstTpetraMV2(vec); }

    static const Tpetra::MultiVector<SC,LO,GO,NO>&          MV2TpetraMV(const MultiVector& vec)         { return Utilities::MV2TpetraMV(vec); }
    static       Tpetra::MultiVector<SC,LO,GO,NO>&          MV2NonConstTpetraMV(MultiVector& vec)       { return Utilities::MV2NonConstTpetraMV(vec); }

    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2TpetraCrs(RCP<const Matrix> Op)          { return Utilities::Op2TpetraCrs(Op); }
    static RCP<      Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraCrs(RCP<Matrix> Op)        { return Utilities::Op2NonConstTpetraCrs(Op); }

    static const Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2TpetraCrs(const Matrix& Op)              { return Utilities::Op2TpetraCrs(Op); }
    static       Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2NonConstTpetraCrs(Matrix& Op)            { return Utilities::Op2NonConstTpetraCrs(Op); }

    static RCP<const Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2TpetraRow(RCP<const Matrix> Op)          { return Utilities::Op2TpetraRow(Op); }
    static RCP<      Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraRow(RCP<Matrix> Op)        { return Utilities::Op2NonConstTpetraRow(Op); }

    static const RCP<const Tpetra::Map<LO, GO, NO> >        Map2TpetraMap(const Map& map)               { return Utilities::Map2TpetraMap(map); }
#endif
    static RCP<Xpetra::Matrix<SC,LO,GO,NO> >                Crs2Op(RCP<CrsMatrix> Op)                   { return Utilities::Crs2Op(Op); }

    static RCP<Vector> GetMatrixDiagonal(const Matrix& A) {
      const auto rowMap = A.getRowMap();
      auto diag = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap,true);

      A.getLocalDiagCopy(*diag);

      return diag;
    }
    static RCP<Vector> GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = Teuchos::ScalarTraits<SC>::eps()*100, const bool doLumped=false);

    static RCP<Vector> GetLumpedMatrixDiagonal(Matrix const &A, const bool doReciprocal=false, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100, Scalar tolReplacement = Teuchos::ScalarTraits<Scalar>::zero(), const bool replaceSingleEntryRowWithZero = false, const bool useAverageAbsDiagVal = false) {
      return UtilitiesBase::GetLumpedMatrixDiagonal(A, doReciprocal, tol, tolReplacement, replaceSingleEntryRowWithZero, useAverageAbsDiagVal);
    }
    static RCP<Vector> GetMatrixOverlappedDiagonal(const Matrix& A) {
      return UtilitiesBase::GetMatrixOverlappedDiagonal(A);
    }
    static RCP<Vector> GetInverse(RCP<const Vector> v, Magnitude tol = Teuchos::ScalarTraits<SC>::eps()*100, SC tolReplacement = Teuchos::ScalarTraits<SC>::zero()) {
      return UtilitiesBase::GetInverse(v,tol,tolReplacement);
    }
    static Teuchos::Array<Magnitude> ResidualNorm(const Operator& Op, const MultiVector& X, const MultiVector& RHS) {
      return UtilitiesBase::ResidualNorm(Op,X,RHS);
    }
    static RCP<MultiVector> Residual(const Operator& Op, const MultiVector& X, const MultiVector& RHS) {
      return UtilitiesBase::Residual(Op,X,RHS);
    }
    static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream& os) {
      return UtilitiesBase::MakeFancy(os);
    }
    static void SetRandomSeed(const Teuchos::Comm<int> &comm) {
      UtilitiesBase::SetRandomSeed(comm);
    }

    // todo: move this to UtilitiesBase::kokkos
    static Kokkos::View<bool*, typename Node::device_type> DetectDirichletRows(const Matrix& A, const Magnitude& tol = Teuchos::ScalarTraits<SC>::zero(), const bool count_twos_as_dirichlet=false);

    static Kokkos::View<bool*, typename Node::device_type> DetectDirichletCols(const Matrix& A, const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows);

    static void FindNonZeros(const typename Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dual_view_type::t_dev_const_um vals,
                             Kokkos::View<bool*, typename Node::device_type> nonzeros);

    static void DetectDirichletColsAndDomains(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
                                              const Kokkos::View<bool*, typename Node::device_type> & dirichletRows,
                                              Kokkos::View<bool*, typename Node::device_type> dirichletCols,
                                              Kokkos::View<bool*, typename Node::device_type> dirichletDomain);

    static void ZeroDirichletRows(RCP<Matrix>& A, const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows, SC replaceWith=Teuchos::ScalarTraits<SC>::zero());

    static void ZeroDirichletRows(RCP<MultiVector>& X, const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows, SC replaceWith=Teuchos::ScalarTraits<SC>::zero());

    static void ZeroDirichletCols(RCP<Matrix>& A, const Kokkos::View<const bool*, typename Node::device_type>& dirichletCols, SC replaceWith=Teuchos::ScalarTraits<SC>::zero());

    static void ApplyRowSumCriterion(const Matrix& A,
                                     const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol,
                                     Kokkos::View<bool*, typename NO::device_type> & dirichletRows);

    static RCP<MultiVector> RealValuedToScalarMultiVector(RCP<RealValuedMultiVector> X);

    static Scalar PowerMethod(const Matrix& A, bool scaleByDiag = true, LO niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      return UtilitiesBase::PowerMethod(A,scaleByDiag,niters,tolerance,verbose,seed);
    }

    static Scalar PowerMethod(const Matrix& A, const Teuchos::RCP<Vector> &invDiag, LO niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      return UtilitiesBase::PowerMethod(A, invDiag, niters, tolerance, verbose, seed);
    }

    static void MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const SC>& scalingVector, bool doInverse = true,
        bool doFillComplete = true, bool doOptimizeStorage = true) {
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
        MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;

      default:
        throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled.");
#ifndef __NVCC__ //prevent nvcc warning
        break;
#endif
      }
    }

    // TODO This is the <double,int,int> specialization
    static void MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
        bool doFillComplete, bool doOptimizeStorage) {
  #ifdef HAVE_MUELU_TPETRA
  #ifdef HAVE_MUELU_TPETRA_INST_INT_INT
      try {
        Tpetra::CrsMatrix<SC,LO,GO,NO>& tpOp = Op2NonConstTpetraCrs(Op);

        const RCP<const Tpetra::Map<LO,GO,NO> > rowMap    = tpOp.getRowMap();
        const RCP<const Tpetra::Map<LO,GO,NO> > domainMap = tpOp.getDomainMap();
        const RCP<const Tpetra::Map<LO,GO,NO> > rangeMap  = tpOp.getRangeMap();

        size_t maxRowSize = tpOp.getLocalMaxNumRowEntries();
        if (maxRowSize == Teuchos::as<size_t>(-1)) // hasn't been determined yet
          maxRowSize = 20;

        std::vector<SC> scaledVals(maxRowSize);
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
              scaledVals[j] = vals[j]*scalingVector[i];

            if (nnz > 0) {
	      Teuchos::ArrayView<const LocalOrdinal> cols_view(cols.data(), nnz);
              Teuchos::ArrayView<const SC> valview(&scaledVals[0], nnz);
              tpOp.replaceLocalValues(i, cols_view, valview);
            }
          } //for (size_t i=0; ...

        } else {
          typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::global_inds_host_view_type cols;
          typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type vals;

          for (size_t i = 0; i < rowMap->getLocalNumElements(); ++i) {
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
	      Teuchos::ArrayView<const GlobalOrdinal> cols_view(cols.data(), nnz);
              Teuchos::ArrayView<const SC> valview(&scaledVals[0], nnz);
              tpOp.replaceGlobalValues(gid, cols_view, valview);
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

    static void MyOldScaleMatrix_Epetra(Matrix& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector, bool doFillComplete, bool doOptimizeStorage) {
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
    static RCP<Matrix> Transpose(Matrix& Op, bool /* optimizeTranspose */ = false, const std::string & label = std::string(),const Teuchos::RCP<Teuchos::ParameterList> &params=Teuchos::null) {
      switch (Op.getRowMap()->lib()) {
      case Xpetra::UseTpetra:
      {
#ifdef HAVE_MUELU_TPETRA
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
        try {
          const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpetraOp = Utilities::Op2TpetraCrs(Op);

          // Compute the transpose A of the Tpetra matrix tpetraOp.
          RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
          Tpetra::RowMatrixTransposer<SC,LO,GO,NO> transposer(rcpFromRef(tpetraOp),label);
          {
            using Teuchos::ParameterList;
            using Teuchos::rcp;
            RCP<ParameterList> transposeParams = params.is_null () ?
              rcp (new ParameterList) :
              rcp (new ParameterList (*params));
            transposeParams->set ("sort", false);
            A = transposer.createTranspose(transposeParams);
          }

          RCP<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO> > AA   = rcp(new Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>(A));
          RCP<CrsMatrix>                             AAA  = rcp_implicit_cast<CrsMatrix>(AA);
          RCP<CrsMatrixWrap>                         AAAA = rcp(new CrsMatrixWrap(AAA));

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
#ifndef __NVCC__ //prevent nvcc warning
        break;
#endif
      }
      case Xpetra::UseEpetra:
      {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
        Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("ZZ Entire Transpose"));
        // Epetra case
        Epetra_CrsMatrix& epetraOp = Utilities::Op2NonConstEpetraCrs(Op);
        EpetraExt::RowMatrix_Transpose transposer;
        Epetra_CrsMatrix * A = dynamic_cast<Epetra_CrsMatrix*>(&transposer(epetraOp));
        transposer.ReleaseTranspose(); // So we can keep A in Muelu...

        RCP<Epetra_CrsMatrix> rcpA(A);
        RCP<Xpetra::EpetraCrsMatrixT<GO,NO> >  AA   = rcp(new Xpetra::EpetraCrsMatrixT<GO,NO> (rcpA));
        RCP<CrsMatrix>                         AAA  = rcp_implicit_cast<CrsMatrix>(AA);
        RCP<CrsMatrixWrap>                     AAAA = rcp( new CrsMatrixWrap(AAA));
        AAAA->fillComplete(Op.getRangeMap(), Op.getDomainMap());

        return AAAA;
#else
        throw Exceptions::RuntimeError("Epetra (Err. 2)");
#endif
#ifndef __NVCC__ //prevent nvcc warning
        break;
#endif
      }
      default:
        throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be transposed.");
#ifndef __NVCC__ //prevent nvcc warning
        break;
#endif
      }

#ifndef __NVCC__ //prevent nvcc warning
      return Teuchos::null;
#endif
    }

    /*! @brief Extract coordinates from parameter list and return them in a Xpetra::MultiVector
    */
    static RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO> > ExtractCoordinatesFromParameterList(ParameterList& paramList) {
      RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO> > coordinates = Teuchos::null;

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
      typedef Tpetra::MultiVector<float, LO,GO,NO> tfMV;
      RCP<tfMV> floatCoords = Teuchos::null;
#endif

      // define Tpetra::MultiVector type with Scalar=double only if
      // * ETI is turned off, since then the compiler will instantiate it automatically OR
      // * Tpetra is instantiated on Scalar=double
      typedef Tpetra::MultiVector<SC,LO,GO,NO> tdMV;
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
        coordinates = Teuchos::rcp(new Xpetra::TpetraMultiVector<SC,LO,GO,NO>(doubleCoords));
        TEUCHOS_TEST_FOR_EXCEPT(doubleCoords->getNumVectors() != coordinates->getNumVectors());
      }
#endif // Tpetra instantiated on GO=int and EpetraNode
#endif // endif HAVE_TPETRA

#if defined(HAVE_MUELU_EPETRA)
      RCP<Epetra_MultiVector> doubleEpCoords;
      if (paramList.isType<RCP<Epetra_MultiVector> >("Coordinates")) {
        doubleEpCoords = paramList.get<RCP<Epetra_MultiVector> >("Coordinates");
        paramList.remove("Coordinates");
        RCP<Xpetra::EpetraMultiVectorT<GO,NO> > epCoordinates = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO,NO>(doubleEpCoords));
        coordinates = rcp_dynamic_cast<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO> >(epCoordinates);
        TEUCHOS_TEST_FOR_EXCEPT(doubleEpCoords->NumVectors() != Teuchos::as<int>(coordinates->getNumVectors()));
      }
#endif
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(coordinates));
      return coordinates;
    }

    /*! Perform a Cuthill-McKee (CM) or Reverse Cuthill-McKee (RCM) ordering of the local component of the matrix
      Kokkos-Kernels has an RCM implementation, so we reverse that here if we call CM.
     */
    static RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > ReverseCuthillMcKee(const Matrix &Op);

    /*! Perform a Cuthill-McKee (CM) or Reverse Cuthill-McKee (RCM) ordering of the local component of the matrix
      Kokkos-Kernels has an RCM implementation, so we reverse that here if we call CM.
    */
    static RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > CuthillMcKee(const Matrix &Op);

    static void ApplyOAZToMatrixRows(RCP<Matrix>& A, const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows);

  }; // class Utilities (specialization SC=double LO=GO=int)



  // Useful Kokkos conversions
  template < class View, unsigned AppendValue >
  struct AppendTrait {
    // static_assert(false, "Error: NOT a Kokkos::View");
  };

  template < class MT, unsigned T >
  struct CombineMemoryTraits {
    // empty
  };

  template < unsigned U, unsigned T>
  struct CombineMemoryTraits<Kokkos::MemoryTraits<U>, T> {
    typedef Kokkos::MemoryTraits<U|T> type;
  };

  template < class DataType, unsigned T, class... Pack >
  struct AppendTrait< Kokkos::View< DataType, Pack... >, T> {
    typedef Kokkos::View< DataType, Pack... > view_type;
    using type = Kokkos::View< DataType, typename view_type::array_layout, typename view_type::device_type, typename CombineMemoryTraits<typename view_type::memory_traits,T>::type >;
  };

} //namespace MueLu

#define MUELU_UTILITIES_KOKKOS_SHORT

#endif // MUELU_UTILITIES_KOKKOS_DECL_HPP
