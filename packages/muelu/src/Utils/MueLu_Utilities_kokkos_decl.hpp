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
#if defined(HAVE_MUELU_KOKKOS_REFACTOR)

#include <string>

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>

#ifdef HAVE_MUELU_EPETRA
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Xpetra_EpetraCrsMatrix_fwd.hpp>
#include <Xpetra_EpetraMultiVector_fwd.hpp>
#endif

#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_MUELU_TPETRA
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Xpetra_TpetraCrsMatrix_fwd.hpp>
#include <Xpetra_TpetraMultiVector_fwd.hpp>
#endif


namespace MueLu {

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
  class Utils_kokkos {
#undef MUELU_UTILITIES_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

#ifdef HAVE_MUELU_EPETRA
    //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
    // @{
    static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<MultiVector> const Vec)     { return Utils::MV2EpetraMV(Vec); }
    static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<MultiVector> Vec)   { return Utils::MV2NonConstEpetraMV(Vec); }

    static const Epetra_MultiVector&                        MV2EpetraMV(const MultiVector& Vec)         { return Utils::MV2EpetraMV(Vec); }
    static       Epetra_MultiVector&                        MV2NonConstEpetraMV(MultiVector& Vec)       { return Utils::MV2NonConstEpetraMV(Vec); }

    static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Matrix> Op)          { return Utils::Op2EpetraCrs(Op); }
    static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Matrix> Op)        { return Utils::Op2NonConstEpetraCrs(Op); }

    static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Matrix& Op)              { return Utils::Op2EpetraCrs(Op); }
    static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Matrix& Op)            { return Utils::Op2NonConstEpetraCrs(Op); }

    static const Epetra_Map&                                Map2EpetraMap(const Map& map)               { return Utils::Map2EpetraMap(map); }
    // @}
#endif

#ifdef HAVE_MUELU_TPETRA
    //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
    // @{
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> >     MV2TpetraMV(RCP<MultiVector> const Vec)     { return Utils::MV2TpetraMV(Vec); }
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV(RCP<MultiVector> Vec)   { return Utils::MV2NonConstTpetraMV(Vec); }
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV2(MultiVector& Vec)      { return Utils::MV2NonConstTpetraMV2(Vec); }

    static const Tpetra::MultiVector<SC,LO,GO,NO>&          MV2TpetraMV(const MultiVector& Vec)         { return Utils::MV2TpetraMV(Vec); }
    static       Tpetra::MultiVector<SC,LO,GO,NO>&          MV2NonConstTpetraMV(MultiVector& Vec)       { return Utils::MV2NonConstTpetraMV(Vec); }

    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2TpetraCrs(RCP<const Matrix> Op)          { return Utils::Op2TpetraCrs(Op); }
    static RCP<      Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraCrs(RCP<Matrix> Op)        { return Utils::Op2NonConstTpetraCrs(Op); }

    static const Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2TpetraCrs(const Matrix& Op)              { return Utils::Op2TpetraCrs(Op); }
    static       Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2NonConstTpetraCrs(Matrix& Op)            { return Utils::Op2NonConstTpetraCrs(Op); }

    static RCP<const Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2TpetraRow(RCP<const Matrix> Op)          { return Utils::Op2TpetraRow(Op); }
    static RCP<      Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraRow(RCP<Matrix> Op)        { return Utils::Op2NonConstTpetraRow(Op); }

    static const RCP<const Tpetra::Map<LO, GO, NO> >        Map2TpetraMap(const Map& map)               { return Utils::Map2TpetraMap(map); }
#endif

    static RCP<Xpetra::Matrix<SC,LO,GO,NO> >                Crs2Op(RCP<CrsMatrix> Op)                   { return Utils::Crs2Op(Op); }

    /*! @brief Helper function to do matrix-matrix multiply

    Returns C = AB.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
    */
    static RCP<Matrix> Multiply(const Matrix & A,
                                bool transposeA,
                                const Matrix & B,
                                bool transposeB,
                                //Teuchos::FancyOStream &fos = *(Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout))),
                                Teuchos::FancyOStream &fos,
                                bool callFillCompleteOnResult = true,
                                bool doOptimizeStorage        = true,
                                const std::string & label     = std::string()){
      return Utils::Multiply(A, transposeA, B, transposeB, Teuchos::null, fos, callFillCompleteOnResult, doOptimizeStorage,label);
    }

    static RCP<Matrix> Jacobi(Scalar omega,
                              const Vector& D,
                              const Matrix& A,
                              const Matrix& B,
                              RCP<Matrix> C_in,
                              Teuchos::FancyOStream &fos,
                              const std::string & label     = std::string()) {
      return Utils::Jacobi(omega, D, A, B, C_in, fos, label);
    }


    /*! @brief Helper function to do matrix-matrix multiply

    Returns C = AB.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B

    @param C_in advanced usage. Use Teuchos::null by default.
           When C_in is available, its memory is reused to build
           This is useful in the case of multiple solve: if the pattern of C does not change, we can keep the memory and pattern of previous C matrix (C_in)
           C_in is modified in place and is not valid after the call.
    @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
    @param doOptimizedStorage if true, optimize storage
    */
    static RCP<Matrix> Multiply(const Matrix& A,
                                bool transposeA,
                                const Matrix& B,
                                bool transposeB,
                                RCP<Matrix> C_in,
                                //Teuchos::FancyOStream &fos = *(Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)))
                                Teuchos::FancyOStream &fos,
                                bool callFillCompleteOnResult = true,
                                bool doOptimizeStorage        = true,
                                const std::string & label     = std::string()) {
      return Utils::Multiply(A, transposeA, B, transposeB, C_in, fos, callFillCompleteOnResult, doOptimizeStorage, label);
    }

    /*! @brief Helper function to do matrix-matrix multiply "in-place"

    Returns RCP to non-constant Xpetra::BlockedCrsMatrix.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param doOptimizeStorage if true, the resulting matrix should be fillComplete'd
    */
    static RCP<BlockedCrsMatrix> TwoMatrixMultiplyBlock(BlockedCrsMatrix& A, bool transposeA,
                                                        BlockedCrsMatrix& B, bool transposeB,
                                                        Teuchos::FancyOStream& fos,
                                                        bool doFillComplete    = true,
                                                        bool doOptimizeStorage = true) {
      Utils::TwoMatrixMultiplyBlock(A, transposeA, B, transposeB, fos, doFillComplete, doOptimizeStorage);
    }

    /*! @brief Extract Matrix Diagonal

    Returns Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetMatrixDiagonal(const Matrix& A); // FIXME

    /*! @brief Extract Matrix Diagonal

    Returns inverse of the Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = Teuchos::ScalarTraits<SC>::eps()*100); // FIXME



    /*! @brief Extract Matrix Diagonal of lumped matrix

    Returns Matrix diagonal of lumped matrix in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetLumpedMatrixDiagonal(const Matrix& A); // FIXME

    /*! @brief Extract Overlapped Matrix Diagonal

    Returns overlapped Matrix diagonal in ArrayRCP.

    The local overlapped diagonal has an entry for each index in A's column map.
    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixOverlappedDiagonal(const Matrix& A); // FIXME

    /*! @brief Left scale matrix by an arbitrary vector.

    Algorithmically, this left scales a matrix by a diagonal matrix.
    The inverse of a diagonal matrix can also be applied.

    @param Op matrix to be scaled
    @param scalingVector vector that represents diagonal matrix
    @doInverse Indicates whether the inverse of the diagonal matrix should be applied.  (Default is to use inverse.)
    */
    static void ScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doInverse = true); // FIXME


    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
    static Teuchos::Array<Magnitude> ResidualNorm(const Operator& Op, const MultiVector& X, const MultiVector& RHS) {
      return Utils::ResidualNorm(Op, X, RHS);
    }

    static RCP<MultiVector> Residual(const Operator& Op, const MultiVector& X, const MultiVector& RHS) {
      return Utils::Residual(Op, X, RHS);
    }

    // NOTE:
    // A better place for the Read/Write function is probably Xpetra

    //! Read/Write methods
    //@{
    /*! @brief Save map to file. */
    static void Write(const std::string& fileName, const Map& M)                { Utils::Write(fileName, M); }

    /*! @brief Save vector to file in Matrix Market format.  */
    static void Write(const std::string& fileName, const MultiVector& Vec)      { Utils::Write(fileName, Vec); }

    /*! @brief Save matrix to file in Matrix Market format. */
    static void Write(const std::string& fileName, const Matrix& Op)            { Utils::Write(fileName, Op); }

    //! @brief Read matrix from file in Matrix Market or binary format.
    static Teuchos::RCP<Matrix> Read(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm, bool binary = false) {
      return Utils::Read(fileName, lib, comm, binary);
    }

    /*! @brief Read matrix from file in Matrix Market or binary format.

        If only rowMap is specified, then it is used for the domainMap and rangeMap, as well.
    */
    static Teuchos::RCP<Matrix> Read(const std::string&   fileName,
                                     const RCP<const Map> rowMap,
                                           RCP<const Map> colMap           = Teuchos::null,
                                     const RCP<const Map> domainMap        = Teuchos::null,
                                     const RCP<const Map> rangeMap         = Teuchos::null,
                                     const bool           callFillComplete = true,
                                     const bool           binary           = false,
                                     const bool           tolerant         = false,
                                     const bool           debug            = false)
    { return Utils::Read(fileName, rowMap, colMap, domainMap, rangeMap, callFillComplete, binary, tolerant, debug); }
    //@}

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
    static SC PowerMethod(const Matrix& A, bool scaleByDiag = true,
                          LO niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      return Utils::PowerMethod(A, scaleByDiag, niters, tolerance, verbose, seed);
    }

    static void MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const SC>& scalingVector, bool doInverse = true,
                                 bool doFillComplete = true, bool doOptimizeStorage = true); // FIXME

    static void MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                                        bool doFillComplete, bool doOptimizeStorage); // FIXME

    static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream& os)       { return Utils::MakeFancy(os); }

    /*! @brief Squared distance between two rows in a multivector

       Used for coordinate vectors.
    */
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType Distance2(const MultiVector& v, LocalOrdinal i0, LocalOrdinal i1) { return Utils::Distance2(v, i0, i1); }

    /*! @brief Detect Dirichlet rows

        @param[in] A matrix
        @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.

        @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
    */
    static Teuchos::ArrayRCP<const bool> DetectDirichletRows(const Matrix& A, const Magnitude& tol = Teuchos::ScalarTraits<SC>::zero()); // FIXME

    /*! @brief Set seed for random number generator.

      Distribute the seeds evenly in [1,INT_MAX-1].  This guarantees nothing
      about where random number streams on difference processes will intersect.
      This does avoid overflow situations in parallel when multiplying by a PID.
      It also avoids the pathological case of having the *same* random number stream
      on each process.
    */

    static void SetRandomSeed(const Teuchos::Comm<int> &comm)                       { return Utils::SetRandomSeed(comm); }

    static void findDirichletRows(Teuchos::RCP<Matrix> A,
                                  std::vector<LO>& dirichletRows);          // FIXME
    static void findDirichletCols(Teuchos::RCP<Matrix> A,
                                  std::vector<LO>& dirichletRows,
                                  std::vector<LO>& dirichletCols);          // FIXME
    static void Apply_BCsToMatrixRows(Teuchos::RCP<Matrix>& A,
                                      std::vector<LO>& dirichletRows);          // FIXME
    static void Apply_BCsToMatrixCols(Teuchos::RCP<Matrix>& A,
                                      std::vector<LO>& dirichletCols);          // FIXME
    static void Remove_Zeroed_Rows(Teuchos::RCP<Matrix>& A, double tol=1.0e-14);            // FIXME

  }; // class Utils

  /*!
    @class Utils2
    @brief MueLu utility class.

    Separate class for utilities that need a specialization for Epetra.
  */
  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = KokkosClassic::DefaultNode::DefaultNodeType>
  class Utils2_kokkos {
#include "MueLu_UseShortNames.hpp"

  public:

    /*! @brief Transpose a Xpetra::Matrix

    Note: Currently, an error is thrown if the matrix isn't a Tpetra::CrsMatrix or Epetra_CrsMatrix.
    In principle, however, we could allow any Epetra_RowMatrix because the Epetra transposer does.
    */
    static RCP<Matrix> Transpose(Matrix& Op, bool optimizeTranspose = false, const std::string & label = std::string()) {
      return Utils2::Transpose(Op, optimizeTranspose, label);
    }

    /*! @brief Helper function to calculate B = alpha*A + beta*B.

    @param A      left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha  scalar multiplier for A
    @param B      right matrix operand
    @param beta   scalar multiplier for B

    @return sum in B.

    Note that B does not have to be fill-completed.
    */
    static void TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta) {
      Utils2::TwoMatrixAdd(A, transposeA, alpha, B, beta);
    }

    /*! @brief Helper function to calculate C = alpha*A + beta*B.

    @param A          left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha      scalar multiplier for A, defaults to 1.0
    @param B          right matrix operand
    @param transposeB indicate whether to use transpose of B
    @param beta       scalar multiplier for B, defaults to 1.0
    @param C          resulting sum

    It is up to the caller to ensure that the resulting matrix sum is fillComplete'd.
    */
    static void TwoMatrixAdd(const Matrix& A, bool transposeA, const SC& alpha,
                             const Matrix& B, bool transposeB, const SC& beta,
                             RCP<Matrix>& C,  Teuchos::FancyOStream &fos, bool AHasFixedNnzPerRow = false) {
      Utils2::TwoMatrixAdd(A, transposeA, alpha, B, transposeB, beta, C, fos, AHasFixedNnzPerRow);
    }

    static RCP<MultiVector> ReadMultiVector (const std::string& fileName, const RCP<const Map>& map) {
      return Utils2::ReadMultiVector(fileName, map);
    }
    static RCP<const Map>   ReadMap         (const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm) {
      return Utils2::ReadMap(fileName, lib, comm);
    }

  }; // class Utils2_kokkos

} //namespace MueLu

#define MUELU_UTILITIES_KOKKOS_SHORT

#endif // #if defined(HAVE_MUELU_KOKKOS_REFACTOR)

#endif // MUELU_UTILITIES_KOKKOS_DECL_HPP
