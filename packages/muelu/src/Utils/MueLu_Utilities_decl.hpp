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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_UTILITIES_DECL_HPP
#define MUELU_UTILITIES_DECL_HPP

#include <unistd.h> //necessary for "sleep" function in debugging methods

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ScalarTraits.hpp>

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>

#ifdef HAVE_MUELU_EPETRA
namespace Xpetra {
  class EpetraCrsMatrix; // TODO: replace by include of _fwd.hpp
  //  class
}

// needed because of inlined function
//TODO: remove inline function?
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#endif

#include "MueLu_Exceptions.hpp"

#ifdef HAVE_MUELU_EPETRAEXT
class Epetra_CrsMatrix;
class Epetra_MultiVector;
#endif

#ifdef HAVE_MUELU_TPETRA
#include <Xpetra_TpetraMultiVector_fwd.hpp>
#include <Xpetra_TpetraCrsMatrix_fwd.hpp>

namespace Tpetra {
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node> class MultiVector;
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps> class CrsMatrix;
}

#endif

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
  //defined after Utils class
  template<typename SC,typename LO,typename GO,typename NO, typename LMO>
  RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB);
#endif

  /*!
    @class Utils
    @brief MueLu utility class.

    This class provides a number of static helper methods.  Some are temporary and will eventually
    go away, while others should be moved to Xpetra.
  */
  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = Kokkos::DefaultNode::DefaultNodeType,
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
  class Utils {
#undef MUELU_UTILITIES_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
#ifdef HAVE_MUELU_EPETRA
    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static RCP<const Epetra_MultiVector> MV2EpetraMV(RCP<MultiVector> const Vec);

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static RCP<Epetra_MultiVector> MV2NonConstEpetraMV(RCP<MultiVector> Vec);

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static Epetra_MultiVector& MV2NonConstEpetraMV(MultiVector &Vec);

    static Epetra_MultiVector const& MV2EpetraMV(MultiVector const &Vec);

    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Xpetra::Matrix.
    static RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<const Matrix> Op);

    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Xpetra::Matrix.
    static RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Matrix> Op);
#endif

#ifdef HAVE_MUELU_TPETRA
    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > MV2TpetraMV(RCP<MultiVector> const Vec);

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV(RCP<MultiVector> Vec);

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static Tpetra::MultiVector<SC,LO,GO,NO> & MV2NonConstTpetraMV(MultiVector &Vec);

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV2(MultiVector &Vec);

    static Tpetra::MultiVector<SC,LO,GO,NO>  const& MV2TpetraMV(MultiVector const &Vec);
    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Xpetra::Matrix.
    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2TpetraCrs(RCP<Matrix> Op);

    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Xpetra::Matrix.
    static RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2NonConstTpetraCrs(RCP<Matrix> Op);

#endif

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
                                bool callFillCompleteOnResult = true,
                                bool doOptimizeStorage        = true,
                                bool allowMLMultiply          = true) {
      return Utils<SC,LO,GO,NO,LMO>::Multiply(A, transposeA, B, transposeB, Teuchos::null, callFillCompleteOnResult, doOptimizeStorage, allowMLMultiply);

    }


    /*! @brief Helper function to do matrix-matrix multiply

    Returns C = AB.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd

    @param C_in advanced usage. Use Teuchos::null by default.
           When C_in is available, its memory is reused to build
           This is useful in the case of multiple solve: if the pattern of C does not change, we can keep the memory and pattern of previous C matrix (C_in)
           C_in is modified in place and is not valid after the call.

           ML MxM multiply does not reuse the pattern of C_in. If a C_in matrix is provided, then it is ignored.
           This can create a memory penalty if both the useless C_in and the new C are present in memory at the same time
           => Do not enable ML MxM at the same time as the option "reuse pattern".
    */
    static RCP<Matrix> Multiply(const Matrix & A,
                                bool transposeA,
                                const Matrix & B,
                                bool transposeB,
                                RCP<Matrix> C_in,
                                bool callFillCompleteOnResult = true,
                                bool doOptimizeStorage        = true,
                                bool allowMLMultiply          = true);

#ifdef HAVE_MUELU_EPETRAEXT
    // Michael Gee's MLMultiply
    static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
                                                     const Epetra_CrsMatrix& epB);
#endif //ifdef HAVE_MUELU_EPETRAEXT

    /*! @brief Helper function to do matrix-matrix multiply "in-place"

    Returns RCP to non-constant Xpetra::BlockedCrsMatrix.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param doOptimizeStorage if true, the resulting matrix should be fillComplete'd
    */
    static RCP<BlockedCrsMatrix> TwoMatrixMultiplyBlock(RCP<BlockedCrsMatrix> const &A, bool transposeA,
                                                        RCP<BlockedCrsMatrix> const &B, bool transposeB,
                                                        bool doFillComplete    = true,
                                                        bool doOptimizeStorage = true);

    /*! @brief Helper function to calculate B = alpha*A + beta*B.

    @param A      left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha  scalar multiplier for A
    @param B      right matrix operand
    @param beta   scalar multiplier for B

    @return sum in B.

    Note that B does not have to be fill-completed.
    */

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

    /*! @brief Get Matrix Diagonal
     */
    static RCP<Matrix> BuildMatrixDiagonal(RCP<Matrix> const &A);

    /*! @brief Extract Matrix Diagonal

    Returns Matrix diagonal in ArrayRCP.

    Note -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetMatrixDiagonal(const Matrix &A);

    /*! @brief Extract Matrix Diagonal of lumped matrix

    Returns Matrix diagonal of lumped matrix in ArrayRCP.

    Note -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetLumpedMatrixDiagonal(const Matrix &A);

    /*! @brief Extract Overlapped Matrix Diagonal

    Returns overlapped Matrix diagonal in ArrayRCP.

    The local overlapped diagonal has an entry for each index in A's column map.
    Note -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixOverlappedDiagonal(const Matrix &A);

    /*! @brief Left scale matrix by an arbitrary vector.

    Algorithmically, this left scales a matrix by a diagonal matrix.
    The inverse of a diagonal matrix can also be applied.

    @param Op matrix to be scaled
    @param scalingVector vector that represents diagonal matrix
    @doInverse Indicates whether the inverse of the diagonal matrix should be applied.  (Default is to use inverse.)
    */
    static void ScaleMatrix(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse=true);

#ifdef UNUSED // and does not work with SC=complex
    /*! @brief Get reciprocal of Matrix diagonal
     */

    static RCP<Matrix> BuildMatrixInverseDiagonal(RCP<Matrix> const &A);
#endif

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
    static Teuchos::Array<Magnitude>
    ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS);

    static RCP<MultiVector> Residual(Matrix const &Op, MultiVector const &X, MultiVector const &RHS);

    /*! @brief Save matrix to file in Matrix Market format.
     TODO Move this to Xpetra?
    */
   static void Write(std::string const & fileName, Matrix const & Op); //Write

    /*! @brief Save vector to file in Matrix Market format.
     TODO Move this to Xpetra?
    */
   static void Write(std::string const & fileName, const MultiVector& x); // Write

    /*! @brief Save map to file in Matrix Market format.
     TODO Move this to Xpetra?
    */
   static void Write(std::string const & fileName, const Map& M); // Write

   //! @brief Read matrix from file in Matrix Market format.
   static Teuchos::RCP<Matrix> Read(std::string const & fileName, Xpetra::UnderlyingLib lib, RCP<const Teuchos::Comm<int> > const &comm);


    /*! @brief Read vector from file in Matrix Market format.
     TODO Move this to Xpetra?
    */
    static RCP<MultiVector> Read(std::string const & fileName,const RCP<const Map> &map);


    static void PauseForDebugger();


    /*! @brief Simple transpose for Tpetra::CrsMatrix types

    Note:  This is very inefficient, as it inserts one entry at a time.
    */
#ifdef HAVE_MUELU_TPETRA
    static RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > simple_Transpose(RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > const &A);
#endif // HAVE_MUELU_TPETRA

#ifdef HAVE_MUELU_EPETRAEXT
    /*! @brief Simple transpose for Epetra_CrsMatrix types

    Note:  This is very inefficient, as it inserts one entry at a time.
    */
    static RCP<Epetra_CrsMatrix> simple_EpetraTranspose(RCP<const Epetra_CrsMatrix> const &A);
#endif

    /*! @brief Power method.

    @param A matrix
    @param scaleByDiag if true, estimate the largest eigenvalue of \f$ D^; A \f$.
    @param niters maximum number of iterations
    @param tolerance stopping tolerance
    @verbose if true, print iteration information

    (Shamelessly grabbed from tpetra/examples.)
    */
    static Scalar PowerMethod(Matrix const &A, bool scaleByDiag=true,
                              LO niters=10, Magnitude tolerance=1e-2, bool verbose=false, unsigned int seed = 123);

    static void MyOldScaleMatrix(RCP<Matrix> &Op, Teuchos::ArrayRCP<const SC> scalingVector, bool doInverse=true,
                                 bool doFillComplete=true,
                                 bool doOptimizeStorage=true);

    static Teuchos::ArrayRCP<double> CoalesceCoordinates(Teuchos::ArrayRCP<double> coord, LocalOrdinal blksize);

    static void MyOldScaleMatrix_Tpetra(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector,
                                        bool doFillComplete, bool doOptimizeStorage);

    static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream & os);

    /*! @brief Squared distance between two rows in a multivector

       Used for coordinate vectors.
    */
    static typename Teuchos::ScalarTraits<SC>::magnitudeType Distance2(const MultiVector& v, LO i0, LO i1);

    /*! @brief Detect Dirichlet rows

        @param[in] A matrix
        @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.

        @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
    */
    static Teuchos::ArrayRCP<const bool> DetectDirichletRows(Matrix const &A, typename Teuchos::ScalarTraits<SC>::magnitudeType const &tol=Teuchos::ScalarTraits<SC>::zero());

    /*! @brief print matrix info
    */
    static std::string PrintMatrixInfo(const Matrix& A, const std::string& msgTag, RCP<const ParameterList> params = Teuchos::null);

  }; // class Utils

#ifdef HAVE_MUELU_EPETRA
  //This non-member templated function exists so that the matrix-matrix multiply will compile if Epetra, Tpetra, and ML are enabled.
  template<class SC,class LO,class GO,class NO, class LMO>
  RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
    return Teuchos::null;
  }

  typedef Kokkos::DefaultNode::DefaultNodeType KDNT;
  typedef Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps KDKSO;

  //specialization for the case of ScalarType=double and LocalOrdinal=GlobalOrdinal=int
  template<>
  inline RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<double,int,int,KDNT,KDKSO > (RCP<Epetra_CrsMatrix> &epAB) {
    RCP<Xpetra::EpetraCrsMatrix> tmpC1 = rcp(new Xpetra::EpetraCrsMatrix(epAB));
    RCP<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> > tmpC2 = rcp_implicit_cast<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> >(tmpC1);
    RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO> > tmpC3 = rcp(new Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO>(tmpC2));
    return tmpC3;
  }
#endif

  //! Little helper function to convert non-string types to strings
  template<class T>
  std::string toString(T const &what) {
    std::ostringstream buf; buf << what;
    return buf.str();
  }

  //RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<double,int,int,KDNT,KDKSO > (RCP<Epetra_CrsMatrix> epAB)


  /*!
    @class Utils2
    @brief MueLu utility class.

    Separate class for utilities that need a specialization for Epetra.
  */
  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = Kokkos::DefaultNode::DefaultNodeType,
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
  class Utils2 {

#include "MueLu_UseShortNames.hpp"

  public:

    /*! @brief Transpose a Xpetra::Matrix

    Note: Currently, an error is thrown if the matrix isn't a Tpetra::CrsMatrix or Epetra_CrsMatrix.
    In principle, however, we could allow any Epetra_RowMatrix because the Epetra transposer does.
    */
    static RCP<Matrix> Transpose(RCP<Matrix> const &Op, bool const & optimizeTranspose=false);

    //! Scale an Epetra matrix.
    static void MyOldScaleMatrix_Epetra(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doFillComplete, bool doOptimizeStorage);

    //! @brief Add two .
    static void TwoMatrixAdd(RCP<Matrix> const &A, bool transposeA, SC alpha, RCP<Matrix> &B, SC beta);

    //! Add two .
    static void TwoMatrixAdd(RCP<Matrix> const &A, bool const &transposeA, SC const &alpha,
                             RCP<Matrix> const &B, bool const &transposeB, SC const &beta,
                             RCP<Matrix> &C, bool const &AHasFixedNnzPerRow=false);
  }; // class Utils2

  // specialization Utils2 for SC=double, LO=GO=int
  template<>
  class Utils2<double,int,int>//, Kokkos::DefaultNode::DefaultNodeType,
               //Kokkos::DefaultKernels<double,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps >
  {
    typedef double SC;
    typedef int LO;
    typedef int GO;
    typedef Kokkos::DefaultNode::DefaultNodeType NO;
    typedef Kokkos::DefaultKernels<double,int,NO>::SparseOps LMO;
    typedef Xpetra::Matrix<double,int,int,NO,LMO> Matrix;

  public:

    static RCP<Matrix> Transpose(RCP<Matrix> const &Op, bool const & optimizeTranspose=false);
    static void MyOldScaleMatrix_Epetra(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doFillComplete, bool doOptimizeStorage);
    static void TwoMatrixAdd(RCP<Matrix> const &A, bool transposeA, SC alpha, RCP<Matrix> &B, SC beta);
    static void TwoMatrixAdd(RCP<Matrix> const &A, bool const &transposeA, SC const &alpha,
                             RCP<Matrix> const &B, bool const &transposeB, SC const &beta,
                             RCP<Matrix> &C, bool const &AHasFixedNnzPerRow=false);
  }; //specialization to Scalar=double


} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DECL_HPP
