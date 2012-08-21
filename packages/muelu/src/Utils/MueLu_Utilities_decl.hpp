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

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ScalarTraits.hpp>

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_OperatorFactory_fwd.hpp>
#include <Xpetra_CrsOperator_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_BlockedCrsOperator_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

#ifdef HAVE_MUELU_EPETRA
namespace Xpetra {
  class EpetraCrsMatrix; // TODO: replace by include of _fwd.hpp
  //  class 
}

// needed because of inlined function
//TODO: remove inline function?
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_CrsOperator.hpp>

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
RCP<Xpetra::CrsOperator<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator(RCP<Epetra_CrsMatrix> &epAB);
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
    static RCP<const Epetra_MultiVector> MV2EpetraMV(RCP<MultiVector> const Vec); //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static RCP<Epetra_MultiVector> MV2NonConstEpetraMV(RCP<MultiVector> Vec); //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static Epetra_MultiVector& MV2NonConstEpetraMV(MultiVector &Vec); //MV2EpetraMV

    static Epetra_MultiVector const& MV2EpetraMV(MultiVector const &Vec); //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Xpetra::Operator.
   static RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<Operator> Op); //Op2EpetraCrs

    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Xpetra::Operator.
   static RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Operator> Op); //Op2NonConstEpetraCrs
#endif

#ifdef HAVE_MUELU_TPETRA
    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > MV2TpetraMV(RCP<MultiVector> const Vec); //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV(RCP<MultiVector> Vec); //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static Tpetra::MultiVector<SC,LO,GO,NO> & MV2NonConstTpetraMV(MultiVector &Vec); //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV2(MultiVector &Vec); //MV2TpetraMV

    static Tpetra::MultiVector<SC,LO,GO,NO>  const& MV2TpetraMV(MultiVector const &Vec); //MV2TpetraMV
    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Xpetra::Operator.
    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2TpetraCrs(RCP<Operator> Op); //Op2TpetraCrs

    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Xpetra::Operator.
   static RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2NonConstTpetraCrs(RCP<Operator> Op); //Op2NonConstTpetraCrs

#endif

    /*! @brief Helper function to do matrix-matrix multiply "in-place"

      Returns RCP to non-constant Xpetra::Operator.

      @param A left matrix
      @param transposeA if true, use the transpose of A
      @param B right matrix
      @param transposeB if true, use the transpose of B
      @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
    */
   static RCP<Operator> TwoMatrixMultiply(RCP<Operator> const &A, bool transposeA,
                                          RCP<Operator> const &B, bool transposeB,
                                          bool doFillComplete=true,
                                          bool doOptimizeStorage=true); //TwoMatrixMultiply()

#ifdef HAVE_MUELU_EPETRAEXT
   // Michael Gee's MLMultiply
   static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
            const Epetra_CrsMatrix& epB);

#endif //ifdef HAVE_MUELU_EPETRAEXT

   /*! @brief Helper function to do matrix-matrix multiply "in-place"

     Returns RCP to non-constant Xpetra::BlockedCrsOperator.

     @param A left matrix
     @param transposeA if true, use the transpose of A
     @param B right matrix
     @param transposeB if true, use the transpose of B
     @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
   */
  static RCP<BlockedCrsOperator> TwoMatrixMultiplyBlock(RCP<BlockedCrsOperator> const &A, bool transposeA,
                                         RCP<BlockedCrsOperator> const &B, bool transposeB,
                                         bool doFillComplete=true,
                                         bool doOptimizeStorage=true); // TwoMatrixMultiplyBlock

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

    static void MatrixPrint(RCP<Operator> const &Op);

    static void MatrixPrint(RCP<Operator> const &Op, std::string const &label);

    /*! @brief Get Operator Diagonal
     */
    static RCP<Operator> BuildMatrixDiagonal(RCP<Operator> const &A); //BuildMatrixDiagonal()

    /*! @brief Extract Operator Diagonal

        Returns Operator diagonal in ArrayRCP.

        Note -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetMatrixDiagonal(RCP<Operator> const &A); //GetMatrixDiagonal

    /*! @brief Left scale matrix by an arbitrary vector.

       Algorithmically, this left scales a matrix by a diagonal matrix.
       The inverse of a diagonal matrix can also be applied.

       @param Op matrix to be scaled
       @param scalingVector vector that represents diagonal matrix
       @doInverse Indicates whether the inverse of the diagonal matrix should be applied.  (Default is to use inverse.)
     */
   static void ScaleMatrix(RCP<Operator> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse=true); //ScaleMatrix()

#ifdef UNUSED // and does not work with SC=complex
    /*! @brief Get reciprocal of Operator diagonal
     */

   static RCP<Operator> BuildMatrixInverseDiagonal(RCP<Operator> const &A); //BuildMatrixInverseDiagonal()
#endif

   typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Operator const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Operator const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
   static Teuchos::Array<Magnitude>
   ResidualNorm(Operator const &Op, MultiVector const &X, MultiVector const &RHS);
    
    static RCP<MultiVector> Residual(Operator const &Op, MultiVector const &X, MultiVector const &RHS);

   /*! @brief Save matrix to file in Matrix Market format.

     TODO Move this to Xpetra?
   */
   static void Write(std::string const & fileName, Operator const & Op); //Write

#include <unistd.h>

   static void PauseForDebugger(); //PauseForDebugger


   /*! @brief Simple transpose for Tpetra::CrsMatrix types

      Note:  This is very inefficient, as it inserts one entry at a time.
   */
#ifdef HAVE_MUELU_TPETRA
   static RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > simple_Transpose(RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > const &A); //simple_Transpose
#endif // HAVE_MUELU_TPETRA

#ifdef HAVE_MUELU_EPETRAEXT
   /*! @brief Simple transpose for Epetra_CrsMatrix types

      Note:  This is very inefficient, as it inserts one entry at a time.
   */
   static RCP<Epetra_CrsMatrix> simple_EpetraTranspose(RCP<const Epetra_CrsMatrix> const &A); //simple_Transpose
#endif

    /*! @brief Power method.

      @param A matrix
      @param scaleByDiag if true, estimate the largest eigenvalue of \f$ D^; A \f$.
      @param niters maximum number of iterations
      @param tolerance stopping tolerance
      @verbose if true, print iteration information
      
      (Shamelessly grabbed from tpetra/examples.)
    */
    static Scalar PowerMethod(Operator const &A, bool scaleByDiag=true,
                              LO niters=10, Magnitude tolerance=1e-2, bool verbose=false, unsigned int seed = 123); //PowerMethod

   static void MyOldScaleMatrix(RCP<Operator> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse=true,
                                bool doFillComplete=true,
                                bool doOptimizeStorage=true); //ScaleMatrix()

   static void MyOldScaleMatrix_Tpetra(RCP<Operator> &Op, Teuchos::ArrayRCP<SC> const &scalingVector,
                               bool doFillComplete, bool doOptimizeStorage);

   static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream & os);

  }; // class Utils

#ifdef HAVE_MUELU_EPETRA
//This non-member templated function exists so that the matrix-matrix multiply will compile if Epetra, Tpetra, and ML are enabled.
template<class SC,class LO,class GO,class NO, class LMO>
RCP<Xpetra::CrsOperator<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator(RCP<Epetra_CrsMatrix> &epAB) {
   TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
   return Teuchos::null;
}

typedef Kokkos::DefaultNode::DefaultNodeType KDNT;
typedef Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps KDKSO;

//specialization for the case of ScalarType=double and LocalOrdinal=GlobalOrdinal=int
template<>
inline RCP<Xpetra::CrsOperator<double,int,int,KDNT,KDKSO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator<double,int,int,KDNT,KDKSO > (RCP<Epetra_CrsMatrix> &epAB) {
  RCP<Xpetra::EpetraCrsMatrix> tmpC1 = rcp(new Xpetra::EpetraCrsMatrix(epAB));
  RCP<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> > tmpC2 = rcp_implicit_cast<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> >(tmpC1);
  RCP<Xpetra::CrsOperator<double,int,int,KDNT,KDKSO> > tmpC3 = rcp(new Xpetra::CrsOperator<double,int,int,KDNT,KDKSO>(tmpC2));
  return tmpC3;
}
#endif

//! Little helper function to convert non-string types to strings
template<class T>
std::string toString(T const &what) {
  std::ostringstream buf; buf << what;
  return buf.str();
}

//RCP<Xpetra::CrsOperator<double,int,int,KDNT,KDKSO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator<double,int,int,KDNT,KDKSO > (RCP<Epetra_CrsMatrix> epAB)


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

    /*! @brief Transpose a Xpetra::Operator

      Note: Currently, an error is thrown if the matrix isn't a Tpetra::CrsMatrix or Epetra_CrsMatrix.
      In principle, however, we could allow any Epetra_RowMatrix because the Epetra transposer does.
    */
   static RCP<Operator> Transpose(RCP<Operator> const &Op, bool const & optimizeTranspose=false); //Transpose

   //! Scale an Epetra matrix.
   static void MyOldScaleMatrix_Epetra(RCP<Operator> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doFillComplete, bool doOptimizeStorage);

   //! @brief Add two Operators.
   static void TwoMatrixAdd(RCP<Operator> const &A, bool transposeA, SC alpha, RCP<Operator> &B, SC beta);

   //! Add two Operators.
   static void TwoMatrixAdd(RCP<Operator> const &A, bool const &transposeA, SC const &alpha,
                                     RCP<Operator> const &B, bool const &transposeB, SC const &beta,
                                     RCP<Operator> &C); //TwoMatrixAdd()
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
   typedef Xpetra::Operator<double,int,int,NO,LMO> Operator;

public:

   static RCP<Operator> Transpose(RCP<Operator> const &Op, bool const & optimizeTranspose=false); //Transpose
   static void MyOldScaleMatrix_Epetra(RCP<Operator> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doFillComplete, bool doOptimizeStorage);
   static void TwoMatrixAdd(RCP<Operator> const &A, bool transposeA, SC alpha, RCP<Operator> &B, SC beta);
   static void TwoMatrixAdd(RCP<Operator> const &A, bool const &transposeA, SC const &alpha,
                                     RCP<Operator> const &B, bool const &transposeB, SC const &beta,
                                     RCP<Operator> &C);
  }; //specialization to Scalar=double


} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DECL_HPP
