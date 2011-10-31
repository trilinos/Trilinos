#ifndef MUELU_UTILITIES_HPP
#define MUELU_UTILITIES_HPP

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Utils.hpp>

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_OperatorFactory.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include "Epetra_RowMatrixTransposer.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#endif

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Memory.hpp"

#ifdef HAVE_MUELU_EPETRAEXT
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#ifdef HAVE_MUELU_TPETRA
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include "TpetraExt_MatrixMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer_decl.hpp"
#include "Tpetra_RowMatrixTransposer_def.hpp"
#endif

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

#ifdef HAVE_MUELU_ML
#include "ml_operator.h"
#include "ml_epetra_utils.h"
#endif

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
  using Xpetra::EpetraCrsMatrix;   // TODO: mv in Xpetra_UseShortNamesScalar
  using Xpetra::EpetraMultiVector;
#endif

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
    
#include "MueLu_UseShortNames.hpp"

  public:
#ifdef HAVE_MUELU_EPETRA
    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static RCP<const Epetra_MultiVector> MV2EpetraMV(RCP<MultiVector> const Vec) ; //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static RCP<Epetra_MultiVector> MV2NonConstEpetraMV(RCP<MultiVector> Vec) ; //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static Epetra_MultiVector& MV2NonConstEpetraMV(MultiVector &Vec) ; //MV2EpetraMV

    static Epetra_MultiVector const& MV2EpetraMV(MultiVector const &Vec) ; //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Xpetra::Operator.
   static RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<Operator> Op) ; //Op2EpetraCrs


    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Xpetra::Operator.
   static RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Operator> Op) ; //Op2NonConstEpetraCrs
#endif

#ifdef HAVE_MUELU_TPETRA
    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > MV2TpetraMV(RCP<MultiVector> const Vec) ; //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV(RCP<MultiVector> Vec) ; //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static Tpetra::MultiVector<SC,LO,GO,NO> & MV2NonConstTpetraMV(MultiVector &Vec) ; //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV2(MultiVector &Vec) ; //MV2TpetraMV

    static Tpetra::MultiVector<SC,LO,GO,NO>  const& MV2TpetraMV(MultiVector const &Vec) ; //MV2TpetraMV
    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Xpetra::Operator.
    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2TpetraCrs(RCP<Operator> Op) ; //Op2TpetraCrs

    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Xpetra::Operator.
   static RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2NonConstTpetraCrs(RCP<Operator> Op) ; //Op2NonConstTpetraCrs

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
                                          bool doOptimizeStorage=true)
    ; //TwoMatrixMultiply()

#ifdef HAVE_MUELU_EPETRAEXT
   // Michael Gee's MLMultiply
   static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(RCP<Epetra_CrsMatrix> epA,
            RCP<Epetra_CrsMatrix> epB)
    ;
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
                                         bool doOptimizeStorage=true)
  ; // TwoMatrixMultiplyBlock

    /*! @brief Helper function to calculate B = alpha*A + beta*B.

      @param A      left matrix operand
      @param transposeA indicate whether to use transpose of A
      @param alpha  scalar multiplier for A
      @param B      right matrix operand
      @param beta   scalar multiplier for B

      @return sum in B.

      Note that B does not have to be fill-completed.
    */
   static void TwoMatrixAdd(RCP<Operator> const &A, bool transposeA, SC alpha, RCP<Operator> &B, SC beta)
   ; //TwoMatrixAdd()

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
   static void TwoMatrixAdd(RCP<Operator> const &A, bool const &transposeA, SC const &alpha,
                                     RCP<Operator> const &B, bool const &transposeB, SC const &beta,
                                     RCP<Operator> &C)
   ; //TwoMatrixAdd()

    static void MatrixPrint(RCP<Operator> const &Op) ;

    static void MatrixPrint(RCP<Operator> const &Op, std::string const &label) ;

    /*! @brief Get Operator Diagonal
     */
   static RCP<Operator> BuildMatrixDiagonal(RCP<Operator> const &A)
    ;