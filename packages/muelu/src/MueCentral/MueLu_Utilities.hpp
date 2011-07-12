#ifndef MUELU_UTILITIES_HPP
#define MUELU_UTILITIES_HPP

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Cthulhu_Map.hpp>
#include <Cthulhu_CrsMatrix.hpp>
#include <Cthulhu_OperatorFactory.hpp>
#include <Cthulhu_Vector.hpp>
#include <Cthulhu_VectorFactory.hpp>
#include <Cthulhu_MultiVectorFactory.hpp>
#ifdef HAVE_MUELU_EPETRA
#include <Cthulhu_EpetraCrsMatrix.hpp>
#include <Cthulhu_EpetraVector.hpp>
#include <Cthulhu_EpetraMultiVector.hpp>
#endif

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Memory.hpp"

#ifdef HAVE_MUELU_EPETRAEXT
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#ifdef HAVE_MUELU_TPETRA
#include <Cthulhu_TpetraCrsMatrix.hpp>
#include <Cthulhu_TpetraVector.hpp>
#include <Cthulhu_TpetraMultiVector.hpp>
#include "Tpetra_MatrixMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer_decl.hpp"
#include "Tpetra_RowMatrixTransposer_def.hpp"
#endif

namespace MueLu {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
#ifdef HAVE_MUELU_EPETRA
  using Cthulhu::EpetraCrsMatrix;   // TODO: mv in Cthulhu_UseShortNamesScalar
  using Cthulhu::EpetraMultiVector;
#endif

/*!
  @class Utils
  @brief MueLu utility class.

  This class provides a number of static helper methods.  Some are temporary and will eventually
  go away, while others should be moved to Cthulhu.
*/
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
  class Utils {
    
#include "MueLu_UseShortNames.hpp"

  public:
#ifdef HAVE_MUELU_EPETRAEXT
    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Cthulhu::MultiVector.
    static RCP<const Epetra_MultiVector> MV2EpetraMV(RCP<MultiVector> const Vec) {
      //rcp<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
      RCP<const EpetraMultiVector > tmpVec;
      tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
      if (tmpVec == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::MultiVector to Cthulhu::EpetraMultiVector failed"));
      RCP<const Epetra_MultiVector> epVec = tmpVec->getEpetra_MultiVector();
      return epVec;
    } //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Cthulhu::MultiVector.
    static RCP<Epetra_MultiVector> MV2NonConstEpetraMV(RCP<MultiVector> Vec) {
      RCP<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
      if (tmpVec == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::MultiVector to Cthulhu::EpetraMultiVector failed"));
      RCP<Epetra_MultiVector> epVec = tmpVec->getEpetra_MultiVector();
      return epVec;
    } //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Cthulhu::MultiVector.
    static Epetra_MultiVector& MV2NonConstEpetraMV(MultiVector &Vec) {
      EpetraMultiVector const &tmpVec = dynamic_cast<EpetraMultiVector const&>(Vec);
      RCP<Epetra_MultiVector> epVec = tmpVec.getEpetra_MultiVector();
      return *epVec;
    } //MV2EpetraMV

    static Epetra_MultiVector const& MV2EpetraMV(MultiVector const &Vec) {
      EpetraMultiVector const &tmpVec = dynamic_cast<EpetraMultiVector const&>(Vec);
      RCP<Epetra_MultiVector const> epVec = tmpVec.getEpetra_MultiVector();
      return *epVec;
    } //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Cthulhu::Operator.
   static RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<Operator> Op) {
      RCP<const Epetra_CrsMatrix> A;
      // Get the underlying Epetra Mtx
      RCP<const CrsOperator> crsOp = rcp_dynamic_cast<const CrsOperator>(Op);
      if (crsOp == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::Operator to Cthulhu::CrsOperator failed"));
      RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
      const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::CrsMatrix to Cthulhu::EpetraCrsMatrix failed"));
      A = tmp_ECrsMtx->getEpetra_CrsMatrix();
      return A;
    } //Op2EpetraCrs


    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Cthulhu::Operator.
   static RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Operator> Op) {
      RCP<Epetra_CrsMatrix> A;
      // Get the underlying Epetra Mtx
      RCP<const CrsOperator> crsOp = rcp_dynamic_cast<const CrsOperator>(Op);
      if (crsOp == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::Operator to Cthulhu::CrsOperator failed"));
      RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
      const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::CrsMatrix to Cthulhu::EpetraCrsMatrix failed"));
      A = tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
      return A;
    } //Op2NonConstEpetraCrs
#endif

#ifdef HAVE_MUELU_TPETRA
    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Cthulhu::MultiVector.
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > MV2TpetraMV(RCP<MultiVector> const Vec) {
      //rcp<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
      RCP<const TpetraMultiVector > tmpVec;
      tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
      if (tmpVec == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::MultiVector to Cthulhu::TpetraMultiVector failed"));
      RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec->getTpetra_MultiVector();
      return tpVec;
    } //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Cthulhu::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV(RCP<MultiVector> Vec) {
      RCP<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
      if (tmpVec == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::MultiVector to Cthulhu::TpetraMultiVector failed"));
      RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec->getTpetra_MultiVector();
      return tpVec;
    } //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Cthulhu::MultiVector.
    static Tpetra::MultiVector<SC,LO,GO,NO> & MV2NonConstTpetraMV(MultiVector &Vec) {
      TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
      RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec.getTpetra_MultiVector();
      return *tpVec;
    } //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Cthulhu::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV2(MultiVector &Vec) {
      TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
      RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec.getTpetra_MultiVector();
      return tpVec;
    } //MV2TpetraMV

    static Tpetra::MultiVector<SC,LO,GO,NO>  const& MV2TpetraMV(MultiVector const &Vec) {
      TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
      RCP<Tpetra::MultiVector<SC,LO,GO,NO>  const> tpVec = tmpVec.getTpetra_MultiVector();
      return *tpVec;
    } //MV2TpetraMV
    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Cthulhu::Operator.
    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2TpetraCrs(RCP<Operator> Op) {
     RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A;
      // Get the underlying Tpetra Mtx
      RCP<const CrsOperator> crsOp = rcp_dynamic_cast<const CrsOperator>(Op);
      if (crsOp == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::Operator to Cthulhu::CrsOperator failed"));
      RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
      const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::CrsMatrix to Cthulhu::TpetraCrsMatrix failed"));
      A = tmp_ECrsMtx->getTpetra_CrsMatrix();
      return A;
    } //Op2TpetraCrs

    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Cthulhu::Operator.
   static RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2NonConstTpetraCrs(RCP<Operator> Op) {
      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A;
      // Get the underlying Tpetra Mtx
      RCP<const CrsOperator> crsOp = rcp_dynamic_cast<const CrsOperator>(Op);
      if (crsOp == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::Operator to Cthulhu::CrsOperator failed"));
      RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
      const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Cthulhu::CrsMatrix to Cthulhu::TpetraCrsMatrix failed"));
      A = tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
      return A;
    } //Op2NonConstTpetraCrs

#endif

    /*! @brief Helper function to do matrix-matrix multiply "in-place"

      Returns RCP to non-constant Cthulhu::Operator.

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
    {
      //FIXME 30 is likely a big overestimate
      RCP<Operator> C;
      if(transposeA) C = OperatorFactory::Build(A->getDomainMap(), 1);
      else C = OperatorFactory::Build(A->getRowMap(), 1);

      if (!A->isFillComplete())
        throw(Exceptions::RuntimeError("A is not fill-completed"));
      if (!B->isFillComplete())
        throw(Exceptions::RuntimeError("B is not fill-completed"));

      if (C->getRowMap()->lib() == Cthulhu::UseEpetra) {
#ifdef HAVE_MUELU_EPETRAEXT
        RCP<const Epetra_CrsMatrix> epA = Op2EpetraCrs(A);
        RCP<const Epetra_CrsMatrix> epB = Op2EpetraCrs(B);
        RCP<Epetra_CrsMatrix>       epC = Op2NonConstEpetraCrs(C);
        
        int i = EpetraExt::MatrixMatrix::Multiply(*epA,transposeA,*epB,transposeB,*epC,doFillComplete);
        
        if (i != 0) {
          std::ostringstream buf;
          buf << i;
        std::string msg = "EpetraExt::MatrixMatrix::Multiply return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
        }
#else
        throw(Exceptions::RuntimeError("MueLu must be compiled with EpetraExt."));
#endif
      } else if(C->getRowMap()->lib() == Cthulhu::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpA = Op2TpetraCrs(A);
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpB = Op2TpetraCrs(B);
        RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >       tpC = Op2NonConstTpetraCrs(C);
        
        if (!doOptimizeStorage) {
          Tpetra::MatrixMatrix::Multiply(*tpA,transposeA,*tpB,transposeB,*tpC,false);
          tpC->fillComplete((transposeB) ? tpB->getRangeMap() : tpB->getDomainMap(),
                             (transposeA) ? tpA->getDomainMap() : tpA->getRangeMap(),
                            Tpetra::DoNotOptimizeStorage);
        } else {
          Tpetra::MatrixMatrix::Multiply(*tpA,transposeA,*tpB,transposeB,*tpC,doFillComplete);
        }
#else
        throw(Exceptions::RuntimeError("MueLu must be compiled with Tpetra."));
#endif
      }

      return C;
    } //TwoMatrixMultiply()

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
   {
      if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
        throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
      }

      if (A->getRowMap()->lib() == Cthulhu::UseEpetra) {
#ifdef HAVE_MUELU_EPETRAEXT
        RCP<const Epetra_CrsMatrix> epA = Op2EpetraCrs(A);
        RCP<Epetra_CrsMatrix> epB = Op2NonConstEpetraCrs(B);
        
        //FIXME is there a bug if beta=0?
        int i = EpetraExt::MatrixMatrix::Add(*epA,transposeA,(double)alpha,*epB,(double)beta);

        if (i != 0) {
          std::ostringstream buf;
          buf << i;
          std::string msg = "EpetraExt::MatrixMatrix::Add return value of " + buf.str();
          throw(Exceptions::RuntimeError(msg));
        }
#else
      throw(Exceptions::RuntimeError("MueLu must be compile with EpetraExt."));
#endif
      } else if(A->getRowMap()->lib() == Cthulhu::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpA = Op2TpetraCrs(A);
        RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpB = Op2NonConstTpetraCrs(B);
        
        Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, beta);
#else
        throw(Exceptions::RuntimeError("MueLu must be compiled with Tpetra."));
#endif
      }

   } //TwoMatrixAdd()

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
   {
      if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
        throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
      }
      if (C==Teuchos::null)
        //FIXME 5 is a complete guess as to the #nonzeros per row
        C = rcp( new CrsOperator(A->getRowMap(), 5) );

      if (C->getRowMap()->lib() == Cthulhu::UseEpetra) {
#ifdef HAVE_MUELU_EPETRAEXT
        RCP<const Epetra_CrsMatrix> epA = Op2EpetraCrs(A);
        RCP<const Epetra_CrsMatrix> epB = Op2EpetraCrs(B);
        RCP<Epetra_CrsMatrix>       epC = Op2NonConstEpetraCrs(C);
        Epetra_CrsMatrix* ref2epC = &*epC; //to avoid a compiler error...

        //FIXME is there a bug if beta=0?
        int i = EpetraExt::MatrixMatrix::Add(*epA,transposeA,(double)alpha,*epB,transposeB,(double)beta,ref2epC);

        if (i != 0) {
          std::ostringstream buf;
          buf << i;
          std::string msg = "EpetraExt::MatrixMatrix::Add return value of " + buf.str();
          throw(Exceptions::RuntimeError(msg));
        }
#else
      throw(Exceptions::RuntimeError("MueLu must be compile with EpetraExt."));
#endif
      } else if(C->getRowMap()->lib() == Cthulhu::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpA = Op2TpetraCrs(A);
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpB = Op2TpetraCrs(B);
        RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >       tpC = Op2NonConstTpetraCrs(C);

        Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, transposeB, beta, tpC);
#else
        throw(Exceptions::RuntimeError("MueLu must be compile with Tpetra."));
#endif
      }

   } //TwoMatrixAdd()

    static void MatrixPrint(RCP<Operator> const &Op) {
      std::string label = "unlabeled operator";
      MatrixPrint(Op, label);
    }

    static void MatrixPrint(RCP<Operator> const &Op, std::string const &label) {
#ifdef HAVE_MUELU_EPETRA 
      RCP<const Epetra_CrsMatrix> epOp = Op2EpetraCrs(Op);
      int mypid = epOp->RowMap().Comm().MyPID();
      if (mypid == 0)
        std::cout << "\n===============\n" << label << "\n==============" << std::endl;

      if (mypid == 0) std::cout << "\n   -- row map -- \n" << std::endl;
      std::cout << epOp->RowMap() << std::endl;
      sleep(1);
      epOp->RowMap().Comm().Barrier();

      if (mypid == 0) std::cout << "\n   -- column map -- \n" << std::endl;
      std::cout << epOp->ColMap() << std::endl;
      sleep(1);
      epOp->RowMap().Comm().Barrier();

      std::cout << *epOp << std::endl;
#endif
    }

    /*! @brief Get Operator Diagonal
     */
   static RCP<Operator> BuildMatrixDiagonal(RCP<Operator> const &A)
    {
      const RCP<const Map> rowmap = A->getRowMap();
      //std::vector<SC> diag(A->getNodeNumRows());
      std::vector<SC> diag(rowmap->getNodeNumElements());
      Teuchos::ArrayView<const LO> cols;
      Teuchos::ArrayView<const SC> vals;
      //for (size_t i=0; i<A->getNodeNumRows(); ++i) {
      for (size_t i=0; i<rowmap->getNodeNumElements(); ++i) {
        A->getLocalRowView(i,cols,vals);
        //for (Teuchos::ArrayView<const LO>::size_type j=0; j<cols.size(); j++) {
        for (size_t j=0; j<cols.size(); j++) {
          //TODO this will break down if diagonal entry is not present
          //if (!(cols[j] > i)) {  //JG says this will work ... maybe
          if (cols[j] == i) {
            diag[i] = vals[j];
            break;
          }
        }
      }

      Teuchos::RCP< Operator > D = Teuchos::rcp( new CrsOperator(rowmap, 1) );
      std::vector<LO> diagInd(1);
      Teuchos::ArrayView<GO> iv(&diagInd[0],1);
      //for (size_t i=0; i< A->getNodeNumRows(); ++i) {
      for (size_t i=0; i< rowmap->getNodeNumElements(); ++i) {
        Teuchos::ArrayView<SC> av(&diag[i],1);
        diagInd[0] = rowmap->getGlobalElement(i);
        D->insertGlobalValues(i,iv,av);
      }
      D->fillComplete();
      //MatrixPrint(D);

      return D;

    } //BuildMatrixDiagonal()

    /*! @brief Extract Operator Diagonal

        Returns Operator diagonal in ArrayRCP.
    */
    static Teuchos::ArrayRCP<SC> GetMatrixDiagonal(RCP<Operator> const &A)
    {
      const RCP<const Map> rowmap = A->getRowMap();
      size_t locSize = rowmap->getNodeNumElements();
      Teuchos::ArrayRCP<SC> diag(locSize);
      Teuchos::ArrayView<const LO> cols;
      Teuchos::ArrayView<const SC> vals;
      for (size_t i=0; i<locSize; ++i) {
        A->getLocalRowView(i,cols,vals);
        for (LO j=0; j<cols.size(); ++j) {
          //TODO this will break down if diagonal entry is not present
          //if (!(cols[j] > i))   //JG says this will work ... maybe
          if (Teuchos::as<size_t>(cols[j]) == i) {
            diag[i] = vals[j];
            break;
          }
        }
      }
      //for (int i=0; i<locSize; ++i) std::cout << "diag[" << i << "] = " << diag[i] << std::endl;
      return diag;
    } //GetMatrixDiagonal

    /*! @brief Left scale matrix by an arbitrary vector.

       Algorithmically, this left scales a matrix by a diagonal matrix.
       The inverse of a diagonal matrix can also be applied.

       @param Op matrix to be scaled
       @param scalingVector vector that represents diagonal matrix
       @doInverse Indicates whether the inverse of the diagonal matrix should be applied.  (Default is to use inverse.)
     */
   static void ScaleMatrix(RCP<Operator> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse=true)
   {
      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpOp;
      try {
        tpOp = Op2NonConstTpetraCrs(Op);
      }
      catch(...) {
        throw(Exceptions::RuntimeError("Sorry, haven't implemented matrix scaling for epetra"));
      }
      Tpetra::Vector<SC,LO,GO,NO> x(tpOp->getRowMap(),scalingVector());
      if(doInverse){
        Tpetra::Vector<SC,LO,GO,NO> xi(tpOp->getRowMap());
        xi.reciprocal(x);
        tpOp->leftScale(xi);
      }
      else
        tpOp->leftScale(x);
   } //ScaleMatrix()

    /*! @brief Get reciprocal of Operator diagonal
     */

   static RCP<Operator> BuildMatrixInverseDiagonal(RCP<Operator> const &A)
    {
      const RCP<const Map> rowmap = A->getRowMap();
      //std::vector<SC> diag(A->getNodeNumRows());
      std::vector<SC> diag(rowmap->getNodeNumElements());
      Teuchos::ArrayView<const LO> cols;
      Teuchos::ArrayView<const SC> vals;
      //for (size_t i=0; i<A->getNodeNumRows(); ++i) {
      LO rowmapLocalSize = (LO) rowmap->getNodeNumElements();
      for (LO i=0; i<rowmapLocalSize; ++i) {
        A->getLocalRowView(i,cols,vals);
        for (LO j=0; j<cols.size(); ++j) {
          //TODO this will break down if diagonal entry is not present
          if (cols[j] == i) {
            diag[i] = 1 / vals[j];
            break;
          }
        }
      }

      RCP< Operator > D = Teuchos::rcp( new CrsOperator(rowmap, 1) );
      std::vector<GO> diagInd(1);
      Teuchos::ArrayView<GO> iv(&diagInd[0],1);
      //for (size_t i=0; i< A->getNodeNumRows(); ++i) {

      for (size_t i=0; i< rowmap->getNodeNumElements(); ++i) {
        Teuchos::ArrayView<SC> av(&diag[i],1);
        diagInd[0] = rowmap->getGlobalElement(i);
        D->insertGlobalValues(rowmap->getGlobalElement(i),iv,av); //TODO is this expensive?
      }
      D->fillComplete();

      return D;

    } //BuildMatrixInverseDiagonal()

   typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Operator const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Operator const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
   static Teuchos::Array<Magnitude>
   ResidualNorm(Operator const &Op, MultiVector const &X, MultiVector const &RHS)
   {
     //if (X.getNumVectors() != RHS.getNumVectors())
     //  throw(Exceptions::RuntimeError("Number of solution vectors != number of right-hand sides"));
     //const size_t numVecs = X.getNumVectors();
     const size_t numVecs = 1;
     RCP<MultiVector> RES = Residual(Op,X,RHS);
     Teuchos::Array<Magnitude> norms(numVecs);
     RES->norm2(norms);
     return norms;
   }

   static RCP<MultiVector> Residual(Operator const &Op, MultiVector const &X, MultiVector const &RHS)
   {
     SC one = 1.0;
     SC negone = -1.0;
     //if (X.getNumVectors() != RHS.getNumVectors())
     //  throw(Exceptions::RuntimeError("Number of solution vectors != number of right-hand sides"));
     //const size_t numVecs = X.getNumVectors();
     const size_t numVecs = 1;
     RCP<MultiVector> RES = MultiVectorFactory::Build(Op.getRangeMap(),numVecs);
     Op.apply(X,*RES,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
     RES->update(one,RHS,negone);
     return RES;
   }

   /*! @brief Save matrix to file in Matrix Market format.

     TODO Move this to Cthulhu?
   */
   static void Write(std::string const & fileName, Operator const & Op) {
    CrsOperator const & crsOp = dynamic_cast<CrsOperator const &>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp.getCrsMatrix();
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    const RCP<const TpetraCrsMatrix> &tmp_TCrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx != Teuchos::null) {

      RCP<const Epetra_CrsMatrix> A = tmp_ECrsMtx->getEpetra_CrsMatrix();
      int rv = EpetraExt::RowMatrixToMatrixMarketFile(fileName.c_str(), *A);
      if (rv != 0) {
        std::ostringstream buf;
        buf << rv;
      std::string msg = "EpetraExt::RowMatrixToMatrixMarketFile return value of " + buf.str();
      throw(Exceptions::RuntimeError(msg));
      }

    } else if (tmp_TCrsMtx != Teuchos::null) {

      RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A = tmp_TCrsMtx->getTpetra_CrsMatrix();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >::writeSparseFile(fileName,A);

    } else {

      throw(Exceptions::BadCast("Could not cast to EpetraCrsMatrix or TpetraCrsMatrix in matrix writing"));
    }

   } //Write

#include <unistd.h>

   static void PauseForDebugger()
   {
     RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
   
     int mypid = comm->getRank();
   
     for (int i = 0; i <comm->getSize() ; i++) {
       if (i == mypid ) {
         char buf[80];
         char hostname[80];
         gethostname(hostname, sizeof(hostname));
         int pid = getpid();
         sprintf(buf, "Host: %s\tMPI rank: %d,\tPID: %d\n\tattach %d\n\tcontinue\n",
             hostname, mypid, pid, pid);
         printf("%s\n",buf);
         fflush(stdout);
         sleep(1);
       }
     }
   
     if (mypid == 0) {
       printf( "** Enter a character to continue > "); fflush(stdout);
       char go = ' ';
       scanf("%c",&go);
     }
     comm->barrier();
   } //PauseForDebugger

    //! @brief Transpose a Cthulhu::Operator
   static RCP<Operator> Transpose(RCP<Operator> const &Op, bool const & optimizeTranspose=false)
   {
     RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpetraOp = Op2TpetraCrs(Op); //TODO JJH try/catch this
     //     Tpetra::RowMatrixTransposer<SC,LO,GO,NO,LMO> transposer(*tpetraOp); //more than meets the eye
     //     RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A = transposer.createTranspose(optimizeTranspose ? Tpetra::DoOptimizeStorage : Tpetra::DoNotOptimizeStorage); //couldn't have just used a bool...
     RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A=simple_Transpose(tpetraOp);
     RCP<TpetraCrsMatrix> AA = rcp(new TpetraCrsMatrix(A) );
     RCP<CrsMatrix> AAA = Teuchos::rcp_implicit_cast<CrsMatrix>(AA);
     RCP<CrsOperator> AAAA = rcp( new CrsOperator(AAA) );
     return AAAA;
   } //Transpose


   static RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > simple_Transpose(RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > const &A)
   {
      LocalOrdinal N=A->getNodeNumRows();
      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > AT=rcp(new Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>(A->getDomainMap(),0));
      const RCP<const Tpetra::Map<LO,GO,NO> > & rowMap=A->getRowMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > & colMap=A->getColMap();

      for(LO i=0;i<N;i++){
        GO grid= rowMap->getGlobalElement(i);
        Teuchos::ArrayRCP<const LO> indices;
        Teuchos::ArrayRCP<const SC> vals;
        A->getLocalRowView(i,indices,vals);
        for(LO j=0;j<indices.size();j++){
          GO gcid=colMap->getGlobalElement(indices[j]);
          AT->insertGlobalValues(gcid,Teuchos::tuple(grid),Teuchos::tuple(vals[j]));
        }
      }
      AT->fillComplete(A->getRangeMap(),A->getDomainMap());
      
      return AT;
    } //simple_Tranpose

    /*! @brief Power method.

      @param A matrix
      @param scaleByDiag if true, estimate the largest eigenvalue of \f$ D^{-1} A \f$.
      @param niters maximum number of iterations
      @param tolerance stopping tolerance
      @verbose if true, print iteration information
      
      (Shamelessly grabbed from tpetra/examples.)
    */
    static Scalar PowerMethod(Operator const &A, bool scaleByDiag=true,
                              LO niters=10, Magnitude tolerance=1e-2, bool verbose=false)
    {
      if ( !(A.getRangeMap()->isSameAs(*(A.getDomainMap()))) ) {
        throw(Exceptions::Incompatible("Utils::PowerMethod: operator must have domain and range maps that are equivalent."));
      }
      // create three vectors, fill z with random numbers
      RCP<MultiVector> q = MultiVectorFactory::Build(A.getRangeMap(),1);
      RCP<MultiVector> r = MultiVectorFactory::Build(A.getRangeMap(),1);
      RCP<MultiVector> z = MultiVectorFactory::Build(A.getRangeMap(),1);
      z->randomize();
      
      // typedef Teuchos::ScalarTraits<SC> ST;

      Teuchos::Array<Magnitude> norms(1);
  
      //std::vector<Scalar> lambda(1);
      //lambda[0] = 0.0;
      Scalar lambda=0.0;
      Magnitude residual = 0.0;
      // power iteration
      Teuchos::ArrayView<Scalar> avLambda(&lambda,1);
      RCP<Vector> diagVec,oneOverDiagonal;
      if (scaleByDiag) {
        diagVec = VectorFactory::Build(A.getRowMap());
        A.getLocalDiagCopy(*diagVec);
        oneOverDiagonal = VectorFactory::Build(A.getRowMap());
        oneOverDiagonal->reciprocal(*diagVec);
      }
      for (int iter = 0; iter < niters; ++iter) {
        z->norm2(norms);                               // Compute 2-norm of z
        q->update(1.0/norms[0],*z,0.);                 // Set q = z / normz
        A.apply(*q, *z);                               // Compute z = A*q
        if (scaleByDiag) z->elementWiseMultiply(1.0, *oneOverDiagonal, *z, 0.0);
        q->dot(*z,avLambda);                            // Approximate maximum eigenvalue: lamba = dot(q,z)
        if ( iter % 100 == 0 || iter + 1 == niters ) {
          r->update(1.0, *z, -lambda, *q, 0.0);         // Compute A*q - lambda*q
          r->norm2(norms);
          residual = Teuchos::ScalarTraits<Scalar>::magnitude(norms[0] / lambda);
          if (verbose) {
            std::cout << "Iter = " << iter
                      << "  Lambda = " << lambda
                      << "  Residual of A*q - lambda*q = " << residual
                      << std::endl;
          }
        }
        if (residual < tolerance) {
          break;
        }
      }
      return lambda;
    } //PowerMethod

   static void MyOldScaleMatrix(RCP<Operator> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse=true,
                                bool doFillComplete=true,
                                bool doOptimizeStorage=true)
   {

      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpOp;
      try {
        tpOp = Op2NonConstTpetraCrs(Op);
      }
      catch(...) {
        throw(Exceptions::RuntimeError("Sorry, haven't implemented matrix scaling for epetra"));
      }

      const RCP<const Tpetra::Map<LO,GO,NO> > rowMap = tpOp->getRowMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > domainMap = tpOp->getDomainMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > rangeMap = tpOp->getRangeMap();
      Teuchos::ArrayView<const LO> cols;
      Teuchos::ArrayView<const SC> vals;
      size_t maxRowSize = tpOp->getNodeMaxNumRowEntries();
      if (maxRowSize==(size_t)-1) //hasn't been determined yet
        maxRowSize=20;
      std::vector<SC> scaledVals(maxRowSize);
      if (tpOp->isFillComplete()) {
        std::cout << "In MyOldScale, resuming fill" << std::endl;
        tpOp->resumeFill();
      }

      Teuchos::ArrayRCP<SC> sv(scalingVector.size());
      if (doInverse) {
        for (int i=0; i<scalingVector.size(); ++i)
          sv[i] = 1.0 / scalingVector[i];
      } else {
        for (int i=0; i<scalingVector.size(); ++i)
          sv[i] = scalingVector[i];
      }

      if (Op->isLocallyIndexed() == true) {
        for (size_t i=0; i<rowMap->getNodeNumElements(); ++i) {
          tpOp->getLocalRowView(i,cols,vals);
          size_t nnz = tpOp->getNumEntriesInLocalRow(i);
          if (nnz>maxRowSize) {
            maxRowSize=nnz;
            scaledVals.resize(maxRowSize);
          }
          for (size_t j=0; j<nnz; ++j) {
            scaledVals[j] = vals[j]*sv[i];
          }
          if (nnz>0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0],nnz);
            tpOp->replaceLocalValues(i,cols,valview);
          }
        } //for (size_t i=0; ...
      } else {
        for (size_t i=0; i<rowMap->getNodeNumElements(); ++i) {
          GO gid = rowMap->getGlobalElement(i);
          tpOp->getGlobalRowView(gid,cols,vals);
          size_t nnz = tpOp->getNumEntriesInGlobalRow(gid);
          if (nnz>maxRowSize) {
            maxRowSize=nnz;
            scaledVals.resize(maxRowSize);
          }
          for (size_t j=0; j<nnz; ++j) {
            scaledVals[j] = vals[j]*sv[i]; //FIXME i or gid?
          }
          if (nnz>0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0],nnz);
            tpOp->replaceGlobalValues(gid,cols,valview);
          }
        } //for (size_t i=0; ...
      }

      if (doFillComplete) {
        if (domainMap == Teuchos::null || rangeMap == Teuchos::null)
          throw(Exceptions::RuntimeError("In Utils::Scaling: cannot fillComplete because the domain and/or range map hasn't been defined"));
        if (doOptimizeStorage)
          tpOp->fillComplete(domainMap,rangeMap,Tpetra::DoOptimizeStorage);
        else
          tpOp->fillComplete(domainMap,rangeMap,Tpetra::DoNotOptimizeStorage);
      }

   } //ScaleMatrix()

  }; // class

} //namespace MueLu
#define MUELU_UTILITIES_SHORT

#endif //ifndef MUELU_UTILITIES_HPP
