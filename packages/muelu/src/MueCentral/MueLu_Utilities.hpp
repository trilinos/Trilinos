#ifndef MUELU_UTILITIES_HPP
#define MUELU_UTILITIES_HPP

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Comm.hpp>

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
      @param B right matrix
      @param transposeA if true, use the transpose of A
      @param transposeB if true, use the transpose of B
    */
   static RCP<Operator> TwoMatrixMultiply(RCP<Operator> const &A, RCP<Operator> const &B,
                                          bool transposeA=false, bool transposeB=false) //TODO: modify definition to respect definition of Epetra/Tpetra::MatrixMatrix::Multiply (order of input args)
    {
      //FIXME 30 is likely a big overestimate
      RCP<Operator> C = OperatorFactory::Build(A->getRowMap(), 30);

      if (!A->isFillComplete())
        throw(Exceptions::RuntimeError("A is not fill-completed"));
      if (!B->isFillComplete())
        throw(Exceptions::RuntimeError("B is not fill-completed"));

      if (C->getRowMap()->lib() == Cthulhu::UseEpetra) {
#ifdef HAVE_MUELU_EPETRAEXT
        RCP<const Epetra_CrsMatrix> epA = Op2EpetraCrs(A);
        RCP<const Epetra_CrsMatrix> epB = Op2EpetraCrs(B);
        RCP<Epetra_CrsMatrix>       epC = Op2NonConstEpetraCrs(C);
        
        int i = EpetraExt::MatrixMatrix::Multiply(*epA,transposeA,*epB,transposeB,*epC);
        
        if (i != 0) {
          std::ostringstream buf;
          buf << i;
        std::string msg = "EpetraExt::MatrixMatrix::Multiply return value of " + buf.str();
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
        
        Tpetra::MatrixMatrix::Multiply(*tpA,transposeA,*tpB,transposeB,*tpC);
#else
        throw(Exceptions::RuntimeError("MueLu must be compile with Tpetra."));
#endif
      }

      return C;
    } //TwoMatrixMultiply()

    /*! @brief Helper function to calculate alpha*A + beta*B.

      @param A (required) left matrix operand
      @param B (required) right matrix operand
      @param alpha (optional) scalar multiplier for A, defaults to 1.0
      @param beta  (optional) scalar multiplier for B, defaults to 1.0

      @return Teuchos::RCP to non-constant Cthulhu::Operator.

      Note that the returned matrix is fill-complete'd.
    */
   static RCP<Operator> TwoMatrixAdd(RCP<Operator> const &A, RCP<Operator> const &B, SC alpha=1.0, SC beta=1.0)
    {
      //FIXME 30 is a complete guess as to the #nonzeros per row
      RCP<Operator> C = rcp( new CrsOperator(A->getRowMap(), 30) );

      if (C->getRowMap()->lib() == Cthulhu::UseEpetra) {
#ifdef HAVE_MUELU_EPETRAEXT
        RCP<const Epetra_CrsMatrix> epA = Op2EpetraCrs(A);
        RCP<const Epetra_CrsMatrix> epB = Op2EpetraCrs(B);
        RCP<Epetra_CrsMatrix>       epC = Op2NonConstEpetraCrs(C);
        Epetra_CrsMatrix* ref2epC = &*epC; //to avoid a compiler error...
        
        //FIXME is there a bug if beta=0?
        int i = EpetraExt::MatrixMatrix::Add(*epA,false,(double)alpha,*epB,false,(double)beta,ref2epC);

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
        
        Tpetra::MatrixMatrix::Add(*tpA, false, alpha, *tpB, false, beta, tpC);
#else
        throw(Exceptions::RuntimeError("MueLu must be compile with Tpetra."));
#endif
      }

      C->fillComplete(A->getDomainMap(),A->getRangeMap()); //TODO: should be done outside or well documented
      return C;
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
      std::vector<LO> diagInd(1);
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
     RCP<MultiVector> RES = MultiVectorFactory::Build(Op.getRowMap(),numVecs);
     Op.multiply(X,*RES,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
     RES->update(one,RHS,negone);
     return RES;
   }

#include <unistd.h>
   /*@brief Utility for pausing execution to attach debugger.
   
     WARNING -- This utility checks for the existence of a file.  This could be *very* expensive
                if there are many simultaneous processes doing this.
   */
   static void BreakForDebugger(const Teuchos::Comm<GO> &Comm)
   {
     // Pulled directly from ML.
     // print out some junk for debugging
     // LAM/MPI has some difficulties related to environmental variables.
     // The problem is that LAM/MPI uses ssh to log in, and on some
     // machine "export ML_BREAK_FOR_DEBUGGER" does not work. So, we
     // allow two options:
     // 1.) export MUELU_BREAK_FOR_DEBUGGER=1
     // 2.) create a file in the executable directory, called MueLu_debug_now
   
     char * str = (char *) getenv("MUELU_BREAK_FOR_DEBUGGER");
     GO i = 0, j = 0;
     char buf[80];
     char go = ' ';
     char hostname[80];
     if (str != NULL) i++;
   
     int mypid = Comm.getRank();
     int root = 0;
     char buff;
     FILE * MueLu_capture_flag=0;

     if (mypid == 0) {
       MueLu_capture_flag = fopen("MueLu_debug_now","r");
       if (MueLu_capture_flag) {
         buff = 'T';
         fclose(MueLu_capture_flag);
       } else 
         buff = 'F';
     }
     //Comm.broadcast(root,1,buff); //FIXME

     if (buff == 'T')
     {
       if (mypid  == 0) std::cout << "Host and Process Ids for tasks" << std::endl;
       for (i = 0; i <Comm.getSize() ; i++) {
         if (i == mypid ) {
       gethostname(hostname, sizeof(hostname));
       LO pid = getpid();
       sprintf(buf, "Host: %s\tMPI rank: %d\tPID: %d\n\tattach %d\n\tcontinue\n",
           hostname, mypid, pid, pid);
       printf("%s\n",buf);
       fflush(stdout);
       sleep(1);
         }
       }
        if(mypid == 0) {
          printf("\n");
          printf("** Pausing because environment variable MUELU_BREAK_FOR_DEBUGGER has been set,\n");
          puts("** or file MueLu_debug_now has been created");
          printf("**\n");
          printf("** You may now attach debugger to the processes listed above.\n");
          printf( "**\n");
          printf( "** Enter a character to continue > "); fflush(stdout);
          scanf("%c",&go);
        }
        Comm.barrier();
      }
   
   } //BreakForDebugger()


   /*! @brief Save matrix to file in Matrix Market format.

     TODO Move this to Cthulhu?
   */
   static void Write(std::string const & fileName, RCP<Operator> const & Op) {
    RCP<const CrsOperator> crsOp = rcp_dynamic_cast<const CrsOperator>(Op);
    if (crsOp == Teuchos::null)
      throw(Exceptions::BadCast("Only Cthulhu::CrsOperator is supported in matrix writing"));
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
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

  }; // class

} //namespace MueLu

#define MUELU_UTILITIES_SHORT

#endif //ifndef MUELU_UTILITIES_HPP
