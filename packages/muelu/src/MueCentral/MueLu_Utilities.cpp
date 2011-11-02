#include <MueLu_Utilities_def.hpp>

namespace MueLu {

  RCP<Xpetra::Operator<double, int, int> > Utils2<double, int, int>::Transpose(RCP<Xpetra::Operator<double, int, int> > const &Op, bool const & optimizeTranspose)
  {
   typedef Xpetra::Operator<double,int,int> Operator;
   typedef double Scalar;
   typedef int LocalOrdinal;
   typedef int GlobalOrdinal;
   typedef Kokkos::DefaultNode::DefaultNodeType Node;
   typedef Kokkos::DefaultKernels<double,int,NO>::SparseOps LocalMatOps;

#ifdef HAVE_MUELU_EPETRAEXT
    std::string TorE = "epetra";
#else
    std::string TorE = "tpetra";
#endif

#ifdef HAVE_MUELU_EPETRAEXT
    RCP<Epetra_CrsMatrix> epetraOp;
    try {
      epetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstEpetraCrs(Op);
    }
    catch (...) {
      TorE = "tpetra";
    }
#endif

#ifdef HAVE_MUELU_TPETRA
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpetraOp;
    if (TorE=="tpetra") {
      try {
        tpetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2TpetraCrs(Op);
      }
      catch (...) {
        throw(Exceptions::RuntimeError("Utils::Transpose: Can only transpose Crs matrices"));
      }
    } //if
#endif

    if (TorE == "tpetra") {
#ifdef HAVE_MUELU_TPETRA
      //     Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> transposer(*tpetraOp); //more than meets the eye
      //     RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A = transposer.createTranspose(optimizeTranspose ? Tpetra::DoOptimizeStorage : Tpetra::DoNotOptimizeStorage); //couldn't have just used a bool...
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A = Utils<SC,LO,GO>::simple_Transpose(tpetraOp);
      RCP<Xpetra::TpetraCrsMatrix<SC> > AA = rcp(new Xpetra::TpetraCrsMatrix<SC>(A) );
      RCP<Xpetra::CrsMatrix<SC> > AAA = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
      RCP<Xpetra::CrsOperator<SC> > AAAA = rcp( new Xpetra::CrsOperator<SC> (AAA) );
      AAAA->fillComplete(Op->getRangeMap(),Op->getDomainMap());
      return AAAA;
#else
      throw(Exceptions::RuntimeError("Tpetra"));
#endif
    } else {
#ifdef HAVE_MUELU_EPETRAEXT
      //epetra case
      Epetra_RowMatrixTransposer et(&*epetraOp);
      Epetra_CrsMatrix *A;
      int rv = et.CreateTranspose(false,A);
      if (rv != 0) {
        std::ostringstream buf;
        buf << rv;
        std::string msg = "Utils::Transpose: Epetra::RowMatrixTransposer returned value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }

      RCP<Epetra_CrsMatrix> rcpA(A);
      //RCP<Epetra_CrsMatrix> rcpA = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::simple_EpetraTranspose(epetraOp);
      RCP<EpetraCrsMatrix> AA = rcp(new EpetraCrsMatrix(rcpA) );
      RCP<Xpetra::CrsMatrix<SC> > AAA = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
      RCP<Xpetra::CrsOperator<SC> > AAAA = rcp( new Xpetra::CrsOperator<SC>(AAA) );
      AAAA->fillComplete(Op->getRangeMap(),Op->getDomainMap());
      return AAAA;
#else
      throw(Exceptions::RuntimeError("Epetra (Err. 2)"));
#endif
    }
     
  } //Transpose

}
