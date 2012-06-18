#ifndef MUELU_SCHURCOMPLEMENTFACTORY_DEF_HPP_
#define MUELU_SCHURCOMPLEMENTFACTORY_DEF_HPP_

#include <Xpetra_BlockedCrsOperator.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_OperatorFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_HierarchyHelpers.hpp"

#include "MueLu_SchurComplementFactory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SchurComplementFactory(Teuchos::RCP<const FactoryBase> Afact, Scalar omega/**size_t row, size_t col, LocalOrdinal blksize*/)
: AFact_(Afact), omega_(omega)/*, row_(row), col_(col), blksize_(blksize)*/
  { }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~SchurComplementFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
  currentLevel.DeclareInput("A",AFact_.get(),this);
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
  FactManager_ = FactManager;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const
{
  FactoryMonitor  m(*this, "Building Schur Complement", currentLevel);
  Teuchos::RCP<Operator> A = currentLevel.Get<RCP<Operator> >("A", AFact_.get());


  RCP<BlockedCrsOperator> bA = Teuchos::rcp_dynamic_cast<BlockedCrsOperator>(A);
  TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::BlockedPFactory::Build: input matrix A is not of type BlockedCrsOperator! error.");

  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixClass;
  typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOperatorClass;

  Teuchos::RCP<CrsMatrixClass> A00 = bA->getMatrix(0,0);
  Teuchos::RCP<CrsMatrixClass> A01 = bA->getMatrix(0,1);
  Teuchos::RCP<CrsMatrixClass> A10 = bA->getMatrix(1,0);
  Teuchos::RCP<CrsMatrixClass> A11 = bA->getMatrix(1,1);

  Teuchos::RCP<CrsOperatorClass> Op00 = Teuchos::rcp(new CrsOperatorClass(A00));
  Teuchos::RCP<CrsOperatorClass> Op01 = Teuchos::rcp(new CrsOperatorClass(A01));
  Teuchos::RCP<CrsOperatorClass> Op10 = Teuchos::rcp(new CrsOperatorClass(A10));
  Teuchos::RCP<CrsOperatorClass> Op11 = Teuchos::rcp(new CrsOperatorClass(A11));

  Teuchos::RCP<Operator> F = Teuchos::rcp_dynamic_cast<Operator>(Op00);
  Teuchos::RCP<Operator> G = Teuchos::rcp_dynamic_cast<Operator>(Op01);
  Teuchos::RCP<Operator> D = Teuchos::rcp_dynamic_cast<Operator>(Op10);
  Teuchos::RCP<Operator> Z = Teuchos::rcp_dynamic_cast<Operator>(Op11);


  //*************
  //Version with MyOldScaleMatrix
  //DA: in MueLu::Utils there is a left scale method but it works with an ArrayRCP
  //DA: This method applies the inverse by default and get the data from the Operator, therefore the reciprocal step is not needed.

  Teuchos::ArrayRCP<SC> AdiagFinv = Utils::GetMatrixDiagonal(F);

  //copy the value of G so we can do the left scale. To do so, a new Operator is built and them the original values of G are added.
  RCP<Operator> FhatinvG = OperatorFactory::Build(G->getRowMap(), G->getGlobalMaxNumRowEntries());
  Utils2::TwoMatrixAdd(G, false, 1.0, FhatinvG, 0.0);
  //  this->GetOStream(Runtime0, 0) << "======= G after build =======" << std::endl;
  //  Utils::MatrixPrint(G_);
  //DA: This step is necessary, otherwise the Left scale does not work. why?
  FhatinvG->fillComplete(G->getDomainMap(),G->getRowMap());
  //  this->GetOStream(Runtime0, 0) << "FhatinvG after build" << std::endl;
  //  Utils::MatrixPrint(FhatinvG);
  // Here the first boolean is true, it means we multiply by the inverse of the diagonal of F_
  Utils::MyOldScaleMatrix(FhatinvG,AdiagFinv,true,false,false);
  //  this->GetOStream(Runtime0, 0) << "======= FhatinvG =======" << std::endl;
  //  Utils::MatrixPrint(FhatinvG);

  RCP<Operator> DFhatinvG = Utils::TwoMatrixMultiply(D,false,FhatinvG,false);
  //  this->GetOStream(Runtime0, 0) << "======= DFhatinvG =======" << std::endl;
  //  Utils::MatrixPrint(DFhatinvG);
  RCP<Operator> S;
  Utils2::TwoMatrixAdd(Z,false,1.0,DFhatinvG,false,-1.0/omega_,S);
  S->fillComplete();
  //  this->GetOStream(Runtime0, 0) << "======= S =======" << std::endl;
  //  Utils::MatrixPrint(S_);
  //  this->GetOStream(Runtime0, 0) << "Omega = " << omega_<< std::endl;

  {
    //The Schur complement must be set in the new factory manager
    //this->GetOStream(Runtime0, 0) << "SchurComplementFactoy::Build: Setting Schur Complement";
    //DA: New version: Instead of using Factory Managers, try to set S
    //SetFactoryManager currentSFM  (rcpFromRef(currentLevel),   FactManager_);


    currentLevel.Set("A", S, this);
  }
  //this->GetOStream(Runtime0, 0) << " ... Done." << std::endl;
}

} // namespace MueLu

#endif /* MUELU_SCHURCOMPLEMENTFACTORY_DEF_HPP_ */
