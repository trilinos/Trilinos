#ifndef MUELU_TRANSPFACTORY_DEF_HPP
#define MUELU_TRANSPFACTORY_DEF_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>

#include <Xpetra_Operator_fwd.hpp>

#include "MueLu_TransPFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TransPFactory(RCP<FactoryBase> PFact)
    : PFact_(PFact)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TransPFactory() {}
 
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    coarseLevel.DeclareInput("P", PFact_.get(), this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const {
    return BuildR(fineLevel,coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildR(Level & fineLevel, Level & coarseLevel) const {

    std::ostringstream buf; buf << coarseLevel.GetLevelID();
    RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("TransPFactory::OldBuildR_"+buf.str()));
    timer->start(true);

    Teuchos::OSTab tab(this->getOStream());
    Teuchos::ParameterList matrixList;
    RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", PFact_.get());

    //doesn't work -- bug in EpetraExt?
    //RCP<Operator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>("Identity",P->getRangeMap(),matrixList);
    //      RCP<CrsOperator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>("Identity",P->getDomainMap(),matrixList);
    //RCP<Operator> R = Utils::TwoMatrixMultiply(P,true,I,false); //doesn't work -- bug in EpetraExt?
    //      RCP<Operator> R = Utils::TwoMatrixMultiply(I,false,P,true);

    RCP<Operator> R = Utils2::Transpose(P,true);

    coarseLevel.Set("R", R, this);

    timer->stop();
    MemUtils::ReportTimeAndMemory(*timer, *(P->getRowMap()->getComm()));

  } //BuildR

} //namespace MueLu

#endif // MUELU_TRANSPFACTORY_DEF_HPP
