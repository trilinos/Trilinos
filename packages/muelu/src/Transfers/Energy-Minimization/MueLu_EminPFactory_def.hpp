#ifndef MUELU_EMINPFACTORY_DEF_HPP
#define MUELU_EMINPFACTORY_DEF_HPP

#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_EminPFactory_decl.hpp"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_PatternFactory.hpp"

#include "MueLu_SteepestDescentSolver.hpp"
#include "MueLu_CGSolver.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Constraint_fwd.hpp"

#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::EminPFactory(RCP<const FactoryBase> InitialPFact, RCP<const FactoryBase> ConstrFact, RCP<const FactoryBase> AFact)
  : initialPFact_(InitialPFact), constrFact_(ConstrFact), AFact_(AFact)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~EminPFactory()
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    RCP<const FactoryBase> initialPFact = initialPFact_;
    if (initialPFact == Teuchos::null)
      initialPFact  = coarseLevel.GetFactoryManager()->GetFactory("Ptent");

    fineLevel.  DeclareInput("A",           AFact_.get(), this);
    fineLevel.  DeclareInput("Nullspace",   fineLevel.GetFactoryManager()->GetFactory("Nullspace").get(), this);
    coarseLevel.DeclareInput("P",           initialPFact.get(), this);
    coarseLevel.DeclareInput("Constraint",  constrFact_.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level& coarseLevel) const {
    BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Prolongator minimization", coarseLevel);

    RCP<Matrix>      A          = fineLevel.  Get< RCP<Matrix> >("A", AFact_.get());
    RCP<MultiVector> B          = fineLevel.  Get< RCP<MultiVector> >("Nullspace", fineLevel.GetFactoryManager()->GetFactory("Nullspace").get());
    RCP<Matrix>      P0         = coarseLevel.Get< RCP<Matrix> >("P", initialPFact_.get());
    RCP<Constraint>  X          = coarseLevel.Get< RCP<Constraint> >("Constraint", constrFact_.get());

    RCP<Matrix>      P;
    // SteepestDescentSolver EminSolver(1, 1.0);
    CGSolver EminSolver(3);
    EminSolver.Iterate(*A, *X, *P0, *B, P);

    coarseLevel.Set("P", P, this);
  }

} // namespace MueLu

#endif // MUELU_EMINPFACTORY_DEF_HPP
