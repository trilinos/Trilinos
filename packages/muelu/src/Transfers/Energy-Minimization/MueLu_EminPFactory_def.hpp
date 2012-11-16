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
  EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::EminPFactory()
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~EminPFactory()
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel,   "A");
    Input(fineLevel,   "Nullspace");
    Input(coarseLevel, "P");
    Input(coarseLevel, "Constraint");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level& coarseLevel) const {
    BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Prolongator minimization", coarseLevel);

    RCP<Matrix>      A          = Get< RCP<Matrix> >     (fineLevel,   "A");
    RCP<MultiVector> B          = Get< RCP<MultiVector> >(fineLevel,   "Nullspace");
    RCP<Matrix>      P0         = Get< RCP<Matrix> >     (coarseLevel, "P");
    RCP<Constraint>  X          = Get< RCP<Constraint> > (coarseLevel, "Constraint");

    RCP<Matrix>      P;
    CGSolver EminSolver(3);
    EminSolver.Iterate(*A, *X, *P0, *B, P);

    Set(coarseLevel, "P", P);
  }

} // namespace MueLu

#endif // MUELU_EMINPFACTORY_DEF_HPP
