#ifndef MUELU_EMINPFACTORY_DEF_HPP
#define MUELU_EMINPFACTORY_DEF_HPP

#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_EminPFactory_decl.hpp"

#include "MueLu_Utilities.hpp"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_PatternFactory.hpp"

#include "MueLu_SteepestDescentSolver.hpp"
#include "MueLu_CGSolver.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Constraint_fwd.hpp"

#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",                 Teuchos::null, "Generating factory for the matrix A used during internal iterations");
    validParamList->set< RCP<const FactoryBase> >("P",                 Teuchos::null, "Generating factory for the initial guess");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",         Teuchos::null, "Generating factory for the nullspace");
    validParamList->set< RCP<const FactoryBase> >("Constraint",        Teuchos::null, "Generating factory for constraints");
    validParamList->set< int >                   ("Niterations",                   3, "Number of iterations of the internal iterative method");
    validParamList->set< int >                   ("Reuse Niterations",             1, "Number of iterations of the internal iterative method");
    validParamList->set< RCP<Matrix> >           ("P0",                Teuchos::null, "Initial guess at P");
    validParamList->set< bool >                  ("Keep P0",                   false, "Keep an initial P0 (for reuse)");
    validParamList->set< RCP<Constraint> >       ("Constraint0",       Teuchos::null, "Initial Constraint");
    validParamList->set< bool >                  ("Keep Constraint0",          false, "Keep an initial Constraint (for reuse)");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel,   "A");
    Input(fineLevel,   "Nullspace");
    // FIXME: we should not request for P or Constraint if we reuse stuff
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

    const ParameterList & pL = GetParameterList();

    // Set keep flags
    if (pL.isParameter("Keep P0") && pL.get<bool>("Keep P0"))
      coarseLevel.Keep("P0",this);
    if (pL.isParameter("Keep Constraint0") && pL.get<bool>("Keep Constraint0"))
      coarseLevel.Keep("Constraint0",this);

    // Reuse
    int Niterations;

    // Get A, B
    RCP<Matrix>      A = Get< RCP<Matrix> >     (fineLevel,   "A");
    RCP<MultiVector> B = Get< RCP<MultiVector> >(fineLevel,   "Nullspace");

    // Get P0 or make P
    RCP<Matrix>      P0;
    if (coarseLevel.IsAvailable("P0", this)) {
      P0          = coarseLevel.Get<RCP<Matrix> >("P0", this);
      Niterations = pL.get<int>("Reuse Niterations");
      GetOStream(Runtime0, 0) << "EminPFactory: Reusing P0"<<std::endl;

    } else {
      P0          = Get< RCP<Matrix> >(coarseLevel, "P");
      Niterations = pL.get<int>("Niterations");
    }

    // Get Constraint0 or make Constraint
    RCP<Constraint> X;
    if (coarseLevel.IsAvailable("Constraint0", this)) {
      X = coarseLevel.Get<RCP<Constraint> >("Constraint0", this);
      GetOStream(Runtime0, 0) << "EminPFactory: Reusing Constraint0"<<std::endl;

    } else {
      X = Get< RCP<Constraint> > (coarseLevel, "Constraint");
    }


    RCP<Matrix> P;
    CGSolver EminSolver(Niterations);
    EminSolver.Iterate(*A, *X, *P0, *B, P);

    Set(coarseLevel, "Constraint0", X);
    Set(coarseLevel, "P",           P);
    Set(coarseLevel, "P0",          P);

    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    GetOStream(Statistics0,0) << Utils::PrintMatrixInfo(*P, "P", params);
  }

} // namespace MueLu

#endif // MUELU_EMINPFACTORY_DEF_HPP
