// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.

#include <Xpetra_UseShortNamesScalar.hpp>

#ifdef MUELU_HIERARCHY_SHORT
typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Hierarchy;
#endif

#ifdef MUELU_SAPFACTORY_SHORT
typedef MueLu::SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SaPFactory;
#endif

#ifdef MUELU_PGPFACTORY_SHORT
typedef MueLu::PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> PgPFactory;
#endif

/*#ifdef MUELU_GENERICPRFACTORY_SHORT
typedef MueLu::GenericPRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> GenericPRFactory;
#endif*/

#ifdef MUELU_GENERICRFACTORY_SHORT
typedef MueLu::GenericRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> GenericRFactory;
#endif

#ifdef MUELU_TRANSPFACTORY_SHORT
typedef MueLu::TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TransPFactory;
#endif

#ifdef MUELU_RAPFACTORY_SHORT
typedef MueLu::RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> RAPFactory;
#endif

#ifdef MUELU_SMOOTHERPROTOTYPE_SHORT
typedef MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherPrototype;
#endif

#ifdef MUELU_FAKESMOOTHERPROTOTYPE_SHORT
typedef MueLu::FakeSmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> FakeSmootherPrototype;
#endif

#ifdef MUELU_SMOOTHERBASE_SHORT
typedef MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherBase;
#endif

#ifdef MUELU_SMOOTHERFACTORY_SHORT
typedef MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherFactory;
#endif

#ifdef MUELU_TENTATIVEPFACTORY_SHORT
typedef MueLu::TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TentativePFactory;
#endif

#ifdef MUELU_SMOOTHER_SHORT
typedef MueLu::Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Smoother;
#endif

#ifdef MUELU_UTILITIES_SHORT
typedef MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Utils;
#endif

#ifdef MUELU_GAUSSSEIDELSMOOTHER_SHORT
typedef MueLu::GaussSeidelSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> GaussSeidelSmoother;
#endif

#ifdef MUELU_IFPACK2_SMOOTHER_SHORT
typedef MueLu::Ifpack2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Ifpack2Smoother;
#endif

#ifdef MUELU_AMESOS2_SMOOTHER_SHORT
typedef MueLu::Amesos2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Amesos2Smoother;
#endif

#ifdef MUELU_DIRECT_SOLVER_SHORT
typedef MueLu::DirectSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> DirectSolver;
#endif

#ifdef MUELU_TRILINOS_SMOOTHER_SHORT
typedef MueLu::TrilinosSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> TrilinosSmoother;
#endif

#ifdef MUELU_MERGED_SMOOTHER_SHORT
typedef MueLu::MergedSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> MergedSmoother;
#endif

#ifdef MUELU_COALESCEDROPFACTORY_SHORT
typedef MueLu::CoalesceDropFactory<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CoalesceDropFactory;
#endif

#ifdef MUELU_PREDROPFUNCTIONBASECLASS_SHORT
typedef MueLu::PreDropFunctionBaseClass<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PreDropFunctionBaseClass;
#endif

#ifdef MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
typedef MueLu::PreDropFunctionConstVal<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PreDropFunctionConstVal;
#endif

/*#ifdef MUELU_DEFAULTFACTORYHANDLER_SHORT
typedef MueLu::DefaultFactoryHandler<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> DefaultFactoryHandler;
#endif*/

#ifdef MUELU_NULLSPACEFACTORY_SHORT
typedef MueLu::NullspaceFactory<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> NullspaceFactory;
#endif

#ifdef MUELU_FACTORYMANAGER_SHORT
typedef MueLu::FactoryManager<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> FactoryManager;
#endif

#ifdef MUELU_THRESHOLDAFILTERFACTORY_SHORT
typedef MueLu::ThresholdAFilterFactory<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ThresholdAFilterFactory;
#endif
