// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.

#include <Xpetra_UseShortNamesScalar.hpp>

#ifdef MUELU_HIERARCHY_SHORT
typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>           Hierarchy;
#endif

#ifdef MUELU_SAPFACTORY_SHORT
typedef MueLu::SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      SaPFactory;
#endif

#ifdef MUELU_PRFACTORY_SHORT
typedef MueLu::PRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      PRFactory;
#endif

#ifdef MUELU_GENERICPRFACTORY_SHORT
typedef MueLu::GenericPRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> GenericPRFactory;
#endif

#ifdef MUELU_PFACTORY_SHORT
typedef MueLu::PFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      PFactory;
#endif

#ifdef MUELU_RFACTORY_SHORT
typedef MueLu::RFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      RFactory;
#endif

#ifdef MUELU_TRANSPFACTORY_SHORT
typedef MueLu::TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>   TransPFactory;
#endif

#ifdef MUELU_RAPFACTORY_SHORT
typedef MueLu::RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      RAPFactory;
#endif

#ifdef MUELU_SMOOTHERPROTOTYPE_SHORT
typedef MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherPrototype;
#endif

#ifdef MUELU_SMOOTHERBASE_SHORT
typedef MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>     SmootherBase;
#endif

#ifdef MUELU_SMOOTHERFACTORYBASE_SHORT
typedef MueLu::SmootherFactoryBase<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherFactoryBase;
#endif

#ifdef MUELU_SMOOTHERFACTORY_SHORT
typedef MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherFactory;
#endif

#ifdef MUELU_TENTATIVEPFACTORY_SHORT
typedef MueLu::TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TentativePFactory;
#endif

#ifdef MUELU_TWOLEVELFACTORY_SHORT
typedef MueLu::TwoLevelFactoryBase<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TwoLevelFactoryBase;
#endif

#ifdef MUELU_SMOOTHER_SHORT
typedef MueLu::Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>        Smoother;
#endif

#ifdef MUELU_SALEVEL_SHORT
typedef MueLu::SaLevel<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>               SaLevel;
#endif

#ifdef MUELU_UTILITIES_SHORT
typedef MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>               Utils;
#endif

#ifdef MUELU_GAUSSSEIDEL_SHORT
typedef MueLu::GaussSeidel<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>          GaussSeidel;
#endif

#ifdef MUELU_IFPACK_SMOOTHER_SHORT
typedef MueLu::IfpackSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       IfpackSmoother;
#endif

#ifdef MUELU_IFPACK2_SMOOTHER_SHORT
typedef MueLu::Ifpack2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       Ifpack2Smoother;
#endif

#ifdef MUELU_AMESOS_SMOOTHER_SHORT
typedef MueLu::AmesosSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       AmesosSmoother;
#endif

#ifdef MUELU_AMESOS2_SMOOTHER_SHORT
typedef MueLu::Amesos2Smoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       Amesos2Smoother;
#endif

#ifdef MUELU_MERGED_SMOOTHER_SHORT
typedef MueLu::MergedSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       MergedSmoother;
#endif

#ifdef MUELU_UCAGGREGATIONCOMMHELPER_SHORT
typedef MueLu::UCAggregationCommHelper<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       UCAggregationCommHelper;
#endif

#ifdef MUELU_COALESCEDROPFACTORY_SHORT
typedef MueLu::CoalesceDropFactory<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CoalesceDropFactory;
#endif

#ifdef MUELU_UCAGGREGATIONFACTORY_SHORT
typedef MueLu::UCAggregationFactory<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       UCAggregationFactory;
#endif

#ifdef MUELU_LOCALAGGREGATIONFACTORY_SHORT
typedef MueLu::LocalAggregationFactory<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       LocalAggregationFactory;
#endif

#ifdef MUELU_DEFAULTFACTORYHANDLER_SHORT
typedef MueLu::DefaultFactoryHandler<Scalar, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       DefaultFactoryHandler;
#endif

