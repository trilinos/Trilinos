// New definition of types using the types ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.

#include <Cthulhu_UseShortNamesScalar.hpp>

#ifdef MUELU_LEVEL_SHORT
typedef MueLu::Level<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>               Level;
#endif

#ifdef MUELU_HIERARCHY_SHORT
typedef MueLu::Hierarchy<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>           Hierarchy;
#endif

#ifdef MUELU_SAPFACTORY_SHORT
typedef MueLu::SaPFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      SaPFactory;
#endif

#ifdef MUELU_PRFACTORY_SHORT
typedef MueLu::PRFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      PRFactory;
#endif

#ifdef MUELU_GENERICPRFACTORY_SHORT
typedef MueLu::GenericPRFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> GenericPRFactory;
#endif

#ifdef MUELU_PFACTORY_SHORT
typedef MueLu::PFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      PFactory;
#endif

#ifdef MUELU_RFACTORY_SHORT
typedef MueLu::RFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      RFactory;
#endif

#ifdef MUELU_TRANSPFACTORY_SHORT
typedef MueLu::TransPFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>   TransPFactory;
#endif

#ifdef MUELU_RAPFACTORY_SHORT
typedef MueLu::RAPFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>      RAPFactory;
#endif

#ifdef MUELU_SMOOTHERPROTOTYPE_SHORT
typedef MueLu::SmootherPrototype<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherPrototype;
#endif

#ifdef MUELU_SMOOTHERBASE_SHORT
typedef MueLu::SmootherBase<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>     SmootherBase;
#endif

#ifdef MUELU_SMOOTHERFACTORYBASE_SHORT
typedef MueLu::SmootherFactoryBase<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherFactoryBase;
#endif

#ifdef MUELU_SMOOTHERFACTORY_SHORT
typedef MueLu::SmootherFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherFactory;
#endif

#ifdef MUELU_TENTATIVEPFACTORY_SHORT
typedef MueLu::TentativePFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TentativePFactory;
#endif

#ifdef MUELU_OPERATORFACTORY_SHORT
typedef MueLu::OperatorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OperatorFactory;
#endif

#ifdef MUELU_SMOOTHER_SHORT
typedef MueLu::Smoother<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>        Smoother;
#endif

#ifdef MUELU_SALEVEL_SHORT
typedef MueLu::SaLevel<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>               SaLevel;
#endif

#ifdef MUELU_UTILITIES_SHORT
typedef MueLu::Utils<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>               Utils;
#endif

#ifdef MUELU_GAUSSSEIDEL_SHORT
typedef MueLu::GaussSeidel<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>          GaussSeidel;
#endif

#ifdef MUELU_IFPACK_SMOOTHER_SHORT
typedef MueLu::IfpackSmoother<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       IfpackSmoother;
#endif

#ifdef MUELU_AMESOS_SMOOTHER_SHORT
typedef MueLu::AmesosSmoother<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       AmesosSmoother;
#endif

#ifdef MUELU_UCAGGREGATIONCOMMHELPER_SHORT
typedef MueLu::UCAggregationCommHelper<ScalarType, LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>       UCAggregationCommHelper;
#endif
