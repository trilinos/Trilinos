// Type definitions for templated classes (generally graph-related) that do not require a scalar.

#include <Cthulhu_UseShortNamesOrdinal.hpp>

#ifdef MUELU_GRAPH_SHORT
typedef MueLu::Graph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>                          Graph;
#endif

#ifdef MUELU_AGGREGATES_SHORT
typedef MueLu::Aggregates<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>               Aggregates;
#endif
