#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_LocalAggregationAlgorithm_def.hpp"

#ifdef HAVE_MUELU_EI_DOUBLE_INT_INT
template class MueLu::LocalAggregationAlgorithm<int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;
#else
#error
#endif

#ifdef HAVE_MUELU_EI_DOUBLE_INT_LONGLONGINT
# ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::LocalAggregationAlgorithm<int, long long int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;
#else
# warning To compile MueLu with 'long long int' support, please turn on HAVE_TEUCHOS_LONG_LONG_INT
# endif
#endif
