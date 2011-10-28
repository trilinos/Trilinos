#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION

#include "MueLu_LeftoverAggregationAlgorithm_def.hpp"

template class MueLu::LeftoverAggregationAlgorithm<int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::LeftoverAggregationAlgorithm<int, long long int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;
#endif // HAVE_TEUCHOS_LONG_LONG_INT

#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
