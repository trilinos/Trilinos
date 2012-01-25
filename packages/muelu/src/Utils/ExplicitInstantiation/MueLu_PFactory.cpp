// #include "MueLu_ExplicitInstantiation.hpp"

// #include "MueLu_PFactory_def.hpp"

// #ifdef HAVE_MUELU_INST_DOUBLE_INT_INT
// template class MueLu::PFactory<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;
// #else
// #error
// #endif

// #ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
// # ifdef HAVE_TEUCHOS_LONG_LONG_INT
// template class MueLu::PFactory<double, int, long long int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;
// #else
// # warning To compile MueLu with 'long long int' support, please turn on HAVE_TEUCHOS_LONG_LONG_INT
// # endif
// #endif
