
// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"
#include "Kokkos_CuspSparseOps.cuh"

#define INSTANTIATE_CUSPSPARSEGRAPH_ORDINAL(ORDINAL) \
        template class Kokkos::CrsGraph< ORDINAL, Kokkos::ThrustGPUNode, Kokkos::CuspSparseOps<void,ORDINAL,Kokkos::ThrustGPUNode> >;

#define INSTANTIATE_CUSPSPARSEOPS_ORDINAL_SCALAR(ORDINAL,SCALAR) \
        template class Kokkos::CrsMatrix< SCALAR, ORDINAL, Kokkos::ThrustGPUNode, Kokkos::CuspSparseOps<void,ORDINAL,Kokkos::ThrustGPUNode> >; \
        template class Kokkos::CuspSparseOps<SCALAR,ORDINAL,Kokkos::ThrustGPUNode>;

typedef short int ShortInt;

INSTANTIATE_CUSPSPARSEGRAPH_ORDINAL(int)
INSTANTIATE_CUSPSPARSEGRAPH_ORDINAL(ShortInt)

#ifdef HAVE_KOKKOS_CUDA_FLOAT
INSTANTIATE_CUSPSPARSEOPS_ORDINAL_SCALAR(int,float)
INSTANTIATE_CUSPSPARSEOPS_ORDINAL_SCALAR(ShortInt,float)
#endif

#ifdef HAVE_KOKKOS_CUDA_DOUBLE
INSTANTIATE_CUSPSPARSEOPS_ORDINAL_SCALAR(int,double)
INSTANTIATE_CUSPSPARSEOPS_ORDINAL_SCALAR(ShortInt,double)
#endif

INSTANTIATE_CUSPSPARSEOPS_ORDINAL_SCALAR(int,int)
INSTANTIATE_CUSPSPARSEOPS_ORDINAL_SCALAR(ShortInt,int)
