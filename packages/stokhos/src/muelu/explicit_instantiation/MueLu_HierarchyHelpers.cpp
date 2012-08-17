#include "MueLu_ExplicitInstantiation.hpp"
#include "Stokhos_ConfigDefs.h"

#if defined(HAVE_STOKHOS_MUELU) && defined(HAVE_MUELU_EXPLICIT_INSTANTIATION) && defined(HAVE_STOKHOS_SACADO)

// Sacado headers must be included first so that overloaded operators
// are defined in the muelu template code
#include "Stokhos_Sacado.hpp"
#include "MueLu_HierarchyHelpers_def.hpp"

typedef Stokhos::StandardStorage<int,double> Storage;
typedef Sacado::PCE::OrthogPoly<double,Storage> pce_type;
template class MueLu::TopRAPFactory<pce_type, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;
template class MueLu::TopSmootherFactory<pce_type, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;

#endif
