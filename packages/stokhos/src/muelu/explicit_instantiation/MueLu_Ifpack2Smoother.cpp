#include "MueLu_ConfigDefs.hpp"
#include "MueLu_ExplicitInstantiation.hpp"
#include "Stokhos_ConfigDefs.h"

#if defined(HAVE_STOKHOS_MUELU) && defined(HAVE_MUELU_EXPLICIT_INSTANTIATION) && defined(HAVE_STOKHOS_SACADO) && defined(HAVE_MUELU_IFPACK2)

// Sacado headers must be included first so that overloaded operators
// are defined in the muelu template code
#include "Stokhos_Sacado.hpp"
#include "MueLu_Ifpack2Smoother_def.hpp"

typedef Stokhos::StandardStorage<int,double> Storage;
typedef Sacado::PCE::OrthogPoly<double,Storage> pce_type;
template class MueLu::Ifpack2Smoother<pce_type, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;

#endif

