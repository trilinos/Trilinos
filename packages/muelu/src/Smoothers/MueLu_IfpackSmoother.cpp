#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_IFPACK
#include "MueLu_IfpackSmoother.hpp"

namespace MueLu {

  //
  template <>
  RCP<MueLu::SmootherPrototype<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> > GetIfpackSmoother<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps>(std::string const & type, Teuchos::ParameterList const & paramList, int const &overlap) { 
    return rcp( new IfpackSmoother(type, paramList) );
  }

}

#endif // HAVE_MUELU_IFPACK
