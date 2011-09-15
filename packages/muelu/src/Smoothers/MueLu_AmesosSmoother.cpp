#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_AMESOS

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

#include "MueLu_AmesosSmoother.hpp"

namespace MueLu {

  //
  template <>
  RCP<MueLu::SmootherPrototype<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> > GetAmesosSmoother<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps>(std::string const & type, Teuchos::ParameterList const & paramList) { 
    return rcp( new AmesosSmoother(type, paramList) );
  }

}

#endif // HAVE_MUELU_AMESOS
