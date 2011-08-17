#ifndef FAKEKOKKOS_SERIALNODE_HPP_
#define FAKEKOKKOS_SERIALNODE_HPP_

#include <Teuchos_ParameterList.hpp>

namespace Kokkos {

  class SerialNode {
    public:
      SerialNode(Teuchos::ParameterList &pl) {}
  };

} // end of Kokkos namespace

#endif
