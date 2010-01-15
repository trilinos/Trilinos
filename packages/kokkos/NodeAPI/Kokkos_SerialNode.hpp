#ifndef KOKKOS_SERIALNODE_HPP_
#define KOKKOS_SERIALNODE_HPP_

#include <Teuchos_ParameterList.hpp>
#include <Kokkos_StandardNodeMemoryModel.hpp>
#include "Kokkos_NodeHelpers.hpp"

namespace Kokkos {

class SerialNode : public StandardNodeMemoryModel {
  public:
    SerialNode(Teuchos::ParameterList &pl) {}

    template <class WDP>
    static void parallel_for(int beg, int end, WDP wd) {
      for (int i=beg; i != end; ++i) {
        wd.execute(i);
      }
    }

    template <class WDP>
    static typename WDP::ReductionType
    parallel_reduce(int begin, int end, WDP wd) {
      typename WDP::ReductionType result = wd.identity();
      for (int i=begin; i != end; ++i) {
        result = wd.reduce( result, wd.generate(i) );
      }
      return result;
    }

};

template <> class ArrayOfViewsHelper<SerialNode> : public ArrayOfViewsHelperTrivialImpl<SerialNode> {};

}

#endif
