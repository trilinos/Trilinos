#ifndef KOKKOS_SERIALNODE_HPP_
#define KOKKOS_SERIALNODE_HPP_

#include <Teuchos_ParameterList.hpp>
#include <Kokkos_StandardNodeMemoryModel.hpp>
#include "Kokkos_NodeHelpers.hpp"

namespace Kokkos {

  /** \brief %Kokkos node interface for a serial, CPU node.
      \ingroup kokkos_node_api
   */
  class SerialNode : public StandardNodeMemoryModel {
    public:
      /*! \brief Default constructor, accepts a parameter list but reads no parameters. */
      SerialNode(Teuchos::ParameterList &pl) {}

      //! \begin parallel for skeleton, with a trivial serial implementation. See \ref kokkos_node_api "Kokkos Node API"
      template <class WDP>
      static void parallel_for(int beg, int end, WDP wd) {
        for (int i=beg; i != end; ++i) {
          wd.execute(i);
        }
      }

      //! \begin parallel reduction skeleton, with a trivial serial implementation. See \ref kokkos_node_api "Kokkos Node API"
      template <class WDP>
      static typename WDP::ReductionType
      parallel_reduce(int begin, int end, WDP wd) {
        typename WDP::ReductionType result = wd.identity();
        for (int i=begin; i != end; ++i) {
          result = wd.reduce( result, wd.generate(i) );
        }
        return result;
      }

      //! \begin No-op for SerialNode.
      inline void sync() const {};

  };

  template <> class ArrayOfViewsHelper<SerialNode> : public ArrayOfViewsHelperTrivialImpl<SerialNode> {};

} // end of Kokkos namespace

#endif
