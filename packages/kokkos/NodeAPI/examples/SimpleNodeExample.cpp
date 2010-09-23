#include "Kokkos_NodeExampleKernels.hpp"
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_NodeHelpers.hpp>

/** \file SimpleNodeExample.cpp
    \brief A file illustrating some basic usage of the \ref kokkos_node_api "Kokkos Node API"
 */

int main() {
  typedef Kokkos::DefaultNode::DefaultNodeType NODE;
  const int VEC_LENGTH = 100;

  Teuchos::RCP<NODE> node = Kokkos::DefaultNode::getDefaultNode();
  Teuchos::ArrayRCP<int> x = node->allocBuffer<int>( VEC_LENGTH );

  KokkosExamples::initVec( node, x );
  int ret = KokkosExamples::reduceVec( node, x );
  std::cout << "Result is " << ret << std::endl;
  if (ret == (VEC_LENGTH-1)*VEC_LENGTH/2) std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;
}
