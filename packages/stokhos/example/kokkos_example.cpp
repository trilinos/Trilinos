
#include "Kokkos_DefaultArithmetic.hpp"

int main() {

  typedef Kokkos::ThrustGPUNode Node;
  typedef Kokkos::MultiVector<float,Node> MV;
  typedef Kokkos::DefaultArithmetic<MV> DA;

  Teuchos::ParameterList gpuParams;
  Teuchos::RCP<Node> node = Teuchos::rcp(new Node(gpuParams));

  MV A(node);
  A.initializeValues(5, 5, node->allocBuffer<float>(25), 5);
  DA::Init(A, 1.345);

  return 0;
}
