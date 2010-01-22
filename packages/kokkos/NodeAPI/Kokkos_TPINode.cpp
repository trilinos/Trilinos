#include "Kokkos_TPINode.hpp"
#include <iostream>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TestForException.hpp>

namespace Kokkos {

  TPINode::TPINode(Teuchos::ParameterList &plist) {
    using std::cout;
    using std::cerr;
    using std::endl;

    curNumThreads_ = plist.get<int>("Num Threads", 0);
    int verbose = plist.get<int>("Verbose",0);
    TEST_FOR_EXCEPTION(curNumThreads_ < 0, std::runtime_error, 
        "TPINode::TPINode(): invalid ""Num Threads"" specification.");
    if (verbose) {
      cout << "TPINode initializing with numThreads == " << curNumThreads_ << std::endl;
    }
    init(curNumThreads_);
  }

  void TPINode::init(int numThreads) {
    if (curNumThreads_ >= 1) {
      TPI_Finalize();
    }
    curNumThreads_ = numThreads;
    if (curNumThreads_ >= 1) {
      TPI_Init(curNumThreads_);
    }
  }

  TPINode::~TPINode()
  {
    if (curNumThreads_ >= 1) {
      TPI_Finalize();
    }
  }

}
