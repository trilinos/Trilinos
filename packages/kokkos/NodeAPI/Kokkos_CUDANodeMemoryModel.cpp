#include "Kokkos_CUDANodeMemoryModel.hpp"

namespace Kokkos {

  CUDANodeMemoryModel::CUDANodeMemoryModel() 
  : numCopiesD2H_(0)
  , numCopiesH2D_(0)
  , numCopiesD2D_(0)
  , bytesCopiedD2H_(0)
  , bytesCopiedH2D_(0)
  , bytesCopiedD2D_(0)
  {
  }

  void CUDANodeMemoryModel::printStatistics() const {
    // finish
  }

}
