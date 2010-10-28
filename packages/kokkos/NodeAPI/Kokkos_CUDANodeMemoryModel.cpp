#include "Kokkos_CUDANodeMemoryModel.hpp"
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_FancyOStream.hpp>

namespace Kokkos {

  CUDANodeMemoryModel::CUDANodeMemoryModel() 
  : allocSize_(0)
  {
    clearStatistics();
  }

  void CUDANodeMemoryModel::clearStatistics() {
    numCopiesD2H_ = 0;
    numCopiesH2D_ = 0;
    numCopiesD2D_ = 0;
    bytesCopiedD2H_ = 0;
    bytesCopiedH2D_ = 0;
    bytesCopiedD2D_ = 0;
  }

  void CUDANodeMemoryModel::printStatistics(const RCP< Teuchos::FancyOStream > &os) const {
    using std::setw;
    using std::endl;
    *os << Teuchos::typeName(*this) << " memory transfer statistics" << endl
        << setw(3) << ""      << setw(4) << "" << setw(14) << "Num copies"  << setw(4) << "" << setw(14) << "Bytes copied"  << endl
        << setw(3) << "D2H"   << setw(4) << "" << setw(14) << numCopiesD2H_ << setw(4) << "" << setw(14) << bytesCopiedD2H_ << endl
        << setw(3) << "H2D"   << setw(4) << "" << setw(14) << numCopiesH2D_ << setw(4) << "" << setw(14) << bytesCopiedH2D_ << endl
        << setw(3) << "D2D"   << setw(4) << "" << setw(14) << numCopiesD2D_ << setw(4) << "" << setw(14) << bytesCopiedD2D_ << endl;
  }

}
