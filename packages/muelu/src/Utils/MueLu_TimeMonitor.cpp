#include "MueLu_TimeMonitor.hpp"

namespace MueLu {

  // TODO: this function can be templated (T=double).
  ArrayRCP<double> ReduceMaxMinAvg(double localValue, Teuchos::Comm<int> const &comm, int rootNode) {
    ArrayRCP<double> r = ArrayRCP<double>(3, localValue);

#ifdef HAVE_MPI
    double & maxTime = r[0], & minTime = r[1], & avgTime = r[2];

    // Note: workaround because reduce() is not implemented in Teuchos::Comm
    const Teuchos::MpiComm<int> & mpiComm = dynamic_cast<const Teuchos::MpiComm<int>& >(comm);
    MPI_Comm rawMpiComm = (*mpiComm.getRawMpiComm())();
    //

    // DEBUG std::cout << comm.getRank() << ": " << localValue << std::endl;

    int ntimers=1;
    MPI_Reduce(&localValue, &maxTime, ntimers, MPI_DOUBLE, MPI_MAX, rootNode, rawMpiComm);
    MPI_Reduce(&localValue, &minTime, ntimers, MPI_DOUBLE, MPI_MIN, rootNode, rawMpiComm);
    MPI_Reduce(&localValue, &avgTime, ntimers, MPI_DOUBLE, MPI_SUM, rootNode, rawMpiComm); avgTime /= comm.getSize();
#endif // HAVE_MPI
      
    return r;
  }

} // namespace MueLu
