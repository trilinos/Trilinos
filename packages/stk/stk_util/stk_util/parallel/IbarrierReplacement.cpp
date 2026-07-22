#include "IbarrierReplacement.hpp"
#include "stk_util/parallel/CouplingVersions.hpp"

// On some long-running SM tests, the IntelMPI gets confused about 
// how many Ibarriers are in progress at any given time (the correct
// answer is 1, but it thinks there are several hundred thousand).
// For IntelMPI, use a replacement Ibarrier with a specially chosen
// tag to work around this bug.
#ifdef I_MPI_VERSION
  #define USE_IBARRIER_REPLACEMENT true
#else
  #define USE_IBARRIER_REPLACEMENT false
#endif

namespace stk {
namespace impl {

namespace impl {
int sum_op(const int& lhs, const int& rhs)
{
  return lhs + rhs;
}
}

IbarrierReplacement::IbarrierReplacement(MPI_Comm comm, int tag) :
  m_comm(comm)
{
  if (useReplacementIbarrier()) {
    m_allreduce = std::make_shared<IAllreduceReplacement<int>>(comm, MPI_INT, tag);
  }
}


void IbarrierReplacement::startBarrier()
{
  if (useReplacementIbarrier()) {
    m_allreduce->startReduction(&impl::sum_op, m_input_data, m_output_data);
  } else {
    STK_ThrowRequireMsg(m_ibarrier_request == MPI_REQUEST_NULL, 
                    "Cannot start a new IbarrierReplacement before finishing the old one");
    MPI_Ibarrier(m_comm, &m_ibarrier_request);
  }
}

bool IbarrierReplacement::progressBarrier()
{
  int flag;
  if (useReplacementIbarrier()) {
    flag = m_allreduce->progressReduction();
  } else {
    MPI_Test(&m_ibarrier_request, &flag, MPI_STATUS_IGNORE);
  }

  return flag;
}

void IbarrierReplacement::finishBarrier()
{
  if (useReplacementIbarrier()) {
    m_allreduce->finishReduction();
  } else {
    MPI_Wait(&m_ibarrier_request, MPI_STATUS_IGNORE);
  }
}

bool IbarrierReplacement::useReplacementIbarrier() const
{
  return USE_IBARRIER_REPLACEMENT;
}

}
}
