#ifndef stk_util_parallel_IbarrierReplacement
#define stk_util_parallel_IbarrierReplacement

#include "IallreduceReplacement.hpp"
#include <memory>

namespace stk {
namespace impl {

namespace impl {
int sum_op(const int& lhs, const int& rhs);
}

class IbarrierReplacement
{
  public:
    explicit IbarrierReplacement(MPI_Comm comm, int tag);

    void startBarrier();

    bool progressBarrier();

    void finishBarrier();

  private:
    bool useReplacementIbarrier() const;

    MPI_Comm m_comm;
    std::shared_ptr<IAllreduceReplacement<int>> m_allreduce;
    std::vector<int> m_input_data = {0};
    std::vector<int> m_output_data = {0};
    MPI_Request m_ibarrier_request = MPI_REQUEST_NULL;
};


}
}


#endif