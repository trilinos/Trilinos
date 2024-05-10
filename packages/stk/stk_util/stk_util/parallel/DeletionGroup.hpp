#ifndef stk_util_parallel_DeletionGroup
#define stk_util_parallel_DeletionGroup

#include <vector>
#include "stk_util/parallel/MPITag.hpp"
#include "IbarrierReplacement.hpp"

namespace stk {
namespace impl {

class DeletionGroup
{
  public:
    using EntryCountInt = unsigned long long int;
    MPI_Datatype EntryCountIntDatatype = MPI_UNSIGNED_LONG_LONG;

    DeletionGroup(MPI_Comm comm, int barrier_tag);

    ~DeletionGroup();

    MPI_Comm get_comm() const { return m_comm; }

    void insert_tag(int tag) { m_tags.push_back(tag); }

    const std::vector<int>& get_tags() const { return m_tags; }

    // returns true if the MPI_Ibarrier is *semantically* in progress (ie. if start_barrier() has
    // been called but finish_barrier() has not)
    bool is_barrier_in_progress() const { return m_barrierSemanticallyInProgress; }

    void start_barrier(int entryCount);

    void finish_barrier();

    EntryCountInt get_entry_count_barrier_start() const;

    // use this function to trip the MPI progress engine.
    // You still must call finish_barrier() to complete
    // the operation
    void test_barrier();


  private:
   MPI_Comm m_comm;
   bool m_barrierSemanticallyInProgress;
   bool m_barrierActuallyInProgress;
   EntryCountInt m_entryCount;
   std::vector<int> m_tags;
   IbarrierReplacement m_ibarrier;

#ifndef NDEBUG
    MPI_Request m_req_debug;
    std::vector<EntryCountInt> m_entryCounts;
#endif
};

}
}

#endif
