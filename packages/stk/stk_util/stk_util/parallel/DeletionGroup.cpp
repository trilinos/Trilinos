#include <iostream>
#include "stk_util/parallel/DeletionGroup.hpp"

namespace stk {
namespace impl {

DeletionGroup::DeletionGroup(MPI_Comm comm, int barrier_tag)
    : m_comm(comm),
      m_barrierSemanticallyInProgress(false),
      m_barrierActuallyInProgress(false),
      m_entryCount(0),
      m_ibarrier(comm, barrier_tag)
#ifndef NDEBUG
      ,
      m_req_debug(MPI_REQUEST_NULL)
#endif
{
#ifndef NDEBUG
  int commSize, myRank;
  MPI_Comm_rank(comm, &myRank);
  MPI_Comm_size(comm, &commSize);
  if (myRank == 0)
    m_entryCounts.resize(commSize);
#endif
}

DeletionGroup::~DeletionGroup()
{
  if (m_barrierActuallyInProgress)
    std::cerr << "Error: DeletionGroup destructed without calling finish_barrier() first. "
              << "MPI may terminate abnormally or produce undefined behavior" << std::endl;
}


void DeletionGroup::start_barrier(int entryCount)
{ 
  assert(!m_barrierSemanticallyInProgress);

  m_ibarrier.startBarrier();

  m_barrierActuallyInProgress     = true;
  m_barrierSemanticallyInProgress = true;
  m_entryCount = entryCount;


#ifndef NDEBUG
  MPI_Igather(&m_entryCount, 1, EntryCountIntDatatype, m_entryCounts.data(), 1,
              EntryCountIntDatatype, 0, m_comm, &m_req_debug); 
#endif
}

void DeletionGroup::finish_barrier()
{
  assert(m_barrierSemanticallyInProgress);
  if (m_barrierActuallyInProgress) {
    m_ibarrier.finishBarrier();
    m_barrierActuallyInProgress     = false;
    m_barrierSemanticallyInProgress = false;
  }

#ifndef NDEBUG
  int myRank;
  MPI_Comm_rank(m_comm, &myRank);
  MPI_Wait(&m_req_debug, MPI_STATUS_IGNORE);

  if (myRank == 0) {
    for (unsigned int i=1; i < m_entryCounts.size(); ++i) {
      assert(m_entryCounts[i] == m_entryCounts[0]);
    }
  }
#endif
}

DeletionGroup::EntryCountInt DeletionGroup::get_entry_count_barrier_start() const
{ 
  assert(m_barrierSemanticallyInProgress);
  return m_entryCount;
}


void DeletionGroup::test_barrier()
{
  assert(m_barrierSemanticallyInProgress);
  if (!m_barrierActuallyInProgress)
    return;

  m_barrierActuallyInProgress = !(m_ibarrier.progressBarrier());
}

}
}
