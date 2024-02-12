/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_coupling/SyncInfo.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_util/parallel/CommSparse.hpp>

namespace stk
{
namespace coupling
{

void SyncInfo::pack(stk::CommBuffer & b) const
{
  // pack and unpack calls must be in the same order
  b.pack(m_name);
  m_vals.pack(b);
}

void SyncInfo::unpack(stk::CommBuffer & b)
{
  // pack and unpack calls must be in the same order
  b.unpack(m_name);
  m_vals.unpack(b);
}

SyncInfo SyncInfo::exchange(const SplitComms & splitComms, int otherColor) const
{
  MPI_Comm global = splitComms.get_pairwise_comm(otherColor);
  MPI_Comm local = splitComms.get_split_comm();

  int globalRank = -1;
  MPI_Comm_rank(global, &globalRank);

  int localRank = -1;
  MPI_Comm_rank(local, &localRank);

  PairwiseRanks rootRanks = splitComms.get_pairwise_root_ranks(otherColor);
  int myRootProc = rootRanks.localColorRoot;
  int otherRootProc = rootRanks.otherColorRoot;

  if (globalRank == myRootProc)
  {
    STK_ThrowRequireMsg(localRank == 0, "SyncInfo object " << m_name << ": Local rank on the root proc is not zero");
  }

  SyncInfo recvInfo;

  { // First exchange between rootA <-> rootB
    stk::CommSparse comm(global);

    stk::pack_and_communicate(comm, [&]() {
      if (globalRank == myRootProc)
      {
        pack(comm.send_buffer(otherRootProc));
      }
    });

    if (globalRank == myRootProc)
    {
      recvInfo.unpack(comm.recv_buffer(otherRootProc));
    }
  }

  { // Then broadcast on the local communicator
    stk::CommBroadcast comm(local, 0);
    stk::pack_and_communicate(comm, [&](){
      if (localRank == 0) {
        recvInfo.pack(comm.send_buffer());
      }
    });

    if (localRank != 0)
    {
      recvInfo.unpack(comm.recv_buffer());
    }
  }

  return recvInfo;
}

SyncInfo::ColorToSyncInfoMap SyncInfo::exchange(const SplitComms & splitComms) const
{
  ColorToSyncInfoMap otherInfos;
  std::vector<int> otherColors = splitComms.get_other_colors();
  for(int otherColor : otherColors) {
    otherInfos[otherColor] = exchange(splitComms, otherColor);
  }
  return otherInfos;
}

}
}
