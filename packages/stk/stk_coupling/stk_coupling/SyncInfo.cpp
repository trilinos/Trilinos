/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_coupling/SyncInfo.hpp>
#include <stk_coupling/CommSplitting.hpp>
#include <stk_util/parallel/CommSparse.hpp>

namespace stk
{
namespace coupling
{

void SyncInfo::pack(stk::CommBuffer & b) const
{
  // pack and unpack calls must be in the same order
  b.pack(m_name);
  b.pack(bvals);
  b.pack(ivals);
  b.pack(dvals);
  b.pack(svals);
}

void SyncInfo::unpack(stk::CommBuffer & b)
{
  // pack and unpack calls must be in the same order
  b.unpack(m_name);
  b.unpack(bvals);
  b.unpack(ivals);
  b.unpack(dvals);
  b.unpack(svals);
}

SyncInfo
SyncInfo::exchange(stk::ParallelMachine global, stk::ParallelMachine local)
{

  int globalRank = -1;
  MPI_Comm_rank(global, &globalRank);

  int localRank = -1;
  MPI_Comm_rank(local, &localRank);

  int myRootProc;
  int otherRootProc;
  std::tie(myRootProc, otherRootProc) = stk::coupling::calc_my_root_and_other_root_ranks(global, local);

  if (globalRank == myRootProc)
  {
    ThrowRequireMsg(localRank == 0, "SyncInfo object " << m_name << ": Local rank on the root proc is not zero");
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

    if (localRank == 0)
    {
      recvInfo.pack(comm.send_buffer());
    }

    comm.allocate_buffer();

    if (localRank == 0)
    {
      recvInfo.pack(comm.send_buffer());
    }

    comm.communicate();

    if (localRank != 0)
    {
      recvInfo.unpack(comm.recv_buffer());
    }
  }

  return recvInfo;
}

}
}
