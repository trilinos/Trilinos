/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_coupling/OldSyncInfo.hpp>
#include <stk_coupling/OldCommSplitting.hpp>
#include <stk_util/parallel/CommSparse.hpp>

#ifndef STK_HIDE_DEPRECATED_CODE  // remove October 2022

namespace stk
{
namespace coupling
{

void OldSyncInfo::pack(stk::CommBuffer & b) const
{
  // pack and unpack calls must be in the same order
  b.pack(svals);
  b.pack(dvals);
  b.pack(ivals);
}

void OldSyncInfo::unpack(stk::CommBuffer & b)
{
  // pack and unpack calls must be in the same order
  b.unpack(svals);
  b.unpack(dvals);
  b.unpack(ivals);
}

OldSyncInfo
OldSyncInfo::exchange(stk::ParallelMachine global, stk::ParallelMachine local)
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
    ThrowRequireMsg(localRank == 0, "OldSyncInfo object Local rank on the root proc is not zero");
  }

  OldSyncInfo recvInfo;

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
#endif