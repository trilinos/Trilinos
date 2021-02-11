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
#include <stk_coupling/Utils.hpp>
#include <stk_util/parallel/CommSparse.hpp>

namespace stk
{
namespace coupling
{

void SyncInfo::pack(stk::CommBuffer & b) const
{
  // pack and unpack calls must be in the same order
  b.pack(bvals);
  b.pack(ivals);
  b.pack(dvals);
  b.pack(svals);
}

void SyncInfo::unpack(stk::CommBuffer & b)
{
  // pack and unpack calls must be in the same order
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
    ThrowRequireMsg(localRank == 0, "Local rank on the root proc is not zero");
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

void check_sync_mode_consistency(const SyncInfo & myInfo,
                                 const SyncInfo & otherInfo)
{
  ThrowRequireMsg(myInfo.has_value<int>(stk::coupling::TimeSyncMode), "myInfo doesn't contain value for '"<<stk::coupling::TimeSyncMode);   
  stk::coupling::SyncMode myMode = static_cast<stk::coupling::SyncMode>(myInfo.get_value<int>(stk::coupling::TimeSyncMode));
  if (myMode == stk::coupling::Minimum || myMode == stk::coupling::Receive) {
    ThrowRequireMsg(otherInfo.has_value<int>(stk::coupling::TimeSyncMode), "otherInfo doesn't contain value for '"<<stk::coupling::TimeSyncMode);   
  }
  bool otherInfoHasMode = otherInfo.has_value<int>(stk::coupling::TimeSyncMode);
  if (otherInfoHasMode) {
    stk::coupling::SyncMode otherMode = static_cast<stk::coupling::SyncMode>(otherInfo.get_value<int>(stk::coupling::TimeSyncMode));
    bool consistent = (myMode==stk::coupling::Minimum && otherMode==stk::coupling::Minimum) ||
                    ((myMode==stk::coupling::Send) && (otherMode==stk::coupling::Receive)) ||
                    ((myMode==stk::coupling::Receive) && (otherMode==stk::coupling::Send));
    ThrowRequireMsg(consistent, "Inconsistent TimeSyncMode (my mode="<<myMode<<", other mode="<<otherMode
                               <<"). Required to both be Minimum, or one Send and one Receive.");
  }
  else {
    ThrowRequireMsg(myMode == stk::coupling::Send, "Other info has no TimeSyncMode value, which is incompatible with myMode="<<myMode);
  }
}

double choose_value(const SyncInfo & myInfo,
                    const SyncInfo & otherInfo,
                    const std::string & parameterName,
                    stk::coupling::SyncMode syncMode)
{
  ThrowRequireMsg(myInfo.has_value<double>(parameterName), "sync_value: myInfo doesn't contain "<<parameterName);
  ThrowRequireMsg(otherInfo.has_value<double>(parameterName), "sync_value: otherInfo doesn't contain "<<parameterName);

  double myValue = myInfo.get_value<double>(parameterName);
  double otherValue = otherInfo.get_value<double>(parameterName);

  if (syncMode == stk::coupling::Minimum) {
    return std::min(myValue, otherValue);
  }

  return (syncMode == stk::coupling::Receive) ? otherValue : myValue;
}

}
}
