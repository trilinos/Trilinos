/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/Version.hpp>
#include <stk_coupling/Version.hpp>
#include <stk_coupling/OldCommSplitting.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/command_line/CommandLineParser.hpp>
#include <stk_util/util/SortAndUnique.hpp>

#include <stdexcept>
#include <string>
#include <vector>
#include <functional>
#include <limits>
#include <cctype>

namespace stk
{
namespace coupling
{

MPI_Comm split_comm(MPI_Comm parentCommunicator, int color)
{
  MPI_Comm splitComm;
  int myrank;
  MPI_Comm_rank(parentCommunicator, &myrank);
  MPI_Comm_split(parentCommunicator, color, myrank, &splitComm);
  return splitComm;
}

bool has_split_comm(MPI_Comm global, MPI_Comm local)
{
  int result = 0;
  MPI_Comm_compare(global, local, &result);
  return result == MPI_UNEQUAL;
}

std::pair<int, int> calc_my_root_and_other_root_ranks(MPI_Comm global, MPI_Comm local)
{
  int num_global_procs = -1;
  MPI_Comm_size(global, &num_global_procs);

  int num_local_procs = -1;
  MPI_Comm_size(local, &num_local_procs);

  if (num_global_procs == num_local_procs)
  {
    return std::make_pair(0, 0);
  }

  int globalRank = -1;
  MPI_Comm_rank(global, &globalRank);

  int localRank = -1;
  MPI_Comm_rank(local, &localRank);

  int my_root_rank_ = -1;
  int other_root_rank_ = -1;

  if (localRank == 0)
  {
    int ind;
    std::vector<int> glbv(num_local_procs);
    glbv[0] = globalRank;
    my_root_rank_ = globalRank;
    for (int i = 1; i < num_local_procs; ++i)
    {
      MPI_Status err;

      MPI_Recv(&ind, 1, MPI_INT, i, 1, local, &err);
      MPI_Recv(glbv.data()+ind, 1, MPI_INT, i, 0, local, &err);
    }

    if (localRank == globalRank)
    {
      int i;
      bool comm_is_non_contiguous = false;
      for (i = 0; i < num_local_procs; ++i)
      {
        if (i != glbv[i])
        {
          comm_is_non_contiguous = true;
          break;
        }
      }
      other_root_rank_ = (comm_is_non_contiguous) ? glbv[i] - 1 : num_local_procs;
    }
    else
    {
      other_root_rank_ = 0;
    }
  }
  else
  {
    MPI_Send(&localRank, 1, MPI_INT, 0, 1, local);
    MPI_Send(&globalRank, 1, MPI_INT, 0, 0, local);
  }

  MPI_Bcast(&my_root_rank_, 1, MPI_INT, 0, local);
  MPI_Bcast(&other_root_rank_, 1, MPI_INT, 0, local);
  return {my_root_rank_, other_root_rank_};
}

}
}
