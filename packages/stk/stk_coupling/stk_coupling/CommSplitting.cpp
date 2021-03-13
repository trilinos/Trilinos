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
#include <stk_coupling/CommSplitting.hpp>
#include <stk_coupling/impl_VersionUtils.hpp>
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

SplitComms::SplitComms(MPI_Comm parentComm, int localColor)
  : m_parentComm(parentComm), m_splitComm(MPI_COMM_NULL), m_localColor(localColor),
    m_isInitialized(true)
{
  m_compatibilityMode = impl::get_coupling_compatibility_mode(parentComm);
  ThrowRequireMsg(m_compatibilityMode != impl::Incompatible,
    "Incompatible stk_coupling versions detected." << std::endl
    << "Local stk_coupling version " << stk::coupling::version() << ", within STK version: " << stk::version_string());

  int localRank = stk::parallel_machine_rank(m_parentComm);
  MPI_Comm_split(m_parentComm, m_localColor, localRank, &m_splitComm);

  initialize();
}

MPI_Comm SplitComms::get_split_comm() const
{
  ThrowRequireMsg(m_isInitialized, "SplitComms has not been initialized.");
  return m_splitComm;
}

MPI_Comm SplitComms::get_parent_comm() const
{
  ThrowRequireMsg(m_isInitialized, "SplitComms has not been initialized.");
  return m_parentComm;
}


void SplitComms::setup_pairwise_comms()
{
  int localRank = stk::parallel_machine_rank(m_parentComm);
  int numRanks = stk::parallel_machine_size(m_parentComm);
  std::vector<int> allColors(numRanks);
  MPI_Allgather(&m_localColor, 1, MPI_INT, allColors.data(), 1, MPI_INT, m_parentComm);

  stk::util::sort_and_unique(allColors);

  if (allColors.size() == 2) {
    int otherColor = (m_localColor == allColors[0]) ? allColors[1] : allColors[0];
    m_pairwiseComms[otherColor] = m_parentComm;
    return;
  }

  for(unsigned i = 0; i < allColors.size(); i++) {
    for(unsigned j = i + 1; j < allColors.size(); j++) {
      int color1 = allColors[i];
      int color2 = allColors[j];

      MPI_Comm pairwiseComm;
      bool inComm = (m_localColor == color1 || m_localColor == color2);
      int newColor = inComm ? 0 : MPI_UNDEFINED;

      MPI_Comm_split(m_parentComm, newColor, localRank, &pairwiseComm);

      if(m_localColor == color1) {
        m_pairwiseComms[color2] = pairwiseComm;
      }
      else if(m_localColor == color2) {
        m_pairwiseComms[color1] = pairwiseComm;
      }
    }
  }
}

void SplitComms::fill_other_colors()
{
  m_otherColors.clear();
  for (auto& item : m_pairwiseComms) {
    m_otherColors.push_back(item.first);
  }
}

void SplitComms::compute_pairwise_root_ranks()
{
  for (int otherColor : m_otherColors) {
    std::pair<int, int> rootRanks = calc_my_root_and_other_root_ranks(m_pairwiseComms[otherColor], m_splitComm);
    m_rootRanks[otherColor] = PairwiseRanks{rootRanks.first, rootRanks.second};
  }
}

bool SplitComms::is_coupling_version_deprecated() const
{
  return (impl::Deprecated == m_compatibilityMode);
}

void SplitComms::initialize()
{
  setup_pairwise_comms();
  fill_other_colors();
  compute_pairwise_root_ranks();
}

MPI_Comm SplitComms::get_pairwise_comm(int otherColor) const
{
  ThrowRequireMsg(m_isInitialized, "SplitComms has not been initialized.");
  auto iter = m_pairwiseComms.find(otherColor);
  ThrowRequireMsg(iter != m_pairwiseComms.end(), "SplitComms with color " << m_localColor <<
                                                 " has no pairwise communicator with color: " << otherColor);
  return iter->second;
}

const std::vector<int>& SplitComms::get_other_colors() const
{
  ThrowRequireMsg(m_isInitialized, "SplitComms has not been initialized.");
  return m_otherColors;
}

PairwiseRanks SplitComms::get_pairwise_root_ranks(int otherColor) const
{
  ThrowRequireMsg(m_isInitialized, "SplitComms has not been initialized.");
  auto otherColorIter = m_rootRanks.find(otherColor);
  ThrowRequireMsg(otherColorIter != m_rootRanks.end(), "SplitComms with color " << m_localColor <<
                                                 " has no pairwise communicator with color: " << otherColor);

  return otherColorIter->second;
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
