/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/Version.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/command_line/CommandLineParser.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/parallel/CouplingVersions.hpp>

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
namespace impl 
{

SplitCommsImpl::SplitCommsImpl()
{}

SplitCommsImpl::SplitCommsImpl(MPI_Comm parentComm, int localColor)
  : m_parentComm(parentComm), m_splitComm(MPI_COMM_NULL), m_localColor(localColor),
    m_finalizationDestructor([this](){free_comms_from_destructor();}),
    m_isInitialized(true)
{
  stk::util::set_coupling_version(parentComm);
  int localRank = stk::parallel_machine_rank(m_parentComm);
  MPI_Comm_split(m_parentComm, m_localColor, localRank, &m_splitComm);

  initialize();
}

SplitCommsImpl::~SplitCommsImpl()
{
  m_finalizationDestructor.destructor();
}


MPI_Comm SplitCommsImpl::get_split_comm() const
{
  STK_ThrowRequireMsg(m_isInitialized, "SplitCommsImpl has not been initialized.");
  return m_splitComm;
}

MPI_Comm SplitCommsImpl::get_parent_comm() const
{
  STK_ThrowRequireMsg(m_isInitialized, "SplitCommsImpl has not been initialized.");
  return m_parentComm;
}


void SplitCommsImpl::setup_pairwise_comms()
{
  int localRank = stk::parallel_machine_rank(m_parentComm);
  int numRanks = stk::parallel_machine_size(m_parentComm);
  std::vector<int> allColors(numRanks);
  MPI_Allgather(&m_localColor, 1, MPI_INT, allColors.data(), 1, MPI_INT, m_parentComm);

  stk::util::sort_and_unique(allColors);

  if (allColors.size() == 2) {
    int otherColor = (m_localColor == allColors[0]) ? allColors[1] : allColors[0];
    //TODO: duplicating the comm is the correct thing to do, but it
    //      breaks some tests
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

void SplitCommsImpl::fill_other_colors()
{
  m_otherColors.clear();
  for (auto& item : m_pairwiseComms) {
    m_otherColors.push_back(item.first);
  }
}

void SplitCommsImpl::compute_all_pairwise_root_ranks()
{
  for (int otherColor : m_otherColors) {
    m_rootRanks[otherColor] = compute_pairwise_root_ranks(m_pairwiseComms[otherColor]);
  }
}

PairwiseRanks SplitCommsImpl::compute_pairwise_root_ranks(MPI_Comm pairwiseComm)
{
  const int localSize = stk::parallel_machine_size(m_splitComm);
  const int globalRank = stk::parallel_machine_rank(pairwiseComm);

  std::vector<int> globalRanks(localSize);
  MPI_Allgather(&globalRank, 1, MPI_INT, globalRanks.data(), 1, MPI_INT, m_splitComm);

  int localRootRank = globalRanks[0];
  int otherRootRank = compute_other_root_rank(globalRanks);

  return PairwiseRanks{localRootRank, otherRootRank};
}

int SplitCommsImpl::compute_other_root_rank(const std::vector<int>& globalRanks)
{
  int otherRootRankCandidate = 0;
  for (int rank : globalRanks) {
    if (otherRootRankCandidate != rank) break;
    otherRootRankCandidate++;
  }
  return otherRootRankCandidate;
}

bool SplitCommsImpl::is_initialized() const
{
  return m_isInitialized;
}

void SplitCommsImpl::free_comms()
{
  STK_ThrowRequireMsg(m_isInitialized, "SplitCommsImpl has not been initialized."); 
  STK_ThrowRequireMsg(!m_haveFreedComms, "SplitCommsImpl has already freed the comms");
  STK_ThrowRequireMsg(!m_freeCommsInDestructor, std::string("SplitCommsImpl is going to free the comms in the destructor. ") +
                                            "Call set_free_comms_in_destructor(false) if you want to manage the memory manually");
  free_comms_impl();
}

void SplitCommsImpl::set_free_comms_in_destructor(bool willFree)
{
  m_freeCommsInDestructor = willFree;
}

bool SplitCommsImpl::get_free_comms_in_destructor() const
{
  return m_freeCommsInDestructor;
}

void SplitCommsImpl::initialize()
{
  setup_pairwise_comms();
  fill_other_colors();
  compute_all_pairwise_root_ranks();
}

MPI_Comm SplitCommsImpl::get_pairwise_comm(int otherColor) const
{
  STK_ThrowRequireMsg(m_isInitialized, "SplitCommsImpl has not been initialized.");
  auto iter = m_pairwiseComms.find(otherColor);
  STK_ThrowRequireMsg(iter != m_pairwiseComms.end(), "SplitCommsImpl with color " << m_localColor <<
                                                 " has no pairwise communicator with color: " << otherColor);
  return iter->second;
}

const std::vector<int>& SplitCommsImpl::get_other_colors() const
{
  STK_ThrowRequireMsg(m_isInitialized, "SplitCommsImpl has not been initialized.");
  return m_otherColors;
}

int SplitCommsImpl::get_local_color() const
{  
  STK_ThrowRequireMsg(m_isInitialized, "SplitCommsImpl has not been initialized.");
  return m_localColor;
}


PairwiseRanks SplitCommsImpl::get_pairwise_root_ranks(int otherColor) const
{
  STK_ThrowRequireMsg(m_isInitialized, "SplitCommsImpl has not been initialized.");
  auto otherColorIter = m_rootRanks.find(otherColor);
  STK_ThrowRequireMsg(otherColorIter != m_rootRanks.end(), "SplitCommsImpl with color " << m_localColor <<
                                                 " has no pairwise communicator with color: " << otherColor);

  return otherColorIter->second;
}


void SplitCommsImpl::free_comms_from_destructor()
{
  if (!get_free_comms_in_destructor())
    return;

  if (!m_isInitialized) {
    std::cerr << "Warning: Cannot free communicators of uninitialized SplitCommsImpl" << std::endl;
    return;
  }

  if (m_haveFreedComms) {
    std::cerr << "Warning: SplitCommsImpl destructor is trying to free communicators that have aready been freed" << std::endl;
    return;
  }

  free_comms_impl();
}

void SplitCommsImpl::free_comms_impl()
{
  int isMPIFinalized;
  MPI_Finalized(&isMPIFinalized);
  if (isMPIFinalized) {
    std::cerr << "Attempting to free commuicators after MPI was finalized."
              << "  The communicators will not be freed, and you have leaked memory" << std::endl;
    return;
  }

  MPI_Comm_free(&m_splitComm);

  for (auto& item : m_pairwiseComms) {
    MPI_Comm pairComm = item.second;
    if (pairComm != m_parentComm) {
      MPI_Comm_free(&pairComm);
    }
  }

  m_haveFreedComms = true;
}

}

bool are_comms_unequal(MPI_Comm global, MPI_Comm local)
{
  int result = 0;
  MPI_Comm_compare(global, local, &result);
  return result == MPI_UNEQUAL;
}


}
}
