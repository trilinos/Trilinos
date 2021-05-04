/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_COUPLING_SPLITCOMMS_HPP
#define STK_COUPLING_SPLITCOMMS_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/impl_VersionUtils.hpp>
#include <string>
#include <vector>
#include <utility>
#include <map>

namespace stk
{
namespace coupling
{

struct PairwiseRanks
{
  int localColorRoot;
  int otherColorRoot;
};

class SplitComms
{
public:
  SplitComms() = default;
  SplitComms(MPI_Comm parentComm, int localColor);

  MPI_Comm get_split_comm() const;

  MPI_Comm get_parent_comm() const;

  MPI_Comm get_pairwise_comm(int otherColor) const;

  const std::vector<int>& get_other_colors() const;

  PairwiseRanks get_pairwise_root_ranks(int otherColor) const;

  bool is_coupling_version_deprecated() const;

  bool is_initialized() const;

private:
  void setup_pairwise_comms();
  void fill_other_colors();
  void compute_all_pairwise_root_ranks();
  void initialize();

  PairwiseRanks compute_pairwise_root_ranks(MPI_Comm pairwiseComm);
  int compute_other_root_rank(const std::vector<int>& globalRanks);

  MPI_Comm m_parentComm;
  MPI_Comm m_splitComm;

  int m_localColor;

  std::map<int, MPI_Comm> m_pairwiseComms;
  std::vector<int> m_otherColors;
  std::map<int, PairwiseRanks> m_rootRanks;
  bool m_isInitialized = false;
  impl::CouplingCompatibilityMode m_compatibilityMode;
};

}
}

#endif /* STK_COUPLING_SPLITCOMMS_HPP */
