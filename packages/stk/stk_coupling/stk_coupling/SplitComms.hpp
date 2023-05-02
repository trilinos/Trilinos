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
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <memory>
#include "stk_util/parallel/MPIFinalizationCallback.hpp"

namespace stk
{
namespace coupling
{

struct PairwiseRanks
{
  int localColorRoot;
  int otherColorRoot;
};

namespace impl {

class SplitCommsImpl
{
  public:
  SplitCommsImpl();

  SplitCommsImpl(MPI_Comm parentComm, int localColor);

  ~SplitCommsImpl();

  MPI_Comm get_split_comm() const;

  MPI_Comm get_parent_comm() const;

  MPI_Comm get_pairwise_comm(int otherColor) const;

  const std::vector<int>& get_other_colors() const;

  int get_local_color() const;

  PairwiseRanks get_pairwise_root_ranks(int otherColor) const;

  bool is_initialized() const;

  void free_comms();

  void set_free_comms_in_destructor(bool shouldFree);

  bool get_free_comms_in_destructor() const;


private:
  void setup_pairwise_comms();
  void fill_other_colors();
  void compute_all_pairwise_root_ranks();
  void initialize();
  void free_comms_from_destructor();
  void free_comms_impl();

  PairwiseRanks compute_pairwise_root_ranks(MPI_Comm pairwiseComm);
  int compute_other_root_rank(const std::vector<int>& globalRanks);

  MPI_Comm m_parentComm;
  MPI_Comm m_splitComm;
  int m_localColor;

  std::map<int, MPI_Comm> m_pairwiseComms;
  std::vector<int> m_otherColors;
  std::map<int, PairwiseRanks> m_rootRanks;
  MPIFinalizationCallback m_finalizationDestructor;
  bool m_isInitialized = false;
  bool m_freeCommsInDestructor = false;
  bool m_haveFreedComms = false;
};

}

class SplitComms
{
  public:
    SplitComms() :
      m_impl(std::make_shared<impl::SplitCommsImpl>())
    {}

    SplitComms(MPI_Comm parentComm, int localColor) :
      m_impl(std::make_shared<impl::SplitCommsImpl>(parentComm, localColor))
    {}

    MPI_Comm get_split_comm() const { return m_impl->get_split_comm(); }

    MPI_Comm get_parent_comm() const { return m_impl->get_parent_comm(); }

    MPI_Comm get_pairwise_comm(int otherColor) const { return m_impl->get_pairwise_comm(otherColor); }

    const std::vector<int>& get_other_colors() const { return m_impl->get_other_colors(); }

    int get_local_color() const { return m_impl->get_local_color(); }

    PairwiseRanks get_pairwise_root_ranks(int otherColor) const { return m_impl->get_pairwise_root_ranks(otherColor); }

    bool is_initialized() const { return m_impl->is_initialized(); }

    void free_comms() { return m_impl->free_comms(); }

    void set_free_comms_in_destructor(bool shouldFree) { m_impl->set_free_comms_in_destructor(shouldFree); }

    bool get_free_comms_in_destructor() const { return m_impl->get_free_comms_in_destructor(); }

  private:
    std::shared_ptr<impl::SplitCommsImpl> m_impl;
};

bool are_comms_unequal(MPI_Comm global, MPI_Comm local);

}
}

#endif /* STK_COUPLING_SPLITCOMMS_HPP */
