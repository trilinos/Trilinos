/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>

#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <iostream>

using namespace stk_classic::diag;

static const int box_count = 100;

namespace stk_classic {
namespace percept {
namespace unit_tests {
/// PLATFORM_NOTE gcc4.3
/// If DIM here is set to 2, the code failes to compile due to dependence of Oct tree code on assuming Dim = 3 in bounding boxes
#define DIM 3

#ifndef REDS
STKUNIT_UNIT_TEST(search, test1)
{
  typedef stk_classic::search::ident::IdentProc<uint64_t,unsigned> IdentProc;
  typedef stk_classic::search::box::PointBoundingBox<IdentProc,float,DIM> BoundingPoint;
  typedef stk_classic::search::box::AxisAlignedBoundingBox<IdentProc,float,DIM> BoundingBox;
  typedef std::vector<std::pair<IdentProc, IdentProc> > IdentProcRelation;

  stk_classic::ParallelMachine  comm = MPI_COMM_WORLD;
  //stk_classic::diag::WriterThrowSafe _write_throw_safe(dw());

  //!dw().m(LOG_SEARCH) << "Use case 1" << stk_classic::diag::push << stk_classic::diag::dendl;

  int parallel_rank = stk_classic::parallel_machine_rank(comm);
  //  int parallel_size = stk_classic::parallel_machine_size(comm);

  std::vector<BoundingBox> domain_vector;

  for (int i = parallel_rank*box_count; i < (parallel_rank + 1)*box_count; ++i) {
#if DIM == 3
    float box[6] = {(float)i,0.0f,0.0f,i+1.0f,1.0f,1.0f };
#else
    float box[4] = {(float)i,0.0f,0.0f,i+1.0f};
#endif

    BoundingBox   domain;
    domain.key.ident = i;
    domain.set_box(box);

    domain_vector.push_back(domain);
  }

  std::vector<BoundingPoint> range_vector;
  for (int i = parallel_rank*box_count; i < (parallel_rank + 1)*box_count; ++i) {
#if DIM == 3
    float center[3] = {i +0.5f,0.5f, 0.5f};
#else
    float center[2] = {i +0.5f,0.5f};
#endif

    BoundingPoint   p;
    p.key.ident = i;
    p.set_center(center);

    range_vector.push_back(p);
  }

  //dw().m(LOG_SEARCH) << "range  " << range_vector << dendl;
  //dw().m(LOG_SEARCH) << "domain " << domain_vector << dendl;

  stk_classic::search::FactoryOrder order;
  order.m_communicator = comm;
  order.m_algorithm = stk_classic::search::FactoryOrder::BIHTREE;

  //dw().m(LOG_SEARCH) << "Search algorithm " << order.m_algorithm << dendl;

  IdentProcRelation relation;

  stk_classic::search::coarse_search(relation, domain_vector, range_vector, order);

  if (0)
  {
    //for (unsigned i = 0; i < relation.size(); i++)
    for (unsigned i = 0; i < 10; i++)
    {
      std::cout << "relation[ " << i << "]= {" << relation[i].first << "} --> { " << relation[i].second << "}" << std::endl;
    }
  }
  //dw().m(LOG_SEARCH) << "relation " << relation << dendl;

  //dw().m(LOG_SEARCH) << stk_classic::diag::pop;
}

#endif
}
}
}
