
#include <unit_tests/stk_utest_macros.hpp>

#include <mpi.h>

#include <stk_util/diag/Timer.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/diag/IdentProc.hpp>


typedef stk::search::ident::IdentProc<uint64_t,unsigned> IdentProc;
typedef std::vector<std::pair<IdentProc, IdentProc> > IdentProcRelation;
typedef stk::search::box::AxisAlignedBoundingBox<IdentProc,float,3> BoundingVolume;

const int END = 25;

using namespace stk::search;
using namespace stk::diag;
using namespace use_case;

namespace stk_search_coarse_test {

void check_results(const IdentProcRelation & relation, int domain_size)
{
  std::vector<unsigned> count;
  count.resize(domain_size);

  for (unsigned i = 0; i<relation.size(); ++i) {
    count[relation[i].first.ident]++;
  }

  int extra_intersections = 0;
  for (unsigned i=0; i<count.size(); i++) {
    STKUNIT_ASSERT(count[i] >= 64);
    extra_intersections += count[i] - 64;
  }

  //dw() << "\ttotal extra intersections: " << extra_intersections
  //     << " average extra intersections: " << static_cast<float>(extra_intersections) / domain_size << dendl;
}

void testOctree()
{


  unsigned current_id = 0;

  std::vector<BoundingVolume> domain;
  for (int i = 2; i<END-2; ++i) {
    for (int j = 2; j<END-2; ++j) {
      for (int k = 2; k<END-2; ++k) {
        BoundingVolume box;
        box.key.ident = current_id++;
        box.box[0] = i-1.5; box.box[1] = j-1.5; box.box[2] = k-1.5;
        box.box[3] = i+0.5; box.box[4] = j+0.5; box.box[5] = k+0.5;
        domain.push_back(box);
      }
    }
  }

  std::vector<BoundingVolume> range;
  for (int i = 0; i<END; ++i) {
    for (int j = 0; j<END; ++j) {
      for (int k = 0; k<END; ++k) {
        BoundingVolume box;
        box.key.ident = current_id++;
        box.box[0] = i-1; box.box[1] = j-1; box.box[2] = k-1;
        box.box[3] = i+1; box.box[4] = j+1; box.box[5] = k+1;
        range.push_back(box);
      }
    }
  }

  //dw() << dendl;
  //dw() << "Octree 64 expected intersections per test. " << domain.size() << " tests." << dendl;

  FactoryOrder order;
  order.m_algorithm = FactoryOrder::OCTREE;

  IdentProcRelation relation;

  std::ostringstream strout;
  strout << "Octree ";
  stk::diag::Timer timer_tree(strout.str(), use_case::TIMER_SEARCH, use_case::timer());

  {
    stk::diag::TimeBlock timer_tree__(timer_tree);
    coarse_search(relation,range,domain,order);
  }
  check_results(relation,domain.size());

}

void testBihtree()
{
  unsigned current_id = 0;

  std::vector<BoundingVolume> domain;
  for (int i = 2; i<END-2; ++i) {
    for (int j = 2; j<END-2; ++j) {
      for (int k = 2; k<END-2; ++k) {
        BoundingVolume box;
        box.key.ident = current_id++;
        box.box[0] = i-1.5; box.box[1] = j-1.5; box.box[2] = k-1.5;
        box.box[3] = i+0.5; box.box[4] = j+0.5; box.box[5] = k+0.5;
        domain.push_back(box);
      }
    }
  }

  std::vector<BoundingVolume> range;
  for (int i = 0; i<END; ++i) {
    for (int j = 0; j<END; ++j) {
      for (int k = 0; k<END; ++k) {
        BoundingVolume box;
        box.key.ident = current_id++;
        box.box[0] = i-1; box.box[1] = j-1; box.box[2] = k-1;
        box.box[3] = i+1; box.box[4] = j+1; box.box[5] = k+1;
        range.push_back(box);
      }
    }
  }

  //dw() << dendl;
  //dw() << "Bihtree 64 expected intersections per test. " << domain.size() << " tests" << dendl;

  FactoryOrder order;
  order.m_algorithm = FactoryOrder::BIHTREE;

  IdentProcRelation relation;

  stk::diag::Timer timer_tree("Bihtree ", use_case::TIMER_SEARCH, use_case::timer());

  {
    stk::diag::TimeBlock timer_tree__(timer_tree);
    coarse_search(relation,range,domain,order);
  }
  check_results(relation,domain.size());

}

} // namespace stk_search_coarse_test

STKUNIT_UNIT_TEST(UnitTestingOfSearchCoarse, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_search_coarse_test::testOctree ();
  stk_search_coarse_test::testBihtree();
}

