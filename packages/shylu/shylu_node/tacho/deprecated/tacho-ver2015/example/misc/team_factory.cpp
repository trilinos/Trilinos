#include <Kokkos_Core.hpp>
#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#include <Kokkos_Qthread.hpp>
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp>

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#include "util.hpp"

#include "sequential_for.hpp"
#include "parallel_for.hpp"

#include "team_factory.hpp"

using namespace std;
using namespace Example;

// #define USE_SEQUENTIAL_FOR
#ifdef USE_SEQUENTIAL_FOR
// Debug Interface
// ===============
typedef TeamPolicy PolicyType;
typedef SequentialFor ForType;

template<typename OrdinalType, typename MemberType> 
using LoopType = TeamThreadLoopRegion<OrdinalType,MemberType>;
#else
// Kokkos Interface
// ================
typedef Kokkos::Serial space_type;  // threads space does not work in this example
typedef space_type ExecSpace;

typedef Kokkos::Experimental::TaskPolicy<space_type> PolicyType;
typedef ParallelFor ForType;

template<typename OrdinalType, typename MemberType> 
using LoopType = Kokkos::Impl::TeamThreadRangeBoundariesStruct<OrdinalType,MemberType>;
#endif

typedef TeamFactory<PolicyType,LoopType> TeamFactoryType;

int main (int argc, char *argv[]) {
  const int nthreads = 16;
  ExecSpace::initialize(nthreads);

  ForType(TeamFactoryType::createThreadLoopRegion(PolicyType::member_null(), 1, 10),
          [&](int i) {
            cout << " i = " << i << endl;
          });

  ExecSpace::finalize();

  return 0;
}
