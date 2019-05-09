#include <Kokkos_Core.hpp> 
#include "util.hpp"

#include "sequential_for.hpp"

using namespace std;
using namespace Example;

typedef int ordinal_type;

int main (int argc, char *argv[]) {
  Kokkos::initialize();
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  { 
    const ordinal_type begin = 1, end = 11;
    const TeamThreadMember member;  
    SequentialFor(TeamThreadLoop(member, begin, end),
                  [&](ordinal_type i) {
                    cout << " i = " << i << endl;
                  });
  }
  
  Kokkos::finalize();
  
  return 0;
}
