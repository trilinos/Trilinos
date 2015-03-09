#include "util.hpp"

#include "sequential_for.hpp"
#include "team_factory.hpp"

using namespace std;

using namespace Example;

// Testing
typedef TeamFactory<TeamPolicy,TeamThreadLoopRegion> TeamFactoryType;

int main (int argc, char *argv[]) {
  SequentialFor(TeamFactoryType::createThreadLoopRegion(TeamPolicy::member_type(), 1, 10),
                [&](int i) {
                  cout << " i = " << i << endl;
                });

  return 0;
}
