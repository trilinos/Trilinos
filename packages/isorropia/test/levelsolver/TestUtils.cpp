#include "TestUtils.h"

int parseNumThreads(const char *arg, std::vector<int> &numThreads) 
{
  int maxT, minT, modT;
  if (sscanf(arg,"%d:%d+%d",&minT,&maxT,&modT) == 3) {
    for (int nt=minT; nt<=maxT; nt += modT) numThreads.push_back(nt);
  }
  else if (sscanf(arg,"%d:%d*%d",&minT,&maxT,&modT) == 3) {
    for (int nt=minT; nt<=maxT; nt *= modT) numThreads.push_back(nt);
  }
  else if (sscanf(arg,"%d",&minT) == 1) {
    if (minT > 0) numThreads.push_back(minT);
    arg = strchr(arg,',');
    while (arg != NULL) {
      if (sscanf(arg+1,"%d", &minT) != 1) {
        break;
      }
      if (minT > 0) numThreads.push_back(minT);
      arg = strchr(arg+1,',');
    }
  }
  else {
    return -1;
  }
  return 0;
}
