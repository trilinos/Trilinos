#include <cstdlib>

namespace SEAMS {
  class Stats {
  public:
    Stats();
  
    void newsample(int);
    double mean() const;
    double deviation() const;
    double variance() const;

  private:
    size_t Numnums;
    double Mean;
    double StdDev;
  };
}
