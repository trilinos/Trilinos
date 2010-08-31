#ifndef __TSQR_TimeStats_hpp
#define __TSQR_TimeStats_hpp

#include <cstddef> // size_t

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class TimeStats
  /// \brief Collect running statistics on a particular timing
  class TimeStats {
  public:
    TimeStats();
    ~TimeStats() {}

    void init ();
    void update (const double curTime);

    double min() const { return min_; }
    double max() const { return max_; }
    double mean() const { return mean_; }
    double total() const { return total_; }
    size_t count() const { return count_; }

  private:
    double min_, max_, mean_, total_;
    size_t count_;
  };

} // namespace TSQR

#endif // __TSQR_TimeStats_hpp
