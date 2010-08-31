#ifndef __TSQR_TimeStats_hpp
#define __TSQR_TimeStats_hpp

#include <cstddef> // size_t

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class TimeStats
  /// \brief Collect running statistics
  ///
  /// TimeStats collects running statistics on a particular timing,
  /// which is collected count() times.  When you get a new timing,
  /// call update().
  class TimeStats {
  public:
    TimeStats();
    ~TimeStats() {}

    /// Reset the statistics
    ///
    void init ();

    /// Add a new data point and update the running statistics
    ///
    void update (const double curTime);

    /// Minimum value seen thus far (+Inf if no data has been
    /// collected)
    double min() const { return min_; }

    /// Maximum value seen thus far (-Inf if no data has been
    /// collected)
    double max() const { return max_; }

    /// Arithmetic mean thus far (0 if no data has been collected)
    ///
    double mean() const { return mean_; }

    /// Total thus far (0 if no data has been collected)
    ///
    double total() const { return total_; }

    /// Count of data points collected thus far
    ///
    size_t count() const { return count_; }

  private:
    double min_, max_, mean_, total_;
    size_t count_;
  };

} // namespace TSQR

#endif // __TSQR_TimeStats_hpp
