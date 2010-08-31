#include <Tsqr_TimeStats.hpp>
#include <limits>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  TimeStats::TimeStats() { init(); }

  void
  TimeStats::init () {
    min_ = std::numeric_limits< double >::infinity();
    max_ = -std::numeric_limits< double >::infinity();
    mean_ = double (0);
    total_ = double (0);
    count_ = size_t (0);
  }

  void
  TimeStats::update (const double curTime) {
    total_ += curTime;
    count_++;

    if (curTime < min_)
      min_ = curTime;
    if (curTime > max_)
      max_ = curTime;

    // Mean(1:n) = ((n-1) / n) * Mean(1:n-1) + x(n) / n
    const double scale = double(count_ - 1) / double(count_);
    mean_ = scale * mean_ + curTime / double(count_);
  }

} // namespace TSQR
