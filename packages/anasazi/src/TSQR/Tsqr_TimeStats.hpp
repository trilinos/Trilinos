#ifndef __TSQR_TimeStats_hpp
#define __TSQR_TimeStats_hpp

// #include <cstddef> // size_t
#include <ostream>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // Forward declaration
  template< class Scalar >
  class MessengerBase;

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

    /// Print to out
    /// 
    /// \param out [in/out] Output stream to which to print
    /// \param humanReadable [in] Whether to print in a format easy
    ///   for humans to read, or easy for automatic parsing
    /// \param label [in] If not humanReadable, then print this string
    ///   as a row identifier at the beginning of the row
    /// \param labelLabel [in] If not humanReadable, then use this
    ///   as the column header for the "label" (first) column
    /// \param printHeaders [in] If not humanReadable, then print
    ///   column headers, preceded by a "%" so that the parser will
    ///   ignore the line
    void
    print (std::ostream& out, 
	   const bool humanReadable,
	   const std::string& label,
	   const std::string& labelLabel,
	   const bool printHeaders) const;

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
    int count() const { return count_; }

    /// Produce global time statistics out of all the local ones.
    ///
    static TimeStats
    globalTimeStats (MessengerBase< double >* const comm,
		     const TimeStats& localStats);

  private:
    double min_, max_, mean_, total_;
    int count_;
  };

  

} // namespace TSQR

#endif // __TSQR_TimeStats_hpp
