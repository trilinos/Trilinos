#include <Tsqr_TimeStats.hpp>
#include <Tsqr_MessengerBase.hpp>
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
    count_ = int (0);
  }

  void
  TimeStats::update (const double curTime) {
    total_ += curTime;
    count_++;

    if (curTime < min_)
      min_ = curTime;
    if (curTime > max_)
      max_ = curTime;

    // Mean(1:n) = ((n-1) / n) * Mean(1:n-1) + x(n) / n.
    //
    // Casting int to double is exact.
    const double scale = double(count_ - 1) / double(count_);
    mean_ = scale * mean_ + curTime / double(count_);
  }

  void
  TimeStats::print (std::ostream& out, 
		    const bool humanReadable,
		    const std::string& label,
		    const std::string& labelLabel,
		    const bool printHeaders) const
  {
    using std::endl;

    if (humanReadable)
      {
	const char prefix[] = "-- ";
	out << label << ":" << endl;
	if (count() == 0)
	  out << prefix << "No values collected" << endl;
	else if (count() == 1)
	  out << prefix << "One value collected: " << min() << endl;
	else
	  {
	    out << prefix << "Count: " << count() << endl
		<< prefix << "Min:   " << min() << endl
		<< prefix << "Mean:  " << mean() << endl
		<< prefix << "Max:   " << max() << endl
		<< prefix << "Total: " << total() << endl;
	  }
      }
    else
      {
	// "%" identifies this line as a "comment" line to filter out.
	// First print field labels on one line, then print field
	// values on the next line.
	if (printHeaders)
	  {
	    out << "%" << labelLabel
		<< "," << "count"
		<< "," << "min"
		<< "," << "mean"
		<< "," << "max"
		<< "," << "total" 
		<< endl;
	  }
	out << label
	    << "," << count() 
	    << "," << min() 
	    << "," << mean()
	    << "," << max()
	    << "," << total()
	    << endl;
      }
  }

  TimeStats
  TimeStats::globalTimeStats (const Teuchos::RCP< MessengerBase< double > >& comm,
			      const TimeStats& localStats)
  {
    TimeStats globalStats;

    // Casting int to double is exact.
    const double localCount = static_cast<double> (localStats.count());
    const double minCount = comm->globalMin (localCount);
    const double maxCount = comm->globalMax (localCount);
    if (minCount != maxCount)
      throw std::logic_error ("Global stats don\'t make sense, because counts differ");
    globalStats.count_ = localStats.count();
      
    // Casting int to double is exact.
    const double P = static_cast<double> (comm->size());
    globalStats.mean_ = comm->globalSum (localStats.mean() / P);

    globalStats.min_ = comm->globalMin (localStats.min());
    globalStats.max_ = comm->globalMax (localStats.max());
    // Note that this is not the sum of the totals of all the
    // processes, but rather the "global total."  I've chosen to
    // define that as the max of the totals.
    globalStats.total_ = comm->globalMax (localStats.total());

    return globalStats;
  }


} // namespace TSQR
