/*******************************************************************
 *	STATS.C: Statistics routines:
 *
 *	newsample(n):	Add a new sample to the mean/average totals.
 *	mean():	        Returns the mean of the samples.
 *	deviation():	Returns the standard deviation of the sample.
 *
 ********************************************************************/

#include "apr_stats.h"
#include <cmath>

namespace SEAMS {
  Stats::Stats() :
    Numnums(0), Mean(0.0), StdDev(0.0)
  {}

  void Stats::newsample (int n)
  {
    double TMean;

    // See Knuth, TAOCP vol 2, 3rd edition, page 232
    TMean = Mean;
    Numnums++;
    Mean = TMean + (n - TMean) / Numnums;

    if (Numnums > 1)
      StdDev += (n - TMean) * (n - Mean);
  }

  double Stats::mean () const
  {
    return Mean;
  }

  double Stats::variance () const
  {
    return (Numnums > 1) ? StdDev/(Numnums-1) : 0.0;
  }

  double Stats::deviation (void) const
  {
    return std::sqrt(variance());
  }
}

