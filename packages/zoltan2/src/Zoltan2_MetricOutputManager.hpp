// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_MetricOutputManager.hpp

    \brief Metric output manager for Zoltan2 (timing, memory use, etc.)

    Derived from DebugManager
    
*/
#ifndef _ZOLTAN2_METRICOUTPUT_MANAGER_CPP_
#define _ZOLTAN2_METRICOUTPUT_MANAGER_CPP_

#include <string>
#include <iostream>
#include "Teuchos_oblackholestream.hpp"

namespace Zoltan2
{
/*! Zoltan2::MetricOutputManager
    \brief Methods to output timing, memory use and other vals.

*/

template <typename T>
  class MetricOutputManager
{
    public:

    MetricOutputManager ( int rank, bool doPrinting, std::ostream &Os, bool metricsOn)
    {
      me_ = rank;
      iPrint_ = doPrinting;
      os_ = &Os;
      wePrint_ = metricsOn;
    }

    virtual ~MetricOutputManager() {};

    inline void setRank(int p)
    {
      me_ = p;
    }

    inline void setIPrint(bool p)
    {
      iPrint_ = p;
    }

    inline void setOStream(std::ostream &os)
    {
      os_ = &os;
    };

    inline std::ostream *getOStream() const { return os_; }

    inline bool getMetricsOn() const { return wePrint_; }

    inline void print(const std::string &msg, const std::string units, 
      const T val){
      if (iPrint_)
        *os_ << me_ << ": " << msg << ": " << val << " " << units << std::endl;
    }

    inline void print(const char *msg, const char *units, T val){
      if (iPrint_)
        *os_ << me_ << ": " << msg << ": " << val << " " << units << std::endl;
    }

private:

    int me_;
    std::ostream *os_;
    bool iPrint_;    // I am printing metrics
    bool wePrint_;   // at least one process is printing metrics
};

} //namespace Zoltan2

#endif
