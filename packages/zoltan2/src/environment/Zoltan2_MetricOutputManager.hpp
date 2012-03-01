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
/*! \brief MetricOutputManager handles output of timing, memory use and other values

   The template parameter is the data type of the value that is being measured.
   An Environment has a MetricOutputManager for timing and one for
   memory profiling.

   Parameters governing may be output:

     \li \c timing_level
     \li \c timing_procs
     \li \c timing_output_stream
     \li \c timing_output_file
     \li \c memory_profiling_level
     \li \c memory_profiling_procs
     \li \c memory_profiling_output_stream
     \li \c memory_profiling_output_file

   For more information see at their definitions in
   createAllParameters() in Zoltan2_Parameters.cpp.

*/

template <typename T>
  class MetricOutputManager
{
    public:

    /*! \brief Constructor
     *   \param rank         the MPI rank of this process.
     *   \param doPrinting  true if this process outputs messages.
     *   \param messageOs      the output stream for messages.
     *   \param metricsOn   true if any process outputs messages.
     */

    MetricOutputManager ( int rank, bool doPrinting, std::ostream &Os, bool metricsOn)
    {
      me_ = rank;
      iPrint_ = doPrinting;
      os_ = &Os;
      wePrint_ = metricsOn;
    }

    /*! \brief Destructor
     */
    virtual ~MetricOutputManager() {};

    /*! \brief Return the output stream for debug/status messages.
     */
    inline std::ostream *getOStream() const { return os_; }

    /*! \brief Return true if any process outputs metrics.
     */
    inline bool getMetricsOn() const { return wePrint_; }

    /*! \brief Print a line of information
     *   \param msg  a message to go with the information
     *   \param units  the units in which the value is measured
     *   \param val  the value to print
     */
    inline void print(const std::string &msg, const std::string &units, 
      const T val){
      if (iPrint_)
        *os_ << me_ << ": " << msg << ": " << val << " " << units << std::endl;
    }

    /*! \brief A version of print that takes char to avoid the cost
     *                of converting char * to string.
     *   \param msg  a message to go with the information
     *   \param units  the units in which the value is measured
     *   \param val  the value to print
     */
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
