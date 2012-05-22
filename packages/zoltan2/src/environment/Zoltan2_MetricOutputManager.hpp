// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_MetricOutputManager.hpp
    \brief Defines the MetricOutputManager class.
*/
#ifndef _ZOLTAN2_METRICOUTPUT_MANAGER_CPP_
#define _ZOLTAN2_METRICOUTPUT_MANAGER_CPP_

#include <string>
#include <iostream>
#include <fstream>

namespace Zoltan2
{
/*! \brief MetricOutputManager handles output of profiling messages.

   The template parameter is the data type of the value that is being measured.
   An Environment has a MetricOutputManager for timing and one for
   memory profiling.

   Parameters governing may be output:

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

    /*! \brief Constructor for output to a file.
     *   \param rank         the MPI rank of this process.
     *   \param doPrinting  true if this process outputs messages.
     *   \param messageOs      the output stream for messages.
     *   \param metricsOn   true if any process outputs messages.
     */

    MetricOutputManager (int rank, bool doPrinting, std::ofstream &Os, 
      bool metricsOn): me_(rank), 
        os_(static_cast<std::ostream *>(&Os)), osFile_(&Os),
        iPrint_(doPrinting), wePrint_(metricsOn){}


    /*! \brief Constructor
     *   \param rank         the MPI rank of this process.
     *   \param doPrinting  true if this process outputs messages.
     *   \param messageOs      the output stream for messages.
     *   \param metricsOn   true if any process outputs messages.
     */

    MetricOutputManager (int rank, bool doPrinting, std::ostream &Os, 
      bool metricsOn): me_(rank), os_(&Os), osFile_(NULL), 
        iPrint_(doPrinting), wePrint_(metricsOn){}

    /*! \brief Destructor
     */
    ~MetricOutputManager() 
    {
      os_->flush();
      if (osFile_)
        osFile_->close();
    };

    /*! \brief Return the output stream for  messages.
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
    std::ofstream *osFile_;
    bool iPrint_;    // I am printing metrics
    bool wePrint_;   // at least one process is printing metrics
};

} //namespace Zoltan2

#endif
