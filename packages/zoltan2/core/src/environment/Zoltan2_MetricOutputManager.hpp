// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
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
   An Environment has a MetricOutputManager for memory profiling.

   Parameters governing memory output:

     \li \c memory_procs
     \li \c memory_output_stream
     \li \c memory_output_file

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
     *   \param units     name of units for metrics printed
     *   \param width     field width for metric value
     *
     *  It is the caller's responsibility to close and free the
     *  ofstream.
     */

    MetricOutputManager (int rank, bool doPrinting, std::ofstream &Os, 
      bool metricsOn, std::string units=std::string(""), int width=10): 
        me_(rank), 
        os_(static_cast<std::ostream *>(&Os)), osFile_(&Os),
        iPrint_(doPrinting), wePrint_(metricsOn), 
        units_(units), width_(width){}


    /*! \brief Constructor for output to a stream.
     *   \param rank         the MPI rank of this process.
     *   \param doPrinting  true if this process outputs messages.
     *   \param messageOs      the output stream for messages.
     *   \param metricsOn   true if any process outputs messages.
     *   \param units     name of units for metrics printed
     *   \param width     field width for metric value
     */

    MetricOutputManager (int rank, bool doPrinting, std::ostream &Os, 
      bool metricsOn, std::string units=std::string(""), int width=10): 
        me_(rank), os_(&Os), osFile_(NULL), 
        iPrint_(doPrinting), wePrint_(metricsOn),
        units_(units), width_(width){}

    /*! \brief Destructor
     */
    ~MetricOutputManager() {}

    /*! \brief Return the output stream for  messages.
     */
    inline std::ostream *getOStream() const { return os_; }

    /*! \brief Return true if any process outputs metrics.
     */
    inline bool getMetricsOn() const { return wePrint_; }

    /*! \brief Print a line of information
     *   \param msg  a message to go with the information
     *   \param val  the value to print
     */
    inline void print(const std::string &msg, const T val){
      if (iPrint_){
        if (os_){
          os_->width(10); os_->fill('*');
          *os_ << me_ << ": " << val << " " << units_ << ", " << msg << std::endl;
          os_->width(0); os_->fill(' ');
        }
        else{
          osFile_->width(10); osFile_->fill('*');
          *osFile_ << me_ << ": " << val << " " << units_ << ", " << msg << std::endl;
          osFile_->width(0); osFile_->fill(' ');
        }
      }
    }

    /*! \brief A version of print that takes char to avoid the cost
     *                of converting char * to string.
     *   \param msg  a message to go with the information
     *   \param units  the units in which the value is measured
     *   \param val  the value to print
     */
    inline void print(const char *msg, T val){
      if (iPrint_){
        if (os_){
          os_->width(10); os_->fill('*');
          *os_ << me_ << ": " << val << " " << units_ << ", " << msg << std::endl;
          os_->width(0); os_->fill(' ');
        }
        else{
          osFile_->width(10); osFile_->fill('*');
          *osFile_ << me_ << ": " << val << " " << units_ << ", " << msg << std::endl;
          osFile_->width(0); osFile_->fill(' ');
        }
      }
    }

private:

    int me_;
    std::ostream *os_;
    std::ofstream *osFile_;
    bool iPrint_;    // I am printing metrics
    bool wePrint_;   // at least one process is printing metrics
    std::string units_;
    int width_;  
};

} //namespace Zoltan2

#endif
