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

/*! \file Zoltan2_MetricBase.hpp
 *  \this is a base class, currently for MetricValues and GraphMetricValues
 */

#ifndef ZOLTAN2_METRICBASE_HPP
#define ZOLTAN2_METRICBASE_HPP

#include <Zoltan2_Standards.hpp>

namespace Zoltan2{

/*! \base class for the metric classes. */
template <typename scalar_t>
  class MetricBase {

private:
	/*! \memory allocation is forced by all constructors of this class. */
    void resetValues( int memCount ){
      scalar_t *tmp = new scalar_t [memCount];
      memset(tmp, 0, sizeof(scalar_t) * memCount);
      values_ = arcp(tmp, 0, memCount, true);
    }
    std::string metricName_;

protected:
    ArrayRCP<scalar_t> values_;

public:

    /*! \brief Constructor */
    MetricBase(int memCount, std::string mname) :
      values_(), metricName_(mname) {
    	resetValues(memCount);
        }

    /*! \brief Constructor */
    MetricBase(int memCount) :
      values_(), metricName_("unset") {
    	resetValues(memCount);
        }

    /*! \brief Get the name of the item measured. */
    const std::string &getName() const { return metricName_; }

    /*! \brief Set or reset the name. */
    void setName(std::string name) { metricName_ = name;}

    /*! \hasMetricValue.  */
    bool hasMetricValue(const std::string & metric_name) const {
    	const std::set<std::string> & metrics = getMetrics();
    	return metrics.find(metric_name) != metrics.end();
    }

    /*! \abstract getMetrics.  */
    virtual const std::set<std::string> & getMetrics() const = 0;

    /*! \abstract printHeader.  */
    // MDM - it may be desirable to switch to this form - but not clear which will be more elegant yet
    // right now we don't enforce printHeader for derived classes but it is implemented individually (in utility files)
    //virtual void printHeader(std::ostream &os) = 0;

    /*! \abstract printLine. */
    virtual void printLine(std::ostream &os) const = 0;

    /// \abstract Return a metric value specified by name
    ///
    /// @param metric_name Name of metric to return
    /// @param[out] value metric value returned by reference
    ///
    /// @return Returns a boolean indicated whether or not the metric was returned
    virtual scalar_t getMetricValue(const std::string & metric_name) const = 0;
};

} //namespace Zoltan2
#endif
