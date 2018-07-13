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

/*! \file Zoltan2_BaseClassMetrics.hpp
 *  \this is a base class, currently for MetricValues and GraphMetricValues
 */

#ifndef ZOLTAN2_BASEMETRICVALUES_HPP
#define ZOLTAN2_BASEMETRICVALUES_HPP

#include <Zoltan2_Standards.hpp>

#define UNKNOWN_METRICS_TYPE_NAME "UnknownMetricClass"  // Unknown would be error
#define METRICS_UNSET_STRING "unset"

namespace Zoltan2{

/*! \base class BaseClassMetrics for the metric classes. */
template <typename scalar_t>
  class BaseClassMetrics {

private:
/*! \memory allocation is forced by all constructors of this class. */
void resetValues( int memCount ){
  scalar_t *tmp = new scalar_t [memCount];
  memset(tmp, 0, sizeof(scalar_t) * memCount);
  values_ = arcp(tmp, 0, memCount, true);
}

/*! \name of this metric */
std::string metricName_;

/*! \array of values which is synchronized to constructed memCount in length
 * and the getMetrics string.
 */
ArrayRCP<scalar_t> values_;

protected:

/*! \access to getting values_ */
scalar_t getValue(int enumIndex) const { return values_[enumIndex]; }

/*! \access to setting _values */
void setValue(int enumIndex, scalar_t value) { values_[enumIndex] = value; }

public:

/*! \brief Constructor - for compiling but not used */
BaseClassMetrics() :
  metricName_(METRICS_UNSET_STRING), values_() {
}

/*! \brief Constructor */
BaseClassMetrics(int memCount, std::string mname) :
  metricName_(mname), values_() {
  resetValues(memCount);
}

/*! \brief Constructor */
BaseClassMetrics(int memCount) :
  metricName_(METRICS_UNSET_STRING), values_() {
  resetValues(memCount);
}
    
/*! \virtual Deconstructor */
virtual ~BaseClassMetrics() {}

/*! \abstract printLine. Not abstract so that we can generically support
 * stl containers like maps.
*/
virtual void printLine(std::ostream &os) const {};

/*! \abstract getMetrics. Forces declaration of a static string
 * list of the different metric types
 */
virtual const std::vector<std::string> & getMetrics() const
  { return static_metricNames_; }

/*! \brief Get the class type of the metric. */
virtual const std::string & getMetricType() const
  { return static_unknown_metricTypeName_; }

/*! \brief Get the name of the item measured. */
const std::string &getName() const { return metricName_; }

/*! \brief Set or reset the name. */
void setName(std::string name) { metricName_ = name;}

/*! \hasMetricValue.  */
bool hasMetricValue(const std::string & metric_name) const {
  return( convertMetricNameToIndex( metric_name ) != getMetrics().size() );
}

/*! \ return a metric value specified by name */
scalar_t getMetricValue(const std::string & metric_name) const
{
  size_t metricIndex = convertMetricNameToIndex(metric_name);
  if( metricIndex == getMetrics().size() )
    return 0.0; // throw an error
  return values_[metricIndex];
}

/*! \ set a metric value specified by name */
void setMetricValue(const std::string & metric_name, scalar_t value) const
{
  size_t metricIndex = convertMetricNameToIndex(metric_name);
  if( metricIndex != getMetrics().size() )
    values_[metricIndex] = value;
}

/*! \utility function converts the name to an index. */
size_t convertMetricNameToIndex(const std::string & metric_name) const
{
  const std::vector<std::string> & metricNames = getMetrics();
  size_t metricIndex = std::find(metricNames.begin(), metricNames.end(),
    metric_name) - metricNames.begin();
  return metricIndex; // this can return metricNames.size() if not found
}

/*! \setup a static string name indicating my class name. This stub name exists
 * so that this base class is not virtual to resolve problems with using metrics
 * with stl. It should never be used.
 */
static std::string static_unknown_metricTypeName_;

/*! \setup a static vector of strings. Non virtual so that we can generically
 * support stl containers like maps.
 */
static std::vector<std::string> static_metricNames_;

/*! \This is a list of all possible types - it is used to generate a
 * 'was not used' message, if that's you want.
 */
static std::vector<std::string> static_allMetricNames_;
};

/*! \static class name for string - used to identify by parameter lists. This
 * name should never be used and allows us to be not abstract - so that we can
 * generically support stl containers like maps.
 */
template <typename scalar_t>
std::string BaseClassMetrics<scalar_t>::static_unknown_metricTypeName_ =
  UNKNOWN_METRICS_TYPE_NAME;

/*! \synchronize this with the enum list. Empty list allows us to be not
 * abstract - so that we can generically support stl containers like maps.
 */
template <typename scalar_t>
std::vector<std::string> BaseClassMetrics<scalar_t>::static_metricNames_ = {};

} // namespace Zoltan2
#endif
