// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_GraphMetrics.hpp
 */

#ifndef ZOLTAN2_GRAPHCMETRICS_HPP
#define ZOLTAN2_GRAPHCMETRICS_HPP

#include <Zoltan2_BaseClassMetrics.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_MetricUtility.hpp>

#define GRAPH_METRICS_TYPE_NAME "GraphMetrics"

namespace Zoltan2{

/*! \GraphMetrics class containing the metrics for one measurable item. */
template <typename scalar_t>
  class GraphMetrics : public BaseClassMetrics<scalar_t> {

public:
/*! \brief Constructor */
GraphMetrics(std::string mname) : BaseClassMetrics<scalar_t>(static_metricNames_.size(), mname) {}

/*! \brief Constructor */
GraphMetrics() : BaseClassMetrics<scalar_t>(static_metricNames_.size()) {}

/*! \brief Get the class type of the metric. */
virtual const std::string & getMetricType() const { return GraphMetrics<scalar_t>::static_metricTypeName_; }

/*! \brief Print a standard header */
static void printHeader(std::ostream &os);

/*! \brief Print a standard line of data that fits under the header. */
virtual void printLine(std::ostream &os) const;

/*! \brief Set the global sum.  */
void setGlobalSum(scalar_t x) { this->setMetricValue("global sum", x);}

/*! \brief Set the global maximum across parts.  */
void setGlobalMax(scalar_t x) { this->setValue("global maximum", x);}

/*! \brief Get the global sum of edge cuts for all parts. */
scalar_t getGlobalSum() const { return this->getMetricValue("global sum");}

/*! \brief Get the global maximum of edge cuts per part across all parts. */
scalar_t getGlobalMax() const { return this->getMetricValue("global maximum");}

/*! \this method is enforced by the base class. */
virtual const std::vector<std::string> & getMetrics() const { return GraphMetrics<scalar_t>::static_metricNames_; }

/*! \setup a static vector of strings. */
static std::vector<std::string> static_metricNames_;

/*! \setup a static string name indicating my class name. */
static std::string static_metricTypeName_;
};  // end class

/*! \static class name for string - used to identify by parameter lists. */
template <typename scalar_t>
std::string GraphMetrics<scalar_t>::static_metricTypeName_ = GRAPH_METRICS_TYPE_NAME;

/*! \lists all metrics for this class. */
template <typename scalar_t>
std::vector<std::string> GraphMetrics<scalar_t>::static_metricNames_ = {
  "global sum",
  "global maximum"
};

template <typename scalar_t>
  void GraphMetrics<scalar_t>::printHeader(std::ostream &os)
{
  os << std::setw(20) << " ";
  os << std::setw(12) << "total" << std::setw(12) << "max";
  os << std::endl;
}

template <typename scalar_t>
  void GraphMetrics<scalar_t>::printLine(std::ostream &os) const
{
  std::string label( this->getName() );

  os << std::setw(20) << label;
  os << std::setw(12) << std::setprecision(4) << this->getMetricValue("global sum");
  os << std::setw(12) << std::setprecision(4) << this->getMetricValue("global maximum");

  os << std::endl;
}

} //namespace Zoltan2
#endif
