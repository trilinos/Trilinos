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

/*! \file Zoltan2_ComparisonHelper.hpp
 *  \brief Store and compare solution sets from different problems.
 */

#pragma once

#include "Zoltan2_TestHelpers.hpp"
#include "Zoltan2_MetricAnalyzer.hpp"
#include <Zoltan2_Typedefs.hpp>
#include <Zoltan2_EvaluateFactory.hpp>
#include <Zoltan2_ProblemFactory.hpp>
#include <AdapterForTests.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Time.hpp>

#include <sstream>
#include <string>
#include <map>
#include <iostream>

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::Time;
using std::string;
using std::map;
using std::pair;
using Teuchos::reduceAll;
using namespace Zoltan2_TestingFramework;

/*! \brief A class used to save problem solutions and timers.
 */

class ComparisonSource
{
public:
  /* \brief Add a timer by name to the comparison sources timers map.
   * \param name is the name of the timer to be defined
   */
  void addTimer(const std::string &name)
  {
    timers.insert(std::pair<const std::string &, RCP<Time> >(name,rcp(new Time(name))));
    timers[name]->enable();
  }

  void printTimers() 
  {
    for(auto it = timers.begin(); it != timers.end(); ++it) {
      std::cout << it->first << " " << it->second->totalElapsedTime() 
                << std::endl;
    }
  }

  // TODO:  Add method to print a timer summary:  max/min/avg over all procs
  
  std::map<const std::string, RCP<Time> > timers;

  RCP<EvaluateFactory> evaluateFactory;
  RCP<ProblemFactory> problemFactory;
  RCP<AdapterFactory> adapterFactory;
};

/*! \brief A class for comparing solutions, metrics, and timing data of Zoltan2 problems.
 */
class ComparisonHelper
{
public:
  /* \brief Compare the solutions, metrics or timers of two Zoltan2 solutions.
   * \param pList is a parameter list defining the comparison
   * \param comm is the process communicator
   */
  bool Compare(const ParameterList &pList, const RCP<const Comm<int> > &comm);
  
  /* \brief Add a new source by name to the comparison source map.
   * \param name is the name of the new source
   * \param source a problem source that to be used for comparison to another source
   */
  void AddSource(const string &name, RCP<ComparisonSource> source);
  
  /* \brief Return the total number of saved sources.
   */
  size_t getNumberOfSources() const
  {
    return this->sources.size();
  }
  
private:
  map<const string,RCP<const ComparisonSource> > sources;

  /* \brief Method called to compare two solutions
   * \param p1 is the name of problem 1
   * \param p2 is the name of problem 2
   * \param comm is the process communicator
   */
  bool CompareSolutions(const string &p1,
                        const string &p2,
                        const RCP<const Comm<int> > &comm);
  
  /* \brief Safely get parts list by adapter type
   * \param problemFactory is the ProblemFactory
   */
  const zpart_t * getPartListView(RCP<ProblemFactory> problemFactory) const;

  /* \brief Method called to compare two paritioning solutions
   * \param sourceA is a ptr to problem A's comparison source
   * \param sourceB is a ptr to problem B's comparison source
   * \param comm is the process communicator
   */
  bool ComparePartitionSolutions(const ComparisonSource * sourceA,
                                 const ComparisonSource * sourceB,
                                 const RCP<const Comm<int> > &comm);
  
  /* \brief Method called to compare two coloring solutions
   * \param sourceA is a ptr to problem A's comparison source
   * \param sourceB is a ptr to problem B's comparison source
   * \param comm is the process communicator
   */
  bool CompareColoringSolutions(const ComparisonSource * sourceA,
                                const ComparisonSource * sourceB,
                                const RCP<const Comm<int> > &comm);
  
  /* \brief Method called to compare two ordering solutions
   * \param sourceA is a ptr to problem A's comparison source
   * \param sourceB is a ptr to problem B's comparison source
   * \param comm is the process communicator
   */
  bool CompareOrderingSolutions(const ComparisonSource * sourceA,
                                const ComparisonSource * sourceB,
                                const RCP<const Comm<int> > &comm);
  
  /* \brief Safely get metric info by adapter type.
   * \param problemFactory is the ProblemFactory
   * \param metricInfo will be filled
   * \param metricsPlist are the parameters to read
   */
  void loadMetricInfo(RCP<EvaluateFactory> problemFactory,
                 std::vector<MetricAnalyzerInfo> & metricInfo,
                 const ParameterList &metricsPlist);

  /* \brief Method called to compare the metrics/timers of two problems.
   * \param metricsPlist is a parameter list defining the comparison
   * \param comm is the process communicator
   */
  bool CompareMetrics(const ParameterList &metricsPlist,
                      const RCP<const Comm<int> > &comm);
  
  /* \brief Method that compares two metrics and returns a pass/fail message.
   * \param[in] comm is the process communicator
   * \param[in] metric is the metric to be compared to a reference metric
   * \param[in] ref_metric is the reference metric for comparison
   * \param[in] metricPlist is the parameter list defining the metric tolerances
   * \param[out] msg is a returned pass/fail message
   *
   * \return boolean value indicated pass/fail status
   */
  static bool
  metricComparisonTest(const RCP<const Comm<int> > &comm,
                       const MetricAnalyzerInfo & metric,
                       const MetricAnalyzerInfo &ref_metric,
                       const Teuchos::ParameterList & metricPlist,
                       std::ostringstream &msg);
  
  /* \brief Method that compares two timers and returns a pass/fail message.
   * \param[in] comm is the process communicator
   * \param[in] time is the timer data to be compared to a reference metric
   * \param[in] ref_time is the reference timer for comparison
   * \param[in] metricPlist is the parameter list defining the timer tolerances
   * \param[out] msg is a returned pass/fail message
   *
   * \return boolean value indicated pass/fail status
   */
  static bool
  timerComparisonTest(const RCP<const Comm<int> > &comm,
                      const double time,
                      const double ref_time,
                      const Teuchos::ParameterList & metricPlist,
                      std::ostringstream &msg);
  
  /* \brief Method for inserting data from all timers to a map of clocked times
   * param[in] timers a map of timers
   *
   * \return a map with clocked times from timers
   */
  static std::map<const string, const double>
  timerDataToMap(const map<const std::string, RCP<Time> > &timers);
  
  
  /* \brief Method for extracting all methods to compare from a parameter list
   * param[in] plist a parameter list defining 1 or more metric/timer comparisons
   *
   * \return a queue of metric comparison definitions
   */
  static std::queue<ParameterList>
  getMetricsToCompare(const ParameterList & pList);
  
  static void
  reduceWithMessage(const RCP<const Comm<int> > &comm,
                    const std::string &msg_in,
                    int &local_status, std::ostringstream &msg);
  
};


void ComparisonHelper::AddSource(const string &name,
                                 RCP<ComparisonSource> source)
{
  typedef std::pair<const string &, RCP<const ComparisonSource> > pair_t;
  this->sources.insert(pair_t(name, source));
}

bool ComparisonHelper::Compare(const ParameterList &pList,
                               const RCP<const Comm<int> > &comm)
{
  if(pList.isParameter("A") && pList.isParameter("B")) {
    // comparing solutions
    string pA = pList.get<string>("A");
    if(this->sources.find(pA) == this->sources.end())
    {
      std::cout << "\nProblem: " + pA + ", was not saved for comparison.";
      std::cout << "\nThis typically indicates that an error ";
      std::cout << "occurred while running the problem.";
      std::cout << "\nSolution comparison FAILED." << std::endl;
      return false;
    }
    
    string pB = pList.get<string>("B");
    if(this->sources.find(pB) == this->sources.end()) {
      std::cout << "\nProblem: " + pB + ", was not saved for comparison.";
      std::cout << "\nThis typically indicates that an error ";
      std::cout << "occurred while running the problem.";
      std::cout << "\nSolution comparison FAILED." << std::endl;
      return false;
    }
    
    bool bResult = this->CompareSolutions(pA, pB, comm);
    return bResult;
  }
  else if (pList.isParameter("Problem") && pList.isParameter("Reference")) {
    // comparing metrics/timers
    string prb = pList.get<string>("Problem");
    if(this->sources.find(prb) == this->sources.end()) {
      std::cout << "\nProblem: " + prb + ", was not saved for comparison.";
      std::cout << "\nThis typically indicates that an error ";
      std::cout << "occurred while running the problem.";
      std::cout << "\nMetric comparison FAILED." << std::endl;
      return false;
    }
    
    string ref = pList.get<string>("Reference");
    if(this->sources.find(ref) == this->sources.end()) {
      std::cout << "\nReference: " + ref + ", was not saved for comparison.";
      std::cout << "\nThis typically indicates that an error ";
      std::cout << "occurred while running the problem.";
      std::cout << "\nMetric comparison FAILED." << std::endl;
      return false;
    }
    
    bool bResult = this->CompareMetrics(pList, comm);
    return bResult;
  }
  else if (pList.isParameter("A") || pList.isParameter("B"))
  {
    if(comm->getRank() == 0)
    {
      std::cout << "Problem A or Problem B is not specified -- check input.";
      std::cout <<"\nSolution comparison FAILED." << std::endl;
    }
  }
  else if (pList.isParameter("Problem") || pList.isParameter("Reference")) {
    if(comm->getRank() == 0) {
      std::cout << "Problem or reference is not specified -- check input.";
      std::cout <<"\nMetric comparison FAILED." << std::endl;
    }
  }
  else {
    if (comm->getRank() == 0) {
      std::cout << "ComparisonHelper did not understand how to read the xml. ";
      std::cout << "Test FAILED." << std::endl;
    }
  }
  return false;
}

bool ComparisonHelper::CompareSolutions(const string &p1,
                                        const string &p2,
                                        const RCP<const Comm<int> > &comm)
{
  if(comm->getRank() == 0) printf( "\nComparing: %s and %s\n",
    p1.c_str(), p2.c_str());
  auto A = this->sources[p1];
  auto B = this->sources[p2];
  if(A->problemFactory->getProblemName() != B->problemFactory->getProblemName()) {
    std::cout << "Problem A and B are of a different kind and cannot be compared.";
    std::cout <<"\nSolution comparison FAILED." << std::endl;
  }
  else {
    if(A->problemFactory->getProblemName() == "partitioning") {
      return this->ComparePartitionSolutions(A.getRawPtr(), B.getRawPtr(), comm);
    }
    else if(A->problemFactory->getProblemName() == "coloring") {
      return this->CompareColoringSolutions(A.getRawPtr(), B.getRawPtr(), comm);
    }
    else if(A->problemFactory->getProblemName() == "ordering"){
      return this->CompareOrderingSolutions(A.getRawPtr(), B.getRawPtr(), comm);
    }
    else {
      std::cout << "Problem kind: " << A->problemFactory->getProblemName() <<
        " not recognized.  Check spelling.";
      std::cout <<"\nSolution comparison FAILED." << std::endl;
    }
  }
  return false;
}

void ComparisonHelper::reduceWithMessage(const RCP<const Comm<int> > &comm,
                                         const std::string &msg_in,
                                         int &local_status,
                                         std::ostringstream &msg) {
  comm->barrier();
  int global_buff;
  Teuchos::Ptr<int> global(&global_buff);
  reduceAll<int,int>(*comm.get(), Teuchos::EReductionType::REDUCE_MAX,
    local_status , global);
  local_status = *global;
  if (local_status == 1) {
    msg << msg_in;
  }
}

const zpart_t * ComparisonHelper::getPartListView(
  RCP<ProblemFactory> problemFactory) const {
  #define GET_PROBLEM_PARTS(adapterClass)                                  \
      return (rcp_dynamic_cast<PartitioningProblem<adapterClass>>(         \
        problemFactory->getProblem()))->getSolution().getPartListView();
  Z2_TEST_UPCAST(problemFactory->getAdapterType(), GET_PROBLEM_PARTS)
}

bool ComparisonHelper::ComparePartitionSolutions(const ComparisonSource * sourceA,
                                                 const ComparisonSource * sourceB,
                                                 const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  std::ostringstream status;
  int failed = 0;

  if(sourceA->adapterFactory->getMainAdapter()->getLocalNumIDs()
    != sourceB->adapterFactory->getMainAdapter()->getLocalNumIDs()) {
      failed = 1;
  }
  
  ComparisonHelper::reduceWithMessage(comm,
                                  "Number of parts in Solution A != Solution B. \
                                  Partitioning solution comparison FAILED.",
                                  failed, status);
        
  if (!failed) {
    for(size_t i = 0;
      i < sourceA->adapterFactory->getMainAdapter()->getLocalNumIDs(); i++) {
      if(!failed && getPartListView(sourceA->problemFactory)[i] !=
        getPartListView(sourceB->problemFactory)[i]) {
        failed = 1;
        ComparisonHelper::reduceWithMessage(comm,
          "Solution sets A and B have different values for getPartListView(). "
          "Solution comparison FAILED.", failed, status);
      }
    }
  }
  
  if(!failed) {
    status << "Solution sets A and B are the same. ";
    status << "Solution set comparison PASSED.";
  }
  
  if(rank == 0) {
    std::cout << status.str() << std::endl;
  }
  return (failed == 0);
}


bool ComparisonHelper::CompareColoringSolutions(const ComparisonSource * sourceA,
                                                const ComparisonSource * sourceB,
                                                const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  std::ostringstream status;
  int failed = 0;

  // TO DO - implement coloring comparison
  /*
  if(sourceA->problemFactory->getNumColors()
    != sourceB->problemFactory->getNumColors()) {
      failed = 1;
  }
  
  ComparisonHelper::reduceWithMessage(comm,
                              "Number of colors for Solution A != Solution B. \
                              Coloring solution comparison FAILED.",
                              failed, status);
        
  if (!failed) {
    if(sourceA->problemFactory->getColorsSize()
      != sourceB->problemFactory->getColorsSize()) {
        failed = 1;
    }
    ComparisonHelper::reduceWithMessage(comm,
                              "Size of colors array for Solution A != Solution B. \
                              Coloring solution comparison FAILED.",
                              failed,
                              status);
  }

  if (!failed) {
    for(size_t i = 0; i < sourceA->problemFactory->getColorsSize(); i++) {
      if (sourceA->problemFactory->getColors()[i] !=
        sourceB->problemFactory->getColors()[i]) {
          failed = 1;              // fail
      }
    }
    ComparisonHelper::reduceWithMessage(comm,
                                        "Coloring solution comparison FAILED.",
                                        failed,
                                        status);
  }
  */
  
  if (!failed) {
    status << "Solution sets A and B are the same. ";
    status << "Solution set comparison PASSED.";
  }
  
  if (rank == 0) {
    std::cout << status.str() << std::endl;
  }
  return (failed == 0);
}

bool ComparisonHelper::CompareOrderingSolutions(const ComparisonSource * sourceA,
                                                const ComparisonSource * sourceB,
                                                const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  std::ostringstream status;
  int failed = 0;
  
  // TO DO - implement ordering comparison

  if(!failed) {
    status << "Solution sets A and B are the same. ";
    status << "Solution set comparison PASSED.";
  }
  
  if(rank == 0) {
    std::cout << status.str() << std::endl;
  }
  return (failed == 0);
}

// Utility function for safe type conversion of adapter
void ComparisonHelper::loadMetricInfo(RCP<EvaluateFactory> evaluateFactory,
  std::vector<MetricAnalyzerInfo> & metricInfo,
  const ParameterList &metricsPlist) {

  #define LOAD_METRIC_INFO(adapterClass, metricAnalyzerClass)                 \
    RCP<EvaluateBaseClass<adapterClass>> pCast =                              \
      rcp_dynamic_cast<EvaluateBaseClass<adapterClass>>(evaluateFactory->getEvaluateClass());  \
    if(pCast == Teuchos::null) throw std::logic_error(                        \
      "Bad evaluate class cast in loadMetricInfo!"  );                        \
      metricAnalyzerClass analyzer(pCast);                                    \
      analyzer.LoadMetricInfo(metricInfo, metricsPlist.sublist("Metrics"));

  #define LOAD_METRIC_INFO_PARTITIONING(adapterClass)                      \
    LOAD_METRIC_INFO(adapterClass, MetricAnalyzerEvaluatePartition<adapterClass>)

  #define LOAD_METRIC_INFO_ORDERING(adapterClass)                          \
    LOAD_METRIC_INFO(adapterClass, MetricAnalyzerEvaluateOrdering<adapterClass>)

  if(evaluateFactory->getProblemName() == "partitioning") {
    Z2_TEST_UPCAST(evaluateFactory->getAdapterType(), LOAD_METRIC_INFO_PARTITIONING)
  }
  else if(evaluateFactory->getProblemName() == "ordering") {
    Z2_TEST_UPCAST(evaluateFactory->getAdapterType(), LOAD_METRIC_INFO_ORDERING)
  }
  else {
    throw std::logic_error(
      "loadMetricInfo not implemented for this problem type!"  );
  }
}

// compare metrics
bool ComparisonHelper::CompareMetrics(const ParameterList &metricsPlist, const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  
  //get sources for problema nd reference
  const string prb_name = metricsPlist.get<string>("Problem");
  const string ref_name = metricsPlist.get<string>("Reference");
  if(rank == 0) {
    std::cout  << "\nMetric/Timer comparison of: " << prb_name << " and ";
    std::cout << ref_name <<" (reference source)\n";
  }
  
  // get sources
  RCP<const ComparisonSource> sourcePrb = this->sources[prb_name];
  RCP<const ComparisonSource> sourceRef = this->sources[ref_name];

  // get timing data
  std::map< const string, const double> prb_timers = this->timerDataToMap(sourcePrb->timers);
  std::map< const string, const double> ref_timers = this->timerDataToMap(sourceRef->timers);

  // get all of the metrics to be tested
  std::queue<ParameterList> metrics = ComparisonHelper::getMetricsToCompare(metricsPlist);

  // run comparison
  int all_tests_pass = 1;
  string metric_name;
  
  while(!metrics.empty()) {
    // print their names...
    std::ostringstream msg;
    metric_name = metrics.front().name();

    if (metric_name == "Metrics") { // special key word means compare the metrics list
      std::vector<MetricAnalyzerInfo> metricInfoSetPrb;
      std::vector<MetricAnalyzerInfo> metricInfoSetRef;

      loadMetricInfo(sourcePrb.get()->evaluateFactory, metricInfoSetPrb, metricsPlist);
      loadMetricInfo(sourceRef.get()->evaluateFactory, metricInfoSetRef, metricsPlist);

      // there is some redundancy here because the metric info holds both the questions and the results
      // this happened because I wanted to reuse the MetricAnalyzer code for loading metric checks or comparisons
      // we can iterate over either to get the questions
      for (size_t n = 0; n < metricInfoSetPrb.size(); ++n) {
        if(!ComparisonHelper::metricComparisonTest(comm, metricInfoSetPrb[n], metricInfoSetRef[n], metrics.front(), msg)) {
          all_tests_pass = 0;
        }
        if(rank == 0) {
          std::cout << msg.str() << std::endl;
        }
      }
    }
    else if(prb_timers.find(metric_name) != prb_timers.end() && ref_timers.find(metric_name) != ref_timers.end()) {
      if(rank == 0) std::cout << "\ncomparing timer: " << metric_name << std::endl;
      if(!ComparisonHelper::timerComparisonTest(comm,
                                                prb_timers.at(metric_name),
                                                ref_timers.at(metric_name),
                                                metrics.front(), msg)) {
        all_tests_pass = 0;
        if (rank == 0) {
          std::cout << "timer comparison test caused a FAILED event." << std::endl;
        }
      }
      
      if(rank == 0) {
        std::cout << msg.str() << std::endl;
      }
    }

    metrics.pop();
  }
  
  if(rank == 0) {
    if(all_tests_pass == 1) {
      std::cout << "\nAll metric/timer comparisons PASSED." << std::endl;
    }
    else {
      std::cout << "\nMetric/timer metric comparisons FAILED." << std::endl;
    }
  }

  return (all_tests_pass == 1);
}

std::map<const string, const double> ComparisonHelper::timerDataToMap(const map<const std::string, RCP<Time> > &timers)
{
  typedef std::pair<const string,const double> pair_t;
  std::map<const string, const double> time_data;
  for (auto &i : timers) {
    time_data.insert(pair_t(i.first, i.second->totalElapsedTime()));
  }
  return time_data;
}

bool ComparisonHelper::metricComparisonTest(const RCP<const Comm<int> > &comm,
                                       const MetricAnalyzerInfo & metric,
                                       const MetricAnalyzerInfo & ref_metric,
                                       const Teuchos::ParameterList & metricPlist,
                                       std::ostringstream &msg)
{
  // run a comparison of min and max against a given metric
  // return an error message on failure
  bool pass = true;
  string test_name = metricPlist.name() + " test";
  double ref_value = ref_metric.theValue;
  double value = metric.theValue;

  if (ref_value == 0) {
    throw std::logic_error( "The parameter list had a 0 value for the reference value so a percentage cannot be calculated." );
  }
  double percentRatio = value / ref_value;
  
  // want to reduce value to max value for all procs
  
  if (ref_metric.bFoundLowerBound) {
    double min = ref_metric.lowerValue;
    if (percentRatio < min) {
      msg << test_name << " FAILED: " << ref_metric.parameterDescription << ": " << value << " is " << percentRatio << " percent of the reference value " << ref_value << ", which less than specified allowable minimum percent, " << min << ".\n";
      pass = false;
    }
    else {
      msg << test_name << " PASSED: " << ref_metric.parameterDescription << ": " << value << " is " << percentRatio << " percent of the reference value " << ref_value << ", which is greater than specified allowable minimum percent, " << min << ".\n";
    }
  }
  
  if (ref_metric.bFoundUpperBound) {
    double max = ref_metric.upperValue;
    if (percentRatio > max) {
      msg << test_name << " FAILED: " << ref_metric.parameterDescription << ": " << value << " is " << percentRatio << " percent of the reference value " << ref_value << ", which is greater than specified allowable maximum percent, " << max << ".\n";
      pass = false;
    }
    else {
      msg << test_name << " PASSED: " << ref_metric.parameterDescription << ": " << value << " is " << percentRatio << " percent of the reference value " << ref_value << ", which is less than specified allowable maximum percent, " << max << ".\n";
    }
  }

  return pass;
}
// BDD, to do: print metrics even for pass
//             reduce max metric to process 0
//             print only on process 0 --- duh.
bool ComparisonHelper::timerComparisonTest(const RCP<const Comm<int> > &comm,
                                           const double time,
                                           const double ref_time,
                                           const Teuchos::ParameterList & metricPlist,
                                           std::ostringstream &msg)
{
  // Reduce time from test
  double global_time;
  Teuchos::Ptr<double> global(&global_time);
  comm->barrier();
  reduceAll<int, double>(*comm.get(),Teuchos::EReductionType::REDUCE_MAX,time,global);
  
  // Reduce time from reference
  double global_ref_time;
  Teuchos::Ptr<double> globalRef(&global_ref_time);
  comm->barrier();
  reduceAll<int, double>(*comm.get(),Teuchos::EReductionType::REDUCE_MAX,ref_time,globalRef);
  
  // run a comparison of min and max against a given metric
  // return an error message on failure
  bool pass = true;
  string test_name = metricPlist.name() + " test";

  if (metricPlist.isParameter("lower")) {
    double min = metricPlist.get<double>("lower")*global_ref_time;
    
    if (global_time < min) {
      msg << test_name << " FAILED: Minimum time, "
      << time <<
      "[s], less than specified allowable minimum time, " << min <<"[s]"<< ".\n";
      pass = false;
    }
    else {
      msg << test_name << " PASSED: Minimum time, "
      << time <<
      "[s], greater than specified allowable minimum time, " << min <<"[s]"<< ".\n";
    }
  }
  
  if (metricPlist.isParameter("upper" ) && pass != false) {
    double max = metricPlist.get<double>("upper") * global_ref_time;
    if (global_time > max) {
      msg << test_name << " FAILED: Maximum time, "
      << global_time <<
      "[s], greater than specified allowable maximum time, " << max <<"[s]"<< ".\n";
      pass = false;
    }
    else {
      msg << test_name << " PASSED: Maximum time, "
      << global_time <<
      "[s], less than specified allowable maximum time, " << max <<"[s]"<< ".\n";
    }
  }
  
  return pass;
}

std::queue<ParameterList> ComparisonHelper::getMetricsToCompare(const ParameterList &pList)
{
  // extract all of the metrics to be tested
  std::queue<ParameterList> metrics;
  for(auto it = pList.begin(); it != pList.end(); ++it) {
    if (pList.isSublist(it->first)) {
      metrics.push(pList.sublist(it->first));
    }
  }
  return metrics;
}

