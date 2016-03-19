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
#include <Zoltan2_Typedefs.hpp>
#include <AdapterForTests.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>
#include <Zoltan2_OrderingProblem.hpp>

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
using std::cout;
using std::endl;
using std::string;
using std::map;
using std::pair;
using std::ostringstream;
using Teuchos::reduceAll;
using namespace Zoltan2_TestingFramework;

/*! \brief A class used to save problem solutions and timers.
 */
class ComparisonSource
{
  
public:
  
  /*! \brief Destructor.
   */
  ~ComparisonSource()
  {
    if(adapter_kind == "XpetraCrsGraph")
      delete reinterpret_cast<xcrsGraph_adapter *>(adapter.getRawPtr())->getCoordinateInput();
    if(adapter_kind == "XpetraCrsMatrix")
      delete reinterpret_cast<xcrsMatrix_adapter *>(adapter.getRawPtr())->getCoordinateInput();
  }
  /* \brief Add a timer by name to the comparison sources timers map.
   * \param name is the name of the timer to be defined
   */
  void addTimer(const std::string &name)
  {
    timers.insert(std::pair<const std::string &, RCP<Time> >(name,rcp(new Time(name))));
    timers[name]->enable();
  }
  
  RCP<Zoltan2::EvaluatePartition<basic_id_t> > metricObject;
  RCP<base_problem_t> problem;
  RCP<basic_id_t> adapter;
  string problem_kind;
  string adapter_kind;
  std::map<const std::string, RCP<Time> > timers;
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
  void Compare(const ParameterList &pList, const RCP<const Comm<int> > &comm);
  
  /* \brief Add a new source by name to the comparison source map.
   * \param name is the name of the new source
   * \param source a problem source that to be used for comparison to another source
   */
  void AddSource(const string &name, ComparisonSource * source);
  
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
  void CompareSolutions(const string &p1,
                        const string &p2,
                        const RCP<const Comm<int> > &comm);
  
  /* \brief Method called to compare two paritioning solutions
   * \param sourceA is a ptr to problem A's comparison source
   * \param sourceB is a ptr to problem B's comparison source
   * \param comm is the process communicator
   */
  void ComparePartitionSolutions(const ComparisonSource * sourceA,
                                 const ComparisonSource * sourceB,
                                 const RCP<const Comm<int> > &comm);
  
  /* \brief Method called to compare two coloring solutions
   * \param sourceA is a ptr to problem A's comparison source
   * \param sourceB is a ptr to problem B's comparison source
   * \param comm is the process communicator
   */
  void CompareColoringSolutions(const ComparisonSource * sourceA,
                                const ComparisonSource * sourceB,
                                const RCP<const Comm<int> > &comm);
  
  /* \brief Method called to compare two ordering solutions
   * \param sourceA is a ptr to problem A's comparison source
   * \param sourceB is a ptr to problem B's comparison source
   * \param comm is the process communicator
   */
  void CompareOrderingSolutions(const ComparisonSource * sourceA,
                                const ComparisonSource * sourceB,
                                const RCP<const Comm<int> > &comm);
  
  /* \brief Method called to compare the metrics/timers of two problems.
   * \param metricsPlist is a parameter list defining the comparison
   * \param comm is the process communicator
   */
  void CompareMetrics(const ParameterList &metricsPlist,
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
                       const metric_t & metric,
                       const metric_t &ref_metric,
                       const Teuchos::ParameterList & metricPlist,
                       ostringstream &msg);
  
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
                      ostringstream &msg);
  
  /* \brief Method for inserting an array of metrics into a map
   * param[in] metrics an array of metric objects
   *
   * \return a map with metrics assigned to keys assigned by name
   */
  static std::map<const string, const metric_t>
  metricArrayToMap(const ArrayRCP<const metric_t> &metrics);
  
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
  reduceWithMessage(const RCP<const Comm<int> > &comm,const std::string &msg_in,
                    int &local_status, std::ostringstream &msg);
  
};


void ComparisonHelper::AddSource(const string &name, ComparisonSource * source)
{
  typedef std::pair<const string &, RCP<const ComparisonSource> > pair_t;
  this->sources.insert(pair_t(name, RCP<ComparisonSource>(source)));
}

void ComparisonHelper::Compare(const ParameterList &pList, const RCP<const Comm<int> > &comm)
{
  if(pList.isParameter("A") && pList.isParameter("B"))
  {
    // comparing solutions
    
    string pA = pList.get<string>("A");
    if(this->sources.find(pA) == this->sources.end())
    {
      cout << "\nProblem: " + pA + ", was not saved for comparison.";
      cout << "\nThis typically indicates that an error occured while running the problem.";
      cout << "\nSolution comparison FAILED." << endl;
      return;
    }
    
    string pB = pList.get<string>("B");
    if(this->sources.find(pB) == this->sources.end())
    {
      cout << "\nProblem: " + pB + ", was not saved for comparison.";
      cout << "\nThis typically indicates that an error occured while running the problem.";
      cout << "\nSolution comparison FAILED." << endl;
      
      return;
    }
    
    this->CompareSolutions(pA, pB, comm);
  }else if (pList.isParameter("Problem") && pList.isParameter("Reference"))
  {
    // comparing metrics/timers
    string prb = pList.get<string>("Problem");
    if(this->sources.find(prb) == this->sources.end())
    {
      cout << "\nProblem: " + prb + ", was not saved for comparison.";
      cout << "\nThis typically indicates that an error occured while running the problem.";
      cout << "\nMetric comparison FAILED." << endl;
      return;
    }
    
    string ref = pList.get<string>("Reference");
    if(this->sources.find(ref) == this->sources.end())
    {
      cout << "\nReference: " + ref + ", was not saved for comparison.";
      cout << "\nThis typically indicates that an error occured while running the problem.";
      cout << "\nMetric comparison FAILED." << endl;
      return;
    }
    
    this->CompareMetrics(pList,
                         comm);
    
  }else if (pList.isParameter("A") || pList.isParameter("B"))
  {
    if(comm->getRank() == 0)
    {
      cout << "Problem A or Problem B is not specified -- check input.";
      cout <<"\nSolution comparison FAILED." << endl;
    }
    
  }else if (pList.isParameter("Problem") || pList.isParameter("Reference"))
  {
    if(comm->getRank() == 0)
    {
      cout << "Problem or reference is not specified -- check input.";
      cout <<"\nMetric comparison FAILED." << endl;
    }
    
  }
  
}

void ComparisonHelper::CompareSolutions(const string &p1,
                                        const string &p2,
                                        const RCP<const Comm<int> > &comm)
{
  
  if(comm->getRank() == 0) printf("\nComparing: %s and %s\n",p1.c_str(),p2.c_str());
  auto A = this->sources[p1];
  auto B = this->sources[p2];
  if(A->problem_kind != B->problem_kind)
  {
    cout << "Problem A and B are of a different kind and cannot be compared.";
    cout <<"\nSolution comparison FAILED." << endl;
  }else{
    
    if(A->problem_kind == "partitioning")
    {
      this->ComparePartitionSolutions(A.getRawPtr(), B.getRawPtr(), comm);
      
    }else if(A->problem_kind == "coloring")
    {
      this->CompareColoringSolutions(A.getRawPtr(), B.getRawPtr(), comm);
      
    }else if(A->problem_kind == "ordering"){
      
      this->CompareOrderingSolutions(A.getRawPtr(), B.getRawPtr(), comm);
      
    }else{
      cout << "Problem kind not recognized.  Check spelling.";
      cout <<"\nSolution comparison FAILED." << endl;
    }
  }
  
}

void
ComparisonHelper::reduceWithMessage(const RCP<const Comm<int> > &comm, const std::string &msg_in,
                                    int &local_status, std::ostringstream &msg)
{
  comm->barrier();
  int global_buff;
  Teuchos::Ptr<int> global(&global_buff);
  reduceAll<int,int>(*comm.get(), Teuchos::EReductionType::REDUCE_MAX, local_status , global);
  
  local_status = *global;
  if(local_status == 1)
  {
    msg << msg_in;
  }

}

void ComparisonHelper::ComparePartitionSolutions(const ComparisonSource * sourceA,
                                                 const ComparisonSource * sourceB,
                                                 const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  ostringstream status;
  int failed = 0;

  if(!sourceA->problem.getRawPtr()){ failed = 1;}
  ComparisonHelper::reduceWithMessage(comm,
                                      "Solution A is NULL. Solution comparison FAILED.",
                                      failed,
                                      status);
  
  if(!failed && !sourceB->problem.getRawPtr()){ failed = 1;}
  ComparisonHelper::reduceWithMessage(comm,
                                      "Solution B is NULL. Solution comparison FAILED.",
                                      failed,
                                      status);

  if(!failed)
  {
    //    typedef Zoltan2::PartitioningSolution<basic_id_t> partitioning_solution_t; // BDD unused
    // have some solutions lets compare them
    if(partitioning_problem_t * problem_a = reinterpret_cast<partitioning_problem_t *>(sourceA->problem.getRawPtr()))
    {
      if(partitioning_problem_t * problem_b = reinterpret_cast<partitioning_problem_t *>(sourceB->problem.getRawPtr()))
      {
        auto solution_a = problem_a->getSolution();
        auto solution_b = problem_b->getSolution();
        
        if(sourceA->adapter->getLocalNumIDs() != sourceB->adapter->getLocalNumIDs()){failed = 1;}
        ComparisonHelper::reduceWithMessage(comm,
                                            "Number of parts in Solution A != Solution B. \
                                            Partitioning solution comparison FAILED.",
                                            failed,
                                            status);
        
        if(!failed)
        {
          for(size_t i = 0; i < sourceA->adapter->getLocalNumIDs(); i++)
          {
            if(solution_a.getPartListView()[i] != solution_b.getPartListView()[i])
            {
              if(!failed){ failed = 1; }
            }
          }
          
          ComparisonHelper::reduceWithMessage(comm,
                                              "Partitioning solution comparison FAILED.",
                                              failed,
                                              status);
        }
      }else{
        failed = 1;
        ComparisonHelper::reduceWithMessage(comm,
                                            "Solution sets A and B are from different problem types. \
                                            Solution comparison FAILED.",
                                            failed,
                                            status);
      }
      
    }else{
        failed = 1;
        ComparisonHelper::reduceWithMessage(comm,
                                            "Could not cast solution A to valid problem type.  \
                                            Solution comparison FAILED.",
                                            failed,
                                            status);
    }
  }
  
  if(!failed)
  {
    status << "Solution sets A and B are the same. ";
    status << "Solution set comparison PASSED.";
  }
  

  if(rank == 0)
  {
    cout << status.str() << endl;
  }
  
}


void ComparisonHelper::CompareColoringSolutions(const ComparisonSource * sourceA,
                                                const ComparisonSource * sourceB,
                                                const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  ostringstream status;
  int failed = 0;
  
  if(!sourceA->problem.getRawPtr())
  {
    failed = 1;
  }
  ComparisonHelper::reduceWithMessage(comm,
                                      "Solution A is NULL. Solution comparison FAILED.",
                                      failed,
                                      status);
  
  if(!failed && !sourceB->problem.getRawPtr())
  {
    failed = 1;
  }
  ComparisonHelper::reduceWithMessage(comm,
                                      "Solution B is NULL. Solution comparison FAILED.",
                                      failed,
                                      status);
  
  if(!failed)
  {
    // have some solutions lets compare them
    typedef Zoltan2::ColoringProblem<basic_id_t> coloring_problem_t; //BDD unused
    // have some solutions lets compare them
    if(coloring_problem_t * problem_a = reinterpret_cast<coloring_problem_t *>(sourceA->problem.getRawPtr()))
    {
      if(coloring_problem_t * problem_b = reinterpret_cast<coloring_problem_t *>(sourceB->problem.getRawPtr()))
      {
        auto solution_a = problem_a->getSolution();
        auto solution_b = problem_b->getSolution();
        
        if(solution_a->getNumColors() != solution_b->getNumColors())
        {
          failed = 1;
        }
        ComparisonHelper::reduceWithMessage(comm,
                                            "Number of colors for Solution A != Solution B. \
                                            Coloring solution comparison FAILED.",
                                            failed,
                                            status);
        
        if(!failed)
        {
          if(solution_a->getColorsSize() != solution_b->getColorsSize())
          {
            failed = 1;
          }
          ComparisonHelper::reduceWithMessage(comm,
                                              "Size of colors array for Solution A != Solution B. \
                                              Coloring solution comparison FAILED.",
                                              failed,
                                              status);
          
        }
        
        if(!failed)
        {
          for(size_t i = 0; i < solution_a->getColorsSize(); i++)
          {
            if(solution_a->getColors()[i] != solution_b->getColors()[i])
            {
              // fail
              if(!failed) failed = 1;
            }
          }
          ComparisonHelper::reduceWithMessage(comm,
                                              "Coloring solution comparison FAILED.",
                                              failed,
                                              status);
        }
      }else{
        failed = 1;
        ComparisonHelper::reduceWithMessage(comm,
                                            "Solution sets A and B are from different problem types. \
                                            Solution comparison FAILED.",
                                            failed,
                                            status);
      }
      
    }else{
        failed = 1;
        ComparisonHelper::reduceWithMessage(comm,
                                            "Could not cast solution A to valid problem type.  \
                                            Solution comparison FAILED.",
                                            failed,
                                            status);
    }
  }
  
  if(!failed)
  {
    status << "Solution sets A and B are the same. ";
    status << "Solution set comparison PASSED.";
  }
  
  if(rank == 0)
  {
    cout << status.str() << endl;
  }
  
}

void ComparisonHelper::CompareOrderingSolutions(const ComparisonSource * sourceA,
                                                const ComparisonSource * sourceB,
                                                const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  ostringstream status;
  int failed = 0;
  
  if(!sourceA->problem.getRawPtr()){ failed = 1;}
  ComparisonHelper::reduceWithMessage(comm,
                                      "Solution A is NULL. Solution comparison FAILED.",
                                      failed,
                                      status);
  
  if(!failed && !sourceB->problem.getRawPtr()){ failed = 1;}
  ComparisonHelper::reduceWithMessage(comm,
                                      "Solution B is NULL. Solution comparison FAILED.",
                                      failed,
                                      status);
  
  //  if(!failed) //BDD, finish implementation when ordering problem metrics defined
  //  {
  //    // have some solutions lets compare them
  //    typedef Zoltan2::OrderingProblem<basic_id_t> ordering_problem_t;
  //    // have some solutions lets compare them
  //    if(ordering_problem_t * problem_a = reinterpret_cast<ordering_problem_t *>(sourceA->problem.getRawPtr()))
  //    {
  //      if(ordering_problem_t * problem_b = reinterpret_cast<ordering_problem_t *>(sourceB->problem.getRawPtr()))
  //      {
  //
  //      }else{
  //        status << "Solution sets A and B are from different problem types. ";
  //        status << "Solution comparison FAILED.";
  //        failed = true;
  //      }
  //
  //
  //    }else{
  //      if(rank == 0)
  //      {
  //        status << "Could not cast solution A to valid problem type. ";
  //        status << "Solution comparison FAILED.";
  //      }
  //    }
  //  }
  
  
  if(!failed)
  {
    status << "Solution sets A and B are the same. ";
    status << "Solution set comparison PASSED.";
  }
  
  if(rank == 0)
  {
    cout << status.str() << endl;
  }
  
}

// compare metrics
void ComparisonHelper::CompareMetrics(const ParameterList &metricsPlist,
                                      const RCP<const Comm<int> > &comm)
{
  
  int rank = comm->getRank();
  
  //get sources for problema nd reference
  const string prb_name = metricsPlist.get<string>("Problem");
  const string ref_name = metricsPlist.get<string>("Reference");
  if(rank == 0)
  {
    cout  << "\nMetric/Timer comparison of: " << prb_name << " and ";
    cout << ref_name <<" (reference source)\n";
  }
  
  // get sources
  RCP<const ComparisonSource> sourcePrb = this->sources[prb_name];
  RCP<const ComparisonSource> sourceRef = this->sources[ref_name];
  
  // get metric objects
  auto metricObjectPrb = sourcePrb.get()->metricObject.get();
  auto metricObjectRef = sourceRef.get()->metricObject.get();
  
  // get metrics
  std::map<const string, const metric_t> prb_metrics = this->metricArrayToMap
    (metricObjectPrb->getMetrics());
  std::map<const string, const metric_t> ref_metrics = this->metricArrayToMap
    (metricObjectRef->getMetrics());
  
  // get timing data
  std::map< const string, const double> prb_timers = this->timerDataToMap(sourcePrb->timers);
  std::map< const string, const double> ref_timers = this->timerDataToMap(sourceRef->timers);
  
  // get all of the metrics to be tested
  std::queue<ParameterList> metrics = ComparisonHelper::getMetricsToCompare(metricsPlist);
  
  // run comparison
  int all_tests_pass = 1;
  string metric_name;
  while(!metrics.empty())
  {
    // print their names...
    ostringstream msg;
    metric_name = metrics.front().name();
    
    if(prb_metrics.find(metric_name) != prb_metrics.end() &&
       ref_metrics.find(metric_name) != ref_metrics.end())
    {
      if(rank == 0) cout << "\ncomparing metric: " << metric_name << endl;
      if(!ComparisonHelper::metricComparisonTest(comm,
                                                 prb_metrics[metric_name],
                                                 ref_metrics[metric_name],
                                                 metrics.front(), msg))
      {
        all_tests_pass = 0;
      }
      if(rank == 0) cout << msg.str() << endl;
      
    }
    else if(prb_timers.find(metric_name) != prb_timers.end() &&
            ref_timers.find(metric_name) != ref_timers.end())
    {
      if(rank == 0) cout << "\ncomparing timer: " << metric_name << endl;
      if(!ComparisonHelper::timerComparisonTest(comm,
                                                prb_timers.at(metric_name),
                                                ref_timers.at(metric_name),
                                                metrics.front(), msg))
      {
        all_tests_pass = 0;
      }
      
      if(rank == 0) cout << msg.str() << endl;
    }
    
    metrics.pop();
  }
  
  
  if(rank == 0)
  {
    if(all_tests_pass == 1) cout << "\nAll metric/timer comparisons PASSED." << endl;
    else cout << "\nMetric/timer metric comparisons FAILED." << endl;
  }
}

std::map<const string, const metric_t>
ComparisonHelper::metricArrayToMap(const ArrayRCP<const metric_t> &metrics)
{
  typedef std::pair<const string,const metric_t> pair_t;
  std::map<const string, const metric_t> metric_map;
  ArrayRCP<const metric_t>::size_type idx;
  for(idx = 0; idx < metrics.size(); idx++)
  {
    metric_map.insert(pair_t(metrics[idx].getName(),metrics[idx]));
  }
  
  return metric_map;
  
}

std::map<const string, const double>
ComparisonHelper::timerDataToMap(const map<const std::string, RCP<Time> > &timers)
{
  typedef std::pair<const string,const double> pair_t;
  std::map<const string, const double> time_data;
  for(auto &i : timers)
  {
    time_data.insert(pair_t(i.first, i.second->totalElapsedTime()));
  }
  
  return time_data;
}

bool
ComparisonHelper::metricComparisonTest(const RCP<const Comm<int> > &comm,
                                       const Zoltan2::MetricValues<zscalar_t> & metric,
                                       const Zoltan2::MetricValues<zscalar_t> & ref_metric,
                                       const Teuchos::ParameterList & metricPlist,
                                       ostringstream &msg)
{
  // run a comparison of min and max agains a given metric
  // return an error message on failure
  bool pass = true;
  string test_name = metricPlist.name() + " test";
  double ref_value = ref_metric.getMaxImbalance();
  double value = metric.getMaxImbalance();
  
  // want to reduce value to max value for all procs
  
  if (metricPlist.isParameter("lower"))
  {
    double min = metricPlist.get<double>("lower")*ref_value;
    
    if(value < min)
    {
      msg << test_name << " FAILED: imbalance per part, "
      << value << ", less than specified allowable minimum, " << min << ".\n";
      pass = false;
    }else{
      msg << test_name << " PASSED: imbalance per part, "
      << value << ", greater than specified allowable minimum, " << min << ".\n";
    }
  }
  
  if(metricPlist.isParameter("upper" ) && pass != false) {
    
    double max = metricPlist.get<double>("upper") * ref_value;
    if (value > max)
    {
      msg << test_name << " FAILED: imbalance per part, "
      << value << ", greater than specified allowable maximum, " << max << ".\n";
      pass = false;
    }else{
      msg << test_name << " PASSED: imbalance per part, "
      << value << ", less than specified allowable maximum, " << max << ".\n";
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
                                           ostringstream &msg)
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
  
  
  // run a comparison of min and max agains a given metric
  // return an error message on failure
  bool pass = true;
  string test_name = metricPlist.name() + " test";
  if (metricPlist.isParameter("lower"))
  {
    double min = metricPlist.get<double>("lower")*global_ref_time;
    
    if(global_time < min)
    {
      msg << test_name << " FAILED: Minimum time, "
      << time <<
      "[s], less than specified allowable minimum time, " << min <<"[s]"<< ".\n";
      pass = false;
    }else{
      msg << test_name << " PASSED: Minimum time, "
      << time <<
      "[s], greater than specified allowable minimum time, " << min <<"[s]"<< ".\n";
    }
  }
  
  if(metricPlist.isParameter("upper" ) && pass != false) {
    
    double max = metricPlist.get<double>("upper") * global_ref_time;
    if (global_time > max)
    {
      msg << test_name << " FAILED: Maximum time, "
      << global_time <<
      "[s], greater than specified allowable maximum time, " << max <<"[s]"<< ".\n";
      pass = false;
    }else{
      msg << test_name << " PASSED: Maximum time, "
      << global_time <<
      "[s], less than specified allowable maximum time, " << max <<"[s]"<< ".\n";
    }
    
  }
  
  return pass;
}

std::queue<ParameterList>
ComparisonHelper::getMetricsToCompare(const ParameterList &pList)
{
  // extract all of the metrics to be testd
  std::queue<ParameterList> metrics;
  for(auto it = pList.begin(); it != pList.end(); ++it)
  {
    if(pList.isSublist(it->first))
    {
      metrics.push(pList.sublist(it->first));
    }
  }
  
  return metrics;
}

