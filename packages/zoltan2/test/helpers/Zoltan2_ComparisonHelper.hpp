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



class ComparisonSource
{
public:
  
  typedef AdapterForTests::base_adapter_t base_t;
  typedef AdapterForTests::basic_id_t basic_id_t; // basic_identifier_type
  typedef AdapterForTests::xpetra_mv_adapter xpetra_mv_t; // xpetra_mv_type
  typedef AdapterForTests::xcrsGraph_adapter xcrsGraph_t;
  typedef AdapterForTests::xcrsMatrix_adapter xcrsMatrix_t;
  typedef AdapterForTests::basic_vector_adapter basic_vector_t;
  
  typedef Zoltan2::Problem<base_t> problem_t;
  typedef Zoltan2::PartitioningProblem<base_t> partioning_problem_t; // base abstract type
  typedef Zoltan2::PartitioningProblem<basic_id_t> basic_problem_t; // basic id problem type
  typedef Zoltan2::PartitioningProblem<xpetra_mv_t> xpetra_mv_problem_t; // xpetra_mb problem type
  typedef Zoltan2::PartitioningProblem<xcrsGraph_t> xcrsGraph_problem_t; // xpetra_mb problem type
  typedef Zoltan2::PartitioningProblem<xcrsMatrix_t> xcrsMatrix_problem_t; // xpetra_mb problem type
  typedef Zoltan2::PartitioningProblem<basic_vector_t> basicVector_problem_t; // xpetra_mb problem type
  

  ~ComparisonSource()
  {
    if(adapter_kind == "XpetraCrsGraph")
    delete reinterpret_cast<xcrsGraph_t *>(adapter.getRawPtr())->getCoordinateInput();
    if(adapter_kind == "XpetraCrsMatrix")
    delete reinterpret_cast<xcrsMatrix_t *>(adapter.getRawPtr())->getCoordinateInput();
  }
  
  void addTimer(const std::string &name)
  {
    timers.insert(std::pair<const std::string &, RCP<Time> >(name,rcp(new Time(name))));
    timers[name]->enable();
  }
  
  RCP<basic_problem_t> problem;
  RCP<basic_id_t> adapter;
  string problem_kind;
  string adapter_kind;
  std::map<const std::string, RCP<Time> > timers;
};


class ComparisonHelper
{
public:
  
  typedef AdapterForTests::base_adapter_t base_t;
  typedef AdapterForTests::basic_id_t basic_id_t; // basic_identifier_type
  typedef AdapterForTests::xpetra_mv_adapter xpetra_mv_t; // xpetra_mv_type
  typedef AdapterForTests::xcrsGraph_adapter xcrsGraph_t;
  typedef AdapterForTests::xcrsMatrix_adapter xcrsMatrix_t;
  typedef AdapterForTests::basic_vector_adapter basic_vector_t;
  
  typedef Zoltan2::Problem<base_t> problem_t;
  typedef Zoltan2::PartitioningProblem<base_t> partioning_problem_t; // base abstract type
  typedef Zoltan2::PartitioningProblem<basic_id_t> basic_problem_t; // basic id problem type
  typedef Zoltan2::PartitioningProblem<xpetra_mv_t> xpetra_mv_problem_t; // xpetra_mb problem type
  typedef Zoltan2::PartitioningProblem<xcrsGraph_t> xcrsGraph_problem_t; // xpetra_mb problem type
  typedef Zoltan2::PartitioningProblem<xcrsMatrix_t> xcrsMatrix_problem_t; // xpetra_mb problem type
  typedef Zoltan2::PartitioningProblem<basic_vector_t> basicVector_problem_t; // xpetra_mb problem type
  
  typedef const Zoltan2::MetricValues<zscalar_t> metric_t;
  
  void Compare(const ParameterList &pList, const RCP<const Comm<int> > &comm);
  
  void AddSource(const string &name, ComparisonSource * source);
  
  size_t getNumberOfSources() const
  {
    return this->sources.size();
  }
  
private:
  map<const string,RCP<const ComparisonSource> > sources;
  
  
  // Solution comparisons
  void CompareSolutions(const string &p1,
                        const string &p2,
                        const RCP<const Comm<int> > &comm);
  
  
  void ComparePartitionSolutions(const ComparisonSource * sourceA,
                                 const ComparisonSource * sourceB,
                                 const RCP<const Comm<int> > &comm);
  
  void CompareColoringSolutions(const ComparisonSource * sourceA,
                                const ComparisonSource * sourceB,
                                const RCP<const Comm<int> > &comm);
  
  void CompareOrderingSolutions(const ComparisonSource * sourceA,
                                const ComparisonSource * sourceB,
                                const RCP<const Comm<int> > &comm);
  
  // metric comparisons
  void CompareMetrics(const ParameterList &metricsPlist,
                      const RCP<const Comm<int> > &comm);
  
  static bool
  metricComparisonTest(const metric_t & metric,
                       const metric_t &ref_metric,
                       const Teuchos::ParameterList & metricPlist,
                       ostringstream &msg);
  
  static bool
  timerComparisonTest(const double time,
                      const double ref_time,
                      const Teuchos::ParameterList & metricPlist,
                      ostringstream &msg);
  
  static std::map<const string, const metric_t>
  metricArrayToMap(const ArrayRCP<const metric_t> &metrics);
  
  static std::map<const string, const double>
  timerDataToMap(const map<const std::string, RCP<Time> > &timers);
  
  static std::queue<ParameterList>
  getMetricsToCompare(const ParameterList & pList);
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
      cout << "\nThis typically indicates that an error occured while running the probelm.";
      cout << "\nSolution comparison FAILED." << endl;
      return;
    }
    
    string pB = pList.get<string>("B");
    if(this->sources.find(pB) == this->sources.end())
    {
      cout << "\nProblem: " + pB + ", was not saved for comparison.";
      cout << "\nThis typically indicates that an error occured while running the probelm.";
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
      cout << "\nThis typically indicates that an error occured while running the probelm.";
      cout << "\nMetric comparison FAILED." << endl;
      return;
    }
    
    string ref = pList.get<string>("Reference");
    if(this->sources.find(ref) == this->sources.end())
    {
      cout << "\nReference: " + ref + ", was not saved for comparison.";
      cout << "\nThis typically indicates that an error occured while running the probelm.";
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

void ComparisonHelper::ComparePartitionSolutions(const ComparisonSource * sourceA,
                                                 const ComparisonSource * sourceB,
                                                 const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  ostringstream status;
  bool failed = false;
  
  if(!sourceA->problem.getRawPtr())
  {
    status << "Solution A is NULL. Solution comparison FAILED.";
    failed = true;
  }
  if(!failed && !sourceB->problem.getRawPtr())
  {
    status << "Solution B is NULL. Solution comparison FAILED.";
    failed = true;
  }
  
  if(!failed)
  {
    typedef Zoltan2::PartitioningSolution<basic_id_t> partitioning_solution_t;
    // have some solutions lets compare them
    if(basic_problem_t * problem_a = reinterpret_cast<basic_problem_t *>(sourceA->problem.getRawPtr()))
    {
      if(basic_problem_t * problem_b = reinterpret_cast<basic_problem_t *>(sourceB->problem.getRawPtr()))
      {
        auto solution_a = problem_a->getSolution();
        auto solution_b = problem_b->getSolution();
        
        if(sourceA->adapter->getLocalNumIDs() != sourceB->adapter->getLocalNumIDs())
        {
          status << "Number of partitions in Solution A != Solution B. ";
          status <<"Partitioning solution comparison FAILED.";
          failed = true;
        }
        
        if(rank == 0)
        {
          fprintf(stdout, "Parts A: %lu, Parts B: %lu\n", sourceA->adapter->getLocalNumIDs(),
                  sourceB->adapter->getLocalNumIDs());
        }
        
        if(!failed)
        for(size_t i = 0; i < sourceA->adapter->getLocalNumIDs(); i++)
        {
          if(solution_a.getPartListView()[i] != solution_b.getPartListView()[i])
          {
            // fail
            //              status << "Partition_A[" <<  i << "] = " << solution_a.getPartListView()[i];
            //              status << ", and Partition_B[" <<  i << "] = " << solution_b.getPartListView()[i];
            
            if(!failed){
              status <<"Partitioning solution comparison FAILED.";
              failed = true;
            }
          }
        }
        
      }else{
        status << "Solution sets A and B are from different problem types. ";
        status << "Solution comparison FAILED.";
        failed = true;
      }
      
    }else{
      if(rank == 0)
      {
        status << "Could not cast solution A to valid problem type. ";
        status << "Solution comparison FAILED.";
      }
    }
  }
  
  if(!failed)
  {
    status << "Solution sets A and B are the same. ";
    status << "Solution set comparison PASSED.";
  }
  
  if(rank == 0)
  {
    cout << "Rank " << rank <<": ";
    cout << status.str() << endl;
  }
  
}

void ComparisonHelper::CompareColoringSolutions(const ComparisonSource * sourceA,
                                                const ComparisonSource * sourceB,
                                                const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  ostringstream status;
  bool failed = false;
  
  if(!sourceA->problem.getRawPtr())
  {
    status << "Solution A is NULL. Solution comparison FAILED.";
    failed = true;
  }
  if(!failed && !sourceB->problem.getRawPtr())
  {
    status << "Solution B is NULL. Solution comparison FAILED.";
    failed = true;
  }
  
  if(!failed)
  {
    // have some solutions lets compare them
    typedef Zoltan2::ColoringProblem<basic_id_t> coloring_problem_t;
    // have some solutions lets compare them
    if(coloring_problem_t * problem_a = reinterpret_cast<coloring_problem_t *>(sourceA->problem.getRawPtr()))
    {
      if(coloring_problem_t * problem_b = reinterpret_cast<coloring_problem_t *>(sourceB->problem.getRawPtr()))
      {
        auto solution_a = problem_a->getSolution();
        auto solution_b = problem_b->getSolution();
        
        if(solution_a->getNumColors() != solution_b->getNumColors())
        {
          status << "Number of colors for Solution A != Solution B. ";
          status <<"Coloring solution comparison FAILED.";
          failed = true;
        }
        
        if(!failed)
        if(solution_a->getColorsSize() != solution_b->getColorsSize())
        {
          status << "Size of colors array for Solution A != Solution B. ";
          status <<"Coloring solution comparison FAILED.";
          failed = true;
        }
        
        if(!failed)
        for(size_t i = 0; i < solution_a->getColorsSize(); i++)
        {
          if(solution_a->getColors()[i] != solution_b->getColors()[i])
          {
            // fail
            status << "Colors_A[" <<  i << "] = " << solution_a->getColors()[i];
            status << ", and Colors_A[" <<  i << "] = " << solution_b->getColors()[i];
            status <<"\n Coloring solution comparison FAILED." <<"\n";
            if(!failed) failed = true;
          }
        }
      }else{
        status << "Solution sets A and B are from different problem types. ";
        status << "Solution comparison FAILED.";
        failed = true;
      }
      
    }else{
      if(rank == 0)
      {
        status << "Could not cast solution A to valid problem type. ";
        status << "Solution comparison FAILED.";
      }
    }
  }
  
  if(!failed)
  {
    status << "Solution sets A and B are the same. ";
    status << "Solution set comparison PASSED.";
  }
  
  cout << "Rank " << rank <<": ";
  cout << status.str() << endl;
  
}

void ComparisonHelper::CompareOrderingSolutions(const ComparisonSource * sourceA,
                                                const ComparisonSource * sourceB,
                                                const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  ostringstream status;
  bool failed = false;
  
  if(!sourceA->problem.getRawPtr())
  {
    status << "Solution A is NULL. Solution comparison FAILED.";
    failed = true;
  }
  if(!failed && !sourceB->problem.getRawPtr())
  {
    status << "Solution B is NULL. Solution comparison FAILED.";
    failed = true;
  }
  
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
  
  cout << "Rank " << rank <<": ";
  cout << status.str() << endl;
  
}

// compare metrics
void ComparisonHelper::CompareMetrics(const ParameterList &metricsPlist,
                                      const RCP<const Comm<int> > &comm)
{
  
  typedef std::pair<const string,double> pair_t;
  
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
  
  // get problems
  auto problem = sourcePrb.get()->problem.get();
  auto reference = sourceRef.get()->problem.get();
  
  // get metrics
  std::map<const string, const metric_t> prb_metrics = this->metricArrayToMap(problem->getMetrics());
  std::map<const string, const metric_t> ref_metrics = this->metricArrayToMap(reference->getMetrics());

  // get timing data
  std::map< const string, const double> prb_timers = this->timerDataToMap(sourcePrb->timers);
  std::map< const string, const double> ref_timers = this->timerDataToMap(sourceRef->timers);

  
//  if(rank == 0)
//  {
//    cout << "Have the following timing data for the problem:" << endl;
//    for(auto &i : prb_timers) cout << i.first <<" = " << i.second << endl;
//    
//    cout << "\nHave the following timing data for the reference:" << endl;
//    for(auto &i : ref_timers) cout << i.first <<" = " << i.second << endl;
//    cout << endl;
//
//  }
  
  // get all of the metrics to be tested
  std::queue<ParameterList> metrics = ComparisonHelper::getMetricsToCompare(metricsPlist);
  
  // run comparison
  bool all_tests_pass = true;
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
      if(!ComparisonHelper::metricComparisonTest(prb_metrics[metric_name],
                                                 ref_metrics[metric_name],
                                                 metrics.front(), msg))
      {
        all_tests_pass = false;
      }
      cout << msg.str() << endl;
      
    }
    else if(prb_timers.find(metric_name) != prb_timers.end() &&
            ref_timers.find(metric_name) != ref_timers.end())
    {
      if(rank == 0) cout << "\ncomparing timer: " << metric_name << endl;
      if(!ComparisonHelper::timerComparisonTest(prb_timers.at(metric_name),
                                                ref_timers.at(metric_name),
                                                metrics.front(), msg))
      {
        all_tests_pass = false;
      }
      cout << msg.str() << endl;
    }
    
    metrics.pop();
  }
  
  
  if(rank == 0)
  {
    if(all_tests_pass) cout << "\nAll metric/timer comparisons PASSED." << endl;
    else cout << "\nMetric/timer metric comparisons FAILED." << endl;
  }
}

std::map<const string, const ComparisonHelper::metric_t>
ComparisonHelper::metricArrayToMap(const ArrayRCP<const ComparisonHelper::metric_t> &metrics)
{
  typedef std::pair<const string,const metric_t> pair_t;
  std::map<const string, const metric_t> metric_map;
  ArrayRCP<const ComparisonHelper::metric_t>::size_type idx;
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
ComparisonHelper::metricComparisonTest(const Zoltan2::MetricValues<zscalar_t> & metric,
                                       const Zoltan2::MetricValues<zscalar_t> & ref_metric,
                                       const Teuchos::ParameterList & metricPlist,
                                       ostringstream &msg)
{
  // run a comparison of min and max agains a given metric
  // return an error message on failure
  bool pass = true;
  string test_name = metricPlist.name() + " test";
  double ref_value = ref_metric.getMaxImbalance()/ref_metric.getAvgImbalance();
  double value = metric.getMaxImbalance()/metric.getAvgImbalance();
  
  if (metricPlist.isParameter("lower"))
  {
    double min = metricPlist.get<double>("lower")*ref_value;
    
    if(value < min)
    {
      msg << test_name << " FAILED: Minimum imbalance per part, "
      << value << ", less than specified allowable minimum, " << min;
      pass = false;
    }
  }
  
  if(metricPlist.isParameter("upper" ) && pass != false) {
    
    double max = metricPlist.get<double>("upper") * ref_value;
    if (value > max)
    {
      msg << test_name << " FAILED: Maximum imbalance per part, "
      << value << ", greater than specified allowable maximum, " << max;
      pass = false;
    }
    
  }
  
  if(pass){
    msg << test_name << " PASSED.";
    pass = true;
  }
  
  return pass;
}

bool ComparisonHelper::timerComparisonTest(const double time,
                                           const double ref_time,
                                           const Teuchos::ParameterList & metricPlist,
                                           ostringstream &msg)
{
  // run a comparison of min and max agains a given metric
  // return an error message on failure
  bool pass = true;
  string test_name = metricPlist.name() + " test";
  if (metricPlist.isParameter("lower"))
  {
    double min = metricPlist.get<double>("lower")*ref_time;
    
    if(time < min)
    {
      msg << test_name << " FAILED: Minimum time, "
      << time <<
      ", less than specified allowable minimum time, " << min;
      pass = false;
    }
  }
  
  if(metricPlist.isParameter("upper" ) && pass != false) {
    
    double max = metricPlist.get<double>("upper") * ref_time;
    if (time > max)
    {
      msg << test_name << " FAILED: Maximum time, "
      << time <<
      ", greater than specified allowable maximum time, " << max;
      pass = false;
    }
    
  }
  
  if(pass){
    msg << test_name << " PASSED.";
    pass = true;
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

