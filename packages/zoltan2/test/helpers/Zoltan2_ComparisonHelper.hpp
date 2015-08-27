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

#include <sstream>
#include <string>
#include <map>
#include <iostream>

using Teuchos::Comm;
using Teuchos::RCP;

using std::cout;
using std::endl;
using std::string;
using std::map;
using std::pair;
using std::ostringstream;

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

class ComparisonSource
{
public:
  ~ComparisonSource()
  {
      if(adapter_kind == "XpetraCrsGraph")
        delete reinterpret_cast<xcrsGraph_t *>(adapter.getRawPtr())->getCoordinateInput();
      if(adapter_kind == "XpetraCrsMatrix")
        delete reinterpret_cast<xcrsMatrix_t *>(adapter.getRawPtr())->getCoordinateInput();
  }
  
  RCP<basic_problem_t> problem;
  RCP<basic_id_t> adapter;
  string problem_kind;
  string adapter_kind;
};

typedef std::pair<const string, RCP<const ComparisonSource> > pair_t;

class ComparisonHelper
{
public:
  
  void CompareSolutions(const string &p1,
                        const string &p2,
                        const RCP<const Comm<int> > &comm);
  
  void AddSource(const string &name, ComparisonSource * source);
  
  size_t getNumberOfSources() const
  {
    return this->sources.size();
  }
  
private:
  map<const string,RCP<const ComparisonSource> > sources;
  
  void ComparePartitionSolutions(const ComparisonSource * sourceA,
                                 const ComparisonSource * sourceB,
                                 const RCP<const Comm<int> > &comm);
  
  void CompareColoringSolutions(const ComparisonSource * sourceA,
                                const ComparisonSource * sourceB,
                                const RCP<const Comm<int> > &comm);
  
  void CompareOrderingSolutions(const ComparisonSource * sourceA,
                                const ComparisonSource * sourceB,
                                const RCP<const Comm<int> > &comm);
  
};


void ComparisonHelper::AddSource(const string &name, ComparisonSource * source)
{
  this->sources.insert(pair_t(name, RCP<ComparisonSource>(source)));
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
    
    if(A->problem_kind == "Partitioning")
    {
      this->ComparePartitionSolutions(A.getRawPtr(), B.getRawPtr(), comm);
      
    }else if(A->problem_kind == "Partitioning")
    {
      this->CompareColoringSolutions(A.getRawPtr(), B.getRawPtr(), comm);
      
    }else if(A->problem_kind == "Partitioning"){
      
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
  
  if(!failed)
  {
    // have some solutions lets compare them
    typedef Zoltan2::OrderingProblem<basic_id_t> ordering_problem_t;
    // have some solutions lets compare them
    if(ordering_problem_t * problem_a = reinterpret_cast<ordering_problem_t *>(sourceA->problem.getRawPtr()))
    {
      if(ordering_problem_t * problem_b = reinterpret_cast<ordering_problem_t *>(sourceB->problem.getRawPtr()))
      {
        
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

