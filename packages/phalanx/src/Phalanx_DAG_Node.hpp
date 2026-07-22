// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_DAG_NODE_HPP
#define PHALANX_DAG_NODE_HPP

#include "Phalanx_Evaluator.hpp"
#include <unordered_set>
#include <chrono>

namespace PHX {

  // class PHX::Evaluator;

  //! Colors for DAG DFS algorithm
  enum class Color {WHITE,GREY,BLACK};

  //! DAG Node wrapper for graph algorithms (DFS and topological sort).
  template<typename Traits>
  class DagNode {
    
    //! Evaluator index.
    int index_;
    //! Current node color.
    PHX::Color c_;
    //! Node data.
    Teuchos::RCP<PHX::Evaluator<Traits>> e_;
    //! DAG depth first search discovery time.
    int discovery_time_;
    //! DAG depth first search final time.
    int final_time_;
    //! Out-edge node adjacency indices. Use a set to avoid duplicate edges.
    std::unordered_set<int> adjacencies_;
    //! Execution time
    std::chrono::duration<double> exec_time_;
    //! Graph analysis node start time
    std::chrono::duration<double> t_start_;
    //! Graph analysis node finish time
    std::chrono::duration<double> t_finish_;

  public:
    DagNode(const int index,
	    const Teuchos::RCP<PHX::Evaluator<Traits>>& e) :
      index_(index), 
      c_(PHX::Color::WHITE),
      e_(e),
      discovery_time_(-1),
      final_time_(-1),
      exec_time_(0.0),
      t_start_(0.0),
      t_finish_(0.0)
    {}
    DagNode(const DagNode<Traits>& ) = default;
    DagNode(DagNode<Traits>&& ) = default;
    DagNode<Traits>& operator=(const DagNode<Traits>& ) = default;
    int index() const {return index_;}
    Teuchos::RCP<const PHX::Evaluator<Traits>> get() const {return e_;}
    Teuchos::RCP<PHX::Evaluator<Traits>> getNonConst() const {return e_;}
    void setColor(const PHX::Color& c) {c_ = c;}
    PHX::Color color() const {return c_;}
    int discoveryTime() const {return discovery_time_;}
    void setDiscoveryTime(int dt) {discovery_time_ = dt;}
    int finalTime() const {return final_time_;}
    void setFinalTime(int ft) {final_time_ = ft;}
    void resetDfsParams(const PHX::Color c = PHX::Color::WHITE) {
      c_ = c;
      discovery_time_ = -1;
      final_time_ = -1;
      adjacencies_.clear();
    }
    void addAdjacency(const int& node_index)
    {adjacencies_.insert(node_index);}
    const std::unordered_set<int>& adjacencies() const
    {return adjacencies_;}
    void setExecutionTime(const std::chrono::duration<double>& exec_time)
    {exec_time_ = exec_time;}
    void sumIntoExecutionTime(const std::chrono::duration<double>& exec_time)
    {exec_time_ += exec_time;}
    const std::chrono::duration<double>& executionTime() const
    {return exec_time_;}
    void setStartTime(const std::chrono::duration<double>& t)
    {t_start_ = t;}
    const std::chrono::duration<double>& startTime() const
    {return t_start_;}
    void setFinishTime(const std::chrono::duration<double>& t)
    {t_finish_ = t;}
    const std::chrono::duration<double>& finishTime() const
    {return t_finish_;}
  };

}

#endif
