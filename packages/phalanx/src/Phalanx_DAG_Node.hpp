#ifndef PHALANX_DAG_NODE_HPP
#define PHALANX_DAG_NODE_HPP

#include "Phalanx_Evaluator.hpp"
#include <unordered_set>

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
    //! DAG discovery time.
    int discovery_time_;
    //! DAG final time.
    int final_time_;
    //! Out-edge node adjacency indices. Use a set to avoid duplicate edges.
    std::unordered_set<int> adjacencies_;

  public:
    DagNode(const int index,
	    const Teuchos::RCP<PHX::Evaluator<Traits>>& e) :
      index_(index), 
      c_(PHX::Color::WHITE),
      e_(e),
      discovery_time_(-1),
      final_time_(-1)
    {}
    DagNode(const DagNode<Traits>& ) = default;
    DagNode(DagNode<Traits>&& ) = default;
    DagNode<Traits>& operator=(const DagNode<Traits>& ) = default;
    int index() const {return index_;}
    Teuchos::RCP<const PHX::Evaluator<Traits>> get() const {return e_;}
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
  };

}

#endif
