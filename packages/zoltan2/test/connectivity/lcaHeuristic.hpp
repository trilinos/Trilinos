#ifndef __lca_heuristic_hpp__
#define __lca_heuristic_hpp__
#include "edgeDS.hpp"
#include"graph.h"
#include<queue>

void bfs(graph* g, int root, int* parents, int* levels, edge_set& nonTreeEdges){
  std::queue<int> frontier;
  frontier.push(root);
  parents[root] = root;
  levels[root] = 0;
  while(!frontier.empty()){
    int curr_node = frontier.front();
    frontier.pop();
    
    int out_degree = out_degree(g, curr_node);
    int* outs = out_vertices(g, curr_node);
    for(int i = 0; i < out_degree; ++i){
      int neighbor = outs[i];
      if(parents[neighbor] < 0){
        parents[neighbor] = curr_node;
        levels[neighbor] = levels[curr_node] + 1;
        frontier.push(neighbor);
      } else if(neighbor != parents[curr_node]){
        edge e(curr_node, neighbor);
        e.validate();
        nonTreeEdges.insert(e);
      }
    }
  }
}

int lca(graph* g, int* levels, int* parents, int curr_node, int neighbor, edge_map<int>& visitedEdges){
  int p_curr_node = parents[curr_node];
  int p_neighbor = parents[neighbor];
  
  if(p_curr_node == neighbor){
    edge e(curr_node, neighbor);
    e.validate();
    visitedEdges[e] += 1;
    return -1;
  }

  if(p_neighbor == curr_node){
    edge e(curr_node, neighbor);
    e.validate();
    visitedEdges[e]+=1;
  }
 
  edge edge1(curr_node, p_curr_node);
  edge1.validate();
  edge edge2(neighbor, p_neighbor);
  edge2.validate();
  visitedEdges[edge1]+=1;
  visitedEdges[edge2]+=1;

  if(levels[p_curr_node] < levels[p_neighbor]){
    edge e(p_neighbor, parents[p_neighbor]);
    e.validate();
    visitedEdges[e] += 1;
    p_neighbor = parents[p_neighbor];
  } else if(levels[p_curr_node] > levels[p_neighbor]){
    edge e(p_curr_node, parents[p_curr_node]);
    e.validate();
    visitedEdges[e]+=1;
    p_curr_node = parents[p_curr_node];
  }

  while(p_curr_node != p_neighbor){
    edge edge1(p_curr_node, parents[p_curr_node]);
    edge edge2(p_neighbor, parents[p_neighbor]);
    edge1.validate();
    edge2.validate();
    visitedEdges[edge1] += 1;
    visitedEdges[edge2] += 1;
    p_curr_node = parents[p_curr_node];
    p_neighbor = parents[p_neighbor];
  }
  return p_curr_node;
}

void findPotentialArtPts(graph* g, int*& potential_art_pts, edge_map<int>& visited_edges){
  //create the arrays to hold the output of the bfs
  int* parents = new int[g->n];
  int* levels = new int[g->n];
  potential_art_pts = new int[g->n];
  for(int i = 0; i < g->n; i++){
    parents[i] = -1;
    levels[i] = -1;
    potential_art_pts[i] = 0;
  }

  edge_set nontree_edges;
  
  bfs(g,0,parents,levels,nontree_edges);
  
  for(edge_set::iterator itr = nontree_edges.begin(); itr != nontree_edges.end(); itr++){
    int endpt1= itr->getU();
    int endpt2 = itr->getV();
    int art_pt = lca(g, levels, parents, endpt1, endpt2, visited_edges);
    if(art_pt != -1) potential_art_pts[art_pt]++;
    edge e(endpt1,endpt2);
    e.validate();
    visited_edges[e]++;
  }

  for(edge_map<int>::iterator itr = visited_edges.begin(); itr != visited_edges.end(); itr++){
    if(itr->second == 0){
      potential_art_pts[itr->first.getU()]++;
      potential_art_pts[itr->first.getV()]++;
    }
  }
}

#endif
