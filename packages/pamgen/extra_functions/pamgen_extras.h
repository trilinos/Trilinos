// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <list>
#include <vector>
#include <string>
#include <assert.h>
class topo_entity{
public:
  topo_entity(){
    local_id = -1;
    global_id = -1;
    owned = true;
  };
  void add_node(long long the_val,long long * global_nids){
    local_node_ids.push_back(the_val);
    sorted_local_node_ids.push_back(the_val);
    sorted_global_node_ids.push_back(global_nids[the_val-1]);
  }
  void sort(){
    sorted_local_node_ids.sort();
    sorted_global_node_ids.sort();
  }
  ~topo_entity(){};
  std::list <long long > local_node_ids;
  std::list <long long > sorted_local_node_ids;
  std::list <long long > sorted_global_node_ids;
  long long local_id;
  long long global_id;
  bool owned;
};

/*******************************************************************************/
inline bool compare_sorted_global_node_ids ( topo_entity* const x,  topo_entity* const y )
/*******************************************************************************/
{
  assert(x->sorted_global_node_ids.size() == y->sorted_global_node_ids.size());
  if(x->sorted_global_node_ids < y->sorted_global_node_ids)return true;    
  return false;
}

void calc_global_ids(std::vector < topo_entity * > eof_vec,
		     long long **comm_node_ids,
		     long long * node_comm_proc_ids,
		     long long * node_cmap_node_cnts,
		     int num_node_comm_maps,
		     int rank,
		     std::string fname_string);

void calc_global_node_ids(long long * globalNodeIds,
			  bool * nodeIsOwned,
			  long long numNodes,
			  long long num_node_comm_maps,
			  long long * node_cmap_node_cnts,
			  long long * node_comm_proc_ids,
			  long long * * comm_node_ids,
			  int rank);
