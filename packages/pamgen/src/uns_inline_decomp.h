// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef uns_inline_decompH
#define uns_inline_decompH

#include <list>
#include <cstdlib>
#include <vector>
#include "topology_enum.h"
#include "inline_geometries.h"

namespace PAMGEN_NEVADA {

class Tel{
public:
  Tel(){visits = 0;
  global_id = -1;
  real_element = false;
  periodic_minj = false;
  periodic_maxj = false;};
  ~Tel(){};
  bool real_element;
  bool periodic_minj;
  bool periodic_maxj;
  long long visits;
  long long global_id;
  std::list < std::pair < long long , Topo_Loc > > conn_connections;
  std::list <long long>proc_neighbors;
};


class Partition{
public:
  
  Partition(long long ls, long long is, long long js, long long ks, long long le, long long ie, long long je, long long ke, InlineDecompositionType idt,long long rcuts[]){
    high = NULL;
    low = NULL;
    lows[0] = is;
    lows[1] = js;
    lows[2] = ks;
    lows[3] = ls;
    
    highs[0] = ie;
    highs[1] = je;
    highs[2] = ke;
    highs[3] = le;

    for(long long i = 0; i < 3; i ++)remaining_cuts[i] = rcuts[i];

    centroid = (((double)(ks+ke))/2.0)*(((double)(ks+ke))/2.0) + (((double)(js+je))/2.0)*(((double)(js+je))/2.0) + (((double)(is+ie))/2.0)*(((double)(is+ie))/2.0);
    numels =(ie-is)*(je-js)*(ke-ks);
    proc_id = -1;
    split_value = -1;
    split_direction = 0;
    unique_id = partition_count;
    partition_count++;
    inline_decomposition_type = idt;
  };

  Partition(long long is, long long js, long long ks, long long ie, long long je, long long ke, InlineDecompositionType idt,long long rcuts[]){
    high = NULL;
    low = NULL;
    lows[0] = is;
    lows[1] = js;
    lows[2] = ks;
    lows[3] = 0;
    
    highs[0] = ie;
    highs[1] = je;
    highs[2] = ke;
    highs[3] = 1;

    for(long long i = 0; i < 3; i ++)remaining_cuts[i] = rcuts[i];

    centroid = (((double)(ks+ke))/2.0)*(((double)(ks+ke))/2.0) + (((double)(js+je))/2.0)*(((double)(js+je))/2.0) + (((double)(is+ie))/2.0)*(((double)(is+ie))/2.0);
    numels =(ie-is)*(je-js)*(ke-ks);
    proc_id = -1;
    split_value = -1;
    split_direction = 0;
    unique_id = partition_count;
    partition_count++;
    inline_decomposition_type = idt;
  };
  
  ~Partition(){
  }

  void empty(){
    if(high)high->empty();
    if(low)low->empty();
    delete high;
    high = NULL;
    delete low;
    low = NULL;
  };

  InlineDecompositionType inline_decomposition_type;

  static long long partition_count;

  void Processor_Partition(std::vector  < Partition *> & ,long long inc_nels []);

  long long Element_Proc(long long []);
  Partition * high;
  Partition * low;
  unsigned split_direction;
  long long split_value;
  long long lows[4];
  long long highs[4];
  long long remaining_cuts[3];
  double centroid;
  long long numels;
  long long proc_id;
  long long unique_id;

  void print();
};


}// end namespace
#endif
