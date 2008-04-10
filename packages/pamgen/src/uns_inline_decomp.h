#ifndef uns_inline_decompH
#define uns_inline_decompH

#include <list>
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
  int visits;
  int global_id;
  std::list < std::pair < int , Topo_Loc > > conn_connections;
  std::list <int>proc_neighbors;
};


class Partition{
public:
  
  Partition(int ls, int is, int js, int ks, int le, int ie, int je, int ke, InlineDecompositionType idt,int rcuts[]){
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

    for(int i = 0; i < 3; i ++)remaining_cuts[i] = rcuts[i];

    centroid = (((double)(ks+ke))/2.0)*(((double)(ks+ke))/2.0) + (((double)(js+je))/2.0)*(((double)(js+je))/2.0) + (((double)(is+ie))/2.0)*(((double)(is+ie))/2.0);
    numels =(ie-is)*(je-js)*(ke-ks);
    proc_id = -1;
    split_value = -1;
    split_direction = -1;
    unique_id = partition_count;
    partition_count++;
    inline_decomposition_type = idt;
  };

  Partition(int is, int js, int ks, int ie, int je, int ke, InlineDecompositionType idt,int rcuts[]){
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

    for(int i = 0; i < 3; i ++)remaining_cuts[i] = rcuts[i];

    centroid = (((double)(ks+ke))/2.0)*(((double)(ks+ke))/2.0) + (((double)(js+je))/2.0)*(((double)(js+je))/2.0) + (((double)(is+ie))/2.0)*(((double)(is+ie))/2.0);
    numels =(ie-is)*(je-js)*(ke-ks);
    proc_id = -1;
    split_value = -1;
    split_direction = -1;
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

  static int partition_count;

  void Processor_Partition(std::vector  < Partition *> & ,int inc_nels []);

  int Element_Proc(int []);
  Partition * high;
  Partition * low;
  int split_direction;
  int split_value;
  int lows[4];
  int highs[4];
  int remaining_cuts[3];
  double centroid;
  int numels;
  int proc_id;
  int unique_id;

  void print();
};


}// end namespace
#endif
