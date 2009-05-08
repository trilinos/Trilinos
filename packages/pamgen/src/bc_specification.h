#ifndef bc_specificationH
#define bc_specificationH

#include "StrLoopLimits.h"
#include "topology_enum.h"


class PG_BC_Specification{
public:
  PG_BC_Specification(long long the_id,Topo_Loc the_loc, bool is_block_boundary, long long the_block_id){
  id = the_id;
  location = the_loc;
  block_boundary_set = is_block_boundary;
  block_id = the_block_id;
  };

  ~PG_BC_Specification(){};

  long long id;
  long long block_id;
  Topo_Loc location;
  PAMGEN_NEVADA::LoopLimits limits;
  bool block_boundary_set;
};


#endif
