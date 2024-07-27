// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef pamgen_bc_specificationH
#define pamgen_bc_specificationH

#include "PamgenStrLoopLimits.h"
#include "topology_enum.h"
#include <vector>


class PG_BC_Spec_Loc{
public:
  long long block_id;
  Topo_Loc location;
  PAMGEN_NEVADA::LoopLimits limits;
  bool block_boundary_set;
};

class PG_BC_Specification{
public:
  PG_BC_Specification(long long the_id,Topo_Loc the_loc, bool is_block_boundary, long long the_block_id){
    PG_BC_Spec_Loc sl;
    
    id = the_id;
    sl.location = the_loc;
    sl.block_boundary_set = is_block_boundary;
    sl.block_id = the_block_id;
    the_locs.push_back(sl);
  };
  
  void addEntry(Topo_Loc the_loc, bool is_block_boundary, long long the_block_id){
    PG_BC_Spec_Loc sl;
    sl.location = the_loc;
    sl.block_boundary_set = is_block_boundary;
    sl.block_id = the_block_id;
    the_locs.push_back(sl);
  };

  long long id;

  ~PG_BC_Specification(){};
  std::vector<PG_BC_Spec_Loc>the_locs;

};


#endif
