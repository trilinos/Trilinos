// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef legacy_inline_mesh_descH
#define legacy_inline_mesh_descH
#include "inline_mesh_desc.h"


namespace PAMGEN_NEVADA {

class Legacy_Inline_Mesh_Desc : public Inline_Mesh_Desc
{
public:
  Legacy_Inline_Mesh_Desc(){
  };

  Legacy_Inline_Mesh_Desc(long long dim){
    dimension = dim;
  };


  virtual ~Legacy_Inline_Mesh_Desc(){};
  virtual void calculateSize(long long & total_el_count, 
			     long long & total_node_count, 
			     long long & total_edge_count);
  virtual long long Set_Up();
  virtual long long Calc_Coord_Vectors();
};


class Cartesian_Inline_Mesh_Desc : public Legacy_Inline_Mesh_Desc
{
public:
  Cartesian_Inline_Mesh_Desc(long long dim){dimension = dim;};
  virtual ~Cartesian_Inline_Mesh_Desc(){};
  virtual void Populate_Coords(double * coords,   
		       std::vector<long long> & global_node_vector, 
		       std::map <long long, long long> & global_node_map,
		       long long num_nodes);
};

class Cylindrical_Inline_Mesh_Desc : public Legacy_Inline_Mesh_Desc
{
public:
  Cylindrical_Inline_Mesh_Desc(long long dim){dimension = dim;};
  virtual ~Cylindrical_Inline_Mesh_Desc(){};
  virtual void Populate_Coords(double * coords,   
		       std::vector<long long> & global_node_vector, 
		       std::map <long long, long long> & global_node_map,
		       long long num_nodes);
};

class Spherical_Inline_Mesh_Desc : public Legacy_Inline_Mesh_Desc
{
public:
  Spherical_Inline_Mesh_Desc(long long dim){dimension = dim;};
  virtual ~Spherical_Inline_Mesh_Desc(){};
  virtual void Populate_Coords(double * coords,   
		       std::vector<long long> & global_node_vector, 
		       std::map <long long, long long> & global_node_map,
		       long long num_nodes);
};
} //end namespace
#endif
