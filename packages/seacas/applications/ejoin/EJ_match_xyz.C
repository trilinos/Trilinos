// Copyright(C) 2010 Sandia Corporation.
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#include <cmath>
#include <float.h>
#include "EJ_match_xyz.h"
#include "EJ_vector3d.h"
#include "EJ_index_sort.h"
#include "smart_assert.h"
#include "EJ_mapping.h"
#include <Ioss_SubSystem.h>

namespace {
  template <typename INT>
  void do_matching(std::vector<INT> &i_inrange, const RealVector &i_coord, size_t i_offset,
		   std::vector<INT> &j_inrange, const RealVector &j_coord, size_t j_offset,
		   double epsilon, int XYZ, std::vector<INT> &local_node_map);

  double max3(double x, double y, double z)
  {
    double max = x;
    if (y > max)
      max = y;
    if (z > max)
      max = z;
    return max;
  }
  
  void find_range(std::vector<double> &coord, vector3d &min, vector3d &max)
  {
    if (!coord.empty()) {
      min.set(coord[0], coord[1], coord[2]);
      max = min;
      for (size_t i=3; i < coord.size(); i+=3) {
	if (min.x > coord[i+0]) min.x = coord[i+0];
	if (min.y > coord[i+1]) min.y = coord[i+1];
	if (min.z > coord[i+2]) min.z = coord[i+2];

	if (max.x < coord[i+0]) max.x = coord[i+0];
	if (max.y < coord[i+1]) max.y = coord[i+1];
	if (max.z < coord[i+2]) max.z = coord[i+2];
      }
    } else {
      min.set(0,0,0);
      max = min;
    }
  }

  template <typename INT>
  void find_in_range(std::vector<double> coord, vector3d &min, vector3d &max, std::vector<INT> &in_range)
  {
    if (!coord.empty()) {
      for (size_t i=0; i < coord.size(); i+=3) {
	if (coord[i+0] > min.x  && coord[i+0] < max.x &&
	    coord[i+1] > min.y  && coord[i+1] < max.y &&
	    coord[i+2] > min.z  && coord[i+2] < max.z) {
	  in_range.push_back(i/3);
	}
      }
    }
  }
}

template <typename INT>
void match_node_xyz(RegionVector &part_mesh,
		    double tolerance,
		    std::vector<INT> &global_node_map, std::vector<INT> &local_node_map)
{
  // See if any omitted element blocks...
  bool has_omissions = false;
  for (size_t p = 0; p < part_mesh.size(); p++) {
    if (part_mesh[p]->get_property("block_omission_count").get_int() > 0) {
      has_omissions = true;
      break;
    }
  }

  if (!has_omissions) {
    for (size_t i=0; i < local_node_map.size(); i++) {
      local_node_map[i] = i;
    }
  } else {
    std::vector<INT> dummy;
    eliminate_omitted_nodes(part_mesh, dummy, local_node_map);

    // The local_node_map is not quite in the correct format after the
    // call to 'eliminate_omitted_nodes'.  We need all non-omitted
    // nodes to have local_node_map[i] == i.
    for (size_t i=0; i < local_node_map.size(); i++) {
      if (local_node_map[i] >= 0) local_node_map[i] = i;
    }
  }
  
  size_t part_count = part_mesh.size();
  enum {X=0, Y=1, Z=2};

  for (size_t ip=0; ip < part_count; ip++) {
    vector3d i_max;
    vector3d i_min;
    std::vector<double> i_coord;
    Ioss::NodeBlock *inb = part_mesh[ip]->get_node_blocks()[0];
    inb->get_field_data("mesh_model_coordinates", i_coord);
    find_range(i_coord, i_min, i_max);
    
    size_t i_offset = part_mesh[ip]->get_property("node_offset").get_int();

    for (size_t jp=ip+1; jp < part_count; jp++) {
      vector3d j_max;
      vector3d j_min;
      std::vector<double> j_coord;
      Ioss::NodeBlock *jnb = part_mesh[jp]->get_node_blocks()[0];
      jnb->get_field_data("mesh_model_coordinates", j_coord);
      find_range(j_coord, j_min, j_max);

      size_t j_offset = part_mesh[jp]->get_property("node_offset").get_int();

      // See if the ranges overlap...
      vector3d max;
      vector3d min;
      max.x = std::min(i_max.x, j_max.x);
      max.y = std::min(i_max.y, j_max.y);
      max.z = std::min(i_max.z, j_max.z);

      min.x = std::max(i_min.x, j_min.x);
      min.y = std::max(i_min.y, j_min.y);
      min.z = std::max(i_min.z, j_min.z);

      double delta[3];
      int XYZ = X;
      delta[XYZ] = max.x - min.x;
      delta[Y] = max.y - min.y;
      if (delta[Y] > delta[XYZ])
	XYZ = Y;
      delta[Z] = max.z - min.z;
      if (delta[Z] > delta[XYZ])
	XYZ = Z;

      double epsilon = (delta[X] + delta[Y] + delta[Z]) / 1.0e3;
      if (epsilon < 0.0) {
	std::cout << "Parts " << ip << " and " << jp << " do not overlap.\n";
	continue;
      }

      min -= epsilon;
      max += epsilon;

      if (tolerance >= 0.0) epsilon = tolerance;

      std::vector<INT> j_inrange;
      std::vector<INT> i_inrange;

      find_in_range(j_coord, min, max, j_inrange);
      find_in_range(i_coord, min, max, i_inrange);

      // Index sort all nodes on the coordinate range with the maximum delta.
      index_coord_sort(i_coord, i_inrange, XYZ);
      index_coord_sort(j_coord, j_inrange, XYZ);

      if (i_inrange.size() < j_inrange.size()) {
	do_matching(i_inrange, i_coord, i_offset,
		    j_inrange, j_coord, j_offset,
		    epsilon, XYZ, local_node_map);
      } else {
	do_matching(j_inrange, j_coord, j_offset,
		    i_inrange, i_coord, i_offset,
		    epsilon, XYZ, local_node_map);
      }
    }
  }

  // Build the global and local maps...
  size_t j = 1;
  for (size_t i=0; i < local_node_map.size(); i++) {
    if (local_node_map[i] == (INT)i) {
      global_node_map.push_back(j);
      local_node_map[i] = j-1;
      j++;
    } else if (local_node_map[i] >= 0) {
      local_node_map[i] = local_node_map[local_node_map[i]];
    }
  }
}
template void match_node_xyz(RegionVector &part_mesh,
			     double tolerance,
			     std::vector<int> &global_node_map, std::vector<int> &local_node_map);
template void match_node_xyz(RegionVector &part_mesh,
			     double tolerance,
			     std::vector<int64_t> &global_node_map, std::vector<int64_t> &local_node_map);
  
namespace {
  template <typename INT>
  void do_matching(std::vector<INT> &i_inrange, const RealVector &i_coord, size_t i_offset,
		   std::vector<INT> &j_inrange, const RealVector &j_coord, size_t j_offset,
		   double epsilon, int XYZ, std::vector<INT> &local_node_map)
  {
    INT j2beg = 0;
      
    INT match = 0;
    INT compare = 0;

    double g_dismin =  FLT_MAX;
    double dismax = -FLT_MAX;

    for (size_t i=0; i < i_inrange.size(); i++) {
      INT ii = i_inrange[i];
      if (local_node_map[ii+i_offset] < 0)
	continue;

      double dismin =  FLT_MAX;
      double dmin = FLT_MAX;
      INT node_dmin = -1;

      for (size_t j = j2beg; j < j_inrange.size(); j++) {
	compare++;
	INT  jj = j_inrange[j];
	if (jj < 0 || local_node_map[jj+j_offset] < 0)
	  continue;

	if (i_coord[3*ii+XYZ]-epsilon > j_coord[3*jj+XYZ]) {
	  j2beg = j;
	  continue;
	}

	//... Since we are sorted on coordinate X|Y|Z,
	//    if set 'j' X|Y|Z greater than set 'i' X|Y|Z+eps, go to next 'i' X1|Y1|Z1 coord.
	if (j_coord[3*jj+XYZ]-epsilon > i_coord[3*ii+XYZ])
	  break;

	double distance = max3(std::fabs(j_coord[3*jj+0] - i_coord[3*ii+0]),
			       std::fabs(j_coord[3*jj+1] - i_coord[3*ii+1]),
			       std::fabs(j_coord[3*jj+2] - i_coord[3*ii+2]));
	  
	if (float(distance) <= epsilon) {
	  if (distance < dmin) {
	    dmin = distance;
	    node_dmin = j;
	  }
	} else {
	  if (distance < dismin) dismin = distance;
	}
	if (distance == 0.0)
	  break;
      }

      if (dmin <= epsilon) {
	INT jnod = j_inrange[node_dmin] + j_offset;
	INT inod = ii+i_offset;
	match++;
	if (dmin > dismax) dismax = dmin;
	j_inrange[node_dmin] *= -1;
	SMART_ASSERT(jnod < (INT)local_node_map.size());
	if (inod < jnod)
	  local_node_map[jnod] = inod;
	else
	  local_node_map[inod] = jnod;
      } else {
	if (dismin < g_dismin)
	  g_dismin = dismin;
      }
    }
    std::cout << "\nNumber of nodes matched                   = " << match << "\n";
    std::cout << "Number of comparisons                     = " << compare  << "\n";
    std::cout << "Tolerance used for matching               = " << epsilon << "\n";
    if (dismax > -FLT_MAX)
      std::cout << "Maximum distance between matched nodes    = " << dismax << "\n";
    if (g_dismin < FLT_MAX)
      std::cout << "Minimum distance between nonmatched nodes = " << g_dismin << "\n";
    std::cout << "\n";
  }
}
