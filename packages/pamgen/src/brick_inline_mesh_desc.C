// @HEADER
// ***************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// ***************************************************************************
// @HEADER

#include "inline_mesh_desc.h"
#include "brick_inline_mesh_desc.h"
#include "uns_inline_decomp.h"
#include <sstream>
#include "pamgen_fudges.h"
#include <math.h>

#include <Kokkos_Core.hpp>


namespace PAMGEN_NEVADA {

/*****************************************************************************/
void Brick_Inline_Mesh_Desc::calculateSize(long long & total_el_count,
									 long long & total_node_count,
									 long long & total_edge_count)
/*****************************************************************************/
{
  total_el_count =0L;
  total_node_count = 0L;
  total_edge_count = 0L;

  if(dimension == 3){
    for(long long k = 0; k < inline_b[2]; k ++){
      for(long long j = 0; j < inline_b[1]; j ++){
	for(long long i = 0; i < inline_b[0]; i ++){
	  total_el_count +=
	    (long long)interval[0][i]*
	    (long long)interval[1][j]*
	    (long long)interval[2][k];
	}
      }
    }
  }
  else{
    for(long long j = 0; j < inline_b[1]; j ++){
      for(long long i = 0; i < inline_b[0]; i ++){
	total_el_count +=
	  (long long)interval[0][i]*
	  (long long)interval[1][j];
      }
    }
  }

  long long nodes_temp[3];
  nodes_temp[0] = 0;
  nodes_temp[1] = 0;
  nodes_temp[2] = 0;

  for(long long i = 0; i < inline_b[0]; i ++){
    nodes_temp[0] += interval[0][i];
  }
  for(long long j = 0; j < inline_b[1]; j ++){
    nodes_temp[1] += interval[1][j];
  }
  if(dimension == 3){
    for(long long k = 0; k < inline_b[2]; k ++){
      nodes_temp[2] += interval[2][k];
    }
  }

  total_node_count = (nodes_temp[0]+1L)*(nodes_temp[1]+1L);

  if(dimension == 3){
    total_node_count *= nodes_temp[2] + 1L;
  }

  if(dimension == 3){
    total_edge_count =
      (nodes_temp[0]+1L)*(nodes_temp[1]+1L)*(nodes_temp[2]+0L)+
      (nodes_temp[0]+1L)*(nodes_temp[1]+0L)*(nodes_temp[2]+1L)+
      (nodes_temp[0]+0L)*(nodes_temp[1]+1L)*(nodes_temp[2]+1L);
  }
  else{
    total_edge_count =
      (nodes_temp[0]+0L)*(nodes_temp[1]+1L)+
      (nodes_temp[0]+1L)*(nodes_temp[1]+0L);
  }
}

/*****************************************************************************/
  std::string  Brick_Inline_Mesh_Desc::Calc_Intervals()
/*****************************************************************************/
{
  std::string error_string;

  for(long long axis = 0; axis < dimension; axis ++){
    for(long long i = 0; i < inline_b[axis]; i ++){
      if((first_size[axis][i] > 0.) && (last_size[axis][i] == 0.)){
	last_size[axis][i] = first_size[axis][i];
      }
      if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
	double xave = (first_size[axis][i]+last_size[axis][i])/2.;
	double k = a_lot_positive + (block_dist[axis][i]/xave);
	long long ktil =(long long)k;
	if(ktil < 1)ktil = 1;
	double delxtil = block_dist[axis][i]/(double)ktil;
	double s = xave-delxtil;
	interval[axis][i] = ktil;
	first_size[axis][i]-=s;
      last_size[axis][i] -=s;
      }
    }
  }
  return error_string;
}

/*****************************************************************************/
long long Brick_Inline_Mesh_Desc::Set_Up()
/*****************************************************************************/
{
  for(long long axis = 0 ; axis < 3; axis ++){
    a_inline_n[axis] = new long long [inline_b[axis]];
    c_inline_n[axis] = new long long [inline_b[axis]+1];
  }

  if(dimension != 3){
    nel_tot[2] = 1;
    a_inline_n[2][0] = 1;
    c_inline_n[2][0] = 0;
    c_inline_n[2][1] = 1;
  }

  for(long long axis = 0 ; axis < dimension; axis ++){
    for(long long i = 0; i < inline_b[axis]; i ++){
      a_inline_n[axis][i] = interval[axis][i];
      c_inline_n[axis][i] = 0;
      nel_tot[axis] += a_inline_n[axis][i];
      if(i)c_inline_n[axis][i] = c_inline_n[axis][i-1]+a_inline_n[axis][i-1];
      c_block_dist[axis][i] = inline_gmin[axis];//inner radius
      if(i)c_block_dist[axis][i] = c_block_dist[axis][i-1]+block_dist[axis][i-1];
    }
    c_inline_n[axis][inline_b[axis]] = c_inline_n[axis][inline_b[axis] - 1]+a_inline_n[axis][inline_b[axis] - 1];
    c_block_dist[axis][inline_b[axis]] = c_block_dist[axis][inline_b[axis] - 1]+block_dist[axis][inline_b[axis] - 1];
  }


  cum_block_totals = new long long[numBlocks()];
  els_in_block = new long long[numBlocks()];

  long long bl_ct = 0;
  for(long long k = 0; k < inline_b[2]; k ++){
    for(long long j = 0; j < inline_b[1]; j ++){
      for(long long i = 0; i < inline_b[0]; i ++){
	els_in_block[bl_ct] = a_inline_n[0][i]*a_inline_n[1][j]*a_inline_n[2][k];
	cum_block_totals[bl_ct]=0;
        if(bl_ct){
	  cum_block_totals[bl_ct] = cum_block_totals[bl_ct-1]+els_in_block[bl_ct-1];
        }
        bl_ct ++;
      }
    }
  }

  // this tolerance is dangerous
  if(fabs(c_block_dist[1][inline_b[1]]-360.0)< 1.0) periodic_j = true;


  for(long long axis = 0; axis < dimension; axis ++){
    inline_gmax[axis] = inline_gmin[axis];
    for(long long i = 0; i < inline_b[axis]; i ++){
      inline_gmax[axis] += block_dist[axis][i];
    }
  }

  return 0;
}

/****************************************************************************/
void Brick_Inline_Mesh_Desc::Populate_Coords(double * coords,
					    std::vector<long long> & global_node_vector,
					    std::map <long long, long long> & global_node_map,
					    long long num_nodes)

/****************************************************************************/
{

  // Create host and device views for coordinate data
  // Note: The input coords array uses structure-of-arrays layout (all x, then all y)
  // but we need array-of-structures layout for the 2D view (x,y for node 0, x,y for node 1, etc.)
  HostView2D<double> coords_host("coords_host", num_nodes, dimension);

  // Copy data from flattened structure-of-arrays to array-of-structures layout
  for (long long i = 0; i < num_nodes; i++) {
    for (long long axis = 0; axis < dimension; axis++) {
      coords_host(i, axis) = coords[i + axis * num_nodes];
    }
  }

  // Create device view and copy data
  auto coords_device = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), coords_host);
  Kokkos::deep_copy(coords_device, coords_host);

  // Create host views for node data first
  View1D<long long> global_node_vector_host("global_node_vector_host", global_node_vector.size());
  View1D<long long> global_node_map_keys_host("global_node_map_keys_host", global_node_map.size());
  View1D<long long> global_node_map_values_host("global_node_map_values_host", global_node_map.size());

  // Fill host views with data
  for (size_t i = 0; i < global_node_vector.size(); ++i) {
    global_node_vector_host(i) = global_node_vector[i];
  }

  size_t idx = 0;
  for (const auto& pair : global_node_map) {
    global_node_map_keys_host(idx) = pair.first;
    global_node_map_values_host(idx) = pair.second;
    idx++;
  }

  // Create device views and copy data
  auto global_node_vector_view = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), global_node_vector_host);
  auto global_node_map_keys = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), global_node_map_keys_host);
  auto global_node_map_values = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), global_node_map_values_host);

  // Execute device computation
  Populate_Coords_Device(coords_device, global_node_vector_view,
                        global_node_map_keys, global_node_map_values, num_nodes);

  // Copy results back to host memory
  Kokkos::deep_copy(coords_host, coords_device);

  // Copy back from array-of-structures layout to structure-of-arrays layout
  for (long long i = 0; i < num_nodes; i++) {
    for (long long axis = 0; axis < dimension; axis++) {
      coords[i + axis * num_nodes] = coords_host(i, axis);
    }
  }
}


/****************************************************************************/
void Brick_Inline_Mesh_Desc::Populate_Coords_Device(View2D<double> coords,
                                                         const View1D<long long> global_node_vector,
                                                         const View1D<long long> global_node_map_keys,
                                                         const View1D<long long> global_node_map_values,
                                                         long long num_nodes)
/****************************************************************************/
{
  Kokkos::parallel_for("BrickPopulateCoords", num_nodes,
    KOKKOS_LAMBDA(const size_t gnv) {
      long long the_node = global_node_vector(gnv);
      long long global_ind[3];

      // Calculate global indices
      global_ind[2] = the_node / knstride;
      global_ind[1] = (the_node - global_ind[2] * knstride) / jnstride;
      global_ind[0] = the_node - global_ind[2] * knstride - global_ind[1] * jnstride;

      // Find the local node using map lookup
      long long the_local_node = get_map_entry_device(global_node_map_keys,
                                                     global_node_map_values,
                                                     the_node,
                                                     global_node_map_keys.size());

      // Set coordinates
      for (long long axis = 0; axis < dimension; axis++) {
        coords(the_local_node, axis) = IJKcoors[axis][global_ind[axis]];
      }
    });

  // Ensure all device operations are complete
  Kokkos::fence();
}


}//end namespace PAMGEN_NEVADA