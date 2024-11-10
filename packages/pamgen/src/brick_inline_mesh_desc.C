// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "inline_mesh_desc.h"
#include "brick_inline_mesh_desc.h"
#include "uns_inline_decomp.h"
#include <sstream>
#include "pamgen_fudges.h"
#include <math.h>


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
  long long global_ind[3];
  for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
    long long the_node = global_node_vector[gnv];
    global_ind[2] = the_node/knstride;
    global_ind[1] = (the_node-global_ind[2]*knstride)/jnstride;
    global_ind[0] = the_node - global_ind[2]*knstride-global_ind[1]*jnstride;
    long long the_local_node = get_map_entry(global_node_map,the_node);
    for(long long axis = 0; axis < dimension; axis ++){
      coords[the_local_node+axis*num_nodes]= IJKcoors[axis][global_ind[axis]];
    }
  }
}


}//end namespace PAMGEN_NEVADA
