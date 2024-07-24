// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <math.h>
#include "uns_inline_decomp.h"
#include "inline_mesh_desc.h"
#include "legacy_inline_mesh_desc.h"

namespace PAMGEN_NEVADA {

/*****************************************************************************/
void Legacy_Inline_Mesh_Desc::calculateSize(long long & total_el_count, 
					    long long & total_node_count, 
					    long long & total_edge_count)
  /*****************************************************************************/
{
  total_el_count = 
    (long long) inline_b[0] *
    (long long) inline_b[1] *
    (long long) inline_b[2] * 
    (long long) inline_n[0] * 
    (long long) inline_n[1] * 
    (long long) inline_n[2];

  if(dimension == 3){
  total_node_count = 
    ((long long)inline_b[0] * (long long)inline_n[0] +1L)*
    ((long long)inline_b[1] * (long long)inline_n[1] +1L)*
    ((long long)inline_b[2] * (long long)inline_n[2] +1L);
  }
  else{
  total_node_count = 
    ((long long)inline_b[0] * (long long)inline_n[0] +1L)*
    ((long long)inline_b[1] * (long long)inline_n[1] +1L);
  }

  if(dimension == 3){
  total_edge_count = 
    ((long long)inline_b[0] * (long long)inline_n[0] +1L)* 
    ((long long)inline_b[1] * (long long)inline_n[1])    *
    ((long long)inline_b[2] * (long long)inline_n[2] +1L)
    +
    ((long long)inline_b[1] * (long long)inline_n[1] +1L)* 
    ((long long)inline_b[0] * (long long)inline_n[0])    *
    ((long long)inline_b[2] * (long long)inline_n[2] +1L)
    +
    ((long long)inline_b[0] * (long long)inline_n[0] +1L)* 
    ((long long)inline_b[2] * (long long)inline_n[2])    *
    ((long long)inline_b[1] * (long long)inline_n[1] +1L);
  }
  else{
 total_edge_count = 
    ((long long)inline_b[0] * (long long)inline_n[0] +1L)* 
    ((long long)inline_b[1] * (long long)inline_n[1]) 
    +
    ((long long)inline_b[0] * (long long)inline_n[0])* 
    ((long long)inline_b[1] * (long long)inline_n[1] +1L)  ;
  }
}


/*****************************************************************************/
long long Legacy_Inline_Mesh_Desc::Set_Up()
/*****************************************************************************/
{
  a_inline_n[0] = new  long long [inline_b[0]];
  a_inline_n[1] = new  long long [inline_b[1]];
  a_inline_n[2] = new  long long [inline_b[2]];
  c_inline_n[0] = new  long long [inline_b[0]+1];
  c_inline_n[1] = new  long long [inline_b[1]+1];
  c_inline_n[2] = new  long long [inline_b[2]+1];

  c_block_dist[0] = new double[inline_b[0] + 1];
  c_block_dist[1] = new double[inline_b[1] + 1];
  c_block_dist[2] = new double[inline_b[2] + 1];
  
  block_dist[0] = new double[inline_b[0]];
  block_dist[1] = new double[inline_b[1]];
  block_dist[2] = new double[inline_b[2]];


  for(long long i = 0; i < inline_b[0]; i ++){
    a_inline_n[0][i] = inline_n[0];
    c_inline_n[0][i] = 0;
    nel_tot[0] += a_inline_n[0][i];
    if(i)c_inline_n[0][i] = c_inline_n[0][i-1]+a_inline_n[0][i-1];
    block_dist[0][i] =(inline_gmax[0] - inline_gmin[0])/(double)inline_b[0];
    c_block_dist[0][i] = inline_gmin[0];//inner radius
    if(i)c_block_dist[0][i] = c_block_dist[0][i-1]+block_dist[0][i-1];
  }
  c_inline_n[0][inline_b[0]] = c_inline_n[0][inline_b[0] - 1]+a_inline_n[0][inline_b[0] - 1];
  c_block_dist[0][inline_b[0]] = c_block_dist[0][inline_b[0] - 1]+block_dist[0][inline_b[0] - 1];

  for(long long i = 0; i < inline_b[1]; i ++){
    a_inline_n[1][i] = inline_n[1];
    c_inline_n[1][i] = 0;
    nel_tot[1] += a_inline_n[1][i];
    if(i)c_inline_n[1][i] = c_inline_n[1][i-1]+a_inline_n[1][i-1];
    block_dist[1][i] =(inline_gmax[1] - inline_gmin[1])/(double)inline_b[1];
    c_block_dist[1][i] = inline_gmin[1];
    if(i)c_block_dist[1][i] = c_block_dist[1][i-1]+block_dist[1][i-1];
  }
  c_inline_n[1][inline_b[1]] = c_inline_n[1][inline_b[1]-1]+a_inline_n[1][inline_b[1]-1];
  c_block_dist[1][inline_b[1]] = c_block_dist[1][inline_b[1] - 1]+block_dist[1][inline_b[1] - 1];


  for(long long i = 0; i < inline_b[2]; i ++){
    a_inline_n[2][i] = inline_n[2];
    c_inline_n[2][i] = 0;
    nel_tot[2] += a_inline_n[2][i];
    if(i)c_inline_n[2][i] = c_inline_n[2][i-1]+a_inline_n[2][i-1];
    block_dist[2][i] =(inline_gmax[2] - inline_gmin[2])/(double)inline_b[2];
    c_block_dist[2][i] = inline_gmin[2];
    if(i)c_block_dist[2][i] = c_block_dist[2][i-1]+block_dist[2][i-1];
  }
  c_inline_n[2][inline_b[2]] = c_inline_n[2][inline_b[2]-1]+a_inline_n[2][inline_b[2]-1];
  c_block_dist[2][inline_b[2]] = c_block_dist[2][inline_b[2] - 1]+block_dist[2][inline_b[2] - 1];

  cum_block_totals = new long long[inline_b[0]*inline_b[1]*inline_b[2]];
  els_in_block = new long long[inline_b[0]*inline_b[1]*inline_b[2]];

   long long bl_ct = 0;
  for(long long k = 0; k < inline_b[2]; k ++){
    for(long long j = 0; j < inline_b[1]; j ++){
      for(long long i = 0; i < inline_b[0]; i ++){
        els_in_block[bl_ct] = a_inline_n[0][i]*a_inline_n[1][j]*a_inline_n[2][k];
        cum_block_totals[bl_ct]=0;
        if(bl_ct){
	  cum_block_totals[bl_ct] = cum_block_totals[bl_ct-1]+els_in_block[bl_ct-1];;
        }
        bl_ct ++;
      }
    }
  }
  return 0;
}


/****************************************************************************/
long long Legacy_Inline_Mesh_Desc::Calc_Coord_Vectors()
/****************************************************************************/
{
   long long nnx = nel_tot[0]+1;
   long long nny = nel_tot[1]+1;
  double xdelta = inline_gmax[0]-inline_gmin[0];
  double ydelta = inline_gmax[1]-inline_gmin[1];
  IJKcoors[0] = new double[nnx];
  IJKcoors[1] = new double[nny];

   long long nct = 0;
  for(long long i = 0; i < inline_b[0]; i ++){
    for(long long j = 0; j < a_inline_n[0][i]; j ++){
      IJKcoors[0][nct] = c_block_dist[0][i]+j*block_dist[0][i]/(double)a_inline_n[0][i];
      IJKcoors[0][nct+1] = c_block_dist[0][i]+(j+1)*block_dist[0][i]/(double)a_inline_n[0][i];
      nct ++;
    }
  }

  nct = 0;
  for(long long i = 0; i < inline_b[1]; i ++){
    for(long long j = 0; j < a_inline_n[1][i]; j ++){
      IJKcoors[1][nct] = c_block_dist[1][i]+j*block_dist[1][i]/(double)a_inline_n[1][i];
      IJKcoors[1][nct+1] = c_block_dist[1][i]+(j+1)*block_dist[1][i]/(double)a_inline_n[1][i];
      nct ++;
    }
  }

  //   for(long long i = 0; i < nnx; i++)IJKcoors[0][i] = inline_gmin[0]+(double)i*dx;
  //   for(long long i = 0; i < nny; i++)IJKcoors[1][i] = inline_gmin[1]+(double)i*dy;
  if(Element_Density_Functions[0])Element_Density_Functions[0]->Integrate(inline_gmin[0],inline_gmax[0], error_stream);
  if(!error_stream.str().empty()){return 1;}
  if(Element_Density_Functions[1])Element_Density_Functions[1]->Integrate(inline_gmin[1],inline_gmax[1], error_stream);
  if(!error_stream.str().empty()){return 1;}
  if(Element_Density_Functions[0]){
    for(long long ict = 0; ict < nnx; ict ++){
      double factor = (IJKcoors[0][ict]-inline_gmin[0])/xdelta;
      double interpolant =  Element_Density_Functions[0]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      double new_coord = inline_gmin[0]+interpolant*xdelta;
      IJKcoors[0][ict] = new_coord;
    }
  }
  if(Element_Density_Functions[1]){
    for(long long ict = 0; ict < nny; ict ++){
      double factor = (IJKcoors[1][ict]-inline_gmin[1])/ydelta;
      double interpolant =  Element_Density_Functions[1]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      double new_coord = inline_gmin[1]+interpolant*ydelta;
      IJKcoors[1][ict] = new_coord;
    }
  }
  if(dimension == 3){
   long long nnz = nel_tot[2]+1;
  double zdelta = inline_gmax[2]-inline_gmin[2];
  IJKcoors[2] = new double[nnz];

  nct = 0;
  for(int i = 0; i < inline_b[2]; i ++){
    for(long long j = 0; j < a_inline_n[2][i]; j ++){
      IJKcoors[2][nct] = c_block_dist[2][i]+j*block_dist[2][i]/(double)a_inline_n[2][i];
      IJKcoors[2][nct+1] = c_block_dist[2][i]+(j+1)*block_dist[2][i]/(double)a_inline_n[2][i];
      nct ++;
    }
  }

  //   for(long long i = 0; i < nnz; i++)IJKcoors[2][i] = inline_gmin[2]+(double)i*dz;
  if(Element_Density_Functions[2])Element_Density_Functions[2]->Integrate(inline_gmin[2],inline_gmax[2],error_stream);
  if(!error_stream.str().empty()){return 1;}
  if(Element_Density_Functions[2]){
    for(long long ict = 0; ict < nnz; ict ++){
      double factor = (IJKcoors[2][ict]-inline_gmin[2])/zdelta;
      double interpolant =  Element_Density_Functions[2]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      double new_coord = inline_gmin[2]+interpolant*zdelta;
      IJKcoors[2][ict] = new_coord;
    }
  }
  }
  return 0;
}

/****************************************************************************/
void Cylindrical_Inline_Mesh_Desc::Populate_Coords(double * coords,   
				    std::vector<long long> & global_node_vector,                             
				    std::map <long long, long long> & global_node_map,
				     long long num_nodes)

/****************************************************************************/
{
  if(dimension == 3){
  double deg_to_rad = M_PI/180.0;
    for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
       long long the_node = global_node_vector[gnv];
       long long global_k = the_node/knstride;
       long long global_j = (the_node-global_k*knstride)/jnstride;
       long long global_i = the_node - global_k*knstride-global_j*jnstride;
       long long the_local_node = get_map_entry(global_node_map,the_node);
      coords[the_local_node+0*num_nodes]= IJKcoors[0][global_i]*cos(IJKcoors[1][global_j]*deg_to_rad);
      coords[the_local_node+1*num_nodes]= IJKcoors[0][global_i]*sin(IJKcoors[1][global_j]*deg_to_rad);
      coords[the_local_node+2*num_nodes]= IJKcoors[2][global_k];
    }
  }
}

/****************************************************************************/
void Spherical_Inline_Mesh_Desc::Populate_Coords(double * coords,   
				    std::vector<long long> & global_node_vector,                             
				    std::map <long long, long long> & global_node_map,
				     long long num_nodes)

/****************************************************************************/
{
  double deg_to_rad = M_PI/180.0;

  for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
     long long the_node = global_node_vector[gnv];
     long long global_k = the_node/knstride;
     long long global_j = (the_node-global_k*knstride)/jnstride;
     long long global_i = the_node - global_k*knstride-global_j*jnstride;
     long long the_local_node = get_map_entry(global_node_map,the_node);
    
    if(dimension == 2){
    coords[the_local_node+0*num_nodes]= IJKcoors[0][global_i]*cos(IJKcoors[1][global_j]*deg_to_rad);
    coords[the_local_node+1*num_nodes]= IJKcoors[0][global_i]*sin(IJKcoors[1][global_j]*deg_to_rad);
    }
    else{
    coords[the_local_node+0*num_nodes]= IJKcoors[0][global_i]*cos(IJKcoors[1][global_j]*deg_to_rad);
    coords[the_local_node+1*num_nodes]= IJKcoors[0][global_i]*sin(IJKcoors[1][global_j]*deg_to_rad)*cos(IJKcoors[2][global_k]*deg_to_rad);
    coords[the_local_node+2*num_nodes]= IJKcoors[0][global_i]*sin(IJKcoors[1][global_j]*deg_to_rad)*sin(IJKcoors[2][global_k]*deg_to_rad);
    }
  }
}

/****************************************************************************/
void Cartesian_Inline_Mesh_Desc::Populate_Coords(double * coords,   
				    std::vector<long long> & global_node_vector,                             
				    std::map <long long, long long> & global_node_map,
				     long long num_nodes)

/****************************************************************************/
{
  for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
     long long the_node = global_node_vector[gnv];
     long long global_k = the_node/knstride;
     long long global_j = (the_node-global_k*knstride)/jnstride;
     long long global_i = the_node - global_k*knstride-global_j*jnstride;
     long long the_local_node = get_map_entry(global_node_map,the_node);
    coords[the_local_node+0*num_nodes] = IJKcoors[0][global_i];
    coords[the_local_node+1*num_nodes] = IJKcoors[1][global_j];
  if(dimension == 3){
    coords[the_local_node+2*num_nodes] = IJKcoors[2][global_k];
  }
  }
}
}//end namespace PAMGEN_NEVADA
