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
  a_inline_nx = new  long long [inline_b[0]];
  a_inline_ny = new  long long [inline_b[1]];
  a_inline_nz = new  long long [inline_b[2]];
  c_inline_nx = new  long long [inline_b[0]+1];
  c_inline_ny = new  long long [inline_b[1]+1];
  c_inline_nz = new  long long [inline_b[2]+1];

  c_block_dist[0] = new double[inline_b[0] + 1];
  c_block_dist[1] = new double[inline_b[1] + 1];
  c_block_dist[2] = new double[inline_b[2] + 1];
  
  block_dist[0] = new double[inline_b[0]];
  block_dist[1] = new double[inline_b[1]];
  block_dist[2] = new double[inline_b[2]];


  for(long long i = 0; i < inline_b[0]; i ++){
    a_inline_nx[i] = inline_n[0];
    c_inline_nx[i] = 0;
    nelx_tot += a_inline_nx[i];
    if(i)c_inline_nx[i] = c_inline_nx[i-1]+a_inline_nx[i-1];
    block_dist[0][i] =(inline_gmax[0] - inline_gmin[0])/(double)inline_b[0];
    c_block_dist[0][i] = inline_gmin[0];//inner radius
    if(i)c_block_dist[0][i] = c_block_dist[0][i-1]+block_dist[0][i-1];
  }
  c_inline_nx[inline_b[0]] = c_inline_nx[inline_b[0] - 1]+a_inline_nx[inline_b[0] - 1];
  c_block_dist[0][inline_b[0]] = c_block_dist[0][inline_b[0] - 1]+block_dist[0][inline_b[0] - 1];

  for(long long i = 0; i < inline_b[1]; i ++){
    a_inline_ny[i] = inline_n[1];
    c_inline_ny[i] = 0;
    nely_tot += a_inline_ny[i];
    if(i)c_inline_ny[i] = c_inline_ny[i-1]+a_inline_ny[i-1];
    block_dist[1][i] =(inline_gmax[1] - inline_gmin[1])/(double)inline_b[1];
    c_block_dist[1][i] = inline_gmin[1];
    if(i)c_block_dist[1][i] = c_block_dist[1][i-1]+block_dist[1][i-1];
  }
  c_inline_ny[inline_b[1]] = c_inline_ny[inline_b[1]-1]+a_inline_ny[inline_b[1]-1];
  c_block_dist[1][inline_b[1]] = c_block_dist[1][inline_b[1] - 1]+block_dist[1][inline_b[1] - 1];


  for(long long i = 0; i < inline_b[2]; i ++){
    a_inline_nz[i] = inline_n[2];
    c_inline_nz[i] = 0;
    nelz_tot += a_inline_nz[i];
    if(i)c_inline_nz[i] = c_inline_nz[i-1]+a_inline_nz[i-1];
    block_dist[2][i] =(inline_gmax[2] - inline_gmin[2])/(double)inline_b[2];
    c_block_dist[2][i] = inline_gmin[2];
    if(i)c_block_dist[2][i] = c_block_dist[2][i-1]+block_dist[2][i-1];
  }
  c_inline_nz[inline_b[2]] = c_inline_nz[inline_b[2]-1]+a_inline_nz[inline_b[2]-1];
  c_block_dist[2][inline_b[2]] = c_block_dist[2][inline_b[2] - 1]+block_dist[2][inline_b[2] - 1];

  cum_block_totals = new long long[inline_b[0]*inline_b[1]*inline_b[2]];
  els_in_block = new long long[inline_b[0]*inline_b[1]*inline_b[2]];

   long long bl_ct = 0;
  for(long long k = 0; k < inline_b[2]; k ++){
    for(long long j = 0; j < inline_b[1]; j ++){
      for(long long i = 0; i < inline_b[0]; i ++){
        els_in_block[bl_ct] = a_inline_nx[i]*a_inline_ny[j]*a_inline_nz[k];
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
   long long nnx = nelx_tot+1;
   long long nny = nely_tot+1;
  double xdelta = inline_gmax[0]-inline_gmin[0];
  double ydelta = inline_gmax[1]-inline_gmin[1];
  Icoors = new double[nnx];
  Jcoors = new double[nny];

   long long nct = 0;
  for(long long i = 0; i < inline_b[0]; i ++){
    for(long long j = 0; j < a_inline_nx[i]; j ++){
      Icoors[nct] = c_block_dist[0][i]+j*block_dist[0][i]/(double)a_inline_nx[i];
      Icoors[nct+1] = c_block_dist[0][i]+(j+1)*block_dist[0][i]/(double)a_inline_nx[i];
      nct ++;
    }
  }

  nct = 0;
  for(long long i = 0; i < inline_b[1]; i ++){
    for(long long j = 0; j < a_inline_ny[i]; j ++){
      Jcoors[nct] = c_block_dist[1][i]+j*block_dist[1][i]/(double)a_inline_ny[i];
      Jcoors[nct+1] = c_block_dist[1][i]+(j+1)*block_dist[1][i]/(double)a_inline_ny[i];
      nct ++;
    }
  }

  //   for(long long i = 0; i < nnx; i++)Icoors[i] = inline_gmin[0]+(double)i*dx;
  //   for(long long i = 0; i < nny; i++)Jcoors[i] = inline_gmin[1]+(double)i*dy;
  if(Element_Density_Functions[0])Element_Density_Functions[0]->Integrate(inline_gmin[0],inline_gmax[0], error_stream);
  if(!error_stream.str().empty()){return 1;}
  if(Element_Density_Functions[1])Element_Density_Functions[1]->Integrate(inline_gmin[1],inline_gmax[1], error_stream);
  if(!error_stream.str().empty()){return 1;}
  if(Element_Density_Functions[0]){
    for(long long ict = 0; ict < nnx; ict ++){
      double factor = (Icoors[ict]-inline_gmin[0])/xdelta;
      double interpolant =  Element_Density_Functions[0]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      double new_coord = inline_gmin[0]+interpolant*xdelta;
      Icoors[ict] = new_coord;
    }
  }
  if(Element_Density_Functions[1]){
    for(long long ict = 0; ict < nny; ict ++){
      double factor = (Jcoors[ict]-inline_gmin[1])/ydelta;
      double interpolant =  Element_Density_Functions[1]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      double new_coord = inline_gmin[1]+interpolant*ydelta;
      Jcoors[ict] = new_coord;
    }
  }
  if(dimension == 3){
   long long nnz = nelz_tot+1;
  double zdelta = inline_gmax[2]-inline_gmin[2];
  Kcoors = new double[nnz];

  nct = 0;
  for(int i = 0; i < inline_b[2]; i ++){
    for(long long j = 0; j < a_inline_nz[i]; j ++){
      Kcoors[nct] = c_block_dist[2][i]+j*block_dist[2][i]/(double)a_inline_nz[i];
      Kcoors[nct+1] = c_block_dist[2][i]+(j+1)*block_dist[2][i]/(double)a_inline_nz[i];
      nct ++;
    }
  }

  //   for(long long i = 0; i < nnz; i++)Kcoors[i] = inline_gmin[2]+(double)i*dz;
  if(Element_Density_Functions[2])Element_Density_Functions[2]->Integrate(inline_gmin[2],inline_gmax[2],error_stream);
  if(!error_stream.str().empty()){return 1;}
  if(Element_Density_Functions[2]){
    for(long long ict = 0; ict < nnz; ict ++){
      double factor = (Kcoors[ict]-inline_gmin[2])/zdelta;
      double interpolant =  Element_Density_Functions[2]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      double new_coord = inline_gmin[2]+interpolant*zdelta;
      Kcoors[ict] = new_coord;
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
      coords[the_local_node+0*num_nodes]= Icoors[global_i]*cos(Jcoors[global_j]*deg_to_rad);
      coords[the_local_node+1*num_nodes]= Icoors[global_i]*sin(Jcoors[global_j]*deg_to_rad);
      coords[the_local_node+2*num_nodes]= Kcoors[global_k];
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
    coords[the_local_node+0*num_nodes]= Icoors[global_i]*cos(Jcoors[global_j]*deg_to_rad);
    coords[the_local_node+1*num_nodes]= Icoors[global_i]*sin(Jcoors[global_j]*deg_to_rad);
    }
    else{
    coords[the_local_node+0*num_nodes]= Icoors[global_i]*cos(Jcoors[global_j]*deg_to_rad);
    coords[the_local_node+1*num_nodes]= Icoors[global_i]*sin(Jcoors[global_j]*deg_to_rad)*cos(Kcoors[global_k]*deg_to_rad);
    coords[the_local_node+2*num_nodes]= Icoors[global_i]*sin(Jcoors[global_j]*deg_to_rad)*sin(Kcoors[global_k]*deg_to_rad);
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
    coords[the_local_node+0*num_nodes] = Icoors[global_i];
    coords[the_local_node+1*num_nodes] = Jcoors[global_j];
  if(dimension == 3){
    coords[the_local_node+2*num_nodes] = Kcoors[global_k];
  }
  }
}
}//end namespace PAMGEN_NEVADA
