#include "inline_mesh_desc.h"
#include "radial_inline_mesh_desc.h"
#include "uns_inline_decomp.h"
#include "pamgen_fudges.h"
#include "Vector.h"

namespace PAMGEN_NEVADA {

/*****************************************************************************/
void Radial_Inline_Mesh_Desc::calculateSize(long long & total_el_count, 
						  long long & total_node_count, 
						  long long & total_edge_count)
/*****************************************************************************/
{
  total_el_count   = 0L; 
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

  long long j_add = 1L;
  if(periodic_j)j_add = 0L;
  total_node_count = (nodes_temp[0]+1L)*(nodes_temp[1]+j_add);

  if(dimension == 3){
    total_node_count *= nodes_temp[2] + 1L;
  }

  if(dimension == 3){
    total_edge_count =
      (nodes_temp[0]+0L)*(nodes_temp[1]+j_add)*(nodes_temp[2]+1L)+
      (nodes_temp[0]+1L)*(nodes_temp[1]+0L)*(nodes_temp[2]+1L)+
      (nodes_temp[0]+1L)*(nodes_temp[1]+j_add)*(nodes_temp[2]+0L);
  }  
  else{
    total_edge_count =
      (nodes_temp[0]+0L)*(nodes_temp[1]+j_add)+
      (nodes_temp[0]+1L)*(nodes_temp[1]+0L);
  }
}

/*****************************************************************************/
std::string Radial_Inline_Mesh_Desc::Calc_Intervals()
/*****************************************************************************/
{
  std::string errorString;
  for(long long axis = 0; axis < dimension;axis ++){
    for(long long i = 0; i < inline_b[axis]; i ++){
      if((first_size[axis][i] > 0.) && (last_size[axis][i] == 0.)){
	last_size[axis][i] = first_size[axis][i];
      }
      if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
	Real xave = (first_size[axis][i]+last_size[axis][i])/2.;
	Real k = a_lot_positive + (block_dist[axis][i]/xave);
	long long ktil =(long long)k;
	if(ktil < 1)ktil = 1;
	Real delxtil = block_dist[axis][i]/(Real)ktil;
	Real s = xave-delxtil;
	interval[axis][i] = ktil;
	first_size[axis][i]-=s;
	last_size[axis][i] -=s;
      }
    }
  }
  return errorString;
}


/*****************************************************************************/
long long Radial_Inline_Mesh_Desc::Set_Up()
/*****************************************************************************/
{
  a_inline_n[0] = new long long [inline_b[0]];
  a_inline_n[1] = new long long [inline_b[1]];
  a_inline_n[2] = new long long [inline_b[2]];
  c_inline_n[0] = new long long [inline_b[0]+1];
  c_inline_n[1] = new long long [inline_b[1]+1];
  c_inline_n[2] = new long long [inline_b[2]+1];
  
  for(long long i = 0; i < inline_b[0]; i ++){
    a_inline_n[0][i] = inline_n[0];
    a_inline_n[0][i] = interval[0][i];
    c_inline_n[0][i] = 0;
    nel_tot[0] += a_inline_n[0][i];
    if(i)c_inline_n[0][i] = c_inline_n[0][i-1]+a_inline_n[0][i-1];
    c_block_dist[0][i] = inline_gmin[0];//inner radius
    if(i)c_block_dist[0][i] = c_block_dist[0][i-1]+block_dist[0][i-1];
  }
  c_inline_n[0][inline_b[0]] = c_inline_n[0][inline_b[0] - 1]+a_inline_n[0][inline_b[0] - 1];
  c_block_dist[0][inline_b[0]] = c_block_dist[0][inline_b[0] - 1]+block_dist[0][inline_b[0] - 1];

  for(long long i = 0; i < inline_b[1]; i ++){
    a_inline_n[1][i] = inline_n[1];
    a_inline_n[1][i] = interval[1][i];
    c_inline_n[1][i] = 0;
    nel_tot[1] += a_inline_n[1][i];
    if(i)c_inline_n[1][i] = c_inline_n[1][i-1]+a_inline_n[1][i-1];
    c_block_dist[1][i] = inline_gmin[1];
    if(i)c_block_dist[1][i] = c_block_dist[1][i-1]+block_dist[1][i-1];
  }
  c_inline_n[1][inline_b[1]] = c_inline_n[1][inline_b[1]-1]+a_inline_n[1][inline_b[1]-1];
  c_block_dist[1][inline_b[1]] = c_block_dist[1][inline_b[1] - 1]+block_dist[1][inline_b[1] - 1];

  if(dimension == 3){
    for(long long i = 0; i < inline_b[2]; i ++){
      a_inline_n[2][i] = inline_n[2];
      a_inline_n[2][i] = interval[2][i];
      
      c_inline_n[2][i] = 0;
      nel_tot[2] += a_inline_n[2][i];
      if(i)c_inline_n[2][i] = c_inline_n[2][i-1]+a_inline_n[2][i-1];
      c_block_dist[2][i] = inline_gmin[2];
      if(i)c_block_dist[2][i] = c_block_dist[2][i-1]+block_dist[2][i-1];
    }
    c_inline_n[2][inline_b[2]] = c_inline_n[2][inline_b[2]-1]+a_inline_n[2][inline_b[2]-1];
    c_block_dist[2][inline_b[2]] = c_block_dist[2][inline_b[2] - 1]+block_dist[2][inline_b[2] - 1];
  }
  else{
    nel_tot[2] = 1;
    a_inline_n[2][0] = 1;
    c_inline_n[2][0] = 0;
    c_inline_n[2][1] = 1;
  }

  cum_block_totals = new long long[inline_b[0]*inline_b[1]*inline_b[2]];
  els_in_block = new long long[inline_b[0]*inline_b[1]*inline_b[2]];

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
  double total_angle = c_block_dist[1][inline_b[1]]-inline_gmin[1];
  if(fabs(total_angle - 360.0)< 1.0) periodic_j = true;
  
  
  inline_gmax[0] = inline_gmin[0];
  inline_gmax[1] = inline_gmin[1];    
  inline_gmax[2] = inline_gmin[2];
  for(long long i = 0; i < inline_b[0]; i ++){
    inline_gmax[0] += block_dist[0][i];
  }
  for(long long i = 0; i < inline_b[1]; i ++){
    inline_gmax[1] += block_dist[1][i];
  }
  if(dimension == 3){
    for(long long i = 0; i < inline_b[2]; i ++){
      inline_gmax[2] += block_dist[2][i];
    }
  }

  if(enforce_periodic){
    Real total_theta = c_block_dist[1][inline_b[1]];
   if(total_theta != 90. && total_theta != 180. && total_theta != 360.){
      error_stream << "Radial_Inline_Mesh_Desc::Set_Up(...): "
		   << "ENFORCE PERIODIC requires the extent of the mesh in theta to be 90, 180.0, or 360.0 degrees.";
      return 1;
    }
    //must have 90/180/360 degrees and even numbers of elements
    long long mod = (long long)(c_block_dist[1][inline_b[1]]/90.);
    if(nel_tot[1] % (mod*2) != 0){
      error_stream << "Radial_Inline_Mesh_Desc::Set_Up(...): "
		   << "ENFORCE PERIODIC Requires an even number of elements in ech 90 degree quadrant.";
      return 1;
    }
    if(inline_b[1] != 1){
      error_stream << "Radial_Inline_Mesh_Desc::Set_Up(...): "
		   << "ENFORCE PERIODIC Requires a single block in the circumferential direction.";
      return 1;
    }
  }
  return 0;
}

/****************************************************************************/
long long Radial_Inline_Mesh_Desc::Calc_Coord_Vectors()
/****************************************************************************/
{

  for(long long axis = 0; axis < dimension; axis ++){
    IJKcoors[axis] = new double[nel_tot[axis]+1];

    long long nct = 0;
    for(long long i = 0; i < inline_b[axis]; i ++){
      double sum = 0.;
      for(long long j = 0; j < a_inline_n[axis][i]; j ++){
	if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
	  IJKcoors[axis][nct] = c_block_dist[axis][i]+sum;
	  sum += first_size[axis][i];
	  if(interval[axis][i]-1) sum += (double)j*(last_size[axis][i]-first_size[axis][i])/((double)interval[axis][i]-1);
	  IJKcoors[axis][nct+1] = c_block_dist[axis][i+1];
	}
	else{
	  IJKcoors[axis][nct] = c_block_dist[axis][i]+j*block_dist[axis][i]/(double)a_inline_n[axis][i];
	  IJKcoors[axis][nct+1] = c_block_dist[axis][i]+(j+1)*block_dist[axis][i]/(double)a_inline_n[axis][i];
	}
	nct ++;
      }
    }
    if(Element_Density_Functions[axis]){
      Element_Density_Functions[axis]->Integrate(inline_gmin[axis],inline_gmax[axis], error_stream);
      if(!error_stream.str().empty()){return 1;}
      double delta = inline_gmax[axis]-inline_gmin[axis];
      for(long long ict = 0; ict < (nel_tot[axis]+1); ict ++){
	double factor = (IJKcoors[axis][ict]-inline_gmin[axis])/delta;
	double interpolant =  Element_Density_Functions[axis]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
	double new_coord = inline_gmin[axis]+interpolant*delta;
	IJKcoors[axis][ict] = new_coord;
      }
    }
  }

  return 0;
}

/****************************************************************************/
void Radial_Inline_Mesh_Desc::Populate_Coords(Real * coords,   
						    std::vector<long long> & global_node_vector,                             
						    std::map <long long, long long> & global_node_map,
						    long long num_nodes)

/****************************************************************************/
{
  Real deg_to_rad = M_PI/180.0;
  Real total_theta = c_block_dist[1][inline_b[1]];
  long long global_ind[3];
  for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
    long long the_node = global_node_vector[gnv];
    global_ind[2] = the_node/knstride;
    global_ind[1] = (the_node-global_ind[2]*knstride)/jnstride;
    global_ind[0] = the_node - global_ind[2]*knstride-global_ind[1]*jnstride;

    long long the_local_node = get_map_entry(global_node_map,the_node);
    coords[the_local_node+0*num_nodes]= IJKcoors[0][global_ind[0]]*cos(IJKcoors[1][global_ind[1]]*deg_to_rad);
    coords[the_local_node+1*num_nodes]= IJKcoors[0][global_ind[0]]*sin(IJKcoors[1][global_ind[1]]*deg_to_rad);
    if(dimension == 3){
      coords[the_local_node+2*num_nodes]= IJKcoors[2][global_ind[2]];
    }
    
    if(enforce_periodic){
      Vector tv = calc_coords_periodic(total_theta, global_ind[0], global_ind[1], global_ind[2]);
      coords[the_local_node+0*num_nodes]= tv.X();
      coords[the_local_node+1*num_nodes]= tv.Y();
      if(dimension == 3){
	coords[the_local_node+2*num_nodes]= tv.Z();
      }
    }
  }
}

/****************************************************************************/
Vector Radial_Inline_Mesh_Desc::calc_coords_periodic(double total_theta,
						     long long i, 
						     long long j, 
						     long long k)
/****************************************************************************/
{
  // this function is used if ENFORCE PERIODIC is requested
  // it calculates all coordinates in the first 45 degrees of the domain
  // and then transforms them to the appropriate octant
  long long per=0;
  if(total_theta == 90.)per = nel_tot[1]/2;
  if(total_theta == 180.)per = nel_tot[1]/4;
  if(total_theta == 360.)per = nel_tot[1]/8;
  
  long long jmod = j%per;
  long long jmult = j/per;
  if(jmult %2 == 0){//even
    
  }
  else{
    jmod = per - jmod;
  }
  double xval,yval,zval;
  double deg_to_rad = M_PI/180.0;

  xval = IJKcoors[0][i]*cos(IJKcoors[1][jmod]*deg_to_rad);
  yval = IJKcoors[0][i]*sin(IJKcoors[1][jmod]*deg_to_rad);
  if(jmod == per){
    xval = IJKcoors[0][i]*cos(IJKcoors[1][jmod]*deg_to_rad);
    yval = xval;
  }

  zval = 0.0;
  if(dimension == 3){
    zval = IJKcoors[2][k];
  }

  Vector res = Vector( xval,yval,zval);

  //transforming back to original quadrant coordinates
  switch(jmult){
  case  0:{
    break;
  }
  case 1:{
    res = Vector(res.Y(),res.X(),res.Z());
    break;
  }
  case 2:{
    res =  Vector(-res.Y(),res.X(),res.Z());
    break;
  }
  case 3:{
    res = Vector(-res.X(),res.Y(),res.Z());
    break;
  }
  case 4:{
    res = Vector(-res.X(),-res.Y(),res.Z());
    break;
  }
  case 5:{
    res = Vector(-res.Y(),-res.X(),res.Z());
    break;
  }
  case 6:{
    res = Vector(res.Y(),-res.X(),res.Z());
    break;
  }
  case 7:{
    res = Vector(res.X(),-res.Y(),res.Z());
    break;
  }
  }
  return res;
  
}

}// end namespace PAMGEN_NEVADA
