#include "inline_mesh_desc.h"
#include "radial_inline_mesh_desc.h"
#include "uns_inline_decomp.h"
#include "../asrc/fudges.h"
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
    for(int k = 0; k < inline_bz; k ++){
      for(int j = 0; j < inline_by; j ++){
	for(int i = 0; i < inline_bx; i ++){
	  total_el_count += 
	    (long long)interval[0][i]*
	    (long long)interval[1][j]*
	    (long long)interval[2][k];
	}
      }
    }
  }
  else{
    for(int j = 0; j < inline_by; j ++){
      for(int i = 0; i < inline_bx; i ++){
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

  for(int i = 0; i < inline_bx; i ++){
    nodes_temp[0] += interval[0][i];
  }
  for(int j = 0; j < inline_by; j ++){
    nodes_temp[1] += interval[1][j];
  }
  if(dimension == 3){
    for(int k = 0; k < inline_bz; k ++){
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
void Radial_Inline_Mesh_Desc::Calc_Intervals()
/*****************************************************************************/
{
  for(int i = 0; i < inline_bx; i ++){
    int axis = 0;
    if((first_size[axis][i] > 0.) && (last_size[axis][i] == 0.)){
      last_size[axis][i] = first_size[axis][i];
    }
    if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
      Real xave = (first_size[axis][i]+last_size[axis][i])/2.;
      Real k = a_lot_positive + (block_dist[axis][i]/xave);
      int ktil =(int)k;
      if(ktil < 1)ktil = 1;
      Real delxtil = block_dist[axis][i]/(Real)ktil;
      Real s = xave-delxtil;
      interval[axis][i] = ktil;
      first_size[axis][i]-=s;
      last_size[axis][i] -=s;
    }
  }
  for(int i = 0; i < inline_by; i ++){
    int axis = 1;
    if((first_size[axis][i] > 0.) && (last_size[axis][i] == 0.)){
      last_size[axis][i] = first_size[axis][i];
    }
    if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
      Real xave = (first_size[axis][i]+last_size[axis][i])/2.;
      Real k = a_lot_positive + (block_dist[axis][i]/xave);
      int ktil =(int)k;
      if(ktil < 1)ktil = 1;
      Real delxtil = block_dist[axis][i]/(Real)ktil;
      Real s = xave-delxtil;
      interval[axis][i] = ktil;
      first_size[axis][i]-=s;
      last_size[axis][i] -=s;
    }
  }

  if(dimension == 3){
    for(int i = 0; i < inline_bz; i ++){
      int axis = 2;
      if((first_size[axis][i] > 0.) && (last_size[axis][i] == 0.)){
	last_size[axis][i] = first_size[axis][i];
      }
      if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
	Real xave = (first_size[axis][i]+last_size[axis][i])/2.;
	Real k = a_lot_positive + (block_dist[axis][i]/xave);
	int ktil =(int)k;
	if(ktil < 1)ktil = 1;
	Real delxtil = block_dist[axis][i]/(Real)ktil;
	Real s = xave-delxtil;
	interval[axis][i] = ktil;
	first_size[axis][i]-=s;
	last_size[axis][i] -=s;
      }
    }
  }
}


/*****************************************************************************/
int Radial_Inline_Mesh_Desc::Set_Up()
/*****************************************************************************/
{
  a_inline_nx = new int [inline_bx];
  a_inline_ny = new int [inline_by];
  a_inline_nz = new int [inline_bz];
  c_inline_nx = new int [inline_bx+1];
  c_inline_ny = new int [inline_by+1];
  c_inline_nz = new int [inline_bz+1];
  
  for(int i = 0; i < inline_bx; i ++){
    a_inline_nx[i] = inline_nx;
    a_inline_nx[i] = interval[0][i];
    c_inline_nx[i] = 0;
    nelx_tot += a_inline_nx[i];
    if(i)c_inline_nx[i] = c_inline_nx[i-1]+a_inline_nx[i-1];
    c_block_dist[0][i] = inline_gminx;//inner radius
    if(i)c_block_dist[0][i] = c_block_dist[0][i-1]+block_dist[0][i-1];
  }
  c_inline_nx[inline_bx] = c_inline_nx[inline_bx - 1]+a_inline_nx[inline_bx - 1];
  c_block_dist[0][inline_bx] = c_block_dist[0][inline_bx - 1]+block_dist[0][inline_bx - 1];

  for(int i = 0; i < inline_by; i ++){
    a_inline_ny[i] = inline_ny;
    a_inline_ny[i] = interval[1][i];
    c_inline_ny[i] = 0;
    nely_tot += a_inline_ny[i];
    if(i)c_inline_ny[i] = c_inline_ny[i-1]+a_inline_ny[i-1];
    c_block_dist[1][i] = inline_gminy;
    if(i)c_block_dist[1][i] = c_block_dist[1][i-1]+block_dist[1][i-1];
  }
  c_inline_ny[inline_by] = c_inline_ny[inline_by-1]+a_inline_ny[inline_by-1];
  c_block_dist[1][inline_by] = c_block_dist[1][inline_by - 1]+block_dist[1][inline_by - 1];

  if(dimension == 3){
    for(int i = 0; i < inline_bz; i ++){
      a_inline_nz[i] = inline_nz;
      a_inline_nz[i] = interval[2][i];
      
      c_inline_nz[i] = 0;
      nelz_tot += a_inline_nz[i];
      if(i)c_inline_nz[i] = c_inline_nz[i-1]+a_inline_nz[i-1];
      c_block_dist[2][i] = inline_gminz;
      if(i)c_block_dist[2][i] = c_block_dist[2][i-1]+block_dist[2][i-1];
    }
    c_inline_nz[inline_bz] = c_inline_nz[inline_bz-1]+a_inline_nz[inline_bz-1];
    c_block_dist[2][inline_bz] = c_block_dist[2][inline_bz - 1]+block_dist[2][inline_bz - 1];
  }
  else{
    nelz_tot = 1;
    a_inline_nz[0] = 1;
    c_inline_nz[0] = 0;
    c_inline_nz[1] = 1;
  }

  cum_block_totals = new int[inline_bx*inline_by*inline_bz];
  els_in_block = new int[inline_bx*inline_by*inline_bz];

  int bl_ct = 0;
  for(int k = 0; k < inline_bz; k ++){
    for(int j = 0; j < inline_by; j ++){
      for(int i = 0; i < inline_bx; i ++){
        els_in_block[bl_ct] = a_inline_nx[i]*a_inline_ny[j]*a_inline_nz[k];
        cum_block_totals[bl_ct]=0;
        if(bl_ct){
	  cum_block_totals[bl_ct] = cum_block_totals[bl_ct-1]+els_in_block[bl_ct-1];
        }
        bl_ct ++;
      }
    }
  }

  // this tolerance is dangerous
  double total_angle = c_block_dist[1][inline_by]-inline_gminy;
  if(fabs(total_angle - 360.0)< 1.0) periodic_j = true;
  
  
  inline_gmaxx = inline_gminx;
  inline_gmaxy = inline_gminy;    
  inline_gmaxz = inline_gminz;
  for(int i = 0; i < inline_bx; i ++){
    inline_gmaxx += block_dist[0][i];
  }
  for(int i = 0; i < inline_by; i ++){
    inline_gmaxy += block_dist[1][i];
  }
  if(dimension == 3){
    for(int i = 0; i < inline_bz; i ++){
      inline_gmaxz += block_dist[2][i];
    }
  }

  if(enforce_periodic){
    Real total_theta = c_block_dist[1][inline_by];
   if(total_theta != 90. && total_theta != 180. && total_theta != 360.){
      error_stream << "Radial_Inline_Mesh_Desc::Set_Up(...): "
		   << "ENFORCE PERIODIC requires the extent of the mesh in theta to be 90, 180.0, or 360.0 degrees.";
      return 1;
    }
    //must have 90/180/360 degrees and even numbers of elements
    int mod = (int)(c_block_dist[1][inline_by]/90.);
    if(nely_tot % (mod*2) != 0){
      error_stream << "Radial_Inline_Mesh_Desc::Set_Up(...): "
		   << "ENFORCE PERIODIC Requires an even number of elements in ech 90 degree quadrant.";
      return 1;
    }
    if(inline_by != 1){
      error_stream << "Radial_Inline_Mesh_Desc::Set_Up(...): "
		   << "ENFORCE PERIODIC Requires a single block in the circumferential direction.";
      return 1;
    }
  }
  return 0;
}

/****************************************************************************/
int Radial_Inline_Mesh_Desc::Calc_Coord_Vectors()
/****************************************************************************/
{
  int nnx = nelx_tot+1;
  int nny = nely_tot+1;
  Real xdelta = inline_gmaxx-inline_gminx;
  Real ydelta = inline_gmaxy-inline_gminy;
  Icoors = new Real[nnx];
  Jcoors = new Real[nny];

  int nct = 0;
  for(int i = 0; i < inline_bx; i ++){
    int axis = 0;
    Real sum = 0.;
    for(int j = 0; j < a_inline_nx[i]; j ++){
      if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
        Icoors[nct] = c_block_dist[axis][i]+sum;
        sum += first_size[axis][i];
        if(interval[axis][i]-1) sum += (Real)j*(last_size[axis][i]-first_size[axis][i])/((Real)interval[axis][i]-1);
        Icoors[nct+1] = c_block_dist[axis][i+1];
      }
      else{
        Icoors[nct] = c_block_dist[0][i]+j*block_dist[0][i]/(Real)a_inline_nx[i];
        Icoors[nct+1] = c_block_dist[0][i]+(j+1)*block_dist[0][i]/(Real)a_inline_nx[i];
      }
      nct ++;
    }
  }

  nct = 0;
  for(int i = 0; i < inline_by; i ++){
    int axis = 1;
    Real sum = 0.;
    for(int j = 0; j < a_inline_ny[i]; j ++){
      if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
        Jcoors[nct] = c_block_dist[axis][i]+sum;
        sum += first_size[axis][i];
        if(interval[axis][i]-1) sum += (Real)j*(last_size[axis][i]-first_size[axis][i])/((Real)interval[axis][i]-1);
        Jcoors[nct+1] = c_block_dist[axis][i+1];
      }
      else{
        Jcoors[nct] = c_block_dist[1][i]+j*block_dist[1][i]/(Real)a_inline_ny[i];
        Jcoors[nct+1] = c_block_dist[1][i]+(j+1)*block_dist[1][i]/(Real)a_inline_ny[i];
      }
      nct ++;
    }
  }
  if(Element_Density_Functions[0])Element_Density_Functions[0]->Integrate(inline_gminx,inline_gmaxx, error_stream);
  if(!error_stream.str().empty()){return 1;}
  if(Element_Density_Functions[1])Element_Density_Functions[1]->Integrate(inline_gminy,inline_gmaxy, error_stream);
  if(!error_stream.str().empty()){return 1;}
  if(Element_Density_Functions[0]){
    for(int ict = 0; ict < nnx; ict ++){
      Real factor = (Icoors[ict]-inline_gminx)/xdelta;
      Real interpolant =  Element_Density_Functions[0]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      Real new_coord = inline_gminx+interpolant*xdelta;
      Icoors[ict] = new_coord;
    }
  }
  if(Element_Density_Functions[1]){
    for(int ict = 0; ict < nny; ict ++){
      Real factor = (Jcoors[ict]-inline_gminy)/ydelta;
      Real interpolant =  Element_Density_Functions[1]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      Real new_coord = inline_gminy+interpolant*ydelta;
      Jcoors[ict] = new_coord;
    }
  }

  if(dimension == 3){
    int nnz = nelz_tot+1;
    Real zdelta = inline_gmaxz-inline_gminz;
    Kcoors = new Real[nnz];
    
    nct = 0;
    for(int i = 0; i < inline_bz; i ++){
      int axis = 2;
      Real sum = 0.;
      for(int j = 0; j < a_inline_nz[i]; j ++){
	if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
	  Kcoors[nct] = c_block_dist[axis][i]+sum;
	  sum += first_size[axis][i];
          if(interval[axis][i]-1) sum += (Real)j*(last_size[axis][i]-first_size[axis][i])/((Real)interval[axis][i]-1);
	  Kcoors[nct+1] = c_block_dist[axis][i+1];
	}
	else{
	  Kcoors[nct] = c_block_dist[2][i]+j*block_dist[2][i]/(Real)a_inline_nz[i];
	  Kcoors[nct+1] = c_block_dist[2][i]+(j+1)*block_dist[2][i]/(Real)a_inline_nz[i];
	}
	nct ++;
      }
    }
    
    if(Element_Density_Functions[2])Element_Density_Functions[2]->Integrate(inline_gminz,inline_gmaxz, error_stream);
    if(!error_stream.str().empty()){return 1;}
    if(Element_Density_Functions[2]){
      for(int ict = 0; ict < nnz; ict ++){
	Real factor = (Kcoors[ict]-inline_gminz)/zdelta;
	Real interpolant =  Element_Density_Functions[2]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
	Real new_coord = inline_gminz+interpolant*zdelta;
	Kcoors[ict] = new_coord;
      }
    }
  }
  return 0;
}

/****************************************************************************/
void Radial_Inline_Mesh_Desc::Populate_Coords(Real * coords,   
						    std::vector<int> & global_node_vector,                             
						    std::map <int, int> & global_node_map,
						    int num_nodes)

/****************************************************************************/
{
  Real deg_to_rad = M_PI/180.0;
  Real total_theta = c_block_dist[1][inline_by];
  for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
    int the_node = global_node_vector[gnv];
    int global_k = the_node/knstride;
    int global_j = (the_node-global_k*knstride)/jnstride;
    int global_i = the_node - global_k*knstride-global_j*jnstride;
    int the_local_node = get_map_entry(global_node_map,the_node);
    coords[the_local_node+0*num_nodes]= Icoors[global_i]*cos(Jcoors[global_j]*deg_to_rad);
    coords[the_local_node+1*num_nodes]= Icoors[global_i]*sin(Jcoors[global_j]*deg_to_rad);
    if(dimension == 3){
      coords[the_local_node+2*num_nodes]= Kcoors[global_k];
    }
    
    if(enforce_periodic){
      Vector tv = calc_coords_periodic(total_theta, global_i, global_j, global_k);
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
						     int i, 
						     int j, 
						     int k)
/****************************************************************************/
{
  // this function is used if ENFORCE PERIODIC is requested
  // it calculates all coordinates in the first 45 degrees of the domain
  // and then transforms them to the appropriate octant
  int per=0;
  if(total_theta == 90.)per = nely_tot/2;
  if(total_theta == 180.)per = nely_tot/4;
  if(total_theta == 360.)per = nely_tot/8;
  
  int jmod = j%per;
  int jmult = j/per;
  if(jmult %2 == 0){//even
    
  }
  else{
    jmod = per - jmod;
  }
  double xval,yval,zval;
  double deg_to_rad = M_PI/180.0;

  xval = Icoors[i]*cos(Jcoors[jmod]*deg_to_rad);
  yval = Icoors[i]*sin(Jcoors[jmod]*deg_to_rad);
  if(jmod == per){
    xval = Icoors[i]*cos(Jcoors[jmod]*deg_to_rad);
    yval = xval;
  }

  zval = 0.0;
  if(dimension == 3){
    zval = Kcoors[k];
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
