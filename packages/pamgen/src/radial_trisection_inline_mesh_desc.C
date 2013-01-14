#include "uns_inline_decomp.h"
#include "inline_mesh_desc.h"
#include "radial_inline_mesh_desc.h"
#include "radial_trisection_inline_mesh_desc.h"
#include <algorithm>
#include "calc_decomp_cuts.h"
#include "pamgen_fudges.h"
#include "Random.h"
#include <stdlib.h>


namespace PAMGEN_NEVADA {

bool part_compare_size4(const Partition *a, const Partition *b) {
  if(a->numels < b->numels)
    return true;
  if(a->numels == b->numels)
    return (a->unique_id < b->unique_id);
  return false;
}

bool part_compare_size2(const Partition *a, const Partition *b) {
  if((a->highs[3] - a->lows[3]) < (b->highs[3] - b->lows[3]))return true;
  if(a->numels < b->numels)
    return true;
  if(a->numels == b->numels)
    return (a->unique_id < b->unique_id);
  return false;
}

bool part_compare_size3(const Partition *a, const Partition *b) {
  if((a->highs[3] - a->lows[3]) > (b->highs[3] - b->lows[3]))return true;
  if(a->numels > b->numels)
    return true;
  if(a->numels == b->numels)
    return (a->unique_id > b->unique_id);
  return false;
}


/*****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::Calc_Intervals()
/*****************************************************************************/
{
  if(transition_radius <= 0.){
    transition_radius = block_dist[0][0]/(sqrt(2.));
  }
 
 for(long long i = 0; i < inline_bx; i ++){
    long long axis = 0;
    if((first_size[axis][i] > 0.) && (last_size[axis][i] == 0.)){
      last_size[axis][i] = first_size[axis][i];
    }
    if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
      Real xave = (first_size[axis][i]+last_size[axis][i])/2.;
      Real xdist = block_dist[axis][i];
      if(i == 0)xdist -= transition_radius;
      Real k = a_lot_positive + (xdist/xave);
      long long ktil =(long long)k;
      if(ktil < 1)ktil = 1;
      Real delxtil = xdist/(Real)ktil;
      Real s = xave-delxtil;
      interval[axis][i] = ktil;
      first_size[axis][i]-=s;
      last_size[axis][i] -=s;
    }
  }
  for(long long i = 0; i < inline_by; i ++){
    long long axis = 1;
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
  if(dimension == 3){
    for(long long i = 0; i < inline_bz; i ++){
      long long axis = 2;
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
}

/*****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::calculateSize(long long & total_el_count, 
							     long long & total_node_count, 
							     long long & total_edge_count)
/*****************************************************************************/
{
  total_el_count = 0L;
  total_node_count = 0L;
  total_edge_count = 0L;

  long long nodes_temp[3];
  nodes_temp[0] = 0;
  nodes_temp[1] = 0;
  nodes_temp[2] = 0;

  for(long long i = 0; i < inline_bx; i ++){
    nodes_temp[0] += interval[0][i];
  }
  for(long long j = 0; j < inline_by; j ++){
    nodes_temp[1] += interval[1][j];
  }

  if(dimension == 3){

    for(long long k = 0; k < inline_bz; k ++){
      nodes_temp[2] += interval[2][k];
    }
  }
  else{
    nodes_temp[2] = 1L;
  }
  
  long long total_theta = nodes_temp[1];
  long long local_div = total_theta/(2L*(long long)trisection_blocks);

  total_el_count =0L; 
  
  total_el_count = nodes_temp[0]*nodes_temp[1]*nodes_temp[2];
  total_el_count += local_div*local_div*trisection_blocks*nodes_temp[2];

  long long j_add = 1L;
  if(periodic_j)j_add = 0L;

  total_node_count = (nodes_temp[0]+1L)*(nodes_temp[1]+j_add);

  if(periodic_j){
    total_node_count += (long long)(trisection_blocks)*((local_div-1L)*(local_div +0L));
    total_node_count += 1L;
  }
  else{
    total_node_count += (local_div + 0L)*(local_div + 0L);
    if(trisection_blocks > 1){
      total_node_count += (long long)(trisection_blocks-1L)*((local_div-1L)*(local_div +0L));
    }
  }

  if(dimension == 3){
    total_node_count *= nodes_temp[2] + 1L;
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
  
  
  if(periodic_j){
    if(dimension == 3){
      total_edge_count += (long long)(trisection_blocks) * 
	((local_div*(local_div - 0L)*(nodes_temp[2]+1L))+ 
	 (local_div*(local_div - 1L)*(nodes_temp[2]+1L))+ 
	 (local_div*(local_div - 1L)*(nodes_temp[2]+0L))); 
      total_edge_count += nodes_temp[2];
    }    
    else{
      total_edge_count += (long long)(trisection_blocks) * 
	((local_div*(local_div - 0L))+ 
	 (local_div*(local_div - 1L))); 
    }
  }
  else{
    // the first trisection block
    if(dimension == 3){
      total_edge_count += 2L*local_div*local_div*(nodes_temp[2]+1L);
      total_edge_count += local_div*local_div*(nodes_temp[2]);
    }
    else{
      total_edge_count += 2L*local_div*local_div;
    }
    
    if(trisection_blocks > 1){
      if(dimension == 3){
	total_edge_count += (long long)(trisection_blocks-1L) * 
	  ((local_div*(local_div - 0L)*(nodes_temp[2]+1L))+ 
	   (local_div*(local_div - 1L)*(nodes_temp[2]+1L))+ 
	   (local_div*(local_div - 1L)*(nodes_temp[2]+0L))); 
      }
      else{
	total_edge_count += (long long)(trisection_blocks-1L) * 
	  ((local_div*(local_div - 0L))+ 
	   (local_div*(local_div - 1L))); 
      }
    }
  }
}

/*****************************************************************************/
long long Radial_Trisection_Inline_Mesh_Desc::Set_Up()
/*****************************************************************************/
{

  inline_gminx = 0.;//always zero for trisection


  a_inline_nx = new long long [inline_bx];
  a_inline_ny = new long long [inline_by];
  a_inline_nz = new long long [inline_bz];
  c_inline_nx = new long long [inline_bx+1];
  c_inline_ny = new long long [inline_by+1];
  c_inline_nz = new long long [inline_bz+1];

  for(long long i = 0; i < inline_bx; i ++){
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

  for(long long i = 0; i < inline_by; i ++){
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
    for(long long i = 0; i < inline_bz; i ++){
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

  // this tolerance is dangerous
  double total_angle = c_block_dist[1][inline_by] - inline_gminy;;
  if(fabs(total_angle-360.0)< 1.0) periodic_j = true;
  
  
  inline_gmaxx = inline_gminx;
  inline_gmaxy = inline_gminy;    
  inline_gmaxz = inline_gminz;
  for(long long i = 0; i < inline_bx; i ++){
    inline_gmaxx += block_dist[0][i];
  }
  for(long long i = 0; i < inline_by; i ++){
    inline_gmaxy += block_dist[1][i];
  }
  if(dimension == 3){
    for(long long i = 0; i < inline_bz; i ++){
      inline_gmaxz += block_dist[2][i];
    }
  }
  
  long long total_theta = c_inline_ny[inline_by];
  div = total_theta/(2*trisection_blocks);
  mod = total_theta % (2*trisection_blocks);

  if(mod != 0){
    error_stream << "Radial_Trisection_Inline_Mesh_Desc::Set_Up"
		 << "The number of elements in the circumferential direction " << total_theta << "\n"
		 << "must be divisible by twice the number of trisection_blocks " << trisection_blocks << "." << "\n";
    return 1;
  }


  cum_block_totals = new long long[inline_bx*inline_by*inline_bz];
  els_in_block = new long long[inline_bx*inline_by*inline_bz];


  long long bl_ct = 0;
  for(long long k = 0; k < inline_bz; k ++){
    for(long long j = 0; j < inline_by; j ++){
      for(long long i = 0; i < inline_bx; i ++){
	bl_ct = i + j * inline_bx + k * inline_bx*inline_by;
        els_in_block[bl_ct] = a_inline_nx[i]*a_inline_ny[j]*a_inline_nz[k];
      }
    }
  }

  bl_ct = 0;
  long long total_element_count = 0;
  for(long long k = 0; k < inline_bz; k ++){

    for(long long j = 0; j < inline_by; j ++){
      for(long long i = 0; i < 1; i ++){
	bl_ct = i + j * inline_bx + k * inline_bx*inline_by;
 
	cum_block_totals[bl_ct]=total_element_count;
	total_element_count += els_in_block[bl_ct];

      }
    }

    for(long long j = 0; j < inline_by; j ++){
      for(long long i = 1; i < inline_bx; i ++){
	bl_ct = i + j * inline_bx + k * inline_bx*inline_by;

	cum_block_totals[bl_ct]=total_element_count;
	total_element_count += els_in_block[bl_ct];
      }
    }
    total_element_count += trisection_blocks*a_inline_nz[k]*div*div;

  }


  if(enforce_periodic){
    Real ctotal_theta = c_block_dist[1][inline_by];
    if(ctotal_theta != 90. && ctotal_theta != 180. && ctotal_theta != 360.){
      error_stream << "Radial_Inline_Mesh_Desc::Set_Up(...): "
	  << "ENFORCE PERIODIC requires the extent of the mesh in theta to be 90, 180.0, or 360.0 degrees.";
      return 1;
    }
    //must have 90/180/360 degrees and even numbers of elements
    long long num_90_sections = (long long)(c_block_dist[1][inline_by]/90.);
    if(trisection_blocks != num_90_sections){
      error_stream << "Radial_Inline_Mesh_Desc::Set_Up(...): "
	  << "ENFORCE PERIODIC Requires a single trisection block in each quadrant.";
      return 1;
    }
  }
  return 0;
}

/****************************************************************************/
long long Radial_Trisection_Inline_Mesh_Desc::Calc_Coord_Vectors()
/****************************************************************************/
{
  long long nnx = nelx_tot+1;
  long long nny = nely_tot+1;
  Real xdelta = inline_gmaxx-inline_gminx;
  Real ydelta = inline_gmaxy-inline_gminy;
  Icoors = new Real[nnx];
  Jcoors = new Real[nny];

  long long nct = 0;
  for(long long i = 0; i < inline_bx; i ++){
    long long axis = 0;
    Real sum = 0.;
    for(long long j = 0; j < a_inline_nx[i]; j ++){
      if(i == 0){
	if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
	  Icoors[nct] = c_block_dist[axis][i]+transition_radius + sum;
	  sum += first_size[axis][i];
          if(interval[axis][i]-1) sum += (Real)j*(last_size[axis][i]-first_size[axis][i])/((Real)interval[axis][i]-1);
	  Icoors[nct+1] = c_block_dist[axis][i+1];
	}
	else{
	  Icoors[nct] = c_block_dist[0][i] + transition_radius + j*(block_dist[0][i]-transition_radius)/(Real)a_inline_nx[i];
	  Icoors[nct+1] = c_block_dist[0][i]+ transition_radius + (j+1)*(block_dist[0][i]-transition_radius)/(Real)a_inline_nx[i];
	}
      }
      else{
	if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
	  Icoors[nct] = c_block_dist[axis][i] + sum;
	  sum += first_size[axis][i];
          if(interval[axis][i]-1) sum += (Real)j*(last_size[axis][i]-first_size[axis][i])/((Real)interval[axis][i]-1);
	  Icoors[nct+1] = c_block_dist[axis][i+1];
	}
	else{
	  Icoors[nct] = c_block_dist[0][i]+j*block_dist[0][i]/(Real)a_inline_nx[i];
	  Icoors[nct+1] = c_block_dist[0][i]+(j+1)*block_dist[0][i]/(Real)a_inline_nx[i];
	}
      }
      nct ++;
    }
  }

  nct = 0;
  for(long long i = 0; i < inline_by; i ++){
    long long axis = 1;
    Real sum = 0.;
    for(long long j = 0; j < a_inline_ny[i]; j ++){
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
  std::stringstream terror_stream;
  if(Element_Density_Functions[0])Element_Density_Functions[0]->Integrate(inline_gminx,inline_gmaxx,terror_stream);
  if(!terror_stream.str().empty()){error_stream << terror_stream.str();return 1;}
  if(Element_Density_Functions[1])Element_Density_Functions[1]->Integrate(inline_gminy,inline_gmaxy,terror_stream);
  if(!terror_stream.str().empty()){error_stream << terror_stream.str();return 1;}
  if(Element_Density_Functions[0]){
    for(long long ict = 0; ict < nnx; ict ++){
      Real factor = (Icoors[ict]-inline_gminx)/xdelta;
      Real interpolant =  Element_Density_Functions[0]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      Real new_coord = inline_gminx+interpolant*xdelta;
      Icoors[ict] = new_coord;
    }
  }
  if(Element_Density_Functions[1]){
    for(long long ict = 0; ict < nny; ict ++){
      Real factor = (Jcoors[ict]-inline_gminy)/ydelta;
      Real interpolant =  Element_Density_Functions[1]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
      Real new_coord = inline_gminy+interpolant*ydelta;
      Jcoors[ict] = new_coord;
    }
  }
  if(dimension == 3){
    long long nnz = nelz_tot+1;
    Real zdelta = inline_gmaxz-inline_gminz;
    Kcoors = new Real[nnz];
    
    nct = 0;
    for(long long i = 0; i < inline_bz; i ++){
      long long axis = 2;
      Real sum = 0.;
      for(long long j = 0; j < a_inline_nz[i]; j ++){
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
    
    if(Element_Density_Functions[2])Element_Density_Functions[2]->Integrate(inline_gminz,inline_gmaxz,terror_stream);
    if(!terror_stream.str().empty()){error_stream << terror_stream.str();return 1;}
    if(Element_Density_Functions[2]){
      for(long long ict = 0; ict < nnz; ict ++){
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
Vector Radial_Trisection_Inline_Mesh_Desc::calc_coords_periodic_trisect_blocks(double total_theta,
									       long long nl,
									       long long ni, 
									       long long nj, 
									       long long nk,
									       Quad_Patch ** quads)
/****************************************************************************/
{
  //This function is called if ENFORCE PERIODICITY is requested
  // it calculates the coordinats based on the lower i,j section of the first 
  // trisection block (the first octant, first 45 degrees) 
  // and then transforms the result back to the original octant
  long long per = 0;
  if(total_theta == 90.)per = nely_tot/2;
  if(total_theta == 180.)per = nely_tot/4;
  if(total_theta == 360.)per = nely_tot/8;
  
  long long jmod = nj%per;
  long long jmult = nj/per;
  if(jmult %2 == 0){//even
    
  }
  else{
    jmod = per - jmod;
  }
  double xval,yval,zval;
  double deg_to_rad = M_PI/180.0;

  Real angle = total_theta *((Real)jmod/(Real)(c_inline_ny[inline_by]));
  Real outer_radius = c_block_dist[0][1];
  
  // need inner radius.
  // this is a vile trick. By getting the node number of the node on the border with the trisection
  // blocks you can then get the l,i,j,k for that node and use the quads to get the inner radius for 
  // this node. The outer radius is just the outer radius for block 0.
  
  long long border_nn = get_node_number_from_l_i_j_k(nl,0,jmod,nk);
  
  long long bl,bi,bj,bk;
  get_l_i_j_k_from_node_number(border_nn,bl,bi,bj,bk);
  long long maxi = div+1;
  long long maxj = div+1;
  Vector loc = quads[bl]->Interpolate(Vector((Real)bi/(Real)maxi,(Real)bj/(Real)maxj,0.));
  
  Real param = (Icoors[ni]-transition_radius)/(block_dist[0][0]-transition_radius);
  
  xval = loc.x  + param*( outer_radius*cos(deg_to_rad*angle) -loc.x);
  yval = loc.y  + param*( outer_radius*sin(deg_to_rad*angle) -loc.y);
  zval = 0.0;
  if(dimension == 3){
    zval = Kcoors[nk];
  }
  
  if(jmod == per){
    yval = xval;
  }


  Vector res = Vector( xval,yval,zval);

  //transforming back to original octant
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

/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::Populate_Coords(Real * coords,   
							 std::vector<long long> & global_node_vector,                             
							 std::map <long long, long long> & global_node_map,
							 long long num_nodes)

/****************************************************************************/
{
  long long num_nodes_per_element = 4;

  if(dimension == 3){
    num_nodes_per_element = 8;
  }
 

  Real deg_to_rad = M_PI/180.0;

  Real total_theta = c_block_dist[1][inline_by]-inline_gminy;
  Real dtheta = total_theta/(Real)(2*trisection_blocks);

  Quad_Patch ** quads;
  quads = new Quad_Patch *[trisection_blocks];

  for(long long i = 0; i < trisection_blocks; i ++){
   quads[i] = new Quad_Patch(Vector(0,0,0),
			      Vector(transition_radius*cos(deg_to_rad*(inline_gminy + dtheta*(Real)(2*i+0))), 
				     transition_radius*sin(deg_to_rad*(inline_gminy + dtheta*(Real)(2*i+0))), 0.),
			      Vector(transition_radius*cos(deg_to_rad*(inline_gminy + dtheta*(Real)(2*i+1))), 
				     transition_radius*sin(deg_to_rad*(inline_gminy + dtheta*(Real)(2*i+1))), 0.),
			      Vector(transition_radius*cos(deg_to_rad*(inline_gminy + dtheta*(Real)(2*i+2))), 
				     transition_radius*sin(deg_to_rad*(inline_gminy + dtheta*(Real)(2*i+2))), 0.));
  }
  
  for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
    long long l,i,j,k;
    long long the_node = global_node_vector[gnv];
    get_l_i_j_k_from_node_number(the_node,l,i,j,k);

    long long the_local_node = get_map_entry(global_node_map,the_node);
    if(l == trisection_blocks){
      coords[the_local_node+0*num_nodes]= Icoors[i]*cos(Jcoors[j]*deg_to_rad);
      coords[the_local_node+1*num_nodes]= Icoors[i]*sin(Jcoors[j]*deg_to_rad);
      if(dimension == 3){
	coords[the_local_node+2*num_nodes]= Kcoors[k];
      }
      //enforce_periodicity on outer blocks
      if(enforce_periodic){
	Vector tv = calc_coords_periodic(total_theta, i, j, k);
	coords[the_local_node+0*num_nodes]= tv.X();
	coords[the_local_node+1*num_nodes]= tv.Y();
	if(dimension == 3){
	  coords[the_local_node+2*num_nodes]= tv.Z();
	}
      }      
    }
    else{
      long long maxi = div+1;
      long long maxj = div+1;
      if(dimension == 3){
	coords[the_local_node+2*num_nodes]= Kcoors[k];
      }
      Vector loc = quads[l]->Interpolate(Vector((Real)i/(Real)maxi,(Real)j/(Real)maxj,0.));
      coords[the_local_node+0*num_nodes] = loc.x;
      coords[the_local_node+1*num_nodes] = loc.y;
      //Evaluate each coord based on the first octant (45.0 degrees) 
      // and then transform them back
      if(enforce_periodic){
	if(l == 0){
	  if(j > i){
	    loc = quads[l]->Interpolate(Vector((Real)j/(Real)maxi,(Real)i/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = loc.y;
	    coords[the_local_node+1*num_nodes] = loc.x;
	  }
	  else if (j == i){
	    loc = quads[l]->Interpolate(Vector((Real)j/(Real)maxi,(Real)i/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = loc.y;
	    coords[the_local_node+1*num_nodes] = loc.y;
	  }
	}
	if(l == 1){
	  long long tj = i;
	  long long ti = j;
	  if(tj > ti){
	    loc = quads[0]->Interpolate(Vector((Real)tj/(Real)maxi,(Real)ti/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = -loc.y;
	    coords[the_local_node+1*num_nodes] = loc.x;
	  }
	  else if (tj == ti){
	    loc = quads[0]->Interpolate(Vector((Real)tj/(Real)maxi,(Real)ti/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = -loc.y;
	    coords[the_local_node+1*num_nodes] = loc.y;
	  }
	  else{
	    loc = quads[0]->Interpolate(Vector((Real)ti/(Real)maxi,(Real)tj/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = -loc.x;
	    coords[the_local_node+1*num_nodes] = loc.y;
	  }
	}

	if(l == 2){
	  long long tj = j;
	  long long ti = i;
	  if(tj > ti){
	    loc = quads[0]->Interpolate(Vector((Real)tj/(Real)maxi,(Real)ti/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = -loc.y;
	    coords[the_local_node+1*num_nodes] = -loc.x;
	  }
	  else if (tj == ti){
	    loc = quads[0]->Interpolate(Vector((Real)tj/(Real)maxi,(Real)ti/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = -loc.y;
	    coords[the_local_node+1*num_nodes] = -loc.y;
	  }
	  else{
	    loc = quads[0]->Interpolate(Vector((Real)ti/(Real)maxi,(Real)tj/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = -loc.x;
	    coords[the_local_node+1*num_nodes] = -loc.y;
	  }
	}
	if(l == 3){
	  long long tj = i;
	  long long ti = j;
	  if(tj > ti){
	    loc = quads[0]->Interpolate(Vector((Real)tj/(Real)maxi,(Real)ti/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = loc.y;
	    coords[the_local_node+1*num_nodes] = -loc.x;
	  }
	  else if (tj == ti){
	    loc = quads[0]->Interpolate(Vector((Real)tj/(Real)maxi,(Real)ti/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = loc.y;
	    coords[the_local_node+1*num_nodes] = -loc.y;
	  }
	  else{
	    loc = quads[0]->Interpolate(Vector((Real)ti/(Real)maxi,(Real)tj/(Real)maxj,0.));
	    coords[the_local_node+0*num_nodes] = loc.x;
	    coords[the_local_node+1*num_nodes] = -loc.y;
	  }
	}
	
      }
    }
    
  }


  //try looping over elements block wise with block 0 holding the trisection blocks
  for(unsigned bct = 0; bct < numBlocks();bct += blockKstride()){
    for(unsigned elct = 0;elct < element_block_lists[bct].size();elct++){
      long long the_el = element_block_lists[bct][elct];
      long long tl;
      long long ti;
      long long tj;
      long long tk;
      long long nn[8];
      long long ln[8];

      get_l_i_j_k_from_element_number(the_el,tl,ti,tj,tk);
      nn[0] = get_node_number_from_l_i_j_k(tl, ti + 0, tj + 0, tk + 0);
      nn[1] = get_node_number_from_l_i_j_k(tl, ti + 1, tj + 0, tk + 0);
      nn[2] = get_node_number_from_l_i_j_k(tl, ti + 1, tj + 1, tk + 0);
      nn[3] = get_node_number_from_l_i_j_k(tl, ti + 0, tj + 1, tk + 0);

      ln[0] = get_map_entry(global_node_map,nn[0]);
      ln[1] = get_map_entry(global_node_map,nn[1]);
      ln[2] = get_map_entry(global_node_map,nn[2]);
      ln[3] = get_map_entry(global_node_map,nn[3]);

      if(dimension == 3){
	nn[4] = get_node_number_from_l_i_j_k(tl, ti + 0, tj + 0, tk + 1);
	nn[5] = get_node_number_from_l_i_j_k(tl, ti + 1, tj + 0, tk + 1);
	nn[6] = get_node_number_from_l_i_j_k(tl, ti + 1, tj + 1, tk + 1);
	nn[7] = get_node_number_from_l_i_j_k(tl, ti + 0, tj + 1, tk + 1);

	ln[4] = get_map_entry(global_node_map,nn[4]);
	ln[5] = get_map_entry(global_node_map,nn[5]);
	ln[6] = get_map_entry(global_node_map,nn[6]);
	ln[7] = get_map_entry(global_node_map,nn[7]);
      }

      for(long long i = 0; i < num_nodes_per_element; i ++){
	long long nl;
	long long ni;
	long long nj;
	long long nk;
	get_l_i_j_k_from_node_number(nn[i],nl,ni,nj,nk);
	if((nl == trisection_blocks) && (ni > 0)){// these are block 0 nodes that need to be repositioned

	  Real angle = total_theta *((Real)nj/(Real)(c_inline_ny[inline_by]));
	  angle += inline_gminy;
	  Real outer_radius = c_block_dist[0][1];

	  // need inner radius.
	  // this is a vile trick. By getting the node number of the node on the border with the trisection
	  // blocks you can then get the l,i,j,k for that node and use the quads to get the inner radius for 
	  // this node. The outer radius is just the outer radius for block 0.
	  
	  long long border_nn = get_node_number_from_l_i_j_k(nl,0,nj,nk);

	  long long bl,bi,bj,bk;
	  get_l_i_j_k_from_node_number(border_nn,bl,bi,bj,bk);
	  long long maxi = div+1;
	  long long maxj = div+1;
	  Vector loc = quads[bl]->Interpolate(Vector((Real)bi/(Real)maxi,(Real)bj/(Real)maxj,0.));

	  Real param = (Icoors[ni]-transition_radius)/(block_dist[0][0]-transition_radius);

	  coords[ln[i]+0*num_nodes] = loc.x  + param*( outer_radius*cos(deg_to_rad*angle) -loc.x);
	  coords[ln[i]+1*num_nodes] = loc.y  + param*( outer_radius*sin(deg_to_rad*angle) -loc.y);        

	  if(enforce_periodic){
	    Vector res = calc_coords_periodic_trisect_blocks( total_theta,
							      nl,
							      ni, 
							      nj, 
							      nk,
							      quads);
	    coords[ln[i]+0*num_nodes] = res.X();
	    coords[ln[i]+1*num_nodes] = res.Y();
	    if(dimension == 3){  
	      coords[ln[i]+2*num_nodes] = res.Z();
	    }	    
	    
	  }//enforce periodic
 
	}
      }
    }
  }
  // clean up

  for(long long i = 0; i < trisection_blocks; i ++){
    delete quads[i];
  }
  delete [] quads;
}


/*****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::setStrides()
/*****************************************************************************/
{
  instride = 1;
  jnstride = nelx_tot+1;
  
  if(inline_geometry_type == RADIAL && periodic_j){
    instride = 1;
    jnstride = nelx_tot+1;
  }
  
  iestride = 1;
  jestride = nelx_tot;
  kestride = (nelx_tot)*(nely_tot);
  for(long long i = 0; i < trisection_blocks; i ++){
    kestride += div * div;
  }

  knstride = 0;
  // sum of trisection squares
  for(long long i = 0; i < trisection_blocks; i ++){
    knstride += (div+1)*(div+1);
    if(i)knstride -= div+1;
  }
  if(periodic_j){
    knstride -=div;
  }
  nn_center = knstride;
  // add rainbow
  knstride += (nely_tot+1)*(nelx_tot+1-1);
  if(periodic_j){
    knstride -= (nelx_tot+1-1);
  }

  //cumulative nodes in trisection blocks
  tri_block_cum_nn = new long long [trisection_blocks+1];
  long long cnode_count = 0;
  for(long long i = 0; i < trisection_blocks; i ++){
    tri_block_cum_nn[i] = cnode_count;
    cnode_count += (div+1)*(div+1);
    if(i)cnode_count -= div+1;
  }
  tri_block_cum_nn[trisection_blocks] = nn_center;
}


//! Partitions all elements by a recursive bisection into 
//! rectangular chunks. Alternative decompositions are possible.
//! This one is simple to implement and easy to run on all processors.
//! Zoltan could be used for this. The elements on the processor are
//! packed into the stl list.
/****************************************************************************/
Partition * Radial_Trisection_Inline_Mesh_Desc::Decompose(std::list <long long> & global_el_ids,
							  long long & error_code)
  /****************************************************************************/
{
  std::vector <Partition *> sorted_partition_list;
  
  if(inline_decomposition_type == SEQUENTIAL){
    long long total = kestride * nelz_tot;
    long long num_per_proc = total/num_processors;
    long long remainder = total - num_per_proc*num_processors;
    long long my_start = my_rank * num_per_proc;
    long long my_end = my_start + num_per_proc;
    if(my_rank == num_processors-1)my_end +=remainder;

    for(long long mtotal = my_start; mtotal < my_end; mtotal ++){
      global_el_ids.push_back(mtotal);
    }
    return base_partition;
  }
  else if(inline_decomposition_type == RANDOM){
    long long total = kestride * nelz_tot;
    for(long long i = 0; i < total; i ++){
      SRANDOM(i);
      long long rand_num = RANDOM();
      long long proc = rand_num%num_processors;
      if(proc == my_rank)global_el_ids.push_back(i);
    }
    return base_partition;
  }
  else if(inline_decomposition_type == BISECTION){

    long long remaining_cuts[3];
    long long decomp_result = 0;
    if(dimension == 3){
      //fix i/r/x at 1 to prevent r cuts
      decomp_result = dom_decomp_3d(1,nely_tot,nelz_tot,num_processors,
					&(inc_nels[0]),&(inc_nels[1]),&(inc_nels[2]));
    }
    else{
      decomp_result = dom_decomp_2d(1,nely_tot,num_processors,
					&(inc_nels[0]),&(inc_nels[1]));
    }
    
    if(decomp_result != 0){
      error_stream << "Terminating from Inline_Mesh_Desc::Decompose, ";
      error_stream << "non-zero return value from dom_decomp_2/3d ";
      error_stream << decomp_result;
      error_code = 1;
      return NULL;
    }

    if(dimension == 3){
      remaining_cuts[0] = 1;
      remaining_cuts[1] = inc_nels[1];
      remaining_cuts[2] = inc_nels[2];
      inc_nels[0] = nelx_tot/inc_nels[0];
      inc_nels[1] = nely_tot/inc_nels[1];
      inc_nels[2] = nelz_tot/inc_nels[2];
    }
    else{
      remaining_cuts[0] = 1;
      remaining_cuts[1] = inc_nels[1];
      remaining_cuts[2] = 1;
      inc_nels[0] = nelx_tot/inc_nels[0];
      inc_nels[1] = nely_tot/inc_nels[1];
      inc_nels[2] = 1;
    }


      info_stream << "Using BISECTION LAYOUT decomposition.\n";
      info_stream << "Number of elements/segment in directions I/X/R \t\t" << inc_nels[0] << "\n";
      info_stream << "Number of elements/segment in directions J/Y/THETA \t" << inc_nels[1] << "\n";
      info_stream << "Number of elements/segment in directions K/Z/PHI \t" << inc_nels[2] << "\n";
      info_stream << "Number of mesh segments in directions I/X/R \t\t" << remaining_cuts[0] << "\n";
      info_stream << "Number of mesh segments in directions J/Y/THETA \t" << remaining_cuts[1] << "\n";
      info_stream << "Number of mesh segments in directions K/Z/PHI \t" << remaining_cuts[2] << "\n";

      info_stream << "Using  BISECTION decomposition." << std::endl
		  << " No cuts will be made in radial direction." << std::endl;
	
    base_partition  = new Partition(0,0,0,
				    0,1,
				    nelx_tot,nely_tot,nelz_tot,
				    inline_decomposition_type,remaining_cuts);

    sorted_partition_list.push_back(base_partition);
    Partition* biggest;
    
    while(sorted_partition_list.size() < num_processors){
      //  get first entry in list and the partition object it points to
      biggest = sorted_partition_list.back();
      //  remove the entry from the sorted_list
      sorted_partition_list.pop_back();
      
      //  bisect the partition object and
      //  add the children of the partition object to the sorted list
      biggest->Processor_Partition(sorted_partition_list, inc_nels);
      //  re-sort the list
      std::sort(sorted_partition_list.begin(),sorted_partition_list.end(),part_compare_size4);
    }

    std::sort(sorted_partition_list.begin(),sorted_partition_list.end(),part_compare_size2);


    //Assign proc id numbers to the centroid sorted list entries
    std::vector < Partition * > :: iterator citer;
    long long proc_cnt = 0;
    Partition *my_part = NULL;
    for(citer = sorted_partition_list.begin();citer != sorted_partition_list.end();citer ++,proc_cnt++){
      (*citer)->proc_id = proc_cnt;
      if(proc_cnt == my_rank){
	my_part = (*citer);
      }
    }

    // we have split up the mesh using only the non-trisection portion of the mesh.
    // the assignment of elements must be expanded to include the trisection elements
    //Then loop over the elements for the current processor and push them into
    //the global el_ids list
    
    for(long long k = my_part->lows[2]; k < my_part->highs[2]; k++){
      for(long long j = my_part->lows[1]; j < my_part->highs[1]; j++){
	for(long long i = my_part->lows[0]; i < my_part->highs[0]; i++){
	  long long elnum = get_element_number_from_l_i_j_k(trisection_blocks,i,j,k);	  
	  global_el_ids.push_back(elnum);
	}
	//add in the trisection blocks here if controlled.
	if(((j/div) %2) == 0){
	  long long l = j/(2*div);
	  long long jj = j % (2*div);
	  for(long long ii = 0; ii < div; ii ++){
	    long long elnum =  get_element_number_from_l_i_j_k(l,ii,jj,k);
	    global_el_ids.push_back(elnum);
	  }
	}
      }
    }
    

  }
  else if(inline_decomposition_type == PROCESSOR_LAYOUT){
    if(num_processors != (long long)(inline_nprocs[0]*inline_nprocs[1]*inline_nprocs[2])){
      error_stream << "Inline_Mesh_Desc::Decompose"
	  << "The specified inline processor layout " << "\n"
	  << inline_nprocs[0] << " X/R " 
	  << inline_nprocs[1] << " Y/THETA " 
	  << inline_nprocs[2] << " Z " << "\n"
	  << " does not correspond to the number of processors " << num_processors << "." ;
      error_code = 1;
      return NULL;
    }
    if(inline_nprocs[0] != 1){
      std::ostringstream oss;
      error_stream << "The number of processors in the X/R direction for a trisection mesh is limited to 1\n"
	  << "The value " << "\n"
	  << inline_nprocs[0] << " X/R " 
	  << " is unacceptable.";
      error_code = 1;
      return NULL;
    }
    long long remaining_cuts[3];
    remaining_cuts[0] = 1;
    remaining_cuts[1] = inline_nprocs[1];
    if(dimension == 3){
      remaining_cuts[2] = inline_nprocs[2];
    }
    else{
      remaining_cuts[2] = 1;
    }
    
    inc_nels[0] = nelx_tot/inline_nprocs[0];
    inc_nels[1] = nely_tot/inline_nprocs[1];
    if(dimension == 3){
      inc_nels[2] = nelz_tot/inline_nprocs[2];
    }
    else{
      inc_nels[2] = 1;
    }
    
    if(inc_nels[0] == 0 ){
      error_stream << "Radial_Trisection_inline_Mesh_Desc::Decompose"
	  << " Value for numprocs specified in I direction " << inline_nprocs[0] 
	  << " is greater than the number of elements in I " << nelx_tot << ".";
      error_code = 1;
      return NULL;
    }

    if(inc_nels[1] == 0 ){
      error_stream << "Radial_Trisection_inline_Mesh_Desc::Decompose"
	  << " Value for numprocs specified in J direction " << inline_nprocs[1] 
	  << " is greater than the number of elements in J " << nely_tot << ".";
      error_code = 1;
      return NULL;
    }

    if(inc_nels[2] == 0 ){
      error_stream << "Radial_Trisection_inline_Mesh_Desc::Decompose"
	  << " Value for numprocs specified in K direction " << inline_nprocs[2] 
		   << " is greater than the number of elements in K " << nelz_tot << ".";
      error_code = 1;
      return NULL;
    }
	
    info_stream << "Using PROCESSOR LAYOUT decomposition.\n";
    info_stream << "Number of elements/segment in directions I/X/R \t\t" << inc_nels[0] << "\n";
    info_stream << "Number of elements/segment in directions J/Y/THETA \t" << inc_nels[1] << "\n";
    info_stream << "Number of elements/segment in directions K/Z/PHI \t" << inc_nels[2] << "\n";
    info_stream << "Number of mesh segments in directions I/X/R \t\t" << remaining_cuts[0] << "\n";
    info_stream << "Number of mesh segments in directions J/Y/THETA \t" << remaining_cuts[1] << "\n";
    info_stream << "Number of mesh segments in directions K/Z/PHI \t" << remaining_cuts[2] << "\n";

    base_partition  = new Partition(0,0,0,
				    0,1,
				    nelx_tot,nely_tot,nelz_tot,
				    inline_decomposition_type,remaining_cuts);

    sorted_partition_list.push_back(base_partition);
    Partition* biggest;
    
    while(sorted_partition_list.size() < num_processors){
      //  get first entry in list and the partition object it points to
      biggest = sorted_partition_list.back();
      //  remove the entry from the sorted_list
      sorted_partition_list.pop_back();
      
      //  bisect the partition object and
      //  add the children of the partition object to the sorted list
      biggest->Processor_Partition(sorted_partition_list, inc_nels);
      //  re-sort the list
      std::sort(sorted_partition_list.begin(),sorted_partition_list.end(),part_compare_size4);
    }

    std::sort(sorted_partition_list.begin(),sorted_partition_list.end(),part_compare_size2);


    //Assign proc id numbers to the centroid sorted list entries
    std::vector < Partition * > :: iterator citer;
    long long proc_cnt = 0;
    Partition *my_part = NULL;
    for(citer = sorted_partition_list.begin();citer != sorted_partition_list.end();citer ++,proc_cnt++){
      (*citer)->proc_id = proc_cnt;
      if(proc_cnt == my_rank){
	my_part = (*citer);
      }
    }

    // we have split up the mesh using only the non-trisection portion of the mesh.
    // the assignment of elements must be expanded to include the trisection elements
    //Then loop over the elements for the current processor and push them into
    //the global el_ids list
    
    for(long long k = my_part->lows[2]; k < my_part->highs[2]; k++){
      for(long long j = my_part->lows[1]; j < my_part->highs[1]; j++){
	for(long long i = my_part->lows[0]; i < my_part->highs[0]; i++){
	  long long elnum = get_element_number_from_l_i_j_k(trisection_blocks,i,j,k);	  
	  global_el_ids.push_back(elnum);
	}
	//add in the trisection blocks here if controlled.
	if(((j/div) %2) == 0){
	  long long l = j/(2*div);
	  long long jj = j % (2*div);
	  for(long long ii = 0; ii < div; ii ++){
	    long long elnum =  get_element_number_from_l_i_j_k(l,ii,jj,k);
	    global_el_ids.push_back(elnum);
	  }
	}
      }
    }
    
  }
  global_el_ids.sort();
  return base_partition;
}

/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::get_l_and_remainder_from_elno(long long el,
									     long long & l, 
									     long long & remainder)
/****************************************************************************/
{
  remainder = el;
  long long c_block_size = 0;
  for(long long i = 0; i < trisection_blocks; i ++){
     c_block_size += div * div;
    if(c_block_size > el){
      l = i;
      return;
    }
    remainder -=  div * div;
  }
  l = trisection_blocks;
  return;
}

/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::get_l_and_remainder_from_node_number(long long nn,
									     long long & l, 
									     long long & remainder)
/****************************************************************************/
{
  l = 0;
  remainder = nn;
  for(long long i = 0; i < trisection_blocks; i ++){
    if(tri_block_cum_nn[i+1] > nn){
      l = i;
      remainder = nn-tri_block_cum_nn[i];
      return;
    }
  }

  
  return;
}

/****************************************************************************/
long long Radial_Trisection_Inline_Mesh_Desc::get_node_number_from_l_i_j_k( long long  l,
								     long long  i,
								     long long  j,
								     long long  k)
/****************************************************************************/
{
  long long nno = k*knstride;
  if(l == trisection_blocks){//rainbow
    if(i == 0){
      long long soff = (j-1)/(div);// this is actually an element calculation
      long long loff = soff/2;
      if(soff%2==0){//even
	long long iloc = div;
	long long jloc = j-soff*div;
	nno +=get_node_number_from_l_i_j_k(loff,iloc,jloc,0);
      }
      else{//odd
	long long iloc = div-(j-soff*div);
	long long jloc = div;
	nno +=get_node_number_from_l_i_j_k(loff,iloc,jloc,0);
      }
    }
    else{
      if((j == nely_tot) && periodic_j){
	nno+= nn_center;
	nno += i-1;
      }
      else{
	nno+= nn_center;
	nno += i-1;
	nno += j*nelx_tot;
      }
    }
  }
  else{
    if (l == 0){
      nno += i + j * (div+1);
    }
    else if((l == (trisection_blocks-1)) && periodic_j && (i == 0)){//360 case
      long long loff = 0;
      long long ioff = j;
      long long joff = 0;
      nno += get_node_number_from_l_i_j_k(loff,ioff,joff,0);
    }
    else if((l == (trisection_blocks-1)) && periodic_j && (j != 0)){//360 case
	nno += tri_block_cum_nn[l];
	nno += (i-1) + (j-1) * (div);
    }
    else{
      if (j == 0){
	nno += get_node_number_from_l_i_j_k(l-1,0,i,0);
      }
      else{
	nno += tri_block_cum_nn[l];
	nno += i + (j-1) * (div+1);
      }
    }
  }
  return nno;
}

/****************************************************************************/
long long Radial_Trisection_Inline_Mesh_Desc::get_element_number_from_l_i_j_k( long long  l,
								     long long  i,
								     long long  j,
								     long long  k)
/****************************************************************************/
{
  //adding bounds checking to deal with neighbor calcultion
  // return -1 if off mesh
  // adjust and recurse if i or j throw element onto adjacent block
  if(k < 0) return -1;
  if(k >= nelz_tot) return -1;
  if(l == trisection_blocks){
    
    if(periodic_j){
      //MOD for for 360
      if(j < 0)j = nely_tot-1;
      if(j >= nely_tot)j = 0;
    }
    else{
      if(j < 0) return -1;
      if(j >= nely_tot)return -1;
    }

    if(i < 0){
      long long soff = (j)/(div);
      long long loff = soff/2;
      if(soff%2==0){//even
	long long iloc = div-1;
	long long jloc = j-soff*div;
	return get_element_number_from_l_i_j_k(loff,iloc,jloc,k);
      }
      else{//odd
	long long iloc = (div-1)-(j-soff*div);
	long long jloc = div-1;
	return get_element_number_from_l_i_j_k(loff,iloc,jloc,k);
      }
    }
    if(i >= nelx_tot)return -1;
  }
  else{
    //i
    if(i < 0){
      if(l+1 == trisection_blocks){
	if(periodic_j){
	  //mod for 360
	  long long loff = 0;
	  long long ioff = j;
	  long long joff = 0;
	  return get_element_number_from_l_i_j_k(loff,ioff,joff,k);
	}
	else{
	  return -1;
	}
      }
      long long loff = l+1;
      long long joff = 0;
      long long ioff = j;
      return get_element_number_from_l_i_j_k(loff,ioff,joff,k);
    }
    if(i >= div){
      long long loff = trisection_blocks;
      long long ioff = 0;
      long long joff = j + l*2*div;
      return get_element_number_from_l_i_j_k(loff,ioff,joff,k);
    }
    if(j < 0){
      if(l == 0){
	if(periodic_j){
	  //mod for 360
	  long long loff = trisection_blocks-1;
	  long long ioff = 0;
	  long long joff = i;
	  return get_element_number_from_l_i_j_k(loff,ioff,joff,k);
	}
	else{
	  return -1;
	}
      }
      long long loff = l - 1;
      long long ioff = 0;
      long long joff = i;
      return get_element_number_from_l_i_j_k(loff,ioff,joff,k);
    }
    if(j >= div){
      long long loff = trisection_blocks;
      long long ioff = 0;
      long long joff = (l+1)*2*div - 1 - i;
      return get_element_number_from_l_i_j_k(loff,ioff,joff,k);
    }
    //j
  }

  long long elno;
  elno = k*kestride;
  for(long long ict = 0; ict < l; ict ++){
    elno += div * div;
  }
  if(l == trisection_blocks){
    elno += i + j * jestride;
  }
  else{
    elno += i + j * div;
  }

  return elno;
}

//! A utility function to build up required bookkeeping objects.
/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::get_l_i_j_k_from_element_number(long long el,
								     long long & l,
								     long long & i,
								     long long & j,
								     long long & k)
/****************************************************************************/
{
  k = el/kestride;
  long long lel = el-k*kestride;
  long long remainder;
  get_l_and_remainder_from_elno(lel,l,remainder);
  if(l == trisection_blocks){
    j = (remainder)/(jestride);
    i = remainder - j*(jestride);
  }
  else{
    j = remainder/(div);
    i = remainder - j * (div);
  }
}

//! A utility function to build up required bookkeeping objects.
/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::get_l_i_j_k_from_node_number(long long nn,
									    long long & l,
									    long long & i,
									    long long & j,
									    long long & k)
/****************************************************************************/
{
  k = nn/knstride;
  long long remainder = nn - k*knstride;  
  if (remainder/nn_center){
    if(periodic_j){
      l = trisection_blocks;
      remainder -= nn_center;
      j = remainder / nelx_tot;
      i = remainder - j*nelx_tot + 1;
    }
    else{
      l = trisection_blocks;
      remainder -= nn_center;
      j = remainder / nelx_tot;
      i = remainder - j*nelx_tot + 1;
    }
  }
  else{
    get_l_and_remainder_from_node_number(remainder, l, remainder);
    if(periodic_j && (l == (trisection_blocks - 1))){
      j = remainder /(div);
      i = remainder - j * (div);
      j++;
      i++;
    }
    else{
      j = remainder /(div+1);
      i = remainder - j * (div+1);
      if(l)j++;
    }
  }
}

/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::Build_Global_Lists(std::list <long long> & element_list,
                                       std::vector <long long> & element_vector,
                                       std::list <long long> & global_node_list,
                                       std::vector <long long> & global_node_vector,
                                       std::map <long long, long long> & global_node_map,
                                       std::map <long long, long long> & global_element_map)
/****************************************************************************/
{
  //can muck with element_list and element vector here

  element_block_lists = new std::vector <long long> [numBlocks()];

  std::vector < std::pair < long long, long long > > * tblist;
  tblist = new std::vector < std::pair < long long , long long > > [numBlocks()];

  std::list <long long> ::iterator lit;
  for(lit = element_list.begin();lit != element_list.end();lit ++){
    long long the_element = *lit;
    // These are the indices of the element in the entire domain
    long long global_l;
    long long global_i;
    long long global_j;
    long long global_k;

    get_l_i_j_k_from_element_number(the_element,global_l,global_i,global_j,global_k);

    //allow only blocks in i/x/r direction
    long long block_k  = get_block_index(global_k,inline_bz,c_inline_nz);
    long long block_j = 0;
    if(global_l == trisection_blocks){
      block_j = get_block_index(global_j,inline_by,c_inline_ny);
    }
    long long block_i = 0;
    if(global_l == trisection_blocks){
      block_i = get_block_index(global_i,inline_bx,c_inline_nx);
    }   

    // This is the ordinal number of the block the element resides in
    if(block_i == 0)block_j = 0;
    long long local_block = block_i + block_j*(inline_bx-1)+ block_k*(blockKstride());

    {
      long long global_block_el_id = GetBlockBasedGlobalID(the_element, local_block);
      std::pair <long long ,long long > the_pair(global_block_el_id,the_element);
      tblist[local_block].push_back(the_pair);
    }
    long long nn;

      nn=get_node_number_from_l_i_j_k(global_l, global_i + 0, global_j + 0, global_k + 0); global_node_list.push_back(nn);
      nn=get_node_number_from_l_i_j_k(global_l, global_i + 1, global_j + 0, global_k + 0); global_node_list.push_back(nn);
      nn=get_node_number_from_l_i_j_k(global_l, global_i + 1, global_j + 1, global_k + 0); global_node_list.push_back(nn);
      nn=get_node_number_from_l_i_j_k(global_l, global_i + 0, global_j + 1, global_k + 0); global_node_list.push_back(nn);

      if(dimension == 3){
	nn=get_node_number_from_l_i_j_k(global_l, global_i + 0, global_j + 0, global_k + 1); global_node_list.push_back(nn);
	nn=get_node_number_from_l_i_j_k(global_l, global_i + 1, global_j + 0, global_k + 1); global_node_list.push_back(nn);
	nn=get_node_number_from_l_i_j_k(global_l, global_i + 1, global_j + 1, global_k + 1); global_node_list.push_back(nn);
	nn=get_node_number_from_l_i_j_k(global_l, global_i + 0, global_j + 1, global_k + 1); global_node_list.push_back(nn);
      }
  }

  for(long long ict = 0; ict < numBlocks(); ict++){
    std::sort(tblist[ict].begin(),tblist[ict].end());
  }

  for(long long ict = 0; ict < numBlocks(); ict++){
   std::vector < std::pair < long long, long long > > :: iterator tit;
   for(tit = tblist[ict].begin(); tit != tblist[ict].end(); tit ++){
     element_block_lists[ict].push_back((*tit).second);
   }
  }

  for(long long ict = 0; ict < numBlocks(); ict++){
    std::vector <long long> :: iterator vlit;
    for(vlit = element_block_lists[ict].begin();vlit != element_block_lists[ict].end();vlit++){
      element_vector.push_back(*vlit);
    }
  }
  
  global_node_list.sort();
  global_node_list.unique();


  // Create the global_node_map
  std::list <long long> ::iterator nit;
  long long total = 0;
  for(nit = global_node_list.begin();nit != global_node_list.end();nit ++,total++){
    global_node_vector.push_back(*nit);
    global_node_map[*nit] = total;
  }

  // Create the global_element_map
  long long total_element_count = 0;
  for(long long bct = 0; bct < numBlocks();bct ++ ){
    for(unsigned elct = 0;elct < element_block_lists[bct].size();elct++,total_element_count++){
      long long the_el = element_block_lists[bct][elct];
      global_element_map[the_el] = total_element_count;
    }
  }

  delete[]tblist;
}

/****************************************************************************/
LoopLimits Radial_Trisection_Inline_Mesh_Desc::get_tri_block_limits(bool is_block_bc, 
								    long long input_id, 
								    long long l_value, 
								    Topo_Loc the_location, 
								    Topo_Loc & new_location,
								    long long nodeset_plus_1)
  /****************************************************************************/
{
  LoopLimits ll;


  long long bid = input_id-1;
  long long kind = bid/(inline_bx*inline_by);
  long long jind = (bid-kind*(inline_bx * inline_by))/inline_bx;
  long long iind = bid - jind *(inline_bx) - kind*(inline_bx * inline_by);
     

  new_location = the_location;
  switch(the_location) {
  
  case ALL_NODES:{
    ll = getLimits(the_location,0,div+nodeset_plus_1,0,div+nodeset_plus_1,0,nelz_tot+nodeset_plus_1,div+nodeset_plus_1,div+nodeset_plus_1);
    break;
  }
  case Z_AXIS:{
    ll = getLimits(EDGE4,0,div+nodeset_plus_1,0,div+nodeset_plus_1,0,nelz_tot+nodeset_plus_1,div+nodeset_plus_1,div+nodeset_plus_1);
    break;
  }
  case PLUS_I:{
    //nop
    break;
  }
  case MINUS_J:{
    if(l_value==0){
      if(is_block_bc){
      ll = getLimits(the_location,
		     0,div+nodeset_plus_1,
		     0,div+nodeset_plus_1,
		     c_inline_nz[kind],c_inline_nz[kind+1]+nodeset_plus_1,
		     div+nodeset_plus_1,div+nodeset_plus_1);
      }
      else{
	ll = getLimits(the_location,
		       0,div+nodeset_plus_1,
		       0,div+nodeset_plus_1,
		       0,nelz_tot+nodeset_plus_1,
		       div+nodeset_plus_1,div+nodeset_plus_1);
      }
    }
    
    break;
  }
  case PLUS_J:{
    if(l_value==(trisection_blocks-1)){// need to change the Topo_Loc for this one
      if(is_block_bc){
      ll = getLimits(MINUS_I,
		     0,div+nodeset_plus_1,
		     0,div+nodeset_plus_1,
		     c_inline_nz[kind],c_inline_nz[kind+1]+nodeset_plus_1,
		     div+nodeset_plus_1,div+nodeset_plus_1);
      new_location = MINUS_I;
      }
      else{
	ll = getLimits(MINUS_I,
		       0,div+nodeset_plus_1,
		       0,div+nodeset_plus_1,
		       0,nelz_tot+nodeset_plus_1,
		       div+nodeset_plus_1,div+nodeset_plus_1);
	new_location = MINUS_I;
      }
    }
    break;
  }
  case MINUS_K:{
    if(is_block_bc){
    ll = getLimits(the_location,
		   0,div+nodeset_plus_1,
		   0,div+nodeset_plus_1,
		   c_inline_nz[kind],c_inline_nz[kind+1]+nodeset_plus_1,
		   div+nodeset_plus_1,div+nodeset_plus_1);
    }
    else{
      ll = getLimits(the_location,
		     0,div+nodeset_plus_1,
		     0,div+nodeset_plus_1,
		     0,nelz_tot+nodeset_plus_1,
		     div+nodeset_plus_1,div+nodeset_plus_1);
    }
    break;
  }
  case PLUS_K:{
    if(is_block_bc){
      ll = getLimits(the_location,
		     0,div+nodeset_plus_1,
		     0,div+nodeset_plus_1,
		     c_inline_nz[kind],c_inline_nz[kind+1]+nodeset_plus_1,
		     div+nodeset_plus_1,div+nodeset_plus_1);
    }
    else{
      ll = getLimits(the_location,
		     0,div+nodeset_plus_1,
		     0,div+nodeset_plus_1,
		     0,nelz_tot+nodeset_plus_1,
		     div+nodeset_plus_1,div+nodeset_plus_1);
    }
    break;		     
  }
  default:
    
    break;
  }   
  return ll;
}

//! Reads in/creates the serial component of the unstructured mesh in parallel
/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::Calc_Serial_Component(Partition * my_part,
                                          std::vector <long long> & element_vector,
                                          std::list <long long> & global_node_list,
                                          std::vector<long long> & global_node_vector,
                                          std::map <long long, long long> & global_node_map,
                                          std::map <long long, long long> & global_element_map)
/****************************************************************************/
{
  //NODESETS
  // Nodesets are simple because we are numbering nodes across the entire domain
  // sequentially in i,j,k
  // The loop limits are the index bounds that allow traversal of the nodes of interest.
  //DMHMOD

  if(nodeset_list.size() > 0)nodeset_vectors = new std::vector <long long> [nodeset_list.size()];

  std::list < PG_BC_Specification *> ::iterator setit;
  long long nsct = 0;
  for(setit = nodeset_list.begin(); setit != nodeset_list.end();setit++,nsct ++){
    LoopLimits ll = (*setit)->limits;
    Topo_Loc the_location = (*setit)->location;

    std::list < long long > nodes_vector;

    if(dimension == 3){
      for ( long long _nk_ = ll.ks; _nk_ < ll.ke; _nk_ ++){ 
	for ( long long _nj_ = ll.js; _nj_ < ll.je; _nj_ ++){ 
	  for ( long long _ni_ = ll.is; _ni_ < ll.ie; _ni_ ++) {
	    long long global_node_id = get_node_number_from_l_i_j_k(trisection_blocks,_ni_,_nj_,_nk_);
	    for(unsigned the_nct = 0; the_nct < global_node_vector.size();the_nct++){
	      if(global_node_id == global_node_vector[the_nct]){
		nodes_vector.push_back(global_node_id);
	      }
	    } 
	  }
	}
      }
    }
    else{
      long long _nk_ = 0;
      for ( long long _nj_ = ll.js; _nj_ < ll.je; _nj_ ++){ 
	for ( long long _ni_ = ll.is; _ni_ < ll.ie; _ni_ ++) {
	  long long global_node_id = get_node_number_from_l_i_j_k(trisection_blocks,_ni_,_nj_,_nk_);
	  for(unsigned the_nct = 0; the_nct < global_node_vector.size();the_nct++){
	    if(global_node_id == global_node_vector[the_nct]){
	      nodes_vector.push_back(global_node_id);
	    }
	  }
	}
      }
    }
    if((!(*setit)->block_boundary_set) || 
       ((*setit)->block_boundary_set && (((*setit)->block_id-1)%inline_bx == 0))){
      // need to do trisect blocks now
      for(long long i = 0; i < trisection_blocks; i ++){
	Topo_Loc new_location;
	long long is_nodes = 1;
	LoopLimits tll = get_tri_block_limits((*setit)->block_boundary_set,(*setit)->block_id,i,the_location,new_location,is_nodes);

	if(dimension == 3){
	  for ( long long _nk_ = tll.ks; _nk_ < tll.ke; _nk_ ++){ 
	    for ( long long _nj_ = tll.js; _nj_ < tll.je; _nj_ ++){ 
	      for ( long long _ni_ = tll.is; _ni_ < tll.ie; _ni_ ++) {
		
		long long global_node_id = get_node_number_from_l_i_j_k(i,_ni_,_nj_,_nk_);  
		for(unsigned the_nct = 0; the_nct < global_node_vector.size();the_nct++){
		  if(global_node_id == global_node_vector[the_nct]){
		    nodes_vector.push_back(global_node_id);
		  }
		}
	      }
	    }
	  }
	}
	else{ 
	  long long _nk_ = 0;
	  for ( long long _nj_ = tll.js; _nj_ < tll.je; _nj_ ++){ 
	    for ( long long _ni_ = tll.is; _ni_ < tll.ie; _ni_ ++) {
	      
	      long long global_node_id = get_node_number_from_l_i_j_k(i,_ni_,_nj_,_nk_);  
	      for(unsigned the_nct = 0; the_nct < global_node_vector.size();the_nct++){
		if(global_node_id == global_node_vector[the_nct]){
		  nodes_vector.push_back(global_node_id);
		}
	      }
	    }
	  }  
	}
      }
    }


    nodes_vector.sort();
    nodes_vector.unique();

    std::list < long long > :: iterator nit;
    for(nit = nodes_vector.begin(); nit != nodes_vector.end(); nit ++){
      for(unsigned the_nct = 0; the_nct < global_node_vector.size();the_nct++){
        if((*nit) == global_node_vector[the_nct])nodeset_vectors[nsct].push_back((*nit));
      }
    }
  }
  // END OF NODESETS

  //SIDESETS
  // Sidesets are defined by the element id and the id of the face to which the set
  // applies.

  // For the purposes of calculation the sideset definition process starts by 
  // considering all elements to be numbered sequentially in i then j then k
  // across all blocks. This is not actually the case since exodus requires elements
  // to be numbered sequentially within a block.

  // If we consider the elements to be numbered sequentially across blocks in i,j, and k
  // then we can use a loop limits type call to get the i,j,and k index limits of the
  // elements in the sideset. We can then loop over these limits to calculate the 
  // index of the element if it were numbered across blocks. It is then left to 
  // calculate the index of this element if it were numbered sequentially within blocks.
  // This is done by calculating the global i,j,and k indices of the element within the 
  // entire domain, calculating the indices of the block within which the element resides,
  // calculating the ordinal number of the block in which the element resides, and the 
  // local i,j,and k indices of the element within the block it resides.

  //  These values are combined to calculate the index of the element that corresponds
  // to numbering the elements sequentially within blocks. 
  nsct = 0;
  if(sideset_list.size() > 0)sideset_vectors = new std::vector < std::pair <long long ,Topo_Loc > > [sideset_list.size()];  

  for(setit = sideset_list.begin(); setit != sideset_list.end();setit++,nsct ++){
     Topo_Loc the_location = (*setit)->location;
    //Sidesets allowed only on faces of block, not on edges or corners

    LoopLimits ll = (*setit)->limits;
	if(dimension == 3){
	  for ( long long _nk_ = ll.ks; _nk_ < ll.ke; _nk_ ++){ 
	    for ( long long _nj_ = ll.js; _nj_ < ll.je; _nj_ ++){ 
	      for ( long long _ni_ = ll.is; _ni_ < ll.ie; _ni_ ++) {
		long long elnumber = get_element_number_from_l_i_j_k(trisection_blocks,_ni_,_nj_,_nk_);
		long long the_proc_id = Element_Proc(elnumber);
		if(the_proc_id == my_rank){// add only if on this proc
		  std::pair <long long ,Topo_Loc > el_loc_pair(elnumber,the_location);
		  sideset_vectors[nsct].push_back(el_loc_pair);
		}
	      }
	    }
	  }
	}
	else{
	  long long _nk_ = 0;
	  for ( long long _nj_ = ll.js; _nj_ < ll.je; _nj_ ++){ 
	    for ( long long _ni_ = ll.is; _ni_ < ll.ie; _ni_ ++) {
	      long long elnumber = get_element_number_from_l_i_j_k(trisection_blocks,_ni_,_nj_,_nk_);
	      long long the_proc_id = Element_Proc(elnumber);
	      if(the_proc_id == my_rank){// add only if on this proc
		std::pair <long long ,Topo_Loc > el_loc_pair(elnumber,the_location);
		sideset_vectors[nsct].push_back(el_loc_pair);
	      }
	    }
	  }
	}
	
	if((!(*setit)->block_boundary_set) || 
	   ((*setit)->block_boundary_set && (((*setit)->block_id-1)%inline_bx == 0))){
	// need to do trisect blocks now
	for(long long i = 0; i < trisection_blocks; i ++){
	  Topo_Loc new_location;
	  long long is_nodes = 0;
	  LoopLimits tll = get_tri_block_limits((*setit)->block_boundary_set,(*setit)->block_id,i,the_location,new_location,is_nodes);
	  
	  if(dimension == 3){
	    for ( long long _nk_ = tll.ks; _nk_ < tll.ke; _nk_ ++){ 
	      for ( long long _nj_ = tll.js; _nj_ < tll.je; _nj_ ++){ 
		for ( long long _ni_ = tll.is; _ni_ < tll.ie; _ni_ ++) {
		  long long elnumber = get_element_number_from_l_i_j_k(i,_ni_,_nj_,_nk_);
		  long long the_proc_id = Element_Proc(elnumber);
		  if(the_proc_id == my_rank){// add only if on this proc	  
		    std::pair <long long ,Topo_Loc > el_loc_pair(elnumber,new_location);
		    sideset_vectors[nsct].push_back(el_loc_pair);
		  }
		}
	      }
	    }
	  }
	  else{
	    long long _nk_ = 0;
	    for ( long long _nj_ = tll.js; _nj_ < tll.je; _nj_ ++){ 
	      for ( long long _ni_ = tll.is; _ni_ < tll.ie; _ni_ ++) {
		long long elnumber = get_element_number_from_l_i_j_k(i,_ni_,_nj_,_nk_);
		long long the_proc_id = Element_Proc(elnumber);
		if(the_proc_id == my_rank){// add only if on this proc	  
		  std::pair <long long ,Topo_Loc > el_loc_pair(elnumber,new_location);
		  sideset_vectors[nsct].push_back(el_loc_pair);
		}
	      }
	    }
	  }
	}
	}
  }
  //END of SIDESETS
}

/****************************************************************************/
long long Radial_Trisection_Inline_Mesh_Desc::Rename_Block_BC_Sets()
/****************************************************************************/
{
  std::list < PG_BC_Specification * > ::iterator setit;
  for(setit = nodeset_list.begin(); setit != nodeset_list.end();setit++){
    if((*setit)->block_boundary_set){
      long long bid = (*setit)->block_id-1;
      long long kind = bid/(blockKstride());
      long long jord = (bid-kind*blockKstride());
      //interior
      {
	long long jind = 0;
	if(jord && (inline_bx-1)){
	  jind = (jord-1)/(inline_bx-1);
	}
	long long iind = jord;
	if (jind){
	  iind =1+ ((jord-1)-jind*(inline_bx-1));
	}
	(*setit)->block_id = 1 + iind + jind*inline_bx + kind * inline_bx*inline_by;
      }
    }
  }
  for(setit = sideset_list.begin(); setit != sideset_list.end();setit++){
    if((*setit)->block_boundary_set){
      long long bid = (*setit)->block_id-1;
      long long kind = bid/(blockKstride());
      long long jord = (bid-kind*blockKstride());
      //interior
      {
	long long jind = 0;
	if(jord && (inline_bx-1)){
	  jind = (jord-1)/(inline_bx-1);
	}
	long long iind = jord;
	if (jind){
	  iind =1+ ((jord-1)-jind*(inline_bx-1));
	}
	(*setit)->block_id = 1 + iind + jind*inline_bx + kind * inline_bx*inline_by;
      }

    }
  }
  return 0;
}


/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::Populate_Connectivity(long long * const * conn_array,                                         
								     std::map <long long, long long> & global_node_map
)
/****************************************************************************/
{
  long long num_nodes_per_element = 4;
    
  if(dimension == 3){
    num_nodes_per_element = 8;
  }

  //Nodes are ordered across the entire block domain.
  //To calculate connectivity for a given element in a given block.
  //It is necessary to calculate the global element indices in i,j,k
  //that identify that element in the global space.
  long long total_element_count = 0;
  for(long long bct = 0; bct < numBlocks(); bct ++ ){
    long long * conn = conn_array[bct];
    //build connectivity for each element
    //nodes are numbered 1-tot_num_nodes+1 across all blocks
    //incrementing i fastest
    for(unsigned elct = 0;elct < element_block_lists[bct].size();elct++,total_element_count++){
      long long the_el = element_block_lists[bct][elct];
      long long l;
      long long i;
      long long j;
      long long k;

      get_l_i_j_k_from_element_number(the_el,l,i,j,k);
      long long nn0 = get_node_number_from_l_i_j_k(l, i + 0, j + 0, k + 0);
      long long nn1 = get_node_number_from_l_i_j_k(l, i + 1, j + 0, k + 0);
      long long nn2 = get_node_number_from_l_i_j_k(l, i + 1, j + 1, k + 0);
      long long nn3 = get_node_number_from_l_i_j_k(l, i + 0, j + 1, k + 0);

      conn[elct*num_nodes_per_element + 0 + 0] = get_map_entry(global_node_map,nn0)+1;
      conn[elct*num_nodes_per_element + 1 + 0] = get_map_entry(global_node_map,nn1)+1;
      conn[elct*num_nodes_per_element + 2 + 0] = get_map_entry(global_node_map,nn2)+1;
      conn[elct*num_nodes_per_element + 3 + 0] = get_map_entry(global_node_map,nn3)+1;
      if(dimension == 3){
	long long nn4 = get_node_number_from_l_i_j_k(l, i + 0, j + 0, k + 1);
	long long nn5 = get_node_number_from_l_i_j_k(l, i + 1, j + 0, k + 1);
	long long nn6 = get_node_number_from_l_i_j_k(l, i + 1, j + 1, k + 1);
	long long nn7 = get_node_number_from_l_i_j_k(l, i + 0, j + 1, k + 1);
	conn[elct*num_nodes_per_element + 0 + 4] = get_map_entry(global_node_map,nn4)+1;
	conn[elct*num_nodes_per_element + 1 + 4] = get_map_entry(global_node_map,nn5)+1;
	conn[elct*num_nodes_per_element + 2 + 4] = get_map_entry(global_node_map,nn6)+1;
	conn[elct*num_nodes_per_element + 3 + 4] = get_map_entry(global_node_map,nn7)+1;
      }
    }
  }
}

/****************************************************************************/
long long Radial_Trisection_Inline_Mesh_Desc::Populate_Sideset_Info(std::map <long long, long long> & global_element_map,
								     std::map <long long, long long> & global_node_map,
								     long long * const * side_set_elements,
								     long long * const * side_set_faces,
								     long long * const * side_set_nodes,
								     long long * const * side_set_node_counter)
/****************************************************************************/
{
  long long num_nodes_per_face = 2;
  if(dimension == 3){
    num_nodes_per_face = 4;
  }

  long long nsct = 0;
   std::list < PG_BC_Specification *> ::iterator setit;

  for(setit = sideset_list.begin(); setit != sideset_list.end();setit++,nsct ++){

    long long * the_elements = side_set_elements[nsct];
    long long * the_faces = side_set_faces[nsct];
    long long * the_nodes = side_set_nodes[nsct];
    long long * the_node_counter = side_set_node_counter[nsct];

    //Sidesets allowed only on faces of block, not on edges or corners


    for(unsigned elct = 0; elct < sideset_vectors[nsct].size();elct ++){
      long long the_element = sideset_vectors[nsct][elct].first;
      Topo_Loc the_location =  sideset_vectors[nsct][elct].second;
      // These are the indices of the element in the entire domain
      long long gl_l;
      long long gl_k;
      long long gl_j;
      long long gl_i;

      get_l_i_j_k_from_element_number(the_element,gl_l,gl_i,gl_j,gl_k);

      //The result of the calculation is the id of the element in the block as found in the connectivity array
      the_elements[elct] = get_map_entry(global_element_map,the_element)+1;
      the_faces[elct] = topo_loc_to_exo_face[the_location];
      the_node_counter[elct] = num_nodes_per_face;
      //Since the nodes are numbered across blocks, the global indices are the actual indices of the node
      //It is only required to permute the indices to collect the appropriate nodes for the indicated face and 
      //calculated the node index using the nodal stride values for the entire mesh and the indices

      // adjust for periodicity in j 
      long long glj_plus1 = gl_j + 1;
      if(periodic_j){
        if(glj_plus1 == nely_tot){
          glj_plus1 = 0;
        } 
      }

      long long gli_plus0 = gl_i;
      long long gli_plus1 = gl_i+1;
      long long glj_plus0 = gl_j;
      long long glk_plus0 = gl_k;
      long long glk_plus1 = gl_k+1;
      
      switch(the_location) {
	
      case MINUS_I:{
	if(dimension == 2){
	  the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus1, glk_plus0));
	  the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus0, glk_plus0));
	}
	else{
	  the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus0, glk_plus0));
	  the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus0, glk_plus1));
	  the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus1, glk_plus1));
	  the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus1, glk_plus0));
	}
	break;
      }
      case PLUS_I:{
	if(dimension == 2){
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus0, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus1, glk_plus0));
	}
	else{
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus0, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus1, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus1, glk_plus1));
	the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus0, glk_plus1));
	}
	break;
      }
      case MINUS_J:{
      if(dimension == 2){
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus0, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus0, glk_plus0));
      }
      else{
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus0, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus0, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus0, glk_plus1));
	the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus0, glk_plus1));
      }
	break;
      }
      case PLUS_J:{
      if(dimension == 2){
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus1, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus1, glk_plus0));
      }
      else{
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus1, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus1, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus1, glk_plus1));
	the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus1, glk_plus1));
      }
	break;
      }
      case MINUS_K:{
      if(dimension == 2){
      }
      else{
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus0, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus1, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus1, glk_plus0));
	the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus0, glk_plus0));
      }
      break;
      }
      case PLUS_K:{
	if(dimension == 2){
	}
	else{
	  the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus0, glk_plus1));
	  the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus0, glk_plus1));
	  the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus1, glj_plus1, glk_plus1));
	  the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,get_node_number_from_l_i_j_k(gl_l, gli_plus0, glj_plus1, glk_plus1));
	}
	break;		     
      }
      default:
	error_stream << "Inline_Mesh_Desc::Read_mesh(): "
		     << "Sideset applied to unknown Topo Location.";
	return 1;
	break;
      }   
    }
  }
  return 0;
  //END of SIDESETS
}

/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::Populate_Map_and_Global_Element_List(long long * the_map, 
							 long long * global_element_numbers)
/****************************************************************************/
{
  long long total_count = 0;
  for(long long bct = 0; bct < numBlocks();bct ++ ){
    for(unsigned ect = 0; ect < element_block_lists[bct].size();ect ++){
      long long the_el = element_block_lists[bct][ect];
      long long elid = GetBlockBasedGlobalID(the_el,bct);
      the_map[total_count] = elid + 1;
      global_element_numbers[total_count] = elid + 1;
      total_count ++;
    }
  }
}
/****************************************************************************/
long long Radial_Trisection_Inline_Mesh_Desc::GetBlockBasedGlobalID(long long the_el,long long bct)
/****************************************************************************/
{
  long long l;
  long long Ig;
  long long Jg;
  long long Kg;
  
  get_l_i_j_k_from_element_number(the_el,l,Ig,Jg,Kg);
  
  //global element indices
  
  
  long long Kbg = get_block_index(Kg,inline_bz,c_inline_nz);
  long long Jbg = 0;
  if(l == trisection_blocks){
    Jbg = get_block_index(Jg,inline_by,c_inline_ny);
  }
  long long Ibg = 0;
  if(l == trisection_blocks){
    Ibg = get_block_index(Ig,inline_bx,c_inline_nx);
  }
  
  //ordinal of the block
  long long mesh_block_id = Ibg + Jbg*inline_bx + Kbg*inline_bx*inline_by; 
  
  
  //indices inside the block
  long long Kblock = Kg-c_inline_nz[Kbg];
  long long Jblock = Jg-c_inline_ny[Jbg];
  long long Iblock = Ig-c_inline_nx[Ibg];
  
  long long elid = 0;
  if(l == trisection_blocks && (Ibg > 0)){
    elid = l*a_inline_nz[Kbg] * div * div;
    elid += cum_block_totals[mesh_block_id] + Iblock + Jblock*a_inline_nx[Ibg] + Kblock*a_inline_nx[Ibg]*a_inline_ny[Jbg];
  }
  else if (l == trisection_blocks){
    elid = cum_block_totals[mesh_block_id];
    elid += l*a_inline_nz[Kbg] * div * div;
    elid += Iblock + Jblock*a_inline_nx[Ibg] + Kblock*a_inline_nx[Ibg]*a_inline_ny[Jbg];
  }
  else{
    elid = cum_block_totals[mesh_block_id];
    elid += l*a_inline_nz[Kbg]*div * div;
    elid += Ig + Jg * div + Kblock*div*div;
  }

  return elid;
}

//! Queries which processor an element lies on.
//! Calls the recursive Partition::Element_Proc function.
/****************************************************************************/
long long Radial_Trisection_Inline_Mesh_Desc::Element_Proc(long long global_element_id)
/****************************************************************************/
{
  long long proc = 0;
  if(inline_decomposition_type == SEQUENTIAL){
    long long total = kestride * nelz_tot;
    long long num_per_proc = total/num_processors;
    proc = global_element_id/num_per_proc;
    if(proc >= num_processors)proc = num_processors-1;
  }
  else if(inline_decomposition_type == RANDOM){
    SRANDOM(global_element_id);
    long long rand_num = RANDOM();
    proc = rand_num%num_processors;
  }
  else if((inline_decomposition_type == BISECTION) || (inline_decomposition_type == PROCESSOR_LAYOUT)){
    long long l,i,j,k;
    get_l_i_j_k_from_element_number(global_element_id,l,i,j,k);
    long long ginds[4];

    ginds[0] = i;
    ginds[1] = j;
    ginds[2] = k;
    ginds[3] = l;
    if(l != trisection_blocks){
      ginds[0] = 0;//anywhere in this row would do
      ginds[1] = j+l*div*2;// this is the critical one
      ginds[3] = trisection_blocks;
    }
    return base_partition->Element_Proc(ginds);
  }

  return proc;
}


/****************************************************************************/
void Radial_Trisection_Inline_Mesh_Desc::getGlobal_Element_Block_Totals(long long * totals_array)
/****************************************************************************/
{
  //zero out array
  for(long long bct = 0; bct < numBlocks();bct ++ )totals_array[bct] = 0;
  
  //accumulate
  for(long long kct = 0; kct < inline_bz; kct ++){
    for(long long jct = 0; jct < inline_by; jct ++){
      for(long long ict = 1; ict < inline_bx; ict ++){
	long long block_num = kct*blockKstride();
	block_num += jct*(inline_bx-1);
	block_num += ict;
	totals_array[block_num] += a_inline_nx[ict]*a_inline_ny[jct]*a_inline_nz[kct];
      }
    }
  }
  for(long long kct = 0; kct < inline_bz; kct ++){
    for(long long ict = 0; ict < 1; ict ++){
      long long block_num = kct*blockKstride();
      totals_array[block_num] += a_inline_nx[ict]*c_inline_ny[inline_by]*a_inline_nz[kct];
      totals_array[block_num] += trisection_blocks*div*div*a_inline_nz[kct];
    }
  }
}


/****************************************************************************/
void  Radial_Trisection_Inline_Mesh_Desc::Calc_Parallel_Info(                                          
					  std::vector <long long> & element_vector,
					  std::vector<long long> & global_node_vector,   
					  std::map <long long, long long> & global_node_map,                          
					  std::list   <long long> & internal_node_list,	
					  std::list   <long long> & border_nodes_list,
					  std::list   <long long> & internal_element_list,
					  std::list   <long long> & border_elements_list,
					  std::list   <long long> & node_proc_id_list,
					  std::list   <long long> & element_proc_id_list,
					  std::vector <long long> & node_neighbor_vector,
					  std::list <long long> * &  boundary_node_list,
					  std::vector <long long> & element_neighbor_vector,
					  std::list <std::pair <long long ,Topo_Loc > > * & boundary_element_list)
  /****************************************************************************/
{

  Tel *el_array = NULL;
  if(element_vector.size()>0){
    el_array = new Tel[element_vector.size()];
    for(unsigned long long gev = 0; gev < element_vector.size(); gev ++){
      el_array[gev].global_id = element_vector[gev];
      el_array[gev].real_element = true;
    }
  }
  Tel * node_array = NULL;
  if(global_node_vector.size()>0){
    node_array = new Tel[global_node_vector.size()];
    
    for(unsigned long long gnv = 0;gnv < global_node_vector.size();gnv ++){
      node_array[gnv].global_id = global_node_vector[gnv];
    }
  }

  // walk all local elements
  // peer through their neighboring faces and ask if that neighbor is on my processor
  // if it is do nothing
  // if it is not
  //   I am a border
  //   all nodes on that face are border nodes 
  //   expand node_array and el_array objects with neighbor proc ids and topology directions

  // walk all local elements

  long long nfn = 2;
  long long nfaces = 4;
  long long nen = 1;
  
  if(dimension == 3){
    nfn = 4;
    nfaces = 6;
    nen = 2;
  }




  for(unsigned long long gev = 0; gev < element_vector.size(); gev ++){
    long long face_nodes_array[4];
    long long my_id = el_array[gev].global_id;
    long long ll,li,lj,lk;
    get_l_i_j_k_from_element_number(my_id,ll,li,lj,lk);
    el_array[gev].visits ++;
    
    for(long long face_count = 0; face_count < nfaces; face_count ++){
      Topo_Loc tl = (Topo_Loc)face_count;
      long long neighbor = get_neighbor(tl,ll,li,lj,lk);
      if(neighbor >= 0){
	long long neighbor_proc_id = Element_Proc(neighbor);
	if(neighbor_proc_id != my_rank){
	  std::pair < long long,Topo_Loc> conn_pair(neighbor_proc_id,tl);
	  el_array[gev].conn_connections.push_back(conn_pair);

	  el_array[gev].visits ++;
	  get_face_nodes(tl,my_id,face_nodes_array);
     
	  for(long long fnc = 0; fnc < nfn; fnc ++){
	    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].visits ++;
	    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].proc_neighbors.push_back(neighbor_proc_id);
	  }
	}
      }
    }


    // need to do edges and vertex neighbors for nodes
    long long edge_start = EDGE4;
    long long edge_end = EDGE8;
    if(dimension == 3){
      edge_start = EDGE0;
      edge_end = VERTEX0;
    }
    for(long long edge_count = edge_start; edge_count < edge_end; edge_count ++){
      Topo_Loc tl = (Topo_Loc)edge_count;
      long long neighbor = get_neighbor(tl,ll,li,lj,lk);
      if(neighbor >= 0){
	long long neighbor_proc_id = Element_Proc(neighbor);
	if(neighbor_proc_id != my_rank){
	  get_face_nodes(tl,my_id,face_nodes_array);
	  
	  for(long long fnc = 0; fnc < nen; fnc ++){
	    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].visits ++;
	    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].proc_neighbors.push_back(neighbor_proc_id);
	  }
	}
      }
    }
    if(ll < trisection_blocks){
      Topo_Loc ttl = (Topo_Loc)EDGE4;
      if(li == 0 && lj == 0){// means we are at an apex element of transition block
	// run through the possible apex elements
	for(long long tll = 0; tll < trisection_blocks; tll ++){
	  long long tneighbor = get_element_number_from_l_i_j_k(tll,0,0,lk);	  
	  if(tneighbor >=0){
	    long long tneighbor_proc_id = Element_Proc(tneighbor);

	    if(tneighbor_proc_id != my_rank){

	      get_face_nodes(ttl,my_id,face_nodes_array);
	      for(long long fnc = 0; fnc < nen; fnc ++){
		node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].visits ++;
		node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].proc_neighbors.push_back(tneighbor_proc_id);
	      }
	    }
	  }
	}
      }
    }
    if(dimension == 3){
      // need to do vertices and vertex neighbors for nodes
      for(long long vertex_count = EDGE11; vertex_count < NUM_TOPO_CONNECTIONS; vertex_count ++){
	Topo_Loc tl = (Topo_Loc)vertex_count;
	long long neighbor = get_neighbor(tl,ll,li,lj,lk);
	if(neighbor >= 0){
	  long long neighbor_proc_id = Element_Proc(neighbor);
	  if(neighbor_proc_id != my_rank){
	    get_face_nodes(tl,my_id,face_nodes_array);
	    for(long long fnc = 0; fnc < 1; fnc ++){
	      node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].visits ++;
	      node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].proc_neighbors.push_back(neighbor_proc_id);
	    }
	  }
	}
      }
    if(ll < trisection_blocks){// loop over the apex elements associated  with the two cornere nodes
      if(li == 0 && lj == 0){// means we are at an apex element of trisection block
	Topo_Loc tta[2];
	tta[0] = (Topo_Loc)VERTEX0;//bottom of apex edge
	tta[1] = (Topo_Loc)VERTEX4;//top of apex edge
	for(long long ict = 0; ict < 2; ict ++){
	  Topo_Loc ttl = tta[ict];
	  long long tkstart = lk;
	  long long tkend = lk + 1;
	  if(ict == 0)tkstart --;// lower node use lower and this level elements
	  if(ict == 1)tkend ++;// upper node use this and upper level elements
	  // run through the possible apex elements
	  for(long long tll = 0; tll < trisection_blocks; tll ++){
	    for(long long tlk = tkstart; tlk < tkend; tlk ++){
	      long long tneighbor = get_element_number_from_l_i_j_k(tll,0,0,tlk);
	      if(tneighbor >=0){
		long long tneighbor_proc_id = Element_Proc(tneighbor);		
		if(tneighbor_proc_id != my_rank){
		  
		  get_face_nodes(ttl,my_id,face_nodes_array);
		  for(long long fnc = 0; fnc < 1; fnc ++){
		    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].visits ++;
		    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].proc_neighbors.push_back(tneighbor_proc_id);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    }
  }

  for(unsigned long long i = 0; i < element_vector.size();i ++){
    if(el_array[i].visits > 1){
      // loop over all conn_connections
      std::list < std::pair < long long , Topo_Loc > > ::iterator conit;
      for(conit  = el_array[i].conn_connections.begin();
          conit != el_array[i].conn_connections.end();
          conit ++){
        element_proc_id_list.push_back((*conit).first);
      }
    }
  }
  // sort and uniq element_proc_id_list
  element_proc_id_list.sort();
  element_proc_id_list.unique();

  if(element_proc_id_list.size()){
    boundary_element_list = new std::list < std::pair <long long ,Topo_Loc > > [element_proc_id_list.size()];
  }

  std::map <long long,long long> element_neighbor_proc_map; //key is proc_id value is ordinal
  std::list <long long> ::iterator listit;
  long long the_count = 0;
  for(listit = element_proc_id_list.begin(); listit != element_proc_id_list.end(); listit++,the_count ++){
    element_neighbor_proc_map[*listit] = the_count;
    element_neighbor_vector.push_back(*listit);
  }

  // now populate the maps

  for(unsigned long long i = 0; i < element_vector.size();i ++){
    long long the_element = element_vector[i];
    if(el_array[i].visits == 1){
      internal_element_list.push_back(the_element);
    }
    if(el_array[i].visits > 1){
      // loop over all conn_connections
      
      border_elements_list.push_back(the_element);
      std::list < std::pair < long long , Topo_Loc > > ::iterator conit;
      for(conit  = el_array[i].conn_connections.begin();
          conit != el_array[i].conn_connections.end();
          conit ++){
	
        long long index = get_map_entry(element_neighbor_proc_map,(*conit).first);
        boundary_element_list[index].push_back( std::pair <long long ,Topo_Loc >(the_element,(*conit).second));
      }
    }
  }

  border_elements_list.sort();
  border_elements_list.unique();


  for(unsigned long long gnv = 0;gnv < global_node_vector.size();gnv ++){
    if(node_array[gnv].visits > 0){
      // loop over all conn_connections
      std::list < long long > ::iterator conit;
      for(conit  = node_array[gnv].proc_neighbors.begin();
	  conit != node_array[gnv].proc_neighbors.end();
	  conit ++){
	node_proc_id_list.push_back((*conit));
      }
    }
  }
    
  node_proc_id_list.sort();
  node_proc_id_list.unique();

  std::map <long long,long long> node_neighbor_proc_map; //key is proc_id value is ordinal
  std::list <long long> ::iterator nlistit;
  the_count = 0;
  for(nlistit = node_proc_id_list.begin(); nlistit != node_proc_id_list.end(); nlistit++,the_count ++){
    node_neighbor_proc_map[*nlistit] = the_count;
    node_neighbor_vector.push_back(*nlistit);
  }

  if(node_proc_id_list.size()){
    boundary_node_list = new std::list <long long> [node_proc_id_list.size()];
  }



  //node array needs global_id!!!!
  for(unsigned long long i = 0;i < global_node_vector.size();i ++){
    if(node_array[i].visits == 0){
      long long the_node = node_array[i].global_id;
      internal_node_list.push_back(the_node);
    }
    else if(node_array[i].visits > 0){
      long long the_node = node_array[i].global_id;
      // loop over all conn_connections
      std::list < long long > ::iterator conit;
      for(conit  = node_array[i].proc_neighbors.begin();
	  conit != node_array[i].proc_neighbors.end();
	  conit ++){
	long long index = get_map_entry(node_neighbor_proc_map,(*conit));
	boundary_node_list[index].push_back(the_node);
	border_nodes_list.push_back(the_node);
      }
    }
  }
  // sort the boundary_node_list
  for(unsigned long long i = 0; i < node_proc_id_list.size();i++){
    boundary_node_list[i].sort();
    boundary_node_list[i].unique();    
  }
  border_nodes_list.sort();
  border_nodes_list.unique();

  
  
  
  delete [] el_array;
  delete[] node_array;

  //check number node comm_maps
}

/****************************************************************************/
long long Radial_Trisection_Inline_Mesh_Desc::GlobalNumElements()
/****************************************************************************/
{
  long long total_count = 0;
  //accumulate
  for(long long kct = 0; kct < inline_bz; kct ++){
    for(long long jct = 0; jct < inline_by; jct ++){
      for(long long ict = 1; ict < inline_bx; ict ++){
	total_count += a_inline_nx[ict]*a_inline_ny[jct]*a_inline_nz[kct];
      }
    }
  }
  for(long long kct = 0; kct < inline_bz; kct ++){
    for(long long ict = 0; ict < 1; ict ++){
      total_count += a_inline_nx[ict]*c_inline_ny[inline_by]*a_inline_nz[kct];
      total_count += trisection_blocks*div*div*a_inline_nz[kct];
    }
  }
  return total_count;
}
}//end namespace PAMGEN_NEVADA
