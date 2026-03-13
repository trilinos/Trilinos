// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "inline_mesh_desc.h"
#include "radial_inline_mesh_desc.h"
#include "uns_inline_decomp.h"
#include "pamgen_fudges.h"
#include "Vector.h"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace PAMGEN_NEVADA {

  /*****************************************************************************/
  void Radial_Inline_Mesh_Desc::calculateSize(
      long long & total_el_count,
      long long & total_node_count,
      long long & total_edge_count)
  {
    total_el_count   = 0L;
    total_node_count = 0L;
    total_edge_count = 0L;

    if(dimension == 3) {
      for(long long k = 0; k < inline_b[2]; k ++) {
        for(long long j = 0; j < inline_b[1]; j ++) {
          for(long long i = 0; i < inline_b[0]; i ++) {
            total_el_count +=
              (long long)interval[0][i]*
              (long long)interval[1][j]*
              (long long)interval[2][k];
          }
        }
      }
    } else {
      for(long long j = 0; j < inline_b[1]; j ++){
        for(long long i = 0; i < inline_b[0]; i ++){
          total_el_count +=
            (long long)interval[0][i]*
            (long long)interval[1][j];
        }
      }
    }

    long long nodes_temp[3] = {};

    for(long long axis = 0; axis < dimension;axis ++){
      nodes_temp[axis] = 0;
      for(long long i = 0; i < inline_b[axis]; i ++){
        nodes_temp[axis] += interval[axis][i];
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

    for(long long axis = 0; axis < 3;axis ++){
      a_inline_n[axis] = new long long [inline_b[axis]];
      c_inline_n[axis] = new long long [inline_b[axis]+1];
    }

    //Defaults for 3D quantities in 2D mesh
    nel_tot[2] = 1;
    a_inline_n[2][0] = 1;
    c_inline_n[2][0] = 0;
    c_inline_n[2][1] = 1;

    for(long long axis = 0; axis < dimension; axis ++){
      nel_tot[axis] = 0;
      for(long long i = 0; i < inline_b[axis]; i ++){
        a_inline_n[axis][i] = inline_n[axis];
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


    for(long long axis = 0; axis < dimension; axis ++){
      inline_gmax[axis] = inline_gmin[axis];
      for(long long i = 0; i < inline_b[axis]; i ++){
        inline_gmax[axis] += block_dist[axis][i];
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
  void Radial_Inline_Mesh_Desc::Populate_Coords(Real * coords,
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
  Vector Radial_Inline_Mesh_Desc::calc_coords_periodic(
      double total_theta,
      long long i,
      long long j,
      long long k)
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


/****************************************************************************/
KOKKOS_INLINE_FUNCTION
void Radial_Inline_Mesh_Desc::calc_coords_periodic_device(double total_theta,
                                                        long long i, long long j, long long k,
                                                        double& x, double& y, double& z) const
{
  // Device version of periodic coordinate calculation
  long long per = 0;
  if (total_theta == 90.0) per = nel_tot[1] / 2;
  if (total_theta == 180.0) per = nel_tot[1] / 4;
  if (total_theta == 360.0) per = nel_tot[1] / 8;

  long long jmod = j % per;
  long long jmult = j / per;

  double deg_to_rad = M_PI / 180.0;
  double theta = IJKcoors[1][jmod] * deg_to_rad;
  double r = IJKcoors[0][i];

  // Transform to appropriate octant
  if (jmult == 1) {
    x = r * cos(M_PI / 2.0 - theta);
    y = r * sin(M_PI / 2.0 - theta);
  } else if (jmult == 2) {
    x = -r * cos(theta);
    y = r * sin(theta);
  } else if (jmult == 3) {
    x = -r * cos(M_PI / 2.0 - theta);
    y = -r * sin(M_PI / 2.0 - theta);
  } else if (jmult == 4) {
    x = r * cos(theta);
    y = -r * sin(theta);
  } else if (jmult == 5) {
    x = r * cos(theta - M_PI / 2.0);
    y = r * sin(theta - M_PI / 2.0);
  } else if (jmult == 6) {
    x = r * cos(M_PI - theta);
    y = r * sin(M_PI - theta);
  } else if (jmult == 7) {
    x = r * cos(M_PI + theta);
    y = r * sin(M_PI + theta);
  } else {
    x = r * cos(theta);
    y = r * sin(theta);
  }

  if (dimension == 3) {
    z = IJKcoors[2][k];
  } else {
    z = 0.0;
  }
}

/****************************************************************************/
void Radial_Inline_Mesh_Desc::Populate_Coords_Device(View2D<double> coords,
                                                   const View1D<long long> global_node_vector,
                                                   const View1D<long long> global_node_map_keys,
                                                   const View1D<long long> global_node_map_values,
                                                   long long num_nodes)
/****************************************************************************/
{
  double deg_to_rad = M_PI / 180.0;
  double total_theta = c_block_dist[1][inline_b[1]];

  Kokkos::parallel_for("RadialPopulateCoords", global_node_vector.size(),
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

      // Calculate coordinates using radial transformation
      double x = IJKcoors[0][global_ind[0]] * cos(IJKcoors[1][global_ind[1]] * deg_to_rad);
      double y = IJKcoors[0][global_ind[0]] * sin(IJKcoors[1][global_ind[1]] * deg_to_rad);
      double z = 0.0;

      if (dimension == 3) {
        z = IJKcoors[2][global_ind[2]];
      }

      // Apply periodic boundary conditions if enabled
      if (enforce_periodic) {
        calc_coords_periodic_device(total_theta, global_ind[0], global_ind[1], global_ind[2], x, y, z);
      }

      // Set coordinates
      coords(the_local_node, 0) = x;
      coords(the_local_node, 1) = y;
      if (dimension == 3) {
        coords(the_local_node, 2) = z;
      }
    });

  // Ensure all device operations are complete
  Kokkos::fence();
}


}// end namespace PAMGEN_NEVADA
