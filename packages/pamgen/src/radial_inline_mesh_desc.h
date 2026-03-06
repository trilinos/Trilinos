// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef radial_inline_mesh_descH
#define radial_inline_mesh_descH

#include "inline_mesh_desc.h"
#include "Vector.h"
#include "pamgen_kokkos_utils.h"

namespace PAMGEN_NEVADA {

struct Vector;
class Radial_Inline_Mesh_Desc : public Inline_Mesh_Desc
{
public:
  Radial_Inline_Mesh_Desc(long long dim){dimension = dim;};
  Radial_Inline_Mesh_Desc(){};
  virtual ~Radial_Inline_Mesh_Desc(){};
  virtual long long Set_Up();
  virtual std::string Calc_Intervals();
  virtual void calculateSize(long long & total_el_count,
                             long long & total_node_count,
                             long long & total_edge_count);
//   virtual long long Calc_Coord_Vectors();
  Vector calc_coords_periodic(double total_theta,
                              long long i,
                              long long j,
                              long long k);
  virtual void Populate_Coords(double * coords,
                       std::vector<long long> & global_node_vector,
                       std::map <long long, long long> & global_node_map,
                       long long num_nodes);

  // Kokkos version of Populate_Coords
  void Populate_Coords_Device(View2D<double> coords,
                             const View1D<long long> global_node_vector,
                             const View1D<long long> global_node_map_keys,
                             const View1D<long long> global_node_map_values,
                             long long num_nodes);

  // Device function for periodic coordinate calculation
  KOKKOS_INLINE_FUNCTION
  void calc_coords_periodic_device(double total_theta,
                                  long long i, long long j, long long k,
                                  double& x, double& y, double& z) const;
};
}// end namespace
#endif
