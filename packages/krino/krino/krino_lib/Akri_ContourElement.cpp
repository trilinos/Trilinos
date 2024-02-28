// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ContourElement.hpp>
#include <Akri_ContourSubElement.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshHelpers.hpp>
#include <math.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include "Akri_LevelSet.hpp"

namespace krino{

ContourElement::ContourElement( const stk::mesh::BulkData & mesh,
		  stk::mesh::Entity element,
		  const FieldRef coords_field,
                  const FieldRef dist_field,
		  const double iso_dist )
    : my_spatial_dim( mesh.mesh_meta_data().spatial_dimension() ),
    my_coords_master_elem(MasterElementDeterminer::getMasterElement(mesh.bucket(element), coords_field)),
    my_dist_master_elem(MasterElementDeterminer::getMasterElement(mesh.bucket(element), dist_field)),
    my_sign(0),
    my_length_scale(-1.0),
    my_edge_linear_tolerance(-1.0),
    my_edge_nonlinear_tolerance(-1.0)
{ /* %TRACE% */  /* %TRACE% */

  //
  // Gather element's coordinates and distance into ArrayContainer's.
  //

  const int dim = my_spatial_dim;

  const int npe_coords = my_coords_master_elem.get_topology().num_nodes();
  const int npe_dist   = my_dist_master_elem.get_topology().num_nodes();

  my_coords.resize( dim, npe_coords );
  my_dist.resize( npe_dist );

  const stk::mesh::Entity* elem_nodes = mesh.begin_nodes(element);

  for ( int i = 0; i < npe_coords; ++i )
  {
    double * var = field_data<double>( coords_field, elem_nodes[i]);
    for ( int d = 0; d < dim; d++ )
    {
      my_coords(d,i) = var[d];
    }
  }

  for ( int i = 0; i < npe_dist; ++i )
  {
    double * var = field_data<double>( dist_field, elem_nodes[i]);
    my_dist(i) = *var - iso_dist;
  }
}

ContourElement::~ContourElement() {}

std::string
ContourElement::debug_output() const
{ /* %TRACE% */  /* %TRACE% */
  const int dim = my_spatial_dim;
  std::ostringstream os;

  os << "Element description:" << std::endl
     << "  Coordinates topolology is " << coord_topology().name()
     << ", Distance topology is " << dist_topology().name() << std::endl;

  int lnn = 0;
  int coord_counter = 0;
  int dist_counter = 0;
  const unsigned npe_coords = my_coords_master_elem.get_topology().num_nodes();
  for (unsigned node_index=0; node_index<npe_coords; ++node_index)
  {
    os << "  Node #" << lnn++;
    stk::math::Vector3d coords(stk::math::Vector3d::ZERO);
    for ( int d = 0; d < dim; d++ )
      coords[d] = my_coords(d,coord_counter);
    os << ", coords = (" << coords[0] << ","
       << coords[1] << ","
       << coords[2] << ")";
    coord_counter++;

    os << ", dist = " << my_dist(dist_counter++);
    os << std::endl;
  }
  return os.str();
}

double
ContourElement::distance( const stk::math::Vector3d & p_coords ) const
  { /* %TRACE% */  /* %TRACE% */
    double dist;
    my_dist_master_elem.interpolate_point(my_spatial_dim, p_coords.data(), 1, my_dist.ptr(), &dist);
    return dist;
  }

stk::math::Vector3d
ContourElement::coordinates( const stk::math::Vector3d & p_coords ) const
  { /* %TRACE% */  /* %TRACE% */
    stk::math::Vector3d coords(stk::math::Vector3d::ZERO);
    my_coords_master_elem.interpolate_point(my_spatial_dim, p_coords.data(), my_spatial_dim, my_coords.ptr(), coords.data());
    return coords;
  }

double ContourElement::compute_domain_integral(const sierra::ArrayContainer<double, DIM, NINT> & intgPtLocations,
    const sierra::ArrayContainer<double, NINT> & intgWeights,
    const sierra::ArrayContainer<double, NINT> &  determinants)
{
  const int numIntgPts = intgWeights.dimension<0>();
  double domainIntegral = 0.;
  for (int ip = 0; ip < numIntgPts; ++ip)
    domainIntegral += intgWeights(ip) * determinants(ip);
  return domainIntegral;
}

double ContourElement::compute_domain_integral(const int signOfDomain) const // 0=interface, -1=negVol, 1=posVol
{
  if (my_sign != 0)
  {
    if (my_sign == signOfDomain)
    {
      return volume();
    }
    return 0;
  }
  STK_ThrowErrorMsgIf(!my_base_subelement, "\ncompute_subelement_decomposition(...) must be called prior to calling gather_intg_pts().");

  if (signOfDomain == 0)
  {
    return my_base_subelement->compute_area_of_interface();
  }

  return my_base_subelement->compute_relative_signed_volume(signOfDomain) * volume() / my_coords_master_elem.parametric_volume();
}

double ContourElement::compute_area_of_interface() const
{
  return compute_domain_integral(0);
}

double ContourElement::compute_signed_volume(const int signOfVolume) const
{
  return compute_domain_integral(signOfVolume);
}

double
ContourElement::determinant( const stk::math::Vector3d & p_coords ) const
  { /* %TRACE% */  /* %TRACE% */
    double detJ, detJ_error;

    const int num_coord_dofs = coord_topology().num_nodes();

    // temporary data
    static sierra::ArrayContainer<double,DIM,NPE_COORD,NINT> d_shapef_coords;

    // reserve enough space for results
    d_shapef_coords.resize(my_spatial_dim, num_coord_dofs, 1);

    my_coords_master_elem.shape_fcn_deriv(1, p_coords.data(), d_shapef_coords.ptr());
    my_coords_master_elem.determinant(
        my_spatial_dim,       // Number of coordinate dimensions
        1,                    // Number of target points
        num_coord_dofs,       // Number of coord shape functions
        d_shapef_coords.ptr(),// Mesh shape function derivatives
        1,                    // Number of elements
        my_coords.ptr(),      // Mesh coordinate values
        &detJ,                // Determinant of the transformation Jacobian for each element (output)
        &detJ_error );        // Determinant error (output)
    return detJ;
  }

stk::math::Vector3d
ContourElement::distance_gradient( const stk::math::Vector3d & p_coords ) const
  { /* %TRACE% */  /* %TRACE% */
    const int dim = my_spatial_dim;
    const int num_coord_dofs = my_coords.dimension<NPE_COORD>();
    const int num_dist_dofs = my_dist.dimension<NPE_VAR>();

    // temporary data
    static sierra::ArrayContainer<double,DIM,NPE_COORD> d_shapef_coords;
    static sierra::ArrayContainer<double,DIM,NPE_VAR>   d_shapef_dist;
    static sierra::ArrayContainer<double,DIM,NPE_VAR>   grad_op;

    // reserve enough space for results
    d_shapef_coords.resize(dim, num_coord_dofs);
    d_shapef_dist.resize(dim, num_dist_dofs);
    grad_op.resize(dim, num_dist_dofs);
    stk::math::Vector3d grad_dist( stk::math::Vector3d::ZERO );
    double det_J;

    // pointers to array data
    const double * intg_pt_loc_ptr = p_coords.data();
    double * grad_dist_ptr = grad_dist.data();

    double gradop_error;

    my_coords_master_elem.shape_fcn_deriv(1, intg_pt_loc_ptr, d_shapef_coords.ptr());
    my_dist_master_elem.shape_fcn_deriv(1, intg_pt_loc_ptr, d_shapef_dist.ptr());

    my_dist_master_elem.gradient_operator(
        dim,                    // Number of coordinate dimensions
        1,                      // Number of target points
	num_coord_dofs,         // Number of coord shape functions
	d_shapef_coords.ptr(),// Mesh shape function derivatives
	num_dist_dofs,          // Number of dof shape functions
	d_shapef_dist.ptr(),    // Dof shape function derivatives
	1,                      // Number of elements
	my_coords.ptr(),        // Mesh coordinate values
	grad_op.ptr(),  // Gradient operator values (output)
	&det_J,                 // Determinant of the transformation Jacobian for each element (output)
	&gradop_error );                // Gradop error (output)

    my_dist_master_elem.scalar_gradient(
        1,                // Number of target points
	1,                // Number of elements
	grad_op.ptr(),    // Gradient operator values
	&det_J,           // Determinant of the transformation Jacobian for each element
	my_dist.ptr(),    // nodal distance
	grad_dist_ptr );  // Gradient of distance at integration pts (output)
    return( grad_dist );
  }

void
ContourElement::compute_distance_gradient( const sierra::Array<const double,DIM,NINT> & intg_pt_locations,
				    sierra::ArrayContainer<double,DIM,NINT> & grad_dist ) const
  { /* %TRACE% */  /* %TRACE% */
    const int dim          = intg_pt_locations.dimension<DIM>();
    const int num_intg_pts = intg_pt_locations.dimension<NINT>();

    const int num_coord_dofs = my_coords.dimension<NPE_COORD>();
    const int num_dist_dofs = my_dist.dimension<NPE_VAR>();

    // temporary data
    static sierra::ArrayContainer<double,DIM,NPE_COORD,NINT> d_shapef_coords;
    static sierra::ArrayContainer<double,DIM,NPE_VAR,NINT>   d_shapef_dist;
    static sierra::ArrayContainer<double,NINT>               det_J;
    static sierra::ArrayContainer<double,DIM,NPE_VAR,NINT>   grad_op;

    // reserve enough space for results
    d_shapef_coords.resize(dim, num_coord_dofs, num_intg_pts);
    d_shapef_dist.resize(dim, num_dist_dofs, num_intg_pts);
    det_J.resize(num_intg_pts);
    grad_op.resize(dim, num_dist_dofs, num_intg_pts);
    grad_dist.resize(dim, num_intg_pts);

    double gradop_error;

    my_coords_master_elem.shape_fcn_deriv(num_intg_pts, intg_pt_locations.ptr(), d_shapef_coords.ptr());
    my_dist_master_elem.shape_fcn_deriv(num_intg_pts, intg_pt_locations.ptr(), d_shapef_dist.ptr());

    my_dist_master_elem.gradient_operator(
        dim,                    // Number of coordinate dimensions
        num_intg_pts,           // Number of target points
	num_coord_dofs,         // Number of coord shape functions
	d_shapef_coords.ptr(),  // Mesh shape function derivatives
	num_dist_dofs,          // Number of dof shape functions
	d_shapef_dist.ptr(),    // Dof shape function derivatives
	1,                      // Number of elements
	my_coords.ptr(),        // Mesh coordinate values
	grad_op.ptr(),          // Gradient operator values (output)
	det_J.ptr(),            // Determinant of the transformation Jacobian for each element (output)
	&gradop_error );                // Gradop error (output)

    my_dist_master_elem.scalar_gradient(
        num_intg_pts,           // Number of target points
	1,                      // Number of elements
	grad_op.ptr(),          // Gradient operator values
	det_J.ptr(),            // Determinant of the transformation Jacobian for each element
	my_dist.ptr(),          // nodal distance
	grad_dist.ptr() );      // Gradient of distance at integration pts (output)
  }

template<int NDIM, int NNODES>
std::array<stk::math::Vector3d,NNODES> get_nodal_parametric_coords_array(const MasterElement & me)
{
  std::array<stk::math::Vector3d,NNODES> paramCoords;
  const double * paramCoordsArray = me.nodal_parametric_coordinates();
  int paramCoordsCounter = 0;
  for ( int i = 0; i < NNODES; ++i )
  {
    for ( int j = 0; j < NDIM; j++ )
      paramCoords[i][j] = paramCoordsArray[paramCoordsCounter++];
    for ( int j = NDIM; j < 3; j++ )
      paramCoords[i][j] = 0.;
  }
  return paramCoords;
}

template<int NSIDES>
std::array<int,NSIDES> get_side_ids_array()
{
  std::array<int,NSIDES> sideIds;
  for ( int i = 0; i < NSIDES; ++i )
    sideIds[i] = i;
  return sideIds;
}

int ContourElement::compute_conservative_nonlinear_distance_sign(const double snapTol, const unsigned numDist, const double * dist)
{
  // For a nonlinear distance function, estimate that interpolated values within
  // the element lie between min-variation < values < max+variation
  // where min and max are min and max nodal values and varation = max-min.

  // TODO: Could probably be less conservative, by shrinking the range in this nonlinear estimate.
  //       But what is the range of the interpolant in a quadratic element?

  double maxDist = dist[0];
  double minDist = dist[0];
  for ( unsigned n = 1; n < numDist; n++ )
  {
    if (dist[n] < minDist) minDist = dist[n];
    if (dist[n] > maxDist) maxDist = dist[n];
  }

  // Incorporate snapping into range over element
  if (minDist > 0. && minDist < snapTol) minDist = 0.;
  if (maxDist < 0. && -maxDist < snapTol) maxDist = 0.;

  const double variation = maxDist - minDist;

  const double estimatedMinInterpolant = minDist - variation;
  const double estimatedMaxInterpolant = maxDist + variation;

  const bool all_hi = estimatedMinInterpolant > 0.0;
  const bool all_lo = estimatedMaxInterpolant < 0.0;

  if (all_hi || all_lo)
  {
    return LevelSet::sign(minDist);
  }

  return 0;
}

int ContourElement::compute_linear_distance_sign(const double snapTol, const unsigned numDist, const double * dist)
{
  // Here "linear" means that the extrema within the element matches the extrema of the nodal values,
  // which is true for bilinear quads (and hex) as well as truly linear tets and tris.

  bool hasNeg = false;
  bool hasPos = false;
  for ( unsigned n = 0; n < numDist; n++ )
  {
    const double nodeDist = (std::abs(dist[n]) < snapTol) ? 0. : dist[n];
    if (LevelSet::sign(nodeDist) < 0)
    {
      if (hasPos)
        return 0;
      hasNeg = true;
    }
    else
    {
      if (hasNeg)
        return 0;
      hasPos = true;
    }
  }
  STK_ThrowAssert(!(hasNeg && hasPos));
  return (hasNeg ? -1 : 1);
}

int ContourElement::compute_conservative_distance_sign(const double snapTol, const stk::topology distTopology, const double * dist)
{
  // Here "linear" means that the extrema within the element matches the extrema of the nodal values,
  // which is true for bilinear quads (and hex) as well as truly linear tets and tris.
  const unsigned numDist = distTopology.num_nodes();
  if (distTopology.base() == distTopology)
    return compute_linear_distance_sign(snapTol, numDist, dist);
  return compute_conservative_nonlinear_distance_sign(snapTol, numDist, dist);
}

void
ContourElement::compute_subelement_decomposition(const double in_length_scale, const double in_edge_linear_tolerance, const double in_edge_nonlinear_tolerance)
  { /* %TRACE% */  /* %TRACE% */
    stk::topology topology = dist_topology();

    const double snapTol = in_length_scale * in_edge_linear_tolerance;

    my_sign = compute_conservative_distance_sign(snapTol, topology, my_dist.ptr());
    if (my_sign != 0)
    {
      return;
    }

    my_length_scale = in_length_scale;
    my_edge_linear_tolerance = in_edge_linear_tolerance;
    my_edge_nonlinear_tolerance = in_edge_nonlinear_tolerance;

    // Now create our base subelement.  This is the subelement that has the same
    // topology as the distance field.  This involves recursive calls that will continue
    // to divide into subelements until conformal subelements are formed.

    if (stk::topology::HEXAHEDRON_8 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Hex_8>( get_nodal_parametric_coords_array<3,8>(my_dist_master_elem), get_side_ids_array<6>(), this );
    }
    else if (stk::topology::HEXAHEDRON_27 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Hex_27>( get_nodal_parametric_coords_array<3,27>(my_dist_master_elem), get_side_ids_array<6>(), this );
    }
    else if (stk::topology::TETRAHEDRON_4 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Tet_4>( get_nodal_parametric_coords_array<3,4>(my_dist_master_elem), get_side_ids_array<4>(), this );
    }
    else if (stk::topology::TETRAHEDRON_10 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Tet_10>( get_nodal_parametric_coords_array<3,10>(my_dist_master_elem), get_side_ids_array<4>(), this );
    }
    else if (stk::topology::WEDGE_6 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Wedge_6>( get_nodal_parametric_coords_array<3,6>(my_dist_master_elem), get_side_ids_array<5>(), this );
    }
    else if (stk::topology::QUADRILATERAL_4_2D == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Quad_4>( get_nodal_parametric_coords_array<2,4>(my_dist_master_elem), get_side_ids_array<4>(), this );
    }
    else if (stk::topology::QUADRILATERAL_9_2D == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Quad_9>( get_nodal_parametric_coords_array<2,9>(my_dist_master_elem), get_side_ids_array<4>(), this );
    }
    else if (stk::topology::TRIANGLE_3_2D == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Tri_3>( get_nodal_parametric_coords_array<2,3>(my_dist_master_elem), get_side_ids_array<3>(), this );
    }
    else if (stk::topology::TRIANGLE_6_2D == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Tri_6>( get_nodal_parametric_coords_array<2,6>(my_dist_master_elem), get_side_ids_array<3>(), this );
    }
    STK_ThrowErrorMsgIf(!my_base_subelement, "Element with topology " << topology.name() << " not supported.");

    if ( krinolog.shouldPrint(LOG_SUBELEMENT) )
      {
	if ( krinolog.shouldPrint(LOG_DEBUG) )
	  {
	    dump_subelement_details();
	  }
	else
	  {
	    dump_subelement_structure();
	  }
      }
  }

static double tet_volume(const sierra::ArrayContainer<double,DIM,NPE_COORD> & coords)
{
  const double * arrayData = coords.ptr();
  const std::array<stk::math::Vector3d,4> tetNodeCoords = {stk::math::Vector3d(arrayData), stk::math::Vector3d(arrayData+3), stk::math::Vector3d(arrayData+6), stk::math::Vector3d(arrayData+9)};
  return compute_tet_volume(tetNodeCoords);
}

static double tri_volume(const sierra::ArrayContainer<double,DIM,NPE_COORD> & coords)
{
  const double * arrayData = coords.ptr();
  const std::array<stk::math::Vector2d,3> triNodeCoords = {stk::math::Vector2d(arrayData), stk::math::Vector2d(arrayData+2), stk::math::Vector2d(arrayData+4)};
  return compute_tri_volume(triNodeCoords);
}

double
ContourElement::volume() const
  { /* %TRACE% */  /* %TRACE% */
    if (coord_topology() == stk::topology::TETRAHEDRON_4)
      return tet_volume(my_coords);
    if (coord_topology() == stk::topology::TRIANGLE_3_2D)
      return tri_volume(my_coords);

    sierra::Array<const double,DIM,NINT> intg_pt_locations;
    sierra::Array<const double,NINT> intg_weights;
    sierra::ArrayContainer<double,NINT> determinants;

    std_intg_pts( intg_pt_locations, intg_weights, determinants );

    double vol = 0.0;

    const int num_intg_pts = intg_pt_locations.dimension<NINT>();

    for ( int ip = 0; ip < num_intg_pts; ++ip )
      {
	vol += intg_weights(ip) * determinants(ip);
      }
    return vol;
  }

double
ContourElement::elem_size() const
  { /* %TRACE% */  /* %TRACE% */
    const double vol      = volume();
    const double vol_root = 1.0 / my_spatial_dim;
    const double h        = pow( vol, vol_root );
    return h;
  }

void
ContourElement::dump_subelement_structure() const
  { /* %TRACE% */  /* %TRACE% */
    if (my_sign == 0)
    {
      STK_ThrowErrorMsgIf(!my_base_subelement, "\ncompute_subelement_decomposition(...) must be called prior to calling dump_subelement_structure().");
      krinolog << "***********************************************\n";
      krinolog << *this;
      krinolog << "Subelement structure:\n";
      my_base_subelement->dump_structure();
      krinolog << "***********************************************" << stk::diag::dendl;
    }
  }

void
ContourElement::dump_subelement_details() const
  { /* %TRACE% */  /* %TRACE% */
    if (my_sign == 0)
    {
      STK_ThrowErrorMsgIf(!my_base_subelement, "\ncompute_subelement_decomposition(...) must be called prior to calling dump_subelement_details().");
      krinolog << "***********************************************\n";
      krinolog << *this;
      krinolog << "Subelement details:\n";
      my_base_subelement->dump_details();
      krinolog << "***********************************************" << stk::diag::dendl;
    }
  }

void
ContourElement::build_subelement_facets( FacetedSurfaceBase & facets )
{
  if (my_sign == 0)
  {
    STK_ThrowErrorMsgIf(!my_base_subelement, "\ncompute_subelement_decomposition(...) must be called prior to calling build_subelement_facets().");
    my_base_subelement->build_facets( facets );
  }
}

int
ContourElement::gather_intg_pts( const int intg_pt_sign,
			  sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
			  sierra::ArrayContainer<double,NINT> & intg_weights,
			  sierra::ArrayContainer<double,NINT> & determinants ) const
{
  if (my_sign != 0)
  {
    if (my_sign == intg_pt_sign)
    {
      return std_intg_pts(intg_pt_locations, intg_weights, determinants);
    }

    intg_pt_locations.resize(my_spatial_dim,0);
    intg_weights.resize(0u);
    determinants.resize(0u);
    return 0;
  }
  STK_ThrowErrorMsgIf(!my_base_subelement, "\ncompute_subelement_decomposition(...) must be called prior to calling gather_intg_pts().");

  const int num_intg_pts = my_base_subelement->num_intg_pts(intg_pt_sign);

  intg_pt_locations.resize(my_spatial_dim,num_intg_pts);
  intg_weights.resize(num_intg_pts);
  determinants.resize(num_intg_pts);

  my_base_subelement->gather_intg_pts( intg_pt_sign,
				       intg_pt_locations,
				       intg_weights,
				       determinants );

  if ( 0 != intg_pt_sign )
    {
      //
      // Include the determinant from element to real space
      //
      const int num_coord_dofs = coord_topology().num_nodes();

      // temporary data
      static sierra::ArrayContainer<double,DIM,NPE_COORD,NINT> d_shapef_coords;
      static sierra::ArrayContainer<double,NINT>               det_J;

      // reserve enough space for results
      d_shapef_coords.resize(my_spatial_dim, num_coord_dofs, num_intg_pts);
      det_J.resize(num_intg_pts);

      double det_J_error;

      my_coords_master_elem.shape_fcn_deriv(num_intg_pts, intg_pt_locations.ptr(), d_shapef_coords.ptr());
      my_coords_master_elem.determinant(
          my_spatial_dim,           // Number of coordinate dimensions
          num_intg_pts,             // Number of target points
          num_coord_dofs,           // Number of coord shape functions
          d_shapef_coords.ptr(),// Mesh shape function derivatives
          1,                        // Number of elements
          my_coords.ptr(),  // Mesh coordinate values
          det_J.ptr(),              // Determinant of the transformation Jacobian for each element (output)
          &det_J_error );           // Determinant error (output)
      for ( int i=0; i<num_intg_pts; ++i )
	{
	  determinants(i) *= det_J(i);
	}
    }

  return( num_intg_pts );
}

int
ContourElement::std_intg_pts( sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
    sierra::ArrayContainer<double,NINT> & intg_weights,
    sierra::ArrayContainer<double,NINT> & determinants,
    const MasterElement & me ) const
{
  sierra::Array<const double,DIM,NINT> std_intg_pt_locations;
  sierra::Array<const double,NINT> std_intg_weights;
  std_intg_pts(std_intg_pt_locations, std_intg_weights, determinants, me);

  const int num_intg_pts = determinants.dimension<0>();

  intg_pt_locations.resize(my_spatial_dim,num_intg_pts);
  intg_weights.resize(num_intg_pts);

  const double * stdIntgPtLocData = std_intg_pt_locations.ptr();
  double * intgPtLocData = intg_pt_locations.ptr();
  for ( int i=0; i<num_intg_pts*my_spatial_dim; ++i )
    intgPtLocData[i] = stdIntgPtLocData[i];
  for ( int i=0; i<num_intg_pts; ++i )
    intg_weights(i) = std_intg_weights(i);

  return num_intg_pts;
}

int
ContourElement::std_intg_pts( const MasterElement & evalMasterElem,
  const MasterElement & coordMasterElem,
  const sierra::Array<double,DIM,NPE_COORD> & coords,
  sierra::Array<const double,DIM,NINT> & intg_pt_locations,
  sierra::Array<const double,NINT> & intg_weights,
  sierra::ArrayContainer<double,NINT> & det_J)
{
  const int dim = coords.dimension(0);
  const int num_intg_pts         = evalMasterElem.num_intg_pts();
  const double * intg_pt_loc_ptr = evalMasterElem.intg_pt_locations();
  const double * intg_wt_ptr     = evalMasterElem.intg_weights();

  intg_pt_locations.set(intg_pt_loc_ptr,dim,num_intg_pts);
  intg_weights.set(intg_wt_ptr,num_intg_pts);
  det_J.resize( num_intg_pts );

  double det_J_error;

  if ( evalMasterElem.get_topology() == coordMasterElem.get_topology() )
  {
    coordMasterElem.determinant( dim, 1, coords.ptr(), det_J.ptr(), &det_J_error );
  }
  else
  {
    const int num_coord_dofs = coordMasterElem.get_topology().num_nodes();
    sierra::ArrayContainer<double,DIM,NPE_COORD,NINT> d_shapef_coords(dim, num_coord_dofs, num_intg_pts);
    coordMasterElem.shape_fcn_deriv(num_intg_pts, intg_pt_locations.ptr(), d_shapef_coords.ptr());
    coordMasterElem.determinant(
        dim,                    // Number of coordinate dimensions
        num_intg_pts,           // Number of target points
        num_coord_dofs,         // Number of coord shape functions
        d_shapef_coords.ptr(),  // Mesh shape function derivatives
        1,                      // Number of elements
        coords.ptr(),           // Mesh coordinate values
        det_J.ptr(),            // Determinant of the transformation Jacobian for each element (output)
        &det_J_error );         // Determinant error (output)
  }

  return( num_intg_pts );
}

int
ContourElement::std_intg_pts( sierra::Array<const double,DIM,NINT> & intg_pt_locations,
		       sierra::Array<const double,NINT> & intg_weights,
		       sierra::ArrayContainer<double,NINT> & det_J,
		       const MasterElement & me ) const
{
  return std_intg_pts(me, my_coords_master_elem, my_coords, intg_pt_locations, intg_weights, det_J);
}

int
ContourElement::std_side_intg_pts( const int iside,
  sierra::Array<const double,DIM,NINT> & intg_pt_locations,
  sierra::Array<const double,NINT> & intg_weights,
  sierra::ArrayContainer<double,NINT> & det_J,
  const MasterElement & me ) const
{
  const stk::topology coordTopo = my_coords_master_elem.get_topology();
  const stk::topology coordSideTopo = coordTopo.side_topology(iside);
  const auto & coordSideMasterElement = MasterElementDeterminer::getMasterElement(coordSideTopo);
  const int numSideNodes = coordSideTopo.num_nodes();
  const unsigned * const lnn = get_side_node_ordinals(coordTopo, iside );

  sierra::ArrayContainer<double,DIM,NPE_COORD> sideNodeCoords(my_spatial_dim, numSideNodes);  // temp array
  for ( int i = 0; i < numSideNodes; i++ )
    for ( int d = 0; d < my_spatial_dim; d++ )
      sideNodeCoords(d,i) = my_coords(d, lnn[i]);

  const auto & evalSideMasterElement = (me.get_topology() == my_coords_master_elem.get_topology()) ?
      coordSideMasterElement :
      MasterElementDeterminer::getMasterElement(me.get_topology().side_topology(iside));

  return std_intg_pts(evalSideMasterElement, coordSideMasterElement, sideNodeCoords, intg_pt_locations, intg_weights, det_J);
}

} // namespace krino
