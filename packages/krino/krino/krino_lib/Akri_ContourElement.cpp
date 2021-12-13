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

namespace krino{

ContourElement::ContourElement( const stk::mesh::BulkData & mesh,
		  stk::mesh::Entity in_mesh_obj,
		  const FieldRef coords_field,
                  const FieldRef dist_field,
		  const double iso_dist )
    : my_mesh( mesh ),
    my_entity( in_mesh_obj ),
    my_spatial_dim( mesh.mesh_meta_data().spatial_dimension() ),
    my_coords_master_elem(MasterElementDeterminer::getMasterElement(mesh.bucket(in_mesh_obj), coords_field)),
    my_dist_master_elem(MasterElementDeterminer::getMasterElement(mesh.bucket(in_mesh_obj), dist_field)),
    my_coords_field( coords_field ),
    my_dist_field( dist_field ),
    my_length_scale(-1.0),
    my_edge_linear_tolerance(-1.0),
    my_edge_nonlinear_tolerance(-1.0)
{ /* %TRACE% */  /* %TRACE% */

  //
  // Gather obj's coordinates, distance, and velocity into ArrayContainer's.
  //

  const int dim = my_spatial_dim;

  const int npe_coords = my_coords_master_elem.get_topology().num_nodes();
  const int npe_dist   = my_dist_master_elem.get_topology().num_nodes();

  my_coords.resize( dim, npe_coords );
  my_dist.resize( npe_dist );

  const stk::mesh::Entity* elem_nodes = mesh.begin_nodes(my_entity);

  for ( int i = 0; i < npe_coords; ++i )
  {
    double * var = field_data<double>( my_coords_field, elem_nodes[i]);
    for ( int d = 0; d < dim; d++ )
    {
      my_coords(d,i) = var[d];
    }
  }

  for ( int i = 0; i < npe_dist; ++i )
  {
    double * var = field_data<double>( my_dist_field, elem_nodes[i]);
    my_dist(i) = *var - iso_dist;
  }

  // the interpolate master element methods require the transpose?!?
  my_coords_transpose.resize( npe_coords, dim );
  for (int i = 0; i < npe_coords; i++)
  {
    for (int d = 0; d < dim; d++)
    {
      my_coords_transpose(i,d) = my_coords(d,i);
    }
  }
}

ContourElement::~ContourElement() {}

std::ostream &
ContourElement::put( std::ostream& os ) const
{ /* %TRACE% */  /* %TRACE% */
  const int dim = my_spatial_dim;

  os << "Element description:" << std::endl
     << "  Coordinates topolology is " << coord_topology().name()
     << ", Distance topology is " << dist_topology().name() << std::endl;

  int lnn = 0;
  int coord_counter = 0;
  int dist_counter = 0;
  const unsigned num_nodes = my_mesh.num_nodes(my_entity);
  const stk::mesh::Entity* nodes = my_mesh.begin_nodes(my_entity);
  for (unsigned node_index=0; node_index<num_nodes; ++node_index)
  {
    stk::mesh::Entity node = nodes[node_index];
    os << "  Node #" << lnn++ << ", Global Node #" << my_mesh.identifier(node);
    double * var = field_data<double>(my_coords_field, node);
    if ( NULL != var )
    {
      Vector3d coords(Vector3d::ZERO);
      for ( int d = 0; d < dim; d++ )
        coords[d] = my_coords(d,coord_counter);
      os << ", coords = (" << coords[0] << ","
         << coords[1] << ","
         << coords[2] << ")";
      coord_counter++;
    }

    var = field_data<double>(my_dist_field, node);
    if ( NULL != var )
    {
      os << ", dist = " << my_dist(dist_counter++);
    }
    os << std::endl;
  }
  return os ;
}

bool
ContourElement::have_interface_sides() const
{
  ThrowAssert(my_base_subelement);
  return my_base_subelement->have_interface_sides();
}

double
ContourElement::distance( const Vector3d & p_coords ) const
  { /* %TRACE% */  /* %TRACE% */
    double dist;
    my_dist_master_elem.interpolate_point(my_spatial_dim, p_coords.data(), 1, my_dist.ptr(), &dist);
    return dist;
  }

Vector3d
ContourElement::coordinates( const Vector3d & p_coords ) const
  { /* %TRACE% */  /* %TRACE% */
    Vector3d coords(Vector3d::ZERO);
    my_coords_master_elem.interpolate_point(my_spatial_dim, p_coords.data(), my_spatial_dim, my_coords.ptr(), coords.data());
    return coords;
  }

double
ContourElement::determinant( const Vector3d & p_coords ) const
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

Vector3d
ContourElement::continuous_distance_gradient( const Vector3d & p_coords ) const
  { /* %TRACE% */  /* %TRACE% */
    // gather continuous distance gradient
    sierra::ArrayContainer<double,DIM,NPE_VAR> grad_distance;
    const int npe = dist_topology().num_nodes();
    grad_distance.resize( my_spatial_dim, npe );
    stk::mesh::FieldBase* field_ptr = my_mesh.mesh_meta_data().get_field(stk::topology::NODE_RANK, "CONT_GRAD");
    ThrowRequireMsg(nullptr != field_ptr, "Field CONT_GRAD not found.");
    const FieldRef contGradField(field_ptr);
    const stk::mesh::Entity* elem_nodes = my_mesh.begin_nodes(my_entity);
    for ( int i = 0; i < npe; ++i )
    {
      double * var = field_data<double>(contGradField, elem_nodes[i]);
      for(int d = 0; d < my_spatial_dim; ++d)
        grad_distance(d,i) = var[d];
    }

    Vector3d grad_dist(Vector3d::ZERO);
    my_vel_master_elem->interpolate_point(my_spatial_dim, p_coords.data(), my_spatial_dim, grad_distance.ptr(), grad_dist.data());
    return grad_dist;
  }

Vector3d
ContourElement::distance_gradient( const Vector3d & p_coords ) const
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
    Vector3d grad_dist( Vector3d::ZERO );
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

void
ContourElement::compute_subelement_decomposition(const double in_length_scale, const double in_edge_linear_tolerance, const double in_edge_nonlinear_tolerance) const
  { /* %TRACE% */  /* %TRACE% */
    my_length_scale = in_length_scale;
    my_edge_linear_tolerance = in_edge_linear_tolerance;
    my_edge_nonlinear_tolerance = in_edge_nonlinear_tolerance;

    stk::topology topology = dist_topology();
    const int num_sides = topology.num_sides();
    std::vector<int> side_ids(num_sides);
    for (int i=0; i<num_sides; ++i)
      side_ids[i] = i;

    // Generate my_dist_p_coords, a vector of Vector3d's of the parametric coordinates
    // at each node where the distance is defined.  We need this because this is the topology
    // that our base subelement will have.
    const int npe_dist = dist_topology().num_nodes();
    my_dist_p_coords.reserve( npe_dist );
    const double * p_coords_array = my_dist_master_elem.nodal_parametric_coordinates();
    int p_coords_array_counter = 0;

    for ( int i = 0; i < npe_dist; ++i )
      {
	Vector3d p_coords(Vector3d::ZERO);
	for ( int j = 0; j < my_spatial_dim; j++ )
	  p_coords[j] = p_coords_array[p_coords_array_counter++];
	my_dist_p_coords.push_back( p_coords );
      }

    // Now create our base subelement.  This is the subelement that has the same
    // topology as the distance field.  This involves recursive calls that will continue
    // to divide into subelements until conformal subelements are formed.

    if (stk::topology::HEXAHEDRON_8 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Hex_8>( my_dist_p_coords, side_ids, this );
    }
    else if (stk::topology::HEXAHEDRON_27 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Hex_27>( my_dist_p_coords, side_ids, this );
    }
    else if (stk::topology::TETRAHEDRON_4 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Tet_4>( my_dist_p_coords, side_ids, this );
    }
    else if (stk::topology::TETRAHEDRON_10 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Tet_10>( my_dist_p_coords, side_ids, this );
    }
    else if (stk::topology::WEDGE_6 == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Wedge_6>( my_dist_p_coords, side_ids, this );
    }
    else if (stk::topology::QUADRILATERAL_4_2D == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Quad_4>( my_dist_p_coords, side_ids, this );
    }
    else if (stk::topology::QUADRILATERAL_9_2D == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Quad_9>( my_dist_p_coords, side_ids, this );
    }
    else if (stk::topology::TRIANGLE_3_2D == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Tri_3>( my_dist_p_coords, side_ids, this );
    }
    else if (stk::topology::TRIANGLE_6_2D == topology)
    {
      my_base_subelement = std::make_unique<ContourSubElement_Tri_6>( my_dist_p_coords, side_ids, this );
    }
    ThrowErrorMsgIf(!my_base_subelement, "Element with topology " << topology.name() << " not supported.");

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

double
ContourElement::volume() const
  { /* %TRACE% */  /* %TRACE% */

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
ContourElement::average_edge_length() const
{ /* %TRACE% */  /* %TRACE% */
  const stk::topology Top = coord_topology();
  int num_edges = Top.num_edges();

  double sum_edge_lengths = 0.0;

  for ( int edge = 0; edge < num_edges; edge++ )
    {
      const unsigned * const lnn = get_edge_node_ordinals(Top, edge);

      double sqr_length = 0.0;
      for ( int d = 0; d < my_spatial_dim; d++ ) sqr_length += (my_coords(d,lnn[0]) - my_coords(d,lnn[1])) *
							       (my_coords(d,lnn[0]) - my_coords(d,lnn[1]));
      sum_edge_lengths += std::sqrt(sqr_length);
    }

  return sum_edge_lengths/num_edges;
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
    ThrowErrorMsgIf(!my_base_subelement, "\ncompute_subelement_decomposition(...) must be called prior to calling dump_subelement_structure().");

    krinolog << "***********************************************" << stk::diag::dendl;
    krinolog << *this;
    krinolog << "Subelement structure:" << stk::diag::dendl;
    my_base_subelement->dump_structure();
    krinolog << "***********************************************" << stk::diag::dendl;
  }

void
ContourElement::dump_subelement_details() const
  { /* %TRACE% */  /* %TRACE% */
    ThrowErrorMsgIf(!my_base_subelement, "\ncompute_subelement_decomposition(...) must be called prior to calling dump_subelement_details().");

    krinolog << "***********************************************" << stk::diag::dendl;
    krinolog << *this;
    krinolog << "Subelement details:" << stk::diag::dendl;
    my_base_subelement->dump_details();
    krinolog << "***********************************************" << stk::diag::dendl;
  }

int
ContourElement::build_subelement_facets( Faceted_Surface & facets )
{
  ThrowErrorMsgIf(!my_base_subelement, "\ncompute_subelement_decomposition(...) must be called prior to calling build_subelement_facets(...).");
  return my_base_subelement->build_facets( facets );
}

int
ContourElement::gather_intg_pts( const int intg_pt_sign,
			  sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
			  sierra::ArrayContainer<double,NINT> & intg_weights,
			  sierra::ArrayContainer<double,NINT> & determinants,
			  const bool map_to_real_coords )
{
  ThrowErrorMsgIf(0 == intg_pt_sign && !map_to_real_coords,
      "\nSubelement decomposition can currently only provide the surface integration with the overall determinant.");
  ThrowErrorMsgIf(!my_base_subelement,
      "\ncompute_subelement_decomposition(...) must be called prior to calling gather_intg_pts(...).");

  const int num_intg_pts = my_base_subelement->num_intg_pts(intg_pt_sign);

  intg_pt_locations.resize(my_spatial_dim,num_intg_pts);
  intg_weights.resize(num_intg_pts);
  determinants.resize(num_intg_pts);

  my_base_subelement->gather_intg_pts( intg_pt_sign,
				       intg_pt_locations,
				       intg_weights,
				       determinants );

  if ( 0 != intg_pt_sign && map_to_real_coords )
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
ContourElement::std_intg_pts( sierra::Array<const double,DIM,NINT> & intg_pt_locations,
		       sierra::Array<const double,NINT> & intg_weights,
		       sierra::ArrayContainer<double,NINT> & det_J,
		       const MasterElement & me ) const
  {
    const int num_intg_pts       = me.num_intg_pts();
    const double * intg_pt_loc_ptr = me.intg_pt_locations();
    const double * intg_wt_ptr     = me.intg_weights();

    intg_pt_locations.set(intg_pt_loc_ptr,my_spatial_dim,num_intg_pts);
    intg_weights.set(intg_wt_ptr,num_intg_pts);
    det_J.resize( num_intg_pts );

    double det_J_error;

    if ( me.get_topology() == my_coords_master_elem.get_topology() )
      {
	my_coords_master_elem.determinant( my_spatial_dim, 1, my_coords.ptr(), det_J.ptr(), &det_J_error );
      }
    else
      {
	const int num_coord_dofs = coord_topology().num_nodes();
	sierra::ArrayContainer<double,DIM,NPE_COORD,NINT> d_shapef_coords(my_spatial_dim, num_coord_dofs, num_intg_pts);
	my_coords_master_elem.shape_fcn_deriv(num_intg_pts, intg_pt_locations.ptr(), d_shapef_coords.ptr());
	my_coords_master_elem.determinant(
	    my_spatial_dim,         // Number of coordinate dimensions
	    num_intg_pts,           // Number of target points
            num_coord_dofs,         // Number of coord shape functions
            d_shapef_coords.ptr(),// Mesh shape function derivatives
            1,                      // Number of elements
            my_coords.ptr(),        // Mesh coordinate values
            det_J.ptr(),            // Determinant of the transformation Jacobian for each element (output)
            &det_J_error );         // Determinant error (output)
      }

    return( num_intg_pts );
  }

} // namespace krino
