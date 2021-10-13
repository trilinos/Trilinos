// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_ContourElement_h
#define Akri_ContourElement_h

#include <Akri_TypeDefs.hpp>
#include <vector>
#include <set>
#include <map>

#include <Akri_DiagWriter.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_Vec.hpp>

namespace krino {

class ContourSubElement;

class ContourElement {
public:

  ContourElement( const stk::mesh::BulkData & mesh,
	   stk::mesh::Entity mesh_obj,
	   const FieldRef coords_field,
	   const FieldRef dist_field,
	   const double iso_dist = 0.0 );
  ~ContourElement(); // Definition must in implementation file because SubElement is incomplete

  std::ostream & put( std::ostream& os ) const;
  friend std::ostream & operator << ( std::ostream &os , const ContourElement &s ) {
    return s.put(os);
  }

  const MasterElement & coord_master_elem() const { return my_coords_master_elem; }
  const MasterElement & dist_master_elem() const { return my_dist_master_elem; }
  const MasterElement & vel_master_elem() const { return *my_vel_master_elem; }

  const stk::topology coord_topology() const { return my_coords_master_elem.get_topology(); }
  const stk::topology dist_topology() const { return my_dist_master_elem.get_topology(); }
  const stk::topology vel_topology() const { return my_vel_master_elem->get_topology(); }

  bool dist_is_linear() const { return dist_topology() == stk::topology::TRIANGLE_3_2D || dist_topology() == stk::topology::TETRAHEDRON_4; }

  stk::mesh::Entity entity() const { return my_entity; }
  stk::mesh::EntityId elem_id() const { return my_mesh.identifier(my_entity); }
  bool have_interface_sides() const;

  Vector3d coordinates( const Vector3d & p_coords ) const;
  double distance( const Vector3d & p_coords ) const;
  const PointVec & nodal_parametric_coordinates() const { return my_dist_p_coords; }
  double determinant( const Vector3d & p_coords ) const;

  Vector3d continuous_distance_gradient( const Vector3d & p_coords ) const;
  Vector3d distance_gradient( const Vector3d & p_coords ) const;
  void compute_distance_gradient( const sierra::Array<const double,DIM,NINT> & intg_pt_locations,
				  sierra::ArrayContainer<double,DIM,NINT> & grad_dist ) const;

  void dump_subelement_structure( void ) const;
  void dump_subelement_details( void ) const;
  const ContourSubElement * get_base_subelement() const { return my_base_subelement.get(); }

  int build_subelement_facets( Faceted_Surface & facets );

  int gather_intg_pts( const int intg_pt_sign,
		       sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
		       sierra::ArrayContainer<double,NINT> & intg_weights,
		       sierra::ArrayContainer<double,NINT> & determinants,
		       const bool map_to_real_coords );

  int std_intg_pts( sierra::Array<const double,DIM,NINT> & intg_pt_locations,
		    sierra::Array<const double,NINT> & intg_weights,
		    sierra::ArrayContainer<double,NINT> & determinants,
		    const MasterElement & me ) const;
  int std_intg_pts( sierra::Array<const double,DIM,NINT> & intg_pt_locations,
                    sierra::Array<const double,NINT> & intg_weights,
                    sierra::ArrayContainer<double,NINT> & determinants ) const { return std_intg_pts(intg_pt_locations, intg_weights, determinants, my_coords_master_elem); }

  int spatial_dim() const { return my_spatial_dim; }
  double length_scale() const { return my_length_scale; }
  double edge_linear_tolerance() const { return my_edge_linear_tolerance; }
  double edge_nonlinear_tolerance() const { return my_edge_nonlinear_tolerance; }

  double volume() const;
  double elem_size() const;
  double average_edge_length() const;

  void compute_subelement_decomposition(const double length_scale, const double edge_linear_tolerance = 1.e-4, const double edge_nonlinear_tolerance = 1.0e-2) const;

private:
  const stk::mesh::BulkData & my_mesh;
  stk::mesh::Entity my_entity;
  const int my_spatial_dim;

  const MasterElement & my_coords_master_elem;
  const MasterElement & my_dist_master_elem;
  const MasterElement * my_vel_master_elem;

  const FieldRef my_coords_field;
  const FieldRef my_dist_field;
  const FieldRef my_vel_field;

  sierra::ArrayContainer<double,DIM,NPE_COORD> my_coords;
  sierra::ArrayContainer<double,NPE_COORD,DIM> my_coords_transpose;
  sierra::ArrayContainer<double,NPE_VAR> my_dist;
  sierra::ArrayContainer<double,DIM,NPE_VAR> my_vel;
  sierra::ArrayContainer<double,NPE_VAR,DIM> my_vel_transpose;

  mutable PointVec my_dist_p_coords;
  mutable double my_length_scale;
  mutable double my_edge_linear_tolerance;
  mutable double my_edge_nonlinear_tolerance;
  mutable std::unique_ptr<ContourSubElement> my_base_subelement;

  //: Default constructor not allowed
  ContourElement();

};

} // namespace krino

#endif // Akri_ContourElement_h
