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
#include <stk_math/StkVector.hpp>

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

  std::string debug_output() const;
  friend std::ostream & operator << ( std::ostream &os , const ContourElement &s ) {
    os << s.debug_output();
    return os;
  }

  const MasterElement & coord_master_elem() const { return my_coords_master_elem; }
  const MasterElement & dist_master_elem() const { return my_dist_master_elem; }

  const stk::topology coord_topology() const { return my_coords_master_elem.get_topology(); }
  const stk::topology dist_topology() const { return my_dist_master_elem.get_topology(); }

  static double compute_domain_integral(const sierra::ArrayContainer<double, DIM, NINT> & intgPtLocations,
    const sierra::ArrayContainer<double, NINT> & intgWeights,
    const sierra::ArrayContainer<double, NINT> &  determinants);

  double compute_domain_integral(const int signOfDomain) const; // 0=interface, -1=negVol, 1=posVol

  double compute_area_of_interface() const;

  double compute_signed_volume(const int signOfVolume) const;

  static int compute_conservative_nonlinear_distance_sign(const double snapTol, const unsigned numDist, const double * dist);
  static int compute_linear_distance_sign(const double snapTol, const unsigned numDist, const double * dist);
  static int compute_conservative_distance_sign(const double snapTol, const stk::topology distTopology, const double * dist);

  bool dist_is_linear() const { return dist_topology() == stk::topology::TRIANGLE_3_2D || dist_topology() == stk::topology::TETRAHEDRON_4; }

  stk::math::Vector3d coordinates( const stk::math::Vector3d & p_coords ) const;
  double distance( const stk::math::Vector3d & p_coords ) const;
  double determinant( const stk::math::Vector3d & p_coords ) const;

  stk::math::Vector3d distance_gradient( const stk::math::Vector3d & p_coords ) const;
  void compute_distance_gradient( const sierra::Array<const double,DIM,NINT> & intg_pt_locations,
				  sierra::ArrayContainer<double,DIM,NINT> & grad_dist ) const;

  void dump_subelement_structure( void ) const;
  void dump_subelement_details( void ) const;

  void build_subelement_facets( FacetedSurfaceBase & facets );

  int gather_intg_pts( const int intg_pt_sign,
		       sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
		       sierra::ArrayContainer<double,NINT> & intg_weights,
		       sierra::ArrayContainer<double,NINT> & determinants ) const;

  static int std_intg_pts( const MasterElement & evalMasterElem,
      const MasterElement & coordMasterElem,
      const sierra::Array<double,DIM,NPE_COORD> & coords,
      sierra::Array<const double,DIM,NINT> & intg_pt_locations,
      sierra::Array<const double,NINT> & intg_weights,
      sierra::ArrayContainer<double,NINT> & det_J);
  int std_intg_pts( sierra::Array<const double,DIM,NINT> & intg_pt_locations,
		    sierra::Array<const double,NINT> & intg_weights,
		    sierra::ArrayContainer<double,NINT> & determinants,
		    const MasterElement & me ) const;
  int std_intg_pts( sierra::Array<const double,DIM,NINT> & intg_pt_locations,
                    sierra::Array<const double,NINT> & intg_weights,
                    sierra::ArrayContainer<double,NINT> & determinants ) const { return std_intg_pts(intg_pt_locations, intg_weights, determinants, my_coords_master_elem); }
  int std_intg_pts( sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
                    sierra::ArrayContainer<double,NINT> & intg_weights,
                    sierra::ArrayContainer<double,NINT> & determinants,
                    const MasterElement & me ) const;
  int std_intg_pts( sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
                    sierra::ArrayContainer<double,NINT> & intg_weights,
                    sierra::ArrayContainer<double,NINT> & determinants ) const { return std_intg_pts(intg_pt_locations, intg_weights, determinants, my_coords_master_elem); }

  int std_side_intg_pts( const int iside,
      sierra::Array<const double,DIM,NINT> & intg_pt_locations,
      sierra::Array<const double,NINT> & intg_weights,
      sierra::ArrayContainer<double,NINT> & det_J,
      const MasterElement & me ) const;
  int std_side_intg_pts( const int iside,
      sierra::Array<const double,DIM,NINT> & intg_pt_locations,
      sierra::Array<const double,NINT> & intg_weights,
      sierra::ArrayContainer<double,NINT> & det_J) const  { return std_side_intg_pts(iside, intg_pt_locations, intg_weights, det_J, my_coords_master_elem); }

  int spatial_dim() const { return my_spatial_dim; }
  double length_scale() const { return my_length_scale; }
  double edge_linear_tolerance() const { return my_edge_linear_tolerance; }
  double edge_nonlinear_tolerance() const { return my_edge_nonlinear_tolerance; }

  double volume() const;
  double elem_size() const;

  void compute_subelement_decomposition(const double length_scale, const double edge_linear_tolerance = 1.e-4, const double edge_nonlinear_tolerance = 1.0e-2);

private:
  const int my_spatial_dim;

  const MasterElement & my_coords_master_elem;
  const MasterElement & my_dist_master_elem;

  sierra::ArrayContainer<double,DIM,NPE_COORD> my_coords;
  sierra::ArrayContainer<double,NPE_VAR> my_dist;

  int my_sign;
  double my_length_scale;
  double my_edge_linear_tolerance;
  double my_edge_nonlinear_tolerance;
  std::unique_ptr<ContourSubElement> my_base_subelement;

  //: Default constructor not allowed
  ContourElement();

};

} // namespace krino

#endif // Akri_ContourElement_h
