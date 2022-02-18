// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_ContourSubElement_h
#define Akri_ContourSubElement_h

#include <Akri_ContourElement.hpp>
#include <vector>
#include <set>
#include <map>
#include <memory>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Akri_Faceted_Surface.hpp>
#include <Akri_Facet.hpp>
#include <Akri_Vec.hpp>

namespace krino {

template<class T>
bool compare_ptr_by_global_id (const T* i, const T* j) { return (i->global_id() < j->global_id()); }

class ContourSubElement {
public:

  ContourSubElement( const stk::topology topo,
	      const PointVec & coords,
	      const std::vector<int> & side_ids,
	      const ContourElement * owner,
              const int in_subelement_depth,
              const int subelement_sign );

  virtual ~ContourSubElement();

  static bool is_more(const Vector3d & v1, const Vector3d & v2);

  double relative_volume() const;
  double side_relative_area( const int side ) const;
  double parametric_quality() const;
  double physical_quality() const;
  double side_quality(const int side) const;

  int num_intg_pts( const int intg_pt_sign );

  int gather_intg_pts( const int intg_pt_sign,
		       sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
		       sierra::ArrayContainer<double,NINT> & intg_weights,
		       sierra::ArrayContainer<double,NINT> & determinants,
		       int index = 0 );

  int build_facets( Faceted_Surface & facets );

  // default implementation
  virtual int side_facets( Faceted_Surface & facets, int side ) const;

  stk::topology topology() const { return my_master_element.get_topology(); }

  int spatial_dim() const { return my_owner->spatial_dim(); }
  int num_subelements() const { return my_subelements.size(); }

  std::vector<int> & get_side_ids() { return my_side_ids; }
  const std::vector<int> & get_side_ids() const { return my_side_ids; }

  const PointVec & get_coords() const { return my_coords; }
  const ContourElement * owner() const { return my_owner; }
  int subelement_depth() const { return my_subelement_depth; }

  void dump_structure() const;
  void dump_details() const;

  bool have_interface_sides() const;

  virtual std::ostream & put( std::ostream& os ) const;

  friend std::ostream & operator << ( std::ostream &os , const ContourSubElement &s ) {
    return s.put(os);
  }

  static double find_quadratic_crossing( const double d0,
                                       const double d1,
                                       const double d2 );

protected:

  const MasterElement& my_master_element;
  const MasterElement& my_side_master_element;
  int my_num_nodes;
  int my_num_sides;
  PointVec my_coords;
  std::vector<int> my_side_ids;
  std::vector<double> my_dist;
  std::vector< ContourSubElement * > my_subelements;
  const ContourElement * my_owner;
  int my_subelement_depth; // depth down the tree of subelements (0 for the base_subelement)
  int my_sign; // -1 for elements on negative side, +1 for positive side, 0 for non-conformal elements spanning interface

private:
  //: Default constructor not allowed
  ContourSubElement();
};

class ContourSubElement_Quad_4 : public ContourSubElement {
public:
  ContourSubElement_Quad_4( const PointVec & coords,
		     const std::vector<int> & side_ids,
		     const ContourElement * in_owner );
  virtual ~ContourSubElement_Quad_4() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Quad_4();

  int non_conformal_decomposition();
};

class ContourSubElement_Quad_9 : public ContourSubElement {
public:
  ContourSubElement_Quad_9( const PointVec & coords,
		     const std::vector<int> & side_ids,
		     const ContourElement * in_owner );

  virtual ~ContourSubElement_Quad_9() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Quad_9();

  int non_conformal_decomposition();
};

class ContourSubElement_Hex_8 : public ContourSubElement {
public:
  ContourSubElement_Hex_8( const PointVec & coords,
		    const std::vector<int> & side_ids,
		    const ContourElement * in_owner );
  virtual ~ContourSubElement_Hex_8() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Hex_8();

  int subpyramid_non_conformal_decomposition( const int face );
};

class ContourSubElement_Hex_27 : public ContourSubElement {
public:
  ContourSubElement_Hex_27( const PointVec & coords,
		     const std::vector<int> & side_ids,
		     const ContourElement * in_owner );
  virtual ~ContourSubElement_Hex_27() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Hex_27();

  int subpyramid_non_conformal_decomposition( const int face );
};

class ContourSubElement_Wedge_6 : public ContourSubElement {
public:
  ContourSubElement_Wedge_6( const PointVec & coords,
                    const std::vector<int> & side_ids,
                    const ContourElement * in_owner );
  virtual ~ContourSubElement_Wedge_6() {}

private:
  int subpyramid_non_conformal_decomposition( const int face );
};

class ContourSubElement_Tri_3 : public ContourSubElement {
public:
  ContourSubElement_Tri_3( const PointVec & coords,
		    const std::vector<int> & side_ids,
		    const ContourElement * in_owner,
                    const int in_subelement_depth = 0,
                    const int subelement_sign = 0 );
  virtual ~ContourSubElement_Tri_3() {}

  virtual int side_facets( Faceted_Surface & facets, int side ) const;

private:
  //: Default constructor not allowed
  ContourSubElement_Tri_3();

  int conformal_decomposition();

  int process_edge( const int i0,
                    const int i1,
                    const int i2,
                    std::vector<int> & is_on_surf,
                    PointVec & lnodes,
                    const std::vector<double> & ldist );

  bool is_degenerate( const std::vector<int> & edge_node_ids,
                      const int i0, const int i1, const int i2 );
};

class ContourSubElement_Adaptive_Tri_3 : public ContourSubElement {
public:
  ContourSubElement_Adaptive_Tri_3( const PointVec & coords,
                             const std::vector<int> & side_ids,
                             const ContourElement * in_owner,
                             const int in_subelement_depth = 0);
  ContourSubElement_Adaptive_Tri_3( const PointVec & coords,
                             const std::vector<int> & side_ids,
                             const std::vector<int> & edge_age,
                             const ContourElement * in_owner,
                             const int in_subelement_depth = 0);
  virtual ~ContourSubElement_Adaptive_Tri_3() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Adaptive_Tri_3();

  static const int MAX_REFINMENT_LEVELS;

  int non_conformal_decomposition();

  std::vector<int> my_edge_age;
};

class ContourSubElement_Tri_6 : public ContourSubElement {
public:
  ContourSubElement_Tri_6( const PointVec & coords,
                    const std::vector<int> & side_ids,
                    const ContourElement * in_owner,
                    const int in_subelement_depth = 0,
                    const int subelement_sign = 0 );
  virtual ~ContourSubElement_Tri_6() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Tri_6();
};

class ContourSubElement_Tet_4 : public ContourSubElement {
public:
  ContourSubElement_Tet_4( const PointVec & coords,
                    const std::vector<int> & side_ids,
                    const ContourElement * in_owner,
                    const int in_subelement_depth = 0,
                    const int subelement_sign = 0 );
  virtual ~ContourSubElement_Tet_4() {}

  virtual int side_facets( Faceted_Surface & facets, int side ) const;

private:
  //: Default constructor not allowed
  ContourSubElement_Tet_4();

  int conformal_decomposition();

  int process_edge( const int i0,
                    const int i1,
                    const int i2,
                    std::vector<int> & is_on_surf,
                    PointVec & lnodes,
                    const std::vector<double> & ldist );

  bool is_degenerate( const std::vector<int> & edge_node_ids,
                      const int i0, const int i1, const int i2, const int i3 );
};

class ContourSubElement_Adaptive_Tet_4 : public ContourSubElement {
public:
  ContourSubElement_Adaptive_Tet_4( const PointVec & coords,
                             const std::vector<int> & side_ids,
                             const ContourElement * in_owner,
                             const int in_subelement_depth = 0);
  ContourSubElement_Adaptive_Tet_4( const PointVec & coords,
                             const std::vector<int> & side_ids,
                             const std::vector<int> & edge_age,
                             const ContourElement * in_owner,
                             const int in_subelement_depth = 0);
  virtual ~ContourSubElement_Adaptive_Tet_4() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Adaptive_Tet_4();

  static const int MAX_REFINMENT_LEVELS;

  int non_conformal_decomposition();

  std::vector<int> my_edge_age;
};

class ContourSubElement_Tet_10 : public ContourSubElement {
public:
  ContourSubElement_Tet_10( const PointVec & coords,
                     const std::vector<int> & side_ids,
                     const ContourElement * in_owner,
                     const int in_subelement_depth = 0,
                     const int subelement_sign = 0 );
  virtual ~ContourSubElement_Tet_10() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Tet_10();
};

} // namespace krino

#endif // Akri_ContourSubElement_h
