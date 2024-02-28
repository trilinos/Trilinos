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
#include <stk_math/StkVector.hpp>
#include "Akri_MasterElementDeterminer.hpp"
#include "Akri_TopologyData.hpp"

namespace krino {

template<class T>
bool compare_ptr_by_global_id (const T* i, const T* j) { return (i->global_id() < j->global_id()); }

class ContourSubElement {
public:

  ContourSubElement(const ContourElement * owner,
              const int in_subelement_depth,
              const int subelement_sign );

  virtual ~ContourSubElement();

  virtual const MasterElement& get_master_element() const = 0;
  virtual const MasterElement& get_side_master_element() const = 0;
  virtual int get_num_nodes() const = 0;
  virtual int get_num_sides() const = 0;
  virtual const int * get_side_ids() const = 0;
  virtual const double * get_distance_at_nodes() const = 0;
  virtual const stk::math::Vector3d * get_coordinates_at_nodes() const = 0;

  static bool is_more(const stk::math::Vector3d & v1, const stk::math::Vector3d & v2);

  double compute_area_of_interface() const;
  double compute_relative_signed_volume(const int signOfDomain) const;

  double relative_volume() const;
  double parametric_quality() const;
  double physical_quality() const;
  double side_quality(const int side) const;

  int num_intg_pts( const int intg_pt_sign );

  int gather_intg_pts( const int intg_pt_sign,
		       sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
		       sierra::ArrayContainer<double,NINT> & intg_weights,
		       sierra::ArrayContainer<double,NINT> & determinants,
		       int index = 0 );

  int build_facets( FacetedSurfaceBase & facets );

  // default implementation
  virtual int side_facets( FacetedSurfaceBase & facets, int side ) const;
  virtual double side_area( int side ) const;

  stk::topology topology() const { return get_master_element().get_topology(); }

  int spatial_dim() const { return my_owner->spatial_dim(); }
  int num_subelements() const { return my_subelements.size(); }

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

  std::vector< ContourSubElement * > my_subelements;
  const ContourElement * my_owner;
  int my_subelement_depth; // depth down the tree of subelements (0 for the base_subelement)
  int my_sign; // -1 for elements on negative side, +1 for positive side, 0 for non-conformal elements spanning interface

private:
  //: Default constructor not allowed
  ContourSubElement();
};

template<stk::topology::topology_t TOPO>
class MasterElementHolder
{
public:
  static const MasterElement & get_master_element()
  {
    static const MasterElement * me = nullptr;
    if (nullptr == me)
      me = &MasterElementDeterminer::getMasterElement(stk::topology(TOPO));
    return *me;
  }
  static const MasterElement & get_side_master_element()
  {
    static const MasterElement * sideMe = nullptr;
    if (nullptr == sideMe)
      sideMe = &MasterElementDeterminer::getMasterElement(stk::topology(TOPO).side_topology());
    return *sideMe;
  }
};

template<stk::topology::topology_t TOPO>
class ContourSubElementWithTopology : public ContourSubElement {
public:
  static constexpr unsigned NUM_NODES = TopologyData<TOPO>::num_nodes();
  static constexpr unsigned NUM_SIDES = TopologyData<TOPO>::num_sides();

  ContourSubElementWithTopology( const std::array<stk::math::Vector3d,NUM_NODES> & coords,
              const std::array<int,NUM_SIDES> & sideIds,
              const ContourElement * owner,
              const int subelementDepth,
              const int subelementSign )
  : ContourSubElement(owner, subelementDepth, subelementSign),
    myCoords(coords),
    mySideIds(sideIds)
  {
    for ( unsigned i = 0; i < NUM_NODES; i++ )
      myDist[i] = my_owner->distance( coords[i] );
  }

  virtual ~ContourSubElementWithTopology() {}

  virtual int get_num_nodes() const override { return NUM_NODES; }
  virtual int get_num_sides() const override { return NUM_SIDES; }
  virtual const int * get_side_ids() const override { return mySideIds.data(); }
  virtual const double * get_distance_at_nodes() const override { return myDist.data(); }
  virtual const stk::math::Vector3d * get_coordinates_at_nodes() const override { return myCoords.data(); }
  virtual const MasterElement & get_master_element() const override { return MasterElementHolder<TOPO>::get_master_element(); }
  virtual const MasterElement & get_side_master_element() const override { return MasterElementHolder<TOPO>::get_side_master_element(); }

protected:
  std::array<stk::math::Vector3d,NUM_NODES> myCoords;
  std::array<int,NUM_SIDES> mySideIds;
  std::array<double,NUM_NODES> myDist;
};

class ContourSubElement_Quad_4 : public ContourSubElementWithTopology<stk::topology::QUAD_4_2D> {
public:
  ContourSubElement_Quad_4( const std::array<stk::math::Vector3d,4> & coords,
		     const std::array<int,4> & sideIds,
		     const ContourElement * in_owner );
  virtual ~ContourSubElement_Quad_4() {}

private:
  int non_conformal_decomposition();
};

class ContourSubElement_Quad_9 : public ContourSubElementWithTopology<stk::topology::QUAD_9_2D> {
public:
  ContourSubElement_Quad_9( const std::array<stk::math::Vector3d,9> & coords,
		     const std::array<int,4> & sideIds,
		     const ContourElement * in_owner );

  virtual ~ContourSubElement_Quad_9() {}

private:
  int non_conformal_decomposition();
};

class ContourSubElement_Hex_8 : public ContourSubElementWithTopology<stk::topology::HEX_8> {
public:
  ContourSubElement_Hex_8( const std::array<stk::math::Vector3d,8> & coords,
		    const std::array<int,6> & sideIds,
		    const ContourElement * in_owner );
  virtual ~ContourSubElement_Hex_8() {}

private:
  int subpyramid_non_conformal_decomposition( const int face );
};

class ContourSubElement_Hex_27 : public ContourSubElementWithTopology<stk::topology::HEX_27> {
public:
  ContourSubElement_Hex_27( const std::array<stk::math::Vector3d,27> & coords,
		     const std::array<int,6> & sideIds,
		     const ContourElement * in_owner );
  virtual ~ContourSubElement_Hex_27() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Hex_27();

  int subpyramid_non_conformal_decomposition( const int face );
};

class ContourSubElement_Wedge_6 : public ContourSubElementWithTopology<stk::topology::WEDGE_6> {
public:
  ContourSubElement_Wedge_6( const std::array<stk::math::Vector3d,6> & coords,
                    const std::array<int,5> & sideIds,
                    const ContourElement * in_owner );
  virtual ~ContourSubElement_Wedge_6() {}

private:
  int subpyramid_non_conformal_decomposition( const int face );
};

class ContourSubElement_Tri_3 : public ContourSubElementWithTopology<stk::topology::TRI_3_2D> {
public:
  ContourSubElement_Tri_3( const std::array<stk::math::Vector3d,3> & coords,
		    const std::array<int,3> & sideIds,
		    const ContourElement * in_owner,
                    const int in_subelement_depth = 0,
                    const int subelement_sign = 0 );
  virtual ~ContourSubElement_Tri_3() {}

  virtual int side_facets( FacetedSurfaceBase & facets, int side ) const override;
  double side_area( int side ) const override;

private:
  //: Default constructor not allowed
  ContourSubElement_Tri_3();

  int conformal_decomposition();

  int process_edge( const int i0,
                    const int i1,
                    const int i2,
                    std::array<int,6> & is_on_surf,
                    std::array<stk::math::Vector3d,6> & lnodes,
                    const std::array<double,3> & ldist );

  bool is_degenerate( const std::array<int,6> & edge_node_ids,
                      const int i0, const int i1, const int i2 );
};

class ContourSubElement_Adaptive_Tri_3 : public ContourSubElementWithTopology<stk::topology::TRI_3_2D> {
public:
  ContourSubElement_Adaptive_Tri_3( const std::array<stk::math::Vector3d,3> & coords,
                             const std::array<int,3> & sideIds,
                             const ContourElement * in_owner,
                             const int in_subelement_depth = 0);
  ContourSubElement_Adaptive_Tri_3( const std::array<stk::math::Vector3d,3> & coords,
                             const std::array<int,3> & sideIds,
                             const std::array<int,3> & edge_age,
                             const ContourElement * in_owner,
                             const int in_subelement_depth = 0);
  virtual ~ContourSubElement_Adaptive_Tri_3() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Adaptive_Tri_3();

  static const int MAX_REFINMENT_LEVELS;

  int non_conformal_decomposition();

  std::array<int,3> my_edge_age;
};

class ContourSubElement_Tri_6 : public ContourSubElementWithTopology<stk::topology::TRI_6_2D> {
public:
  ContourSubElement_Tri_6( const std::array<stk::math::Vector3d,6> & coords,
                    const std::array<int,3> & sideIds,
                    const ContourElement * in_owner,
                    const int in_subelement_depth = 0,
                    const int subelement_sign = 0 );
  virtual ~ContourSubElement_Tri_6() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Tri_6();
};

class ContourSubElement_Tet_4 : public ContourSubElementWithTopology<stk::topology::TET_4> {
public:
  ContourSubElement_Tet_4( const std::array<stk::math::Vector3d,4> & coords,
                    const std::array<int,4> & sideIds,
                    const ContourElement * in_owner,
                    const int in_subelement_depth = 0,
                    const int subelement_sign = 0 );
  virtual ~ContourSubElement_Tet_4() {}

  virtual int side_facets( FacetedSurfaceBase & facets, int side ) const override;
  virtual double side_area( int side ) const override;

private:
  //: Default constructor not allowed
  ContourSubElement_Tet_4();

  int conformal_decomposition();

  int process_edge( const int i0,
                    const int i1,
                    const int i2,
                    std::array<int,10> & is_on_surf,
                    std::array<stk::math::Vector3d,10> & lnodes,
                    const std::array<double,4> & ldist );

  bool is_degenerate( const std::array<int,10> & edge_node_ids,
                      const int i0, const int i1, const int i2, const int i3 );
};

class ContourSubElement_Adaptive_Tet_4 : public ContourSubElementWithTopology<stk::topology::TET_4> {
public:
  ContourSubElement_Adaptive_Tet_4( const std::array<stk::math::Vector3d,4> & coords,
                             const std::array<int,4> & sideIds,
                             const ContourElement * in_owner,
                             const int in_subelement_depth = 0);
  ContourSubElement_Adaptive_Tet_4( const std::array<stk::math::Vector3d,4> & coords,
                             const std::array<int,4> & sideIds,
                             const std::array<int,6> & edge_age,
                             const ContourElement * in_owner,
                             const int in_subelement_depth = 0);
  virtual ~ContourSubElement_Adaptive_Tet_4() {}

private:
  //: Default constructor not allowed
  ContourSubElement_Adaptive_Tet_4();

  static const int MAX_REFINMENT_LEVELS;

  int non_conformal_decomposition();

  std::array<int,6> my_edge_age;
};

class ContourSubElement_Tet_10 : public ContourSubElementWithTopology<stk::topology::TET_10> {
public:
  ContourSubElement_Tet_10( const std::array<stk::math::Vector3d,10> & coords,
                     const std::array<int,4> & sideIds,
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
