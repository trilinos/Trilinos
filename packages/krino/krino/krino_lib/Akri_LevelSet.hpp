// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_LevelSet_h
#define Akri_LevelSet_h

/**

  Create one of these for every interface you want to capture.

 */
#include <stk_util/diag/Timer.hpp>

#include <Akri_BoundingBox.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_LevelSet_Identifier.hpp>
#include <Akri_Vec.hpp>

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <cmath>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace sierra { namespace Sctl { class Event; } }
namespace krino { class AuxMetaData; }
namespace krino { class IC_Alg; }
namespace krino { class ParallelErrorMessage; }

namespace krino {

enum Redistance_Method
{
  CLOSEST_POINT=0,
  FAST_MARCHING,
  MAX_REDISTANCE_METHOD_TYPE
};

/// Return true if field-data exists for the specified meshobj and field.
bool all_nodes_have_field_data(const stk::mesh::BulkData& stk_bulk, stk::mesh::Entity entity, const stk::mesh::FieldBase& field);

class LevelSet {
friend class LevelSet_Size;
public:
  stk::mesh::MetaData & meta();
  const stk::mesh::MetaData & meta() const;
  stk::mesh::BulkData & mesh();
  const stk::mesh::BulkData & mesh() const;
  AuxMetaData & aux_meta();
  const AuxMetaData & aux_meta() const;

  const std::string & name() const { return my_name; }
  LevelSet_Identifier get_identifier() const {return my_identifier; }
  stk::diag::Timer & get_timer() const { return my_timer; }
  stk::diag::Timer & get_parent_timer() const { return my_parent_timer; }

  // finalize setup
  static void setup(stk::mesh::MetaData & meta);
  static void post_commit_setup(stk::mesh::MetaData & meta);
  virtual void setup();

  static void set_current_coordinates(stk::mesh::MetaData & meta, const FieldRef ref);
  static FieldRef get_current_coordinates(stk::mesh::MetaData & meta);

  void register_fields();

  void advance_semilagrangian(const double deltaTime);

  static void gather_nodal_field(
    const stk::mesh::BulkData& stk_mesh,
    stk::mesh::Entity obj,
    const FieldRef  & field_ref,
    double * field);

  double compute_average_edge_length() const;

  void build_facets_locally(const stk::mesh::Selector & selector);

  void compute_sizes( double & area, double & neg_vol, double & pos_vol, const double distance = 0.0 );
  double gradient_magnitude_error();
  void compute_continuous_gradient() const;

  void compute_distance( stk::mesh::Entity n,
			 const double & deltaTime ) const;
  void estimate_error();

  // hack to dump facet list to exoii databse.
  void facets_exoii();
  void facets_exoii(Faceted_Surface & cs);

  bool elem_on_interface(stk::mesh::Entity e) const;

  void snap_to_mesh() const;
  bool remove_wall_features() const;
  bool simple_remove_wall_features() const;
  bool elem_has_field_data(const FieldRef &myField, const stk::mesh::Entity &elem) const;
  double calc_elem_len_normal(const stk::mesh::Entity &elem, const stk::mesh::Entity &side, const FieldRef &coordinates_field) const;


  //--------------------------------------------------------------------------------
  // queries
  //--------------------------------------------------------------------------------

  const std::string & get_distance_name() const { return my_distance_name; }
  void set_distance_name( const std::string & distance_name ) { my_distance_name = distance_name; }

  const FieldRef & get_distance_field() const { return my_distance_field; }
  void set_distance_field( const FieldRef & ref ) { my_distance_field = ref; }

  const FieldRef & get_old_distance_field() const { return my_old_distance_field; }
  void set_old_distance_field( const FieldRef & ref ) { my_old_distance_field = ref; }

  void set_extension_velocity( const Vector3d & extension_velocity ) { my_extension_velocity = extension_velocity; }
  const Vector3d & get_extension_velocity() const { return my_extension_velocity; }

  const FieldRef & get_coordinates_field() const { return my_coordinates_field; }

  const std::string & get_isovar_name() const { return my_isovar_name; }

  const std::string & get_composite_name() const { return my_composite_name; }
  void set_composite_name( const std::string & composite_name ) {
    my_composite_name = composite_name;
    std::transform(my_composite_name.begin(), my_composite_name.end(), my_composite_name.begin(), ::toupper);
  }

  const FieldRef & get_isovar_field() const { return my_isovar_field; }
  void set_isovar_field( const FieldRef & ref ) { my_isovar_field = ref; }

  double get_time_of_arrival_speed(stk::mesh::Entity elem, ParallelErrorMessage& err) const;

  void set_isovar(const std::string & isovar_name, const double isoval) { my_isovar_name = isovar_name; my_threshold = isoval; trackIsoSurface = true; }
  const double & get_isoval() const { return my_threshold; }

  bool get_reinitialize_every_step() const { return my_needs_reinitialize_every_step; }
  void set_reinitialize_every_step(const bool reinit) { set_keep_IC_surfaces(); my_needs_reinitialize_every_step = true; }

  Redistance_Method get_redistance_method() const { return my_redistance_method; }
  void set_redistance_method( const Redistance_Method type ) { my_redistance_method = type; }
  void set_time_of_arrival_element_speed_field_name( const std::string & time_of_arrival_speed_field_name) { my_time_of_arrival_element_speed_field_name = time_of_arrival_speed_field_name; }
  void set_time_of_arrival_block_speed(const std::string & blockName, const double blockSpeed);
  Faceted_Surface & get_facets() { return *facets; }
  const Faceted_Surface & get_facets() const { return *facets; }

  static std::vector<std::string> the_levelSet_names;
  static LevelSet_Identifier get_identifier(const std::string & name);
  static std::string & get_name(const LevelSet_Identifier identifier) { ThrowAssert(identifier.get() < the_levelSet_names.size()); return the_levelSet_names[identifier.get()]; }

  void narrow_band_multiplier( double multiplier ) { my_narrow_band_multiplier = multiplier; }
  const double & narrow_band_size() const { return my_narrow_band_size; }
  void narrow_band_size( double size ) { my_narrow_band_size = size; } // publicly deprecated, make private

  void max_feature_size( double size ) { my_max_feature_size= size; }
  void use_simple_remove_feature( bool is_simple ) {my_use_simple_remove_feature = is_simple; }

  double max_feature_size() { return my_max_feature_size; }
  bool use_simple_remove_feature() { return my_use_simple_remove_feature; }

  void set_surface_parts_vector();

  void set_ic_offset (const double offset) {
    my_ic_offset = offset;
  }
  void set_ic_scale (const double scale) {
    my_ic_scale = scale;
  }
  void perform_initial_redistance(const bool flag) {
    my_perform_initial_redistance = flag;
  }

  void set_keep_IC_surfaces() { my_keep_IC_surfaces = true; }
  bool get_keep_IC_surfaces() const { return my_keep_IC_surfaces; }

  static bool sign_change( double f1, double f2 ) {
    return ( (f1 < 0.) ? (f2 >= 0.) : (f2 < 0.) ); // GOMA sign convention
    //return ( (f1 > 0.) ? (f2 <= 0.) : (f2 > 0.) ); // Marching cubes sign convention
  }

  static int sign( double f ) {
    return ( (f < 0.) ? -1 : 1 ); // GOMA sign convention
    //return ( (f > 0.) ? 1 : -1 ); // Marching cubes sign convention
  }

  const std::vector<stk::mesh::Part *> & get_compute_surface_distance_parts() const { return my_compute_surface_distance_parts; }
  std::vector<stk::mesh::Part *> & get_compute_surface_distance_parts() { return my_compute_surface_distance_parts; }
  void set_surface_distance(std::vector<stk::mesh::Part *> surfaces, const double in_distance);
  void compute_surface_distance(const double narrowBandSize=0.0, const double farFieldValue=0.0);
  static void initialize(stk::mesh::MetaData & meta, const bool requires_additional_initialization);
  void initialize(const double time = 0.0, const bool requires_additional_initialization = false);
  void redistance();
  void redistance(const stk::mesh::Selector & selector);
  void fast_marching_redistance(const stk::mesh::Selector & selector, const bool compute_time_of_arrival = false);

  void set_initial_volume(const double v) { my_initial_neg_vol = v; }
  double constrained_redistance(const bool use_initial_vol = false);

  void compute_nodal_bbox( const stk::mesh::Selector & selector,
    BoundingBox & node_bbox,
    const Vector3d & displacement = Vector3d::ZERO ) const;

  double find_redistance_correction( const double start_area,
				   const double start_neg_vol,
				   const double start_pos_vol,
				   const int max_iterations = 100,
				   const double tol = 1.e-6 );

  bool has_IC_surfaces();
  BoundingBox get_IC_surface_bounding_box();
  IC_Alg& get_IC_alg();

  static LevelSet & build(stk::mesh::MetaData & in_meta, const std::string & ls_name, stk::diag::Timer & parent_timer);

  virtual ~LevelSet();

private:
  LevelSet(stk::mesh::MetaData & in_meta, const std::string & in_name, stk::diag::Timer & parent_timer);

private:
  stk::mesh::MetaData & my_meta;
  AuxMetaData & my_aux_meta;
  const LevelSet_Identifier my_identifier;
  const std::string my_name;
  mutable stk::diag::Timer my_parent_timer;
  mutable stk::diag::Timer my_timer;

public:
  const unsigned spatial_dimension;

private:

  FieldRef my_coordinates_field;
  FieldRef my_distance_field;
  FieldRef my_old_distance_field;
  FieldRef my_isovar_field;
  FieldRef myTimeOfArrivalElementSpeedField;

  std::string my_distance_name;
  std::string my_isovar_name;

  std::string my_composite_name;

  double my_narrow_band_multiplier;
  double my_narrow_band_size;

  double my_max_feature_size;
  bool my_use_simple_remove_feature = false;


  double my_ic_offset;
  double my_ic_scale;
  bool my_perform_initial_redistance;
  bool my_keep_IC_surfaces;

  double my_threshold;
  Redistance_Method my_redistance_method;
  std::string my_time_of_arrival_element_speed_field_name;
  std::map<std::string, double> myTimeOfArrivalBlockSpeedsByName;
  std::vector<double> myTimeOfArrivalBlockSpeeds;

  std::unique_ptr<IC_Alg> my_IC_alg;

  // vector of previous facets
  std::unique_ptr<Faceted_Surface> facets_old;

    // vector of current facets
  std::unique_ptr<Faceted_Surface> facets;

  Vector3d my_extension_velocity;
  const double epsilon;

  bool trackIsoSurface;

  // used to increment file name for facet exoii database hack
  int my_facetFileIndex;

  double my_initial_neg_vol;

  bool my_needs_reinitialize_every_step;

  std::vector<stk::mesh::Part *> my_compute_surface_distance_parts;
  std::vector<stk::mesh::Part *> my_surface_parts;


  void set_distance(const double & distance)  const;
  void increment_distance(const double increment, const bool enforce_sign = false)  const;
  void scale_distance(const double scale)  const;
  void negate_distance()  const;

  void time_integrate(const double deltaTime);

  void prepare_to_compute_distance( const double & deltaTime, const stk::mesh::Selector & selector );

  void compute_distance_semilagrangian(const double & deltaTime, const stk::mesh::Selector & selector );

  double distance( const Vector3d & x,
		 const int previous_sign,
		 const bool enforce_sign ) const;

  bool compute_time_of_arrival() const { return !my_time_of_arrival_element_speed_field_name.empty() || !myTimeOfArrivalBlockSpeeds.empty(); }
};

class LevelSetManager {
public:
  LevelSetManager() {}
  LevelSetManager(LevelSetManager const&) = delete;
  LevelSetManager& operator=(LevelSetManager const&) = delete;

  static LevelSetManager & get(stk::mesh::MetaData & meta);
  static LevelSetManager & get(const stk::mesh::MetaData & meta);

  int numberLevelSets() const { return my_level_sets.size(); }
  LevelSet & levelSet(const unsigned ordinal) const { ThrowAssert(ordinal < my_level_sets.size()); return *my_level_sets[ordinal]; }
  bool has_levelSet(const std::string & ls_name) const;
  void add(LevelSet * ls) { my_level_sets.emplace_back(ls); }

  std::vector< std::unique_ptr<LevelSet> >::const_iterator begin() const { return my_level_sets.begin(); }
  std::vector< std::unique_ptr<LevelSet> >::const_iterator end() const { return my_level_sets.end(); }

  void set_current_coordinates(FieldRef current_coords) { my_current_coordinates = current_coords; }
  FieldRef get_current_coordinates() const { return my_current_coordinates; }
private:
  std::vector< std::unique_ptr<LevelSet> > my_level_sets;
  FieldRef my_current_coordinates;
};

std::string print_sizes(const LevelSet & ls);

} // namespace krino

#endif // Akri_LevelSet_h
