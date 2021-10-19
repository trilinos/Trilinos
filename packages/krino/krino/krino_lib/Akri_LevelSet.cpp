// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_LevelSet.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Compute_Surface_Distance.hpp>
#include <Akri_ContourElement.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_TypeDefs.hpp>
#include <Akri_Facet.hpp>
#include <Akri_IC_Alg.hpp>
#include <Akri_Fast_Marching.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_ParallelErrorMessage.hpp>

#include <math.h>
#include <map>
#include <iomanip>
#include <cstdio>

#include <Ioss_SubSystem.h>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <Akri_MasterElementDeterminer.hpp>

namespace krino {

bool all_nodes_have_field_data(const stk::mesh::BulkData& stk_bulk, stk::mesh::Entity entity, const stk::mesh::FieldBase& field)
{
  const unsigned nnodes = stk_bulk.num_nodes(entity);
  const stk::mesh::Entity* nodes = stk_bulk.begin_nodes(entity);
  for (unsigned i = 0; i < nnodes; ++i )
  {
    if (field_bytes_per_entity(field, nodes[i]) == 0)
    {
      return false;
    }
  }
  return true;
}

std::vector<std::string> LevelSet::the_levelSet_names;

LevelSet_Identifier LevelSet::get_identifier(const std::string & name)
{
  std::string upper_name = name;
  std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(), ::toupper);
  for (unsigned i=0; i<the_levelSet_names.size(); ++i)
  {
    if (upper_name == the_levelSet_names[i])
    {
      return LevelSet_Identifier(i);
    }
  }

  const unsigned id = the_levelSet_names.size();
  the_levelSet_names.push_back(upper_name);
  return LevelSet_Identifier(id);
}

stk::mesh::BulkData & LevelSet::mesh()
{
  return my_meta.mesh_bulk_data();
}
const stk::mesh::BulkData & LevelSet::mesh() const
{
  return my_meta.mesh_bulk_data();
}

stk::mesh::MetaData & LevelSet::meta()
{
  return my_meta;
}
const stk::mesh::MetaData & LevelSet::meta() const
{
  return my_meta;
}

AuxMetaData & LevelSet::aux_meta()
{
  return my_aux_meta;
}
const AuxMetaData & LevelSet::aux_meta() const
{
  return my_aux_meta;
}

bool LevelSet::has_IC_surfaces()
{
  IC_Alg& ic_alg = get_IC_alg();
  return ic_alg.numberSurfaces() > 0;
}

BoundingBox LevelSet::get_IC_surface_bounding_box()
{
  IC_Alg& ic_alg = get_IC_alg();
  return ic_alg.get_surface_bounding_box();
}

IC_Alg& LevelSet::get_IC_alg()
{
  if (!my_IC_alg)
  {
    my_IC_alg = std::make_unique<IC_Alg>(*this);
  }
  return *my_IC_alg;
}

void LevelSet::set_current_coordinates(stk::mesh::MetaData & meta, const FieldRef ref)
{
  LevelSetManager & ls_manager = LevelSetManager::get(meta);
  ls_manager.set_current_coordinates(ref);
}

FieldRef LevelSet::get_current_coordinates(stk::mesh::MetaData & meta)
{
  LevelSetManager & ls_manager = LevelSetManager::get(meta);
  FieldRef current_coords = ls_manager.get_current_coordinates();
  if (!current_coords.valid())
  {
    const stk::mesh::FieldBase * meta_coords = meta.coordinate_field();
    ThrowRequireMsg(nullptr != meta_coords, "Coordinates must be defined before calling LevelSet::get_current_coordinates().");
    ls_manager.set_current_coordinates(meta_coords);
    current_coords = ls_manager.get_current_coordinates();
  }
  return current_coords;
}

void LevelSet::setup(stk::mesh::MetaData & meta)
{
  const LevelSetManager & region_ls = LevelSetManager::get(meta);
  for (auto&& ls : region_ls)
  {
    ls->setup();
  }
}

void LevelSet::post_commit_setup(stk::mesh::MetaData & meta)
{
  const double max_elem_size = compute_maximum_element_size(meta.mesh_bulk_data());

  const LevelSetManager & region_ls = LevelSetManager::get(meta);
  for (auto&& ls : region_ls)
  {
    ls->setup();
    if (ls->my_narrow_band_multiplier > 0.)
    {
      ls->narrow_band_size(ls->my_narrow_band_multiplier * max_elem_size);
    }
    else if (ls->narrow_band_size() > 0.)
    {
      const double narrow_band = ls->narrow_band_size();
      ThrowErrorMsgIf(!(narrow_band > max_elem_size),
          "Currently, narrow_band_size must be greater than the maximum element size of " << max_elem_size << std::endl
          << "in order to avoid unintentional accuracy degradation. If this feature is needed, please contact krino developers.");
    }
  }
}


void LevelSet::setup(void)
{
  register_fields();

  facets.reset(new Faceted_Surface("current_facets"));
  facets_old.reset(new Faceted_Surface("old_facets"));

  // initializes CDFEM_Support
  if (!meta().is_commit())
  {
    CDFEM_Support::get(my_meta);
  }
}
//--------------------------------------------------------------------------------
void LevelSet::register_fields(void)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::register_fields(void)"); /* %TRACE% */

  if (trackIsoSurface || aux_meta().has_field(stk::topology::NODE_RANK, my_distance_name))
  {
    // non-krino region
    if (trackIsoSurface)
    {
      if (aux_meta().has_field(stk::topology::NODE_RANK, my_isovar_name))
      {
        set_isovar_field( aux_meta().get_field(stk::topology::NODE_RANK, my_isovar_name) );
      }
      else if (aux_meta().has_field(stk::topology::ELEMENT_RANK, my_isovar_name))
      {
        set_isovar_field( aux_meta().get_field(stk::topology::ELEMENT_RANK, my_isovar_name) );
      }
      else
      {
        ThrowErrorMsgIf(
            true, "Isosurface variable '" << my_isovar_name << "' should already be registered.");
      }
    }
    else
    {
      const FieldRef distance_ref = aux_meta().get_field(stk::topology::NODE_RANK, my_distance_name);
      if ( distance_ref.number_of_states() == 1 )
      {
        set_distance_field( distance_ref );
        set_old_distance_field( distance_ref );
      }
      else
      {
        set_distance_field( distance_ref.field_state(stk::mesh::StateNew) );
        set_old_distance_field( distance_ref.field_state(stk::mesh::StateOld) );
      }

      set_isovar_field( my_distance_field );
    }
  }
  else
  {
    if (aux_meta().using_fmwk())
    {
      ThrowRuntimeError("ERROR: field " << my_distance_name << " is not registered for level set " << name() << ".  "
        << "This can be caused by an incorrect Distance Variable.  "
        << "Or, in aria, there could be a conflict between the specified subindex and the named species.  "
        << "If so, try using a subindex greater than the number of species for your level set.");
    }

    // krino region
    const FieldType & type_double  = FieldType::REAL;

    if(krinolog.shouldPrint(LOG_DEBUG))
      krinolog << "KRINO: Registering distance variable with name '" << my_distance_name << "'." << stk::diag::dendl;

    const bool cdfem_is_active = krino::CDFEM_Support::is_active(meta());
    if (cdfem_is_active)
    {
      Phase_Support phase_support = Phase_Support::get(meta());
      for (auto partPtr : meta().get_mesh_parts())
      {
        if (partPtr->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
            phase_support.level_set_is_used_by_nonconformal_part(this, phase_support.find_nonconformal_part(*partPtr)))
        {
          FieldRef distance_ref = aux_meta().register_field( my_distance_name, type_double, stk::topology::NODE_RANK, 1, 1, *partPtr );
          set_old_distance_field( distance_ref );
          set_distance_field( distance_ref );
        }
      }
    }
    else
    {
      FieldRef distance_ref = aux_meta().register_field( my_distance_name, type_double, stk::topology::NODE_RANK, 2, 1, meta().universal_part() );
      set_old_distance_field( FieldRef( distance_ref, stk::mesh::StateOld ) );
      set_distance_field( FieldRef( distance_ref, stk::mesh::StateNew ) );
    }

    set_isovar_field( my_distance_field );
  }

  if (!myTimeOfArrivalBlockSpeedsByName.empty())
  {
    myTimeOfArrivalBlockSpeeds.resize(my_meta.get_parts().size());
    for (auto entry : myTimeOfArrivalBlockSpeedsByName)
    {
      if (my_aux_meta.has_part(entry.first))
      {
        myTimeOfArrivalBlockSpeeds[my_aux_meta.get_part(entry.first).mesh_meta_data_ordinal()] = entry.second;
      }
      else
      {
        ThrowErrorMsgIf(true, "Could not find block " << entry.first << " when setting speed for computing time-of-arrival.");
      }
    }
  }

  if (!my_time_of_arrival_element_speed_field_name.empty())
  {
    const bool hasSpeedField = aux_meta().has_field(stk::topology::ELEMENT_RANK, my_time_of_arrival_element_speed_field_name);
    ThrowErrorMsgIf(!hasSpeedField, "Could not find element speed field " << my_time_of_arrival_element_speed_field_name << " for computing time-of-arrival.");
    myTimeOfArrivalElementSpeedField = aux_meta().get_field(stk::topology::ELEMENT_RANK, my_time_of_arrival_element_speed_field_name);
    ThrowRequireMsg(myTimeOfArrivalBlockSpeeds.empty(), "Speed for time-of-arrival calculation should be specified via element speed or block speed (not both).");
  }
}

//-----------------------------------------------------------------------------------

void
LevelSet::set_time_of_arrival_block_speed(const std::string & blockName, const double blockSpeed)
{
  std::string lowerBlockName = blockName;
  std::transform(lowerBlockName.begin(), lowerBlockName.end(), lowerBlockName.begin(), ::tolower);
  auto entry = myTimeOfArrivalBlockSpeedsByName.find(lowerBlockName);
  ThrowRequireMsg(entry == myTimeOfArrivalBlockSpeedsByName.end(), "Speed for block " << blockName << " specified more than once.");
  myTimeOfArrivalBlockSpeedsByName[lowerBlockName] = blockSpeed;
}

//-----------------------------------------------------------------------------------
void
LevelSet::facets_exoii(void)
{
  /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::facets_exoii(void)"); /* %TRACE% */
  Faceted_Surface & f = *facets;
  facets_exoii(f);
}


void
LevelSet::facets_exoii(Faceted_Surface & cs)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::facets_exoii(Faceted_Surface & cs)"); /* %TRACE% */

  int nfaces = cs.size();
  const int nodes_per_elem = spatial_dimension;
  const int nnodes = nfaces * nodes_per_elem;
  const int nelems = nfaces * 1; //writing out faces as elements

  // construct file name
  ThrowAssert(my_facetFileIndex <= 99999);
  char counterChar[6];
  std::sprintf(counterChar, "%.5d", my_facetFileIndex);
  std::string fn = std::string("facets.") + name() + ".e-s" + std::string(counterChar);

  // Type (ExodusII) is hard-wired at this time.
  Ioss::DatabaseIO *db = Ioss::IOFactory::create("exodusII", fn, Ioss::WRITE_RESULTS);
  Ioss::Region io(db, "FacetRegion");

  // find offsets for consistent global numbering
  int elem_start = 0;
  if ( mesh().parallel_size() > 1 ) {
    std::vector< int > num_faces;
    num_faces.resize( mesh().parallel_size() );

    // Gather everyone's facet list sizes
    // I don't think there is a framework call for this...
    MPI_Allgather( &nfaces, 1, MPI_INT,
		   &(num_faces[0]), 1, MPI_INT, mesh().parallel() );
    for (int i = 0; i< mesh().parallel_rank(); ++i) {
      elem_start += num_faces[i];
    }
  }

  const std::string description = "level set interface facets";
  io.property_add(Ioss::Property("title", description));
  io.begin_mode(Ioss::STATE_DEFINE_MODEL);

  // if we have no elements bail now
  if ( 0 == nelems ) {
    io.end_mode(Ioss::STATE_DEFINE_MODEL);
    my_facetFileIndex++;
    return;
  }

  Ioss::NodeBlock *nb = new Ioss::NodeBlock(db, "nodeblock_1", nnodes, spatial_dimension);
  io.add(nb);

  std::string el_type;
  if (spatial_dimension == 3)
    el_type = "trishell3";
  else
    el_type = "shellline2d2";

  Ioss::ElementBlock *eb = new Ioss::ElementBlock(db, "block_1", el_type, nelems);
  io.add(eb);

  io.end_mode(Ioss::STATE_DEFINE_MODEL);
  io.begin_mode(Ioss::STATE_MODEL);

  const FacetOwningVec & cs_surfaces = cs.get_facets();

  // loop over elements and node to create maps
  {
    std::vector< int > nmap, emap;
    nmap.reserve(nnodes);
    emap.reserve(nelems);

    for ( unsigned n=0, e=0; e<cs_surfaces.size(); ++e ) {
      emap.push_back(elem_start + e + 1);
      for ( int j = 0; j < nodes_per_elem; ++j ) {
	nmap.push_back(nodes_per_elem*elem_start + n + 1);
	++n;
      }
    }
    nb->put_field_data("ids", nmap);
    eb->put_field_data("ids", emap);
  }

  // generate coordinates
  {
    std::vector< double > xyz;
    xyz.reserve(spatial_dimension*nnodes);

    for ( auto&& cs_surface : cs_surfaces ) {
      for ( int j = 0; j < nodes_per_elem; ++j ) {
	const Vector3d & vert = cs_surface->facet_vertex(j);
	xyz.push_back(vert[0]);
	xyz.push_back(vert[1]);
	if (3 == spatial_dimension) xyz.push_back(vert[2]);
      }
    }
    nb->put_field_data("mesh_model_coordinates", xyz);
  }

  // generate connectivity
  {
    std::vector< int > conn;
    conn.reserve(nnodes);
    for ( int n = 0; n < nnodes; ++n ) {
      conn.push_back(nodes_per_elem*elem_start + n+1);
    }
    eb->put_field_data("connectivity", conn);
  }

  io.end_mode(Ioss::STATE_MODEL);
  my_facetFileIndex++;
}

//-----------------------------------------------------------------------------------
void
LevelSet::compute_surface_distance(const double narrowBandSize, const double farFieldValue)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::compute_surface_distance(void)"); /* %TRACE% */

  stk::mesh::Selector surface_selector = selectUnion(my_compute_surface_distance_parts);

  if (narrowBandSize > 0.0)
  {
    my_narrow_band_size = narrowBandSize;
  }

  const stk::mesh::Field<double>& coords = reinterpret_cast<const stk::mesh::Field<double>&>(get_coordinates_field().field());
  const stk::mesh::Field<double>& dist = reinterpret_cast<const stk::mesh::Field<double>&>(get_distance_field().field());

  Compute_Surface_Distance::calculate(
     mesh(),
     get_timer(),
     coords,
     dist,
     surface_selector,
     my_narrow_band_size,
     farFieldValue);

  // output for time 0 is from old
  if (!(my_distance_field == my_old_distance_field))
  {
    stk::mesh::field_copy(my_distance_field, my_old_distance_field);
  }

}

//-----------------------------------------------------------------------------------
void
LevelSet::set_surface_distance(std::vector<stk::mesh::Part *> surfaces, const double in_distance)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::compute_surface_distance(void)"); /* %TRACE% */

  const FieldRef dField = get_distance_field();

  const stk::mesh::Selector selector = stk::mesh::selectField(dField) & selectUnion(surfaces);
  stk::mesh::BucketVector const& buckets = mesh().get_buckets( stk::topology::NODE_RANK, selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const size_t length = b.size();
    double *dist = field_data<double>(dField , b);

    for ( size_t n = 0; n < length; ++n )
    {
      dist[n] = in_distance;
    }
  }
}

//-----------------------------------------------------------------------------------
void
LevelSet::advance_semilagrangian(const double deltaTime)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::advance_semilagrangian(const double deltaTime)"); /* %TRACE% */

  if(trackIsoSurface)
  {
    redistance();
  }
  else
  {
    stk::mesh::Selector selector(my_meta.universal_part());

    // store existing facets in facets_old
    facets->swap( *facets_old );

    // get non-local facets such that we have copies of all old facets
    // within the range of this proc's nodes
    prepare_to_compute_distance( deltaTime, selector );

    // compute nodal distances with semi-lagrangian step
    stk::mesh::field_copy(my_old_distance_field, my_distance_field);
    compute_distance_semilagrangian( deltaTime, selector );

    // build local facet list
    build_facets_locally(my_meta.universal_part());

    // debugging
    if (krinolog.shouldPrint(LOG_FACETS))
    {
      facets_exoii();
    }
  }
}

//-----------------------------------------------------------------------------------
void
LevelSet::initialize(stk::mesh::MetaData & meta, const bool requires_additional_initialization)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::initialize(void)"); /* %TRACE% */

  const LevelSetManager & region_ls = LevelSetManager::get(meta);
  for (auto&& ls : region_ls)
  {
    ls->initialize(0., requires_additional_initialization);
  }
}

//-----------------------------------------------------------------------------------
void
LevelSet::initialize(const double time, const bool requires_additional_initialization)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::initialize(void)"); /* %TRACE% */

  if (trackIsoSurface)
  {
    return;
  }

  if(get_isovar_field().valid()) get_isovar_field().field().sync_to_host();
  if(get_distance_field().valid()) get_distance_field().field().sync_to_host();
  if(get_old_distance_field().valid()) get_old_distance_field().field().sync_to_host();

  if(get_isovar_field().valid()) get_isovar_field().field().modify_on_host();
  if(get_distance_field().valid()) get_distance_field().field().modify_on_host();
  if(get_old_distance_field().valid()) get_old_distance_field().field().modify_on_host();

  if (!my_compute_surface_distance_parts.empty())
  {
    compute_surface_distance();
    return;
  }

  krinolog << "Initializing levelset " << name() << "..." << stk::diag::dendl;

  /* process analytic surfaces */
  if (my_IC_alg)
  {
    my_IC_alg->execute(time, requires_additional_initialization);
  }

  if (compute_time_of_arrival())
  {
    fast_marching_redistance(my_meta.universal_part(), true);
  }
  else if (my_perform_initial_redistance)
  {
    constrained_redistance();
  }

  // Offset initialized LS if requested
  if (my_ic_offset != 0.0) increment_distance(my_ic_offset, false);

  // Scale initialized LS if requested
  if (my_ic_scale != 1.0) scale_distance(my_ic_scale);

  stk::mesh::field_copy(get_distance_field(), get_old_distance_field());
}

double LevelSet::constrained_redistance(const bool use_initial_vol)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::constrained_redistance(const bool use_initial_vol)"); /* %TRACE% */

  // Steps:
  // 1. measure current volume
  // 2. perform regular redistance
  // 3. find offset needed to conserve volume
  // 4. increment nodal distance by this amount

  if (get_isovar_field().valid()) get_isovar_field().field().sync_to_host();
  if (get_distance_field().valid()) get_distance_field().field().sync_to_host();
  if (get_old_distance_field().valid()) get_old_distance_field().field().sync_to_host();

  if (get_isovar_field().valid()) get_isovar_field().field().modify_on_host();
  if (get_distance_field().valid()) get_distance_field().field().modify_on_host();
  if (get_old_distance_field().valid()) get_old_distance_field().field().modify_on_host();

  // measure current area and volumes
  double start_area, start_neg_vol, start_pos_vol;
  compute_sizes( start_area,start_neg_vol, start_pos_vol, 0. );
  if ( 0. == start_neg_vol || 0. == start_area )
  {
    krinolog << "Skipping redistancing operation.  Volume and/or area have zero value.  Area: "
        << start_area << " Volume: " << start_neg_vol << stk::diag::dendl;
    return 0.;
  }

  if (use_initial_vol)
  {
    krinolog << "Performing conserved redistancing..." << stk::diag::dendl;

    if (my_initial_neg_vol <= 0.0)
    {
      my_initial_neg_vol = start_neg_vol;
    }
    start_pos_vol = start_neg_vol + start_pos_vol - my_initial_neg_vol;
    start_neg_vol = my_initial_neg_vol;
  }
  else
  {
    krinolog << "Performing constrained redistancing..." << stk::diag::dendl;
  }

  // perform regular redistance
  redistance();

  krinolog << "Correcting for volume change:" << stk::diag::dendl;

  // find correction needed to conserve volume
  const double correction = find_redistance_correction( start_area, start_neg_vol, start_pos_vol );

  // update nodal distance field
  increment_distance( -correction, true );

  return my_initial_neg_vol;
}

//--------------------------------------------------------------------------------

static std::function<std::pair<double,double>(const double)> build_volume_error_function_with_derivative(LevelSet & ls, const double startingNegVol, const double startingPosVol)
{
  auto volume_error_function_with_derivative = [&ls, startingNegVol, startingPosVol](const double x)
    {
      double area, negVol, posVol;
      ls.compute_sizes( area, negVol, posVol, x );
      const double totVol = startingNegVol+startingPosVol;
      const double relativeError = (negVol - startingNegVol) / totVol;
      krinolog << "  Correction = " << x
            << ", Current volume = " << negVol
            << ", Target volume = " << startingNegVol
            << ", Relative Error = " << std::abs(relativeError)
            << stk::diag::dendl;
      const double derivative = area/totVol;
      return std::make_pair(relativeError, derivative);
    };
  return volume_error_function_with_derivative;
}

double
LevelSet::find_redistance_correction( const double start_area,
                                      const double start_neg_vol,
                                      const double start_pos_vol,
                                      const int max_iterations,
                                      const double tol )
{ /* %TRACE% */  /* %TRACE% */
  auto volume_error_function_with_derivative = build_volume_error_function_with_derivative(*this, start_neg_vol, start_pos_vol);
  const auto result = find_root_newton_raphson(volume_error_function_with_derivative, 0., max_iterations, tol);

  if (!result.first)
  {
    stk::RuntimeWarningAdHoc() << "\nConstrained renormalization failed to converge to root within "
      << max_iterations << " iterations. Continuing with correction " << result.second << "\n";
  }
  return result.second;
}

void
LevelSet::redistance()
{
  redistance(my_meta.universal_part());
}

void
LevelSet::redistance(const stk::mesh::Selector & selector)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::redistance(void)"); /* %TRACE% */

  ThrowErrorMsgIf(!my_time_of_arrival_element_speed_field_name.empty(), "Redistancing a time-of-arrival field will corrupt it.");

  if (get_isovar_field().valid()) get_isovar_field().field().sync_to_host();
  if (get_distance_field().valid()) get_distance_field().field().sync_to_host();
  if (get_old_distance_field().valid()) get_old_distance_field().field().sync_to_host();

  if (get_isovar_field().valid()) get_isovar_field().field().modify_on_host();
  if (get_distance_field().valid()) get_distance_field().field().modify_on_host();
  if (get_old_distance_field().valid()) get_old_distance_field().field().modify_on_host();

  if (FAST_MARCHING == my_redistance_method)
  {
    fast_marching_redistance(selector);
    return;
  }
  ThrowRequire(CLOSEST_POINT == my_redistance_method);

  krinolog << "Redistancing the level set field..." << stk::diag::dendl;

  // our starting point is a nodal variable (like distance or temperature)
  // that needs to be contoured to form the surface
  // after forming the surface, the nodal distance needs to be calculated
  // the newly formed surface should be remain in the vector facets
  build_facets_locally(selector);

  // debugging
  if (krinolog.shouldPrint(LOG_FACETS))
    {
      facets_exoii();
    }

  // swap these facets into facet_old to take advantage of routines
  // that are expecting the facets there
  facets->swap( *facets_old );

  // get non-local facets such that we have copies of all "old" facets
  // within the range of this proc's nodes
  prepare_to_compute_distance( 0., selector );

  // compute nodal distances with semi-lagrangian step
  compute_distance_semilagrangian( 0., selector );

  // swap so that the facets that were formed remain in the vector facets
  facets->swap( *facets_old );

}

double
LevelSet::get_time_of_arrival_speed(stk::mesh::Entity elem, ParallelErrorMessage& err) const
{
  double speed = 1.0;
  if (myTimeOfArrivalBlockSpeeds.empty())
  {
    if (myTimeOfArrivalElementSpeedField.valid())
    {
      const double * speedFieldData = field_data<double>(myTimeOfArrivalElementSpeedField, elem);
      if (!speedFieldData) err << "Missing element speed on element " << mesh().identifier(elem) << "\n";
      else speed = *speedFieldData;

      if (speed <= 0.0)
      {
        err << "Non positive-definite speed " << speed << " found on element " << mesh().identifier(elem) << "\n";
      }
    }
  }
  else
  {
    const stk::mesh::Part & elemPart = find_element_part(mesh(), elem);
    speed = myTimeOfArrivalBlockSpeeds[elemPart.mesh_meta_data_ordinal()];
    ThrowAssert(speed >= 0.0); // Negative speeds should have already been caught and generated error.
    if (speed == 0.0)
      err << "Speed not specified for block " << elemPart.name() << "\n";
  }

  return speed;
}

void
LevelSet::fast_marching_redistance(const stk::mesh::Selector & selector, const bool compute_time_of_arrival)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::fast_marching_redistance(const stk::mesh::Selector & selector)"); /* %TRACE% */

  // Unlike redistance() this method provides an approximate (not exact) distance to the isosurface.
  // On the other hand, this method provide solve the Eikonal equation on the mesh.  This is slightly
  // different than a pure distance function because it provides the distance through domain.  It can't see
  // through walls like the redistance() method does.

  if (compute_time_of_arrival) krinolog << "Initializing the level set field to be the time-of-arrival using a fast marching method..." << stk::diag::dendl;
  else krinolog << "Redistancing the level set field using a fast marching method..." << stk::diag::dendl;
  Fast_Marching fm(*this, selector, get_timer());
  fm.redistance();
}

//--------------------------------------------------------------------------------

void
LevelSet::set_distance( const double & in_distance ) const
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::set_distance( const double & in_distance ) const"); /* %TRACE% */

  stk::mesh::field_fill(in_distance, get_distance_field());
}

void
LevelSet::increment_distance( const double increment, const bool enforce_sign ) const
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::increment_distance( const double & increment ) const"); /* %TRACE% */
  //
  // increment the distance everywhere by the given value
  //

  const FieldRef dField = get_distance_field();

  const stk::mesh::Selector active_field_selector = aux_meta().active_not_ghost_selector() & stk::mesh::selectField(dField);
  stk::mesh::BucketVector const& buckets = mesh().get_buckets( stk::topology::NODE_RANK, active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const size_t length = b.size();

    double *d = field_data<double>(dField , b);

    for (size_t i = 0; i < length; ++i)
    {
      const double change = (enforce_sign && sign_change(d[i],d[i]+increment)) ? -0.5*d[i] : increment;
      d[i] += change;
    }
  }  // end bucket loop
}

void
LevelSet::scale_distance( const double scale) const
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::scale_distance( const double & scale ) const"); /* %TRACE% */
  //
  // increment the distance everywhere by the given value
  //

  const FieldRef dField = get_distance_field();

  const stk::mesh::Selector active_field_selector = aux_meta().active_not_ghost_selector() & stk::mesh::selectField(dField);
  stk::mesh::BucketVector const& buckets = mesh().get_buckets( stk::topology::NODE_RANK, active_field_selector);

  stk::mesh::BucketVector::const_iterator ib = buckets.begin();
  stk::mesh::BucketVector::const_iterator ib_end = buckets.end();

  // iterate nodes, by buckets, and set distance
  for ( ; ib != ib_end ; ++ib )
  {
    const stk::mesh::Bucket & b = **ib;

    const size_t length = b.size();

    double *d = field_data<double>(dField , b);

    for (size_t i = 0; i < length; ++i)
    {
      d[i] *= scale;
    }
  }  // end bucket loop
}

void
LevelSet::negate_distance() const
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::increment_distance( const double & increment ) const"); /* %TRACE% */

  const FieldRef dRef = get_distance_field();
  stk::mesh::field_scale(-1.0, dRef);
}

//-----------------------------------------------------------------------------------

void
LevelSet::compute_continuous_gradient() const
{ /* %TRACE[ON]% */ /* %TRACE% */
  //
  // Compute a mass-lumped continuous distance gradient
  //

  const int dim = mesh().mesh_meta_data().spatial_dimension();

  const FieldRef contGradRef = aux_meta().get_field(stk::topology::NODE_RANK, "CONT_GRAD");
  const FieldRef nodeAreaRef = aux_meta().get_field(stk::topology::NODE_RANK, "NODE_AREA");


  // initialize
  stk::mesh::field_fill(0.0, contGradRef);
  stk::mesh::field_fill(0.0, nodeAreaRef);

  // intg_pt_locations and intg_weights are just wrappers for framework data
  // determinants and grad_distance are arrays containing data that are resized
  // only as needed
  sierra::Array<const double,DIM,NINT> intg_pt_locations;
  sierra::Array<const double,NINT> intg_weights;
  sierra::ArrayContainer<double,NINT> determinants;
  sierra::ArrayContainer<double,DIM,NINT> grad_distance;

  // **************************
  // Not-covered-in-nightly RWH
  // **************************

  const FieldRef xField = get_coordinates_field();
  const FieldRef isoField = get_isovar_field();

  stk::mesh::Selector active_field_selector = stk::mesh::selectField(contGradRef) & aux_meta().active_not_ghost_selector();
  std::vector< stk::mesh::Entity> objs;
  stk::mesh::get_selected_entities( active_field_selector, mesh().buckets( stk::topology::ELEMENT_RANK ), objs );

  for ( auto && elem : objs )
  {
    // create element
    ContourElement ls_elem( mesh(), elem, xField, isoField );

    ls_elem.std_intg_pts( intg_pt_locations, intg_weights, determinants,
                          ls_elem.dist_master_elem() );
    ls_elem.compute_distance_gradient( intg_pt_locations, grad_distance );

    const unsigned num_intg_pts = intg_pt_locations.dimension<NINT>();
    const int npe = ls_elem.dist_topology().num_nodes();

    const double * shape_fcn_ptr = ls_elem.dist_master_elem().shape_fcn();
    const sierra::Array<const double,NPE_VAR,NINT> shape_fcn(shape_fcn_ptr,npe,num_intg_pts);

    const stk::mesh::Entity* elem_nodes = mesh().begin_nodes(elem);

    for ( int i = 0; i < npe; ++i )
    {
      double * cont_grad = field_data<double>(contGradRef, elem_nodes[i]);
      double * area_ptr = field_data<double>(nodeAreaRef, elem_nodes[i]);
      double & area = *area_ptr;

      for (unsigned ip = 0; ip < num_intg_pts; ++ip )
      {
        const double NdV = shape_fcn(i,ip) * intg_weights(ip) * determinants(ip);
        area += NdV;

        for ( int d = 0; d < dim; d++ )
        {
          cont_grad[d] += grad_distance(d,ip) * NdV;
          //krinolog << "grad_distance(" << d << "," << ip << ") = " << grad_distance(d,ip) << stk::diag::dendl;
        }
      }
    }
  }

  // iterate nodes, by buckets, and calculate average gradient

  stk::mesh::BucketVector const& buckets = mesh().get_buckets(stk::topology::NODE_RANK, active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const size_t length = b.size();

    double * cont_grad = field_data<double>(contGradRef, b);
    double * area = field_data<double>(nodeAreaRef, b);

    for (size_t i = 0; i < length; ++i)
    {
      for ( int d = 0; d < dim; d++ )
      {
        cont_grad[dim*i+d] = cont_grad[dim*i+d] / area[i];
      }
    }
  }
}

//-----------------------------------------------------------------------------------

void
LevelSet::compute_nodal_bbox( const stk::mesh::Selector & selector,
    BoundingBox & node_bbox,
    const Vector3d & displacement ) const
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::compute_nodal_bbox( BoundingBox & node_bboxes, const double & deltaTime ) const"); /* %TRACE% */

  // find the local nodal bounding box

  const FieldRef xField = get_coordinates_field();
  const FieldRef dField = get_distance_field();

  const stk::mesh::Selector active_field_selector = aux_meta().active_not_ghost_selector() & selector & stk::mesh::selectField(dField);
  stk::mesh::BucketVector const& buckets = mesh().get_buckets( stk::topology::NODE_RANK, active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const size_t length = b.size();

    double *x = field_data<double>(xField, b);

    for (size_t i = 0; i < length; ++i)
    {

      Vector3d x_bw(Vector3d::ZERO);
      for ( unsigned dim = 0; dim < spatial_dimension; ++dim )
      {
        int index = i*spatial_dimension+dim;
        x_bw[dim] = x[index] - displacement[dim];
      }

      // incrementally size bounding box
      node_bbox.accommodate( x_bw );
    }
  }  // end bucket loop
}

//-----------------------------------------------------------------------------------

void
LevelSet::prepare_to_compute_distance( const double & deltaTime, const stk::mesh::Selector & selector )
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::get_nonlocal_facets( const double & deltaTime )"); /* %TRACE% */

  // Get all of the facets that are within our processor's bounding box
  // To do this,  see if any local facets lie in the nodal bounding box
  // of another proc. if so, send them a copy of those facets

  // First, find the bounding box for each proc that contains all of the
  // nodes on that proc plus the narrow_band size

  BoundingBox node_bbox;
  const Vector3d displacement = deltaTime * get_extension_velocity();
  compute_nodal_bbox( selector, node_bbox, displacement );

  facets_old->prepare_to_compute(node_bbox, my_narrow_band_size);
}

//-----------------------------------------------------------------------------------

void
LevelSet::compute_distance_semilagrangian( const double & deltaTime, const stk::mesh::Selector & selector )
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::compute_distance_semilagrangian( const double & deltaTime ) const"); /* %TRACE% */

  const double h_avg = compute_average_edge_length();

  const FieldRef xField = get_coordinates_field();
  const FieldRef dField = get_distance_field();
  const Vector3d extv = get_extension_velocity();

  const stk::mesh::Selector active_field_selector = aux_meta().active_not_ghost_selector() & selector & stk::mesh::selectField(dField);
  stk::mesh::BucketVector const& buckets = mesh().get_buckets(stk::topology::NODE_RANK, active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const size_t length = b.size();

    double *d = field_data<double>( dField , b);
    double *x = field_data<double>( xField , b);

    // Handle special case of ( deltaTime == 0. ) so that we
    // can handle situation where velocity is not defined at all
    // nodes where distance is defined.  (This is currently a requirement
    // for regular semilagrangian advancement).
    if ( deltaTime == 0. )
    {
      for (size_t i = 0; i < length; ++i)
      {
        Vector3d x_node(Vector3d::ZERO);
        for ( unsigned dim = 0; dim < spatial_dimension; ++dim )
        {
          int index = i*spatial_dimension+dim;
          x_node[dim] = x[index];
        }

        int previous_sign = LevelSet::sign(d[i]);
        // If this is too large, then sharp edges can propagate incorrect signs
        // (even through walls, etc).
        // If this is too small, then a phase can't disappear because the sign
        // preservation will prevent it even if the subelement contouring process
        // neglects it.  So this should be slightly larger than the tolerance in
        // compute_subelement_decomposition.
        bool enforce_sign = (std::abs(d[i]) > 5.e-4*h_avg);

        d[i] = distance( x_node, previous_sign, enforce_sign );
      }
    }
    else
    {
      for (size_t i = 0; i < length; ++i)
      {
        Vector3d x_bw(Vector3d::ZERO);
        for ( unsigned dim = 0; dim < spatial_dimension; ++dim )
        {
          int index = i*spatial_dimension+dim;
          x_bw[dim] = x[index] - extv[dim] * deltaTime;
        }

        int previous_sign = LevelSet::sign(d[i]);

        d[i] = distance( x_bw, previous_sign, false );
      }
    }
  }  // end bucket loop
}

//-----------------------------------------------------------------------------------

void
LevelSet::compute_distance( stk::mesh::Entity n,
			    const double & deltaTime ) const
{ /* %TRACE% */  /* %TRACE% */

  // use the facet cell array to compute the distance to a node n

  const FieldRef xField = get_coordinates_field();
  const FieldRef dField = get_distance_field();
  const Vector3d extv = get_extension_velocity();

  double *x = field_data<double>( xField , n);
  double *d = field_data<double>( dField , n);

  // Handle special case of ( deltaTime == 0. ) so that we
  // can handle situation where velocity is not defined at all
  // nodes where distance is defined.  (This is currently a requirement
  // for regular semilagrangian advancement).
  if ( deltaTime == 0. )
  {
    Vector3d x_node(Vector3d::ZERO);
    for ( unsigned dim = 0; dim < spatial_dimension; ++dim )
    {
      x_node[dim] = x[dim];
    }

    int previous_sign = LevelSet::sign(*d);
    *d = distance( x_node, previous_sign, true );
  }
  else
  {
    Vector3d x_bw(Vector3d::ZERO);
    for ( unsigned dim = 0; dim < spatial_dimension; ++dim )
    {
      x_bw[dim] = x[dim] - extv[dim] * deltaTime;
    }

    int previous_sign = LevelSet::sign(*d);
    *d = distance( x_bw, previous_sign, false );
  }
}

//-----------------------------------------------------------------------------------

double
LevelSet::distance( const Vector3d & x,
		    const int previous_sign,
		    const bool enforce_sign ) const
{ /* %TRACE% */  /* %TRACE% */

  if (enforce_sign)
  {
    return previous_sign * facets_old->point_unsigned_distance(x, my_narrow_band_size, my_narrow_band_size);
  }
  return facets_old->truncated_point_signed_distance(x, my_narrow_band_size, previous_sign*my_narrow_band_size);
}

//-----------------------------------------------------------------------------------

void
LevelSet::snap_to_mesh() const
{ /* %TRACE[ON]% */ Trace trace__("LevelSet::snap_to_mesh()"); /* %TRACE% */
  const double tol = 1.0e-2;

  // Remove sliver subelements by setting nodal values near zero to zero exactly.
  // This should probably be an edge-based check.  But this poses a problem for higher order
  // elements, which we are going to decompose into lower order elements.  So we make it
  // simpler by compare ALL pairs of nodes within the element.  If the crossing between any
  // pair of nodes is is within a relative distance of tol, the distance at the nearest node
  // is set to zero.

  // This seems like a great way to consistently handle degeneracies.
  // One problem, however, is this only handles the near zero's on the original elements.
  // We will generate others as we decompose into non-conformal subelements.  This won't fix
  // those degenerate situations.

  std::vector<double> dist;

  const FieldRef dField = get_distance_field();

  const stk::mesh::Selector active_field_selector = stk::mesh::selectField(dField) & aux_meta().active_locally_owned_selector();
  stk::mesh::BucketVector const& buckets = mesh().get_buckets(stk::topology::ELEMENT_RANK, active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const stk::topology dist_topology = MasterElementDeterminer::get_field_topology(b, dField);
    const int npe_dist = dist_topology.num_nodes();
    dist.resize( npe_dist );

    const size_t length = b.size();
    for (size_t i = 0; i < length; ++i)
    {
      stk::mesh::Entity elem = b[i];

      if (!elem_on_interface(elem)) continue;

      const stk::mesh::Entity* elem_nodes = mesh().begin_nodes(elem);

      for ( int n = 0; n < npe_dist; ++n )
      {
        dist[n] = *field_data<double>(dField, elem_nodes[n]);
      }

      for ( int n = 0; n < npe_dist; ++n )
      {
        for ( int np = n+1; np < npe_dist; ++np )
        {
          if (sign_change(dist[n],dist[np]))
          {
            const double d0 = std::fabs(dist[n]);
            const double d1 = std::fabs(dist[np]);
            if (d0 < d1)
            {
              if (d0 / (d0+d1) < tol)
              {
                *field_data<double>(dField, elem_nodes[n]) = 0.0;
              }
            }
            else
            {
              if (d1 / (d1+d0) < tol)
              {
                *field_data<double>(dField, elem_nodes[np]) = 0.0;
              }
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------------

bool
LevelSet::remove_wall_features() const
{ /* %TRACE[ON]% */ Trace trace__("LevelSet::remove_wall_features()"); /* %TRACE% */

  std::vector<double> dist;
  std::vector<double*> coords;
  int small_feature_removed = false;

  const FieldRef dField = get_distance_field();
  const FieldRef coordinates_field = get_coordinates_field();

  stk::mesh::Selector surface_selector = selectUnion(my_surface_parts);

  const stk::mesh::Selector active_field_selector = stk::mesh::selectField(dField) & aux_meta().active_locally_owned_selector() & surface_selector;
  stk::mesh::BucketVector const& buckets = mesh().get_buckets(meta().side_rank(), active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const stk::topology dist_topology = MasterElementDeterminer::get_field_topology(b, dField);
    const int npe_dist = dist_topology.num_nodes();
    dist.resize( npe_dist );
    coords.resize( npe_dist );

    const size_t length = b.size();
    for (size_t i = 0; i < length; ++i)
    {
      stk::mesh::Entity side = b[i];

      const stk::mesh::Entity* side_nodes = mesh().begin_nodes(side);

      for ( int n = 0; n < npe_dist; ++n )
      {
        dist[n] = *field_data<double>(dField, side_nodes[n]);

        double *coord_data = field_data<double>(coordinates_field,side_nodes[n]);
        coords[n] = coord_data;
      }

      for ( int n = 0; n < npe_dist; ++n )
      {
        const stk::mesh::Entity *elems = mesh().begin_elements(side); //get the element
        const int num_elems = mesh().num_elements(side);

        for (int l = 0; l < num_elems; ++l )
        {
          const stk::mesh::Entity elem = elems[l];

          //make sure we have the right element (active) and has the distance function
          //at every point in the element (to take gradient)

          if(!aux_meta().active_locally_owned_selector()(mesh().bucket(elem))) continue;
          if(!elem_has_field_data(dField, elem)) continue;

          //const double elem_len = calc_elem_len_normal(elem, side, coordinates_field);

          if(std::fabs(dist[n]) > my_max_feature_size) continue;

          ContourElement ls_elem( mesh(), elem, coordinates_field, dField );
          const Vector3d p_coords(1/3., 1/3., 1/3.);
          const Vector3d grad_dist_vec = ls_elem.distance_gradient(p_coords);

          Vector3d face_normal;

          //assume linear tet or tri elements!
          if(spatial_dimension == 2)
          {
            face_normal = {-coords[1][1]+coords[0][1], coords[1][0]-coords[0][0], 0}; //simple 90 deg rotation
            const stk::mesh::Entity* elem_nodes = mesh().begin_nodes(elem);

            int num_elem_nodes = mesh().num_nodes(elem);

            for (int j = 0; j < num_elem_nodes; ++j) //find node in elem not part of side part, confirm normal points into elem
            {
              if(elem_nodes[j] != side_nodes[0] && elem_nodes[j] != side_nodes[spatial_dimension-1])
              {
                double *coord = field_data<double>(coordinates_field,elem_nodes[j]);
                const Vector3d vec_check(coord[0]-coords[0][0], coord[1]-coords[0][1], 0);

                if(Dot(face_normal, vec_check) < 0)
                {
                  face_normal *= -1;
                  break;
                }
              }
            }
          }
          else
          {
            const Vector3d x1 (coords[1][0]-coords[0][0], coords[1][1]-coords[0][1],coords[1][2]-coords[0][2]);
            const Vector3d x2 (coords[2][0]-coords[0][0], coords[2][1]-coords[0][1],coords[2][2]-coords[0][2]);
            face_normal = -1.0*Cross(x1,x2);
          }

          face_normal.unitize();
          const double eta = -dist[n]/(Dot(grad_dist_vec, face_normal));

          if(eta < my_max_feature_size && eta > 0)
          {
            *field_data<double>(dField, side_nodes[n]) = -1.0 * dist[n];
            small_feature_removed = true;
          }
        }
      }
    }
  }

  int global_small_feature_removed = false;
  stk::all_reduce_sum(mesh().parallel(), &small_feature_removed, &global_small_feature_removed, 1);

  if (global_small_feature_removed)
  {
    stk::mesh::communicate_field_data(mesh(), {&dField.field()});
    return true;
  }

  return false;
}
//--------------------------------------------------------------------------------

bool
LevelSet::simple_remove_wall_features() const
{ /* %TRACE[ON]% */ Trace trace__("LevelSet::remove_wall_features()"); /* %TRACE% */

  int small_feature_removed = false;
  const FieldRef dField = get_distance_field();

  stk::mesh::Selector surface_selector = selectUnion(my_surface_parts);

  const stk::mesh::Selector active_field_selector = stk::mesh::selectField(dField) & aux_meta().active_part() & surface_selector;
  stk::mesh::BucketVector const& buckets = mesh().get_buckets(stk::topology::NODE_RANK, active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;
    const size_t length = b.size();
    for (size_t i = 0; i < length; ++i)
    {
      stk::mesh::Entity node = b[i];
      double & dist = *field_data<double>(dField, node);

      if(std::fabs(dist) < my_max_feature_size)
      {
        dist*=-1;
        small_feature_removed = true;
      }
    }
  }

  return small_feature_removed;
}
//--------------------------------------------------------------------------------

void
LevelSet::set_surface_parts_vector()
{
  Phase_Support my_phase_support = Phase_Support::get(meta());

  std::vector<stk::mesh::Part *> conformal_parts =  my_phase_support.get_conformal_parts();

  krinolog << "Removing small features less than " << my_max_feature_size << " on surfaces: ";
  for (auto && mypart : conformal_parts)
  {
    if(mypart->primary_entity_rank() == meta().side_rank() && !my_phase_support.is_interface(mypart))
    {
      my_surface_parts.push_back(mypart);
      krinolog << mypart->name() << " ";
    }
  }

  krinolog << stk::diag::dendl;
}

//--------------------------------------------------------------------------------

bool
LevelSet::elem_has_field_data(const FieldRef &myField, const stk::mesh::Entity &elem) const
{
  const stk::mesh::Entity* elem_nodes = mesh().begin_nodes(elem);
  int num_nodes = mesh().num_nodes(elem);

  for (int k = 0; k < num_nodes; ++k)
  {
    if(!has_field_data(myField, elem_nodes[k])) return false;
  }

 return true;
}
//--------------------------------------------------------------------------------

bool
LevelSet::elem_on_interface(stk::mesh::Entity e) const
{ /* %TRACE% */  /* %TRACE% */

  const FieldRef isoField = get_isovar_field();

  const unsigned nnodes = mesh().num_nodes(e);
  ThrowAssert( 0 < nnodes );
  const stk::mesh::Entity* nodes = mesh().begin_nodes(e);

  // guilty till proven innocent here
  bool on_interface = false;
  bool have_crossing = false;

  double *d = field_data<double>(isoField, nodes[0]);
  double first_value = *d - my_threshold;
  bool inside_narrow_band = fabs(first_value) < 0.5*my_narrow_band_size;

  for (unsigned i = 1; i < nnodes; ++i )
    {
      d = field_data<double>(isoField, nodes[i]);
      if ( NULL == d ) continue;  // account for lower order interpolation
      double value = *d - my_threshold;

      have_crossing |= sign_change(value, first_value);
      inside_narrow_band |= (fabs(value) < 0.5*my_narrow_band_size);
    }

  // It is the user's job to make sure that the narrow band is sufficiently
  // large that we don't to test if this crossing is within the narrow band
  if ( have_crossing )
    on_interface = true;

  return on_interface;
}

//--------------------------------------------------------------------------------
void
LevelSet::compute_sizes( double & area, double & neg_vol, double & pos_vol, const double iso_val )
{ /* %TRACE[ON]% */ /* %TRACE% */

  // initialize
  area = 0.0;
  neg_vol = 0.0;
  pos_vol = 0.0;

  const double h_avg = compute_average_edge_length();

  sierra::ArrayContainer<double,DIM,NINT> intg_pt_locations;
  sierra::ArrayContainer<double,NINT> intg_weights;
  sierra::ArrayContainer<double,NINT> determinants;

  const FieldRef xField = get_coordinates_field();
  const FieldRef isoField = get_isovar_field();
  //const Vector3d extv = get_extension_velocity();

  stk::mesh::Selector active_field_selector = stk::mesh::selectField(isoField) & aux_meta().active_locally_owned_selector();
  std::vector< stk::mesh::Entity> objs;
  stk::mesh::get_selected_entities(active_field_selector, mesh().buckets( stk::topology::ELEMENT_RANK ), objs);

  for ( auto && elem : objs )
  {
    // create element that is decomposed into subelements
    ContourElement ls_elem( mesh(), elem, xField, isoField, iso_val );

    ls_elem.compute_subelement_decomposition(h_avg);

    // get integration point locations, weights, and determinants for surface
    int num_intg_pts = ls_elem.gather_intg_pts( 0,		// interface
                                                intg_pt_locations,// gauss pt locations
                                                intg_weights,	// integration wts at gauss pts
                                                determinants,	// determinant at gauss pts
                                                true );		// map_to_real_coords
    for ( int ip = 0; ip < num_intg_pts; ++ip )
      {
        area += intg_weights(ip) * determinants(ip);
      }

    // get integration point locations, weights, and determinants for negative volume
    num_intg_pts = ls_elem.gather_intg_pts( -1,			// negative side of interface
                                            intg_pt_locations,	// gauss pt locations
                                            intg_weights,		// integration wts at gauss pts
                                            determinants,		// determinant at gauss pts
                                            true );		// map_to_real_coords
    for ( int ip = 0; ip < num_intg_pts; ++ip )
      {
        neg_vol += intg_weights(ip) * determinants(ip);
      }

    // get integration point locations, weights, and determinants for positive volume
    num_intg_pts = ls_elem.gather_intg_pts( 1,			// positive side of interface
                                            intg_pt_locations,	// gauss pt locations
                                            intg_weights,		// integration wts at gauss pts
                                            determinants,		// determinant at gauss pts
                                            true );		// map_to_real_coords
    for ( int ip = 0; ip < num_intg_pts; ++ip )
      {
        pos_vol += intg_weights(ip) * determinants(ip);
      }
  }

  // communicate global sums
  const int vec_length = 3;
  std::vector <double> local_sum( vec_length );
  std::vector <double> global_sum( vec_length );
  local_sum[0] = area;
  local_sum[1] = neg_vol;
  local_sum[2] = pos_vol;

  stk::all_reduce_sum(mesh().parallel(), &local_sum[0], &global_sum[0], vec_length);

  area    = global_sum[0];
  neg_vol = global_sum[1];
  pos_vol = global_sum[2];

}
//--------------------------------------------------------------------------------
double
LevelSet::gradient_magnitude_error(void)
{ /* %TRACE[ON]% */ /* %TRACE% */

  double area = 0., sum_L2 = 0., global_L2 = 0., local_Loo = 0., global_Loo = 0.;

  sierra::ArrayContainer<double,DIM,NINT> intg_pt_locations;
  sierra::ArrayContainer<double,NINT> intg_weights;
  sierra::ArrayContainer<double,NINT> determinants;
  sierra::ArrayContainer<double,DIM,NINT> grad_dist;

  const double h_avg = compute_average_edge_length();

  const FieldRef xField = get_coordinates_field();
  const FieldRef isoField = get_isovar_field();
  //const Vector3d extv = get_extension_velocity();

  stk::mesh::Selector active_field_selector = stk::mesh::selectField(isoField) & aux_meta().active_locally_owned_selector();
  std::vector< stk::mesh::Entity> objs;
  stk::mesh::get_selected_entities( active_field_selector, mesh().buckets( stk::topology::ELEMENT_RANK ), objs );

  for ( auto && elem : objs )
  {
    // create element that is decomposed into subelements
    ContourElement ls_elem( mesh(), elem, xField, isoField );

    ls_elem.compute_subelement_decomposition(h_avg);

    // get integration point locations, weights, and determinants for surface
    int num_intg_pts = ls_elem.gather_intg_pts( 0,		// interface
                                                intg_pt_locations,// gauss pt locations
                                                intg_weights,	// integration wts at gauss pts
                                                determinants,	// determinant at gauss pts
                                                true );		// map_to_real_coords

    ls_elem.compute_distance_gradient( intg_pt_locations, grad_dist );
    for ( int ip = 0; ip < num_intg_pts; ++ip )
    {
      double mag_grad_phi = 0.;
      for ( unsigned dim = 0; dim < spatial_dimension; ++dim ) mag_grad_phi += grad_dist(dim,ip) * grad_dist(dim,ip);
      mag_grad_phi = sqrt(mag_grad_phi);

      sum_L2 += (mag_grad_phi - 1.) * (mag_grad_phi - 1.) * intg_weights(ip) * determinants(ip);
      area += intg_weights(ip) * determinants(ip);

      if ( fabs(mag_grad_phi - 1.) > local_Loo ) local_Loo = fabs(mag_grad_phi - 1.);
    }
  }

  // communicate global norms
  const int vec_length = 2;
  std::vector <double> local_sum( vec_length );
  std::vector <double> global_sum( vec_length );
  local_sum[0] = sum_L2;
  local_sum[1] = area;

  stk::all_reduce_sum(mesh().parallel(), &local_sum[0], &global_sum[0], vec_length);
  stk::all_reduce_max(mesh().parallel(), &local_Loo, &global_Loo, 1);

  if ( global_sum[1] > 0. )
    {
      global_L2 = global_sum[0] / global_sum[1];
    }

  krinolog << "Gradient norm error for " << name() << ": L2 = " << global_L2 << ", Loo = " << global_Loo << stk::diag::dendl;

  // L2 is the standard now, maybe Loo would be better?
  return global_L2;
}
//--------------------------------------------------------------------------------
double
LevelSet::compute_average_edge_length() const
{ /* %TRACE[ON]% */ /* %TRACE% */

  double sum_avg_edge_lengths = 0.0;
  int num_elems = 0;

  const FieldRef xField = get_coordinates_field();
  const FieldRef isoField = get_isovar_field();

  stk::mesh::Selector active_field_selector = stk::mesh::selectField(isoField) & aux_meta().active_locally_owned_selector();
  std::vector< stk::mesh::Entity> objs;
  stk::mesh::get_selected_entities( active_field_selector, mesh().buckets( stk::topology::ELEMENT_RANK ), objs );

  for ( auto && elem : objs )
  {
    ++num_elems;

    // create element
    ContourElement ls_elem( mesh(), elem, xField, isoField );

    sum_avg_edge_lengths += ls_elem.average_edge_length();
  }

  // communicate global sums
  const int vec_length = 2;
  std::vector <double> local_sum( vec_length );
  std::vector <double> global_sum( vec_length );
  local_sum[0] = sum_avg_edge_lengths;
  local_sum[1] = 1.0*num_elems;

  stk::all_reduce_sum(mesh().parallel(), &local_sum[0], &global_sum[0], vec_length);

  const double h_avg = ( global_sum[1] != 0.0 ) ? global_sum[0]/global_sum[1] : 0.0;

  return h_avg;
}

//--------------------------------------------------------------------------------

void
LevelSet::build_facets_locally(const stk::mesh::Selector & selector)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::build_facets_locally(void)"); /* %TRACE% */

  // clear vector of facets
  facets->clear();

  const double h_avg = compute_average_edge_length();

  const FieldRef xField = get_coordinates_field();
  const FieldRef isoField = get_isovar_field();

  stk::mesh::Selector active_field_selector = selector & stk::mesh::selectField(isoField) & aux_meta().active_locally_owned_selector();
  std::vector< stk::mesh::Entity> objs;
  stk::mesh::get_selected_entities( active_field_selector, mesh().buckets( stk::topology::ELEMENT_RANK ), objs );

  for ( auto && elem : objs )
  {
    // create element that is decomposed into subelements
    ContourElement ls_elem( mesh(), elem, xField, isoField );
    ls_elem.compute_subelement_decomposition(h_avg);

    ls_elem.build_subelement_facets( *facets );
  }
}

LevelSet &
LevelSet::build(
    stk::mesh::MetaData & in_meta,
    const std::string & ls_name,
    stk::diag::Timer & parent_timer )
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::build(stk::mesh::MetaData & in_meta, const std::string & ls_name, stk::diag::Timer & parent_timer)"); /* %TRACE% */
  LevelSetManager & region_ls = LevelSetManager::get(in_meta);

  ThrowRequire(!region_ls.has_levelSet(ls_name));
  LevelSet * ls = new LevelSet(in_meta, ls_name, parent_timer);
  region_ls.add(ls);
  return *ls;
}

//--------------------------------------------------------------------------------
LevelSet::LevelSet(
    stk::mesh::MetaData & in_meta,
    const std::string & in_name,
    stk::diag::Timer & parent_timer ) :
    my_meta(in_meta),
    my_aux_meta(AuxMetaData::get(in_meta)),
    my_identifier(LevelSet::get_identifier(in_name)),
    my_name(LevelSet::get_name(my_identifier)),
    my_parent_timer(parent_timer),
    my_timer("LevelSet", parent_timer),
    spatial_dimension(in_meta.spatial_dimension()),
    my_narrow_band_multiplier(0.0),
    my_narrow_band_size(0.0),
    my_max_feature_size(-1.0),
    my_ic_offset(0.0),
    my_ic_scale(1.0),
    my_perform_initial_redistance(false),
    my_keep_IC_surfaces(false),
    my_threshold(0.0),
    my_redistance_method(CLOSEST_POINT),
    epsilon(1.0e-16),
    trackIsoSurface(false),
    my_facetFileIndex(1),
    my_initial_neg_vol(0.0),
    my_needs_reinitialize_every_step(false)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::LevelSet(stk::mesh::MetaData & in_meta, const std::string & ls_name, stk::diag::Timer & parent_timer)"); /* %TRACE% */
  my_coordinates_field = get_current_coordinates(in_meta);

  // default names for distance and velocity
  // line commands are available for overriding these names
  const std::string distName = "D_" + name();
  set_distance_name(distName);
}

LevelSet::~LevelSet()
{
}

//-----------------------------------------------------------------------------------
void
LevelSet::gather_nodal_field(
  const stk::mesh::BulkData& stk_mesh,
  stk::mesh::Entity obj,
  const FieldRef & field,
  double * gathered_field_data )
{ /* %TRACE% */  /* %TRACE% */

  int ncomp_field = field.length();

  // Gather obj's nodal field into a single flat array
  // of dimension (ncomp_field,num_nodes)

  int j = 0;
  const unsigned num_nodes = stk_mesh.num_nodes(obj);
  const stk::mesh::Entity* nodes = stk_mesh.begin_nodes(obj);
  for (unsigned node_index=0; node_index<num_nodes; ++node_index)
  {
    stk::mesh::Entity node = nodes[node_index];
    double * var = field_data<double>(field, node);
    for ( int i = 0; i < ncomp_field; ++i )
    {
      gathered_field_data[j++] = var[i];
    }
  }

  ThrowAssert( (unsigned)j == ncomp_field * stk_mesh.num_nodes(obj));
}
//--------------------------------------------------------------------------------
LevelSetManager &
LevelSetManager::get(stk::mesh::MetaData & meta)
{
  LevelSetManager * mgr = const_cast<LevelSetManager *>(meta.get_attribute<LevelSetManager>());
  if (nullptr == mgr)
  {
    mgr = new LevelSetManager;
    meta.declare_attribute_with_delete<LevelSetManager>(mgr);
  }
  return *mgr;
}
//--------------------------------------------------------------------------------
LevelSetManager &
LevelSetManager::get(const stk::mesh::MetaData & meta)
{
  LevelSetManager * mgr = const_cast<LevelSetManager *>(meta.get_attribute<LevelSetManager>());
  ThrowRequireMsg(nullptr != mgr, "No LevelSetManager found for MetaData.");
  return *mgr;
}
//--------------------------------------------------------------------------------
bool LevelSetManager::has_levelSet(const std::string & ls_name) const
{
  for (auto&& ls : my_level_sets)
  {
    if (ls->name() == ls_name || ls->get_composite_name() == ls_name)
    {
      return true;
    }
  }
  return false;
}
//--------------------------------------------------------------------------------
std::string
print_sizes(const LevelSet & ls)
{ /* %TRACE[ON]% */ Trace trace__("krino::print_sizes(const LevelSet & levelSet)"); /* %TRACE% */

  //
  // find sizes of facets stored in vector of facet descriptions
  //

  const auto & ls_surfaces = ls.get_facets().get_facets();
  unsigned local_facet_num = ls_surfaces.size();
  unsigned global_facet_num = 0;
  stk::all_reduce_sum(ls.mesh().parallel(), &local_facet_num , &global_facet_num, 1);

  // iterate local facets and calc area

  double local_facets_totalArea = 0.0; // total surface area of interface from facets
  double local_facets_maxArea = -1.0; // area of largest facet on interface
  double local_facets_minArea = std::numeric_limits<double>::max();  // area of smallest facet on interface

  // loop over facets
  for ( auto&& surface : ls_surfaces )
  {
    double area = surface->facet_area();
    local_facets_totalArea += area;

    local_facets_maxArea = std::min(area, local_facets_maxArea);
    local_facets_minArea = std::max(area, local_facets_maxArea);
  }

  double global_facets_totalArea = 0.0;
  double global_facets_maxArea = 0.0;
  double global_facets_minArea = 0.0;

  stk::all_reduce_min(ls.mesh().parallel(), &local_facets_minArea, &global_facets_minArea, 1);
  stk::all_reduce_max(ls.mesh().parallel(), &local_facets_maxArea, &global_facets_maxArea, 1);
  stk::all_reduce_sum(ls.mesh().parallel(), &local_facets_totalArea, &global_facets_totalArea, 1);

  std::ostringstream out;

  // facet info
  out << "P" << ls.mesh().parallel_rank() << ": facets for level set '" << ls.name() << "' " << std::endl ;
  out << "\t " << "Global sizes: { " << global_facet_num << " facets on }" << std::endl;
  out << "\t " << "Local sizes: { " << local_facet_num << " facets on }" << std::endl;
  out << "\t Global areas: { Min = " << global_facets_minArea << ", Max = " << global_facets_maxArea
  << ", Total = " << global_facets_totalArea << " }" << std::endl;
  return out.str();
}
//--------------------------------------------------------------------------------

} // namespace krino
