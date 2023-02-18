// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Region.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_AdaptivityInterface.hpp>
#include <Akri_AdaptivityHelpers.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_MeshInputOptions.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_ResultsOutputOptions.hpp>
#include <Akri_Simulation.hpp>
#include <Ioss_Region.h>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <Akri_BoundingSurface.hpp>
#include <Akri_CDMesh_Refinement.hpp>
#include <Akri_Surface_Manager.hpp>
#include <Akri_RefinementInterface.hpp>
#include <Akri_RefinementSupport.hpp>

namespace krino{

//--------------------------------------------------------------------------------

Region::Region(Simulation & owning_simulation, const std::string & regionName)
: my_simulation(owning_simulation),
  my_meta(nullptr),
  my_bulk(nullptr),
  my_name(regionName),
  my_timerRegion(std::string("Region ") + regionName, sierra::Diag::TIMER_REGION, my_simulation.get_timer()),
  my_timerInitialize("Initialize", my_timerRegion),
  my_timerExecute("Execute", my_timerRegion),
  my_timerMeshInput("Mesh input", my_timerRegion),
  my_timerMeshOutput("Mesh output", my_timerRegion),
  my_output_file_index(0),
  my_output_file_created(false)
{ /* %TRACE[ON]% */ Trace trace__("krino::Region::Region()"); /* %TRACE% */
  my_simulation.add_region(this);
  myIOBroker = std::make_unique<stk::io::StkMeshIoBroker>(stk::EnvData::parallel_comm());
  myIOBroker->use_simple_fields();

  std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
  entity_rank_names.push_back("FAMILY_TREE");
  stk_IO().set_rank_name_vector(entity_rank_names);

  my_results_options = std::make_unique<ResultsOutputOptions>();
}
//--------------------------------------------------------------------------------
Region::~Region()
{
}

static bool locally_has_64bit_ids_in_use_for_nodes_or_elements(const stk::mesh::BulkData & mesh)
{
  const uint64_t max32bitId = ~0U;
  for (auto entityRank : {stk::topology::NODE_RANK, stk::topology::ELEMENT_RANK})
    for (auto && bucketPtr : mesh.buckets(entityRank))
      for (auto && entity : *bucketPtr)
        if (mesh.identifier(entity) > max32bitId)
          return true;

  return false;
}

static bool cdfem_mesh_displacements_requested_in_results_fields(const std::string & cdfemMeshDisplacementsFieldName, const std::set<FieldName_OutputName_Pair> & resultsFields)
{
  for (auto && resultsField : resultsFields)
    if (stk::equal_case(resultsField.first, cdfemMeshDisplacementsFieldName))
      return true;
  return false;
}

//--------------------------------------------------------------------------------
void Region::commit()
{ /* %TRACE[ON]% */ Trace trace__("krino::Region::commit()"); /* %TRACE% */

  auto & meta = get_stk_mesh_meta_data();
  LevelSet::setup(meta);
  CDFEM_Support & cdfem_support = CDFEM_Support::get(meta);
  RefinementSupport & refinementSupport = RefinementSupport::get(meta);

  if (krino::CDFEM_Support::is_active(meta))
  {
    if (cdfem_mesh_displacements_requested_in_results_fields(CDFEM_Support::cdfem_mesh_displacements_field_name(), my_results_options->get_nodal_fields()))
      cdfem_support.register_cdfem_mesh_displacements_field();

    if (cdfem_support.get_num_initial_decomposition_cycles() > 1)
      cdfem_support.register_parent_node_ids_field(); // Needed to use piecewise linear location on previously cut edges

    cdfem_support.setup_fields();
  }

  auto & active_part = AuxMetaData::get(meta).active_part();
  stk::mesh::BulkData::AutomaticAuraOption auto_aura_option = stk::mesh::BulkData::NO_AUTO_AURA;

  if (refinementSupport.get_initial_refinement_levels() > 0 || refinementSupport.get_interface_maximum_refinement_level() > 0 ||
      (krino::CDFEM_Support::is_active(meta) && cdfem_support.get_post_cdfem_refinement_levels() > 0))
  {
    auto_aura_option = stk::mesh::BulkData::AUTO_AURA;
    RefinementInterface & refinement = create_refinement(meta, refinementSupport.get_use_percept(), my_timerExecute);
    refinementSupport.set_non_interface_conforming_refinement(refinement);
  }

  if (cdfem_support.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
  {
    auto_aura_option = stk::mesh::BulkData::AUTO_AURA;
  }

  if (krino::CDFEM_Support::is_active(meta))
  {
    cdfem_support.add_interpolation_field(cdfem_support.get_coords_field());
    const Surface_Manager & surfaceManager = Surface_Manager::get(meta);
    for(auto&& ls : surfaceManager.get_levelsets())
    {
      cdfem_support.add_interpolation_field(ls->get_distance_field());
    }
    cdfem_support.finalize_fields();
  }

  if (nullptr != my_generated_mesh.get())
  {
    set_generated_mesh_domain();
    my_generated_mesh->populate_mesh(stk::EnvData::parallel_comm(), auto_aura_option);
    my_bulk = &my_generated_mesh->bulk_data();
    stk_IO().set_bulk_data( *my_bulk );
    if (my_generated_mesh->has_flat_boundaries())
      my_generated_mesh->create_domain_sides();
  }
  else
  {
    std::shared_ptr<stk::mesh::BulkData> shared_bulk = stk::mesh::MeshBuilder(stk::EnvData::parallel_comm()).set_aura_option(auto_aura_option).create(std::shared_ptr<stk::mesh::MetaData>(&meta,[](auto ptrWeWontDelete){}));
    my_bulk = shared_bulk.get();
    my_bulk->mesh_meta_data().use_simple_fields();
    stk_IO().set_bulk_data( shared_bulk );
    stk_IO().populate_bulk_data();
  }

  if (AuxMetaData::get(get_stk_mesh_meta_data()).get_assert_32bit_flag())
  {
    const bool has_64bit_ids_in_use = stk::is_true_on_any_proc(my_bulk->parallel(), locally_has_64bit_ids_in_use_for_nodes_or_elements(*my_bulk));
    ThrowErrorMsgIf(has_64bit_ids_in_use, "Option use_32_bit ids is active, but input file uses 64 bit ids.");
  }
  else
  {
    my_bulk->set_large_ids_flag(true);
  }

  if (my_bulk->is_automatic_aura_on() || my_bulk->parallel_size() == 1)
  {
    // Used for element side connectivty checks
    stk::mesh::create_exposed_block_boundary_sides(*my_bulk, meta.universal_part(), {&AuxMetaData::get(meta).exposed_boundary_part()});
  }



  activate_all_entities(*my_bulk, active_part);

  LevelSet::post_commit_setup(meta);
}

namespace {

void zero_error_indicator(stk::mesh::BulkData & mesh)
{
  auto & meta = mesh.mesh_meta_data();
  const auto & refinementSupport = RefinementSupport::get(meta);
  const auto & aux_meta = AuxMetaData::get(meta);
  auto indicator_field =
      aux_meta.get_field(stk::topology::ELEMENT_RANK,
          refinementSupport.get_nonconformal_adapt_indicator_name());

  const auto & buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK,
      meta.locally_owned_part());
  for(auto && b_ptr : buckets)
  {
    auto * data = field_data<double>(indicator_field, *b_ptr);
    auto length = b_ptr->size();
    std::fill(data, data + length, 0.);
  }
}

static void refine_elements_near_interface(RefinementInterface & refinement,
    const RefinementSupport & refinementSupport,
    stk::mesh::BulkData & mesh,
    const int numRefinementSteps)
{
  std::function<void(int)> mark_elements_near_interface =
      [&mesh, &refinementSupport, numRefinementSteps]
       (int num_refinements)
      {
        LevelSet::initialize(mesh.mesh_meta_data());
        if (num_refinements < numRefinementSteps)
        {
          const std::unique_ptr<InterfaceGeometry> interfaceGeometry = create_interface_geometry(mesh.mesh_meta_data());
          krino::mark_possible_cut_elements_for_adaptivity(mesh,
              refinementSupport.get_non_interface_conforming_refinement(),
              *interfaceGeometry,
              refinementSupport,
              num_refinements);
        }
      };

  perform_multilevel_adaptivity(refinement, mesh, mark_elements_near_interface, refinementSupport.get_do_not_refine_or_unrefine_selector());
}

static void refine_interface_elements(RefinementInterface & refinement,
    const RefinementSupport & refinementSupport,
    stk::mesh::BulkData & mesh,
    const int numRefinementSteps)
{
  std::function<void(int)> marker_function =
      [&mesh, &refinementSupport, numRefinementSteps]
       (int num_refinements)
      {
        LevelSet::initialize(mesh.mesh_meta_data());
        if(num_refinements < numRefinementSteps)
        {
          const auto & auxMeta = krino::AuxMetaData::get(mesh.mesh_meta_data());
          const std::unique_ptr<InterfaceGeometry> interfaceGeometry = create_interface_geometry(mesh.mesh_meta_data());
          CDMesh::mark_interface_elements_for_adaptivity(mesh, auxMeta.get_current_coordinates(), refinementSupport, *interfaceGeometry, num_refinements);
        }
      };

  perform_multilevel_adaptivity(refinement, mesh, marker_function, refinementSupport.get_do_not_refine_or_unrefine_selector());
}

static void refine_based_on_indicator_field(RefinementInterface & refinement,
    const krino::RefinementSupport & refinementSupport,
    stk::mesh::BulkData & mesh,
    const int targetCount,
    const int numRefinementSteps)
{
  zero_error_indicator(mesh);
  const auto & indicatorFieldName = refinementSupport.get_nonconformal_adapt_indicator_name();

  std::function<void(int)> marker_function =
        [&refinement, &mesh, numRefinementSteps, targetCount, indicatorFieldName]
         (int num_refinements)
        {
          zero_error_indicator(mesh);
          LevelSet::initialize(mesh.mesh_meta_data());
          if(num_refinements < numRefinementSteps)
          {
            mark_based_on_indicator_field(mesh, refinement, indicatorFieldName, numRefinementSteps,  num_refinements, targetCount);
          }
        };

  perform_multilevel_adaptivity(refinement, mesh, marker_function, refinementSupport.get_do_not_refine_or_unrefine_selector());
}

void do_adaptive_refinement(const krino::RefinementSupport & refinementSupport, stk::mesh::BulkData & mesh)
{
  stk::diag::TimeBlock adapt_timer__(refinementSupport.get_timer());

  const auto targetCount = refinementSupport.get_nonconformal_adapt_target_count();
  auto & refinement = refinementSupport.get_non_interface_conforming_refinement();

  const int numRefinementLevels = refinementSupport.get_interface_maximum_refinement_level();

  if (targetCount > 0)
  {
    refine_based_on_indicator_field(refinement, refinementSupport, mesh, targetCount, numRefinementLevels);
  }
  else
  {
    const bool doNearbyRefinementBeforeInterfaceRefinement = refinementSupport.do_nearby_refinement_before_interface_refinement();
    if (doNearbyRefinementBeforeInterfaceRefinement)
    {
      const int numNearbyRefinementSteps = numRefinementLevels;
      refine_elements_near_interface(refinement, refinementSupport, mesh, numNearbyRefinementSteps);
    }
    const int numInterfaceRefinementSteps = doNearbyRefinementBeforeInterfaceRefinement ? numRefinementLevels : (2*numRefinementLevels); // Make sure refinement completes so that interfacial elements are fully refined
    refine_interface_elements(refinement, refinementSupport, mesh, numInterfaceRefinementSteps);
  }
}

void do_post_adapt_uniform_refinement(const Simulation & simulation, const RefinementSupport & refinementSupport, const AuxMetaData  & auxMeta, stk::mesh::BulkData & mesh)
{
  if ( refinementSupport.get_post_adapt_refinement_levels() > 0 )
  {
    if ( simulation.is_transient() )
    {
      krinolog << "Can not do post-adaptivity uniform mesh refinement for transient problems...";
    }
    else
    {
      krinolog << "Performing " << refinementSupport.get_post_adapt_refinement_levels() << " levels of post-adaptivity uniform mesh refinement..." << std::endl;

      // Doing adaptive refinement with a uniform marker is better than doing uniform refinement here because of how
      // the transition elements are handled.
      auto & refinement = refinementSupport.get_non_interface_conforming_refinement();
      const int num_levels = refinementSupport.get_post_adapt_refinement_levels();
      FieldRef marker_field = refinement.get_marker_field();

      std::function<void(int)> marker_function =
          [&mesh, marker_field, num_levels]
           (int num_refinements)
          {
            LevelSet::initialize(mesh.mesh_meta_data());
            if (num_refinements < num_levels)
            {
              stk::mesh::field_fill(1, marker_field);
            }
          };

      perform_multilevel_adaptivity(refinement, mesh, marker_function, refinementSupport.get_do_not_refine_or_unrefine_selector());
    }
  }
}

void do_post_cdfem_uniform_refinement(const Simulation & simulation, const CDFEM_Support & cdfemSupport, const RefinementSupport & refinementSupport, const AuxMetaData  & auxMeta, stk::mesh::BulkData & mesh)
{
  if ( cdfemSupport.get_post_cdfem_refinement_levels() > 0 )
  {
    const int num_levels = cdfemSupport.get_post_cdfem_refinement_levels();
    if ( simulation.is_transient() )
    {
      krinolog << "Can not do post-cdfem mesh refinement for transient problems...";
    }
    else
    {
      krinolog << "Performing " << num_levels << " levels of post-cdfem mesh refinement..." << std::endl;

      // Doing adaptive refinement with a uniform marker is better than doing uniform refinement here because of how
      // the transition elements are handled.
      auto & refinement = refinementSupport.get_non_interface_conforming_refinement();
      const stk::mesh::Selector & refinement_selector = cdfemSupport.get_post_cdfem_refinement_selector();
      std::function<void(int)> marker_function =
          [&refinement, num_levels, &refinement_selector](int num_refinements)
          {
            mark_selected_elements_for_refinement(refinement, num_refinements, num_levels, refinement_selector);
          };

      CDMesh::get_new_mesh()->delete_cdfem_parent_elements(); // Extreme work-around for the way percept messes up the active part on cdfem parents.
      perform_multilevel_adaptivity(refinement, mesh, marker_function);
    }
  }
}

}

//--------------------------------------------------------------------------------
void Region::initialize()
{ /* %TRACE[ON]% */ Trace trace__("krino::Region::initialize()"); /* %TRACE% */
  stk::diag::TimeBlock timer__(my_timerInitialize);

  const bool cdfem_is_active = krino::CDFEM_Support::is_active(get_stk_mesh_meta_data());
  krino::CDFEM_Support & cdfem_support = krino::CDFEM_Support::get(get_stk_mesh_meta_data());
  const RefinementSupport & refinementSupport = RefinementSupport::get(get_stk_mesh_meta_data());
  const Surface_Manager & surfaceManager = Surface_Manager::get(get_stk_mesh_meta_data());

  auto & bulk = get_stk_mesh_bulk_data();
  const auto & auxMeta = AuxMetaData::get(get_stk_mesh_meta_data());

  const bool doRefinement = refinementSupport.get_initial_refinement_levels() > 0 ||
      refinementSupport.get_interface_maximum_refinement_level() > 0 ||
      refinementSupport.get_post_adapt_refinement_levels() > 0;

  if (doRefinement)
  {
    cdfem_support.set_global_ids_are_NOT_parallel_consistent();

    auto & refinement = refinementSupport.get_non_interface_conforming_refinement();
    if (refinementSupport.get_initial_refinement_levels() > 0)
      refinement.do_uniform_refinement(refinementSupport.get_initial_refinement_levels());

    do_adaptive_refinement(refinementSupport, bulk);

    do_post_adapt_uniform_refinement(my_simulation, refinementSupport, auxMeta, bulk);
  }

  if (cdfem_is_active)
  {
    const int num_init_decomp_cycles = cdfem_support.get_num_initial_decomposition_cycles();
    for (int icycle=0; icycle<num_init_decomp_cycles; ++icycle)
    {
      LevelSet::initialize(get_stk_mesh_meta_data());
      if (icycle >= num_init_decomp_cycles-1)
        LevelSet::clear_initialization_data(get_stk_mesh_meta_data()); // reclaim memory, harmless to call multiple times after we are done initializing
      const std::unique_ptr<InterfaceGeometry> interfaceGeometry = create_interface_geometry(get_stk_mesh_meta_data());
      krino::CDMesh::decompose_mesh(bulk, *interfaceGeometry, my_simulation.get_time_step_count(), {});
    }

    do_post_cdfem_uniform_refinement(my_simulation, cdfem_support, refinementSupport, auxMeta, bulk);

    for (auto && part : krino::Phase_Support::get(get_stk_mesh_meta_data()).get_nonconformal_parts())
    {
      if (stk::io::is_part_io_part(*part)) stk::io::remove_io_part_attribute(*part);
    }
  }
  else
  {
    for(auto&& ls : surfaceManager.get_levelsets())
    {
      ls->initialize(0.);

      if (my_simulation.is_transient())
      {
        // initialize does not end with the facets constructed so manually construct them now
        ls->build_facets_locally(get_stk_mesh_meta_data().universal_part());

        // debugging
        if (krinolog.shouldPrint(LOG_FACETS))
        {
          ls->write_facets();
        }
      }
    }
  }

  LevelSet::clear_initialization_data(get_stk_mesh_meta_data()); // reclaim memory, harmless to call multiple times after we are done initializing

  // Skip output of empty parts
  for (auto && part : get_stk_mesh_meta_data().get_parts())
  {
    if (stk::io::is_part_io_part(*part))
    {
      uint64_t numEntities = stk::mesh::count_selected_entities(*part & auxMeta.active_part(), get_stk_mesh_bulk_data().buckets(part->primary_entity_rank()));
      const uint64_t localNumEntities = numEntities;
      stk::all_reduce_sum(get_stk_mesh_bulk_data().parallel(), &localNumEntities, &numEntities, 1);
      if(numEntities == 0)
      {
        krinolog << "Skipping output of empty part " << part->name() << stk::diag::dendl;
        stk::io::remove_io_part_attribute(*part);
      }
    }
  }
}
//--------------------------------------------------------------------------------

void Region::execute()
{ /* %TRACE[ON]% */ Trace trace__("krino::Region::execute()"); /* %TRACE% */
  stk::diag::TimeBlock timer__(my_timerExecute);

  // This is a hack for now to exit immediately when CDFEM is active
  if (krino::CDFEM_Support::is_active(get_stk_mesh_meta_data()))
  {
    return;
  }

  double deltaTime = time_step();

  const Surface_Manager & surfaceManager = Surface_Manager::get(get_stk_mesh_meta_data());
  for(auto&& ls : surfaceManager.get_levelsets())
  {
    ls->advance_semilagrangian(deltaTime);
  }
}

unsigned Region::spatial_dimension() const { return my_meta->spatial_dimension(); }
const stk::mesh::BulkData& Region::get_stk_mesh_bulk_data() const { ThrowRequire(my_bulk); return *my_bulk; }
stk::mesh::BulkData& Region::get_stk_mesh_bulk_data() { ThrowRequire(my_bulk); return *my_bulk; }
const stk::mesh::MetaData& Region::get_stk_mesh_meta_data() const { ThrowRequire(my_meta); return *my_meta; }
stk::mesh::MetaData& Region::get_stk_mesh_meta_data() { ThrowRequire(my_meta); return *my_meta; }
double Region::time_step() const { return my_simulation.get_time_step(); }

stk::io::StkMeshIoBroker & Region::stk_IO()
{
  ThrowRequire(myIOBroker);
  return *myIOBroker;
}

Ioss::Region * Region::get_input_io_region()
{
  return stk_IO().get_input_io_region().get();
}

std::string Region::name_of_input_mesh() const
{
  return my_input_model_name;
}

void Region::create_output_mesh()
{
  bool output_mesh = my_results_options->get_num_step_increments() != 0;

  if(!output_mesh) return;

  my_results_options->get_scheduler().set_termination_time(my_simulation.get_stop_time());

  my_output_file_index = stk_IO().create_output_mesh(my_results_options->get_filename(), stk::io::WRITE_RESULTS, my_results_options->get_properties());

  Teuchos::RCP<stk::mesh::Selector> active_selector = Teuchos::rcp(new stk::mesh::Selector(AuxMetaData::get(get_stk_mesh_meta_data()).active_part()));
  stk_IO().set_subset_selector(my_output_file_index, active_selector);

  stk_IO().write_output_mesh(my_output_file_index);

  declare_output_variables(my_output_file_index);

  my_output_file_created = true;
}

void
Region::declare_output_variables(size_t result_output_index)
{
  for (auto && outField : my_results_options->get_nodal_fields())
  {
    const std::string & varName = outField.first;
    const std::string & newName = outField.second;
    stk::mesh::FieldBase *theField = get_stk_mesh_meta_data().get_field(stk::topology::NODE_RANK, varName);
    if ( nullptr == theField )
    {
      krinolog << "Sorry, no nodal field by the name " << varName << std::endl;
      krinolog << "Available fields: " << std::endl;
      for (auto && field : get_stk_mesh_meta_data().get_fields(stk::topology::NODE_RANK))
      {
        krinolog << "  " << field->name() << std::endl;
      }
    }
    else
    {
      stk_IO().add_field(result_output_index, *theField, newName);
    }
  }

  for (auto && outField : my_results_options->get_element_fields())
  {
    const std::string & varName = outField.first;
    const std::string & newName = outField.second;
    stk::mesh::FieldBase *theField = get_stk_mesh_meta_data().get_field(stk::topology::ELEMENT_RANK, varName);
    if ( nullptr == theField )
    {
      krinolog << "Sorry, no element field by the name " << varName << std::endl;
      krinolog << "Available fields: " << std::endl;
      for (auto && field : get_stk_mesh_meta_data().get_fields(stk::topology::ELEMENT_RANK))
      {
        krinolog << "  " << field->name() << std::endl;
      }
    }
    else
    {
      stk_IO().add_field(result_output_index, *theField, newName);
    }
  }
}

void Region::process_output(bool forceOutput)
{
  stk::diag::TimeBlock mesh_output_timeblock(my_timerMeshOutput);

  if(my_results_options->get_num_step_increments() == 0) return;

  stk::util::Scheduler & scheduler = my_results_options->get_scheduler();

  if(forceOutput) scheduler.set_force_schedule();

  const int stepCounter = my_simulation.get_time_step_count();
  const double currentTime = my_simulation.get_current_time();

  const bool doOutput = scheduler.is_it_time(currentTime, stepCounter);

  if(!doOutput) return;

  // A couple of krino tests that just compute surface distances never call Region::initialize()
  // where we normally create the output mesh.
  if(!my_output_file_created) create_output_mesh();

  stk_IO().begin_output_step(my_output_file_index, currentTime);
  stk_IO().write_defined_output_fields(my_output_file_index);
  stk_IO().end_output_step(my_output_file_index);
}

void
Region::associate_input_mesh(const std::string & model_name, bool assert_32bit_ids, bool force_64bit_ids)
{
  stk::diag::TimeBlock mesh_input_timeblock(my_timerMeshInput);

  MeshInputOptions * db_options = MeshInputOptions::get(model_name);
  ThrowRequireMsg(db_options != nullptr, "No finite element model found with name " << model_name);
  my_input_model_name = model_name;

  ThrowRequire(db_options->is_valid());
  if (db_options->use_generated_mesh())
  {
    stk::topology generated_mesh_element_type = db_options->get_generated_mesh_element_type();
    std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
    entity_rank_names.push_back("FAMILY_TREE");
    my_generated_mesh = std::make_unique<BoundingBoxMesh>(generated_mesh_element_type,entity_rank_names);
    my_meta = &my_generated_mesh->meta_data();
    stk::mesh::Field<double> & coordsField = my_meta->declare_field<double>(stk::topology::NODE_RANK, "coordinates", 1);
    stk::mesh::put_field_on_mesh(coordsField, my_meta->universal_part(), generated_mesh_element_type.dimension(), nullptr);
  }
  else
  {
    const std::string directory = sierra::Env::working_directory().c_str();
    const std::string filename(Ioss::Utils::local_filename(db_options->get_filename(), db_options->get_filetype(), directory));

    stk::io::StkMeshIoBroker & stk_io = stk_IO();

    stk_io.property_add(Ioss::Property("MAXIMUM_NAME_LENGTH", 180));

    if (!db_options->get_decomposition_method().empty())
    {
      stk_io.property_add(Ioss::Property("DECOMPOSITION_METHOD", Ioss::Utils::uppercase(db_options->get_decomposition_method())));
    }

    stk_io.add_mesh_database(db_options->get_filename(), stk::io::READ_MESH);
    stk_io.create_input_mesh();
    my_meta = &stk_io.meta_data();
    AuxMetaData::create(*my_meta);
  }

  {
    stk::mesh::Selector activeSelector = AuxMetaData::get(*my_meta).active_part();
    stk_IO().set_active_selector(activeSelector);
  }

  if (assert_32bit_ids)
  {
    AuxMetaData::get(*my_meta).set_assert_32bit_flag();
  }
  if (!force_64bit_ids)
  {
    AuxMetaData::get(*my_meta).clear_force_64bit_flag();
  }
}

void
Region::set_generated_mesh_domain()
{
  MeshInputOptions * db_options = MeshInputOptions::get(my_input_model_name);
  ThrowRequire(db_options != nullptr);
  my_generated_mesh->set_mesh_structure_type(db_options->get_generated_mesh_structure_type());
  const double mesh_size = db_options->get_generated_mesh_size();
  if (db_options->get_generated_mesh_domain_type() == MeshInputOptions::GENERATED_MESH_FOR_SPECIFIED_DOMAIN)
  {
    typename BoundingBoxMesh::BoundingBoxType bbox;
    const std::vector<double> & domain = db_options->get_generated_mesh_domain();
    if (db_options->get_generated_mesh_spatial_dimension() == 2)
    {
      ThrowRequire(domain.size() == 4);
      const Vector3d min(domain[0], domain[1], 0.);
      const Vector3d max(domain[2], domain[3], 0.);
      ThrowRequireMsg(max[0]>min[0] && max[1]>min[1], "Invalid domain specified.");
      bbox = BoundingBoxMesh::BoundingBoxType(min, max);
    }
    else
    {
      ThrowRequire(domain.size() == 6);
      const Vector3d min(domain[0], domain[1], domain[2]);
      const Vector3d max(domain[3], domain[4], domain[5]);
      ThrowRequireMsg(max[0]>min[0] && max[1]>min[1] && max[2]>min[2], "Invalid domain specified.");
      bbox = BoundingBoxMesh::BoundingBoxType(min, max);
    }
    my_generated_mesh->set_domain(bbox, mesh_size);
  }
  else
  {
    typename BoundingBoxMesh::BoundingBoxType domain_bbox;
    const Surface_Manager & surfaceManager = Surface_Manager::get(get_stk_mesh_meta_data());
    for(auto&& ls : surfaceManager.get_levelsets())
    {
      if (ls->has_IC_surfaces())
      {
        const BoundingBox ls_IC_surf_bbox = ls->get_IC_surface_bounding_box();
        domain_bbox.accommodate(ls_IC_surf_bbox);
      }
    }
    for(auto&& boundingSurface : surfaceManager.get_bounding_surfaces())
    {
      const BoundingBox surfaceBBox = boundingSurface->surface().get_bounding_box();
      domain_bbox.accommodate(surfaceBBox);
    }

    if (db_options->get_generated_mesh_spatial_dimension() == 2)
    {
      Vector3d min = domain_bbox.get_min();
      Vector3d max = domain_bbox.get_max();
      min[2] = 0.;
      max[2] = 0.;
      domain_bbox = typename BoundingBoxMesh::BoundingBoxType(min, max);
    }
    my_generated_mesh->set_domain(domain_bbox, mesh_size, 1);
  }
}

} // namespace krino

