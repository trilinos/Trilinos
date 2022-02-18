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
#include <Akri_LevelSet.hpp>
#include <Akri_LevelSetInterfaceGeometry.hpp>
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
  my_output_file_created(false),
  my_initial_refinement_levels(0)
{ /* %TRACE[ON]% */ Trace trace__("krino::Region::Region()"); /* %TRACE% */
  my_simulation.add_region(this);
  myIOBroker = std::make_unique<stk::io::StkMeshIoBroker>(stk::EnvData::parallel_comm());

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
  krino::CDFEM_Support & cdfem_support = krino::CDFEM_Support::get(meta);

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

  if (my_initial_refinement_levels > 0 || (krino::CDFEM_Support::is_active(meta) && cdfem_support.get_interface_maximum_refinement_level() > 0) ||
      (krino::CDFEM_Support::is_active(meta) && cdfem_support.get_post_cdfem_refinement_levels() > 0))
  {
    auto_aura_option = stk::mesh::BulkData::AUTO_AURA;
    HAdapt::setup(meta, active_part, my_timerExecute);
  }

  if (cdfem_support.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
  {
    auto_aura_option = stk::mesh::BulkData::AUTO_AURA;
  }

  if (krino::CDFEM_Support::is_active(meta))
  {
    cdfem_support.add_interpolation_field(cdfem_support.get_coords_field());
    const LevelSetManager & region_ls = LevelSetManager::get(meta);
    for(auto&& ls : region_ls)
    {
      cdfem_support.add_interpolation_field(ls->get_distance_field());
    }
    cdfem_support.finalize_fields();

    if(cdfem_support.get_interface_maximum_refinement_level() > 0 ||
       cdfem_support.get_post_cdfem_refinement_levels() > 0)
    {
      cdfem_support.set_nonconformal_hadapt(
          [&meta](const std::string & marker_field, int debug_level)
          {
            HAdapt::do_adaptive_refinement(meta, marker_field);
          });
    }
  }

  if (nullptr != my_generated_mesh.get())
  {
    set_generated_mesh_domain();
    my_generated_mesh->populate_mesh(stk::EnvData::parallel_comm(), auto_aura_option);
    my_bulk = &my_generated_mesh->bulk_data();
    stk_IO().set_bulk_data( Teuchos::rcpFromRef( *my_bulk ) );
    if (my_generated_mesh->has_flat_boundaries())
      my_generated_mesh->create_domain_sides();
  }
  else
  {
    my_bulk = new stk::mesh::BulkData(meta,stk::EnvData::parallel_comm(),auto_aura_option);
    stk_IO().set_bulk_data( Teuchos::rcp( my_bulk ) );
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

void zero_error_indicator(const LevelSetManager & region_ls, stk::mesh::BulkData & mesh)
{
  auto & meta = mesh.mesh_meta_data();
  const auto & cdfem_support = CDFEM_Support::get(meta);
  const auto & aux_meta = AuxMetaData::get(meta);
  auto indicator_field =
      aux_meta.get_field(stk::topology::ELEMENT_RANK,
          cdfem_support.get_nonconformal_adapt_indicator_name());

  const auto & buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK,
      meta.locally_owned_part());
  for(auto && b_ptr : buckets)
  {
    auto * data = field_data<double>(indicator_field, *b_ptr);
    auto length = b_ptr->size();
    std::fill(data, data + length, 0.);
  }
}

void mark_selected_elements(const std::string & marker_field_name, const int current_refinement_level, const AuxMetaData  & auxMeta, const int max_refinement_levels, const stk::mesh::Selector & selector)
{
  FieldRef marker_field = auxMeta.get_field(stk::topology::ELEMENT_RANK, marker_field_name);
  const bool do_adapt = current_refinement_level < max_refinement_levels;
  if (do_adapt)
  {
    stk::mesh::field_fill(1, marker_field, selector);
  }
}

void do_adaptive_refinement(const krino::CDFEM_Support & cdfemSupport, const LevelSetManager & region_ls, stk::mesh::BulkData & mesh)
{
  const auto target_count = cdfemSupport.get_nonconformal_adapt_target_count();
  if(target_count > 0) zero_error_indicator(region_ls, mesh);
  {
    stk::diag::TimeBlock adapt_timer__(cdfemSupport.get_timer_adapt());

    const std::string & marker_name = "refine_field"; //cdfemSupport.get_nonconformal_adapt_marker_name();
    auto & h_adapt = cdfemSupport.get_nonconformal_hadapt();

    const int num_refinement_levels = cdfemSupport.get_interface_maximum_refinement_level();
    const int num_refinement_steps = (target_count > 0) ? num_refinement_levels : 2*num_refinement_levels;

    const auto & indicator_field_name = cdfemSupport.get_nonconformal_adapt_indicator_name();

    const bool do_post_adapt_refinement = cdfemSupport.get_post_adapt_refinement_levels() > 0;
    std::function<void(const std::string &, int)> marker_function =
        [&mesh, &region_ls, num_refinement_steps, do_post_adapt_refinement, target_count, indicator_field_name]
         (const std::string & marker_field_name, int num_refinements)
        {
          const bool do_adapt = num_refinements < num_refinement_steps;
          const bool more_to_do = do_adapt || do_post_adapt_refinement;
          if(target_count > 0) zero_error_indicator(region_ls, mesh);
          LevelSet::initialize(mesh.mesh_meta_data(), more_to_do);
          if(!do_adapt) return;

          if(target_count == 0)
          {
            const LevelSetInterfaceGeometry interfaceGeometry(mesh.mesh_meta_data());
            CDMesh::mark_interface_elements_for_adaptivity(mesh, interfaceGeometry, marker_field_name, num_refinements);
          }
          else
          {
            HAdapt::mark_based_on_indicator_field(mesh,
                marker_field_name,
                indicator_field_name,
                num_refinement_steps,
                num_refinements,
                target_count);
          }
        };

    perform_multilevel_adaptivity(mesh, marker_name, marker_function, h_adapt, cdfem_do_not_refine_or_unrefine_selector(cdfemSupport));
  }
}

void do_post_adapt_uniform_refinement(const Simulation & simulation, const krino::CDFEM_Support & cdfemSupport, const AuxMetaData  & auxMeta, stk::mesh::BulkData & mesh)
{
  if ( cdfemSupport.get_post_adapt_refinement_levels() > 0 )
  {
    if ( simulation.is_transient() )
    {
      krinolog << "Can not do post-adaptivity uniform mesh refinement for transient problems...";
    }
    else
    {
      krinolog << "Performing " << cdfemSupport.get_post_adapt_refinement_levels() << " levels of post-adaptivity uniform mesh refinement..." << std::endl;

      // Doing adaptive refinement with a uniform marker is better than doing uniform refinement here because of how
      // the transition elements are handled.
      const std::string & marker_name = "refine_field";
      auto & h_adapt = cdfemSupport.get_nonconformal_hadapt();
      const int num_levels = cdfemSupport.get_post_adapt_refinement_levels();
      FieldRef marker_field = auxMeta.get_field(stk::topology::ELEMENT_RANK, marker_name);

      std::function<void(const std::string &, int)> marker_function =
          [&mesh, marker_field, num_levels]
           (const std::string & marker_field_name, int num_refinements)
          {
            const bool do_adapt = num_refinements < num_levels;
            LevelSet::initialize(mesh.mesh_meta_data(), do_adapt);
            if (do_adapt)
            {
              stk::mesh::field_fill(1, marker_field);
            }
          };

      perform_multilevel_adaptivity(mesh, marker_name, marker_function, h_adapt, cdfem_do_not_refine_or_unrefine_selector(cdfemSupport));
    }
  }
}

void do_post_cdfem_uniform_refinement(const Simulation & simulation, const krino::CDFEM_Support & cdfemSupport, const AuxMetaData  & auxMeta, stk::mesh::BulkData & mesh)
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
      const std::string & marker_name = "refine_field";
      auto & h_adapt = cdfemSupport.get_nonconformal_hadapt();
      const stk::mesh::Selector & refinement_selector = cdfemSupport.get_post_cdfem_refinement_selector();
      std::function<void(const std::string &, int)> marker_function =
          [&auxMeta, num_levels, &refinement_selector]
           (const std::string & marker_field_name, int num_refinements)
          {
            mark_selected_elements(marker_field_name, num_refinements, auxMeta, num_levels, refinement_selector);
          };

      CDMesh::get_new_mesh()->delete_cdfem_parent_elements(); // Extreme work-around for the way percept messes up the active part on cdfem parents.
      perform_multilevel_adaptivity(mesh, marker_name, marker_function, h_adapt);
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
  const bool adapt_is_active = cdfem_support.get_interface_maximum_refinement_level() > 0 ||
      cdfem_support.get_post_cdfem_refinement_levels() > 0;

  const LevelSetManager & region_ls = LevelSetManager::get(get_stk_mesh_meta_data());

  if (my_initial_refinement_levels > 0)
  {
    cdfem_support.set_global_ids_are_NOT_parallel_consistent();
    HAdapt::do_uniform_refinement(get_stk_mesh_meta_data(), my_initial_refinement_levels);
  }

  if (cdfem_is_active)
  {
    auto & bulk = get_stk_mesh_bulk_data();
    if(adapt_is_active)
    {
      do_adaptive_refinement(cdfem_support, region_ls, bulk);

      const auto & auxMeta = AuxMetaData::get(get_stk_mesh_meta_data());
      do_post_adapt_uniform_refinement(my_simulation, cdfem_support, auxMeta, bulk);

      const LevelSetInterfaceGeometry interfaceGeometry(get_stk_mesh_meta_data());
      krino::CDMesh::decompose_mesh(bulk, interfaceGeometry, my_simulation.get_time_step_count(), {});

      do_post_cdfem_uniform_refinement(my_simulation, cdfem_support, auxMeta, bulk);
    }
    else
    {
      if ( cdfem_support.get_post_adapt_refinement_levels() > 0 )
      {
        krinolog << "Post-adaptivity uniform refinement levels selected without any adaptivity levels.  Not performing additional uniform refinement." << std::endl;
      }

      const int num_init_decomp_cycles = cdfem_support.get_num_initial_decomposition_cycles();
      for (int icycle=0; icycle<num_init_decomp_cycles; ++icycle )
      {
        const bool not_last_init = (icycle < (num_init_decomp_cycles-1));
        LevelSet::initialize(get_stk_mesh_meta_data(), not_last_init);
        const LevelSetInterfaceGeometry interfaceGeometry(get_stk_mesh_meta_data());
        krino::CDMesh::decompose_mesh(bulk, interfaceGeometry, my_simulation.get_time_step_count(), {});
      }
    }

    for (auto && part : krino::Phase_Support::get(get_stk_mesh_meta_data()).get_nonconformal_parts())
    {
      if (stk::io::is_part_io_part(*part)) stk::io::remove_io_part_attribute(*part);
    }
  }
  else
  {
    for(auto&& ls : region_ls)
    {
      ls->initialize(0.);

      if (my_simulation.is_transient())
      {
        // initialize does not end with the facets constructed so manually construct them now
        ls->build_facets_locally(get_stk_mesh_meta_data().universal_part());

        // debugging
        if (krinolog.shouldPrint(LOG_FACETS))
        {
          ls->facets_exoii();
        }
      }
    }
  }

  // Skip output of empty parts
  for (auto && part : get_stk_mesh_meta_data().get_parts())
  {
    if (stk::io::is_part_io_part(*part))
    {
      uint64_t local_num_entities = stk::mesh::count_selected_entities(
          *part, get_stk_mesh_bulk_data().buckets(part->primary_entity_rank()));
      uint64_t global_num_entities = 0;
      stk::all_reduce_sum(
          get_stk_mesh_bulk_data().parallel(), &local_num_entities, &global_num_entities, 1);
      if(global_num_entities == 0)
      {
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

  const LevelSetManager & region_ls = LevelSetManager::get(get_stk_mesh_meta_data());
  for(auto&& ls : region_ls)
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
    stk::mesh::Field<double, stk::mesh::Cartesian> & coords_field = my_meta->declare_field<stk::mesh::Field<double, stk::mesh::Cartesian>>(stk::topology::NODE_RANK, "coordinates", 1);
    stk::mesh::put_field_on_mesh(coords_field, my_meta->universal_part(), generated_mesh_element_type.dimension(), nullptr);
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
    auto & active_part = AuxMetaData::get(*my_meta).active_part();
    stk::mesh::Selector activeSelector = active_part;
    stk_io.set_active_selector(activeSelector);
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
    const LevelSetManager & region_ls = LevelSetManager::get(get_stk_mesh_meta_data());
    for(auto&& ls : region_ls)
    {
      if (ls->has_IC_surfaces())
      {
        BoundingBox ls_IC_surf_bbox = ls->get_IC_surface_bounding_box();
        domain_bbox.accommodate(ls_IC_surf_bbox);
      }
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

