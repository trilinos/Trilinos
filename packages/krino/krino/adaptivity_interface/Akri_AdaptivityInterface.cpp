// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AdaptivityHelpers.hpp>
#include <Akri_AdaptivityInterface.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_RefinementInterface.hpp>

#ifdef __INTEL_COMPILER
#include <adapt/ElementRefinePredicate.hpp> // for ElementRefinePredicate
#include <adapt/RefinerUtil.hpp>
#include <adapt/TransitionElementAdapter.hpp>
#include <adapt/UniformRefiner.hpp> // for UniformRefiner
#include <adapt/UniformRefinerPattern.hpp>
#include <percept/FieldTypes.hpp>
#include <percept/PerceptMesh.hpp>
#include <adapt/markers/Marker.hpp>
#include <adapt/markers/MarkerUsingErrIndFraction.hpp>
#else
// this disables the checking of shadowed variables on GCC only
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <adapt/ElementRefinePredicate.hpp> // for ElementRefinePredicate
#include <adapt/RefinerUtil.hpp>
#include <adapt/TransitionElementAdapter.hpp>
#include <adapt/UniformRefiner.hpp> // for UniformRefiner
#include <adapt/UniformRefinerPattern.hpp>
#include <percept/FieldTypes.hpp>
#include <percept/PerceptMesh.hpp>
#include <adapt/markers/Marker.hpp>
#include <adapt/markers/MarkerUsingErrIndFraction.hpp>
#pragma GCC diagnostic pop
#endif

#include <stddef.h>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/Writer.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/diag/Timer.hpp>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

namespace krino
{

PerceptRefinement & create_percept_refinement(stk::mesh::MetaData & meta, stk::diag::Timer & parentTimer)
{
  auto & active_part = AuxMetaData::get(meta).active_part();
  HAdapt::setup(meta, active_part, parentTimer);
  PerceptRefinement & percepRefinement = PerceptRefinement::create(meta);
  auto adaptiveRefinement = [&meta](const std::string & markerFieldName, int debug_level)
    {
      HAdapt::do_adaptive_refinement(meta, markerFieldName);
    };
  percepRefinement.set_adaptive_refinement_function(adaptiveRefinement);
  auto uniformRefinement = [&meta](const int uniformRefinementLevels)
  {
    HAdapt::do_uniform_refinement(meta, uniformRefinementLevels);
  };
  percepRefinement.set_uniform_refinement_function(uniformRefinement);
  return percepRefinement;
}

RefinementInterface & create_refinement(stk::mesh::MetaData & meta, const bool usePercept, stk::diag::Timer & parentTimer)
{
  if (usePercept)
  {
    return create_percept_refinement(meta, parentTimer);
  }

  return KrinoRefinement::create(meta);
}

namespace
{

class HAdaptImpl
{
public:
  HAdaptImpl(stk::mesh::MetaData & meta)
      : my_meta(meta),
        my_active_part(nullptr),
        my_root_timer(nullptr)
  {
  }
  void setup(stk::mesh::Part & active_part, stk::diag::Timer & root_timer);
  void do_adaptive_refinement(const std::string & marker_field_name);
  void do_initial_uniform_refinement(const int num_levels);

private:
  void get_counts(const stk::mesh::FieldBase & marker_field,
      unsigned & refinement_count,
      unsigned & unrefinement_count) const;
  void check_supported_element_types() const;

  stk::mesh::MetaData & my_meta;
  std::unique_ptr<percept::PerceptMesh> my_pMesh;
  std::unique_ptr<percept::UniformRefinerPatternBase> my_adaptive_refinement_pattern;
  std::unique_ptr<percept::ElementRefinePredicate> my_element_refine_predicate;
  std::unique_ptr<percept::TransitionElementAdapter<percept::ElementRefinePredicate>> my_breaker;
  stk::mesh::Selector my_selector;
  stk::mesh::Part * my_active_part;
  stk::diag::Timer * my_root_timer;
};

HAdaptImpl & get(stk::mesh::MetaData & meta)
{
  HAdaptImpl * hadapt = const_cast<HAdaptImpl *>(meta.get_attribute<HAdaptImpl>());
  if (nullptr == hadapt)
  {
    hadapt = new HAdaptImpl(meta);
    meta.declare_attribute_with_delete<HAdaptImpl>(hadapt);
  }
  return *hadapt;
}

void HAdaptImpl::check_supported_element_types() const
{
  const stk::mesh::PartVector & all_parts = my_meta.get_parts();
  for (size_t i = 0; i < all_parts.size(); ++i)
  {
    const stk::mesh::Part & part = *all_parts[i];
    if (stk::io::is_part_io_part(part) && part.primary_entity_rank() == stk::topology::ELEMENT_RANK)
    {
      if ((2 == my_meta.spatial_dimension() &&
              part.topology() == stk::topology::TRIANGLE_3_2D) ||
          (3 == my_meta.spatial_dimension() &&
              part.topology() == stk::topology::TETRAHEDRON_4) ||
          part.topology() == stk::topology::PARTICLE)
      {
        continue;
      }
      ThrowErrorMsg("Elements in block "
          << part.name() << " have topology " << part.topology().name()
          << " which is currently not supported for adaptive refinement.");
    }
  }
}

void HAdaptImpl::get_counts(const stk::mesh::FieldBase & marker_field,
    unsigned & refinement_count,
    unsigned & unrefinement_count) const
{
  stk::mesh::BulkData & mesh = marker_field.get_mesh();
  stk::mesh::Part & active_part = *my_active_part;
  stk::mesh::Selector field_selector = stk::mesh::selectField(marker_field) & active_part & mesh.mesh_meta_data().locally_owned_part();

  unsigned local_refinement_count = 0;
  unsigned local_unrefinement_count = 0;

  const stk::mesh::BucketVector & buckets =
      mesh.get_buckets(stk::topology::ELEMENT_RANK, field_selector);

  for (stk::mesh::BucketVector::const_iterator ib = buckets.begin(); ib != buckets.end(); ++ib)
  {
    const stk::mesh::Bucket & b = **ib;
    const int length = b.size();

    for (int i = 0; i < length; ++i)
    {
      stk::mesh::Entity elem = b[i];

      const int & marker = *((int *)stk::mesh::field_data(marker_field, elem));
      if (marker < 0)
      {
        ++local_unrefinement_count;
      }
      else if (marker > 0)
      {
        ++local_refinement_count;
      }
    }
  }
  refinement_count = 0;
  unrefinement_count = 0;
  stk::all_reduce_sum(
      mesh.parallel(), (int *)&local_refinement_count, (int *)&refinement_count, 1);
  stk::all_reduce_sum(mesh.parallel(),
      (int *)&local_unrefinement_count,
      (int *)&unrefinement_count,
      1);
}

void HAdaptImpl::setup(stk::mesh::Part & active_part, stk::diag::Timer & root_timer)
{
  /* %TRACE[ON]% */ Trace trace__(
      "void HAdapt::setup(stk::mesh::Part & active_part)"); /* %TRACE% */

  my_active_part = &active_part;
  my_root_timer = &root_timer;

  static stk::diag::Timer timerAdapt_("Adapt", *my_root_timer);
  static stk::diag::Timer timerSetup_("Setup", timerAdapt_);
  stk::diag::TimeBlock tbTimerSetup_(timerSetup_);

  const unsigned nDim = my_meta.spatial_dimension();

  my_pMesh = std::make_unique<percept::PerceptMesh>(&my_meta, nullptr, false);

  my_pMesh->setCoordinatesField(const_cast<stk::mesh::FieldBase*>(my_meta.coordinate_field()));

  my_pMesh->register_and_set_refine_fields();
  my_pMesh->add_field_int("refine_field", stk::topology::NODE_RANK, 0);

  if (2 == nDim)
  {
    my_adaptive_refinement_pattern = std::make_unique<percept::Local_Tri3_Tri3_N_HangingNode>(*my_pMesh);
  }
  else
  {
    my_adaptive_refinement_pattern = std::make_unique<percept::Local_Tet4_Tet4_N_HangingNode>(*my_pMesh);
  }

  my_pMesh->output_active_children_only(true);
  my_pMesh->setProperty("Refiner_skip_side_part_fixes", "true");
}

namespace
{
void limit_coarsening_to_child_elements_with_all_siblings_marked_for_coarsening_and_no_children(
    percept::PerceptMesh & eMesh,
    const stk::mesh::FieldBase & marker_field)
{
  const stk::mesh::BulkData & mesh = *eMesh.get_bulk_data();
  const auto & local_buckets =
      mesh.get_buckets(stk::topology::ELEMENT_RANK, mesh.mesh_meta_data().locally_owned_part());
  stk::mesh::Part &child_part = *mesh.mesh_meta_data().get_part("refine_active_elements_part_3");

  std::vector<stk::mesh::Entity> children;
  std::vector<stk::mesh::Entity> childrenOfChild;

  // Only mark element for COARSEN if all children are marked COARSEN
  for (auto && b_ptr : local_buckets)
  {
    const bool is_child = b_ptr->member(child_part);
    for (auto && elem : *b_ptr)
    {
      eMesh.getChildren(elem, children, false, false);
      if (children.empty())
      {
        if (!is_child)
        {
          int & marker = (*((int *)stk::mesh::field_data( marker_field, elem )));
          if (marker == -1)
          {
            marker = 0;
          }
        }
      }
      else
      {
        bool all_children_marked_COARSEN = true;
        for (auto && child : children)
        {
          int & marker = (*((int *)stk::mesh::field_data( marker_field, child )));

          if (marker == -1)
            eMesh.getChildren(child, childrenOfChild, false, false);
          if (marker != -1 || !childrenOfChild.empty())
          {
            all_children_marked_COARSEN = false;
            break;
          }
        }
        if (!all_children_marked_COARSEN)
        {
          for (auto && child : children)
          {
            int & marker = (*((int *)stk::mesh::field_data( marker_field, child )));
            if (marker == -1)
            {
              marker = 0;
            }
          }
        }
      }
    }
  }
}


void fill_percept_refine_field(percept::PerceptMesh & eMesh,
    const stk::mesh::FieldBase & marker_field)
{
  limit_coarsening_to_child_elements_with_all_siblings_marked_for_coarsening_and_no_children(eMesh, marker_field);

  stk::mesh::BulkData & mesh = *eMesh.get_bulk_data();
  stk::mesh::communicate_field_data(mesh, {&marker_field});

  stk::mesh::Part &parent_part = *mesh.mesh_meta_data().get_part("refine_inactive_elements_part_3");
  stk::mesh::FieldBase & refine_field = *eMesh.m_refine_field;

  for (auto && b_ptr : mesh.get_buckets(stk::topology::ELEMENT_RANK, stk::mesh::selectField(marker_field)))
  {
    const stk::mesh::Bucket & b = *b_ptr;
    const bool is_parent = b.member(parent_part);
    const int length = b.size();
    const unsigned marker_field_length = stk::mesh::field_scalars_per_entity(marker_field, b);
    ThrowRequire(1 == stk::mesh::field_scalars_per_entity(refine_field, b));

    int * marker = (int *)stk::mesh::field_data(marker_field, b);
    int * refine = (int *)stk::mesh::field_data(refine_field, b);
    if (is_parent)
    {
      for (int i = 0; i < length; ++i)
      {
        marker[i*marker_field_length] = 0;
      }
    }

    for (int i = 0; i < length; ++i)
    {
      refine[i] = marker[i*marker_field_length];
    }
  }
}

void delete_partless_faces_and_edges(stk::mesh::BulkData & mesh)
{
  std::vector<stk::mesh::Entity> entities;
  std::vector<stk::mesh::Entity> relatives;
  std::vector<stk::mesh::ConnectivityOrdinal> relative_ordinals;

  mesh.modification_begin();

  for (stk::mesh::EntityRank entity_rank : {stk::topology::FACE_RANK, stk::topology::EDGE_RANK})
  {
    stk::mesh::get_entities(mesh, entity_rank, entities);

    for (auto && entity : entities)
    {
      stk::mesh::Bucket & bucket = mesh.bucket(entity);

      const stk::mesh::PartVector & bucket_parts = bucket.supersets();
      bool bucket_has_same_rank_part =
          std::find_if(bucket_parts.begin(),
              bucket_parts.end(),
              [&](const stk::mesh::Part * bucket_part)
              {
                return (bucket_part->primary_entity_rank() == bucket.entity_rank() &&
                    !stk::mesh::is_auto_declared_part(*bucket_part));
              }) != bucket_parts.end();

      if (!bucket_has_same_rank_part &&
          mesh.num_connectivity(entity, stk::topology::CONSTRAINT_RANK) == 0)
      {
        for (stk::mesh::EntityRank irank = stk::topology::ELEMENT_RANK; irank != entity_rank;
             --irank)
        {
          // Previously this attempted to delete forward or backward and still the list got
          // corrupted,
          // so just copy into vector and delete from there.
          relatives.assign(mesh.begin(entity, irank), mesh.end(entity, irank));
          relative_ordinals.assign(
              mesh.begin_ordinals(entity, irank), mesh.end_ordinals(entity, irank));

          for (size_t irel = 0; irel < relatives.size(); ++irel)
          {
            ThrowRequireMsg(mesh.destroy_relation(relatives[irel], entity, relative_ordinals[irel]),
                "Could not destroy relation between " << mesh.entity_key(relatives[irel]) << " and "
                                                      << mesh.entity_key(entity));
          }
        }
        ThrowRequireMsg(
            mesh.destroy_entity(entity), "Could not destroy entity " << mesh.entity_key(entity));
      }
    }
  }

  mesh.modification_end();
}

stk::mesh::Part & get_refinement_active_part(const stk::mesh::MetaData & meta, stk::mesh::EntityRank rank)
{
  const std::string active_part_name = "refine_active_elements_part_"+std::to_string((int)rank);
  stk::mesh::Part* active_part = meta.get_part(active_part_name);
  ThrowRequireMsg(nullptr != active_part, "Active part not found: " << active_part_name);
  return *active_part;
}

void update_active_inactive_entities(stk::mesh::BulkData & mesh,
    stk::mesh::Part & active_part)
{
  const auto & meta = mesh.mesh_meta_data();
  stk::mesh::PartVector active_part_vec(1, &active_part);
  stk::mesh::PartVector inactive_part_vec;
  stk::mesh::Selector select_locally_owned(meta.locally_owned_part());
  stk::mesh::Selector percept_active_selector = get_refinement_active_part(meta, stk::topology::ELEMENT_RANK);

  mesh.modification_begin();
  std::vector<stk::mesh::Entity> entities;
  for (stk::mesh::EntityRank entity_rank = stk::topology::NODE_RANK;
       entity_rank <= stk::topology::ELEMENT_RANK;
       ++entity_rank)
  {
    stk::mesh::get_selected_entities(select_locally_owned, mesh.buckets(entity_rank), entities);
    for (auto && entity : entities)
    {
      if (percept_active_selector(mesh.bucket(entity)))
        mesh.change_entity_parts(entity, active_part_vec, inactive_part_vec);
      else
        mesh.change_entity_parts(entity, inactive_part_vec, active_part_vec);
    }
  }
  mesh.modification_end();
}

void
fixup_side_permutation(stk::mesh::BulkData & mesh)
{
  mesh.modification_begin();
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  const stk::mesh::BucketVector & buckets = mesh.buckets(stk::topology::ELEMENT_RANK);

  for (auto&& bucket_ptr : buckets)
  {
    for (auto&& elem : *bucket_ptr)
    {
      const stk::mesh::Entity* elem_sides = mesh.begin(elem, meta.side_rank());
      const unsigned num_elem_sides = mesh.num_connectivity(elem, meta.side_rank());

      for (size_t it_side = 0; it_side < num_elem_sides; ++it_side)
      {
        auto relationship = determine_ordinal_and_permutation(mesh, elem, elem_sides[it_side]);
        set_relation_permutation(mesh, elem, elem_sides[it_side], relationship.first, relationship.second);
      }
    }
  }
  mesh.modification_end();
}

}

void HAdaptImpl::do_adaptive_refinement(const std::string & marker_field_name)
{
  /* %TRACE[ON]% */ Trace trace__(
      "void HAdapt::do_adaptive_refinement(const std::string &marker_field)"); /* %TRACE% */

  ThrowAssertMsg(my_root_timer != nullptr, "HAdapt::setup() not called.");
  static stk::diag::Timer timerAdapt_("Adapt", *my_root_timer);
  static stk::diag::Timer timerSetup_("Setup", timerAdapt_);
  static stk::diag::Timer timerFmwkUpdating_("Update Active Part", timerAdapt_);
  static stk::diag::Timer timerRefine_("Refine", timerAdapt_);
  static stk::diag::Timer timerUnrefine_("Unrefine", timerAdapt_);

  stk::diag::TimeBlock tbTimerAdapt_(timerAdapt_);

  constexpr bool debug = false;

  check_supported_element_types();

  stk::mesh::BulkData & bulk = my_meta.mesh_bulk_data();
  if (!my_pMesh->get_bulk_data())
  {
    my_pMesh->set_bulk_data(&bulk);
  }

  const stk::mesh::FieldBase * element_marker_field =
      my_meta.get_field(stk::topology::ELEMENT_RANK, marker_field_name);
  ThrowRequire(element_marker_field);
  my_selector = stk::mesh::selectField(*element_marker_field);

  if (!my_breaker)
  {
    const double tolerance = 0.0; // not used
    stk::mesh::FieldBase * refine_field = my_pMesh->m_refine_field;
    my_element_refine_predicate = std::make_unique<percept::ElementRefinePredicate>(
        *my_pMesh, &my_selector, refine_field, tolerance); // Note that this stores pointer to selector so temporary is not ok
    my_breaker = std::make_unique<percept::TransitionElementAdapter<percept::ElementRefinePredicate>>(
        *my_element_refine_predicate, *my_pMesh, *my_adaptive_refinement_pattern, nullptr, debug);
    my_breaker->setRemoveOldElements(false);
    my_breaker->setAlwaysInitializeNodeRegistry(false);
  }

  my_breaker->initializeRefine();

  unsigned refinement_count = 0;
  unsigned unrefinement_count = 0;

  fill_percept_refine_field(*my_pMesh, *element_marker_field);

  get_counts(*element_marker_field, refinement_count, unrefinement_count);
  krinolog << "Number of elements marked for refinement = " << refinement_count << "\n";
  krinolog << "Number of elements marked for unrefinement = " << unrefinement_count << stk::diag::dendl;

  const unsigned total_marked_count = refinement_count + unrefinement_count;
  if (0 == total_marked_count)
  {
    krinolog << "Adapt: No elements marked for refinement or unrefinement.  "
                "Skipping adaptivity" << stk::diag::dendl;
    return;
  }

  delete_partless_faces_and_edges(bulk);

  if (refinement_count > 0)
  {
    stk::diag::TimeBlock tbTimerRefine_(timerRefine_);
    my_breaker->setAlternateRootTimer(&timerRefine_);
    my_breaker->setModBegEndRootTimer(&timerAdapt_);

    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts);

    krinolog << "Adapt: before refine, mesh has  " << counts[0] << " nodes, " << counts[1]
             << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

    my_breaker->refine();

    stk::mesh::comm_mesh_counts(bulk, counts);

    krinolog << "Adapt: after refine, mesh has  " << counts[0] << " nodes, " << counts[1]
             << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

    my_breaker->setAlternateRootTimer(0);
    my_breaker->setModBegEndRootTimer(0);
  }

  // update refine field to include newly created children
  fill_percept_refine_field(*my_pMesh, *element_marker_field);

  if (unrefinement_count > 0)
  {
    stk::diag::TimeBlock tbTimerUnrefine_(timerUnrefine_);
    my_breaker->setAlternateRootTimer(&timerUnrefine_);
    my_breaker->setModBegEndRootTimer(&timerAdapt_);
    my_breaker->unrefine();
    my_breaker->setAlternateRootTimer(0);
    my_breaker->setModBegEndRootTimer(0);
  }

  {
    stk::diag::TimeBlock tbTimerRebuilding_(timerFmwkUpdating_);
    // The Encore-Percept interface also called induce_nodal_unranked_superset_parts
    // and topology nodeset inducer but I don't believe those are needed here.
    ThrowRequireMsg(my_active_part != nullptr, "Active part not set for krino::HAdapt");
    update_active_inactive_entities(bulk, *my_active_part);
    fixup_side_permutation(bulk);
  }
}

void HAdaptImpl::do_initial_uniform_refinement(const int num_levels)
{
  /* %TRACE[ON]% */ Trace trace__(
      "void HAdapt::do_uniform_refinement()"); /* %TRACE% */

  ThrowAssertMsg(my_root_timer != nullptr, "HAdapt::setup() not called.");
  static stk::diag::Timer timerAdapt_("Adapt", *my_root_timer);
  static stk::diag::Timer timerSetup_("Setup", timerAdapt_);
  static stk::diag::Timer timerFmwkUpdating_("Update Active Part", timerAdapt_);
  static stk::diag::Timer timerRefine_("Refine", timerAdapt_);

  stk::diag::TimeBlock tbTimerAdapt_(timerAdapt_);

  check_supported_element_types();

  stk::mesh::BulkData & bulk = my_meta.mesh_bulk_data();
  if (!my_pMesh->get_bulk_data())
  {
    my_pMesh->set_bulk_data(&bulk);
  }

  percept::UniformRefiner urefine(*my_pMesh, *my_adaptive_refinement_pattern, nullptr);
  urefine.setAlternateRootTimer(&timerAdapt_);

  {
    stk::diag::TimeBlock tbTimerRefine_(timerRefine_);

    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts);

    krinolog << "Adapt: before refine, mesh has  " << counts[0] << " nodes, " << counts[1]
             << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

    for (int level=0; level<num_levels; ++level)
    {
      urefine.doBreak();
    }

    stk::mesh::comm_mesh_counts(bulk, counts);

    krinolog << "Adapt: after refine, mesh has  " << counts[0] << " nodes, " << counts[1]
             << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;
  }

  {
    stk::diag::TimeBlock tbTimerRebuilding_(timerFmwkUpdating_);
    // The Encore-Percept interface also called induce_nodal_unranked_superset_parts
    // and topology nodeset inducer but I don't believe those are needed here.
    ThrowRequireMsg(my_active_part != nullptr, "Active part not set for krino::HAdapt");
    update_active_inactive_entities(bulk, *my_active_part);
    fixup_side_permutation(bulk);
  }
}

} // anonymous namespace

namespace HAdapt
{

void setup(stk::mesh::MetaData & meta, stk::mesh::Part & active_part, stk::diag::Timer & root_timer)
{
  get(meta).setup(active_part, root_timer);
}

void do_adaptive_refinement(stk::mesh::MetaData & meta, const std::string & marker_field_name)
{
  get(meta).do_adaptive_refinement(marker_field_name);
}

void do_uniform_refinement(stk::mesh::MetaData & meta, const int num_levels)
{
  get(meta).do_initial_uniform_refinement(num_levels);
}

}
} // namespace krino
