/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_io/MeshReadWriteUtils.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_io/IossBridge.hpp>

#include <Ioss_SubSystem.h>

#include <stk_util/util/tokenize.hpp>
#include <iostream>
#include <sstream>
#include <cmath>

#include <limits>
#include <assert.h>

namespace {
  struct RestartFieldAttribute {
    RestartFieldAttribute(const std::string &name) : databaseName(name)
    {}
    std::string databaseName;
  };

  void process_surface_entity(Ioss::SideSet *sset, stk::mesh::MetaData &meta)
  {
    assert(sset->type() == Ioss::SIDESET);
    const Ioss::SideBlockContainer& blocks = sset->get_side_blocks();
    stk::io::default_part_processing(blocks, meta);

    stk::mesh::Part* const ss_part = meta.get_part(sset->name());
    assert(ss_part != NULL);

    stk::mesh::Field<double, stk::mesh::ElementNode> *distribution_factors_field = NULL;
    bool surface_df_defined = false; // Has the surface df field been defined yet?

    size_t block_count = sset->block_count();
    for (size_t i=0; i < block_count; i++) {
      Ioss::SideBlock *sb = sset->get_block(i);
      if (stk::io::include_entity(sb)) {
        stk::mesh::Part * const sb_part = meta.get_part(sb->name());
        assert(sb_part != NULL);
        meta.declare_part_subset(*ss_part, *sb_part);

        if (sb->field_exists("distribution_factors")) {
          if (!surface_df_defined) {
            std::string field_name = sset->name() + "_df";
            distribution_factors_field =
                &meta.declare_field<stk::mesh::Field<double, stk::mesh::ElementNode> >(field_name);
            stk::io::set_field_role(*distribution_factors_field, Ioss::Field::MESH);
            stk::io::set_distribution_factor_field(*ss_part, *distribution_factors_field);
            surface_df_defined = true;
          }
          stk::io::set_distribution_factor_field(*sb_part, *distribution_factors_field);
          int side_node_count = sb->topology()->number_nodes();
          stk::mesh::put_field(*distribution_factors_field,
                               stk::io::part_primary_entity_rank(*sb_part),
                               *sb_part, side_node_count);
        }
      }
    }
  }

  size_t get_entities(stk::mesh::Part &part,
                      const stk::mesh::BulkData &bulk,
                      std::vector<stk::mesh::Entity> &entities,
                      const stk::mesh::Selector *anded_selector)
  {
    stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    stk::mesh::EntityRank type = stk::io::part_primary_entity_rank(part);

    stk::mesh::Selector own = meta.locally_owned_part();
    stk::mesh::Selector selector = part & own;
    if (anded_selector) selector &= *anded_selector;

    get_selected_entities(selector, bulk.buckets(type), entities);
    return entities.size();
  }

  bool is_restart_field_on_part(const stk::mesh::FieldBase *field,
				const stk::mesh::EntityRank part_type,
				const stk::mesh::Part &part)
  {
    return (field->attribute<RestartFieldAttribute>() != NULL &&
	    stk::io::is_field_on_part(field,part_type,part));
  }

  std::string get_restart_field_name(const stk::mesh::FieldBase *field)
  {
    ThrowRequire(field != NULL);
    const RestartFieldAttribute *attr = field->attribute<RestartFieldAttribute>();
    ThrowRequire(attr != NULL);
    
    return attr->databaseName;
  }

  void ioss_define_restart_fields(const stk::mesh::Part &part,
				  const stk::mesh::EntityRank part_type,
				  Ioss::GroupingEntity *entity)
  {
    const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

    std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
    while (I != fields.end()) {
      const stk::mesh::FieldBase *f = *I; ++I;
      if (is_restart_field_on_part(f, part_type, part)) {
	// Only add TRANSIENT Fields -- check role; if not present assume transient...
	const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);

	if (role == NULL || *role == Ioss::Field::TRANSIENT) {
	  const stk::mesh::FieldBase::Restriction &res = f->restriction(part_type, part);
	  std::pair<std::string, Ioss::Field::BasicType> field_type;
	  stk::io::get_io_field_type(f, res, &field_type);
	  if (field_type.second != Ioss::Field::INVALID) {
	    size_t entity_size = entity->get_property("entity_count").get_int();
	    const std::string& name = get_restart_field_name(f);
	    entity->field_add(Ioss::Field(name, field_type.second, field_type.first,
					  Ioss::Field::TRANSIENT, entity_size));
	  }
	}
      }
    }
  }

  int ioss_restore_restart_fields(stk::mesh::BulkData &bulk,
				  const stk::mesh::Part &part,
				  const stk::mesh::EntityRank part_type,
				  Ioss::GroupingEntity *io_entity)
  {
    int missing_fields = 0;

    const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

    assert(io_entity != NULL);
    std::vector<stk::mesh::Entity> entity_list;
    bool entity_list_filled=false;

    std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
    while (I != fields.end()) {
      const stk::mesh::FieldBase *f = *I; ++I;
      if (is_restart_field_on_part(f, part_type, part)) {
	// Only add TRANSIENT Fields -- check role; if not present assume transient...
	const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);

	if (role == NULL || *role == Ioss::Field::TRANSIENT) {
	  const stk::mesh::FieldBase::Restriction &res = f->restriction(part_type, part);
	  std::pair<std::string, Ioss::Field::BasicType> field_type;
	  stk::io::get_io_field_type(f, res, &field_type);
	  if (field_type.second != Ioss::Field::INVALID) {
	    const std::string& name = get_restart_field_name(f);

	    // See if field with that name exists on io_entity...
	    if (io_entity->field_exists(name)) {
	      // Restore data...
	      if (!entity_list_filled) {
		stk::io::get_entity_list(io_entity, part_type, bulk, entity_list);
		entity_list_filled=true;
	      }
	      stk::io::field_data_from_ioss(bulk, f, entity_list, io_entity, name);
	    } else {
	      std::cerr  << "ERROR: Could not find restart input field '"
			 << name << "' on '" << io_entity->name() << "'.\n";

	      missing_fields++;
	    }
	  }
	}
      }
    }
    return missing_fields;
  }

}

// ========================================================================
template <typename INT>
void process_surface_entity(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk, INT /*dummy*/)
{
  assert(sset->type() == Ioss::SIDESET);

  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);

  size_t block_count = sset->block_count();
  for (size_t i=0; i < block_count; i++) {
    Ioss::SideBlock *block = sset->get_block(i);
    if (stk::io::include_entity(block)) {
      std::vector<INT> side_ids ;
      std::vector<INT> elem_side ;

      stk::mesh::Part * const sb_part = meta.get_part(block->name());
      stk::mesh::EntityRank elem_rank = stk::mesh::MetaData::ELEMENT_RANK;

      block->get_field_data("ids", side_ids);
      block->get_field_data("element_side", elem_side);

      assert(side_ids.size() * 2 == elem_side.size());
      stk::mesh::PartVector add_parts( 1 , sb_part );

      // Get topology of the sides being defined to see if they
      // are 'faces' or 'edges'.  This is needed since for shell-type
      // elements, (and actually all elements) a sideset can specify either a face or an edge...
      // For a quad shell, sides 1,2 are faces and 3,4,5,6 are edges.
      int par_dimen = block->topology()->parametric_dimension();

      size_t side_count = side_ids.size();
      for(size_t is=0; is<side_count; ++is) {
        stk::mesh::Entity const elem = bulk.get_entity(elem_rank, elem_side[is*2]);

        // If NULL, then the element was probably assigned to an
        // element block that appears in the database, but was
        // subsetted out of the analysis mesh. Only process if
        // non-null.
        if (bulk.is_valid(elem)) {
          // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
          int side_ordinal = elem_side[is*2+1] - 1;

          if (par_dimen == 1) {
            stk::mesh::Entity side = stk::mesh::declare_element_edge(bulk, side_ids[is], elem, side_ordinal);
            bulk.change_entity_parts( side, add_parts );
          }
          else if (par_dimen == 2) {
            stk::mesh::Entity side = stk::mesh::declare_element_side(bulk, side_ids[is], elem, side_ordinal);
            bulk.change_entity_parts( side, add_parts );
          }
        }
      }
    }
  }
}

// ========================================================================
template <typename INT>
void process_surface_entity_df(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk, INT /*dummy*/)
{
  assert(sset->type() == Ioss::SIDESET);

  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);
  bool check_pre_existing = true;

  size_t block_count = sset->block_count();
  for (size_t i=0; i < block_count; i++) {
    Ioss::SideBlock *block = sset->get_block(i);
    if (stk::io::include_entity(block)) {
      std::vector<INT> side_ids ;
      std::vector<INT> elem_side ;

      stk::mesh::Part * const sb_part = meta.get_part(block->name());
      stk::mesh::EntityRank elem_rank = stk::mesh::MetaData::ELEMENT_RANK;

      block->get_field_data("ids", side_ids);
      block->get_field_data("element_side", elem_side);

      assert(side_ids.size() * 2 == elem_side.size());

      // Get topology of the sides being defined to see if they
      // are 'faces' or 'edges'.  This is needed since for shell-type
      // elements, (and actually all elements) a sideset can specify either a face or an edge...
      // For a quad shell, sides 1,2 are faces and 3,4,5,6 are edges.
      int par_dimen = block->topology()->parametric_dimension();

      size_t side_count = side_ids.size();
      std::vector<stk::mesh::Entity> sides(side_count);
      for(size_t is=0; is<side_count; ++is) {
        stk::mesh::Entity const elem = bulk.get_entity(elem_rank, elem_side[is*2]);

        // If NULL, then the element was probably assigned to an
        // element block that appears in the database, but was
        // subsetted out of the analysis mesh. Only process if
        // non-null.
        if (bulk.is_valid(elem)) {
          // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
          int side_ordinal = elem_side[is*2+1] - 1;

          if (par_dimen == 1) {
            stk::mesh::Entity side = stk::mesh::declare_element_edge(bulk, side_ids[is], elem, side_ordinal, NULL, check_pre_existing);
            sides[is] = side;
          }
          else if (par_dimen == 2) {
            stk::mesh::Entity side = stk::mesh::declare_element_side(bulk, side_ids[is], elem, side_ordinal, NULL, check_pre_existing);
            sides[is] = side;
          }
        } else {
          sides[is] = stk::mesh::Entity();
        }
      }

      const stk::mesh::FieldBase *df_field = stk::io::get_distribution_factor_field(*sb_part);
      if (df_field != NULL) {
        stk::io::field_data_from_ioss(bulk, df_field, sides, block, "distribution_factors");
      }

      // Add all attributes as fields.
      // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
      // named components of the 'attribute' field, so add them instead.
      Ioss::NameList names;
      block->field_describe(Ioss::Field::ATTRIBUTE, &names);
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(*I == "attribute" && names.size() > 1)
          continue;
        stk::mesh::FieldBase *field = meta.get_field<stk::mesh::FieldBase> (*I);
        if (field)
          stk::io::field_data_from_ioss(bulk, field, sides, block, *I);
      }
    }
  }
}

void process_surface_entity(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk)
{
  if (stk::io::db_api_int_size(sset) == 4)
    process_surface_entity(sset, bulk, (int)0);
  else
    process_surface_entity(sset, bulk, (int64_t)0);
}

void process_surface_entity_df(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk)
{
  if (stk::io::db_api_int_size(sset) == 4)
    process_surface_entity_df(sset, bulk, (int)0);
  else
    process_surface_entity_df(sset, bulk, (int64_t)0);
}

void process_nodeblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field =
      meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(stk::io::CoordinateFieldName);
  stk::io::set_field_role(coord_field, Ioss::Field::MESH);

  meta.set_coordinate_field(&coord_field);
  
  Ioss::NodeBlock *nb = node_blocks[0];
  stk::mesh::put_field(coord_field, stk::mesh::MetaData::NODE_RANK, meta.universal_part(),
                       meta.spatial_dimension());
  stk::io::define_io_fields(nb, Ioss::Field::ATTRIBUTE, meta.universal_part(), 0);
}

template <typename INT>
void process_nodeblocks(stk::io::MeshData &mesh, INT /*dummy*/)
{
  stk::mesh::BulkData &bulk = mesh.bulk_data();
  // This must be called after the "process_element_blocks" call
  // since there may be nodes that exist in the database that are
  // not part of the analysis mesh due to subsetting of the element
  // blocks.

  // Currently, all nodes found in the finite element mesh are defined
  // as nodes in the stk_mesh database. If some of the element blocks
  // are omitted, then there will be disconnected nodes defined.
  // However, if we only define nodes that are connected to elements,
  // then we risk missing "free" nodes that the user may want to have
  // existing in the model.
  const Ioss::NodeBlockContainer& node_blocks = mesh.input_io_region()->get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  std::vector<INT> ids;
  nb->get_field_data("ids", ids);

  for (size_t i=0; i < ids.size(); i++) {
    stk::mesh::Entity node = bulk.declare_entity(stk::mesh::MetaData::NODE_RANK, ids[i]);
    bulk.set_local_id(node, i);
  }
}

template <typename INT>
void process_node_coords_and_attributes(stk::io::MeshData &mesh, INT /*dummy*/)
{
  stk::mesh::BulkData &bulk = mesh.bulk_data();
  // This must be called after the "process_element_blocks" call
  // since there may be nodes that exist in the database that are
  // not part of the analysis mesh due to subsetting of the element
  // blocks.

  // Currently, all nodes found in the finite element mesh are defined
  // as nodes in the stk_mesh database. If some of the element blocks
  // are omitted, then there will be disconnected nodes defined.
  // However, if we only define nodes that are connected to elements,
  // then we risk missing "free" nodes that the user may want to have
  // existing in the model.
  const Ioss::NodeBlockContainer& node_blocks = mesh.input_io_region()->get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  size_t node_count = nb->get_property("entity_count").get_int();

  std::vector<INT> ids;
  nb->get_field_data("ids", ids);

  std::vector<stk::mesh::Entity> nodes;
  nodes.reserve(node_count);
  for (size_t i=0; i < ids.size(); i++) {
    stk::mesh::Entity node = bulk.get_entity(stk::mesh::MetaData::NODE_RANK, ids[i]);
    nodes.push_back(node);
  }

  // Temporary (2013/04/02) kluge for Salinas porting to stk-based mesh.
  // Salinas uses the "implicit id" which is the ordering of the nodes
  // in the "undecomposed" or "serial" mesh as user-visible ids
  // instead of the "global" ids. If there exists a stk-field with the
  // name "implicit_ids", then populate the field with the correct
  // data.
  stk::mesh::FieldBase *imp_id_field = mesh.meta_data().get_field<stk::mesh::FieldBase> ("implicit_ids");
  if (imp_id_field) {
    stk::io::field_data_from_ioss(bulk, imp_id_field, nodes, nb, "implicit_ids");
  }


  stk::mesh::FieldBase *coord_field = &mesh.get_coordinate_field();
  stk::io::field_data_from_ioss(bulk, coord_field, nodes, nb, "mesh_model_coordinates");

  // Add all attributes as fields.
  // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
  // named components of the 'attribute' field, so add them instead.
  Ioss::NameList names;
  nb->field_describe(Ioss::Field::ATTRIBUTE, &names);
  for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
    if(*I == "attribute" && names.size() > 1)
      continue;
    stk::mesh::FieldBase *field = mesh.meta_data().get_field<stk::mesh::FieldBase> (*I);
    if (field)
      stk::io::field_data_from_ioss(bulk, field, nodes, nb, *I);
  }
}

// ========================================================================
void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  stk::io::default_part_processing(elem_blocks, meta);
}

template <typename INT>
void process_elementblocks(Ioss::Region &region, stk::mesh::BulkData &bulk, INT /*dummy*/)
{
  const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
      it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string &name = entity->name();
      stk::mesh::Part* const part = meta.get_part(name);
      assert(part != NULL);

      stk::topology topo = part->topology();
      if (topo == stk::topology::INVALID_TOPOLOGY) {
        std::ostringstream msg ;
        msg << " INTERNAL_ERROR: Part " << part->name() << " has invalid topology";
        throw std::runtime_error( msg.str() );
      }

      std::vector<INT> elem_ids ;
      std::vector<INT> connectivity ;

      entity->get_field_data("ids", elem_ids);
      entity->get_field_data("connectivity", connectivity);

      size_t element_count = elem_ids.size();
      int nodes_per_elem = topo.num_nodes();

      std::vector<stk::mesh::EntityId> id_vec(nodes_per_elem);

      size_t offset = entity->get_offset();
      for(size_t i=0; i<element_count; ++i) {
        INT *conn = &connectivity[i*nodes_per_elem];
        std::copy(&conn[0], &conn[0+nodes_per_elem], id_vec.begin());
        stk::mesh::Entity element = stk::mesh::declare_element(bulk, *part, elem_ids[i], &id_vec[0]);

        bulk.set_local_id(element, offset + i);
      }
    }
  }
}

template <typename INT>
void process_elem_attributes_and_implicit_ids(Ioss::Region &region, stk::mesh::BulkData &bulk, INT /*dummy*/)
{
  const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
      it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string &name = entity->name();
      stk::mesh::Part* const part = meta.get_part(name);
      assert(part != NULL);

      stk::topology topo = part->topology();
      if (topo == stk::topology::INVALID_TOPOLOGY) {
        std::ostringstream msg ;
        msg << " INTERNAL_ERROR: Part " << part->name() << " has invalid topology";
        throw std::runtime_error( msg.str() );
      }

      // See if we need to get the list of elements...
      bool elements_needed = false;
      Ioss::NameList names;
      entity->field_describe(Ioss::Field::ATTRIBUTE, &names);

      stk::mesh::FieldBase *imp_id_field = meta.get_field<stk::mesh::FieldBase> ("implicit_ids");
      if (imp_id_field) {
	elements_needed = true;
      } else {
	for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
	  if(*I == "attribute" && names.size() > 1)
	    continue;
	  stk::mesh::FieldBase *field = meta.get_field<stk::mesh::FieldBase> (*I);
	  if (field) {
	    elements_needed = true;
	    break;
	  }
	}
      }
      
      if (!elements_needed)
	continue;
      
      std::vector<INT> elem_ids ;
      entity->get_field_data("ids", elem_ids);

      size_t element_count = elem_ids.size();
      std::vector<stk::mesh::Entity> elements;
      elements.reserve(element_count);

      for(size_t i=0; i<element_count; ++i) {
        stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEMENT_RANK, elem_ids[i]);
	elements.push_back(elem);
      }

      // Temporary (2013/04/17) kluge for Salinas porting to stk-based mesh.
      // Salinas uses the "implicit id" which is the ordering of the nodes
      // in the "undecomposed" or "serial" mesh as user-visible ids
      // instead of the "global" ids. If there exists a stk-field with the
      // name "implicit_ids", then populate the field with the correct
      // data.
      if (imp_id_field) {
        stk::io::field_data_from_ioss(bulk, imp_id_field, elements, entity, "implicit_ids");
      }

      // Add all element attributes as fields.
      // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
      // named components of the 'attribute' field, so add them instead.
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(*I == "attribute" && names.size() > 1)
          continue;
        stk::mesh::FieldBase *field = meta.get_field<stk::mesh::FieldBase> (*I);
        if (field)
          stk::io::field_data_from_ioss(bulk, field, elements, entity, *I);
      }
    }
  }
}

// ========================================================================
// ========================================================================
void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  stk::io::default_part_processing(node_sets, meta);

  stk::mesh::Field<double> & distribution_factors_field =
      meta.declare_field<stk::mesh::Field<double> >("distribution_factors");
  stk::io::set_field_role(distribution_factors_field, Ioss::Field::MESH);

  /** \todo REFACTOR How to associate distribution_factors field
   * with the nodeset part if a node is a member of multiple
   * nodesets
   */

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());

      assert(part != NULL);
      assert(entity->field_exists("distribution_factors"));

      stk::io::set_field_role(distribution_factors_field, Ioss::Field::MESH);
      stk::mesh::put_field(distribution_factors_field, stk::mesh::MetaData::NODE_RANK, *part);
    }
  }

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());

      assert(part != NULL);
      assert(entity->field_exists("distribution_factors"));

      std::string nodesetName = part->name();
      std::string nodesetDistFieldName = "distribution_factors_" + nodesetName;

      stk::mesh::Field<double> & distribution_factors_field_per_nodeset =
           meta.declare_field<stk::mesh::Field<double> >(nodesetDistFieldName);

      stk::io::set_field_role(distribution_factors_field_per_nodeset, Ioss::Field::MESH);
      stk::mesh::put_field(distribution_factors_field_per_nodeset, stk::mesh::MetaData::NODE_RANK, *part);
    }
  }
}

// ========================================================================
void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();
  stk::io::default_part_processing(side_sets, meta);

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity(entity, meta);
    }
  }
}

// ========================================================================
template <typename INT>
void process_nodesets(Ioss::Region &region, stk::mesh::BulkData &bulk, INT /*dummy*/)
{
  // Should only process nodes that have already been defined via the element
  // blocks connectivity lists.
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string & name = entity->name();
      stk::mesh::Part* const part = meta.get_part(name);
      assert(part != NULL);
      stk::mesh::PartVector add_parts( 1 , part );

      std::vector<INT> node_ids ;
      size_t node_count = entity->get_field_data("ids", node_ids);

      stk::mesh::EntityRank n_rank = stk::mesh::MetaData::NODE_RANK;
      for(size_t i=0; i<node_count; ++i) {
        stk::mesh::Entity node = bulk.get_entity(n_rank, node_ids[i] );
        if (!bulk.is_valid(node)) {
          node = bulk.declare_entity(n_rank, node_ids[i], add_parts );
        }
        else {
          bulk.change_entity_parts(node, add_parts);
        }
      }
    }
  }
}

// ========================================================================
template <typename INT>
void process_nodesets_df(Ioss::Region &region, stk::mesh::BulkData &bulk, INT /*dummy*/)
{
  // Should only process nodes that have already been defined via the element
  // blocks connectivity lists.
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string & name = entity->name();
      stk::mesh::Part* const part = meta.get_part(name);
      assert(part != NULL);
      stk::mesh::PartVector add_parts( 1 , part );

      std::vector<INT> node_ids ;
      size_t node_count = entity->get_field_data("ids", node_ids);

      std::vector<stk::mesh::Entity> nodes(node_count);
      stk::mesh::EntityRank n_rank = stk::mesh::MetaData::NODE_RANK;
      for(size_t i=0; i<node_count; ++i) {
        nodes[i] = bulk.get_entity(n_rank, node_ids[i] );
        if (!bulk.is_valid(nodes[i])) {
          bulk.declare_entity(n_rank, node_ids[i], add_parts );
        }
      }

      stk::mesh::Field<double> *df_field =
          meta.get_field<stk::mesh::Field<double> >("distribution_factors");

      if (df_field != NULL) {
        stk::io::field_data_from_ioss(bulk, df_field, nodes, entity, "distribution_factors");
      }

      std::string distributionFactorsPerNodesetFieldName = "distribution_factors_" + part->name();

      stk::mesh::Field<double> *df_field_per_nodeset =
                meta.get_field<stk::mesh::Field<double> >(distributionFactorsPerNodesetFieldName);

      if (df_field_per_nodeset != NULL) {
        stk::io::field_data_from_ioss(bulk, df_field_per_nodeset, nodes, entity, "distribution_factors");
      }

      // Add all attributes as fields.
      // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
      // named components of the 'attribute' field, so add them instead.
      Ioss::NameList names;
      entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(*I == "attribute" && names.size() > 1)
          continue;
        stk::mesh::FieldBase *field = meta.get_field<stk::mesh::FieldBase> (*I);
        if (field)
          stk::io::field_data_from_ioss(bulk, field, nodes, entity, *I);
      }
    }
  }
}

// ========================================================================
void process_sidesets(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity(entity, bulk);
    }
  }
}

// ========================================================================
void process_sidesets_df(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity_df(entity, bulk);
    }
  }
}

std::string get_field_name(const stk::mesh::FieldBase &f, Ioss::DatabaseUsage dbUsage)
{
  std::string name = f.name();
  if(dbUsage == Ioss::WRITE_RESTART || dbUsage == Ioss::READ_RESTART) {
    const RestartFieldAttribute *restartAttribute = f.attribute<RestartFieldAttribute>();
    if(restartAttribute !=NULL) {
      name = restartAttribute->databaseName;
    }
  }
  else {
    name = stk::io::get_results_field_name(f);
  }
  return name;
}

// ========================================================================
void put_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
                    stk::mesh::EntityRank part_type,
                    Ioss::GroupingEntity *io_entity,
                    Ioss::Field::RoleType filter_role,
                    const stk::mesh::Selector *anded_selector=NULL)
{
  std::vector<stk::mesh::Entity> entities;
  if (io_entity->type() == Ioss::SIDEBLOCK) {
    // Temporary Kluge to handle sideblocks which contain internally generated sides
    // where the "ids" field on the io_entity doesn't work to get the correct side...
    // NOTE: Could use this method for all entity types, but then need to correctly
    // specify whether shared entities are included/excluded (See IossBridge version).
    size_t num_sides = get_entities(part, bulk, entities, anded_selector);
    if (num_sides != (size_t)io_entity->get_property("entity_count").get_int()) {
      std::ostringstream msg ;
      msg << " INTERNAL_ERROR: Number of sides on part " << part.name() << " (" << num_sides
          << ") does not match number of sides in the associated Ioss SideBlock named "
          << io_entity->name() << " (" << io_entity->get_property("entity_count").get_int()
          << ").";
      throw std::runtime_error( msg.str() );
    }
  } else {
    stk::io::get_entity_list(io_entity, part_type, bulk, entities);
  }

  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);
  const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

  std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const stk::mesh::FieldBase *f = *I; ++I;
    std::string fieldName = get_field_name(*f, io_entity->get_database()->usage());
    stk::io::field_data_to_ioss(bulk, f, entities, io_entity, fieldName, filter_role);
  }
}

namespace stk {
  namespace io {

    MeshData::MeshData()
      : m_communicator_(MPI_COMM_NULL), m_anded_selector(NULL), m_connectivity_map(stk::mesh::ConnectivityMap::default_map()),
	useNodesetForPartNodesFields(true)
    {
      Ioss::Init::Initializer::initialize_ioss();
    }

    MeshData::MeshData(MPI_Comm comm, stk::mesh::ConnectivityMap * connectivity_map)
      : m_communicator_(comm), m_anded_selector(NULL)
      , m_connectivity_map(connectivity_map != NULL ? *connectivity_map : stk::mesh::ConnectivityMap::default_map()),
	useNodesetForPartNodesFields(true)
    {
      Ioss::Init::Initializer::initialize_ioss();
    }

    MeshData::~MeshData()
    {
      if (!Teuchos::is_null(m_output_region))
        stk::io::delete_selector_property(*m_output_region);

      if (!Teuchos::is_null(m_restart_region))
        stk::io::delete_selector_property(*m_restart_region);
    }

    stk::mesh::FieldBase & MeshData::get_coordinate_field()
    {
      stk::mesh::FieldBase * coord_field = meta_data().coordinate_field();
      ThrowRequire( coord_field != NULL);
      return * coord_field;
    }

    void MeshData::set_output_io_region(Teuchos::RCP<Ioss::Region> ioss_output_region)
    {
      m_output_region = ioss_output_region;
    }

    void MeshData::set_input_io_region(Teuchos::RCP<Ioss::Region> ioss_input_region)
    {
      ThrowErrorMsgIf(!Teuchos::is_null(m_input_region),
		      "This MeshData already has an Ioss::Region associated with it.");
      m_input_region = ioss_input_region;
    }

    void MeshData::set_meta_data( Teuchos::RCP<stk::mesh::MetaData> arg_meta_data )
    {
      ThrowErrorMsgIf( !Teuchos::is_null(m_meta_data),
		       "Meta data already initialized" );
      m_meta_data = arg_meta_data;
    }

    void MeshData::set_bulk_data( Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data )
    {
      ThrowErrorMsgIf( !Teuchos::is_null(m_bulk_data),
		       "Bulk data already initialized" );
      m_bulk_data = arg_bulk_data;

      if (Teuchos::is_null(m_meta_data)) {
        set_meta_data(const_cast<stk::mesh::MetaData&>(bulk_data().mesh_meta_data()));
      }

      m_communicator_ = m_bulk_data->parallel();
    }

    bool MeshData::open_mesh_database(const std::string &mesh_filename)
    {
      std::string type = "exodus";
      std::string filename = mesh_filename;

      // See if filename contains a ":" at the beginning of the filename
      // and if the text preceding that filename specifies a valid
      // database type.  If so, set that as the file type and extract
      // the portion following the colon as the filename.
      // If no colon in name, use default type.

      size_t colon = mesh_filename.find(':');
      if (colon != std::string::npos && colon > 0) {
        type = mesh_filename.substr(0, colon);
        filename = mesh_filename.substr(colon+1);
      }
      return open_mesh_database(filename, type);
    }


    bool MeshData::open_mesh_database(const std::string &mesh_filename,
                                      const std::string &mesh_type)
    {
      ThrowErrorMsgIf(!Teuchos::is_null(m_input_database),
		      "This MeshData already has an Ioss::DatabaseIO associated with it.");
      ThrowErrorMsgIf(!Teuchos::is_null(m_input_region),
                      "This MeshData already has an Ioss::Region associated with it.");

      m_input_database = Teuchos::rcp(Ioss::IOFactory::create(mesh_type, mesh_filename,
                                                              Ioss::READ_MODEL, m_communicator_,
                                                              m_property_manager));
      if (Teuchos::is_null(m_input_database) || !m_input_database->ok(true)) {
        std::cerr  << "ERROR: Could not open database '" << mesh_filename
		   << "' of type '" << mesh_type << "'\n";
        return false;
      }
      return true;
    }


    void MeshData::create_ioss_region()
    {
      // If the m_input_region is null, try to create it from
      // the m_input_database. If that is null, throw an error.
      if (Teuchos::is_null(m_input_region)) {
        ThrowErrorMsgIf(Teuchos::is_null(m_input_database),
			"There is no input mesh database associated with this MeshData. Please call open_mesh_database() first.");
        // The Ioss::Region takes control of the m_input_database pointer, so we need to make sure the
        // RCP doesn't retain ownership...
        m_input_region = Teuchos::rcp(new Ioss::Region(m_input_database.release().get(), "input_model"));
      }
    }

    void MeshData::set_rank_name_vector(const std::vector<std::string> &rank_names)
    {
      m_rank_names.clear();
      std::copy(rank_names.begin(), rank_names.end(), std::back_inserter(m_rank_names));
    }

    void MeshData::create_input_mesh()
    {
      if (Teuchos::is_null(m_input_region)) {
        create_ioss_region();
      }

      // See if meta data is null, if so, create a new one...
      if (Teuchos::is_null(m_meta_data)) {
        set_meta_data(Teuchos::rcp( new stk::mesh::MetaData()));
      }

      size_t spatial_dimension = m_input_region->get_property("spatial_dimension").get_int();
      if (m_rank_names.empty()) {
        initialize_spatial_dimension(meta_data(), spatial_dimension, stk::mesh::entity_rank_names());
      } else {
        initialize_spatial_dimension(meta_data(), spatial_dimension, m_rank_names);
      }

      process_nodeblocks(*m_input_region.get(),    meta_data());
      process_elementblocks(*m_input_region.get(), meta_data());
      process_sidesets(*m_input_region.get(),      meta_data());
      process_nodesets(*m_input_region.get(),      meta_data());
    }


    void MeshData::create_output_mesh(const std::string &filename)
    {
      std::string out_filename = filename;
      if (filename.empty()) {
        out_filename = "default_output_mesh";
      } else {
        // These filenames may be coming from the generated options which
        // may have forms similar to: "2x2x1|size:.05|height:-0.1,1"
        // Strip the name at the first "+:|," character:
        std::vector<std::string> tokens;
        stk::util::tokenize(out_filename, "+|:,", tokens);
        out_filename = tokens[0];
      }

      Ioss::DatabaseIO *dbo = Ioss::IOFactory::create("exodusII", out_filename,
                                                      Ioss::WRITE_RESULTS,
                                                      m_communicator_,
                                                      m_property_manager);
      if (dbo == NULL || !dbo->ok()) {
        std::cerr << "ERROR: Could not open results database '" << out_filename
		  << "' of type 'exodusII'\n";
        std::exit(EXIT_FAILURE);
      }

      // NOTE: 'out_region' owns 'dbo' pointer at this time...
      if (!Teuchos::is_null(m_output_region))
        m_output_region = Teuchos::null;
      m_output_region = Teuchos::rcp(new Ioss::Region(dbo, "results_output"));

      create_output_mesh();
    }

    void MeshData::create_output_mesh()
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_output_region),
                       "There is no Output database associated with this Mesh Data.");
      define_output_database();
      write_output_database();
    }

    // ========================================================================
    int MeshData::process_output_request(double time,
                                         const std::set<const stk::mesh::Part*> &exclude)
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_output_region),
		       "There is no Output mesh region associated with this Mesh Data.");
      m_output_region->begin_mode(Ioss::STATE_TRANSIENT);

      int out_step = m_output_region->add_state(time);
      internal_process_output_request(out_step,exclude);

      m_output_region->end_mode(Ioss::STATE_TRANSIENT);

      return out_step;
    }

    // ========================================================================
    int MeshData::process_output_request(double time)
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_output_region),
		       "There is no Output mesh region associated with this Mesh Data.");
      m_output_region->begin_mode(Ioss::STATE_TRANSIENT);

      int out_step = m_output_region->add_state(time);
      internal_process_output_request(out_step,std::set<const stk::mesh::Part*>());

      m_output_region->end_mode(Ioss::STATE_TRANSIENT);

      return out_step;
    }

    void MeshData::populate_mesh(bool delay_field_data_allocation)
    {
      if (!meta_data().is_commit())
        meta_data().commit();

      ThrowErrorMsgIf (Teuchos::is_null(m_input_region),
                       "There is no Input mesh region associated with this Mesh Data.");

      Ioss::Region *region = m_input_region.get();
      ThrowErrorMsgIf (region==NULL,
                       "INTERNAL ERROR: Mesh Input Region pointer is NULL in populate_mesh.");

      // Check if bulk_data is null; if so, create a new one...
      if (Teuchos::is_null(m_bulk_data)) {
        set_bulk_data(Teuchos::rcp( new stk::mesh::BulkData(   meta_data()
							       , region->get_database()->util().communicator()
#ifdef SIERRA_MIGRATION
							       , false
#endif
							       , &m_connectivity_map
							       )));
      }

      if (delay_field_data_allocation) {
        bulk_data().deactivate_field_updating();
      }

      bool i_started_modification_cycle = bulk_data().modification_begin();

      bool ints64bit = db_api_int_size(region) == 8;
      if (ints64bit) {
        int64_t zero = 0;
        process_nodeblocks(*this, zero);
        process_elementblocks(*region, bulk_data(), zero);
        process_nodesets(*region,      bulk_data(), zero);
        process_sidesets(*region,      bulk_data());
      } else {
        int zero = 0;
        process_nodeblocks(*this, zero);
        process_elementblocks(*region, bulk_data(), zero);
        process_nodesets(*region,      bulk_data(), zero);
        process_sidesets(*region,      bulk_data());
      }

      if (i_started_modification_cycle) {
        bulk_data().modification_end();
      }

      if (region->get_property("state_count").get_int() == 0) {
        region->get_database()->release_memory();
      }
    }

    void MeshData::populate_field_data()
    {
      //if field-data has already been allocated, then the allocate_field_data() method
      //is a harmless no-op.
      bulk_data().allocate_field_data();

      Ioss::Region *region = m_input_region.get();
      ThrowErrorMsgIf (region==NULL,
                       "INTERNAL ERROR: Mesh Input Region pointer is NULL in populate_field_data.");

      bool ints64bit = db_api_int_size(region) == 8;
      if (ints64bit) {
          int64_t zero = 0;
          process_node_coords_and_attributes(*this, zero);
          process_elem_attributes_and_implicit_ids(*region, bulk_data(), zero);
          process_nodesets_df(*region,      bulk_data(), zero);
          process_sidesets_df(*region,      bulk_data());
      }
      else {
          int zero = 0;
          process_node_coords_and_attributes(*this, zero);
          process_elem_attributes_and_implicit_ids(*region, bulk_data(), zero);
          process_nodesets_df(*region,      bulk_data(), zero);
          process_sidesets_df(*region,      bulk_data());
      }
    }

    // ========================================================================
    void MeshData::populate_bulk_data()
    {
      if (!meta_data().is_commit())
        meta_data().commit();

      ThrowErrorMsgIf (Teuchos::is_null(m_input_region),
                       "There is no Input mesh region associated with this Mesh Data.");

      Ioss::Region *region = m_input_region.get();
      ThrowErrorMsgIf (region==NULL,
                       "INTERNAL ERROR: Mesh Input Region pointer is NULL in populate_mesh.");

      // Check if bulk_data is null; if so, create a new one...
      if (Teuchos::is_null(m_bulk_data)) {
        set_bulk_data(Teuchos::rcp( new stk::mesh::BulkData(   meta_data()
                                                             , region->get_database()->util().communicator()
#ifdef SIERRA_MIGRATION
                                                             , false
#endif
                                                             , &m_connectivity_map
                                                           )));
      }

      //to preserve behavior for callers of this method, don't do the
      //delay-field-data-allocation optimization.
      //Folks who want the optimization can call the population_mesh/populate_field_data methods separately.
      bool delay_field_data_allocation = false;
      bulk_data().modification_begin();
      populate_mesh(delay_field_data_allocation);
      populate_field_data();
      bulk_data().modification_end();
    }

    void MeshData::internal_process_output_request(int step, const std::set<const stk::mesh::Part*> &exclude)
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_output_region),
		       "There is no Output mesh region associated with this Mesh Data.");

      Ioss::Region *region = m_output_region.get();
      ThrowErrorMsgIf (region==NULL,
                       "INTERNAL ERROR: Mesh Output Region pointer is NULL in internal_process_output_request.");
      region->begin_state(step);

      // Special processing for nodeblock (all nodes in model)...
      put_field_data(bulk_data(), meta_data().universal_part(), stk::mesh::MetaData::NODE_RANK,
                     region->get_node_blocks()[0], Ioss::Field::Field::TRANSIENT,
                     m_anded_selector.get());

      // Now handle all non-nodeblock parts...
      const stk::mesh::PartVector &all_parts = meta_data().get_parts();
      for ( stk::mesh::PartVector::const_iterator
	      ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

        stk::mesh::Part * const part = *ip;

        // Check whether this part should be output to results database.
        if (stk::io::is_part_io_part(*part) && !exclude.count(part)) {
          stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
          // Get Ioss::GroupingEntity corresponding to this part...
          Ioss::GroupingEntity *entity = region->get_entity(part->name());
          if (entity != NULL && entity->type() != Ioss::SIDESET) {
            put_field_data(bulk_data(), *part, rank, entity,
			   Ioss::Field::Field::TRANSIENT, m_anded_selector.get());
          }

          // If rank is != NODE_RANK, then see if any fields are defined on the nodes of this part
          // (should probably do edges and faces also...)
          // Get Ioss::GroupingEntity corresponding to the nodes on this part...
          if (rank != stk::mesh::MetaData::NODE_RANK && use_nodeset_for_part_nodes_fields()) {
            std::string nodes_name = part->name() + "_nodes";
            Ioss::GroupingEntity *node_entity = region->get_entity(nodes_name);
            if (node_entity != NULL) {
              put_field_data(bulk_data(), *part, stk::mesh::MetaData::NODE_RANK, node_entity,
			     Ioss::Field::Field::TRANSIENT, m_anded_selector.get());
            }
          }
        }
      }
      region->end_state(step);
    }

    void MeshData::define_output_database()
    {
      bool sort_stk_parts = false; // used in stk_adapt/stk_percept
      stk::io::define_output_db(*m_output_region.get(), bulk_data(), m_input_region.get(), m_anded_selector.get(),
                                sort_stk_parts, use_nodeset_for_part_nodes_fields());
    }

    void MeshData::write_output_database()
    {
      stk::io::write_output_db(*m_output_region.get(),  bulk_data(), m_anded_selector.get());
    }

    void MeshData::add_restart_field(stk::mesh::FieldBase &field, std::string db_name)
    {
      // NOTE: This could be implemented as a free function; however, I want to keep the option
      //       open of storing this data in a vector/map on the MeshData class instead of as a
      //       field attribute.  For now, use the field attribute approach.
      if (db_name.empty()) {
	db_name = field.name();
      }

      RestartFieldAttribute *my_field_attribute = new RestartFieldAttribute(db_name);
      stk::mesh::MetaData &m = mesh::MetaData::get(field);
      const RestartFieldAttribute *check = m.declare_attribute_with_delete(field, my_field_attribute);
      if ( check != my_field_attribute ) {
	if (check->databaseName != my_field_attribute->databaseName) {
	  std::ostringstream msg ;
	  msg << "ERROR in MeshReadWriteUtils::add_restart_field:"
	      << " Conflicting restart file names for restart field '"
	      << field.name() << "'. Name was previously specified as '"
	      << check->databaseName << "', but is now being specified as '"
	      << my_field_attribute->databaseName << "'.";
	  throw std::runtime_error( msg.str() );
	}
	delete my_field_attribute;
      }
    }

    void MeshData::create_restart_output(const std::string &filename)
    {
      Ioss::DatabaseIO *dbo = Ioss::IOFactory::create("exodusII", filename,
                                                      Ioss::WRITE_RESTART,
                                                      m_communicator_,
                                                      m_property_manager);
      if (dbo == NULL || !dbo->ok()) {
        std::cerr << "ERROR: Could not open restart database '" << filename
		  << "' of type 'exodusII'\n";
        std::exit(EXIT_FAILURE);
      }

      // NOTE: 'out_region' owns 'dbo' pointer at this time...
      if (!Teuchos::is_null(m_restart_region))
        m_restart_region = Teuchos::null;
      m_restart_region = Teuchos::rcp(new Ioss::Region(dbo, "restart_output"));
      create_restart_output();
    }

    void MeshData::create_restart_output()
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_restart_region),
		       "There is no Restart database associated with this Mesh Data.");

      bool sort_stk_parts = false; // used in stk_adapt/stk_percept
      stk::io::define_output_db(*m_restart_region.get(), bulk_data(), m_input_region.get(), m_anded_selector.get(),
                                sort_stk_parts, use_nodeset_for_part_nodes_fields());
      stk::io::write_output_db(*m_restart_region.get(),  bulk_data(), m_anded_selector.get());
    }

    void MeshData::define_restart_fields()
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_restart_region),
		       "There is no Restart database associated with this Mesh Data.");

      Ioss::Region *region = m_restart_region.get();
      region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

      // Special processing for nodeblock (all nodes in model)...
      ioss_define_restart_fields(meta_data().universal_part(),
				 stk::mesh::MetaData::NODE_RANK,
				 region->get_node_blocks()[0]);

      const stk::mesh::PartVector &all_parts = meta_data().get_parts();
      for ( stk::mesh::PartVector::const_iterator
	      ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

	stk::mesh::Part * const part = *ip;

	// Check whether this part should be output to restart database.
	if (stk::io::is_part_io_part(*part)) {
	  stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
	  // Get Ioss::GroupingEntity corresponding to this part...
	  Ioss::GroupingEntity *entity = region->get_entity(part->name());
	  if (entity != NULL) {
	    ioss_define_restart_fields(*part, rank, entity);
	  }

	  // If rank is != NODE_RANK, then see if any fields are defined on the nodes of this part
	  // (should probably do edges and faces also...)
	  // Get Ioss::GroupingEntity corresponding to the nodes on this part...
	  if (rank != stk::mesh::MetaData::NODE_RANK) {
	    Ioss::GroupingEntity *node_entity = NULL;
	    if (use_nodeset_for_part_nodes_fields()) {
	      std::string nodes_name = part->name() + "_nodes";
	      node_entity = region->get_entity(nodes_name);
	    } else {
	      node_entity = region->get_entity("nodeblock_1");
	    }
	    if (node_entity != NULL) {
	      ioss_define_restart_fields(*part, stk::mesh::MetaData::NODE_RANK, node_entity);
	    }
	  }
	}
      }
      region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
    }

    int MeshData::process_restart_output(double time)
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_restart_region),
		       "There is no Restart output associated with this Mesh Data.");

      m_restart_region->begin_mode(Ioss::STATE_TRANSIENT);

      int out_step = m_restart_region->add_state(time);
      internal_process_restart_output(out_step);

      m_restart_region->end_mode(Ioss::STATE_TRANSIENT);

      return out_step;
    }

    void MeshData::internal_process_restart_output(int step)
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_restart_region),
		       "There is no Restart output region associated with this Mesh Data.");

      Ioss::Region *region = m_restart_region.get();
      ThrowErrorMsgIf (region==NULL,
                       "INTERNAL ERROR: Restart output region pointer is NULL in internal_process_restart_output.");
      region->begin_state(step);

      // Special processing for nodeblock (all nodes in model)...
      put_field_data(bulk_data(), meta_data().universal_part(), stk::mesh::MetaData::NODE_RANK,
                     region->get_node_blocks()[0], Ioss::Field::Field::TRANSIENT,
                     m_anded_selector.get());

      // Now handle all non-nodeblock parts...
      const stk::mesh::PartVector &all_parts = meta_data().get_parts();
      for ( stk::mesh::PartVector::const_iterator
	      ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

        stk::mesh::Part * const part = *ip;

        // Check whether this part should be output to results database.
        if (stk::io::is_part_io_part(*part)) {
          stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
          // Get Ioss::GroupingEntity corresponding to this part...
          Ioss::GroupingEntity *entity = region->get_entity(part->name());
          if (entity != NULL && entity->type() != Ioss::SIDESET) {
            put_field_data(bulk_data(), *part, rank, entity,
			   Ioss::Field::Field::TRANSIENT, m_anded_selector.get());
          }

          // If rank is != NODE_RANK, then see if any fields are defined on the nodes of this part
          // (should probably do edges and faces also...)
          // Get Ioss::GroupingEntity corresponding to the nodes on this part...
          if (rank != stk::mesh::MetaData::NODE_RANK && use_nodeset_for_part_nodes_fields()) {
            std::string nodes_name = part->name() + "_nodes";
            Ioss::GroupingEntity *node_entity = region->get_entity(nodes_name);
            if (node_entity != NULL) {
              put_field_data(bulk_data(), *part, stk::mesh::MetaData::NODE_RANK, node_entity,
			     Ioss::Field::Field::TRANSIENT, m_anded_selector.get());
            }
          }
        }
      }
      region->end_state(step);
    }

    namespace {
      // ========================================================================
      // Transfer transient field data from mesh file for io_entity to
      // the corresponding stk_mesh entities If there is a stk_mesh
      // field with the same name as the database field.
      // Assumes that mesh is positioned at the correct state for reading.
      void internal_process_input_request(Ioss::GroupingEntity *io_entity,
                                          stk::mesh::EntityRank entity_rank,
                                          stk::mesh::BulkData &bulk)
      {

        assert(io_entity != NULL);
        std::vector<stk::mesh::Entity> entity_list;
        stk::io::get_entity_list(io_entity, entity_rank, bulk, entity_list);

        const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);

        Ioss::NameList names;
        io_entity->field_describe(Ioss::Field::TRANSIENT, &names);
        for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
          stk::mesh::FieldBase *field = meta.get_field<stk::mesh::FieldBase>(*I);
          if (field) {
            stk::io::field_data_from_ioss(bulk, field, entity_list, io_entity, *I);
          }
        }
      }

      void input_nodeblock_fields(Ioss::Region &region, stk::mesh::BulkData &bulk)
      {
        const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
        assert(node_blocks.size() == 1);

        Ioss::NodeBlock *nb = node_blocks[0];
        internal_process_input_request(nb, stk::mesh::MetaData::NODE_RANK, bulk);
      }

      void input_elementblock_fields(Ioss::Region &region, stk::mesh::BulkData &bulk)
      {
        const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
        for(size_t i=0; i < elem_blocks.size(); i++) {
          if (stk::io::include_entity(elem_blocks[i])) {
            internal_process_input_request(elem_blocks[i], stk::mesh::MetaData::ELEMENT_RANK, bulk);
          }
        }
      }

      void input_nodeset_fields(Ioss::Region &region, stk::mesh::BulkData &bulk)
      {
        const Ioss::NodeSetContainer& nodesets = region.get_nodesets();
        for(size_t i=0; i < nodesets.size(); i++) {
          if (stk::io::include_entity(nodesets[i])) {
            internal_process_input_request(nodesets[i], stk::mesh::MetaData::NODE_RANK, bulk);
          }
        }
      }

      void input_sideset_fields(Ioss::Region &region, stk::mesh::BulkData &bulk)
      {
        const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);
        if (meta.spatial_dimension() <= meta.side_rank())
          return;

        const Ioss::SideSetContainer& side_sets = region.get_sidesets();
        for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
            it != side_sets.end(); ++it) {
          Ioss::SideSet *entity = *it;
          if (stk::io::include_entity(entity)) {
            const Ioss::SideBlockContainer& blocks = entity->get_side_blocks();
            for(size_t i=0; i < blocks.size(); i++) {
              if (stk::io::include_entity(blocks[i])) {
                internal_process_input_request(blocks[i], meta.side_rank(), bulk);
              }
            }
          }
        }
      }

      void define_input_nodeblock_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
        assert(node_blocks.size() == 1);

        Ioss::NodeBlock *nb = node_blocks[0];
        stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT,
                                  meta.universal_part(), stk::mesh::MetaData::NODE_RANK);
      }

      void define_input_elementblock_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
        for(size_t i=0; i < elem_blocks.size(); i++) {
          if (stk::io::include_entity(elem_blocks[i])) {
            stk::mesh::Part* const part = meta.get_part(elem_blocks[i]->name());
            assert(part != NULL);
            stk::io::define_io_fields(elem_blocks[i], Ioss::Field::TRANSIENT,
				      *part, part_primary_entity_rank(*part));
          }
        }
      }

      void define_input_nodeset_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::NodeSetContainer& nodesets = region.get_nodesets();
        for(size_t i=0; i < nodesets.size(); i++) {
          if (stk::io::include_entity(nodesets[i])) {
            stk::mesh::Part* const part = meta.get_part(nodesets[i]->name());
            assert(part != NULL);
            stk::io::define_io_fields(nodesets[i], Ioss::Field::TRANSIENT,
				      *part, part_primary_entity_rank(*part));
          }
        }
      }

      void define_input_sideset_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        if (meta.spatial_dimension() <= meta.side_rank())
          return;

        const Ioss::SideSetContainer& side_sets = region.get_sidesets();
        for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
            it != side_sets.end(); ++it) {
          Ioss::SideSet *entity = *it;
          if (stk::io::include_entity(entity)) {
            const Ioss::SideBlockContainer& blocks = entity->get_side_blocks();
            for(size_t i=0; i < blocks.size(); i++) {
              if (stk::io::include_entity(blocks[i])) {
                stk::mesh::Part* const part = meta.get_part(blocks[i]->name());
                assert(part != NULL);
                stk::io::define_io_fields(blocks[i], Ioss::Field::TRANSIENT,
					  *part, part_primary_entity_rank(*part));
              }
            }
          }
        }
      }

    }

    // ========================================================================
    // Iterate over all Ioss entities in the input mesh database and
    // define a stk_field for all transient fields found.  The stk
    // field will have the same name as the field on the database.
    //
    // Note that all fields found on the database will have a
    // corresponding stk field defined.  If you want just a selected
    // subset of the defined fields, you will need to define the
    // fields manually.
    //
    // To populate the stk field with data from the database, call
    // process_input_request().
    void MeshData::define_input_fields()
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_input_region),
		       "There is no Mesh Input Region associated with this Mesh Data.");
      Ioss::Region *region = m_input_region.get();
      define_input_nodeblock_fields(*region, meta_data());
      define_input_elementblock_fields(*region, meta_data());
      define_input_nodeset_fields(*region, meta_data());
      define_input_sideset_fields(*region, meta_data());
    }

    // ========================================================================
    // Iterate over all fields defined in the stk mesh data structure.
    // If the field has the io_attribute set, then define that field
    // on the corresponding io entity on the output mesh database.
    // The database field will have the same name as the stk field.
    //
    // To export the data to the database, call
    // process_output_request().

    void MeshData::define_output_fields(bool add_all_fields)
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_output_region),
		       "There is no Results output database associated with this Mesh Data.");

      Ioss::Region *region = m_output_region.get();
      region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

      // Special processing for nodeblock (all nodes in model)...
      stk::io::ioss_add_fields(meta_data().universal_part(), stk::mesh::MetaData::NODE_RANK,
			       region->get_node_blocks()[0],
			       Ioss::Field::TRANSIENT, add_all_fields);

      const stk::mesh::PartVector &all_parts = meta_data().get_parts();
      for ( stk::mesh::PartVector::const_iterator
	      ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

	stk::mesh::Part * const part = *ip;

	// Check whether this part should be output to results database.
	if (stk::io::is_part_io_part(*part)) {
	  stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
	  // Get Ioss::GroupingEntity corresponding to this part...
	  Ioss::GroupingEntity *entity = region->get_entity(part->name());
	  if (entity != NULL) {
	    stk::io::ioss_add_fields(*part, rank, entity, Ioss::Field::TRANSIENT, add_all_fields);
	  }

	  // If rank is != NODE_RANK, then see if any fields are defined on the nodes of this part
	  // (should probably do edges and faces also...)
	  // Get Ioss::GroupingEntity corresponding to the nodes on this part...
	  if (rank != stk::mesh::MetaData::NODE_RANK) {
	    Ioss::GroupingEntity *node_entity = NULL;
	    if (use_nodeset_for_part_nodes_fields()) {
	      std::string nodes_name = part->name() + "_nodes";
	      node_entity = region->get_entity(nodes_name);
	    } else {
	      node_entity = region->get_entity("nodeblock_1");
	    }
	    if (node_entity != NULL) {
	      stk::io::ioss_add_fields(*part, stk::mesh::MetaData::NODE_RANK,
				       node_entity, Ioss::Field::TRANSIENT, add_all_fields);
	    }
	  }
	}
      }
      region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
    }
    // ========================================================================
    void MeshData::process_input_request(double time)
    {
      // Find the step on the database with time closest to the requested time...
      Ioss::Region *region = m_input_region.get();
      int step_count = region->get_property("state_count").get_int();
      double delta_min = 1.0e30;
      int    step_min  = 0;
      for (int istep = 0; istep < step_count; istep++) {
        double state_time = region->get_state_time(istep+1);
        double delta = state_time - time;
        if (delta < 0.0) delta = -delta;
        if (delta < delta_min) {
          delta_min = delta;
          step_min  = istep;
          if (delta == 0.0) break;
        }
      }
      // Exodus steps are 1-based;
      process_input_request(step_min+1);
    }

    void MeshData::process_input_request(int step)
    {
      if (step <= 0)
        return;

      Ioss::Region *region = m_input_region.get();
      if (region) {
        // Pick which time index to read into solution field.
        region->begin_state(step);

        input_nodeblock_fields(*region, bulk_data());
        input_elementblock_fields(*region, bulk_data());
        input_nodeset_fields(*region, bulk_data());
        input_sideset_fields(*region, bulk_data());

        region->end_state(step);

      } else {
        std::cerr << "INTERNAL ERROR: Mesh Input Region pointer is NULL in process_input_request.\n";
        std::exit(EXIT_FAILURE);
      }
    }

    double MeshData::process_restart_input(double time)
    {
      // Find the step on the database with time closest to the requested time...
      Ioss::Region *region = m_input_region.get();
      int step_count = region->get_property("state_count").get_int();
      if (step_count == 0) {
        std::ostringstream msg ;
        msg << " ERROR: Restart database '" << region->get_database()->get_filename()
	    << " has no transient data.";
        throw std::runtime_error( msg.str() );
      }
      
      double delta_min = 1.0e30;
      int    step_min  = 0;
      for (int istep = 0; istep < step_count; istep++) {
        double state_time = region->get_state_time(istep+1);
        double delta = state_time - time;
        if (delta < 0.0) delta = -delta;
        if (delta < delta_min) {
          delta_min = delta;
          step_min  = istep;
          if (delta == 0.0) break;
        }
      }
      // Exodus steps are 1-based;
      return process_restart_input(step_min+1);
    }

    double MeshData::process_restart_input(int step)
    {
      if (step <= 0)
        return 0;

      int missing_fields = 0;
      ThrowErrorMsgIf (Teuchos::is_null(m_input_region),
                       "There is no Input mesh/restart region associated with this Mesh Data.");

      Ioss::Region *region = m_input_region.get();

      // Pick which time index to read into solution field.
      region->begin_state(step);

      // Special processing for nodeblock (all nodes in model)...
      missing_fields += ioss_restore_restart_fields(bulk_data(), meta_data().universal_part(),
						    stk::mesh::MetaData::NODE_RANK,
						    region->get_node_blocks()[0]);

      // Now handle all non-nodeblock parts...
      const stk::mesh::PartVector &all_parts = meta_data().get_parts();
      for ( stk::mesh::PartVector::const_iterator
	      ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

        stk::mesh::Part * const part = *ip;

        // Check whether this part should be output to results database.
        if (stk::io::is_part_io_part(*part)) {
          stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
          // Get Ioss::GroupingEntity corresponding to this part...
          Ioss::GroupingEntity *entity = region->get_entity(part->name());
          if (entity != NULL && entity->type() != Ioss::SIDESET) {
            missing_fields += ioss_restore_restart_fields(bulk_data(), *part, rank, entity);
          }

          // If rank is != NODE_RANK, then see if any fields are defined on the nodes of this part
          // (should probably do edges and faces also...)
          // Get Ioss::GroupingEntity corresponding to the nodes on this part...
          if (rank != stk::mesh::MetaData::NODE_RANK) {
	    Ioss::GroupingEntity *node_entity = NULL;
	    if (use_nodeset_for_part_nodes_fields()) {
	      std::string nodes_name = part->name() + "_nodes";
	      node_entity = region->get_entity(nodes_name);
	    } else {
	      node_entity = region->get_entity("nodeblock_1");
	    }
	    if (node_entity != NULL) {
	      missing_fields += ioss_restore_restart_fields(bulk_data(), *part,
							    stk::mesh::MetaData::NODE_RANK, node_entity);
            }
          }
        }
      }

      if (missing_fields > 0) {
        std::ostringstream msg ;
        msg << "ERROR: Restart read at step " << step << " could not find " << missing_fields << " fields.\n";
        throw std::runtime_error( msg.str() );
      }

      region->end_state(step);

      // return time
      return region->get_state_time(step);
    }
  } // namespace io
} // namespace stk
