// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "IossRegionReport.h"
#include "Ioss_Assembly.h"
#include "Ioss_Blob.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_Property.h"
#include "Ioss_SubSystem.h"

namespace ioss_region_report {

  std::string sep1 = " ";
  std::string sep2 = " | ";

  std::ostream &operator<<(std::ostream &os, const Messages &messages)
  {
    for (auto msg : messages.messages)
      os << msg << std::endl;
    return os;
  }

  Messages grouping_entity_report(bool report_transient, const Ioss::GroupingEntity &entity);

  Messages entity_block_report(bool report_transient, const Ioss::EntityBlock &entity);
  Messages entity_set_report(bool report_transient, const Ioss::EntitySet &entity);

  Messages node_block_report(bool report_transient, const Ioss::NodeBlock &entity);
  Messages elem_block_report(bool report_transient, const Ioss::ElementBlock &entity);
  Messages face_block_report(bool report_transient, const Ioss::FaceBlock &entity);
  Messages edge_block_report(bool report_transient, const Ioss::EdgeBlock &entity);

  Messages node_set_report(bool report_transient, const Ioss::NodeSet &entity);
  Messages elem_set_report(bool report_transient, const Ioss::ElementSet &entity);
  Messages face_set_report(bool report_transient, const Ioss::FaceSet &entity);
  Messages edge_set_report(bool report_transient, const Ioss::EdgeSet &entity);

  Messages side_set_report(bool report_transient, const Ioss::SideSet &entity);
  Messages comm_set_report(bool report_transient, const Ioss::CommSet &entity);
  Messages struct_block_report(bool report_transient, const Ioss::StructuredBlock &entity);
  Messages assembly_report(bool report_transient, const Ioss::Assembly &entity);
  Messages blob_report(bool report_transient, const Ioss::Blob &entity);
  Messages coord_frame_report(bool report_transient, const Ioss::CoordinateFrame &entity);

  Messages    property_report(const Ioss::Property &property);
  Messages    field_report(const Ioss::Field &field);
  std::string field_role(const int id);
  std::string property_type(const int id);

  Messages entity_block_report(bool report_transient, const Ioss::EntityBlock &entity)
  {
    Messages msgs;
    msgs += "EntityBlock" + sep2;
    msgs.begin = "\t";
    msgs += grouping_entity_report(report_transient, entity);
    return msgs;
  }

  Messages entity_set_report(bool report_transient, const Ioss::EntitySet &entity)
  {
    Messages msgs;
    msgs += "EntitySet" + sep2;
    msgs.begin = "\t";
    msgs += grouping_entity_report(report_transient, entity);
    return msgs;
  }

  Messages side_set_report(bool report_transient, const Ioss::SideSet &entity)
  {
    Messages msgs;
    msgs += "SideSet" + sep2;
    msgs.begin = "\t";
    msgs += grouping_entity_report(report_transient, entity);
    return msgs;
  }

  Messages comm_set_report(bool report_transient, const Ioss::CommSet &entity)
  {
    Messages msgs;
    msgs += "CommSet" + sep2;
    msgs.begin = "\t";
    msgs += grouping_entity_report(report_transient, entity);
    return msgs;
  }

  Messages struct_block_report(bool report_transient, const Ioss::StructuredBlock &entity)
  {
    Messages msgs;
    msgs += "StructuredBlock" + sep2;
    msgs.begin = "\t";
    msgs += entity_block_report(report_transient, entity);
    return msgs;
  }

  Messages assembly_report(bool report_transient, const Ioss::Assembly &entity)
  {
    Messages msgs;
    msgs += "Assembly" + sep2;
    msgs.begin = "\t";
    msgs += grouping_entity_report(report_transient, entity);
    return msgs;
  }

  Messages blob_report(bool report_transient, const Ioss::Blob &entity)
  {
    Messages msgs;
    msgs += "Blob" + sep2;
    msgs.begin = "\t";
    msgs += grouping_entity_report(report_transient, entity);
    return msgs;
  }

  Messages coord_frame_report(bool report_transient, const Ioss::CoordinateFrame &entity)
  {
    Messages msgs;
    msgs += "CoordinateFrame" + sep2;
    msgs.begin = "\t";
    msgs += sep2 + "Id" + sep1 + std::to_string(entity.id());
    return msgs;
  }

  Messages node_block_report(bool report_transient, const Ioss::NodeBlock &entity)
  {
    Messages msgs;
    msgs += "NodeBlock" + sep2;
    msgs.begin = "\t";
    msgs += entity_block_report(report_transient, entity);
    return msgs;
  }

  Messages elem_block_report(bool report_transient, const Ioss::ElementBlock &entity)
  {
    Messages msgs;
    msgs += "ElemBlock" + sep2;
    msgs.begin = "\t";
    msgs += entity_block_report(report_transient, entity);
    return msgs;
  }

  Messages face_block_report(bool report_transient, const Ioss::FaceBlock &entity)
  {
    Messages msgs;
    msgs += "FaceBlock" + sep2;
    msgs.begin = "\t";
    msgs += entity_block_report(report_transient, entity);
    return msgs;
  }

  Messages edge_block_report(bool report_transient, const Ioss::EdgeBlock &entity)
  {
    Messages msgs;
    msgs += "EdgeBlock" + sep2;
    msgs.begin = "\t";
    msgs += entity_block_report(report_transient, entity);
    return msgs;
  }

  Messages node_set_report(bool report_transient, const Ioss::NodeSet &entity)
  {
    Messages msgs;
    msgs += "NodeSet" + sep2;
    msgs.begin = "\t";
    msgs += entity_set_report(report_transient, entity);
    return msgs;
  }

  Messages elem_set_report(bool report_transient, const Ioss::ElementSet &entity)
  {
    Messages msgs;
    msgs += "ElemSet" + sep2;
    msgs.begin = "\t";
    msgs += entity_set_report(report_transient, entity);
    return msgs;
  }

  Messages face_set_report(bool report_transient, const Ioss::FaceSet &entity)
  {
    Messages msgs;
    msgs += "FaceSet" + sep2;
    msgs.begin = "\t";
    msgs += entity_set_report(report_transient, entity);
    return msgs;
  }

  Messages edge_set_report(bool report_transient, const Ioss::EdgeSet &entity)
  {
    Messages msgs;
    msgs += "EdgeSet" + sep2;
    msgs.begin = "\t";
    msgs += entity_set_report(report_transient, entity);
    return msgs;
  }

  Messages grouping_entity_report(bool report_transient, const Ioss::GroupingEntity &entity)
  {
    auto filename    = entity.get_filename();
    auto state       = entity.get_database()->get_region()->get_current_state();
    auto entity_name = entity.get_property("name").get_string();

    Messages msgs;

    {
      Message msg = "GroupingEntity" + sep1;
      if (!entity_name.empty())
        msg += sep2 + entity_name;
      msg += sep2 + "State" + sep1 + std::to_string(state);
      msgs += msg;
    }

    {
      msgs.begin = "\t";

      // Transient Region Fields
      std::vector<std::string>           names;
      std::vector<Ioss::Field::RoleType> valid_role_types;
      if (report_transient) {
        valid_role_types.push_back(Ioss::Field::RoleType::REDUCTION);
        valid_role_types.push_back(Ioss::Field::RoleType::TRANSIENT);
      }
      else {
        // Properties
        {
          names.clear();
          entity.property_describe(&names);
          for (auto name : names)
            msgs += property_report(entity.get_property(name));
        }

        valid_role_types.push_back(Ioss::Field::RoleType::INTERNAL);
        valid_role_types.push_back(Ioss::Field::RoleType::MESH);
        valid_role_types.push_back(Ioss::Field::RoleType::ATTRIBUTE);
        valid_role_types.push_back(Ioss::Field::RoleType::COMMUNICATION);
        valid_role_types.push_back(Ioss::Field::RoleType::MESH_REDUCTION);
      }

      // Fields
      for (auto role_type : valid_role_types) {
        names.clear();
        entity.field_describe(role_type, &names);
        for (auto name : names)
          msgs += field_report(entity.get_field(name));
      }
    }

    return msgs;
  }

  Messages region_report(const Ioss::Region &region)
  {
    auto filename         = region.get_filename();
    auto state            = region.get_current_state();
    auto region_name      = region.get_property("name").get_string();
    auto report_transient = state > 0 ? true : false;

    Messages msgs;

    {
      Message msg = "Region" + sep1 + filename;
      if (!region_name.empty())
        msg += sep2 + region_name;
      msg += sep2 + "State" + sep1 + std::to_string(state);
      msgs += msg;
    }

    {
      msgs.begin = "\t";
      msgs += grouping_entity_report(report_transient, region);

      for (auto block : region.get_node_blocks())
        msgs += node_block_report(report_transient, *block);
      for (auto block : region.get_element_blocks())
        msgs += elem_block_report(report_transient, *block);
      for (auto block : region.get_face_blocks())
        msgs += face_block_report(report_transient, *block);
      for (auto block : region.get_edge_blocks())
        msgs += edge_block_report(report_transient, *block);

      for (auto set : region.get_nodesets())
        msgs += node_set_report(report_transient, *set);
      for (auto set : region.get_elementsets())
        msgs += elem_set_report(report_transient, *set);
      for (auto set : region.get_facesets())
        msgs += face_set_report(report_transient, *set);
      for (auto set : region.get_edgesets())
        msgs += edge_set_report(report_transient, *set);

      for (auto set : region.get_sidesets())
        msgs += side_set_report(report_transient, *set);
      for (auto set : region.get_commsets())
        msgs += comm_set_report(report_transient, *set);
      for (auto set : region.get_structured_blocks())
        msgs += struct_block_report(report_transient, *set);
      for (auto set : region.get_assemblies())
        msgs += assembly_report(report_transient, *set);
      for (auto set : region.get_blobs())
        msgs += blob_report(report_transient, *set);
      for (auto set : region.get_coordinate_frames())
        msgs += coord_frame_report(report_transient, set);
    }
    return msgs;
  }

  Messages property_report(const Ioss::Property &property)
  {
    Messages msgs;

    {
      Message msg;
      msg += "Property" + sep1 + property.get_name();
      msg += sep2 + "PropertyType" + sep1 + property_type(property.get_type());
      if (property.get_type() == Ioss::Property::REAL) {
        msg += sep2 + "PropertyValue" + sep1 + std::to_string(property.get_real());
      }
      else if (property.get_type() == Ioss::Property::STRING) {
        msg += sep2 + "PropertyValue" + sep1 + property.get_string();
      }
      else if (property.get_type() == Ioss::Property::INTEGER) {
        msg += sep2 + "PropertyValue" + sep1 + std::to_string(property.get_int());
      }
      msgs += msg;
    }

    return msgs;
  }

  std::string property_type(const int id)
  {

    if (id == -1) {
      return std::string("INVALID");
    }

    std::vector<std::string> property_type_string{"REAL",   "INTEGER",     "POINTER",
                                                  "STRING", "VEC_INTEGER", "VEC_DOUBLE"};

    return property_type_string[id];
  }

  Messages field_report(const Ioss::Field &field)
  {
    Messages messages;

    {
      Message msg;
      msg += "Field" + sep1 + field.get_name();
      msg += sep2 + "RoleType" + sep1 + field_role(field.get_role());
      msg += sep2 + "NumItemsInField" + sep1 + std::to_string(field.raw_count());
      msg += sep2 + "SizeInBytes" + sep1 + std::to_string(field.get_size());
      messages += msg;
    }

    return messages;
  }

  std::string field_role(const int id)
  {

    std::vector<std::string> field_role_string{"INTERNAL",      "MESH",           "ATTRIBUTE",
                                               "COMMUNICATION", "MESH_REDUCTION", "REDUCTION",
                                               "TRANSIENT"};

    return field_role_string[id];
  }

} // End namespace ioss_region_report
