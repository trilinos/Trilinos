// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Iodw_Utils.h>
#include <Ioss_Region.h>

#include <Ioss_EdgeBlock.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_FaceBlock.h>
#include <Ioss_NodeBlock.h>

#include <Ioss_EdgeSet.h>
#include <Ioss_ElementSet.h>
#include <Ioss_FaceSet.h>
#include <Ioss_NodeSet.h>

#include <Ioss_DatabaseIO.h>

#include <iomanip>
#include <string>

namespace Iodw {
  namespace Utils {

    void msg(const std::string &m) { std::cerr << m << "\n"; };
    void msg(int m) { std::cerr << m << "\n"; };
    void msg(int64_t m) { std::cerr << m << "\n"; };
    void msg(double m) { std::cerr << m << "\n"; };
    void msg(float m) { std::cerr << m << "\n"; };

    auto print_Property = [](Ioss::Property p) {
      auto type = p.get_type();
      if (type == Ioss::Property::BasicType::STRING) {
        msg(p.get_string());
      }
      else if (type == Ioss::Property::BasicType::INTEGER) {
        msg(p.get_int());
      }
      else if (type == Ioss::Property::BasicType::REAL) {
        msg(p.get_real());
      }
      // else if( type == Ioss::Property::BasicType::POINTER ) { msg(p.get_pointer()); }
      else if (type == Ioss::Property::BasicType::POINTER) {
        msg("pointer");
      }
      else if (type == Ioss::Property::BasicType::INVALID) {
        msg("INVALID");
      }
      else {
        msg("UNKNOWN type");
      }
    };

    auto print_Field = [](Ioss::Field f) {
      size_t size = f.get_size();
      std::cerr << f.raw_count() << "\n";
    };

    void IossToDW::operator()(Ioss::DatabaseIO *dbi)
    {
      auto process_GroupingEntity = [](Ioss::Region *r, Ioss::GroupingEntity *entity) {
        msg("process GroupingEntity");
        std::cerr << "\t" << entity->type_string() << ", " << entity->name() << "\n";
        {
          std::vector<std::string> descr;
          entity->property_describe(&descr);
          msg("\t\tproperties:");
          for (auto name : descr) {
            std::cerr << "\t\t\t" << r->name() << " * " << entity->name() << " * " << name << "\n";
            // print_Property(entity->get_property(name));
          }
        }
        {
          std::vector<std::string> descr;
          entity->field_describe(&descr);
          entity->field_describe(&descr);
          msg("\t\tfields:");
          for (auto name : descr) {
            std::cerr << "\t\t\t" << r->name() << " * " << entity->name() << " * " << name << "\n";
            // print_Field(entity->get_field(name));
          }
        }
      };

      auto process_EntityBlock = [process_GroupingEntity](Ioss::Region *     r,
                                                          Ioss::EntityBlock *entity) {
        msg("process_block");
        process_GroupingEntity(r, entity);
      };

      auto process_EntitySet = [process_GroupingEntity](Ioss::Region *r, Ioss::EntitySet *entity) {
        msg("process_set");
        process_GroupingEntity(r, entity);
      };

      RegionKeys keys;

      std::cerr << "IossToDW::operator(): " << dbi->get_filename() << "\n";
      std::string   region_name("region");
      Ioss::Region *region = new Ioss::Region(dbi, region_name);

      for (auto entity : region->get_edge_blocks())
        process_EntityBlock(region, entity);

      for (auto entity : region->get_element_blocks())
        process_EntityBlock(region, entity);

      for (auto entity : region->get_face_blocks())
        process_EntityBlock(region, entity);

      for (auto entity : region->get_node_blocks())
        process_EntityBlock(region, entity);

      for (auto entity : region->get_edgesets())
        process_EntitySet(region, entity);

      for (auto entity : region->get_elementsets())
        process_EntitySet(region, entity);

      for (auto entity : region->get_facesets())
        process_EntitySet(region, entity);

      for (auto entity : region->get_nodesets())
        process_EntitySet(region, entity);
    }

  } // namespace Utils
} // namespace Iodw
