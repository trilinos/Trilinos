// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

    void msg(const std::string &m) { std::cout << m << std::endl; };
    void msg(int m) { std::cout << m << std::endl; };
    void msg(int64_t m) { std::cout << m << std::endl; };
    void msg(double m) { std::cout << m << std::endl; };
    void msg(float m) { std::cout << m << std::endl; };

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
      std::cout << f.raw_count() << std::endl;
    };

    void IossToDW::operator()(Ioss::DatabaseIO *dbi)
    {
      auto process_GroupingEntity = [](Ioss::Region *r, Ioss::GroupingEntity *entity) {
        msg("process GroupingEntity");
        std::cout << "\t" << entity->type_string() << ", " << entity->name() << std::endl;
        {
          std::vector<std::string> descr;
          entity->property_describe(&descr);
          msg("\t\tproperties:");
          for (auto name : descr) {
            std::cout << "\t\t\t" << r->name() << " * " << entity->name() << " * " << name
                      << std::endl;
            // print_Property(entity->get_property(name));
          }
        }
        {
          std::vector<std::string> descr;
          entity->field_describe(&descr);
          entity->field_describe(&descr);
          msg("\t\tfields:");
          for (auto name : descr) {
            std::cout << "\t\t\t" << r->name() << " * " << entity->name() << " * " << name
                      << std::endl;
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

      std::cout << "IossToDW::operator(): " << dbi->get_filename() << std::endl;
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
