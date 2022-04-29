// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_RegionInterface_h
#define Akri_RegionInterface_h

#include <stk_util/util/ReportHandler.hpp>
#include <memory>

namespace Ioss { class Region; }
namespace stk { namespace diag { class Timer; } }
namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class RegionInterface {

 public:
  static RegionInterface & set_currently_parsed_region(stk::mesh::MetaData & meta,
      const std::string & region_name,
      stk::diag::Timer & region_timer,
      const std::string & name_of_input_mesh,
      Ioss::Region * ioss_region)
  {
    the_currently_parsed_region = std::make_unique<RegionInterface>(meta, region_name, region_timer, name_of_input_mesh, ioss_region);
    return *the_currently_parsed_region;
  }
  static void clear_currently_parsed_region() { the_currently_parsed_region.reset(); }
  static RegionInterface & get_currently_parsed_region() { ThrowAssert(nullptr != the_currently_parsed_region); return *the_currently_parsed_region; }

  stk::mesh::MetaData & get_stk_mesh_meta_data() { return my_meta; }
  const std::string & name() { return my_region_name; }
  stk::diag::Timer & getRegionTimer() const { return my_region_timer; }
  const std::string & name_of_input_mesh() const { return my_name_of_input_mesh; }
  Ioss::Region * get_input_io_region() { return my_ioss_region; }

  // must be public to be used by make_unique
  RegionInterface(stk::mesh::MetaData & meta,
      const std::string & region_name,
      stk::diag::Timer & region_timer,
      const std::string & input_mesh,
      Ioss::Region * ioss_region)
  : my_meta(meta),
    my_region_name(region_name),
    my_region_timer(region_timer),
    my_name_of_input_mesh(input_mesh),
    my_ioss_region(ioss_region) {}

private:
  static std::unique_ptr<RegionInterface> the_currently_parsed_region;
  stk::mesh::MetaData & my_meta;
  const std::string my_region_name;
  stk::diag::Timer & my_region_timer;
  const std::string my_name_of_input_mesh;
  Ioss::Region * my_ioss_region;
};

} // end krino namespace

#endif /* Akri_RegionInterface_h */
