// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Interface_Name_Generator_h
#define Akri_Interface_Name_Generator_h
//

#include <string>

namespace krino {

// Interface_Name_Generator are used in Phase_Support::decompose_blocks() to determine the name of the
// interface sideset within a block so that cdfem death and LS-like problems can have different naming conventions
// that are more user-friendly for each case.
class Interface_Name_Generator
{
public:
  virtual ~Interface_Name_Generator() {}
  virtual std::string interface_name(const std::string & io_part_name, const std::string & phase_name) const = 0;
  virtual std::string interface_superset_name(const std::string & phase_name) const = 0;
  virtual std::string interface_subset_name(const std::string & io_part_name, const std::string & phase_name) const = 0;
};

class LS_Name_Generator : public Interface_Name_Generator
{
public:
  virtual ~LS_Name_Generator() {}
  virtual std::string interface_name(const std::string & io_part_name, const std::string & phase_name) const
  {
    return "surface_"+io_part_name+"_"+phase_name;
  }
  virtual std::string interface_superset_name(const std::string & phase_name) const
  {
    return "surface_"+phase_name;
  }
  virtual std::string interface_subset_name(const std::string & io_part_name, const std::string & phase_name) const
  {
    return "surface_"+io_part_name+"_"+phase_name;
  }
};

class Death_Name_Generator : public Interface_Name_Generator
{
public:
  Death_Name_Generator(std::string spec_name) : death_spec_name(spec_name) {}
  virtual ~Death_Name_Generator() {}
  virtual std::string interface_name(const std::string & io_part_name, const std::string & /*phase_name*/) const
  {
    return "surface_"+io_part_name+"_"+death_spec_name;
  }
  virtual std::string interface_superset_name(const std::string & /*phase_name*/) const
  {
    return "surface_"+death_spec_name;
  }
  virtual std::string interface_subset_name(const std::string & io_part_name, const std::string & /*phase_name*/) const
  {
    return "surface_"+io_part_name+"_"+death_spec_name;
  }
private:
  std::string death_spec_name;
};


} // namespace krino

#endif // Akri_Interface_Name_Generator_h
