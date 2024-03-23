// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_ResultsOutputOptions_h
#define Akri_ResultsOutputOptions_h

#include <string>
#include <utility>
#include <vector>
#include <stk_util/environment/Scheduler.hpp>
#include <stk_util/diag/String.hpp>     // for String
#include "stk_topology/topology.hpp"    // for topology, etc
#include <Ioss_PropertyManager.h>

namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class MetaData; } }

namespace krino {

  typedef std::pair<std::string, std::string> FieldName_OutputName_Pair;

  class ResultsOutputOptions {

  public:

    ResultsOutputOptions() :
      my_scheduler(),
      my_numStepIncrements(0)
    {
      my_scheduler.set_lookahead(1);
      my_filename = "default_filename";
    }

    stk::util::Scheduler & get_scheduler() { return my_scheduler; }

    void set_title(std::string title) { my_title = title; }
    const std::string & get_title() const { return my_title; }

    void set_name(std::string theName) { my_name = theName; }
    const std::string & get_name() const { return my_name; }

    void set_filename(std::string filename) { my_filename = filename; }
    const std::string & get_filename() const { return my_filename; }

    void add_step_increment(int start, int increment) {
      my_scheduler.add_interval(start, increment);
      my_numStepIncrements++;
    }
    unsigned get_num_step_increments() const { return my_numStepIncrements; }

    void add_nodal_field(const std::string & internalName, const std::string & newName)
    {
      FieldName_OutputName_Pair aPair(internalName, newName);
      my_nodal_fields.insert(aPair);
    }
    const std::set<FieldName_OutputName_Pair> & get_nodal_fields() const { return my_nodal_fields; }

    void add_element_field(const std::string & internalName, const std::string & newName)
    {
      FieldName_OutputName_Pair aPair(internalName, newName);
      my_element_fields.insert(aPair);
    }
    const std::set<FieldName_OutputName_Pair> & get_element_fields() const { return my_element_fields; }

    void add_property(const Ioss::Property & property) { my_properties.add(property); }
    Ioss::PropertyManager & get_properties() { return my_properties; }

  private:
    stk::util::Scheduler my_scheduler;
    unsigned my_numStepIncrements;
    std::string my_name;
    std::string my_title;
    std::string my_filename;
    std::set<FieldName_OutputName_Pair> my_nodal_fields;
    std::set<FieldName_OutputName_Pair> my_element_fields;
    Ioss::PropertyManager my_properties;
  };

} // namespace krino

#endif /* Akri_ResultsOutputOptions_h */
