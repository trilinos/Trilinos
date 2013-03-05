/*
 * Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *           
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *                         
 * * Neither the name of Sandia Corporation nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *                                                 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <cstring>
#include <time.h>

#include "EML_CodeTypes.h"
#include "EML_SystemInterface.h"

#include <Ioss_CodeTypes.h>
#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>
#include <Ioss_SurfaceSplit.h>

#include "add_to_log.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#define OUTPUT std::cerr
#define FOUTPUT std::cout

// ========================================================================
namespace {
  int file_info(const std::string& inpfile, const std::string& input_type,
		SystemInterface& interface);

  void output_names(const std::string &type, const Ioss::NameList &fields,
		    Ioss::GroupingEntity *entity)
  {
    OUTPUT << "\n" << type << " variables on exodus data base:\n";
    Ioss::NameList::const_iterator IF;
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      const Ioss::VariableType *var_type = entity->get_field(field_name).raw_storage();
      OUTPUT << "\t" << field_name << "\t" << var_type->name() << std::endl;
    }
  }
}
// ========================================================================

namespace {
  std::string codename;
  std::string version = "1.0";
}

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  time_t begin_time = time(NULL);
  std::string in_type = "exodusII";

  codename = argv[0];
  size_t ind = codename.find_last_of("/", codename.size());
  if (ind != std::string::npos)
    codename = codename.substr(ind+1, codename.size());

  SystemInterface::show_version();
  Ioss::Init::Initializer io;

  SystemInterface interface;
  bool ok = interface.parse_options(argc, argv);

  if (ok) {
    std::string in_file = interface.input_file();
    std::string output_file = interface.output_file();

    OUTPUT << "Input:    '" << in_file  << "', Type: " << in_type  << '\n';
    OUTPUT << "Output:   '" << output_file << "', Type: matlab script\n\n";

    ok = file_info(in_file, in_type, interface);
  }

  std::string success = ok ? "successful" : "unsuccessful";
  OUTPUT << "\n" << codename << " execution " << success << ".\n";
  time_t end_time = time(NULL);
  add_to_log(codename.c_str(), (int)(end_time - begin_time));
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

namespace {
  int file_info(const std::string& inpfile, const std::string& input_type, SystemInterface &interface)
  {
    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART,
						    (MPI_Comm)MPI_COMM_WORLD);
    if (dbi == NULL || !dbi->ok(true)) {
      return false;
    }

    dbi->set_field_separator(interface.field_suffix());
    dbi->set_lower_case_variable_names(false);
    
    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");

    if (interface.list_vars()) {
      StringIdVector types_to_list = interface.vars_to_list();
      for (size_t i=0; i < types_to_list.size(); i++) {
	std::string type = types_to_list[i].first;

	if (type == "all" || type == "global") {
	  Ioss::NameList fields;
	  region.field_describe(Ioss::Field::TRANSIENT, &fields);
	  output_names("Global", fields, &region);
	}
	if (type == "all" || type == "nodal") {
	  Ioss::NameList fields;
	  Ioss::NodeBlock *nb = region.get_node_blocks()[0];
	  nb->field_describe(Ioss::Field::TRANSIENT, &fields);
	  output_names("Nodal", fields, nb);
	}
      }
      return true;
    }

    Ioss::NameList fields;
    StringIdVector global_vars = interface.global_var_names();
    if (global_vars.size() > 0) {
      if (global_vars[0].first == "all") {
	region.field_describe(Ioss::Field::TRANSIENT, &fields);
      } else if (global_vars[0].first == "none") {
	; // do nothing.  This will be used when nodal, element, ... supported
      } else {
	for (size_t i=0; i < global_vars.size(); i++) {
	  std::string field_name = global_vars[i].first;
	  if (region.field_exists(field_name)) {
	    fields.push_back(field_name);
	  } else {
	    OUTPUT << "WARNING: Global variable named '" << field_name << "' does not exist; it will be skipped.\n";
	  }
	}
      }
    } else {
      region.field_describe(Ioss::Field::TRANSIENT, &fields);
    }

    if (fields.size() == 0) {
      OUTPUT << "No variables selected; no output will be written\n";
      return false;
    }
    
    std::ofstream out_stream;
    out_stream.open(interface.output_file().c_str());

    out_stream.setf(std::ios::scientific);
    out_stream.setf(std::ios::showpoint);

    out_stream << "% number of curves\nnvars = " << fields.size()+1 << "\n";

    size_t namelen = 4; // length of 'time'
    Ioss::NameList::const_iterator IF;
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      if (field_name.length() > namelen)
	namelen = field_name.length();
    }

    out_stream << "names= [\n";
    out_stream << "'" << std::left << std::setw(namelen) << "TIME" << "';" << std::endl;

    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      out_stream << "'" << std::left << std::setw(namelen) << field_name << "';" << std::endl;
    }
    out_stream << "];\n";
      
    // Get number of timesteps...
    int num_steps = 0;
    if (region.property_exists("state_count") && region.get_property("state_count").get_int() > 0) {
      num_steps = region.get_property("state_count").get_int();
    } else {
      out_stream << "GENESIS file -- no time steps written\n";
      return false;
    }

    // ========================================================================
    // Calculate min and max times to extract data...
    int st_min = 1;
    int st_max = num_steps;
    
    double tmin = interface.minimum_time();
    double tmax = interface.maximum_time();
    if (tmax == -1.0) tmax = region.get_max_time().second;
    
    if (tmin > 0.0 || tmax != region.get_max_time().second) {
      st_min = 0;

      // Count number of active steps...
      for (int i=1; i <= num_steps; i++) {
	double time = region.get_state_time(i);
	if (time < tmin) {
	  st_min = i;
	}
	if (time <= tmax) {
	  st_max = i;
	}
      }
      st_min++;
      num_steps = st_max - st_min + 1;
    }
    
    // ========================================================================
    // Output time values...

    out_stream << "TIME = zeros(1, " << num_steps << ");\n";
    out_stream << "TIME= [\n";
    for (int i=st_min; i <= st_max; i++) {
      double time = region.get_state_time(i);
      out_stream << " " << time << "\n";
    }
    out_stream << "];\n";

    // ========================================================================
    // Output field values...

    std::vector<double> data(1);
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      int comp_count = region.get_field(field_name).raw_storage()->component_count();
      out_stream << field_name << " = zeros(" << comp_count << ", " << num_steps << ");\n";
      out_stream << field_name << "= [\n";
      for (int istep=st_min; istep <= st_max; istep++) {
	region.begin_state(istep);
	region.get_field_data(field_name, data);
	for (int i=0; i < comp_count; i++) {
	  out_stream << " " << data[i];
	}
	out_stream << ";\n";
	region.end_state(istep);
      }
      out_stream << "];\n";
    }
    OUTPUT << "Wrote data for " << fields.size() << " variables at "
	   << st_max - st_min + 1 << " timesteps.\n";
    return true;
  }
}
