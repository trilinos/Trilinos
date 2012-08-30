// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
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

#include <Ioss_CodeTypes.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Region.h>
#include <Ioss_Utils.h>
#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <tokenize.h>

#include "Ioss_DBUsage.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_Property.h"
#include "Ioss_SerializeIO.h"
#include "Ioss_State.h"
#include "Ioss_SurfaceSplit.h"

namespace {
  void log_field(const char *symbol, const Ioss::GroupingEntity *entity,
		 const Ioss::Field &field, bool single_proc_only,
		 const Ioss::ParallelUtils &util);

#ifndef NDEBUG
  bool is_parallel_consistent(bool single_proc_only, const Ioss::GroupingEntity *ge,
			      const Ioss::Field &field, const Ioss::ParallelUtils &util)
  {
    if (single_proc_only)
      return true;

    std::string ge_name = ge->name();
    std::string field_name = field.get_name();
    unsigned int hash_code = Ioss::Utils::hash(ge_name) + Ioss::Utils::hash(field_name);
    unsigned int max_hash  = util.global_minmax(hash_code, Ioss::ParallelUtils::DO_MAX);
    unsigned int min_hash  = util.global_minmax(hash_code, Ioss::ParallelUtils::DO_MIN);
    if (max_hash != min_hash) {
      std::string errmsg = "Parallel inconsistency detected for field ";
      errmsg += field_name;
      errmsg += " on entity ";
      errmsg += ge_name;
      errmsg += "\n";
      IOSS_WARNING << errmsg;
      return false;
    } else {
      return true;
    }
  }
#endif
}

Ioss::DatabaseIO::DatabaseIO(Ioss::Region* region, const std::string& filename,
			     Ioss::DatabaseUsage db_usage,
			     MPI_Comm communicator,
			     const Ioss::PropertyManager &props)
  : properties(props), commonSideTopology(NULL), DBFilename(filename), dbState(STATE_INVALID),
    isParallel(false), isSerialParallel(false), myProcessor(0), cycleCount(0), overlayCount(0),
    fieldSuffixSeparator('_'), splitType(Ioss::SPLIT_BY_TOPOLOGIES),
    dbUsage(db_usage),dbIntSizeAPI(USE_INT32_API),
    nodeGlobalIdBackwardCompatibility(false), lowerCaseVariableNames(true),
    util_(communicator), region_(region), isInput(is_input_event(db_usage)),
    singleProcOnly(db_usage == WRITE_HISTORY || db_usage == WRITE_HEARTBEAT || Ioss::SerializeIO::isEnabled()),
    doLogging(false)
{
  isParallel  = util_.parallel_size() > 1;
  myProcessor = util_.parallel_rank();

  // Check environment variable IOSS_PROPERTIES. If it exists, parse
  // the contents and add to the 'properties' map.
  
  std::string env_props;
  if (util_.get_environment("IOSS_PROPERTIES", env_props, isParallel)) {
    // env_props string should be of the form
    // "PROP1=VALUE1:PROP2=VALUE2:..."
    std::vector<std::string> prop_val;
    Ioss::tokenize(env_props, ":", prop_val);
    
    for (size_t i=0; i < prop_val.size(); i++) {
      std::vector<std::string> property;
      Ioss::tokenize(prop_val[i], "=", property);
      if (property.size() != 2) {
	std::ostringstream errmsg;
	errmsg << "ERROR: Invalid property specification found in IOSS_PROPERTIES environment variable\n"
	       << "       Found '" << prop_val[i] << "' which is not of the correct PROPERTY=VALUE form";
	IOSS_ERROR(errmsg);
      }
      std::string prop = Ioss::Utils::uppercase(property[0]);
      std::string value = property[1];
      std::string up_value = Ioss::Utils::uppercase(value);
      bool all_digit = value.find_first_not_of("0123456789") == std::string::npos;

      if (myProcessor == 0) 
	std::cerr << "IOSS: Adding property '" << prop << "' with value '" << value << "'\n";
      
      if (all_digit) {
	int int_value = std::strtol(value.c_str(), NULL, 10);
	properties.add(Ioss::Property(prop, int_value));
      }
      else if (up_value == "TRUE" || up_value == "YES") {
	properties.add(Ioss::Property(prop, 1));
      }
      else if (up_value == "FALSE" || up_value == "NO") {
	properties.add(Ioss::Property(prop, 0));
      }
      else {
	properties.add(Ioss::Property(prop, value));
      }
    }
  }
}

Ioss::DatabaseIO::~DatabaseIO()
{
}

int Ioss::DatabaseIO::int_byte_size_api() const
{
  if (dbIntSizeAPI == Ioss::USE_INT32_API) {
    return 4;
  } else {
    return 8;
  }
}

void Ioss::DatabaseIO::set_int_byte_size_api(Ioss::DataSize size) const
{
  dbIntSizeAPI = size; // mutable
}

void Ioss::DatabaseIO::verify_and_log(const Ioss::GroupingEntity *ge, const Ioss::Field& field) const
{
  assert(is_parallel_consistent(singleProcOnly, ge, field, util_));
  if (get_logging()) {
    log_field(">", ge, field, singleProcOnly, util_);
  }
}

// Default versions do nothing...
bool Ioss::DatabaseIO::begin_state(Ioss::Region */* region */, int /* state */, double /* time */)
{
  return true;
}

bool Ioss::DatabaseIO::end_state(Ioss::Region */* region */, int /* state */, double /* time */)
{
  return true;
}

// Utility function that may be used by derived classes.  Determines
// whether all elements in the model have the same face topology.
// This can be used to speed-up certain algorithms since they don't
// have to check each face (or group of faces) individually.
void Ioss::DatabaseIO::set_common_side_topology() const
{
  Ioss::DatabaseIO *new_this = const_cast<Ioss::DatabaseIO*>(this);

  Ioss::ElementBlockContainer element_blocks =
    get_region()->get_element_blocks();
  Ioss::ElementBlockContainer::const_iterator I  = element_blocks.begin();
  Ioss::ElementBlockContainer::const_iterator IE = element_blocks.end();

  while (I != IE) {
    size_t element_count = (*I)->get_property("entity_count").get_int();

    // Check face types.
    if (element_count > 0) {
      if (commonSideTopology != NULL || I == element_blocks.begin()) {
	ElementTopology* side_type = (*I)->topology()->boundary_type();
	if (commonSideTopology == NULL) // First block
	  new_this->commonSideTopology = side_type;
	if (commonSideTopology != side_type) { // Face topologies differ in mesh
	  new_this->commonSideTopology = NULL;
	  return;
	}
      }
    }
    ++I;
  }
}

void Ioss::DatabaseIO::add_information_records(const std::vector<std::string> &info)
{
  informationRecords.reserve(informationRecords.size()+info.size());
  informationRecords.insert(informationRecords.end(), info.begin(), info.end());
}

void Ioss::DatabaseIO::add_information_record(const std::string &info)
{
  informationRecords.push_back(info);
}

void Ioss::DatabaseIO::set_block_omissions(const std::vector<std::string> &omissions)
{
  blockOmissions.assign(omissions.begin(), omissions.end());
  std::sort(blockOmissions.begin(), blockOmissions.end());
}

// Check topology of all sides (face/edges) in model...
void Ioss::DatabaseIO::check_side_topology() const
{
  // The following code creates the sideTopology sets which contain
  // a list of the side topologies in this model.
  //
  // If sideTopology.size() > 1 --> the model has sides with mixed
  // topology (i.e., quads and tris).
  //
  // If sideTopology.size() == 1 --> the model has homogeneous sides
  // and each side is of the topology type 'sideTopology[0]'
  //
  // This is used in other code speed up some tests.

  // Spheres and Circle have no faces/edges, so handle them special...
  bool all_sphere = true;

  if (sideTopology.empty()) {
    // Set contains (parent_element, boundary_topology) pairs...
    std::set<std::pair<const ElementTopology*, const ElementTopology*> > side_topo;

    Ioss::ElementBlockContainer element_blocks =
      get_region()->get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator I;

    for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
      const Ioss::ElementBlock *block = *I;
      const ElementTopology *elem_type = block->topology();
      const ElementTopology *side_type = elem_type->boundary_type();
      if (side_type == NULL) {
	// heterogeneous sides.  Iterate through...
	int size = elem_type->number_boundaries();
	for (int i=1; i <= size; i++) {
	  side_type = elem_type->boundary_type(i);
	  side_topo.insert(std::make_pair(elem_type, side_type));
	  all_sphere = false;
	}
      } else {
	// homogenous sides.
	side_topo.insert(std::make_pair(elem_type, side_type));
	all_sphere = false;
      }
    }
    if (all_sphere) {
      // If we end up here, the model either contains all spheres, or there are no
      // element blocks in the model...
      const ElementTopology *ftopo = ElementTopology::factory("unknown");
      if (element_blocks.empty()) {
	side_topo.insert(std::make_pair(ftopo, ftopo));
      } else {
	const Ioss::ElementBlock *block = *element_blocks.begin();
	side_topo.insert(std::make_pair(block->topology(), ftopo));
      }
    }
    assert(side_topo.size() > 0);
    assert(sideTopology.size() == 0);
    // Copy into the sideTopology container...
    Ioss::DatabaseIO *new_this = const_cast<Ioss::DatabaseIO*>(this);
    std::copy(side_topo.begin(), side_topo.end(),
	      std::back_inserter(new_this->sideTopology));
  }
  assert(sideTopology.size() > 0);
}

#include <sys/time.h>

namespace {
  static struct timeval tp;
  static double initial_time = -1.0;

  void log_field(const char *symbol, const Ioss::GroupingEntity *entity,
		 const Ioss::Field &field, bool single_proc_only,
		 const Ioss::ParallelUtils &util)
  {
    if (initial_time < 0.0) {
      gettimeofday(&tp, NULL);
      initial_time = (double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    }

    std::vector<int64_t> all_sizes;
    if (single_proc_only) {
      all_sizes.push_back(field.get_size());
    } else {
      util.gather(field.get_size(), all_sizes);
    }

    if (util.parallel_rank() == 0 || single_proc_only) {
      std::string name = entity->name();
      std::ostringstream strm;
      gettimeofday(&tp, NULL);
      double time_now = (double)tp.tv_sec+(1.e-6)*tp.tv_usec;
      strm << symbol << " [" << std::fixed << std::setprecision(3)
	   << time_now-initial_time << "]\t";

      int64_t total = 0;
      // Now append each processors size onto the stream...
      std::vector<int64_t>::const_iterator pos = all_sizes.begin();
      for (; pos != all_sizes.end(); ++pos) {
	strm << std::setw(8) << *pos << ":";
	total += *pos;
      }
      if (util.parallel_size() > 1)
	strm << std::setw(8) << total;
      strm << "\t" << name << "/" << field.get_name() << "\n";
      std::cout << strm.str();
    }
  }
}
