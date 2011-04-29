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
#include <Ioss_SubSystem.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Utils.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_Region.h>
#include <Ioss_ElementTopology.h>
#include <string>

#include <set>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <assert.h>

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
			     MPI_Comm communicator)
  : commonSideTopology(NULL), DBFilename(filename), dbState(STATE_INVALID),
    isParallel(false), myProcessor(0), cycleCount(0), overlayCount(0),
    fieldSuffixSeparator('_'), splitType(Ioss::SPLIT_BY_TOPOLOGIES), dbUsage(db_usage),
    surfaceSplitBackwardCompatibility(false), nodeGlobalIdBackwardCompatibility(false),
    util_(communicator), region_(region), isInput(is_input_event(db_usage)),
    singleProcOnly(db_usage == WRITE_HISTORY || db_usage == WRITE_HEARTBEAT || Ioss::SerializeIO::isEnabled()),
    doLogging(false)
{
  isParallel  = util_.parallel_size() > 1;
  myProcessor = util_.parallel_rank();
}

Ioss::DatabaseIO::~DatabaseIO()
{
}

// The get_field and put_field functions are just a wrapper around the
// pure virtual get_field_internal and put_field_internal functions,
// but this lets me add some debug/checking/common code to the
// functions without having to do it in the calling code or in the
// derived classes code.  This also fulfills the hueristic that a
// public interface should not contain pure virtual functions.

int Ioss::DatabaseIO::get_field(const Ioss::Region* reg, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, reg, field, util_));
  if (get_logging()) {
    log_field(">", reg, field, singleProcOnly, util_);
  }
  return get_field_internal(reg, field, data, data_size);
}

int Ioss::DatabaseIO::get_field(const Ioss::NodeBlock* nb, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, nb, field, util_));
  if (get_logging()) {
    log_field(">", nb, field, singleProcOnly, util_);
  }
  return get_field_internal(nb, field, data, data_size);
}

int Ioss::DatabaseIO::get_field(const Ioss::ElementBlock* eb, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, eb, field, util_));
  if (get_logging()) {
    log_field(">", eb, field, singleProcOnly, util_);
  }
  return get_field_internal(eb, field, data, data_size);
}

int Ioss::DatabaseIO::get_field(const Ioss::SideBlock* fb, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, fb, field, util_));
  if (get_logging()) {
    log_field(">", fb, field, singleProcOnly, util_);
  }
  return get_field_internal(fb,  field, data, data_size);
}

int Ioss::DatabaseIO::get_field(const Ioss::NodeSet* ns, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, ns, field, util_));
  if (get_logging()) {
    log_field(">", ns, field, singleProcOnly, util_);
  }
  return get_field_internal(ns, field, data, data_size);
}

int Ioss::DatabaseIO::get_field(const Ioss::SideSet* fs, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, fs, field, util_));
  if (get_logging()) {
    log_field(">", fs, field, singleProcOnly, util_);
  }
  return get_field_internal(fs, field, data, data_size);
}

int Ioss::DatabaseIO::get_field(const Ioss::CommSet* cs, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, cs, field, util_));
  if (get_logging()) {
    log_field(">", cs, field, singleProcOnly, util_);
  }
  return get_field_internal(cs, field, data, data_size);
}


int Ioss::DatabaseIO::put_field(const Ioss::Region* reg, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, reg, field, util_));
  if (get_logging()) {
    log_field("<", reg, field, singleProcOnly, util_);
  }
  return put_field_internal(reg, field, data, data_size);
}

int Ioss::DatabaseIO::put_field(const Ioss::NodeBlock* nb, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, nb, field, util_));
  if (get_logging()) {
    log_field("<", nb, field, singleProcOnly, util_);
  }
  return put_field_internal(nb, field, data, data_size);
}

int Ioss::DatabaseIO::put_field(const Ioss::ElementBlock* eb, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, eb, field, util_));
  if (get_logging()) {
    log_field("<", eb, field, singleProcOnly, util_);
  }
  return put_field_internal(eb, field, data, data_size);
}

int Ioss::DatabaseIO::put_field(const Ioss::SideBlock* fb, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, fb, field, util_));
  if (get_logging()) {
    log_field("<", fb, field, singleProcOnly, util_);
  }
  return put_field_internal(fb, field, data, data_size);
}

int Ioss::DatabaseIO::put_field(const Ioss::NodeSet* ns, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, ns, field, util_));
  if (get_logging()) {
    log_field("<", ns, field, singleProcOnly, util_);
  }
  return put_field_internal(ns, field, data, data_size);
}

int Ioss::DatabaseIO::put_field(const Ioss::SideSet* fs, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, fs, field, util_));
  if (get_logging()) {
    log_field("<", fs, field, singleProcOnly, util_);
  }
  return put_field_internal(fs, field, data, data_size);
}

int Ioss::DatabaseIO::put_field(const Ioss::CommSet* cs, const Ioss::Field& field,
				void *data, size_t data_size) const
{
  assert(is_parallel_consistent(singleProcOnly, cs, field, util_));
  if (get_logging()) {
    log_field("<", cs, field, singleProcOnly, util_);
  }
  return put_field_internal(cs, field, data, data_size);
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
    int element_count = (*I)->get_property("entity_count").get_int();

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
#include <unistd.h>
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

    std::vector<int> all_sizes;
    if (single_proc_only) {
      all_sizes.push_back(static_cast<int>(field.get_size()));
    } else {
      util.gather(static_cast<int>(field.get_size()), all_sizes);
    }

    if (util.parallel_rank() == 0 || single_proc_only) {
      std::string name = entity->name();
      std::ostringstream strm;
      gettimeofday(&tp, NULL);
      double time_now = (double)tp.tv_sec+(1.e-6)*tp.tv_usec;
      strm << symbol << " [" << std::fixed << std::setprecision(3)
	   << time_now-initial_time << "]\t";

      int total = 0;
      // Now append each processors size onto the stream...
      std::vector<int>::const_iterator pos = all_sizes.begin();
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
