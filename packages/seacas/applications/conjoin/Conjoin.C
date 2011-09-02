// Copyright(C) 2009-2010 Sandia Corporation.
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
#include <iostream>
#include <iomanip>
#include <iterator>
#include <set>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>
#include <string>
#include <exception>
#include <stdexcept>

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cfloat>
#include <cmath>
#include <limits>

#include <cstring>
#include <ctime>
#include <sys/times.h>
#include <ctype.h>
#include "smart_assert.h"
#include "adler.h"
#include <sys/utsname.h>

namespace {
template <typename T>
bool approx_equal(T v1, T v2)
{
#if 0
  static const T tolerance = 100.0 * std::numeric_limits<T>::epsilon();
  return std::fabs(v1 - v2) <= std::fabs(v1+v2)*tolerance;
#else
  return (float)v1 == (float)v2;
#endif
}
}

struct NodeInfo {
  NodeInfo() :
    id(0), x(0.0), y(0.0), z(0.0) {}
  NodeInfo(int id_, double x_, double y_, double z_) :
    id(id_), x(x_), y(y_), z(z_) {}
  int id;
  double x, y, z;

  bool operator==(const NodeInfo &other) const
  {
#if 0
    return (!(*this < other) && !(other < *this));
#else
    return id == other.id &&
      approx_equal(x, other.x) &&
      approx_equal(y, other.y) &&
      approx_equal(z, other.z);
#endif
  }

  bool operator!=(const NodeInfo &other) const
  {
    return !(*this == other);
  }

  bool operator<(const NodeInfo &other) const
  {
    if (id < other.id) return true;
    if (id > other.id) return false;
    SMART_ASSERT(id == other.id);
    if (!approx_equal(x, other.x) && x < other.x) return true;
    if (!approx_equal(x, other.x) && x > other.x) return false;
    SMART_ASSERT(approx_equal(x, other.x));
    if (!approx_equal(y, other.y) && y < other.y) return true;
    if (!approx_equal(y, other.y) && y > other.y) return false;
    SMART_ASSERT(approx_equal(y, other.y));
    if (!approx_equal(z, other.z) && z < other.z) return true;
    return false;
  }
};

typedef std::vector<NodeInfo>   GlobalMap;
typedef GlobalMap::iterator     GMapIter;

typedef std::vector<std::pair<int, size_t> > GlobalElemMap;
typedef GlobalElemMap::iterator GElemMapIter;

#include <Internals.h>
#include <ExodusFile.h>
#include <ExodusEntity.h>
#include <SystemInterface.h>
#include <Version.h>
#include <Variables.h>
#include <ObjectType.h>

#include <to_string.h>
#include <exodusII.h>

#if EX_API_VERS_NODOT <= 467
#error "Requires exodusII version 4.68 or later"
#endif

extern void add_to_log(const char *name, int elapsed);
extern double conjoin_timer();

namespace {
  template <class  T>
    struct TimeStepMap {
      TimeStepMap(size_t part, int step, T time)
	: partNumber(part), localStepNumber(step), timeValue(time) {}
      size_t partNumber;
      int localStepNumber;
      T timeValue;
    };

  std::string time_stamp(const std::string &format);
  std::string format_time(double seconds);
  int get_width(int max_value);

  ex_entity_type exodus_object_type(Excn::ObjectType &conjoin_type) {
    switch (conjoin_type) {
    case Excn::EBLK:
      return EX_ELEM_BLOCK;
    case Excn::SSET:
      return EX_SIDE_SET;
    case Excn::NSET:
      return EX_NODE_SET;
    default:
      SMART_ASSERT(1==0 && "Invalid Object Type in exodus_object_type")(conjoin_type);
      return EX_INVALID;
    }
  }

  char **get_name_array(int size, size_t length) {
    char **names = NULL;
    if (size > 0) {
      names = new char* [size];
      for (int i=0; i < size; i++) {
	names[i] = new char [length+1];
	std::memset(names[i], '\0', length+1);
      }
    }
    return names;
  }

  void free_name_array(char **names, int size) {
    for (int i=0; i < size; i++) {
      delete [] names[i];
    }
    delete [] names;
    names = NULL;
  }
}

std::string tsFormat = "[%H:%M:%S] ";

// prototypes

// The main program templated to permit float/double transfer.
template <typename T>
int conjoin(Excn::SystemInterface &interface, T /* dummy */);

namespace {
void compress_white_space(char *str);
void add_info_record(char *info_record, int size);
void put_mesh_summary(const Excn::Mesh& mesh);

template <typename T>
void get_put_qa(int exo_id, int out_id, const std::vector<TimeStepMap<T> > &global_times,
		Excn::SystemInterface &interface);
int  get_put_coordinate_names(int in, int out, int dimensionality);

template <typename T>
int get_put_coordinates(Excn::Mesh& global, size_t part_count,
			std::vector<Excn::Mesh> &local_mesh,
			T dummy);

template <typename T>
int get_coordinates(int id, int dimensionality, size_t num_nodes,
		    std::vector<int> local_node_to_global,
		    size_t part, std::vector<T> &x, std::vector<T> &y,
		    std::vector<T> &z);

StringVector get_exodus_variable_names(int id, ex_entity_type elType,
				       int count);

template <typename T>
void filter_truth_table(int id,
			Excn::Mesh& global,
			std::vector<T> &glob_blocks,
			Excn::Variables &vars,
			const StringIdVector &variable_names);

template <typename T>
void get_truth_table(Excn::Mesh& global, std::vector<T> &glob_blocks,
		     std::vector<Excn::Mesh> &local, Excn::Variables &vars, int debug);

template <typename U>
void create_output_truth_table(const Excn::Mesh& global,
			       std::vector<U> &global_sets,
			       Excn::Variables &vars,
			       std::vector<int> &truth_table);

template <typename T, typename U>
int read_write_master_values(Excn::Variables &vars,
			     const Excn::Mesh& global,
			     std::vector<U> &global_sets,
			     std::vector<Excn::Mesh> &local_mesh,
			     std::vector<std::vector<U> > &local_sets,
			     std::vector<T> &values, 
			     size_t part, int step, int step_out);

void get_variable_params(int id, Excn::Variables &vars,
			 const StringIdVector &variable_list);

void get_put_variable_names(int id, int idout, Excn::Variables &vars, Excn::SystemInterface &interface,
			      int *combined_status_variable_index = NULL);

void build_reverse_node_map(std::vector<Excn::Mesh> &local_mesh,
			    Excn::Mesh *global, size_t part_count,
			    GlobalMap &global_node_map);

void build_reverse_element_map(std::vector<Excn::Mesh> &local_mesh,
			       std::vector<std::vector<Excn::Block> > &blocks,
			       std::vector<Excn::Block> &glob_blocks,
			       Excn::Mesh *global, size_t part_count,
			       GlobalElemMap &global_element_map);

void get_nodesets(int total_node_count,
		  std::vector<Excn::Mesh> &local_mesh,
		  std::vector<std::vector<Excn::NodeSet> > &nodesets,
		  std::vector<Excn::NodeSet> &glob_sets);

void get_element_blocks(const std::vector<Excn::Mesh> &local_mesh,
			const Excn::Mesh &global,
			std::vector<std::vector<Excn::Block> > &blocks,
			std::vector<Excn::Block> &glob_blocks);
template<typename T>
void put_element_blocks(std::vector<Excn::Mesh> &local_mesh,
			std::vector<std::vector<Excn::Block> > &blocks,
			std::vector<Excn::Block> &glob_blocks,
			T single_or_double);

void put_nodesets(std::vector<Excn::NodeSet> &glob_sets);

void get_sideset_metadata(std::vector<Excn::Mesh> &local_mesh,
			  std::vector< std::vector<Excn::SideSet> > &sets,
			  std::vector<Excn::SideSet> &glob_ssets);
void get_put_sidesets(std::vector<Excn::Mesh> &local_mesh,
		      std::vector< std::vector<Excn::SideSet> > &sets,
		      std::vector<Excn::SideSet> &glob_ssets);

template <typename T>
void add_status_variable(int id_out, const Excn::Mesh& global,
			 const std::vector<Excn::Block> &blocks,
			 const std::vector<Excn::Block> &glob_blocks,
			 const Excn::IntVector &local_element_to_global,
			 int step, int variable, T alive_value,
			 int combined_variable_index);

size_t find_max_entity_count(size_t part_count,
			     std::vector<Excn::Mesh> &local_mesh,
			     const Excn::Mesh& global,
			     std::vector<std::vector<Excn::Block> > &blocks,
			     std::vector<std::vector<Excn::NodeSet> > &nodesets,
			     std::vector<std::vector<Excn::SideSet> > &sidesets);
  
int case_compare(const std::string &s1, const std::string &s2);

template <typename T>
void verify_set_position_mapping(const std::string &type,
				 size_t part_count,
				 const std::vector<T> &global_sets,
				 const std::vector<std::vector<T> > &sets)
{
  bool problem = false;
  for (size_t i=0; i < global_sets.size(); i++) {
    int loc_pos = global_sets[i].position_;
    for (size_t p=0; p < part_count; p++) {
      if (global_sets[i].id  != sets[p][loc_pos].id ||
	  global_sets[i].name_ != sets[p][loc_pos].name_ ||
	  sets[p][loc_pos].position_ != (int)i) {
	problem = true;
	std::cerr << "\nMismatch for global " << type << " at position " << i
		  << " and local " << type << " in part " << p+1 << " position " << loc_pos << "\n";
	global_sets[i].dump();
	sets[p][loc_pos].dump();
      }
    }
  }

  if (problem) {
    throw std::runtime_error(type + " mismatch");
  }
}

template <typename T>
void uniqify(std::vector<T> &map) {
  std::sort(map.begin(), map.end());
  map.erase(std::unique(map.begin(), map.end()), map.end());
  // shrink-to-fit...
  std::vector<T>(map).swap(map);
}
}

unsigned int debug_level = 0;
const float FILL_VALUE = FLT_MAX;

using namespace Excn;

int main(int argc, char* argv[])
{
#if defined(__LIBCATAMOUNT__)
  setlinebuf(stderr);
#endif
  try {
  time_t begin_time = time(NULL);
  SystemInterface::show_version();

  SystemInterface interface;
  bool ok = interface.parse_options(argc, argv);

  if (!ok) {
    std::cerr << "\nERROR: Problems parsing command line arguments.\n\n";
    exit(EXIT_FAILURE);
  }

  debug_level = interface.debug();

  if (debug_level & 64)
    ex_opts(EX_VERBOSE|EX_DEBUG);
  else
    ex_opts(0);

  int error = 0;

  if (!ExodusFile::initialize(interface)) {
    std::cerr << "ERROR: Problem initializing input and/or output files.\n";
    exit(EXIT_FAILURE);
  }

  if (ExodusFile::io_word_size() == 4) {
    error = conjoin(interface, (float)0.0);
  } else {
    error = conjoin(interface, (double)0.0);
  }

  ExodusFile::close_all();

  time_t end_time = time(NULL);
  add_to_log(argv[0], (int)(end_time-begin_time));
  return (error);
  }
  catch (std::exception &e) {
    std::cerr << "ERROR: Standard exception: " << e.what() << std::endl;
  }
}

template <typename T>
int conjoin(SystemInterface &interface, T /* dummy */)
{
  SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());

  const T alive = interface.alive_value();
  size_t part_count = interface.inputFiles_.size();

  char* mytitle = new char[MAX_LINE_LENGTH + 1];
  memset(mytitle, '\0', MAX_LINE_LENGTH+1);

  Mesh global;

  std::vector<Mesh> local_mesh(part_count);
  
  // ******************************************************************
  // 1. Read global info

  int error = 0;

  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);

  std::string title0;

  for (size_t p = 0; p < part_count; p++) {
    error += ex_get_init(ExodusFile(p), mytitle,
			 &local_mesh[p].dimensionality,
                         &local_mesh[p].nodeCount,
			 &local_mesh[p].elementCount,
			 &local_mesh[p].blockCount,
			 &local_mesh[p].nodesetCount,
			 &local_mesh[p].sidesetCount);

    local_mesh[p].title = mytitle;

    if (p == 0) {
      global.title = mytitle;
      global.dimensionality = local_mesh[p].count(DIM);
      global.blockCount     = local_mesh[p].count(EBLK);
    } else {
      SMART_ASSERT(global.count(DIM)  == local_mesh[p].count(DIM))  (global.count(DIM)) (local_mesh[p].count(DIM))(p);
      SMART_ASSERT(global.count(EBLK) == local_mesh[p].count(EBLK)) (global.count(EBLK))(local_mesh[p].count(EBLK))(p);
    }

    local_mesh[p].localNodeToGlobal.resize(local_mesh[p].count(NODE));
    local_mesh[p].localElementToGlobal.resize(local_mesh[p].count(ELEM));

  } // end for (p=0..part_count)

  delete [] mytitle;

  if (interface.omit_nodesets()) {
    global.nodesetCount = 0;
  }

  if (interface.omit_sidesets()) {
    global.sidesetCount = 0;
  }

  // Get database times...  Save mapping from time step to part providing that time step.
  std::vector<TimeStepMap<T> > global_times;
  
  T t_min = FLT_MAX;
  for (size_t p = part_count; p > 0; p--) {
    bool used = false;
    int nts = ex_inquire_int(ExodusFile(p-1), EX_INQ_TIME);
    std::vector<T> times(nts);
    ex_get_all_times(ExodusFile(p-1), &times[0]);

    // A database will include all times from step 0 up to the
    // last time that is less than t_min on the following database.
    // Note that we are iterating through the databases from last to first...
    int i = 0;
    for (i=nts; i > 0; i--) {
      if (times[i-1] <t_min) {
	used = true;
	global_times.push_back(TimeStepMap<T>(p-1, i-1, times[i-1]));
      }
    }
    local_mesh[p-1].timestepCount = i;
    t_min = t_min < times[0] ? t_min : times[0];
    if (!used) {
      std::string part = "Part " + ToString(p) + ": ";
      part += interface.inputFiles_[p-1];
      std::cerr << "\nWARNING: " << part
		<< " does not contain any time steps which will be used in conjoined file.\n";
      local_mesh[p-1].isActive = false;
    }
  }
  global.timestepCount = global_times.size();
  std::reverse(global_times.begin(), global_times.end());

  // Need these throughout run, so declare outside of this block...
  // TODO: Add these to the "Mesh" class.
  std::vector<Block> glob_blocks(global.count(EBLK));
  std::vector<std::vector<Block> > blocks(part_count);

  std::vector<SideSet> glob_ssets;
  std::vector<std::vector<SideSet> > sidesets(part_count);

  std::vector<NodeSet> glob_nsets;
  std::vector<std::vector<NodeSet> > nodesets(part_count);
  {
    // Now, build the reverse global node map which permits access of the
    // local id given the global id.
    GlobalMap global_node_map;
    build_reverse_node_map(local_mesh, &global, part_count, global_node_map);

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);

    // ****************************************************************************
    // 5. Get Block information including element attributes
    // must check for zero length blocks
    get_element_blocks(local_mesh, global, blocks, glob_blocks);

    GlobalElemMap global_element_map(global.count(ELEM));
    build_reverse_element_map(local_mesh, blocks,
			      glob_blocks, &global, part_count, global_element_map);

    //
    //    NOTE:  Node set/side set information can be different for each part
    /************************************************************************/
    // 7. Get Side sets
    if (!interface.omit_sidesets()) {
      if (debug_level & 1)
	std::cerr << time_stamp(tsFormat);
      get_sideset_metadata(local_mesh, sidesets, glob_ssets);
      if (global.count(SSET) != glob_ssets.size()) {
	global.sidesetCount = glob_ssets.size();
      }
    }


    /************************************************************************/
    // 6. Get Node sets
    if (!interface.omit_nodesets()) {
      if (debug_level & 1)
	std::cerr << time_stamp(tsFormat);
      get_nodesets(global.count(NODE), local_mesh, nodesets, glob_nsets);
      if (global.count(NSET) != glob_nsets.size()) {
	global.nodesetCount = glob_nsets.size();
      }
    }

    /************************************************************************/
    // Start writing the output file...

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    CommunicationMetaData comm_data;

    // Create the output file...
    ExodusFile::create_output(interface);

    put_mesh_summary(global);
    
    get_put_qa(ExodusFile(0), ExodusFile::output(), global_times, interface);

    Internals exodus(ExodusFile::output(), ExodusFile::max_name_length());

    exodus.write_meta_data(global, glob_blocks, glob_nsets, glob_ssets, comm_data);

    // Output bulk mesh data....
    put_nodesets(glob_nsets);

    // c.2.  Write Global Node Number Map
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    {
      std::vector<int> global_map(global.count(NODE));
      for (size_t i=0; i < global.count(NODE); i++) {
	global_map[i] = global_node_map[i].id;
      }
      error=ex_put_node_num_map(ExodusFile::output(),&global_map[0]);
    }

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    if (global_element_map.size() > 0) {
      std::vector<int> global_map(global.count(ELEM));
      for (size_t i=0; i < global.count(ELEM); i++) {
	global_map[i] = global_element_map[i].first;
      }
      ex_put_elem_num_map(ExodusFile::output(), &global_map[0]);
    }

    T dummy = 0.0;
    put_element_blocks(local_mesh, blocks, glob_blocks, dummy);
    get_put_sidesets(local_mesh, sidesets, glob_ssets);
  }				
  // ************************************************************************
  // 2. Get Coordinate Info.
  {
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);

    error += get_put_coordinates(global, part_count, local_mesh, (T)0.0);

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    std::cerr << "Wrote coordinate information...\n";
  }
  // ####################TRANSIENT DATA SECTION###########################
  // ***********************************************************************
  // 9. Get Variable Information and names

  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);

  //  I. read number of variables for each type. Note that exodusII does not
  //     provide for history variables
  //  NOTE: it is assumed that every part has the same global, nodal,
  //        and element lists

  bool add_n_status = interface.nodal_status_variable() != "NONE";
  bool add_e_status = interface.element_status_variable() != "NONE";
  Variables global_vars(GLOBAL);
  Variables nodal_vars(NODE, add_n_status);
  Variables element_vars(EBLK, add_e_status);
  Variables nodeset_vars(NSET);
  Variables sideset_vars(SSET);
  
  {
    ExodusFile id(0);

    get_variable_params(id, global_vars,  interface.global_var_names());
    get_variable_params(id, nodal_vars,   interface.node_var_names());
    get_variable_params(id, element_vars, interface.elem_var_names());
    get_variable_params(id, nodeset_vars, interface.nset_var_names());
    get_variable_params(id, sideset_vars, interface.sset_var_names());
    
    get_truth_table(global, glob_blocks, local_mesh, element_vars,  4);
    filter_truth_table(id, global, glob_blocks, element_vars, interface.elem_var_names());

    get_truth_table(global, glob_nsets, local_mesh, nodeset_vars, 32);
    filter_truth_table(id, global, glob_nsets,  nodeset_vars, interface.nset_var_names());

    get_truth_table(global, glob_ssets, local_mesh, sideset_vars, 16);
    filter_truth_table(id, global, glob_ssets,  sideset_vars, interface.sset_var_names());
  }

  // There is a slightly tricky situation here. The truthTable block order
  // is based on the ordering of the blocks on the input databases.
  // These blocks may have been reordered on output to make the 'offset'
  // variables line up correctly when the element ids are mapped back to
  // the "overall global" order.  There is not much problem since most
  // calls  outputting block-related items pass the id and don't rely
  // on ordering. However, the truth table is one of the exceptions
  // and we need to reorder the truth table to match the output block
  // order. After this call, we can use the original ordering, so just
  // need a temporary vector here...
  if (global_vars.count(OUT) + nodal_vars.count(OUT) +
      element_vars.count(OUT) + nodeset_vars.count(OUT) +
      sideset_vars.count(OUT) > 0) {

    std::vector<int> elem_truth_table(global.truthTable[EBLK].size());
    create_output_truth_table(global, glob_blocks, element_vars, elem_truth_table);

    ex_put_all_var_param(ExodusFile::output(),
			 global_vars.count(OUT),
			 nodal_vars.count(OUT),
			 element_vars.count(OUT), &elem_truth_table[0],
			 nodeset_vars.count(OUT), &global.truthTable[NSET][0],
			 sideset_vars.count(OUT), &global.truthTable[SSET][0]);
  }

  // II. read/write the variable names
  int combined_status_variable_index = 0;
  {
    ExodusFile id(0);
    get_put_variable_names(id, ExodusFile::output(), global_vars,  interface);
    get_put_variable_names(id, ExodusFile::output(), nodal_vars,   interface);
    get_put_variable_names(id, ExodusFile::output(), element_vars, interface,
			   &combined_status_variable_index);
    get_put_variable_names(id, ExodusFile::output(), nodeset_vars, interface);
    get_put_variable_names(id, ExodusFile::output(), sideset_vars, interface);
  }
  ex_update(ExodusFile::output());

  /**********************************************************************/
  // 10. Get Transient Data
  //     This routine reads in a time dump from an EXODUSII file

  size_t num_time_steps = global.count(TIME);

  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);

  std::vector<T> global_values(global_vars.count(IN));
  std::vector<T> output_global_values(global_vars.count(OUT));

  // Determine maximum number of entities on any part
  int max_ent = find_max_entity_count(part_count, local_mesh, global, blocks, nodesets, sidesets);
  std::vector<T> values(max_ent);
  
  // Stage II.  Extracting transient variable data.
  //            loop over time steps

  // Determine if user wants a subset of timesteps transferred to the output file.
  // Time steps for output file
  double start_time = conjoin_timer();

  // Used for output formatting at end of loop.
  int percent_width = 0;
  int field_width   = 3;
  if (num_time_steps > 100) {
    percent_width = 1;
    field_width   = 5;
  }

  int element_width = get_width(global.count(ELEM));
  int step_width = get_width(num_time_steps);
  int part_width = get_width(part_count+1);
  int loc_step_width = 0;
  for (size_t time_step = 0; time_step < global_times.size(); time_step++) {
    if (loc_step_width < global_times[time_step].localStepNumber+1)
      loc_step_width = global_times[time_step].localStepNumber+1;
  }
  loc_step_width = get_width(loc_step_width);
  
  size_t time_step_out = 0;
  for (size_t time_step = 0; time_step < num_time_steps; time_step++) {
    time_step_out++;

    T time_val = 0;
    size_t p = global_times[time_step].partNumber;
    ExodusFile id(p);

    // read in and write out the time step information
    error += ex_get_time(id, global_times[time_step].localStepNumber+1,   &time_val);
    SMART_ASSERT(time_val == global_times[time_step].timeValue);
    error += ex_put_time(ExodusFile::output(), time_step_out, &time_val);

    if (global_vars.count(OUT) > 0) {
      if (debug_level & 1)
	std::cerr << time_stamp(tsFormat) << "Global Variables...\n";
      error += ex_get_var(id, global_times[time_step].localStepNumber+1, EX_GLOBAL, 0, 0,
			  global_vars.count(), (void*)&global_values[0]);
      // Map ...
      for (int ig=0; ig < global_vars.count(IN); ig++) {
	if (global_vars.index_[ig] > 0) {
	  SMART_ASSERT(ig < (int)global_values.size());
	  output_global_values[global_vars.index_[ig]-1] = global_values[ig];
	}
      }
      error += ex_put_var(ExodusFile::output(),time_step_out, EX_GLOBAL, 1, 0,
			  global_vars.count(OUT), &output_global_values[0]);
      
    }

    // ========================================================================
    // Nodal Values...
    if (debug_level & 1) {
      std::cerr << time_stamp(tsFormat) << "Nodal Variables...\n";
    }

    if (nodal_vars.count(OUT) > 0) {
      int node_count = local_mesh[p].count(NODE);
      std::vector<T> master_nodal_values(global.count(NODE));

      int offset = nodal_vars.addStatus ? 1 : 0;
      for (int i = 0; i < nodal_vars.count(OUT)-offset; i++) { // Last output variable may be status
	for (int j = 0; j < nodal_vars.count(IN); j++) {
	  if (nodal_vars.index_[j]-1 == i) {
	    std::fill(master_nodal_values.begin(), master_nodal_values.end(), 0.0);

	    error += ex_get_var(id, global_times[time_step].localStepNumber+1, EX_NODAL, j+1, 0,
				node_count, &values[0]);
	  
	    // copy values to master nodal value information
	    for (int jj = 0; jj < node_count; jj++) {
	      // Map local nodal value to global location...
	      int nodal_value = local_mesh[p].localNodeToGlobal[jj];
	      master_nodal_values[nodal_value] = values[jj];
	    }
	    error+=ex_put_var(ExodusFile::output(), time_step_out, EX_NODAL, i+1, 0,
			      global.count(NODE), &master_nodal_values[0]);
	    break;
	  }
	}

	// Fill "node_status" variable -- 'alive' for alive; 1-alive for dead.
	// It is the last output variable...
	if (nodal_vars.addStatus) {
	  SMART_ASSERT(alive == 0.0 || alive == 1.0)(alive);
	  std::fill(master_nodal_values.begin(), master_nodal_values.end(), (1.0 - alive));
	  for (int j = 0; j < node_count; j++) {
	    // Map local nodal value to global location...
	    int nodal_value = local_mesh[p].localNodeToGlobal[j];
	    master_nodal_values[nodal_value] = alive;
	  }
	  
	  error+=ex_put_var(ExodusFile::output(), time_step_out, EX_NODAL, nodal_vars.count(OUT), 0,
			    global.count(NODE), &master_nodal_values[0]);
	}
      }
    }

    // ========================================================================
    // Extracting element transient variable data
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat) << "Element Variables...\n";
    
    if (element_vars.count(IN) > 0) {
      read_write_master_values(element_vars, global, glob_blocks, local_mesh, blocks,
			       values, p, global_times[time_step].localStepNumber, time_step_out);
    }
    
    // Add element status variable...
    // Use the output time step for writing data
    if (interface.element_status_variable() != "NONE") {
      add_status_variable(ExodusFile::output(), global, blocks[p], glob_blocks,
			  local_mesh[p].localElementToGlobal,
			  time_step_out,
			  element_vars.index_[element_vars.count(IN)], alive,
			  combined_status_variable_index);
    }

    // ========================================================================
    // Extracting sideset transient variable data
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat) << "Sideset Variables...\n";

    if (sideset_vars.count(IN) > 0) {
      read_write_master_values(sideset_vars, global, glob_ssets, local_mesh, sidesets,
			       values, p, global_times[time_step].localStepNumber, time_step_out);
    }

    // ========================================================================
    // Extracting nodeset transient variable data
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat) << "Nodeset Variables...\n";
    
    if (nodeset_vars.count(IN) > 0) {
      read_write_master_values(nodeset_vars, global, glob_nsets, local_mesh, nodesets,
			       values, p, global_times[time_step].localStepNumber, time_step_out);
    }
    
    // ========================================================================
    ex_update(ExodusFile::output());

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);

    std::cerr << "Step " << std::setw(step_width) << time_step+1 <<  "/" << num_time_steps << ", time "
	      << std::scientific << std::setprecision(4) << time_val
	      << " (Part " << std::setw(part_width) << p+1 << "/" << part_count 
	      << ", step " << std::setw(loc_step_width) 
	      << global_times[time_step].localStepNumber+1 << ")"
	      << "  Active Elem: " << std::setw(element_width) << local_mesh[p].count(ELEM);

    double cur_time = conjoin_timer();
    double elapsed = cur_time - start_time;
    double time_per_step = elapsed / time_step_out;
    double percentage_done = (time_step_out * 100.0) / global.count(TIME);
    double estimated_remaining = time_per_step * (global.count(TIME) - time_step_out);

    std::cerr << "  ["
	      << std::fixed << std::setw(field_width) << std::setprecision(percent_width)
	      << percentage_done << "%, Elapsed=" << format_time(elapsed)
	      << ", ETA=" << format_time(estimated_remaining)
	      << "]\n";
    if (debug_level & 1) {
      std::cerr << "\n";
    }
  }

  /*************************************************************************/
  // EXIT program
  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);
  std::cerr << "******* END *******\n";
  return(error);
}

namespace {
  template <typename T>
  void get_put_qa(int id, int id_out,
		  const std::vector<TimeStepMap<T> > &global_times,
		  Excn::SystemInterface &interface)
  {
    // NOTE: Assuming info and QA records for all parts
    int error = 0;

    int info_string_len = MAX_LINE_LENGTH;

    size_t num_info_records = ex_inquire_int(id,EX_INQ_INFO);
    //    size_t extra_info = global_times.size() + 2 + 1;
    size_t extra_info = 2 + 1;
      
    char **info_records = get_name_array(num_info_records+extra_info, info_string_len);

    if (num_info_records > 0) {
      error += ex_get_info(id, info_records);
    }

    // Add an info record for CONJOIN
    add_info_record(info_records[num_info_records], MAX_LINE_LENGTH);

#if 0
    // Add time/part mapping...
    for (size_t i=0; i < global_times.size(); i++) {
      std::ostringstream os;
      os << "Step " << std::setw(2) << i+1 << ", time "
	 << std::scientific << std::setprecision(4) << global_times[i].timeValue
	 << " (Part " << global_times[i].partNumber+1
	 << ", step "  << global_times[i].localStepNumber+1 << ")  File: "
	 << interface.inputFiles_[global_times[i].partNumber];
      strncpy(info_records[num_info_records+1+i], os.str().c_str(), MAX_LINE_LENGTH);
      info_records[num_info_records+1+i][MAX_LINE_LENGTH] = '\0';
    }
#endif
    
    error += ex_put_info(id_out,num_info_records+extra_info,info_records);

    free_name_array(info_records, num_info_records+extra_info);

    // II. Get and store QA records, if they exist
    struct qa_element {
      char *qa_record[1][4];
    };
    
    int num_qa_records = ex_inquire_int(id, EX_INQ_QA);
    qa_element *qaRecord = new qa_element[num_qa_records+1];
    for (int i=0; i < num_qa_records+1; i++) {
      for (int j=0; j < 4; j++) {
	qaRecord[i].qa_record[0][j] = new char[MAX_STR_LENGTH+1];
      }
    }
    if (num_qa_records) error += ex_get_qa(id, qaRecord[0].qa_record);

    char buffer[MAX_STR_LENGTH+1];

    strncpy(qaRecord[num_qa_records].qa_record[0][0], qainfo[0], MAX_STR_LENGTH); // Code
    strncpy(qaRecord[num_qa_records].qa_record[0][1], qainfo[2], MAX_STR_LENGTH); // Version

    time_t date_time = time(NULL);
    strftime( buffer, MAX_STR_LENGTH, "%Y/%m/%d", localtime(&date_time) );

    strncpy(qaRecord[num_qa_records].qa_record[0][2], buffer, MAX_STR_LENGTH);

    strftime( buffer, MAX_STR_LENGTH, "%H:%M:%S", localtime(&date_time) );
    strncpy(qaRecord[num_qa_records].qa_record[0][3], buffer, MAX_STR_LENGTH);

    error += ex_put_qa(id_out, num_qa_records+1, qaRecord[0].qa_record);

    for (int i=0; i < num_qa_records+1; i++) {
      for (int j=0; j < 4; j++) {
	delete [] qaRecord[i].qa_record[0][j];
      }
    }
    delete [] qaRecord;
  }

  template <typename T>
  int get_put_coordinates(Mesh& global, size_t part_count,
			  std::vector<Mesh> &local_mesh,
			  T /* dummy */)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    std::vector<T> x(global.count(NODE));
    std::vector<T> y(global.count(NODE));
    std::vector<T> z(global.count(NODE));

    if (debug_level & 8) {
      std::fill(x.begin(), x.end(), FILL_VALUE);
      std::fill(y.begin(), y.end(), FILL_VALUE);
      std::fill(z.begin(), z.end(), FILL_VALUE);
    }

    int error = 0;
    for (size_t p = 0; p < part_count; p++) {
      error += get_coordinates(ExodusFile(p), global.count(DIM),
			       local_mesh[p].count(NODE),
			       local_mesh[p].localNodeToGlobal, p, x, y, z);
    }                            // end for p=0..part_count

    // Get Coordinate Names
    // NOTE: Assuming coordinate names should be
    // the same for all files/parts therefore, only one
    // file/part needs to be loaded
    error += get_put_coordinate_names(ExodusFile(0),
				      ExodusFile::output(),
				      global.count(DIM));
    // Write out coordinate information
    error += ex_put_coord(ExodusFile::output(), (void*)&x[0], (void*)&y[0], (void*)&z[0]);
    return error;
  }

  int get_put_coordinate_names(int in, int out, int dimensionality)
  {
    int error = 0;
    char **coordinate_names = get_name_array(dimensionality, ExodusFile::max_name_length());

    error += ex_get_coord_names (in, coordinate_names);
    error += ex_put_coord_names(out,coordinate_names);
    std::cerr << "Wrote coordinate names..." << std::endl;

    free_name_array(coordinate_names, dimensionality);
    return error;
  }

  template <typename T>
  int get_coordinates(int id, int dimensionality, size_t num_nodes,
		      std::vector<int> local_node_to_global,
		      size_t part, std::vector<T> &x, std::vector<T> &y, std::vector<T> &z)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    SMART_ASSERT(local_node_to_global.size() == num_nodes);
    int error = 0;
    std::vector<T> local_x(num_nodes);
    std::vector<T> local_y(num_nodes);
    std::vector<T> local_z(num_nodes);

    error += ex_get_coord(id, &local_x[0], &local_y[0], &local_z[0]);

    // Check for 2D or 3D coordinates

    if (dimensionality  == 3) {
      if (debug_level & 8) {
	for (size_t i = 0; i < num_nodes; i++) {
	  int node = local_node_to_global[i];
	  if (x[node] != FILL_VALUE && y[node] != FILL_VALUE && z[node] != FILL_VALUE) {
	    if (!approx_equal(x[node],local_x[i]) ||
		!approx_equal(y[node],local_y[i]) ||
		!approx_equal(z[node],local_z[i])) {
	      std::cerr << "\nWARNING: Node " << node+1
			<< " has different coordinates in at least two parts.\n"
			<< "         this may indicate that this id has been reused in the current part.\n"
			<< "         cur value = "
			<< std::scientific << std::setprecision(6)
			<< std::setw(14) << x[node] << std::setw(14) << y[node] << std::setw(14) << z[node] << "\n"
			<< "         new value = "
			<< std::setw(14) << local_x[i] << std::setw(14) << local_y[i] << std::setw(14) << local_z[i]
			<< " from part " << part << std::endl;

	    }
	  }
	}
      }

      for (size_t i = 0; i < num_nodes; i++) {

	// The following WILL overwrite x[node],y[node],z[node] if node is the
	// same for different parts
	int node = local_node_to_global[i];
	x[node] = local_x[i];
	y[node] = local_y[i];
	z[node] = local_z[i];

      }
    }
    else {
      if (debug_level & 8) {
	for (size_t i = 0; i < num_nodes; i++) {
	  int node = local_node_to_global[i];
	  if (x[node] != FILL_VALUE && y[node] != FILL_VALUE) {
	    if (!approx_equal(x[node],local_x[i]) ||
		!approx_equal(y[node],local_y[i])) {
	      std::cerr << "\nWARNING: Node " << node+1
			<< " has different coordinates in at least two parts.\n"
			<< "         this may indicate that this id has been reused in the current part.\n"
			<< "         cur value = "
			<< std::scientific << std::setprecision(6)
			<< std::setw(14) << x[node] << std::setw(14) << y[node] << "\n"
			<< "         new value = "
			<< std::setw(14) << local_x[i] << std::setw(14) << local_y[i]
			<< " from part " << part << std::endl;
	    }
	  }
	}
      }
      for (size_t i = 0; i < num_nodes; i++) {

	// The following WILL overwrite x[node],y[node] if node is the same for
	// different parts
	int node = local_node_to_global[i];

	x[node] = local_x[i];
	y[node] = local_y[i];
      }
    }
    return error;
  }

  void get_element_blocks(const std::vector<Mesh> &local_mesh,
			  const Mesh &global,
			  std::vector<std::vector<Block> > &blocks,
			  std::vector<Block> &glob_blocks)
  {
    size_t part_count = local_mesh.size();
    for (size_t ip =0; ip < part_count; ip++)
      blocks[ip].resize(local_mesh[ip].count(EBLK));

    IntVector block_id(global.count(EBLK));

    int error = 0;
    for (size_t p = 0; p < part_count; p++) {
      ExodusFile id(p);

      error += ex_get_ids(id, EX_ELEM_BLOCK, &block_id[0]);

      // Check that the block id ordering is consistent among files...
      if (p > 0) {
	for (size_t b = 0; b < global.count(EBLK); b++) {
	  if (blocks[0][b].id != block_id[b]) {
	    std::cerr << "ERROR: The internal element block id ordering for part "
		      << p << "\n       is not consistent with the ordering for part 0."
		      << std::endl;
	    exit(EXIT_FAILURE);
	  }
	}
      }

      if (debug_level & 4)
	std::cerr << "\nGetting element block info for part " << p
		  << "..." << std::endl;

      for (size_t b = 0; b < global.count(EBLK); b++) {
	if (debug_level & 4)
	  std::cerr << "Block " << b << ", Id = " << block_id[b];

	int temp_elements;
	int temp_node_elements;
	int temp_element_attributes;

	char temp_element_type[MAX_STR_LENGTH+1];

	error += ex_get_elem_block(id, block_id[b], temp_element_type,
				   &temp_elements, &temp_node_elements,
				   &temp_element_attributes);

	std::vector<char> name(ExodusFile::max_name_length()+1);
	ex_get_name(id, EX_ELEM_BLOCK, block_id[b], &name[0]);

	blocks[p][b].id    = block_id[b];
	if (name[0] != '\0')
	  blocks[p][b].name_ = &name[0];
	if (p == 0) {
	  glob_blocks[b].id    = block_id[b];
	  if (name[0] != '\0')
	    glob_blocks[b].name_ = &name[0];
	}

	blocks[p][b].position_         = b;
	if (temp_elements > 0) {
	  blocks[p][b].elementCount    = temp_elements;
	  blocks[p][b].nodesPerElement = temp_node_elements;
	  blocks[p][b].attributeCount  = temp_element_attributes;
	  blocks[p][b].offset_         = temp_elements;
	  strncpy(blocks[p][b].elType,   temp_element_type, MAX_STR_LENGTH);

	  // NOTE: This is not correct, but fixed below.
	  glob_blocks[b].elementCount   += temp_elements;

	  if (glob_blocks[b].nodesPerElement == 0)
	    glob_blocks[b].nodesPerElement = temp_node_elements;
	  else
	    SMART_ASSERT(glob_blocks[b].nodesPerElement == temp_node_elements)
	      (b)(glob_blocks[b].nodesPerElement)(temp_node_elements);

	  if (glob_blocks[b].attributeCount == 0)
	    glob_blocks[b].attributeCount  = temp_element_attributes;
	  else
	    SMART_ASSERT(glob_blocks[b].attributeCount == temp_element_attributes)
	      (b)(glob_blocks[b].attributeCount)(temp_element_attributes);
	  
	  glob_blocks[b].position_  = b;
	  strncpy(glob_blocks[b].elType,   temp_element_type, MAX_STR_LENGTH);
	}

	if (temp_element_attributes > 0 && glob_blocks[b].attributeNames.empty()) {
	  // Get attribute names.  Assume the same on all parts
	  // on which the block exists. 
	  char **names = get_name_array(temp_element_attributes, ExodusFile::max_name_length());

	  ex_get_attr_names(id, EX_ELEM_BLOCK, block_id[b], names);
	  for (int i=0; i < temp_element_attributes; i++) {
	    glob_blocks[b].attributeNames.push_back(std::string(names[i]));
	  }
	  free_name_array(names, temp_element_attributes);
 	}
	if (debug_level & 4) {
	  std::cerr << ", Name = '" << blocks[p][b].name_;
	  std::cerr << "', Elements = " << std::setw(8) << blocks[p][b].entity_count();
	  std::cerr << ", Nodes/element = " << blocks[p][b].nodesPerElement;
	  std::cerr << ", Attributes = " << blocks[p][b].attributeCount << std::endl;
	}
      }
    }	// end for p=0..part_count

    // Convert block_offset from elements/block/part to true offset
    for (size_t p=0; p < part_count; p++) {
      int sum = 0;
      for (size_t b=0; b < global.count(EBLK); b++) {
	int save = blocks[p][b].offset_;
	blocks[p][b].offset_ = sum;
	sum += save;
      }
    }
  }

  template <typename T>
  void put_element_blocks(std::vector<Mesh> &local_mesh,
			  std::vector<std::vector<Block> > &blocks,
			  std::vector<Block> &glob_blocks,
			  T /* dummy */)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    int global_num_blocks = glob_blocks.size();

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    std::cerr << "\nReading and Writing element connectivity & attributes\n";

    for (int b = 0; b < global_num_blocks; b++)  {

      if (debug_level & 4) {
	std::cerr << "\nOutput element block info for...\n";
	std::cerr << "Block " << b << ", Id = " << glob_blocks[b].id;
	std::cerr << ", Name = '" << glob_blocks[b].name_;
	std::cerr << "', Elements = " << std::setw(8) << glob_blocks[b].entity_count();
	std::cerr << ", Nodes/element = " << glob_blocks[b].nodesPerElement;
	std::cerr << ", Attributes = " << glob_blocks[b].attributeCount << std::endl;
      }

      if (debug_level & 4)
	std::cerr << "B" << b << ":\t" << std::flush;


      size_t max_nodes = glob_blocks[b].entity_count();
      max_nodes *= glob_blocks[b].nodesPerElement;
      IntVector block_linkage(max_nodes);

      // Initialize attributes list, if it exists
      std::vector<T> attributes((size_t)glob_blocks[b].attributeCount *
				(size_t)glob_blocks[b].entity_count());

      int error = 0;
      size_t part_count = local_mesh.size();
      for (size_t p = 0; p < part_count; p++) {
	ExodusFile id(p);

	size_t global_pos;
	size_t global_block_pos;

	if (blocks[p][b].entity_count() > 0) {	// non-zero length block
	
	  if (debug_level & 4)
	    std::cerr << "#" << std::flush;
	  size_t maximum_nodes = blocks[p][b].entity_count();
	  maximum_nodes *= blocks[p][b].nodesPerElement;
	  IntVector local_linkage(maximum_nodes);
	
	  int bid = blocks[p][b].id;
	  error = ex_get_conn(id, EX_ELEM_BLOCK, bid, &local_linkage[0], 0, 0);
	  if (error < 0) {
	    std::cerr << "ERROR: Cannot get element block connectivity for block "
		      << bid << " on part " << p << ".\n";
	  }
	  int pos = 0;
	  size_t goffset = glob_blocks[b].offset_;
	  size_t element_count = blocks[p][b].entity_count();
	  size_t boffset = blocks[p][b].offset_;
	  size_t npe = blocks[p][b].nodesPerElement;
	  int* part_loc_elem_to_global = &local_mesh[p].localElementToGlobal[0];
	  int* part_loc_node_to_global = &local_mesh[p].localNodeToGlobal[0];
	
	  for (size_t e = 0; e < element_count; e++) {
	    global_block_pos = part_loc_elem_to_global[(e + boffset)] - goffset;
	    global_pos =  global_block_pos * npe;

	    for (size_t n = 0; n < npe; n++) {
	      size_t node = part_loc_node_to_global[local_linkage[pos++]-1];
	      if (debug_level & 4) {
		SMART_ASSERT(block_linkage[global_pos] == (int)node+1 || block_linkage[global_pos] == 0);
	      }
	      block_linkage[global_pos++] = node+1;
	    }
	  }

	  // Get attributes list,  if it exists
	  if (blocks[p][b].attributeCount > 0) {
	    SMART_ASSERT(blocks[p][b].attributeCount == glob_blocks[b].attributeCount)
	      (p)(b)(blocks[p][b].attributeCount)(glob_blocks[b].attributeCount);
	    
	    size_t max_attr = (size_t)blocks[p][b].entity_count() *
	      blocks[p][b].attributeCount;
	    std::vector<T> local_attr(max_attr);

	    error += ex_get_attr(id, EX_ELEM_BLOCK, blocks[p][b].id, &local_attr[0]);
  	
	    pos = 0;

	    size_t att_count = blocks[p][b].attributeCount;
	    for (size_t e = 0; e < element_count; e++) {
	      // global_pos is global position within this element block...
	      global_block_pos = local_mesh[p].localElementToGlobal[(e + boffset)] - goffset;
	      global_pos = global_block_pos * att_count;
	      for (size_t n = 0; n < att_count; n++) {
		attributes[global_pos++] = local_attr[pos++];
	      }
	    }
	  }
	}	// end if blocks[p][b].entity_count() (non-zero length block)
	else
	  if (debug_level & 4)
	    std::cerr << ".";
      }				// end for p=0..part_count-1

      // Verify that connectivity has been set for all elements...
      for (size_t i=0; i < block_linkage.size(); i++)
	SMART_ASSERT(block_linkage[i] > 0)(block_linkage[i])(i)(b)(glob_blocks[b].id);
      
      // Write out block info
      int id_out = ExodusFile::output();// output file identifier

      if (!block_linkage.empty()) {
	error += ex_put_conn(id_out, EX_ELEM_BLOCK, glob_blocks[b].id, &block_linkage[0], 0, 0);
      }

      // Write out attributes list if it exists
      if (glob_blocks[b].attributeCount > 0) {
	error += ex_put_attr(id_out, EX_ELEM_BLOCK, glob_blocks[b].id, &attributes[0]);
      }				// end for b=0..global_num_blocks-1
      if (debug_level & 4)
	std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }

  void build_reverse_element_map(std::vector<Mesh> &local_mesh,
				 std::vector<std::vector<Block> > &blocks,
				 std::vector<Block> &glob_blocks,
				 Mesh *global, size_t part_count,
				 GlobalElemMap &global_element_map)
  {
    int error = 0;
    // Create the map that maps from a local part element to the
    // global map. This combines the mapping local part element to
    // 'global id' and then 'global id' to global position. The
    // mapping is now a direct lookup instead of a lookup followed by
    // a reverse map.
    //

    // We iterate through each element block a part at a time to
    // determine which elements are in each block. Note that an
    // element will possibly exist in multiple parts, so we need to do
    // the sort/uniqify/shrink on each element block.  The operations
    // for each element block will be:
    // 1. get elements for the block for each part into a vector.
    // 2. sort/uniqify/shrink that vector and copy into the
    // 'global_element_map' in the positions following the previous block.
    // 
    // This is similar to what was done above in a global sense,
    // except it is done an element block at a time. The only need for
    // doing them both is that it lets us detect whether an element id
    // was reused and appears in both blocks.  Currently, just error
    // out if that happens...  In the future, change the id to an
    // unused value and continue....

    size_t tot_size = 0;
    std::vector<GlobalElemMap> global_element_numbers(part_count);
    for (size_t p = 0; p < part_count; p++) {
      ExodusFile id(p);
      global_element_numbers[p].resize(local_mesh[p].count(ELEM));
      IntVector ids(local_mesh[p].count(ELEM));
      error += ex_get_elem_num_map(id, &ids[0]);
      for (size_t i=0; i < local_mesh[p].count(ELEM); i++) {
	global_element_numbers[p][i] = std::make_pair(ids[i], size_t(0));
      }
      tot_size += local_mesh[p].count(ELEM);
    }
    global_element_map.resize(tot_size);

    size_t goffset = 0;
    for (size_t b = 0; b < glob_blocks.size(); b++) {
      size_t block_size = 0;
      for (size_t p = 0; p < part_count; p++) {
	block_size += blocks[p][b].entity_count();
      }
      GlobalElemMap block_element_map(block_size);
      
      size_t poffset = 0;
      for (size_t p = 0; p < part_count; p++) {
	ExodusFile id(p);

	// Get connectivity array for this element block...
	size_t element_count = blocks[p][b].entity_count();
	size_t nodes_per_elem = blocks[p][b].nodesPerElement;
	size_t maximum_nodes = element_count * nodes_per_elem;
	
	IntVector local_linkage(maximum_nodes);

	int bid = blocks[p][b].id;
	error = ex_get_conn(id, EX_ELEM_BLOCK, bid, &local_linkage[0], 0, 0);
	if (error < 0) {
	  std::cerr << "ERROR: Cannot get element block connectivity for block "
		    << bid << " on part " << p << ".\n";
	}

	// Convert connectivity to global node numbers.
	for (size_t i=0; i < element_count * nodes_per_elem; i++) {
	  int local_node = local_linkage[i];
	  int global_node = local_mesh[p].localNodeToGlobal[local_node-1] + 1;
	  local_linkage[i] = global_node;
	}

	// Have element global ids for all elements in this block,
	// and connectivity.  Can now create out "eleminfo" for these elements.
	size_t con_offset = 0;
	size_t boffset = blocks[p][b].offset_;
	for (size_t i=0; i < element_count; i++) {
	  size_t adler_crc = adler(0, &local_linkage[con_offset], nodes_per_elem * sizeof(int));
	  global_element_numbers[p][boffset+i].second = adler_crc;
	  block_element_map[poffset+i] = global_element_numbers[p][boffset+i];
	  con_offset += nodes_per_elem;
	}
	poffset += element_count;
      }

      // Sort, uniqify, shrink 'block_element_map' and the result
      // then contains the list of elements in this block...
      uniqify(block_element_map);
      
      size_t block_total_num_elements = block_element_map.size();
      glob_blocks[b].elementCount = block_total_num_elements;
      glob_blocks[b].offset_ = goffset;
      
      // Copy into the global_element_map...
      std::copy(block_element_map.begin(),
		block_element_map.end(),
		&global_element_map[goffset]);
      goffset += block_total_num_elements;
    }
      
    global->elementCount = goffset;
    global_element_map.resize(goffset);

    size_t max_id = global_element_map[global->elementCount-1].first;
    bool is_contiguous = max_id == global_element_map.size();
    //std::cerr  << "Element id map " << (is_contiguous ? "is" : "is not") << " contiguous.\n";

    // The global_element_map may or may not be globally sorted; however, each
    // block is sorted, so if we do the iteration by blocks, we can
    // use equal_range instead of doing global searches...
    for (size_t b = 0; b < glob_blocks.size(); b++) {
      GElemMapIter gm_begin = global_element_map.begin() + glob_blocks[b].offset_;
      GElemMapIter gm_end   = gm_begin + glob_blocks[b].elementCount;
      GElemMapIter cur_pos  = gm_begin;
      for (size_t p = 0; p < part_count; p++) {
	size_t element_count = blocks[p][b].entity_count();
	size_t boffset = blocks[p][b].offset_;
	for (size_t i = 0; i < element_count; i++) {
	  std::pair<int, size_t> global_element = global_element_numbers[p][boffset+i];
	  
	  if (cur_pos == gm_end || *cur_pos != global_element) {
	    std::pair<GElemMapIter, GElemMapIter> iter =
	      std::equal_range(gm_begin, gm_end, global_element);
	    SMART_ASSERT(iter.first != iter.second);
	    cur_pos = iter.first;
	  }
	  size_t element_value = cur_pos - gm_begin;
	  local_mesh[p].localElementToGlobal[i+boffset] = element_value + glob_blocks[b].offset_;
	  ++cur_pos;
	}				
      }
    }

    // Update the element ids to give a unique, non-repeating set.  If
    // contiguous, then there is nothing to do.  If not contiguous,
    // then need to determine if there are any repeats (id reuse) and
    // if so, generate a new id for the repeated uses.  Note that
    // there is a possibility that elements in two or more element
    // blocks will have the same element id, so we generate a vector
    // containing a pair<id, position_in_global_element_map>, sort it,
    // and then use it to detect duplicates and map them to a new id.

    if (!is_contiguous) {
      std::vector<std::pair<int, int> > id_pos(global_element_map.size());
      for (int i=0; i < global->elementCount; i++) {
	id_pos[i].first  = global_element_map[i].first;
	id_pos[i].second = i;
      }
      std::sort(id_pos.begin(), id_pos.end());
      
      max_id = id_pos[id_pos.size()-1].first;
      // Check again for contiguous ids since we now have a sorted list...
      is_contiguous = max_id == global_element_map.size();
      
      if (!is_contiguous) {
	int repeat_found = 0;
	int id_last = id_pos[0].first;
	SMART_ASSERT(id_last == global_element_map[id_pos[0].second].first);
	
	for (int i=1; i < global->elementCount; i++) {
	  if (id_pos[i].first == id_last) {
	    global_element_map[id_pos[i].second].first = ++max_id;
	    repeat_found++;
	  } else {
	    id_last = id_pos[i].first;
	  }
	}
	if (repeat_found > 0) {
	  std::cerr  << repeat_found << " duplicate element ids were found. Their ids have been renumbered to remove duplicates.\n";
	}
      }
    }
  }

  void build_reverse_node_map(std::vector<Mesh> &local_mesh,
			      Mesh *global, size_t part_count,
			      GlobalMap &global_node_map)
  {
    // Instead of using <set> and <map>, consider using a sorted vector...
    // Append all local node maps to the global node map.
    // Sort the global node map
    // Remove duplicates.
    // Position within map is now the map...
    // When building the local-part node to global id, use binary_search...

    // Global node map and count.
    std::vector<std::vector<NodeInfo> > global_nodes(part_count);

    size_t tot_size = 0;
    for (size_t p = 0; p < part_count; p++) {
      tot_size += local_mesh[p].count(NODE);
      global_nodes[p].resize(local_mesh[p].count(NODE));
    }
    global_node_map.resize(tot_size);

    size_t offset = 0;
    int error = 0;
    for (size_t p = 0; p < part_count; p++) {
      std::vector<double>   x(local_mesh[p].count(NODE));
      std::vector<double>   y(local_mesh[p].count(NODE));
      std::vector<double>   z(local_mesh[p].count(NODE));
      std::vector<int>    nid(local_mesh[p].count(NODE));
      
      ExodusFile id(p);
      error += ex_get_node_num_map(id, &nid[0]);
      error += ex_get_coord(id, &x[0], &y[0], &z[0]);
      for (size_t i=0; i < local_mesh[p].count(NODE); i++) {
	global_nodes[p][i] = NodeInfo(nid[i], x[i], y[i], z[i]);
      }
      std::copy(global_nodes[p].begin(), global_nodes[p].end(),
		&global_node_map[offset]);
      offset += local_mesh[p].count(NODE);
    }

    // Now, sort the global_node_map array and remove duplicates...
    uniqify(global_node_map);
      

    global->nodeCount = global_node_map.size();

    // See whether the node numbers are contiguous.  If so, we can map
    // the nodes back to their original location. Since the nodes are
    // sorted and there are no duplicates, we just need to see if the id
    // at global_node_map.size() == global_node_map.size();
    int max_id = global_node_map[global->nodeCount-1].id;
    bool is_contiguous = max_id == (int)global_node_map.size();
    std::cerr  << "Node map " << (is_contiguous ? "is" : "is not") << " contiguous.\n";

    // Create the map that maps from a local part node to the
    // global map. This combines the mapping local part node to
    // 'global id' and then 'global id' to global position. The
    // mapping is now a direct lookup instead of a lookup followed by
    // a reverse map.
    GMapIter cur_pos = global_node_map.begin();
    for (size_t p = 0; p < part_count; p++) {
      size_t node_count = local_mesh[p].count(NODE);
      for (size_t i = 0; i < node_count; i++) {
	NodeInfo global_node = global_nodes[p][i];

	if (cur_pos == global_node_map.end() || *cur_pos != global_node) {
	  std::pair<GMapIter, GMapIter> iter = std::equal_range(global_node_map.begin(),
								global_node_map.end(),
								global_node);
	  if (iter.first == iter.second) {
	    NodeInfo n = global_node;
	    std::cerr << n.id << "\t" << n.x << "\t" << n.y << "\t" << n.z << "\n";
	    SMART_ASSERT(iter.first != iter.second);
	  }
	  cur_pos = iter.first;
	}
	size_t nodal_value = cur_pos - global_node_map.begin();
	local_mesh[p].localNodeToGlobal[i] = nodal_value;
	++cur_pos;
      }				
    }
    
    // Update the nodal ids to give a unique, non-repeating set.  If contiguous, then
    // there is nothing to do.  If not contiguous, then need to determine if there are any
    // repeats (id reuse) and if so, generate a new id for the repeated uses.
    if (!is_contiguous) {
      bool repeat_found = false;
      int id_last = global_node_map[0].id;
      for (int i=1; i < global->nodeCount; i++) {
	if (global_node_map[i].id == id_last) {
	  global_node_map[i].id = ++max_id;
	  repeat_found = true;
	} else {
	  id_last = global_node_map[i].id;
	}
      }
      if (repeat_found) {
	std::cerr  << "Duplicate node ids were found. Their ids have been renumbered to remove duplicates.\n";
      }
    }
  }

  void get_put_variable_names(int id, int out, Variables &vars, SystemInterface &si,
			      int *combined_status_variable_index)
  {
    if (vars.count(OUT) > 0) {

      char **output_name_list  = get_name_array(vars.count(OUT), ExodusFile::max_name_length());

      int num_input_vars = vars.index_.size();

      char **input_name_list = get_name_array(num_input_vars, ExodusFile::max_name_length());
      ex_get_variable_names(id, vars.type(), num_input_vars, input_name_list);

      /// \todo Check for name collision on status variables.
	if (vars.type() == EX_ELEM_BLOCK || vars.type() == EX_NODAL) {
	  if (vars.type() == EX_ELEM_BLOCK) {
	    std::string status = si.element_status_variable();
	    if (status != "NONE")
	      strcpy(input_name_list[num_input_vars-1], status.c_str());
	  }
	  else if (vars.type() == EX_NODAL) {
	    std::string status = si.nodal_status_variable();
	    if (status != "NONE")
	      strcpy(input_name_list[num_input_vars-1], status.c_str());
	  }
	}

	// Iterate through the 'var_index' and transfer
	// Assume that the number of pointers is limited to
	// the number of results variables
	size_t maxlen = 0;
	for (int i=0; i < num_input_vars; i++) {
	  if (vars.index_[i] > 0) {
	    strcpy(output_name_list[vars.index_[i]-1], input_name_list[i]);
	    if (strlen(input_name_list[i]) > maxlen) {
	      maxlen = strlen(input_name_list[i]);
	    }
	  }
	}
	maxlen += 2;
	// Assume 8 characters for initial tab...
	int width = si.screen_width();
	int nfield = (width-8) / maxlen;
	if (nfield < 1) nfield = 1;

	std::cerr << "Found " << vars.count(OUT) << " " << vars.label() << " variables.\n";
	{
	  int i = 0;
	  int ifld = 1;
	  std::cerr << "\t";
	  while (i < vars.count(OUT)) {
	    std::cerr << std::setw(maxlen) << std::left << output_name_list[i++];
	    if (++ifld > nfield && i < vars.count(OUT)) {
	      std::cerr << "\n\t";
	      ifld = 1;
	    }
	  }
	  std::cerr << "\n\n";
	  std::cerr << std::right; // Reset back to what it was.
	}

	ex_put_variable_names(out, vars.type(), vars.count(OUT), output_name_list);

	// KLUGE: Handle finding combined status index variable here since it is
	// the only place we have the list of output variable names...
	if (vars.type() == EX_ELEM_BLOCK && combined_status_variable_index) {
	  *combined_status_variable_index = 0;
	  std::string comb_stat = si.combined_mesh_status_variable();
	  if (!comb_stat.empty()) {
	    for (int i=0; i <  vars.count(OUT); i++) {
	      if (case_compare(comb_stat, output_name_list[i]) == 0) {
		*combined_status_variable_index = i+1;
		break;
	      }
	    }
	  }
	}
	free_name_array(output_name_list, vars.count(OUT));
	free_name_array(input_name_list, num_input_vars);
    }
  }


  void get_variable_params(int id, Variables &vars, const StringIdVector &variable_list)
  {
    // Determines the number of variables of type 'type()' that will
    // be written to the output database. The 'variable_list' vector
    // specifies a possibly empty list of variable names that the user
    // wants transferred to the output database. If 'variable_list' is
    // empty, then all variables of that type will be transferred; if
    // the 'variable_list' size is 1 and it contains the string 'NONE',
    // then no variables of that type will be transferred; if size is 1
    // and it contains the string 'ALL', then all variables of that type
    // will be transferred.
    //
    // Returns the number of variables which will be output Also creates
    // a 'var_index'.  The var_index is zero-based and of size
    // 'input_variable_count'. If:
    // var_index[i] ==0, variable not written to output database
    // var_index[i] > 0, variable written; is variable 'var_index[i]'

    // If 'type' is ELEMENT or NODE, then reserve space for the 'status' variable.
    int extra = 0;
    if (vars.addStatus)
      extra = 1;

    int num_vars;
    ex_get_variable_param (id, vars.type(),  &num_vars);

    vars.index_.resize(num_vars+extra);

    // Create initial index which defaults to no output...
    std::fill(vars.index_.begin(), vars.index_.end(), 0);
    // ...Except for the status variable (if any)
    if (extra == 1)
      vars.index_[num_vars] = 1;

    // If 'variable_list' is empty or specified 'ALL', then all
    // variables are to be output
    if (variable_list.empty() ||
	(variable_list.size() == 1 && case_compare(variable_list[0].first,"all")==0)) {
      for (size_t i = 0; i < vars.index_.size(); i++) {
	vars.index_[i] = i+1;
      }
      vars.outputCount = num_vars + extra;
      return;
    }

    // Another possibility is user specifies "NONE" for the variable
    // list so no variables will be written.  Just return 0.
    if (variable_list.size() == 1 && case_compare(variable_list[0].first,"none")==0) {
      vars.outputCount = extra;
      return;
    }

    // At this point, the variable_list specifies at least one
    // variable to be output to the database.
    // Get the list that the user entered and the list of files
    // from the input database...
    {
      StringVector exo_names = get_exodus_variable_names(id, vars.type(), num_vars);

      // Iterate 'variable_list' and find position in 'exo_names'.  If
      // not found, there is an error -- set index to -1; otherwise set
      // index to position (1-based) in variable_list.  Others are
      // already set to '0' by default initialization.

      // The variable_list may contain multiple entries for each
      // variable if the user is specifying output only on certain
      // element blocks...
      std::string var_name = "";
      int var_count = 0;
      for (size_t i=0; i < variable_list.size(); i++) {
	if (var_name == variable_list[i].first)
	  continue;
	var_name = variable_list[i].first;
	bool found = false;
	for (size_t j=0; j < exo_names.size() && !found; j++) {
	  if (case_compare(exo_names[j], var_name) == 0) {
	    found = true;
	    vars.index_[j] = ++var_count;
	  }
	}
	if (!found) {
	  std::cerr << "ERROR: Variable '" << variable_list[i].first
		    << "' is not valid." << std::endl;
	  exit(EXIT_FAILURE);
	}
      }
      // Count non-zero entries in var_index;
      int nz_count = 0;
      for (size_t i=0; i < vars.index_.size(); i++) {
	if (vars.index_[i] > 0)
	  nz_count++;
      }
      SMART_ASSERT(nz_count == var_count+extra)(nz_count)(var_count);

      if (vars.addStatus) {
        vars.index_[num_vars] = nz_count; // Already counted above...
      }
      vars.outputCount = nz_count;
      return;
    }
  }


  void put_mesh_summary(const Mesh& mesh)
  {
    // Write out Mesh info

    std::cout << " Title: " << mesh.title.c_str() << "\n\n";
    std::cout << " Number of coordinates per node       =" << std::setw(9)
	      << mesh.count(DIM) << "\n";
    std::cout << " Number of nodes                      =" << std::setw(9)
	      << mesh.count(NODE) << "\n";
    std::cout << " Number of elements                   =" << std::setw(9)
	      << mesh.count(ELEM) << "\n";
    std::cout << " Number of element blocks             =" << std::setw(9)
	      << mesh.count(EBLK) << "\n\n";
    std::cout << " Number of nodal point sets           =" << std::setw(9)
	      << mesh.count(NSET) << "\n";
    std::cout << " Number of element side sets          =" << std::setw(9)
	      << mesh.count(SSET) << "\n\n";
  }

  void get_nodesets(int total_node_count,
		    std::vector<Excn::Mesh> &local_mesh,
		    std::vector<std::vector<NodeSet> > &nodesets,
		    std::vector<NodeSet> &glob_sets)
  {
    // Find number of nodesets in the global model...
    std::set<int> set_ids;
    IntVector ids;

    size_t part_count = local_mesh.size();
    int bad_ns = 0;
    {
      int   ns_count;
      for (size_t p=0; p < part_count; p++) {
	ExodusFile id(p);
	ns_count = ex_inquire_int(id, EX_INQ_NODE_SETS);

	// Get the ids for these
	ids.resize(ns_count);
	ex_get_ids(id, EX_NODE_SET, &ids[0]);

	for (int iset=0; iset < ns_count; iset++) {
	  if (ids[iset] != 0) {
	    set_ids.insert(ids[iset]);
	  } else {
	    bad_ns++;
	  }
	}
      }
    }

    if (bad_ns != 0) {
      std::cerr << "ERROR: There were " << bad_ns
		<< " nodesets (counting all files) which had an id equal to "
		<< "0 which is not allowed.\n";
    }

    if (set_ids.empty())
      return;

    // set_ids now contains the union of all nodeset ids...
    glob_sets.resize(set_ids.size());
    size_t gnset_size = set_ids.size();
    ids.resize(gnset_size);

    {
      std::set<int>::const_iterator I  = set_ids.begin();
      std::set<int>::const_iterator IE = set_ids.end();
      size_t i = 0;
      while (I != IE) {
	glob_sets[i].id = *I++;
	glob_sets[i].position_ = i; i++;
      }
    }
      
    {
      for (size_t p=0; p < part_count; p++) {
	ExodusFile id(p);

	nodesets[p].resize(set_ids.size());

	// Get the ids again so we can map current order back to file order...
	ex_get_ids(id, EX_NODE_SET, &ids[0]);
	int ns_count = ex_inquire_int(id, EX_INQ_NODE_SETS);

	for (int i=0; i < ns_count; i++) {
	  nodesets[p][i].id = ids[i];

	  // Find which global nodeset this id corresponds to...
	  size_t gi = gnset_size;
	  for (size_t j=0; j < gnset_size; j++) {
	    if (ids[i] == glob_sets[j].id) {
	      gi = j;
	      break;
	    }
	  }
	  SMART_ASSERT(gi != gnset_size);
	  glob_sets[gi].position_ = i;
	  nodesets[p][i].position_ = gi;


	  // Get the parameters for this nodeset...
	  ex_get_set_param(id, EX_NODE_SET, nodesets[p][i].id,
			   &nodesets[p][i].nodeCount,
			   &nodesets[p][i].dfCount);

	  std::vector<char> name(ExodusFile::max_name_length()+1);
	  ex_get_name(id, EX_NODE_SET, nodesets[p][i].id, &name[0]);
	  if (name[0] != '\0') {
	    nodesets[p][i].name_ = &name[0];
	    if (glob_sets[gi].name_.empty()) {
	      glob_sets[gi].name_  = &name[0];
	    }
	  }

	  if (debug_level & 32) {
	    std::cerr << "Part " << p+1 << " ";
	    nodesets[p][i].dump();
	  }
	}
      }
    }

    verify_set_position_mapping("nodeset", part_count, glob_sets, nodesets);

    {
      // Now get the nodeset nodes and df.
      // Currently ignore the DF.  Could add if people need it...

      // This is inefficient since the part loop is on
      // the inside...  The other ordering would use more memory...

      IntVector ns_nodes;
      for (size_t ns = 0; ns < set_ids.size(); ns++) {
	
	IntVector glob_ns_nodes(total_node_count+1);
	std::fill(glob_ns_nodes.begin(), glob_ns_nodes.end(), 0);

	size_t lns = glob_sets[ns].position_;
	for (size_t p=0; p < part_count; p++) {
	  ExodusFile id(p);
	
	  glob_sets[ns].name_ = nodesets[p][lns].name_;

	  int size = nodesets[p][lns].entity_count();
	  if (size > 0) {
	    ns_nodes.resize(size);
	    ex_get_set(id, EX_NODE_SET, nodesets[p][lns].id, &ns_nodes[0], 0);

	    // The node ids are in local space -- map to global
	    for (int iset=0; iset < size; iset++) {
	      int global_node = local_mesh[p].localNodeToGlobal[ns_nodes[iset]-1] + 1;
	      glob_ns_nodes[global_node] = 1;
	    }
	  }
	}
	// Count number of nonzero entries and transfer to the
	// output nodeset
	// NOTE: global_node above is 1-based.
	glob_sets[ns].nodeCount = std::accumulate(glob_ns_nodes.begin(), glob_ns_nodes.end(), 0);
	glob_sets[ns].nodeSetNodes.resize(glob_sets[ns].entity_count());
	glob_sets[ns].dfCount = 0;

	size_t j = 0;
	for (int i=1; i <= total_node_count; i++) {
	  if (glob_ns_nodes[i] == 1) {
	    glob_sets[ns].nodeSetNodes[j++] = i;
	    glob_ns_nodes[i] = j;
	  }

	}
	SMART_ASSERT(j == glob_sets[ns].entity_count())(j)(glob_sets[ns].entity_count())(ns);

	// See if we need a map from local nodeset position to global nodeset position
	// Only needed if there are nodeset variables (or dist factors).
	// Assume all files have same number of variables...
	int num_vars;
	{
	  ExodusFile id(0);
	  ex_get_variable_param (id, EX_NODE_SET,  &num_vars);
	}
	if (num_vars > 0) {
	  for (size_t p=0; p < part_count; p++) {
	    ExodusFile id(p);
	    // Get the nodelist, but store it in nodeOrderMap.
	    // global_pos = nodeOrderMap[i]
	    NodeSet &nset = nodesets[p][ns];
	    int nnodes = nset.entity_count();
	    nset.nodeOrderMap.resize(nnodes);
	    ex_get_set(id, EX_NODE_SET, nset.id, &nset.nodeOrderMap[0], NULL);
	    
	    for (int i=0; i < nnodes; i++) {
	      int local_node = nset.nodeOrderMap[i];                       // 1-based
	      int global_node = local_mesh[p].localNodeToGlobal[local_node-1] + 1; // 1-based
	      int global_pos  = glob_ns_nodes[global_node];                // 1-based
	      nset.nodeOrderMap[i] = global_pos - 1;
	    }
#if 0
	    if (debug_level & 32)
	      nset.dump_order();
#endif
	  }
	}
      }
    }
  }


  void put_nodesets(std::vector<NodeSet> &glob_sets)
  {
    int exoid = ExodusFile::output();

    if (debug_level & 32)
      std::cerr << "\nOutput NodeSets:\n";

    for (size_t ns = 0; ns < glob_sets.size(); ns++) {
      ex_put_set(exoid, EX_NODE_SET, glob_sets[ns].id, &glob_sets[ns].nodeSetNodes[0], 0);
      //    ex_put_node_set_dist_fact(exoid, glob_sets[ns].id, &glob_sets[ns].distFactors[0]);
      // Done with the memory; clear out the vector containing the bulk data nodes and distFactors.
      IntVector().swap(glob_sets[ns].nodeSetNodes);
      DistVector().swap(glob_sets[ns].distFactors);

      SMART_ASSERT(glob_sets[ns].nodeSetNodes.empty());
      SMART_ASSERT(glob_sets[ns].distFactors.empty());

      if (debug_level & 32) {
	glob_sets[ns].dump();
      }
    }
  }


  void get_sideset_metadata(std::vector<Excn::Mesh> &local_mesh,
			    std::vector< std::vector<SideSet> > &sets,
			    std::vector<SideSet> &glob_ssets)
  {
    // Find number of sidesets in the global model...
    std::set<int> set_ids;
    IntVector ids;
    
    size_t part_count = local_mesh.size();
    int bad_ss = 0;
    {
      for (size_t p=0; p < part_count; p++) {
	ExodusFile id(p);
	int ss_count = ex_inquire_int(id, EX_INQ_SIDE_SETS);
      
	// Get the ids for these
	ids.resize(ss_count);
	ex_get_ids(id, EX_SIDE_SET, &ids[0]);

	for (int i=0; i < ss_count; i++) {
	  if (ids[i] != 0) {
	    set_ids.insert(ids[i]);
	  } else {
	    bad_ss++;
	  }
	}
      }
    }

    if (bad_ss != 0) {
      std::cerr << "ERROR: There were " << bad_ss
		<< " sidesets (counting all files) which had an id equal to 0 which is not allowed.\n";
    }

    if (set_ids.empty())
      return;

    // set_ids now contains the union of all sideset ids...
    glob_ssets.resize(set_ids.size());
    size_t gsset_size = set_ids.size();
    ids.resize(gsset_size);

    {
      std::set<int>::const_iterator I  = set_ids.begin();
      std::set<int>::const_iterator IE = set_ids.end();
      size_t i = 0;
      while (I != IE) {
	glob_ssets[i].id = *I++;
	glob_ssets[i].position_ = i; i++;
      }
    }
      
    {
      std::vector<char> name(ExodusFile::max_name_length()+1);
      for (size_t p=0; p < part_count; p++) {
	ExodusFile id(p);

	sets[p].resize(set_ids.size());

	// Get the ids again so we can map current order back to file order...
	ex_get_ids(id, EX_SIDE_SET, &ids[0]);

	int ss_count = ex_inquire_int(id, EX_INQ_SIDE_SETS);
	for (int i=0; i < ss_count; i++) {
	  sets[p][i].id = ids[i];

	  // Find which global sideset this id corresponds to...
	  size_t gi = gsset_size;
	  for (size_t j=0; j < gsset_size; j++) {
	    if (ids[i] == glob_ssets[j].id) {
	      gi = j;
	      break;
	    }
	  }
	  SMART_ASSERT(gi != gsset_size);
	  glob_ssets[gi].position_ = i;
	  sets[p][i].position_ = gi;
	  
	  // Get the parameters for this sideset...
	  ex_get_set_param (id, EX_SIDE_SET, sets[p][i].id,
			    &sets[p][i].sideCount,
			    &sets[p][i].dfCount);



	  glob_ssets[gi].sideCount += sets[p][i].entity_count();
	  glob_ssets[gi].dfCount   += sets[p][i].dfCount;

	  ex_get_name(id, EX_SIDE_SET, sets[p][i].id, &name[0]);
	  if (name[0] != '\0') {
	    sets[p][i].name_ = &name[0];
	    glob_ssets[gi].name_ = &name[0];
	  }
	}
      }

      verify_set_position_mapping("sideset", part_count, glob_ssets, sets);

      // See if we need a map from local sideset position to global sideset position
      // Only needed if there are sideeset variables (or dist factors which are currently ignored).
      // Assume all files have same number of variables...
      bool need_sideset_map = false;
      {
	int num_vars;
	ExodusFile id(0);
	ex_get_variable_param (id, EX_SIDE_SET,  &num_vars);
	need_sideset_map = num_vars > 0;
      }

      // Now get the sideset elements/sides. Ignoring DF for now...
      for (size_t ss = 0; ss < set_ids.size(); ss++) {
	// This is maximum possible size; will probably be reduced by duplicate elimination...
	typedef std::vector<std::pair<int,int> > ElemSideMap;
	typedef ElemSideMap::iterator     ESMapIter;
	ElemSideMap elem_side(glob_ssets[ss].sideCount);
	int ss_id = glob_ssets[ss].id;
	size_t offset = 0;
	
	size_t lss = glob_ssets[ss].position_;
	for (size_t p=0; p < part_count; p++) {

	  ExodusFile id(p);
	  sets[p][lss].elems.resize(sets[p][lss].sideCount);
	  sets[p][lss].sides.resize(sets[p][lss].sideCount);
	  ex_get_set(id, EX_SIDE_SET, ss_id,
		     &sets[p][lss].elems[0], &sets[p][lss].sides[0]);

	  // Add these to the elem_side vector...
	  for (int i=0; i < sets[p][lss].sideCount; i++) {
	    int global_elem = local_mesh[p].localElementToGlobal[sets[p][lss].elems[i]-1]+1;
	    elem_side[offset+i] = std::make_pair(global_elem, sets[p][lss].sides[i]);
	  }
	  offset += sets[p][lss].sideCount;
	}

	uniqify(elem_side);

	// Set the output sideset definition...
	glob_ssets[ss].sideCount = elem_side.size();	
	glob_ssets[ss].dfCount   = 0;
	glob_ssets[ss].elems.resize(elem_side.size());	
	glob_ssets[ss].sides.resize(elem_side.size());	

	// Populate the global sideset elements and sides...
	for (size_t i=0; i < elem_side.size(); i++) {
	  glob_ssets[ss].elems[i] = elem_side[i].first;
	  glob_ssets[ss].sides[i] = elem_side[i].second;
	}

	if (need_sideset_map) {
	  // For this sideset in each part, figure out the mapping for
	  // its (elem, side, variable) position into the corresponding
	  // global sideset position...

	  // Try the equal_range searching of elem_side for now.  If
	  // inefficient, fix later...
	  for (size_t p=0; p < part_count; p++) {
	    sets[p][lss].elemOrderMap.resize(sets[p][lss].sideCount);
	    for (int i=0; i < sets[p][lss].sideCount; i++) {
	      int global_elem = local_mesh[p].localElementToGlobal[sets[p][lss].elems[i]-1]+1;
	      std::pair<int,int> es = std::make_pair(global_elem, sets[p][lss].sides[i]);

	      std::pair<ESMapIter, ESMapIter> iter =
		std::equal_range(elem_side.begin(), elem_side.end(), es);
	      SMART_ASSERT(iter.first != iter.second);
	      size_t pos = iter.first - elem_side.begin();
	      sets[p][lss].elemOrderMap[i] = pos;
	    }
	  }
	}
      }

      // Calculate sideset offset
      for (size_t b=0; b < glob_ssets.size(); b++) {
	int sum = 0;
	for (size_t p=0; p < part_count; p++) {
	  sets[p][b].offset_ = sum;
	  sum += sets[p][b].entity_count();

	  if (debug_level & 16) {
	    std::cerr << "Part " << p+1 << " ";
	    sets[p][b].dump();
	  }
	}
      }

      // Free some memory which is no longer needed...
      // (Could move up into sideset loop above)
      for (size_t p=0; p < part_count; p++) {
	for (size_t ss=0; ss < sets[p].size(); ss++) {
	  IntVector().swap(sets[p][ss].elems);
	  IntVector().swap(sets[p][ss].sides);
	  DistVector().swap(sets[p][ss].distFactors);
	  SMART_ASSERT(sets[p][ss].elems.empty());
	  SMART_ASSERT(sets[p][ss].sides.empty());
	  SMART_ASSERT(sets[p][ss].distFactors.empty());
	}
      }
    }
  }

  void get_put_sidesets(std::vector<Excn::Mesh> &local_mesh,
			std::vector< std::vector<SideSet> > &sets,
			std::vector<SideSet> &glob_ssets)
  {
    int exoid = ExodusFile::output();// output file identifier
    for (size_t ss = 0; ss < glob_ssets.size(); ss++) {
      ex_put_set(exoid, EX_SIDE_SET, glob_ssets[ss].id, const_cast<int*>(&glob_ssets[ss].elems[0]),
		 const_cast<int*>(&glob_ssets[ss].sides[0]));
      if (glob_ssets[ss].dfCount > 0) {
	ex_put_set_dist_fact(exoid, EX_SIDE_SET, glob_ssets[ss].id,
			     reinterpret_cast<void*>(&glob_ssets[ss].distFactors[0]));
      }
    }

    for (size_t ss = 0; ss < glob_ssets.size(); ss++) {
      IntVector().swap(glob_ssets[ss].elems);
      IntVector().swap(glob_ssets[ss].sides);
      DistVector().swap(glob_ssets[ss].distFactors);
      SMART_ASSERT(glob_ssets[ss].elems.empty());
      SMART_ASSERT(glob_ssets[ss].sides.empty());
      SMART_ASSERT(glob_ssets[ss].distFactors.empty());
    }
  }

  template <typename T>
  void add_status_variable(int id_out, const Mesh& global,
			   const std::vector<Block> &blocks,
			   const std::vector<Block> &glob_blocks,
			   const Excn::IntVector &local_element_to_global,
			   int step, int variable, T alive,
			   int combined_variable_index)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    std::vector<T> status;
    SMART_ASSERT(alive == 0.0 || alive == 1.0)(alive);

    for (size_t b = 0; b < global.count(EBLK); b++) {
      status.resize(glob_blocks[b].entity_count());
      std::fill(status.begin(), status.end(), (1.0 - alive));
      int boffset = blocks[b].offset_;
      int goffset = glob_blocks[b].offset_;
      int element_count = blocks[b].entity_count();
      for (int e = 0; e < element_count; e++) {
	int global_block_pos = local_element_to_global[(e + boffset)] - goffset;
	status[global_block_pos] = alive;
      }

      // If combining this status variable with a mesh status
      // variable, do that now...
      if (combined_variable_index > 0) {
	std::vector<T> mesh_status(glob_blocks[b].entity_count());
	ex_get_var(id_out, step, EX_ELEM_BLOCK, combined_variable_index,
		   glob_blocks[b].id, glob_blocks[b].entity_count(),
		   &mesh_status[0]);
      
	if (alive == 1.0) {
	  for (size_t i=0; i < glob_blocks[b].entity_count(); i++) {
	    status[i] = status[i] * mesh_status[i];
	  }
	} else {
	  SMART_ASSERT(alive == 0.0);
	  for (size_t i=0; i < glob_blocks[b].entity_count(); i++) {
	    status[i] = 1.0 - ((1.0-status[i]) * (1.0 - mesh_status[i]));
	  }
	}
      }
      ex_put_var(id_out, step, EX_ELEM_BLOCK, variable,
                 glob_blocks[b].id, glob_blocks[b].entity_count(),
                 &status[0]);
    }
  }

  StringVector get_exodus_variable_names(int id, ex_entity_type elType, int var_count)
  {
    // Allocate space for variable names...
    char **name_list = get_name_array(var_count, ExodusFile::max_name_length());
    ex_get_variable_names(id, elType, var_count, name_list);

    StringVector names(var_count);
    for (int j=0; j < var_count; j++) {
      compress_white_space(name_list[j]);
      names[j] = std::string(name_list[j]);
    }

    free_name_array(name_list, var_count);
    return names;
  }

  template <typename T>
  void filter_truth_table(int id,
			  Mesh& global,
			  std::vector<T> &glob_blocks,
			  Variables &vars,
			  const StringIdVector &variable_names)
  {
    // This routine checks the 'variable_names' list to see if the user
    // has restricted the output of certain variables to certain element
    // blocks. If so, then the truth table is modified to match the
    // users request.
    if (variable_names.size() == 0)
      return;

    // Check for a non-zero id entry in the variable_names list which
    // will signify a user-specified element block.
    bool found_it = false;
    for (size_t i=0; i < variable_names.size() && !found_it; i++) {
      if (variable_names[i].second > 0)
	found_it = true;
    }
    if (!found_it)
      return;

    // At this point, know that there is at least one block-restricted
    // variable output specification. For each one, modify the global
    // truth table to match the specification.
    StringVector exo_names = get_exodus_variable_names(id, vars.type(), vars.count());

    std::string var_name = "";
    int out_position = -1;
    int inp_position = -1;
    for (size_t i=0; i < variable_names.size(); i++) {
      if (variable_names[i].second > 0) {
	if (var_name != variable_names[i].first) {
	  var_name = variable_names[i].first;
	  // Find which exodus variable matches this name
	  inp_position = -1;
	  out_position = -1;
	  for (size_t j = 0; j < exo_names.size(); j++) {
	    if (case_compare(exo_names[j], var_name) == 0) {
	      inp_position = j;
	      out_position = vars.index_[j]-1;
	      break;
	    }
	  }
	  SMART_ASSERT(inp_position >= 0)(inp_position);
	  SMART_ASSERT(out_position >= 0)(out_position);

	  // Set all truth table entries for this variable to negative
	  // of current value and then iterate over specified blocks and
	  // set those positive.  This way can make sure that the
	  // variable truly exists for the block that the user specified.
	  for (size_t b=0; b < global.count(vars.objectType); b++) {
	    int truth_table_loc = (b * vars.count(OUT)) + out_position;
	    global.truthTable[vars.objectType][truth_table_loc] *= -1;
	  }
	}
	// Find out which block corresponds to the specified id.
	int block = -1;
	for (size_t b = 0; b < global.count(vars.objectType); b++) {
	  if (glob_blocks[b].id == variable_names[i].second) {
	    block = b;
	    break;
	  }
	}

	if (block == -1) {
	  std::cerr << "ERROR: User-specified block id of "
		    << variable_names[i].second
		    << " for variable '" << variable_names[i].first
		    << "' does not exist.\n";
	  exit(EXIT_FAILURE);
	}

	int truth_table_loc = block * vars.count(OUT) + out_position;
	if (global.truthTable[vars.objectType][truth_table_loc] == 0) {
	  std::cerr << "ERROR: Variable '" << variable_names[i].first
		    << "' does not exist on block " << variable_names[i].second
		    << ".\n";
	  exit(EXIT_FAILURE);
	} else {
	  global.truthTable[vars.objectType][truth_table_loc] = 1;
	}
      }
    }

    // reset truth table values that may be negative
    int output_truth_table_length = vars.count(OUT) * global.count(vars.objectType);
    for (int j = 0; j < output_truth_table_length; j++) {
      if (global.truthTable[vars.objectType][j] < 0)
	global.truthTable[vars.objectType][j] = 0;
    }
  }

  template <typename T>
  void get_truth_table(Mesh& global, std::vector<T> &glob_blocks, std::vector<Mesh> &local,
		       Variables &vars, int debug)
  {
    // read truth table - sum across all parts since many will
    // have values set to zero for zero length blocks the element
    // variable truth table organized as a 2D array:
    // [global.count(EBLK)][num_elem_vars]

    ObjectType object_type = vars.objectType;
    int input_truth_table_length  = vars.count(IN)  * global.count(object_type);
    int output_truth_table_length = vars.count(OUT) * global.count(object_type);

    if (output_truth_table_length) {

      global.truthTable[object_type].resize(output_truth_table_length);
      std::fill(global.truthTable[object_type].begin(), global.truthTable[object_type].end(), 0);

      // For each input exodus file, get it's truth table and fill
      // in the location in the output truth table...

      bool is_sidenodeset = vars.objectType==NSET || vars.objectType==SSET;
      size_t part_count = local.size();
      for (size_t p = 0; p < part_count; p++) {
	ExodusFile id(p);

	if (vars.count(IN) > 0) {
	  local[p].truthTable[object_type].resize(input_truth_table_length);
	  ex_get_truth_table(id, vars.type(), global.count(object_type), vars.count(IN),
			     &local[p].truthTable[object_type][0]);
	}
	for (size_t b=0; b < global.count(object_type); b++) {
	  int bin = b;
	  if (is_sidenodeset) {
	    bin = glob_blocks[b].position_;
	  }

	  for (int j = 0; j < vars.count(IN); j++) {
	    if (vars.index_[j] > 0) {
	      int ki = (bin * vars.count(IN))  + j;
	      int ko = (b   * vars.count(OUT)) + vars.index_[j] - 1;
	      SMART_ASSERT(ko < output_truth_table_length);
	      SMART_ASSERT(ki < input_truth_table_length);
	      global.truthTable[object_type][ko] += local[p].truthTable[object_type][ki];
	    }
	  }
          if (vars.addStatus) {
            int ko = (b * vars.count(OUT)) + vars.count(OUT) - 1;
            SMART_ASSERT(ko < output_truth_table_length);
            global.truthTable[object_type][ko] = 1;
          }
	}
      }

      // reset truth table values that may be greater than 1
      for (int j = 0; j < output_truth_table_length; j++) {
	SMART_ASSERT(global.truthTable[object_type][j] >= 0)(global.truthTable[object_type][j]);
	if (global.truthTable[object_type][j] > 0)
	  global.truthTable[object_type][j] = 1;
      }

      if (debug_level & debug) {
	std::cerr << "Truth table for " << vars.label() << "\t"
		  << vars.count(OUT) << " variables\t"
		  << global.count(object_type) << " sets\n";
	int k = 0;
	for (size_t b=0; b < global.count(object_type); b++) {
	  for (int j = 0; j < vars.count(OUT); j++) {
	    std::cerr << global.truthTable[object_type][k++];
	  }
	  std::cerr << '\n';
	}
      }
    }
  }

  int case_compare(const std::string &s1, const std::string &s2)
  {
    const char *c1 = s1.c_str();
    const char *c2 = s2.c_str();
    for (;;) {
      if (::toupper(*c1) != ::toupper(*c2)) {
	return (::toupper(*c1) - ::toupper(*c2));
      }
      if (*c1 == '\0') {
	return 0;
      }
      c1++;
      c2++;
    }
  }

  void add_info_record(char *info_record, int size)
  {
    // Add 'uname' output to the passed in character string.
    // Maximum size of string is 'size' (not including terminating NULL)
    // This is used as information data in the concatenated results file
    // to help in tracking when/where/... the file was created
    struct utsname sys_info;
    uname(&sys_info);

    std::string info = "CONJOIN: ";
    info += sys_info.nodename;
    info += ", OS: ";
    info += sys_info.sysname;
    info += " ";
    info += sys_info.release;
    info += ", ";
    info += sys_info.version;
    info += ", Machine: ";
    info += sys_info.machine;
    const char *sinfo = info.c_str();
    strncpy(info_record, sinfo, size);
    info_record[size] = '\0';
  }

  inline bool is_whitespace(char c)
  {
    static char white_space[] = {' ', '\t', '\n', '\r', ',', '\0'};
    return (std::strchr(white_space, c) != NULL);
  }

  void compress_white_space(char *str)
  {
    char *ibuf = str;
    char *obuf = str;

    int i = 0;
    int cnt = 0;

    // Don't process an empty string.
    if (str == NULL)
      return;

    // Skip leading...
    while(*ibuf != 0 && is_whitespace(*ibuf))
      ++ibuf;

    while(*ibuf != 0) {
      if(is_whitespace(*ibuf) && cnt > 0)
	ibuf++;
      else {
	if (!is_whitespace(*ibuf))
	  cnt = 0;
	else {
	  *ibuf = ' ';
	  cnt = 1;
	}
	obuf[i++] = *ibuf++;
      }
    }
    obuf[i--] = '\0';

    // Skip trailing whitespace
    while (i > 0 && is_whitespace(obuf[i]))
      obuf[i--] = '\0';
  }

  std::string time_stamp(const std::string &format)
  {
    if (format == "") {
      return std::string("");
    } else {
      const int length=256;
      static char time_string[length];

      time_t calendar_time = time(NULL);
      struct tm *local_time = localtime(&calendar_time);

      int error = strftime(time_string, length, format.c_str(), local_time);
      if (error != 0) {
        time_string[length-1] = (char)NULL;
        return std::string(time_string);
      } else {
        return std::string("[ERROR]");
      }
    }
  }

  std::string format_time(double seconds)
  {
    char suffix = 'u';
    if (seconds > 0.0 && seconds < 1.0) {
      return " <1s";
    }

    if (seconds > 86400) {
      suffix = 'd';
      seconds /= 86400.;
    } else if (seconds > 3600) {
      suffix = 'h';
      seconds /= 3600.;
    } else if (seconds >   60) {
      suffix = 'm';
      seconds /= 60.;
    } else {
      suffix = 's';
    }
    std::ostringstream os;
    os << std::showpoint << std::setprecision(2) << seconds << suffix;
    return os.str();
  }

  int get_width(int max_value)
  {
    // Returns the field width which will accommodate the
    // largest value.
    int width = 0;
    if (max_value >= 10)
      width = int(log10((double)max_value));
    return width + 1;
  }

  template <typename T>
  void map_element_vars(int loffset, int goffset, int entity_count, 
			std::vector<T> &values,
			std::vector<T> &global_values,
			int *part_loc_elem_to_global)
  {
    // copy values to master element value information
    T* local_values = &values[0];
    for (int j = 0; j < entity_count; j++) {
      int global_block_pos = part_loc_elem_to_global[(j + loffset)] - goffset;
      global_values[global_block_pos] = local_values[j];
    }
  }

  template <typename T, typename U>
  void map_sideset_vars(U &, int, int, std::vector<T> &, std::vector<T> &)
  {
    SMART_ASSERT(1==0 && "Internal Error!");
  }

  void map_sideset_vars(Excn::SideSet &local_set, int entity_count, int glob_entity_count,
			std::vector<double> &values,
			std::vector<double> &global_values)
  {
    // copy values to master nodeset value information
    double* local_values = &values[0];
    for (int j = 0; j < entity_count; j++) {
      int global_loc = local_set.elemOrderMap[j];
      SMART_ASSERT(global_loc >= 0 && global_loc < glob_entity_count);
      global_values[global_loc] = local_values[j];
    }
  }

  void map_sideset_vars(Excn::SideSet &local_set, int entity_count, int glob_entity_count,
			std::vector<float> &values,
			std::vector<float> &global_values)
  {
    // copy values to master nodeset value information
    float* local_values = &values[0];
    for (int j = 0; j < entity_count; j++) {
      int global_loc = local_set.elemOrderMap[j];
      SMART_ASSERT(global_loc >= 0 && global_loc < glob_entity_count);
      global_values[global_loc] = local_values[j];
    }
  }

  template <typename T, typename U>
  void map_nodeset_vars(U &, int, int, std::vector<T> &, std::vector<T> &)
  {
    SMART_ASSERT(1==0 && "Internal Error!");
  }

  void map_nodeset_vars(Excn::NodeSet &local_set, int entity_count, int glob_entity_count,
			std::vector<double> &values,
			std::vector<double> &global_values)
  {
    // copy values to master nodeset value information
    double* local_values = &values[0];
    for (int j = 0; j < entity_count; j++) {
      int global_loc = local_set.nodeOrderMap[j];
      SMART_ASSERT(global_loc >= 0 && global_loc < glob_entity_count);
      global_values[global_loc] = local_values[j];
    }
  }

  void map_nodeset_vars(Excn::NodeSet &local_set, int entity_count, int glob_entity_count,
			std::vector<float> &values, 
			std::vector<float> &global_values)
  {
    // copy values to master nodeset value information
    float* local_values = &values[0];
    for (int j = 0; j < entity_count; j++) {
      int global_loc = local_set.nodeOrderMap[j];
      SMART_ASSERT(global_loc >= 0 && global_loc < glob_entity_count);
      global_values[global_loc] = local_values[j];
    }
  }

  template <typename T, typename U>
  int read_write_master_values(Excn::Variables &vars, const Excn::Mesh& global,
			       std::vector<U> &global_sets,
			       std::vector<Excn::Mesh> &local_mesh,
			       std::vector<std::vector<U> > &local_sets,
			       std::vector<T> &values, 
			       size_t p, int time_step, int time_step_out)
  {
    int error = 0;
    bool is_sidenodeset = vars.objectType==NSET || vars.objectType==SSET;

    ExodusFile id(p);
    int id_out = ExodusFile::output();// output file identifier

    size_t max_size = 0;
    for (size_t b = 0; b < global.count(vars.objectType); b++) {
      if (max_size < global_sets[b].entity_count())
	max_size = global_sets[b].entity_count();
    }
    std::vector<T> master_values(max_size);
    
    // Only needed for element, but haven't cleaned this up yet...
    int* part_loc_elem_to_global = &local_mesh[p].localElementToGlobal[0];
	
    for (int i = 0; i < vars.count(IN); i++) {
      if (vars.index_[i] > 0) {
	int ivar = vars.index_[i]-1;
	
	for (size_t b = 0; b < global.count(vars.objectType); b++) {
	  int bin = b;
	  if (is_sidenodeset) {
	    bin = global_sets[b].position_;
	  }
	  
	  int output_truth_table_loc = (b   * vars.count(OUT)) + ivar;
	  int input_truth_table_loc  = (bin * vars.count(IN))  + ivar;
	  if (global.truthTable[vars.objectType][output_truth_table_loc] &&
	      local_sets[p][bin].entity_count() > 0) {
	    
	    int entity_count = local_sets[p][bin].entity_count();
	    
	    if (local_mesh[p].truthTable[vars.objectType][input_truth_table_loc] > 0) {
	      error += ex_get_var(id, time_step+1, exodus_object_type(vars.objectType),
				  i+1, local_sets[p][bin].id, entity_count, &values[0]);
	      
	      switch (vars.objectType) {
	      case EBLK:
		map_element_vars(local_sets[p][bin].offset_, global_sets[b].offset_,
				 entity_count, values, master_values, part_loc_elem_to_global);
		break;
		
	      case SSET:
		map_sideset_vars(local_sets[p][bin], entity_count, global_sets[b].entity_count(),
				 values, master_values);
		break;
		
	      case NSET:
		map_nodeset_vars(local_sets[p][bin], entity_count, global_sets[b].entity_count(),
				 values, master_values);
		break;
	      default:
		break;
	      }
	      
	    }
	    ex_put_var(id_out, time_step_out, exodus_object_type(vars.objectType),
		       ivar+1, global_sets[b].id, global_sets[b].entity_count(),
		       &master_values[0]);
	  }
	}
      }
    }
    return error;
  }

  template <typename U>
  void create_output_truth_table(const Excn::Mesh& global,
				 std::vector<U> &global_sets,
				 Excn::Variables &vars,
				 std::vector<int> &truth_table)
  {
    for (size_t b=0; b < global.count(vars.objectType); b++) {
      int bout = global_sets[b].position_;
      SMART_ASSERT(bout >= 0);
      for (int j = 0; j < vars.count(OUT); j++) {
	int inp_ttable_loc = (b    * vars.count(OUT)) + j;
	int out_ttable_loc = (bout * vars.count(OUT)) + j;
	truth_table[out_ttable_loc] = global.truthTable[vars.objectType][inp_ttable_loc];
      }
    }
  }

  size_t find_max_entity_count(size_t part_count,
			       std::vector<Excn::Mesh> &local_mesh,
			       const Excn::Mesh& global,
			       std::vector<std::vector<Block> > &blocks,
			       std::vector<std::vector<NodeSet> > &nodesets,
			       std::vector<std::vector<SideSet> > &sidesets)
  {
    size_t max_ent = local_mesh[0].count(NODE);
    for (size_t p = 1; p < part_count; p++) {
      if ((size_t)local_mesh[p].count(NODE) > max_ent)
	max_ent = local_mesh[p].count(NODE);
    }
    
    for (size_t p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(EBLK); b++) {
	if (blocks[p][b].entity_count() > max_ent)
	  max_ent = blocks[p][b].entity_count();
      }
    }
    
    // Nodesets...
    for (size_t p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(NSET); b++) {
	if (nodesets[p][b].entity_count() > max_ent)
	  max_ent = nodesets[p][b].entity_count();
      }
    }
    
    // Sidesets...
    for (size_t p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(SSET); b++) {
	if (sidesets[p][b].entity_count() > max_ent)
	  max_ent = sidesets[p][b].entity_count();
      }
    }
    return max_ent;
  }
}
