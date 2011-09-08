/*
 * Copyright(C) 2010 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
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
// concatenates EXODUS/GENESIS output from parallel processors to a single file

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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>

#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <ctype.h>
#include "smart_assert.h"
#include <sys/utsname.h>

//typedef std::vector<int>    IntVector;
typedef std::vector<std::string> StringVector;
//typedef std::vector<std::pair<std::string, int> > StringIdVector;

typedef std::vector<int>    GlobalMap;
typedef GlobalMap::iterator GMapIter;

#include <Internals.h>
#include <ExodusFile.h>
#include <ExodusEntity.h>
#include <SystemInterface.h>
#include <Version.h>
#include <Variables.h>
#include <ObjectType.h>

#include <exodusII.h>

#if EX_API_VERS_NODOT <= 467
#error "Requires exodusII version 4.68 or later"
#endif

extern void add_to_log(const char *name, int elapsed);
extern double epu_timer();
namespace {
  std::string time_stamp(const std::string &format);
  std::string format_time(double seconds);
  int get_width(int max_value);

  ex_entity_type exodus_object_type(Excn::ObjectType &epu_type) {
    switch (epu_type) {
    case Excn::EBLK:
      return EX_ELEM_BLOCK;
    case Excn::SSET:
      return EX_SIDE_SET;
    case Excn::NSET:
      return EX_NODE_SET;
    default:
      SMART_ASSERT(1==0 && "Invalid Object Type in exodus_object_type")(epu_type);
      return EX_INVALID;
    }
  }

  char **get_name_array(int size, int length) {
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
int epu(Excn::SystemInterface &interface, int start_part, int part_count, int cycle,
	T /* dummy */);

namespace {
void compress_white_space(char *str);
void add_info_record(char *info_record, int size);
void put_global_info(const Excn::Mesh& global);
void get_put_coordinate_frames(int id, int id_out);
void get_put_qa(int exo_id, int out_id);
int  get_put_coordinate_names(int in, int out, int dimensionality);

template <typename T>
int get_put_coordinates(Excn::Mesh& global, int part_count,
			std::vector<Excn::Mesh> &local_mesh,
			int **local_node_to_global, T dummy);

template <typename T>
int get_coordinates(int id, int dimensionality, int num_nodes,
		    int **local_node_to_global,
		    int proc, std::vector<T> &x, std::vector<T> &y,
		    std::vector<T> &z);

StringVector get_exodus_variable_names(int id, ex_entity_type elType,
				       int count);

template <typename T>
void filter_truth_table(int id,
			Excn::Mesh& global,
			std::vector<T> &glob_blocks,
			Excn::Variables &vars,
			const Excn::StringIdVector &variable_names);

template <typename T>
void get_truth_table(Excn::Mesh& global, std::vector<T> &glob_blocks,
		     std::vector<Excn::Mesh> &local, Excn::Variables &vars, int debug);

template <typename U>
void create_output_truth_table(const Excn::Mesh& global,
			       std::vector<U> &global_sets,
			       Excn::Variables &vars,
			       std::vector<int> &truth_table);

template <typename T, class U>
void clear_master_values(Excn::Variables &vars, const Excn::Mesh& global,
			 std::vector<U> &glob_sets, T*** master_values);

template <typename T, typename U>
int read_master_values(Excn::Variables &vars, const Excn::Mesh& global,
		       std::vector<U> &global_sets,
		       std::vector<Excn::Mesh> &local_mesh,
		       std::vector<std::vector<U> > &local_sets,
		       T*** master_values, std::vector<T> &values, 
		       int part_count, int step, 
		       int **local_element_to_global);

template <typename T, typename U>
void output_master_values(Excn::Variables &vars, const Excn::Mesh& global,
			  std::vector<U> &glob_sets,
			  T*** master_values, int step);
  
template <typename T, typename U>
void allocate_master_values(Excn::Variables &vars, const Excn::Mesh& global,
			    std::vector<U> &glob_sets, T***& master_values);
  
template <typename T>
void deallocate_master_values(Excn::Variables &vars, const Excn::Mesh& global,
			      T***& master_values);
  
void get_variable_params(int id, Excn::Variables &vars,
			 const Excn::StringIdVector &variable_list);

void get_put_variable_names(int id, int idout, Excn::Variables &vars,
			    Excn::SystemInterface &interface);

void build_reverse_element_map(int **local_element_to_global,
			       const std::vector<Excn::Mesh> &local_mesh,
			       std::vector<std::vector<Excn::Block> > &blocks,
			       std::vector<Excn::Block> &glob_blocks,
			       Excn::Mesh *global, int part_count,
			       GlobalMap &global_element_map,
			       bool map_ids);

void get_nodesets(int part_count, int total_node_count,
		  int **local_node_to_global,
		  std::vector<std::vector<Excn::NodeSet> > &nodesets,
		  std::vector<Excn::NodeSet> &glob_sets);

void build_reverse_node_map(int **local_node_to_global,
			    const std::vector<Excn::Mesh> &local_mesh,
			    Excn::Mesh *global, int part_count,
			    GlobalMap &global_node_map);

void get_element_blocks(int part_count,
			const std::vector<Excn::Mesh> &local_mesh,
			const Excn::Mesh &global,
			std::vector<std::vector<Excn::Block> > &blocks,
			std::vector<Excn::Block> &glob_blocks);
template<typename T>
void put_element_blocks(int part_count, int start_part,
			std::vector<std::vector<Excn::Block> > &blocks,
			std::vector<Excn::Block> &glob_blocks,
			int **local_node_to_global,
			int **local_element_to_global,
			T single_or_double);

void put_nodesets(std::vector<Excn::NodeSet> &glob_sets);

void get_sideset_metadata(int part_count,
			  std::vector< std::vector<Excn::SideSet> > &sets,
			  std::vector<Excn::SideSet> &glob_ssets);
void get_put_sidesets(int part_count, int **local_element_to_global,
		      std::vector< std::vector<Excn::SideSet> > &sets,
		      std::vector<Excn::SideSet> &glob_ssets,
		      Excn::SystemInterface &interface);
template <typename T>
void add_processor_variable(int id_out, int part_count, int start_part,
			    const Excn::Mesh& global,
			    std::vector<std::vector<Excn::Block> > &blocks,
			    const std::vector<Excn::Block> &glob_blocks,
			    int **local_element_to_global,
			    int step, int variable, std::vector<T> &proc);

size_t find_max_entity_count(int part_count,
			     std::vector<Excn::Mesh> &local_mesh,
			     const Excn::Mesh& global,
			     std::vector<std::vector<Excn::Block> > &blocks,
			     std::vector<std::vector<Excn::NodeSet> > &nodesets,
			     std::vector<std::vector<Excn::SideSet> > &sidesets);
  
int case_compare(const std::string &s1, const std::string &s2);
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

  //   1 -- time stamp
  //   2 -- check nodal variable consistency
  //   4 -- Element Blocks
  //   8 -- Nodes
  //  16 -- Sidesets
  //  32 -- Nodesets
  //  64 -- exodus verbose.

  debug_level = interface.debug();

  if (debug_level & 64)
    ex_opts(EX_VERBOSE|EX_DEBUG);
  else
    ex_opts(0);

  int start_part      = interface.start_part();
  int processor_count = interface.processor_count();

  int part_count    = interface.part_count();
  if (part_count <= 1) {
    std::cerr << "INFO: Only one processor or part, no concatenation needed.\n";
    exit(EXIT_SUCCESS);
  }

  if (interface.subcycle() >= 0) {
    start_part = 0;
    int cycles = interface.subcycle();
    if (cycles > 0) {
      // use the specified number of cycles...
      part_count = (processor_count + cycles - 1) / cycles;
    }

    // Sanity check...
    if (part_count < 1) {
      std::cerr << "ERROR: The subcycle specification results in less than 1 part per cycle which is not allowd.\n";
      exit(EXIT_FAILURE);
    }
    interface.subcycle((processor_count+part_count-1) / part_count);

    if (start_part + part_count > processor_count) {
      part_count = processor_count - start_part;
    }
  }

  int error = 0;
  int cycle = 0;
  while (start_part < processor_count) {
    
    if (start_part + part_count > processor_count) {
      part_count = processor_count - start_part;
    }
    
    SMART_ASSERT(part_count > 0);
    SMART_ASSERT(start_part + part_count <= processor_count);

    if (!ExodusFile::initialize(interface, start_part, part_count)) {
      std::cerr << "ERROR: Problem initializing input and/or output files.\n";
      exit(EXIT_FAILURE);
    }

    if (ExodusFile::io_word_size() == 4) {
      error = epu(interface, start_part, part_count, cycle++, (float)0.0);
    } else {
      error = epu(interface, start_part, part_count, cycle++, (double)0.0);
    }

    start_part += part_count;
    ExodusFile::close_all();
    if (!interface.subcycle())
      break;
  }

  time_t end_time = time(NULL);
  add_to_log(argv[0], (int)(end_time-begin_time));
  return (error);
  }
  catch (std::exception &e) {
    std::cerr << "ERROR: Standard exception: " << e.what() << std::endl;
  }
}

template <typename T>
int epu(SystemInterface &interface, int start_part, int part_count, int cycle, T /* dummy */)
{
  SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());

  int p; // file counter p=0..part_count-1
  int j; // generic counter

  char* mytitle = new char[MAX_LINE_LENGTH + 1];
  memset(mytitle, '\0', MAX_LINE_LENGTH+1);

  Mesh global;

  // contains the global node information from each file/processor
  int** local_node_to_global    = new int*[part_count];

  // contains the global element information from each file/processor
  int** local_element_to_global = new int*[part_count];

  std::vector<Mesh> local_mesh(part_count);

  // ******************************************************************
  // 1. Read global info

  int error = 0;

  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);
  std::cerr << "\n**** READ LOCAL (GLOBAL) INFO ****" << std::endl;

  std::string title0;

  // If there are any processors with zero nodes, then
  // that node won't have any nodal variables defined.  We need to
  // find the first processor which has a non-zero node count to use
  // when we look for the nodal variable count.
  int non_zero_node_count = -1;
  for (p = 0; p < part_count; p++) {
    error += ex_get_init(ExodusFile(p), mytitle,
			 &local_mesh[p].dimensionality,
                         &local_mesh[p].nodeCount,
			 &local_mesh[p].elementCount,
			 &local_mesh[p].blockCount,
			 &local_mesh[p].nodesetCount,
			 &local_mesh[p].sidesetCount);

    local_mesh[p].title = mytitle;

    if (local_mesh[p].nodeCount > 0 && non_zero_node_count == -1)
      non_zero_node_count = p;
    
    if (p == 0) {
      global.title = mytitle;
      global.dimensionality = local_mesh[p].dimensionality;
      global.blockCount     = local_mesh[p].count(EBLK);
      global.nodesetCount   = local_mesh[p].count(NSET);
      global.sidesetCount   = local_mesh[p].count(SSET);
    } else {
      SMART_ASSERT(global.dimensionality == local_mesh[p].dimensionality);
      SMART_ASSERT(global.count(EBLK)    == local_mesh[p].count(EBLK));
      SMART_ASSERT(global.count(NSET)    == local_mesh[p].count(NSET));
      SMART_ASSERT(global.count(SSET)    == local_mesh[p].count(SSET));
    }

    local_node_to_global[p]    = new int[local_mesh[p].nodeCount];
    local_element_to_global[p] = new int[local_mesh[p].elementCount];

    // sum required data
    // note that num_blocks is the same for every processor
    global.elementCount += local_mesh[p].elementCount;

  } // end for (p=0..part_count)

  if (non_zero_node_count == -1) {
    non_zero_node_count = 0; // No nodes on entire model...
  }

  delete [] mytitle;

  if (interface.omit_nodesets()) {
    global.nodesetCount = 0;
  }

  if (interface.omit_sidesets()) {
    global.sidesetCount = 0;
  }

  // Need these throughout run, so declare outside of this block...
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
    build_reverse_node_map(local_node_to_global, local_mesh,
			   &global, part_count, global_node_map);

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    std::cerr << "Finished reading/writing Global Info\n";

    // ****************************************************************************
    // 5. Get Block information including element attributes
    // must check for zero length blocks
    get_element_blocks(part_count, local_mesh, global, blocks, glob_blocks);

    GlobalMap global_element_map(global.elementCount);
    build_reverse_element_map(local_element_to_global, local_mesh, blocks,
			      glob_blocks, &global, part_count, global_element_map,
			      interface.map_element_ids());


    //
    //    NOTE:  Node set/side set information can be different for each processor
    /************************************************************************/
    // 7. Get Side sets
    if (!interface.omit_sidesets()) {
      if (debug_level & 1)
	std::cerr << time_stamp(tsFormat);
      std::cerr << "\n**** GET SIDE SETS *****\n";
      get_sideset_metadata(part_count, sidesets, glob_ssets);
      if (global.count(SSET) != glob_ssets.size()) {
	std::cerr << "\nWARNING: Invalid sidesets will not be written to output database.\n";
	global.sidesetCount = glob_ssets.size();
      }
    }


    /************************************************************************/
    // 6. Get Node sets
    if (!interface.omit_nodesets()) {
      if (debug_level & 1)
	std::cerr << time_stamp(tsFormat);
      std::cerr << "\n**** GET NODE SETS *****\n";
      get_nodesets(part_count, global.nodeCount, local_node_to_global, nodesets, glob_nsets);
      if (global.count(NSET) != glob_nsets.size()) {
	std::cerr << "\nWARNING: Invalid nodesets will not be written to output database.\n";
	global.nodesetCount = glob_nsets.size();
      }
    }

    /************************************************************************/
    // Start writing the output file...

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    std::cerr << "\n**** BEGIN WRITING OUTPUT FILE *****\n";
    CommunicationMetaData comm_data;

    // Create the output file...
    ExodusFile::create_output(interface, cycle);

    // Define metadata for model....
    put_global_info(global);

    Internals exodus(ExodusFile::output(), ExodusFile::max_name_length());

    if (interface.append()) {
      bool matches = exodus.check_meta_data(global, glob_blocks, glob_nsets, glob_ssets, comm_data);
      if (!matches) {
	std::cerr << "\n\nERROR: Current mesh dimensions do not match "
		  << "the mesh dimensions in the file being appended to.\n\n";
	exit(EXIT_FAILURE);
      }
    } else {
      exodus.write_meta_data(global, glob_blocks, glob_nsets, glob_ssets, comm_data);

      // Output bulk mesh data....
      put_nodesets(glob_nsets);

      // c.2.  Write Global Node Number Map
      if (debug_level & 1)
	std::cerr << time_stamp(tsFormat);
      std::cerr << "Writing global node number map...\n";
      error=ex_put_node_num_map(ExodusFile::output(),&global_node_map[0]);
    
      if (debug_level & 1)
	std::cerr << time_stamp(tsFormat);
      std::cerr << "Writing out master global elements information...\n";
      if (global_element_map.size() > 0) {
	ex_put_elem_num_map(ExodusFile::output(), &global_element_map[0]);
      }

      // Needed on glory writing to Lustre or we end up with empty maps...
      ex_update(ExodusFile::output());

      T dummy = 0.0;
      put_element_blocks(part_count, start_part, blocks, glob_blocks,
			 local_node_to_global, local_element_to_global, dummy);
    }

    get_put_sidesets(part_count, local_element_to_global, sidesets, glob_ssets, interface);
  }				
  // ************************************************************************
  // 2. Get Coordinate Info.
  if (!interface.append()) {
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    std::cerr << "\n\n**** GET COORDINATE INFO ****\n";

    error += get_put_coordinates(global, part_count, local_mesh,
				 local_node_to_global, (T)0.0);

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    std::cerr << "Wrote coordinate information...\n";
  }
  // ####################TRANSIENT DATA SECTION###########################
  // ***********************************************************************
  // 9. Get Variable Information and names

  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);
  std::cerr << "\n**** GET VARIABLE INFORMATION AND NAMES ****" << std::endl;

  //  I. read number of variables for each type.
  //  NOTE: it is assumed that every processor has the same global, nodal,
  //        and element lists

  Variables global_vars(GLOBAL);
  Variables nodal_vars(NODE);
  Variables element_vars(EBLK);
  Variables nodeset_vars(NSET);
  Variables sideset_vars(SSET);
  
  element_vars.addProcessorId = interface.add_processor_id_field();

  {
    ExodusFile id(non_zero_node_count);
    
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

    if (!interface.append()) {
      ex_put_all_var_param(ExodusFile::output(),
			   global_vars.count(OUT),
			   nodal_vars.count(OUT),
			   element_vars.count(OUT), &elem_truth_table[0],
			   nodeset_vars.count(OUT), &global.truthTable[NSET][0],
			   sideset_vars.count(OUT), &global.truthTable[SSET][0]);
    }
  }

  // II. read/write the variable names
  {
    ExodusFile id(non_zero_node_count);
    get_put_variable_names(id, ExodusFile::output(), global_vars,  interface);
    get_put_variable_names(id, ExodusFile::output(), nodal_vars,   interface);
    get_put_variable_names(id, ExodusFile::output(), element_vars, interface);
    get_put_variable_names(id, ExodusFile::output(), nodeset_vars, interface);
    get_put_variable_names(id, ExodusFile::output(), sideset_vars, interface);
  }
  if (!interface.append())
    ex_update(ExodusFile::output());

  /**********************************************************************/
  // 10. Get Transient Data
  //     This routine reads in a time dump from an EXODUSII file

  int time_step;
  int num_time_steps = 0;

  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);
  std::cerr << "\n**** GET TRANSIENT NODAL, GLOBAL, AND ELEMENT DATA VALUES ****\n";

  // Stage I: Get the number_of_time_steps information

  bool differ = false;
  for (p = 0; p < part_count; p++) {
    ExodusFile id(p);

    int nts = ex_inquire_int(id, EX_INQ_TIME);
    if (p == 0) {
      num_time_steps = nts;
    } else {
      if (nts != num_time_steps) {
	differ = true;
      }
      num_time_steps = num_time_steps < nts ? num_time_steps : nts;
    }
  }
  if (differ) {
    std::cerr << "\nWARNING: The number of time steps is not the same on all input databases.\n"
	      << "         Using minimum count of " << num_time_steps << "\n\n";
  } else {
    std::cerr << "\nNumber of time steps on input databases = " << num_time_steps << "\n\n";
  }

  std::vector<T> global_values(global_vars.count(IN));
  std::vector<T> output_global_values(global_vars.count(OUT));

  T **master_nodal_values   = new T*[nodal_vars.count(OUT)];
  for (int i = 0; i < nodal_vars.count(OUT); i++) {
    master_nodal_values[i] = new T[global.nodeCount];
  }

  // TODO: Handle variables via a class instead of 3-D array.
  T*** master_element_values;
  allocate_master_values(element_vars, global, glob_blocks, master_element_values);
  
  T*** master_sideset_values;
  allocate_master_values(sideset_vars, global, glob_ssets,  master_sideset_values);

  T*** master_nodeset_values;
  allocate_master_values(nodeset_vars, global, glob_nsets,  master_nodeset_values);

  // Determine maximum number of entities on any processor...
  int max_ent = find_max_entity_count(part_count, local_mesh, global, blocks, nodesets, sidesets);
  std::vector<T> values(max_ent);
  
  // Stage II.  Extracting transient variable data.
  //            loop over time steps

  if (num_time_steps == 0 && element_vars.add_processor_id()) {
    // Add a fake timestep with just the processor id information.
    // If adding the processor_id field, do it here...
    T time_val = 0.0;
    ex_put_time(ExodusFile::output(), 1, &time_val);

    std::vector<T> proc;
    add_processor_variable(ExodusFile::output(), part_count, start_part,
			   global, blocks, glob_blocks,
			   local_element_to_global,
			   1, element_vars.index_[element_vars.count()], proc);
  }

  // Determine if user wants a subset of timesteps transferred to the output file.
  int ts_min  = interface.step_min();
  int ts_max  = interface.step_max();
  int ts_step = interface.step_interval();

  if (ts_min == -1 && ts_max == -1) {
    ts_min = num_time_steps;
    ts_max = num_time_steps;
  }
  
  // Time steps for output file
  int time_step_out = 0;
  int time_step_out_init = 0;
  double min_time_to_write = -1.0;

  if (interface.append()) {
    // See how many steps already exist on the output database
    // and the corresponding time.
    int nstep = ex_inquire_int(ExodusFile::output(), EX_INQ_TIME);  
    
    // Get the time corresponding to this step...
    ex_get_time(ExodusFile::output(), nstep, &min_time_to_write);

    time_step_out = nstep;
    time_step_out_init = nstep;
  }

  ts_max = ts_max < num_time_steps ? ts_max : num_time_steps;
  if (ts_min <= ts_max) {
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);
    std::cerr << "\tTransferring step " << ts_min << " to step "
	      << ts_max << " by " << ts_step << std::endl;
  }

  // Determine how many steps will be written...
  int output_steps = (ts_max - ts_min) / ts_step + 1;

  double start_time = epu_timer();

  for (time_step = ts_min-1; time_step < ts_max; time_step+=ts_step) {
    time_step_out++;

    T time_val     = 0.0;
    {
      // read in and write out the time step information
      ExodusFile id(0);

      error += ex_get_time(id,     time_step+1,   &time_val);
      if (time_val <= min_time_to_write)
	continue;
      
      if (min_time_to_write >= 0.0) {
	std::cerr << "\tAppend Mode: Skipping " << time_step - (ts_min-1)
		  << " input steps to align times with already written steps on output file.\n\n";
	min_time_to_write = -1.0;
      }

      error += ex_put_time(ExodusFile::output(), time_step_out, &time_val);

      for (p = 1; p < part_count; p++) {
	ExodusFile idp(p);
	T proc_time_val = 0.0;
	error += ex_get_time(idp, time_step+1, &proc_time_val);
	if (proc_time_val != time_val) {
	  std::cerr << "ERROR: At step " << std::setw(get_width(ts_max+1)) << time_step+1
		    << ", the time on processor " << 0 + start_part << " is "
		    << std::setw(15) << std::scientific << std::setprecision(8)
		    << time_val << " which does not\n       match the time on processor "
		    << p + start_part << " which is " 
		    << std::setw(15) << std::scientific << std::setprecision(8)
		    << proc_time_val << "\n       This usually indicates a corrupt database."
		    << std::endl;
	}
      }

      // NOTE: Assuming that each processor has the exact same global
      // information
      if (global_vars.count(OUT) > 0) {
	if (debug_level & 1)
	  std::cerr << time_stamp(tsFormat) << "Global Variables...\n";
	error += ex_get_var(id, time_step+1, EX_GLOBAL, 0, 0,
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
      
	// Check global variable consistency...
	if (debug_level & 128) {
	  std::vector<T> proc_global_values(global_vars.count(IN));
	  for (p = 1; p < part_count; p++) {
	    ExodusFile idp(p);
	    error += ex_get_var(idp, time_step+1, EX_GLOBAL, 0, 0,
				global_vars.count(IN), (void*)&proc_global_values[0]);
	    for (int ig=0; ig < global_vars.count(IN); ig++) {
	      if (proc_global_values[ig] != global_values[ig]) {
		std::cerr << "At step " << std::setw(get_width(ts_max+1)) << time_step+1
			  << ", Global Variable " << std::setw(get_width(global_vars.count(IN))) << ig+1
			  << ", P" << std::setfill('0') << std::setw(get_width(interface.processor_count()))
			  << 0 + start_part << " = "  << std::setfill(' ')  << std::setw(15)
			  << std::scientific << std::setprecision(8)
			  << global_values[ig]
			  << ", P" << std::setfill('0') << std::setw(get_width(interface.processor_count()))
			  << p + start_part << " = "  << std::setfill(' ')  << std::setw(15)
			  << std::scientific << std::setprecision(8)
			  << proc_global_values[ig]
			  << std::endl;
	      }
	    }
	  }
	}
      }
    }

    // ========================================================================
    // Nodal Values...
    if (debug_level & 1) {
      std::cerr << time_stamp(tsFormat) << "Nodal Variables...\n";
    }
    if (debug_level & 2) {
      for (int i=0; i < nodal_vars.count(OUT); i++) {
	std::fill(&master_nodal_values[i][0],
		  &master_nodal_values[i][global.nodeCount],
		  0.0);
      }
    }

    if (nodal_vars.count(OUT) > 0) {
      for (p = 0; p < part_count; p++) {
	ExodusFile id(p);

	int node_count = local_mesh[p].nodeCount;
	for (int i = 0; i < nodal_vars.count(IN); i++) {
	  if (nodal_vars.index_[i] > 0) {
	    error += ex_get_var(id, time_step+1, EX_NODAL, i+1, 0,
				node_count, &values[0]);
	    
	    int i_out = nodal_vars.index_[i]-1;
	    SMART_ASSERT(i_out < nodal_vars.count(OUT));
	    if (debug_level & 2) {
	      for (j = 0; j < node_count; j++) {
		int nodal_value = local_node_to_global[p][j];
		if (master_nodal_values[i_out][nodal_value] != 0 &&
		    master_nodal_values[i_out][nodal_value] != values[j]) {
		  std::cerr << "Node " << nodal_value
			    << ", old = " << master_nodal_values[i_out][nodal_value]
			    << ", new = " << values[j] << std::endl;
		}
	      }
	    }
	    
	    T *local_nodal_values = &values[0];
	    T *global_nodal_values = &master_nodal_values[i_out][0];
	    if (interface.sum_shared_nodes()) {
	      // sum values into master nodal value information. Note
	      // that for non-shared nodes, this will be the same as a
	      // copy; for shared nodes, it will be a true sum.
	      for (j = 0; j < node_count; j++) {
		// Map local nodal value to global location...
		int nodal_value = local_node_to_global[p][j];
		global_nodal_values[nodal_value] += local_nodal_values[j];
	      }
	    } else {
	      // copy values to master nodal value information
	      for (j = 0; j < node_count; j++) {
		// Map local nodal value to global location...
		int nodal_value = local_node_to_global[p][j];
		global_nodal_values[nodal_value] = local_nodal_values[j];
	      }
	    }
	  }
	}
      }
      // output nodal variable info. for specified time step
      for (int i = 0; i < nodal_vars.count(OUT); i++) {
	error+=ex_put_var(ExodusFile::output(), time_step_out, EX_NODAL, i+1, 0,
			  global.nodeCount, &master_nodal_values[i][0]);
      }
    }
    
    // ========================================================================
    // Extracting element transient variable data
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat) << "Element Variables...\n";

    if (debug_level & 4){
      clear_master_values(element_vars, global, glob_blocks, master_element_values);
    }

    
    if (element_vars.count(IN) > 0) {
      read_master_values(element_vars, global, glob_blocks, local_mesh, blocks,
			 master_element_values, values,
			 part_count, time_step, local_element_to_global);
      output_master_values(element_vars, global, glob_blocks, master_element_values, time_step_out);
    }

    // If adding the processor_id field, do it here...
    // Use the output time step for writing data
    if (element_vars.add_processor_id()) {
      std::vector<T> proc;
      add_processor_variable(ExodusFile::output(), part_count, start_part,
			     global, blocks, glob_blocks,
			     local_element_to_global,
			     time_step_out,
			     element_vars.index_[element_vars.count(IN)], proc);
    }

    // ========================================================================
    // Extracting sideset transient variable data
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat) << "Sideset Variables...\n";

    if (debug_level & 16){
      clear_master_values(sideset_vars, global, glob_ssets, master_sideset_values);
    }

    if (sideset_vars.count(IN) > 0) {
      read_master_values(sideset_vars, global, glob_ssets, local_mesh, sidesets,
			 master_sideset_values, values,
			 part_count, time_step, local_element_to_global);

      output_master_values(sideset_vars, global, glob_ssets, master_sideset_values, time_step_out);
    }

    // ========================================================================
    // Extracting nodeset transient variable data
    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat) << "Nodeset Variables...\n";

    if (debug_level & 32){
      clear_master_values(nodeset_vars, global, glob_nsets, master_nodeset_values);
    }

    if (nodeset_vars.count(IN) > 0) {
      read_master_values(nodeset_vars, global, glob_nsets, local_mesh, nodesets,
			 master_nodeset_values, values,
			 part_count, time_step, local_element_to_global);

      output_master_values(nodeset_vars, global, glob_nsets, master_nodeset_values, time_step_out);
    }

    // ========================================================================
    ex_update(ExodusFile::output());

    if (debug_level & 1)
      std::cerr << time_stamp(tsFormat);

    std::cerr << "Wrote step " << std::setw(2) << time_step+1 << ", time "
	      << std::scientific << std::setprecision(4) << time_val;

    double cur_time = epu_timer();
    double elapsed = cur_time - start_time;
    double time_per_step = elapsed / time_step_out;
    double percentage_done = (time_step_out * 100.0) / output_steps;
    double estimated_remaining = time_per_step * (output_steps - time_step_out);
    std::cerr << "\t\t["
	      << std::fixed << std::setw(5) << std::setprecision(1)
	      << percentage_done << "%, Elapsed=" << format_time(elapsed)
	      << ", ETA=" << format_time(estimated_remaining)
	      << "]\n";
    if (debug_level & 1) {
      std::cerr << "\n";
    }
  }


  for (int n = 0; n < nodal_vars.count(OUT); n++) {
    delete [] master_nodal_values[n];
  }
  delete [] master_nodal_values;
  
  deallocate_master_values(element_vars, global, master_element_values);
  deallocate_master_values(sideset_vars, global, master_sideset_values);
  deallocate_master_values(nodeset_vars, global, master_nodeset_values);

  for (p=0; p < part_count; p++) {
    delete [] local_element_to_global[p];
    delete [] local_node_to_global[p];
  }

  delete [] local_element_to_global;
  delete [] local_node_to_global;

  /*************************************************************************/
  // EXIT program
  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);
  std::cerr << "******* END *******\n";
  return(error);
}

namespace {
  void get_put_coordinate_frames(int id, int id_out)
  {
    // NOTE: Assuming info and QA records for all processors
    // See if there are any coordinate frames...
    int num_frames = ex_inquire_int(id, EX_INQ_COORD_FRAMES);
    if (num_frames <= 0)
      return;

    IntVector           ids(num_frames);
    std::vector<double> coordinates(9*num_frames);
    std::vector<char>   tags(num_frames);

    ex_get_coordinate_frames(id, &num_frames, &ids[0], &coordinates[0], &tags[0]);

    // Now output to the combined file...
    ex_put_coordinate_frames(id_out, num_frames, &ids[0], &coordinates[0], &tags[0]);
  }

  void get_put_qa(int id, int id_out)
  {
    // NOTE: Assuming info and QA records for all processors
    int error = 0;

    // I. Get and store info strings, if they exist
    int num_info_records = ex_inquire_int(id,EX_INQ_INFO);
    char **info_records = new char*[num_info_records+1];
    int info_string_len = MAX_LINE_LENGTH;

    {
      for (int i = 0; i < num_info_records+1; i++) {
	info_records[i] = new char[info_string_len + 1];
	memset(info_records[i], '\0', info_string_len+1);
      }
    }

    if (num_info_records > 0) {
      error += ex_get_info(id, info_records);
    }

    // Add an info record for EPU
    add_info_record(info_records[num_info_records], MAX_LINE_LENGTH);

    error += ex_put_info(id_out,num_info_records+1,info_records);

    {
      for (int i = 0; i < num_info_records+1; i++) {
	delete [] info_records[i];
      }
      delete [] info_records;
    }

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
  int get_put_coordinates(Mesh& global, int part_count,
			  std::vector<Mesh> &local_mesh,
			  int **local_node_to_global, T /* dummy */)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    std::vector<T> x(global.nodeCount);
    std::vector<T> y(global.nodeCount);
    std::vector<T> z(global.nodeCount);

    if (debug_level & 8) {
      std::fill(x.begin(), x.end(), FILL_VALUE);
      std::fill(y.begin(), y.end(), FILL_VALUE);
      std::fill(z.begin(), z.end(), FILL_VALUE);
    }

    int error = 0;
    for (int p = 0; p < part_count; p++) {
      error += get_coordinates(ExodusFile(p), global.dimensionality,
			       local_mesh[p].nodeCount,
			       local_node_to_global, p, x, y, z);
    }                            // end for p=0..part_count

    // Get Coordinate Names
    // NOTE: Assuming coordinate names should be
    // the same for all files/processors therefore, only one
    // file/processor needs to be loaded
    error += get_put_coordinate_names(ExodusFile(0),
				      ExodusFile::output(),
				      global.dimensionality);
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
  int get_coordinates(int id, int dimensionality, int num_nodes,
		      int **local_node_to_global,
		      int proc, std::vector<T> &x, std::vector<T> &y, std::vector<T> &z)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    int error = 0;
    std::vector<T> local_x(num_nodes);
    std::vector<T> local_y(num_nodes);
    std::vector<T> local_z(num_nodes);

    error += ex_get_coord(id, &local_x[0], &local_y[0], &local_z[0]);

    // Check for 2D or 3D coordinates

    if (dimensionality  == 3) {
      if (debug_level & 8) {
	for (int i = 0; i < num_nodes; i++) {
	  int node = local_node_to_global[proc][i];
	  if (x[node] != FILL_VALUE && y[node] != FILL_VALUE && z[node] != FILL_VALUE) {
	    if (x[node] != local_x[i] || y[node] != local_y[i] || z[node] != local_z[i]) {
	      std::cerr << "\nWARNING: Node " << node+1
			<< " has different coordinates in at least two files.\n"
			<< "         cur value = "
			<< std::scientific << std::setprecision(6)
			<< std::setw(14) << x[node] << std::setw(14) << y[node] << std::setw(14) << z[node] << "\n"
			<< "         new value = "
			<< std::setw(14) << local_x[i] << std::setw(14) << local_y[i] << std::setw(14) << local_z[i]
			<< " from processor " << proc << std::endl;

	    }
	  }
	}
      }

      for (int i = 0; i < num_nodes; i++) {

	// The following WILL overwrite x[node],y[node],z[node] if node is the
	// same for different processors
	int node = local_node_to_global[proc][i];
	x[node] = local_x[i];
	y[node] = local_y[i];
	z[node] = local_z[i];

      }
    }
    else {
      if (debug_level & 8) {
	for (int i = 0; i < num_nodes; i++) {
	  int node = local_node_to_global[proc][i];
	  if (x[node] != FILL_VALUE && y[node] != FILL_VALUE) {
	    if (x[node] != local_x[i] || y[node] != local_y[i]) {
	      std::cerr << "\nWARNING: Node " << node+1
			<< " has different coordinates in at least two files.\n"
			<< "         cur value = "
			<< std::scientific << std::setprecision(6)
			<< std::setw(14) << x[node] << std::setw(14) << y[node] << "\n"
			<< "         new value = "
			<< std::setw(14) << local_x[i] << std::setw(14) << local_y[i]
			<< " from processor " << proc << std::endl;
	    }
	  }
	}
      }
      for (int i = 0; i < num_nodes; i++) {

	// The following WILL overwrite x[node],y[node] if node is the same for
	// different processors
	int node = local_node_to_global[proc][i];

	x[node] = local_x[i];
	y[node] = local_y[i];
      }
    }
    return error;
  }

  void get_element_blocks(int part_count,
			  const std::vector<Mesh> &local_mesh,
			  const Mesh &global,
			  std::vector<std::vector<Block> > &blocks,
			  std::vector<Block> &glob_blocks)
  {
    std::cerr << "\n\n**** GET BLOCK INFORMATION (INCL. ELEMENT ATTRIBUTES) ****\n";;

    for (int ip =0; ip < part_count; ip++)
      blocks[ip].resize(local_mesh[ip].count(EBLK));

    std::cerr << "Global block count = " << global.count(EBLK) << std::endl;

    IntVector block_id(global.count(EBLK));

    int error = 0;
    for (int p = 0; p < part_count; p++) {
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
	std::cerr << "\nGetting element block info for processor " << p
		  << "..." << std::endl;
      else {
	if (p == 0)
	  std::cerr << "\nGetting element block info.\n";
      }

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

	std::vector<char> name(Excn::ExodusFile::max_name_length()+1);
	ex_get_name(id, EX_ELEM_BLOCK, block_id[b], &name[0]);

	blocks[p][b].id    = block_id[b];
	if (name[0] != '\0')
	  blocks[p][b].name_ = &name[0];
	if (p == 0) {
	  glob_blocks[b].id    = block_id[b];
	  if (name[0] != '\0')
	    glob_blocks[b].name_ = &name[0];
	}

	if (temp_elements != 0) {
	  blocks[p][b].elementCount    = temp_elements;
	  blocks[p][b].nodesPerElement = temp_node_elements;
	  blocks[p][b].attributeCount  = temp_element_attributes;
	  blocks[p][b].offset_         = temp_elements;
	  blocks[p][b].position_  = b;
	  strncpy(blocks[p][b].elType,   temp_element_type, MAX_STR_LENGTH);

	  glob_blocks[b].elementCount   += temp_elements;
	  glob_blocks[b].nodesPerElement = temp_node_elements;
	  glob_blocks[b].attributeCount  = temp_element_attributes;
	  glob_blocks[b].position_  = b;
	  strncpy(glob_blocks[b].elType,   temp_element_type, MAX_STR_LENGTH);
	}

	if (temp_element_attributes > 0 && glob_blocks[b].attributeNames.empty()) {
	  // Get attribute names.  Assume the same on all processors
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

    // Convert block_offset from elements/block/processor to true offset
    for (int p=0; p < part_count; p++) {
      int sum = 0;
      for (size_t b=0; b < global.count(EBLK); b++) {
	int save = blocks[p][b].offset_;
	blocks[p][b].offset_ = sum;
	sum += save;
      }
    }
  }

  template <typename T>
  void put_element_blocks(int part_count, int start_part, 
			  std::vector<std::vector<Block> > &blocks,
			  std::vector<Block> &glob_blocks,
			  int **local_node_to_global,
			  int **local_element_to_global, T /* dummy */)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    int global_num_blocks = glob_blocks.size();

    int** linkage  = new int*[global_num_blocks];
    T** attributes = new T*[global_num_blocks];

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

      if (max_nodes > 0) {
	linkage[b] = new int[max_nodes];
      } else {
	linkage[b] = NULL;
      }
      int *block_linkage = linkage[b];

      // Initialize attributes list, if it exists
      if (glob_blocks[b].attributeCount > 0) {
	attributes[b]=new T[(size_t)glob_blocks[b].attributeCount *
			    (size_t)glob_blocks[b].entity_count()];
      }

      int error = 0;
      for (int p = 0; p < part_count; p++) {
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
		      << bid << " on part " << p+start_part << ".\n";
	  }
	  int pos = 0;
	  size_t goffset = glob_blocks[b].offset_;
	  size_t element_count = blocks[p][b].entity_count();
	  size_t boffset = blocks[p][b].offset_;
	  size_t npe = blocks[p][b].nodesPerElement;
	  int* proc_loc_elem_to_global = local_element_to_global[p];
	  int* proc_loc_node_to_global = local_node_to_global[p];
	
	  for (size_t e = 0; e < element_count; e++) {
	    global_block_pos = proc_loc_elem_to_global[(e + boffset)] - goffset;
	    global_pos =  global_block_pos * npe;

	    for (size_t n = 0; n < npe; n++) {
	      size_t node = proc_loc_node_to_global[local_linkage[pos++]-1];
	      block_linkage[global_pos++] = node+1;
	    }
	  }

	  // Get attributes list,  if it exists
	  if (blocks[p][b].attributeCount > 0) {
	    size_t max_attr = (size_t)blocks[p][b].entity_count() *
	      blocks[p][b].attributeCount;
	    std::vector<T> local_attr(max_attr);

	    error += ex_get_attr(id, EX_ELEM_BLOCK, blocks[p][b].id, &local_attr[0]);
  	
	    pos = 0;

	    size_t att_count = blocks[p][b].attributeCount;
	    for (size_t e = 0; e < element_count; e++) {
	      // global_pos is global position within this element block...
	      global_block_pos = local_element_to_global[p][(e + boffset)] - goffset;
	      global_pos = global_block_pos * att_count;
	      for (size_t n = 0; n < att_count; n++) {
		attributes[b][global_pos++] = local_attr[pos++];
	      }
	    }
	  }

	}	// end if blocks[p][b].entity_count() (non-zero length block)
	else
	  if (debug_level & 4)
	    std::cerr << ".";
      }				// end for p=0..part_count-1

      // Write out block info
      int id_out = ExodusFile::output();// output file identifier

      if (linkage[b] != NULL) {
	error += ex_put_conn(id_out, EX_ELEM_BLOCK, glob_blocks[b].id, linkage[b], 0, 0);
	delete [] linkage[b];
      }

      // Write out attributes list if it exists
      if (glob_blocks[b].attributeCount > 0) {
	error += ex_put_attr(id_out, EX_ELEM_BLOCK, glob_blocks[b].id,attributes[b]);
	delete [] attributes[b];
      }				// end for b=0..global_num_blocks-1
      if (debug_level & 4)
	std::cerr << std::endl;
    }
    std::cerr << std::endl;
    delete [] linkage;
    delete [] attributes;
  }

  void build_reverse_element_map(int **local_element_to_global,
				 const std::vector<Mesh> &local_mesh,
				 std::vector<std::vector<Block> > &blocks,
				 std::vector<Block> &glob_blocks,
				 Mesh *global, int part_count,
				 GlobalMap &global_element_map, bool map_ids)
  {
    // Global element map and count.
    std::vector<std::vector<int> > global_element_numbers(part_count);

    size_t tot_size = 0;
    for (int p = 0; p < part_count; p++) {
      tot_size += local_mesh[p].elementCount;
      global_element_numbers[p].resize(local_mesh[p].elementCount);
    }
    global_element_map.resize(tot_size);

    {
      size_t error = 0;
      size_t offset = 0;
      for (int p = 0; p < part_count; p++) {
	ExodusFile id(p);
	error += ex_get_elem_num_map(id, &global_element_numbers[p][0]);
	std::copy(global_element_numbers[p].begin(), global_element_numbers[p].end(),
		  &global_element_map[offset]);
	offset += local_mesh[p].elementCount;
      }
    }

    // Now, sort the global_element_map array.
    std::sort(global_element_map.begin(), global_element_map.end());

    // See if any duplicates...
    GlobalMap temp(global_element_map);
    temp.erase(std::unique(temp.begin(), temp.end()), temp.end());

    if (temp.size() != global_element_map.size()) {
      // Duplicates in the element id list...
      // This is not yet handled.  Notify the user and continue for now...
      std::cerr << "\n!!!! POSSIBLE ERROR: There were " << global_element_map.size() - temp.size()
		<< " elements with duplicated ids detected.\n"
		<< "!!!!\tThis may cause problems in the output file.\n\n";
    }

    size_t total_num_elements = global_element_map.size();
    global->elementCount = total_num_elements;

    // See whether the element numbers are contiguous.  If so, we can map
    // the elements back to their original location. Since the elements are
    // sorted and there are no duplicates, we just need to see if the id
    // at global_element_map.size() == global_element_map.size();
    bool is_contiguous = global_element_map[global_element_map.size()-1] == (int)global_element_map.size();
    std::cerr  << "Element id map " << (is_contiguous ? "is" : "is not") << " contiguous.\n";

    // Create the map that maps from a local processor element to the
    // global map. This combines the mapping local processor element to
    // 'global id' and then 'global id' to global position. The
    // mapping is now a direct lookup instead of a lookup followed by
    // a reverse map.
    //
    // If the map is contiguous, then the global_id to global_position map is 1->1
  REMAP:
    if (is_contiguous && map_ids) {
      GMapIter cur_pos = global_element_map.begin();
      size_t element_value;
      for (int p = 0; p < part_count; p++) {
	size_t element_count = local_mesh[p].elementCount;
	for (size_t i = 0; i < element_count; i++) {
	  int global_element = global_element_numbers[p][i];

	  if (cur_pos == global_element_map.end() || *cur_pos != global_element) {
	    std::pair<GMapIter, GMapIter> iter = std::equal_range(global_element_map.begin(),
								  global_element_map.end(),
								  global_element);
	    SMART_ASSERT(iter.first  != iter.second);
	    element_value = iter.first - global_element_map.begin();
	    cur_pos = iter.first;
	  }
	  element_value = cur_pos - global_element_map.begin();
	  local_element_to_global[p][i] = element_value;
	  ++cur_pos;
	}				
      }
    } else {
      IntVector proc_off(part_count);
      std::fill(proc_off.begin(), proc_off.end(), 0);

      size_t gpos = 0;
      for (size_t b = 0; b < glob_blocks.size(); b++) {
	for (int p = 0; p < part_count; p++) {
	  size_t poff = proc_off[p];
	  size_t element_count = blocks[p][b].entity_count();
	  for (size_t e = 0; e < element_count; e++) {
	    local_element_to_global[p][e+poff] = gpos++;
	  }
	  proc_off[p] += element_count;
	}
      }

      IntVector element_map;
      for (int p = 0; p < part_count; p++) {
	size_t element_count = local_mesh[p].elementCount;
	element_map.resize(element_count);
	ex_get_elem_num_map(ExodusFile(p), &element_map[0]);

	for (size_t e = 0; e < element_count; e++) {
	  gpos = local_element_to_global[p][e];
	  global_element_map[gpos] = element_map[e];
	}
      }
    }

    // Now, need to set up the element block offsets.  Within an element
    // block, the ids are contiguous, so to map a global element to its
    // position within the element block, use the formula:
    //
    //    pos_in_block = global_pos - block_offset
    //
    // Since there is the possibility that the element blocks in this
    // file are in a different order than the element blocks in original
    // file, we need to find out the minimum element id in each element
    // block and that is the offset (0-based ids).  Will also get the
    // maximum id to use as a sanity check.
    // Note that ids are contiguous from 0..element_count-1

    for (size_t b = 0; b < glob_blocks.size(); b++) {
      size_t min_id = global_element_map.size();
      size_t max_id = 0;

      for (int p = 0; p < part_count; p++) {
	size_t offset = blocks[p][b].offset_;
	size_t element_count = blocks[p][b].entity_count();
	for (size_t e = 0; e < element_count; e++) {
	  size_t id = local_element_to_global[p][e+offset];
	  min_id = (id < min_id) ? id : min_id;
	  max_id = (id > max_id) ? id : max_id;
	}
      }
      if (glob_blocks[b].entity_count() == 0) {
	min_id = 0;
	max_id = 0;
      } else {
	if (max_id - min_id + 1 != glob_blocks[b].entity_count()) {
	  if (map_ids) {
	    map_ids = false;
	    std::cerr << "WARNING: The element ids are globally contiguous,\n"
		      << "\tbut they are not consistent for element block "
		      << glob_blocks[b].id
		      << ".\n\tRetrying with element id mapping turned off."
		      << std::endl;
	    goto REMAP;
	  } else {
	    std::cerr << "ERROR: The element ids for element block "
		      << glob_blocks[b].id << " are not consistent."
		      << std::endl;
	    std::cerr << "Block " << b << ", Id = " << glob_blocks[b].id
		      << " min/max id = " << min_id+1 << "/" << max_id+1
		      << " size = " << glob_blocks[b].entity_count() << "\n";
	    exit(EXIT_FAILURE);
	  }
	}
      }
      glob_blocks[b].offset_ = min_id;
      if (debug_level & 4) {
	std::cerr << "Block " << b << ", Id = " << glob_blocks[b].id
		  << " min/max id = " << min_id+1 << "/" << max_id+1
		  << " offset = " << glob_blocks[b].offset_ << "\n";
      }
    }
  }

  void build_reverse_node_map(int **local_node_to_global,
			      const std::vector<Mesh> &local_mesh,
			      Mesh *global, int part_count,
			      GlobalMap &global_node_map)
  {
    // Instead of using <set> and <map>, consider using a sorted vector...
    // Append all local node maps to the global node map.
    // Sort the global node map
    // Remove duplicates.
    // Position withing map is now the map...
    // When building the local-proc node to global id, use binary_search...

    // Global node map and count.
    std::vector<std::vector<int> > global_node_numbers(part_count);

    size_t tot_size = 0;
    for (int p = 0; p < part_count; p++) {
      tot_size += local_mesh[p].nodeCount;
      global_node_numbers[p].resize(local_mesh[p].nodeCount);
    }
    global_node_map.resize(tot_size);

    size_t offset = 0;
    size_t error = 0;
    for (int p = 0; p < part_count; p++) {
      //global_node_numbers[p]    = new int[local_mesh[p].nodeCount];
      ExodusFile id(p);
      error += ex_get_node_num_map(id, &global_node_numbers[p][0]);
      std::copy(global_node_numbers[p].begin(), global_node_numbers[p].end(),
		&global_node_map[offset]);
      offset += local_mesh[p].nodeCount;
    }

    // Now, sort the global_node_map array and remove duplicates...
    std::sort(global_node_map.begin(), global_node_map.end());
    global_node_map.erase(std::unique(global_node_map.begin(), global_node_map.end()), global_node_map.end());

    size_t total_num_nodes = global_node_map.size();
    global->nodeCount = total_num_nodes;

    // See whether the node numbers are contiguous.  If so, we can map
    // the nodes back to their original location. Since the nodes are
    // sorted and there are no duplicates, we just need to see if the id
    // at global_node_map.size() == global_node_map.size();
    bool is_contiguous = global_node_map[global_node_map.size()-1] == (int)global_node_map.size();
    std::cerr  << "Node map " << (is_contiguous ? "is" : "is not") << " contiguous.\n";

    // Create the map the maps from a local processor node to the
    // global map. This combines the mapping local processor node to
    // 'global id' and then 'global id' to global position. The
    // mapping is now a direct lookup instead of a lookup followed by
    // a reverse map.
    GMapIter cur_pos = global_node_map.begin();
    size_t nodal_value;
    for (int p = 0; p < part_count; p++) {
      size_t node_count = local_mesh[p].nodeCount;
      for (size_t i = 0; i < node_count; i++) {
	int global_node = global_node_numbers[p][i];

	if (cur_pos == global_node_map.end() || *cur_pos != global_node) {
	  std::pair<GMapIter, GMapIter> iter = std::equal_range(global_node_map.begin(),
								global_node_map.end(),
								global_node);
	  SMART_ASSERT(iter.first  != iter.second);
	  nodal_value = iter.first - global_node_map.begin();
	  cur_pos = iter.first;
	}
	nodal_value = cur_pos - global_node_map.begin();
	local_node_to_global[p][i] = nodal_value;
	++cur_pos;
      }				
    }
  }

  void get_put_variable_names(int id, int out, Variables &vars,
			      Excn::SystemInterface &interface)
  {
    if (vars.count(OUT) > 0) {

      char **output_name_list  = get_name_array(vars.count(OUT), ExodusFile::max_name_length());

      int extra = vars.add_processor_id() ? 1 : 0;
      int num_input_vars = vars.index_.size();

      char **input_name_list = get_name_array(num_input_vars, ExodusFile::max_name_length());
      ex_get_variable_names(id, vars.type(), num_input_vars-extra, input_name_list);

      if (vars.add_processor_id()) {
	// Check that input list of names does not already contain 'processor_id'...
	// If found, create an 'epu'-specific name; don't redo the check; assume ok...
	bool found = false;
	for (int i=0; i < num_input_vars && !found; i++) {
	  if (case_compare(input_name_list[i], "processor_id") == 0) {
	    found = true;
	  }
	}
	if (found) {
	  std::cerr << "\nWARNING: Variable 'processor_id' already exists on database.\n"
		    << "         Adding 'processor_id_epu' instead.\n\n";
	  strcpy(input_name_list[num_input_vars-1], "processor_id_epu");
	} else {
	  strcpy(input_name_list[num_input_vars-1], "processor_id");
	}
      }

      // Iterate through the 'var_index' and transfer
      // Assume that the number of pointers is limited to
      // the number of results variables, plus one
      // extra pointer for optional add_processor_id
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
      int width = interface.screen_width();
      // Assume 8 characters for initial tab...
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

      if (!interface.append()) {
	ex_put_variable_names(out, vars.type(), vars.count(OUT), output_name_list);
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

    // If 'type' is ELEMENT and addProcessorId is true, then reserve
    // space for an additional variable.
    int extra = 0;
    if (vars.type() == EX_ELEM_BLOCK && vars.addProcessorId) extra = 1;

    int num_vars;
    ex_get_variable_param (id, vars.type(),  &num_vars);

    vars.index_.resize(num_vars + extra);

    // Create initial index which defaults to no output...
    std::fill(vars.index_.begin(), vars.index_.end(), 0);
    // ...Except for the processor_id (if any)
    if (vars.addProcessorId)
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
    // list so no variables will be written.  Just return 0 or 1 based
    // on the 'add_processor_id' setting. Default var_index is ok.
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

      if (vars.addProcessorId) {
	vars.index_[num_vars] = nz_count; // Already counted above...
      }
      vars.outputCount = nz_count;
      return;
    }
  }


  void put_global_info(const Mesh& global)
  {
    // Write out Global info

    std::cout << " Title: " << global.title.c_str() << "\n\n";
    std::cout << " Number of coordinates per node       =" << std::setw(9)
	      << global.dimensionality << "\n";
    std::cout << " Number of nodes                      =" << std::setw(9)
	      << global.nodeCount << "\n";
    std::cout << " Number of elements                   =" << std::setw(9)
	      << global.elementCount << "\n";
    std::cout << " Number of element blocks             =" << std::setw(9)
	      << global.count(EBLK) << "\n\n";
    std::cout << " Number of nodal point sets           =" << std::setw(9)
	      << global.count(NSET) << "\n";
    std::cout << " Number of element side sets          =" << std::setw(9)
	      << global.count(SSET) << "\n\n";

    int id_out = ExodusFile::output();
    get_put_qa(ExodusFile(0), id_out);
    get_put_coordinate_frames(ExodusFile(0), id_out);
  }

  void get_nodesets(int part_count, int total_node_count,
		    int **local_node_to_global,
		    std::vector<std::vector<NodeSet> > &nodesets,
		    std::vector<NodeSet> &glob_sets)
  {
    // Find number of nodesets in the global model...
    std::set<int> set_ids;
    IntVector ids;

    int bad_ns = 0;
    {
      int   ns_count;
      for (int p=0; p < part_count; p++) {
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
    {
      for (int p=0; p < part_count; p++) {
	ExodusFile id(p);

	nodesets[p].resize(set_ids.size());
	std::set<int>::const_iterator I  = set_ids.begin();
	std::set<int>::const_iterator IE = set_ids.end();

	// Get the ids again so we can map current order back to file order...
	ex_get_ids(id, EX_NODE_SET, &ids[0]);

	int iset = 0;
	while (I != IE) {

	  nodesets[p][iset].id = *I++;

	  // Get the parameters for this nodeset...
	  ex_get_set_param(id, EX_NODE_SET,
			   nodesets[p][iset].id, &nodesets[p][iset].nodeCount,
			   &nodesets[p][iset].dfCount);

	  std::vector<char> name(Excn::ExodusFile::max_name_length()+1);
	  ex_get_name(id, EX_NODE_SET, nodesets[p][iset].id, &name[0]);
	  if (name[0] != '\0') nodesets[p][iset].name_ = &name[0];

	  nodesets[p][iset].position_ = -1;
	  for (size_t j=0; j < ids.size(); j++) {
	    if (nodesets[p][iset].id == ids[j]) {
	      nodesets[p][iset].position_ = j;
	      break;
	    }
	  }

	  if (debug_level & 32) {
	    std::cerr << "Processor " << p << " ";
	    nodesets[p][iset].dump();
	  }
	  ++iset;
	}
      }
    }

    glob_sets.resize(set_ids.size());
    {
      // Now get the nodeset nodes and df.
      // Currently ignore the DF.  Could add if people need it...

      // This is inefficient since the processor loop is on
      // the inside...  The other ordering would use more memory...

      IntVector ns_nodes;
      for (size_t ns = 0; ns < set_ids.size(); ns++) {

	IntVector glob_ns_nodes(total_node_count+1);
	std::fill(glob_ns_nodes.begin(), glob_ns_nodes.end(), 0);

	for (int p=0; p < part_count; p++) {
	  ExodusFile id(p);
	  if (p == 0) {
	    glob_sets[ns].name_ = nodesets[p][ns].name_;
	    glob_sets[ns].id    = nodesets[p][ns].id;
	  }
	
	  if (nodesets[p][ns].position_ != -1) {
	    if (glob_sets[ns].position_ != -1) {
	      SMART_ASSERT(glob_sets[ns].position_ == nodesets[p][ns].position_);
	    } else {
	      glob_sets[ns].position_ = nodesets[p][ns].position_;
	    }
	  }

	  int size = nodesets[p][ns].entity_count();
	  if (size > 0) {
	    ns_nodes.resize(size);
	    ex_get_set(id, EX_NODE_SET, nodesets[p][ns].id, &ns_nodes[0], 0);

	    // The node ids are in local space -- map to global
	    for (int iset=0; iset < size; iset++) {
	      int global_node = local_node_to_global[p][ns_nodes[iset]-1] + 1;
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
	  for (int p=0; p < part_count; p++) {
	    ExodusFile id(p);
	    // Get the nodelist, but store it in nodeOrderMap.
	    // global_pos = nodeOrderMap[i]
	    NodeSet &nset = nodesets[p][ns];
	    int nnodes = nset.entity_count();
	    nset.nodeOrderMap.resize(nnodes);
	    ex_get_set(id, EX_NODE_SET, nset.id, &nset.nodeOrderMap[0], NULL);
	    
	    for (int i=0; i < nnodes; i++) {
	      int local_node = nset.nodeOrderMap[i];                       // 1-based
	      int global_node = local_node_to_global[p][local_node-1] + 1; // 1-based
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


  void get_sideset_metadata(int part_count,
			    std::vector< std::vector<SideSet> > &sets,
			    std::vector<SideSet> &glob_ssets)
  {
    // Find number of sidesets in the global model...
    std::set<int> set_ids;
    IntVector ids;

    int bad_ss = 0;
    {
      for (int p=0; p < part_count; p++) {
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

    {
      for (int p=0; p < part_count; p++) {
	ExodusFile id(p);

	sets[p].resize(set_ids.size());
	std::set<int>::const_iterator I  = set_ids.begin();
	std::set<int>::const_iterator IE = set_ids.end();

	// Get the ids again so we can map current order back to file order...
	ex_get_ids(id, EX_SIDE_SET, &ids[0]);

	int i = 0;
	while (I != IE) {
	  sets[p][i].id = *I++;

	  // Get the parameters for this sideset...
	  ex_get_set_param (id, EX_SIDE_SET, sets[p][i].id,
			    &sets[p][i].sideCount,
			    &sets[p][i].dfCount);

	  glob_ssets[i].id         = sets[p][i].id;
	  glob_ssets[i].sideCount += sets[p][i].entity_count();
	  glob_ssets[i].dfCount   += sets[p][i].dfCount;

	  std::vector<char> name(Excn::ExodusFile::max_name_length()+1);
	  ex_get_name(id, EX_SIDE_SET, sets[p][i].id, &name[0]);
	  if (name[0] != '\0') {
	    sets[p][i].name_ = &name[0];
	    if (p == 0)
	      glob_ssets[i].name_ = &name[0];
	  }

	  sets[p][i].position_ = -1;
	  for (size_t j=0; j < ids.size(); j++) {
	    if (sets[p][i].id == ids[j]) {
	      sets[p][i].position_ = j;
	      break;
	    }
	  }

	  if (sets[p][i].position_ != -1) {
	    if (glob_ssets[i].position_ != -1) {
	      SMART_ASSERT(glob_ssets[i].position_ == sets[p][i].position_);
	    } else {
	      glob_ssets[i].position_ = sets[p][i].position_;
	    }
	  }

	  i++;
	}
      }

      // Calculate sideset offset
      for (size_t b=0; b < glob_ssets.size(); b++) {
	int sum = 0;
	for (int p=0; p < part_count; p++) {
	  sets[p][b].offset_ = sum;
	  sum += sets[p][b].entity_count();

	  if (debug_level & 16) {
	    std::cerr << "Processor " << p << " ";
	    sets[p][b].dump();
	  }
	}
      }
    }
  }

  void get_put_sidesets(int part_count, int **local_element_to_global,
			std::vector< std::vector<SideSet> > &sets,
			std::vector<SideSet> &glob_ssets,
			Excn::SystemInterface &interface)
  {
    // TODO: See what work is really needed if in append mode...
    
    // Get a temporary vector to maintain the current
    // offset into the glob_ssets for storing sides
    IntVector offset(glob_ssets.size());
    std::fill(offset.begin(), offset.end(), 0);

    IntVector df_offset(glob_ssets.size());
    std::fill(df_offset.begin(), df_offset.end(), 0);

    {
      for (size_t ss = 0; ss < glob_ssets.size(); ss++) {
	glob_ssets[ss].elems.resize(glob_ssets[ss].entity_count());
	glob_ssets[ss].sides.resize(glob_ssets[ss].entity_count());
	glob_ssets[ss].distFactors.resize(glob_ssets[ss].dfCount * ExodusFile::io_word_size());
      }
    }

    // Now get the sideset elements, sides and df.
    for (int p=0; p < part_count; p++) {
      ExodusFile id(p);
      for (size_t ss = 0; ss < glob_ssets.size(); ss++) {

	int size = sets[p][ss].entity_count();
	if (size > 0) {
	  int off = offset[ss];
	  ex_get_set(id, EX_SIDE_SET, sets[p][ss].id,
		     &glob_ssets[ss].elems[off],
		     &glob_ssets[ss].sides[off]);

	  // The element ids are in local space -- map to global
	  for (int i=0; i < size; i++) {
	    int local_elem = glob_ssets[ss].elems[off+i];

	    int global_elem = local_element_to_global[p][local_elem-1];
	    glob_ssets[ss].elems[off+i]=global_elem+1;
	  }
	  offset[ss] += size;
	}

	// Distribution factors...
	int df_size = sets[p][ss].dfCount;
	if (df_size > 0) {
	  int df_off = df_offset[ss] * ExodusFile::io_word_size();
	  ex_get_set_dist_fact(id, EX_SIDE_SET, sets[p][ss].id,
			       &glob_ssets[ss].distFactors[df_off]);
	  df_offset[ss] += df_size;
	}
      }
    }

    if (debug_level & 16)
      std::cerr << "\nOutput SideSets:\n";

    if (debug_level & 16) {
      for (size_t ss = 0; ss < glob_ssets.size(); ss++) {
	glob_ssets[ss].dump();
      }
    }

    if (!interface.append()) {
      // Now write the actual sideset data...
      int exoid = ExodusFile::output();// output file identifier
      for (size_t ss = 0; ss < glob_ssets.size(); ss++) {
	ex_put_set(exoid, EX_SIDE_SET, glob_ssets[ss].id, const_cast<int*>(&glob_ssets[ss].elems[0]),
		   const_cast<int*>(&glob_ssets[ss].sides[0]));
	if (glob_ssets[ss].dfCount > 0) {
	  ex_put_set_dist_fact(exoid, EX_SIDE_SET, glob_ssets[ss].id,
			       reinterpret_cast<void*>(&glob_ssets[ss].distFactors[0]));
	}
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
  void add_processor_variable(int id_out, int part_count, int start_part,
			      const Mesh& global,
			      std::vector<std::vector<Block> > &blocks,
			      const std::vector<Block> &glob_blocks,
			      int **local_element_to_global,
			      int step, int variable, std::vector<T> &proc)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    for (size_t b = 0; b < global.count(EBLK); b++) {
      proc.resize(glob_blocks[b].entity_count());
      for (int p = 0; p < part_count; p++) {
	int boffset = blocks[p][b].offset_;
	int goffset = glob_blocks[b].offset_;
	int element_count = blocks[p][b].entity_count();
	for (int e = 0; e < element_count; e++) {
	  int global_block_pos = local_element_to_global[p][(e + boffset)] - goffset;
	  proc[global_block_pos] = p+start_part;
	}
      }
      ex_put_var(id_out, step, EX_ELEM_BLOCK, variable,
		 glob_blocks[b].id, glob_blocks[b].entity_count(),
		 &proc[0]);
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
	  SMART_ASSERT(inp_position >= 0);
	  SMART_ASSERT(out_position >= 0);

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
    // read truth table - sum across all processors since many will
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
      int part_count = local.size();
      for (int p = 0; p < part_count; p++) {
	ExodusFile id(p);

	if (vars.count(IN) > 0) { // Could be zero if add_processor_id
	  // is the only variable...
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
	  if (vars.addProcessorId) {
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
	std::cerr << "Truth table for " << vars.label() << "\n";
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

    std::string info = "EPU: ";
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
    return (strchr(white_space, c) != NULL);
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
    int width = 1;
    if (max_value >= 10)
      width = int(log10((double)max_value));
    return width + 1;
  }

  template <typename T, typename U>
  void clear_master_values(Excn::Variables &vars, const Excn::Mesh& global,
			   std::vector<U> &glob_sets, T*** master_values)
  {
    for (int i = 0; i < vars.count(IN); i++) {
      if (vars.index_[i] > 0) {
	int ivar = vars.index_[i]-1;
	SMART_ASSERT(ivar < vars.count(OUT));
	// zero out master array
	for (size_t b = 0; b < global.count(vars.objectType); b++) {
	  int output_truth_table_loc = (b * vars.count(OUT)) + ivar;
	  if (global.truthTable[vars.objectType][output_truth_table_loc]) {
	    std::fill(&master_values[ivar][b][0],
		      &master_values[ivar][b][glob_sets[b].entity_count()], 0.0);
	  }
	}
      }
    }
  }

  template <typename T>
  void map_element_vars(int loffset, int goffset, int entity_count, 
			std::vector<T> &values, T *global_values,
			int *proc_loc_elem_to_global)
  {
    // copy values to master element value information
    T* local_values = &values[0];
    for (int j = 0; j < entity_count; j++) {
      int global_block_pos = proc_loc_elem_to_global[(j + loffset)] - goffset;
      global_values[global_block_pos] = local_values[j];
    }
  }

  template <typename T>
  void map_sideset_vars(int loffset, int entity_count, 
			std::vector<T> &values, T *global_values)
  {
    // copy values to master sideset value information
    T* local_values = &values[0];
    for (int j = 0; j < entity_count; j++) {
      global_values[j+loffset] = local_values[j];
    }
  }

  template <typename T, typename U>
  void map_nodeset_vars(U &local_set, int entity_count, int glob_entity_count,
			std::vector<T> &values, T *global_values)
  {
    SMART_ASSERT(1==0 && "Internal Error!");
  }

  void map_nodeset_vars(Excn::NodeSet &local_set, int entity_count, int glob_entity_count,
			std::vector<double> &values, double *global_values)
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
			std::vector<float> &values, float *global_values)
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
  int read_master_values(Excn::Variables &vars, const Excn::Mesh& global,
			 std::vector<U> &global_sets,
			 std::vector<Excn::Mesh> &local_mesh,
			 std::vector<std::vector<U> > &local_sets,
			 T*** master_values, std::vector<T> &values, 
			 int part_count, int time_step,
			 int **local_element_to_global)
  {
    int error = 0;
    bool is_sidenodeset = vars.objectType==NSET || vars.objectType==SSET;

    for (int p = 0; p < part_count; p++) {
      ExodusFile id(p);

      // Only needed for element, but haven't cleaned this up yet...
      int* proc_loc_elem_to_global = local_element_to_global[p];
	
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
		local_sets[p][b].entity_count() > 0) {

	      T* iv_block_mev = master_values[ivar][b];
	      int entity_count = local_sets[p][b].entity_count();

	      if (local_mesh[p].truthTable[vars.objectType][input_truth_table_loc] > 0) {
		error += ex_get_var(id, time_step+1, exodus_object_type(vars.objectType),
				    i+1, local_sets[p][b].id, entity_count, &values[0]);
		    
		switch (vars.objectType) {
		case EBLK:
		    map_element_vars(local_sets[p][b].offset_, global_sets[b].offset_,
				     entity_count, values, iv_block_mev, proc_loc_elem_to_global);
		  break;
		  
		case SSET:
		    map_sideset_vars(local_sets[p][b].offset_, entity_count, values, iv_block_mev);
		  break;
		  
		case NSET:
		    map_nodeset_vars(local_sets[p][b], entity_count, global_sets[b].entity_count(),
				     values, iv_block_mev);
		  break;
		default:
		  break;
		}
	      }
	    }
	  }
	}
      }
    }
    return error;
  }

  template <typename T, typename U>
  void output_master_values(Excn::Variables &vars, const Excn::Mesh& global,
			    std::vector<U> &glob_sets,
			    T*** master_values, int time_step)
  {
    int id_out = ExodusFile::output();// output file identifier
    for (int i = 0; i < vars.count(IN); i++) {
      if (vars.index_[i] > 0) {
	int ivar = vars.index_[i]-1;
	for (size_t b = 0; b < global.count(vars.objectType); b++) {
	  int truth_table_loc = (b * vars.count(OUT)) + ivar;
	
	  if (global.truthTable[vars.objectType][truth_table_loc]) {
	    ex_put_var(id_out, time_step, exodus_object_type(vars.objectType),
		       ivar+1, glob_sets[b].id, glob_sets[b].entity_count(),
		       master_values[ivar][b]);
	  }
	}
      }
    }
  }
  
  template <typename T, typename U>
  void allocate_master_values(Excn::Variables &vars, const Excn::Mesh& global,
			      std::vector<U> &glob_sets, T*** &master_values)
  {
    master_values = new T**[vars.count(OUT)];
    for (int i = 0; i < vars.count(IN); i++) {
      if (vars.index_[i] > 0) {
	int ivar = vars.index_[i]-1;
	master_values[ivar] = new T*[global.count(vars.objectType)];
	for (size_t b = 0; b < global.count(vars.objectType); b++) {
	  int output_truth_table_loc = (b * vars.count(OUT)) + ivar;
	  if (global.truthTable[vars.objectType][output_truth_table_loc] &&
	      glob_sets[b].entity_count() > 0) {
	    master_values[ivar][b] = new T[glob_sets[b].entity_count()];
	  } else {
	    master_values[ivar][b] = 0;
	  }
	}
      }
    }
  }

  template <typename T>
  void deallocate_master_values(Excn::Variables &vars, const Excn::Mesh& global,
				T*** &master_values)
  {
    for (int i=0; i < vars.count(IN); i++) {
      if (vars.index_[i] > 0) {
	int ivar = vars.index_[i]-1;
	for (size_t b = 0; b < global.count(vars.objectType); b++) {
	  delete []  master_values[ivar][b];
	}
	delete [] master_values[ivar];
      }
    }
    delete [] master_values;
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

  size_t find_max_entity_count(int part_count,
			       std::vector<Excn::Mesh> &local_mesh,
			       const Excn::Mesh& global,
			       std::vector<std::vector<Block> > &blocks,
			       std::vector<std::vector<NodeSet> > &nodesets,
			       std::vector<std::vector<SideSet> > &sidesets)
  {
    size_t max_ent = local_mesh[0].nodeCount;
    for (int p = 1; p < part_count; p++) {
      if ((size_t)local_mesh[p].nodeCount > max_ent)
	max_ent = local_mesh[p].nodeCount;
    }
    
    for (int p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(EBLK); b++) {
	if (blocks[p][b].entity_count() > max_ent)
	  max_ent = blocks[p][b].entity_count();
      }
    }
    
    // Nodesets...
    for (int p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(NSET); b++) {
	if (nodesets[p][b].entity_count() > max_ent)
	  max_ent = nodesets[p][b].entity_count();
      }
    }
    
    // Sidesets...
    for (int p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(SSET); b++) {
	if (sidesets[p][b].entity_count() > max_ent)
	  max_ent = sidesets[p][b].entity_count();
      }
    }
    return max_ent;
  }
}
