// ISSUES:
// 1. Does not handle unconnected nodes (not connected to any element)
//
// 2. Sideset distribution factors are klugy and may not fully work in
//    strange cases
//
//
// Copyright(C) 2012
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
#include <exodus/Ioex_Utils.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_SurfaceSplit.h>
#include <Ioss_FileInfo.h>
#include <Ioss_Utils.h>
#include <assert.h>
#include <exodusII_par.h>
#include <par_exo/Iopx_DatabaseIO.h>
#include <par_exo/Iopx_Internals.h>
#include <par_exo/Iopx_DecompositionData.h>
#include <float.h>
#include <stddef.h>
#include <sys/select.h>
#include <time.h>
#include <tokenize.h>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <numeric>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Ioss_CoordinateFrame.h"
#include "Ioss_CommSet.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_EntityBlock.h"
#include "Ioss_EntitySet.h"
#include "Ioss_EntityType.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_Map.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_State.h"
#include "Ioss_VariableType.h"

#include "Ioss_FileInfo.h"
#undef MPICPP

// ========================================================================
// Static internal helper functions
// ========================================================================
namespace {
  inline size_t min(size_t x, size_t y)
  {
    return y ^ ((x^y) & -(x<y));
  }

  inline size_t max(size_t x, size_t y)
  {
    return y ^ ((x^y) & -(x>y));
  }

  const size_t max_line_length   = MAX_LINE_LENGTH;

  const std::string SEP() {return std::string("@");} // Separator for attribute offset storage
  const std::string SCALAR()     {return std::string("scalar");}
  const std::string VECTOR3D()   {return std::string("vector_3d");}
  const std::string SYM_TENSOR() {return std::string("sym_tensor_33");}

  const char *complex_suffix[] = {".re", ".im"};

  int get_parallel_io_mode(const Ioss::PropertyManager &properties)
  {
    static int par_mode = 0;
    static int par_mode_default = EX_PNETCDF;  // Default...
    //    static int par_mode_default = EX_MPIIO; 
    //    static int par_mode_default = EX_MPIPOSIX;

    if (par_mode == 0) {
      if (properties.exists("PARALLEL_IO_MODE")) {
        std::string mode_name = properties.get("PARALLEL_IO_MODE").get_string();
        Ioss::Utils::fixup_name(mode_name);
        if (mode_name == "pnetcdf")
          par_mode = EX_PNETCDF;
        else if (mode_name == "mpiposix")
          par_mode = EX_MPIPOSIX;
        else if (mode_name == "mpiio")
          par_mode = EX_MPIIO;
        else {
          std::ostringstream errmsg;
          errmsg << "ERROR: Unrecognized parallel io mode setting '" << mode_name
              << "' in the Ioss property PARALLEL_IO_MODE.  Valid values are 'pnetcdf', 'mpiposix', or 'mpiio'";
          IOSS_ERROR(errmsg);
        }
      }
      else {
        par_mode = par_mode_default;
      }
    }
    return par_mode;
  }

  void get_connectivity_data(int exoid, void *data, ex_entity_type type, ex_entity_id id, int position)
  {
    int ierr = 0;
    if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
      int64_t* conn[3];
      conn[0] = NULL;
      conn[1] = NULL;
      conn[2] = NULL;
      conn[position] = static_cast<int64_t*>(data);
      assert(1==0 && "Unimplemented fixme");
      ierr = ex_get_conn(exoid, type, id, conn[0], conn[1], conn[2]); // FIXME
    } else {
      int* conn[3];
      conn[0] = NULL;
      conn[1] = NULL;
      conn[2] = NULL;
      conn[position] = static_cast<int*>(data);
      assert(1==0 && "Unimplemented fixme");
      ierr = ex_get_conn(exoid, type, id, conn[0], conn[1], conn[2]); // FIXME
    }
    if (ierr < 0)
      Ioex::exodus_error(exoid, __LINE__, -1);
  }

  void check_variable_consistency(const ex_var_params &exo_params,
                                  int my_processor, const std::string &filename,
                                  const Ioss::ParallelUtils &util);

  template <typename T>
  void compute_internal_border_maps(T* entities, T* internal, size_t count, size_t entity_count)
  {
    for (size_t ij=0; ij < count; ij++)
      internal[ij] = 1;
    for (size_t J=0; J < entity_count; J++)
      internal[entities[J]-1] = 0;

    size_t b = 0;
    for (size_t ij=0; ij < count; ij++) {
      if (internal[ij] == 0) {
        entities[b++] = ij+1;
      }
    }

    size_t k = 0;
    for (size_t ij=0; ij < count; ij++) {
      if (internal[ij] == 1) {
        internal[k++] = ij+1;
      }
    }
  }

  template <typename INT>
  void map_nodeset_id_data(const Ioss::IntVector &owning_processor, Ioss::Int64Vector &owned_nodes,
                           int this_processor,  const INT *ids, size_t ids_size, std::vector<INT> &file_data)
  {
    // Determine which nodes in this nodeset are owned by this processor.
    // Save this mapping in the "owned_nodes" vector for use in
    // mapping nodeset field data (df, transient, attributes, ...) 
    for (size_t i=0; i < ids_size; i++) {
      INT node = ids[i];
      if (owning_processor[node-1] == this_processor) {
        file_data.push_back(ids[i]);
        owned_nodes.push_back(i);
      }
    }
  }

  template <typename T>
  void map_nodeset_data(Ioss::Int64Vector &owned_nodes, const T *data, std::vector<T> &file_data,
                        size_t offset=0, size_t stride=1)
  {
    // Pull out the locally owned nodeset data 
    for (size_t i=0; i < owned_nodes.size(); i++) {
      file_data.push_back(data[stride*owned_nodes[i]+offset]);
    }
  }

  // Ideally, there should only be a single data type for in and out
  // data, but in the node id map mapping, we have an int64_t coming
  // in and either an int or int64_t going out...
  template <typename T, typename U>
  void map_data(const Ioss::IntVector &owning_processor, int this_processor, 
                const T *data, std::vector<U> &file_data, size_t offset=0, size_t stride=1)
  {
    size_t count = owning_processor.size();
    size_t index = offset;
    for (size_t i=0; i < count; i++) {
      if (owning_processor[i] == this_processor) {
        file_data.push_back(data[index]);
      }
      index += stride;
    }
  }

  template <typename INT>
  void map_local_to_global_implicit(INT *data, size_t count, const std::vector<int64_t> &global_implicit_map)
  {
    for (size_t i=0; i < count; i++) {
      data[i] = global_implicit_map[data[i]-1];
    }
  }
}

namespace Iopx {
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
                         Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                         const Ioss::PropertyManager &props) :
    Ioex::DatabaseIO(region, filename, db_usage, communicator, props),
    decomp(NULL), decomp32(NULL), decomp64(NULL),
    metaDataWritten(false)
  {
    if (!is_parallel_consistent()) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Parallel IO cannot be used in an application that is not guaranteeing "
	     << "parallel consistent calls of the get and put field data functions.\n"
	     << "The application created this database with a 'false' setting for the isParallelConsistent property.";
      IOSS_ERROR(errmsg);
    }
  }

  DatabaseIO::~DatabaseIO()
  {
    try {
      delete decomp;
    } catch (...) {
    }
  }

  bool DatabaseIO::ok(bool write_message, std::string *error_msg, int *bad_count) const
  {
    if (fileExists) {
      // File has already been opened at least once...
      return dbState != Ioss::STATE_INVALID;
    }

    // File has not yet been opened, so verify that it is readable or
    // writable.  HOWEVER, *DO NOT* overwrite the file here if it exists
    // and the mode is WRITE.
    int cpu_word_size = sizeof(double);
    int io_word_size  = 0;
    float version;

    int par_mode = get_parallel_io_mode(properties);

    MPI_Info info = MPI_INFO_NULL;
    int exodus_file_ptr = ex_open_par(get_filename().c_str(), EX_READ|par_mode,
				      &cpu_word_size, &io_word_size, &version, util().communicator(), info);

    if (!is_input() && exodus_file_ptr < 0) {
      // File didn't exist above, but this OK if is an output file. See if we can create it...
      int mode = 0;
      if (int_byte_size_api() == 8)
        mode |= EX_ALL_INT64_DB;

      exodus_file_ptr = ex_create_par(get_filename().c_str(), exodusMode|mode|par_mode,
				      &cpu_word_size, &dbRealWordSize, util().communicator(), info);
    }

    // Check for valid exodus_file_ptr (valid >= 0; invalid < 0)
    int global_file_ptr = util().global_minmax(exodus_file_ptr, Ioss::ParallelUtils::DO_MIN);
    if (global_file_ptr < 0 && (write_message || error_msg != NULL || bad_count != NULL)) {
      Ioss::IntVector status;
      util().all_gather(exodus_file_ptr, status);

      if (write_message || error_msg != NULL) {
	// See which processors could not open/create the file...
	std::string open_create = is_input() ? "open input" : "create output";
	bool first = true;
	std::ostringstream errmsg;
	errmsg << "ERROR: Unable to " << open_create << " database '" << get_filename() << "' of type 'exodusII'";
	if (isParallel) {
	  errmsg << "\n\ton processor(s): ";
	  for (int i=0; i < util().parallel_size(); i++) {
	    if (status[i] < 0) {
	      if (!first) errmsg << ", ";
	      errmsg << i;
	      first = false;
	    }
	  }
	}
	if (error_msg != NULL) {
	  *error_msg = errmsg.str();
	}
	if (write_message && myProcessor == 0) {
	  errmsg << "\n";
	  std::cerr << errmsg.str();
	}
      }
      if (bad_count != NULL) {
	for (int i=0; i < util().parallel_size(); i++) {
	  if (status[i] < 0) {
	    (*bad_count)++;
	  }
	}
      }
    }

    if (is_input() && exodus_file_ptr >= 0) {
      // Check byte-size of integers stored on the database...
      int mode = ex_int64_status(exodus_file_ptr) & EX_ALL_INT64_DB;
      if (mode) {
        exodusMode |= mode;
        set_int_byte_size_api(Ioss::USE_INT64_API);
      }
    }

    // Close all open files...
    if (exodus_file_ptr >= 0) {
      ex_close(exodus_file_ptr);
    }
    return global_file_ptr >= 0;
  }

  int DatabaseIO::get_file_pointer() const
  {
    // Returns the file_pointer used to access the file on disk.
    // Checks that the file is open and if not, opens it first.
    if (exodusFilePtr < 0) {
      int cpu_word_size = sizeof(double);
      int io_word_size  = 0;
      float version;

      int mode = exodusMode;
      if (int_byte_size_api() == 8) {
        mode |= EX_ALL_INT64_API;
      }

      int par_mode = get_parallel_io_mode(properties);

      MPI_Info info = MPI_INFO_NULL;
      if (is_input()) {
        exodusFilePtr = ex_open_par(get_filename().c_str(), EX_READ|mode|par_mode,
				    &cpu_word_size, &io_word_size, &version, util().communicator(), info);

      } else {
        if (fileExists) {
          exodusFilePtr = ex_open_par(get_filename().c_str(), EX_WRITE|mode|par_mode,
				      &cpu_word_size, &io_word_size, &version, util().communicator(), info);
        } else {
          // If the first write for this file, create it...
          if (int_byte_size_api() == 8)
            mode |= EX_ALL_INT64_DB;
          exodusFilePtr = ex_create_par(get_filename().c_str(), mode|par_mode,
					&cpu_word_size, &dbRealWordSize, util().communicator(), info);
          if (exodusFilePtr < 0) {
            dbState = Ioss::STATE_INVALID;
            // NOTE: Code will not continue past this call...
            std::ostringstream errmsg;
            errmsg << "ERROR: Cannot create specified file '" << get_filename() << "'";
            IOSS_ERROR(errmsg);
          }
        }
        ex_set_max_name_length(exodusFilePtr, maximumNameLength);
      }

      if (exodusFilePtr < 0) {
        dbState = Ioss::STATE_INVALID;
        fileExists = false;
        // NOTE: Code will not continue past this call...
        std::ostringstream errmsg;
        errmsg << "ERROR: Problem opening specified file '" << get_filename() << "'";
        IOSS_ERROR(errmsg);
      }

      if (is_input()) {
        // Check for maximum name length used on the input file.
        int max_name_length = ex_inquire_int(exodusFilePtr, EX_INQ_DB_MAX_USED_NAME_LENGTH);
        if (max_name_length > maximumNameLength) {
          ex_set_max_name_length(exodusFilePtr, max_name_length);
          maximumNameLength = max_name_length;
        }
      }

      // Check properties handled post-create/open...
      if (properties.exists("COMPRESSION_LEVEL")) {
        int comp_level = properties.get("COMPRESSION_LEVEL").get_int();
        ex_set_option(exodusFilePtr, EX_OPT_COMPRESSION_LEVEL, comp_level);
      }
      if (properties.exists("COMPRESSION_SHUFFLE")) {
        int shuffle = properties.get("COMPRESSION_SHUFFLE").get_int();
        ex_set_option(exodusFilePtr, EX_OPT_COMPRESSION_SHUFFLE, shuffle);
      }

    }
    assert(exodusFilePtr >= 0);
    fileExists = true;
    return exodusFilePtr;
  }

  void DatabaseIO::read_meta_data()
  {
    int exoid = get_file_pointer(); // get_file_pointer() must be called first.

    if (int_byte_size_api() == 8) {
      decomp64 = new DecompositionData<int64_t>(properties, util().communicator());
      decomp = decomp64;
    } else {
      decomp32 = new DecompositionData<int>(properties, util().communicator());
      decomp = decomp32;
    }
    assert(decomp != NULL);
    decomp->decompose_model(exoid);

    read_region();
    get_elemblocks();

    get_step_times();

    get_nodeblocks();
    get_edgeblocks();
    get_faceblocks();

    check_side_topology();

    get_nodesets();
    get_sidesets();
#if 0
    get_edgesets();
    get_facesets();
    get_elemsets();
#endif

    get_commsets();

    handle_groups();

    add_region_fields();

    // This closes the file.  It will be automatically opened the next time the file is
    // accessed and it solves some issues with initial condition
    // data...
    free_file_pointer();
  }

  void DatabaseIO::read_region()
  {
    // Add properties and fields to the 'owning' region.
    // Also defines member variables of this class...
    ex_init_params info;
    int error = ex_get_init_ext(get_file_pointer(), &info);
    if (error < 0)
      Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

    spatialDimension = decomp->spatialDimension;
    nodeCount = decomp->ioss_node_count();
    edgeCount = 0;
    faceCount = 0;
    elementCount = decomp->ioss_elem_count();

    m_groupCount[EX_NODE_BLOCK] = 1;
    m_groupCount[EX_EDGE_BLOCK] = info.num_edge_blk;
    m_groupCount[EX_FACE_BLOCK] = info.num_face_blk;
    m_groupCount[EX_ELEM_BLOCK] = info.num_elem_blk;

    m_groupCount[EX_NODE_SET] = info.num_node_sets;
    m_groupCount[EX_EDGE_SET] = info.num_edge_sets;
    m_groupCount[EX_FACE_SET] = info.num_face_sets;
    m_groupCount[EX_ELEM_SET] = info.num_elem_sets;

    m_groupCount[EX_SIDE_SET] = info.num_side_sets;

    if (nodeCount == 0) {
      IOSS_WARNING << "No nodes were found in the model, file '" << get_filename() << "'";
    } else if (nodeCount < 0) {
      // NOTE: Code will not continue past this call...
      std::ostringstream errmsg;
      errmsg << "ERROR: Negative node count was found in the model\n"
	     << "       File: '" << get_filename() << "'.\n";
      IOSS_ERROR(errmsg);
    }

    if (elementCount == 0) {
      IOSS_WARNING << "No elements were found in the model, file: '" << get_filename() << "'";
    }

    if (elementCount < 0) {
      // NOTE: Code will not continue past this call...
      std::ostringstream errmsg;
      errmsg << "ERROR: Negative element count was found in the model, file: '"
	     << get_filename() << "'";
      IOSS_ERROR(errmsg);
    }

    if (elementCount > 0 && m_groupCount[EX_ELEM_BLOCK] <= 0) {
      // NOTE: Code will not continue past this call...
      std::ostringstream errmsg;
      errmsg << "ERROR: No element blocks were found in the model, file: '" << get_filename() << "'";
      IOSS_ERROR(errmsg);
    }

    Ioss::Region *this_region = get_region();

    // See if any coordinate frames exist on mesh.  If so, define them on region.
    Ioex::add_coordinate_frames(get_file_pointer(), this_region);

    this_region->property_add(Ioss::Property("global_node_count",    (int64_t)decomp->globalNodeCount));
    this_region->property_add(Ioss::Property("global_element_count", (int64_t)decomp->globalElementCount));

    this_region->property_add(Ioss::Property(std::string("title"), info.title));
    this_region->property_add(Ioss::Property(std::string("spatial_dimension"),
                                             spatialDimension));

    // Get information records from database and add to informationRecords...
    int num_info = ex_inquire_int(get_file_pointer(), EX_INQ_INFO);    
    if (num_info > 0) {
      char** info_rec = Ioex::get_exodus_names(num_info, max_line_length); // 'total_lines' pointers to char buffers
      ex_get_info(get_file_pointer(), info_rec);
      for (int i=0; i < num_info; i++) {
        add_information_record(info_rec[i]);
      }
      Ioex::delete_exodus_names(info_rec, num_info);
    }
  }

  void DatabaseIO::get_step_times()
  {
    double last_time = DBL_MAX;
    int timestep_count = 0;
    std::vector<double> tsteps(0);

    {
      timestep_count = ex_inquire_int(get_file_pointer(), EX_INQ_TIME);
      if (timestep_count <= 0)
        return;

      // For an exodusII file, timesteps are global and are stored in the region.
      // Read the timesteps and add to the region
      tsteps.resize(timestep_count);
      int error = ex_get_all_times(get_file_pointer(), TOPTR(tsteps));
      if (error < 0)
	Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // See if the "last_written_time" attribute exists and if it
      // does, check that it matches the largest time in 'tsteps'.
      Ioex::read_last_time_attribute(get_file_pointer(), &last_time);
    }

    // Only add states that are less than or equal to the
    // 'last_time' value which is either DBL_MAX or the value of
    // the last time successfully written to the database and
    // flushed to disk.  This is used to avoid corrupt data arising
    // from a job that crashed during the writing of the last step
    // on the database.  Output a warning message if there is
    // potentially corrupt data on the database...

    Ioss::Region *this_region = get_region();
    for (int i=0; i < timestep_count; i++) {
      if (tsteps[i] <= last_time) {
        this_region->add_state(tsteps[i]*timeScaleFactor);
      } else {
        if (myProcessor == 0) {
          // NOTE: Don't want to warn on all processors if there are
          // corrupt steps on all databases, but this will only print
          // a warning if there is a corrupt step on processor
          // 0... Need better warnings which won't overload in the
          // worst case...
          IOSS_WARNING << "Skipping step " << i+1 << " at time " << tsteps[i]
		       << " in database file\n\t" << get_filename()
		       << ".\nThe data for that step is possibly corrupt.\n";
        }
      }
    }
  }

  const Ioss::Map& DatabaseIO::get_map(ex_entity_type type) const
  {
    switch (type) {
    case EX_NODE_BLOCK:
    case EX_NODE_SET:
      {
	size_t offset = decomp ? decomp->nodeOffset : 0;
	size_t count = decomp ? decomp->nodeCount : nodeCount;
	return get_map(nodeMap, nodeCount,
		       offset, count, EX_NODE_MAP, EX_INQ_NODE_MAP);
      }
    case EX_ELEM_BLOCK:
    case EX_ELEM_SET:
      {
	size_t offset = decomp ? decomp->elementOffset : 0;
	size_t count = decomp ? decomp->elementCount : elementCount;
	return get_map(elemMap, elementCount,
		       offset, count, EX_ELEM_MAP, EX_INQ_ELEM_MAP);
      }

    case EX_FACE_BLOCK:
    case EX_FACE_SET:
      return get_map(faceMap, faceCount,
		     0, 0,
		     EX_FACE_MAP, EX_INQ_FACE_MAP);

    case EX_EDGE_BLOCK:
    case EX_EDGE_SET:
      return get_map(edgeMap, edgeCount,
		     0, 0,
		     EX_EDGE_MAP, EX_INQ_EDGE_MAP);

    default:
      std::ostringstream errmsg;
      errmsg << "INTERNAL ERROR: Invalid map type. "
	     << "Something is wrong in the Iopx::DatabaseIO::get_map() function. "
	     << "Please report.\n";
      IOSS_ERROR(errmsg);
    }      
  }

  const Ioss::Map& DatabaseIO::get_map(Ioss::Map &entity_map,
                                       int64_t entityCount,
                                       int64_t file_offset, int64_t file_count,
                                       ex_entity_type entity_type,
                                       ex_inquiry inquiry_type) const
  {
    // Allocate space for node number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (entity_map.map.empty()) {
      entity_map.map.resize(entityCount+1);

      if (is_input()) {
        Ioss::MapContainer file_data(file_count);
        int error = 0;
        // Check whether there is a "original_global_id_map" map on
        // the database. If so, use it instead of the "node_num_map".
        bool map_read = false;
        int map_count = ex_inquire_int(get_file_pointer(), inquiry_type);
        if (map_count > 0) {
          char **names = Ioex::get_exodus_names(map_count, maximumNameLength);
          int ierr = ex_get_names(get_file_pointer(), entity_type, names);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, -1);

          if (map_count == 1 && Ioss::Utils::case_strcmp(names[0], "original_global_id_map") == 0) {
            if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
              error = ex_get_partial_num_map(get_file_pointer(), entity_type, 1,
					     file_offset+1, file_count, TOPTR(file_data));
            } else {
              // Ioss stores as 64-bit, read as 32-bit and copy over...
              Ioss::IntVector tmp_map(file_count);
              error = ex_get_partial_num_map(get_file_pointer(), entity_type, 1,
					     file_offset+1, file_count, TOPTR(tmp_map));
              std::copy(tmp_map.begin(), tmp_map.end(), file_data.begin());
            }
            if (error >= 0) {
              map_read = true;
            }
          }
	  Ioex::delete_exodus_names(names, map_count);
        }

        if (!map_read) {
          if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
            error = ex_get_partial_id_map(get_file_pointer(), entity_type,
					  file_offset+1, file_count, TOPTR(file_data));
          } else {
            // Ioss stores as 64-bit, read as 32-bit and copy over...
            Ioss::IntVector tmp_map(file_count);
            error = ex_get_partial_id_map(get_file_pointer(), entity_type,
					  file_offset+1, file_count, TOPTR(tmp_map));
            std::copy(tmp_map.begin(), tmp_map.end(), file_data.begin());
          }
        }

        if (error >= 0) {
          if (entity_type == EX_NODE_MAP)
            decomp->communicate_node_data(TOPTR(file_data), &entity_map.map[1], 1);
          else if (entity_type == EX_ELEM_MAP)
            decomp->communicate_element_data(TOPTR(file_data), &entity_map.map[1], 1);
        } else {
          // Clear out the vector...
          Ioss::MapContainer().swap(entity_map.map);
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        }

        // Check for sequential node map.
        // If not, build the reverse G2L node map...
        entity_map.map[0] = -1;
        for (int64_t i=1; i < entityCount+1; i++) {
          if (i != entity_map.map[i]) {
            entity_map.map[0] = 1;
            break;
          }
        }

        entity_map.build_reverse_map(myProcessor);

      } else {
        // Output database; entity_map.map not set yet... Build a default map.
        for (int64_t i=1; i < entityCount+1; i++) {
          entity_map.map[i] = i;
        }
        // Sequential map
        entity_map.map[0] = -1;
      }
    }
    return entity_map;
  }

  void DatabaseIO::get_elemblocks()
  {
    get_blocks(EX_ELEM_BLOCK, 0, "block");
  }

  void DatabaseIO::get_faceblocks()
  {
    //    get_blocks(EX_FACE_BLOCK, 1, "faceblock");
  }

  void DatabaseIO::get_edgeblocks()
  {
    //    get_blocks(EX_EDGE_BLOCK, 2, "edgeblock");
  }

  void DatabaseIO::get_blocks(ex_entity_type entity_type, int rank_offset, const std::string &basename)
  {
    // Attributes of an X block are:  (X = element, face, or edge)
    // -- id
    // -- name
    // -- X type
    // -- number of Xs
    // -- number of attributes per X
    // -- number of nodes per X (derivable from type)
    // -- number of faces per X (derivable from type)
    // -- number of edges per X (derivable from type)

    // In a parallel execution, it is possible that an X block will have
    // no Xs on a particular processor...

    // NOTE: This routine may be called multiple times on a single database.
    //       make sure it is not dependent on being called one time only...

    // Get exodusII X block metadata
    if (m_groupCount[entity_type] == 0)
      return;

    assert(entity_type == EX_ELEM_BLOCK);

    Ioss::Int64Vector X_block_ids(m_groupCount[entity_type]);

    int error;
    {
      if (ex_int64_status(get_file_pointer()) & EX_IDS_INT64_API) {
        error = ex_get_ids(get_file_pointer(), entity_type, TOPTR(X_block_ids));
      } else {
        Ioss::IntVector tmp_set_ids(X_block_ids.size());
        error = ex_get_ids(get_file_pointer(), entity_type, TOPTR(tmp_set_ids));
        if (error >= 0) std::copy(tmp_set_ids.begin(), tmp_set_ids.end(), X_block_ids.begin());
      }
      if (error < 0) {
        Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }

      for (int iblk = 0; iblk < m_groupCount[entity_type]; iblk++) {
        int64_t id = decomp->el_blocks[iblk].id();

        std::string alias = Ioss::Utils::encode_entity_name(basename, id);
	bool db_has_name = false;
        std::string block_name = Ioex::get_entity_name(get_file_pointer(), entity_type, id, basename,
						       maximumNameLength, db_has_name);
	if (get_use_generic_canonical_name()) {
	  std::string temp = block_name;
	  block_name = alias;
	  alias = temp;
	}

        std::string save_type = decomp->el_blocks[iblk].topologyType;
        std::string type = Ioss::Utils::fixup_type(decomp->el_blocks[iblk].topologyType,
                                                   decomp->el_blocks[iblk].nodesPerEntity,
                                                   spatialDimension-rank_offset);

        Ioss::EntityBlock *io_block = NULL;
        if (entity_type == EX_ELEM_BLOCK) {
          Ioss::ElementBlock *eblock = new Ioss::ElementBlock(this, block_name, type, decomp->el_blocks[iblk].ioss_count());
          io_block = eblock;
          io_block->property_add(Ioss::Property("id", id));
	  if (db_has_name) {
	    std::string *db_name = &block_name;
	    if (get_use_generic_canonical_name()) {
	      db_name = &alias;
	    }
	    io_block->property_add(Ioss::Property("db_name", *db_name));
	  }
          get_region()->add(eblock);
#if 0
        } else if (entity_type == EX_FACE_BLOCK) {
          Ioss::FaceBlock *fblock = new Ioss::FaceBlock(this, block_name, type, block.num_entry);
          io_block = fblock;
          io_block->property_add(Ioss::Property("id", id));
	  if (db_has_name) {
	    std::string *db_name = &block_name;
	    if (get_use_generic_canonical_name()) {
	      db_name = &alias;
	    }
	    io_block->property_add(Ioss::Property("db_name", *db_name));
	  }
          get_region()->add(fblock);
        } else if (entity_type == EX_EDGE_BLOCK) {
          Ioss::EdgeBlock *eblock = new Ioss::EdgeBlock(this, block_name, type, block.num_entry);
          io_block = eblock;
          io_block->property_add(Ioss::Property("id", id));
	  if (db_has_name) {
	    std::string *db_name = &block_name;
	    if (get_use_generic_canonical_name()) {
	      db_name = &alias;
	    }
	    io_block->property_add(Ioss::Property("db_name", *db_name));
	  }
          get_region()->add(eblock);
#endif
        } else {
          std::ostringstream errmsg;
          errmsg << "ERROR: Invalid type in get_blocks()";
          IOSS_ERROR(errmsg);
        }

#if 0
        // See which connectivity options were defined for this block.
        // X -> Node is always defined.
        // X -> Face?
        if (block.num_faces_per_entry > 0 && rank_offset < 1) {
          std::string storage = "Real["+Ioss::Utils::to_string(block.num_faces_per_entry)+"]";
          io_block->field_add(Ioss::Field("connectivity_face",
                                          io_block->field_int_type(), storage, Ioss::Field::MESH,
                                          block.num_entry));
        }
        // X -> Edge?
        if (block.num_edges_per_entry > 0 && rank_offset < 2) {
          std::string storage = "Real["+Ioss::Utils::to_string(block.num_edges_per_entry)+"]";
          io_block->field_add(Ioss::Field("connectivity_edge",
                                          io_block->field_int_type(), storage, Ioss::Field::MESH,
                                          block.num_entry));
        }
#endif	

        // Maintain block order on output database...
        io_block->property_add(Ioss::Property("original_block_order", iblk));

        if (save_type != "null" && save_type != "") {
          if (io_block->property_exists("original_topology_type")) {
            if (io_block->get_property("original_topology_type").get_string() != save_type) {
              io_block->property_erase("original_topology_type");
              io_block->property_add(Ioss::Property("original_topology_type", save_type));
            }
          } else {
            if (io_block->get_property("topology_type").get_string() != save_type) {
              // Maintain original X type on output database if possible.
              io_block->property_add(Ioss::Property("original_topology_type", save_type));
            }
          }
        }

        io_block->property_add(Ioss::Property("global_entity_count", (int64_t)decomp->el_blocks[iblk].ioss_count()));

        // See if this block is "omitted" by the calling code.
        // This only affects the generation of surfaces...
        if (!blockOmissions.empty()) {
          std::vector<std::string>::const_iterator I = std::find(blockOmissions.begin(),
                                                                 blockOmissions.end(),
                                                                 block_name);
          if (I != blockOmissions.end()) {
            io_block->property_add(Ioss::Property(std::string("omitted"), 1));
          }
        }
        if (block_name != alias) {
          get_region()->add_alias(block_name, alias);
        }

        // Check for additional variables.
        add_attribute_fields(entity_type, io_block, decomp->el_blocks[iblk].attributeCount, type);
        add_results_fields(entity_type, io_block, iblk);

        if (entity_type == EX_ELEM_BLOCK) {
	  Ioex::add_map_fields(get_file_pointer(), (Ioss::ElementBlock*)io_block,
			       decomp->el_blocks[iblk].ioss_count(), maximumNameLength);
        }
      }
    }
  }

  void DatabaseIO::compute_block_adjacencies() const
  {
    // Add a field to each element block specifying which other element
    // blocks the block is adjacent to (defined as sharing nodes).
    // This is only calculated on request...

    blockAdjacenciesCalculated = true;

    if (m_groupCount[EX_ELEM_BLOCK] == 1) {
      blockAdjacency.resize(m_groupCount[EX_ELEM_BLOCK]);
      blockAdjacency[0].resize(m_groupCount[EX_ELEM_BLOCK]);
      blockAdjacency[0][0] = 0;
      return;
    }

    // Each processor first calculates the adjacencies on their own
    // processor...

    std::vector<int64_t> node_used(nodeCount);
    std::vector<std::vector<int> > inv_con(nodeCount);
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    assert(check_block_order(element_blocks));

    {
      for (int iblk = 0; iblk < m_groupCount[EX_ELEM_BLOCK]; iblk++) {
        Ioss::ElementBlock *eb = element_blocks[iblk];
        int blk_position =  eb->get_property("original_block_order").get_int();
        int64_t id =        eb->get_property("id").get_int();
        int element_nodes = eb->get_property("topology_node_count").get_int();
        int64_t my_element_count = eb->get_property("entity_count").get_int();
        int order = eb->get_property("original_block_order").get_int();
        if (my_element_count > 0) {
          if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
            std::vector<int64_t> conn(my_element_count * element_nodes);
            decomp64->get_block_connectivity(get_file_pointer(), TOPTR(conn), id, order, element_nodes);

            for (int64_t i=0; i < my_element_count * element_nodes; i++) {
              node_used[conn[i]-1] = blk_position+1;
            }
          } else {
            std::vector<int> conn(my_element_count * element_nodes);
            decomp32->get_block_connectivity(get_file_pointer(), TOPTR(conn), id, order, element_nodes);

            for (int64_t i=0; i < my_element_count * element_nodes; i++) {
              node_used[conn[i]-1] = blk_position+1;
            }
          }

          for (int64_t i=0; i < nodeCount; i++) {
            if (node_used[i] == blk_position+1) {
              inv_con[i].push_back(blk_position);
            }
          }
        }
      }
    }

    if (isParallel) {
      // Get contributions from other processors...
      // Get the communication map...
      Ioss::CommSet  *css = get_region()->get_commset("commset_node");
      Ioex::check_non_null(css, "communication map", "commset_node");
      std::vector<std::pair<int,int> > proc_node;
      {
        std::vector<int> entity_processor;
        css->get_field_data("entity_processor", entity_processor);
        proc_node.reserve(entity_processor.size()/2);
        size_t j=0;
        for (size_t i=0; i < entity_processor.size(); j++, i+=2) {
          proc_node.push_back(std::make_pair(entity_processor[i+1], entity_processor[i]));
        }
      }

      // Now sort by increasing processor number.
      std::sort(proc_node.begin() , proc_node.end());

      // Pack the data: global_node_id, bits for each block, ...
      // Use 'int' as basic type...
      size_t id_size  = 1;
      size_t word_size = sizeof(int) * 8;
      size_t bits_size = (m_groupCount[EX_ELEM_BLOCK]+word_size-1)/word_size;
      std::vector<unsigned> send(proc_node.size() * (id_size + bits_size));
      std::vector<unsigned> recv(proc_node.size() * (id_size + bits_size));

      std::vector<int> procs(util().parallel_size());
      size_t offset = 0;
      for (size_t i=0; i < proc_node.size(); i++) {
        int64_t glob_id = proc_node[i].second;
        int proc    = proc_node[i].first;
        procs[proc]++;
        send[offset++] = glob_id;
        int64_t loc_id = nodeMap.global_to_local(glob_id, true) -1;
        for (size_t j = 0; j < inv_con[loc_id].size(); j++) {
          int jblk = inv_con[loc_id][j];
          size_t wrd_off = jblk / word_size;
          size_t bit     = jblk % word_size;
          send[offset+wrd_off] |= (1 << bit);
        }
        offset += bits_size;
      }

      // Count nonzero entries in 'procs' array -- count of
      // sends/receives
      size_t non_zero = util().parallel_size() - std::count(procs.begin(), procs.end(), 0);

      // Post all receives...
      MPI_Request request_null = MPI_REQUEST_NULL;
      std::vector<MPI_Request> request(non_zero, request_null);
      std::vector<MPI_Status>  status(non_zero);

      int result = MPI_SUCCESS;
      size_t req_cnt = 0;
      offset = 0;
      for (int i = 0 ; result == MPI_SUCCESS && i < util().parallel_size(); i++) {
        if (procs[i] > 0) {
          const unsigned size = procs[i] * (id_size + bits_size);
          void * const   recv_buf  = &recv[offset];
          result = MPI_Irecv(recv_buf, size, MPI_INT,
			     i, 10101, util().communicator(), &request[req_cnt]);
          req_cnt++;
          offset += size;
        }
      }
      assert(result != MPI_SUCCESS || non_zero == req_cnt);

      if (result != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Irecv error on processor " << util().parallel_rank()
	       << " in Iopx::DatabaseIO::compute_block_adjacencies";
        std::cerr << errmsg.str();
      }

      int local_error = (MPI_SUCCESS == result) ? 0 : 1 ;
      int global_error = util().global_minmax(local_error, Ioss::ParallelUtils::DO_MAX);

      if (global_error != 0) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Irecv error on some processor "
	       << "in Iopx::DatabaseIO::compute_block_adjacencies";
        IOSS_ERROR(errmsg);
      }

      result = MPI_SUCCESS;
      req_cnt = 0;
      offset = 0;
      for (int i = 0 ; result == MPI_SUCCESS && i < util().parallel_size(); i++) {
        if (procs[i] > 0) {
          const unsigned size = procs[i] * (id_size + bits_size);
          void * const   send_buf  = &send[offset];
          result = MPI_Rsend(send_buf, size, MPI_INT,
			     i, 10101, util().communicator());
          req_cnt++;
          offset += size;
        }
      }
      assert(result != MPI_SUCCESS || non_zero == req_cnt);

      if (result != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Rsend error on processor " << util().parallel_rank()
	       << " in Iopx::DatabaseIO::compute_block_adjacencies";
        std::cerr << errmsg.str();
      }

      local_error = (MPI_SUCCESS == result) ? 0 : 1 ;
      global_error = util().global_minmax(local_error, Ioss::ParallelUtils::DO_MAX);

      if (global_error != 0) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Rsend error on some processor "
	       << "in Iopx::DatabaseIO::compute_block_adjacencies";
        IOSS_ERROR(errmsg);
      }

      result = MPI_Waitall(req_cnt, TOPTR(request), TOPTR(status));

      if (result != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Waitall error on processor " << util().parallel_rank()
	       << " in Iopx::DatabaseIO::compute_block_adjacencies";
        std::cerr << errmsg.str();
      }

      // Unpack the data and update the inv_con arrays for boundary
      // nodes...
      offset = 0;
      for (size_t i=0; i < proc_node.size(); i++) {
        int64_t glob_id = recv[offset++];
        int64_t loc_id = nodeMap.global_to_local(glob_id, true) -1;
        for (int iblk = 0; iblk < m_groupCount[EX_ELEM_BLOCK]; iblk++) {
          size_t wrd_off = iblk / word_size;
          size_t bit     = iblk % word_size;
          if (recv[offset+wrd_off] & (1 << bit)) {
            inv_con[loc_id].push_back(iblk); // May result in duplicates, but that is OK.
          }
        }
        offset += bits_size;
      }
    }

    // Convert from inv_con arrays to block adjacency...
    blockAdjacency.resize(m_groupCount[EX_ELEM_BLOCK]);
    for (int iblk = 0; iblk < m_groupCount[EX_ELEM_BLOCK]; iblk++) {
      blockAdjacency[iblk].resize(m_groupCount[EX_ELEM_BLOCK]);
    }

    for (int64_t i = 0; i < nodeCount; i++) {
      for (size_t j = 0; j < inv_con[i].size(); j++) {
        int jblk = inv_con[i][j];
        for (size_t k = j+1; k < inv_con[i].size(); k++) {
          int kblk = inv_con[i][k];
          blockAdjacency[jblk][kblk] = 1;
          blockAdjacency[kblk][jblk] = 1;
        }
      }
    }

    if (isParallel) {
      // Sync across all processors...
      size_t word_size = sizeof(int) * 8;
      size_t bits_size = (m_groupCount[EX_ELEM_BLOCK]+word_size-1)/word_size;

      std::vector<unsigned> data(m_groupCount[EX_ELEM_BLOCK] * bits_size);
      int64_t offset = 0;
      for (int jblk = 0; jblk < m_groupCount[EX_ELEM_BLOCK]; jblk++) {
        for (int iblk = 0; iblk < m_groupCount[EX_ELEM_BLOCK]; iblk++) {
          if (blockAdjacency[jblk][iblk] == 1) {
            size_t wrd_off = iblk / word_size;
            size_t bit     = iblk % word_size;
            data[offset+wrd_off] |= (1 << bit);
          }
        }
        offset += bits_size;
      }

      std::vector<unsigned> out_data(m_groupCount[EX_ELEM_BLOCK] * bits_size);
      MPI_Allreduce((void*)TOPTR(data), TOPTR(out_data), static_cast<int>(data.size()),
                    MPI_UNSIGNED, MPI_BOR, util().communicator());

      offset = 0;
      for (int jblk = 0; jblk < m_groupCount[EX_ELEM_BLOCK]; jblk++) {
        for (int iblk = 0; iblk < m_groupCount[EX_ELEM_BLOCK]; iblk++) {
          if (blockAdjacency[jblk][iblk] == 0) {
            size_t wrd_off = iblk / word_size;
            size_t bit     = iblk % word_size;
            if (out_data[offset+wrd_off] & (1 << bit)) {
              blockAdjacency[jblk][iblk] = 1;
            }
          }
        }
        offset += bits_size;
      }
    }

    // Make it symmetric...
    for (int iblk = 0; iblk < m_groupCount[EX_ELEM_BLOCK]; iblk++) {
      for (int jblk = iblk; jblk < m_groupCount[EX_ELEM_BLOCK]; jblk++) {
        blockAdjacency[jblk][iblk] = blockAdjacency[iblk][jblk];
      }
    }
  }

  void DatabaseIO::compute_node_status() const
  {
    // Create a field for all nodes in the model indicating
    // the connectivity 'status' of the node.  The status values are:
    // 0 -- node not connected to any elements
    // 1 -- node only connected to omitted elements
    // 2 -- node only connected to active elements
    // 3 -- node at border of active and omitted elements.

    /// \todo Get working for parallel...

    if (nodeConnectivityStatusCalculated)
      return;

    nodeConnectivityStatus.resize(nodeCount);

    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    assert(check_block_order(element_blocks));

    for (int i=0; i < m_groupCount[EX_ELEM_BLOCK]; i++) {
      Ioss::ElementBlock *block = element_blocks[i];
      unsigned char status = 2;
      if (Ioss::Utils::block_is_omitted(block)) {
        status = 1;
      }

      int64_t id = block->get_property("id").get_int();
      int element_nodes = block->get_property("topology_node_count").get_int();
      int64_t my_element_count = block->get_property("entity_count").get_int();
      int order = block->get_property("original_block_order").get_int();
      if (my_element_count > 0) {
        if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
          std::vector<int64_t> conn(my_element_count * element_nodes);
          decomp64->get_block_connectivity(get_file_pointer(), TOPTR(conn), id, order, element_nodes);
          for (int64_t j=0; j < my_element_count * element_nodes; j++) {
            nodeConnectivityStatus[conn[j]-1] |= status;
          }
        } else {
          std::vector<int> conn(my_element_count * element_nodes);
          decomp32->get_block_connectivity(get_file_pointer(), TOPTR(conn), id, order, element_nodes);
          for (int64_t j=0; j < my_element_count * element_nodes; j++) {
            nodeConnectivityStatus[conn[j]-1] |= status;
          }
        }
      }
    }
    nodeConnectivityStatusCalculated = true;
  }

  void DatabaseIO::get_sidesets()
  {
    // This function creates all sidesets (surfaces) for a
    // model.  Note that a sideset contains 1 or more sideblocks
    // which are homogenous (same topology). In serial execution,
    // this is fairly straightforward since there are no null sets and
    // we have all the information we need. (...except see below for
    // surface evolution).
    //
    // However, in a parallel execution, we have the possibility that a
    // side set will have no sides or distribution factors on
    // a particular processor.  We then don't know the block topology of
    // the block(s) contained in this set. We could do some
    // communication and get a good idea of the topologies that are in
    // the set.

    if (m_groupCount[EX_SIDE_SET] > 0) {
      check_side_topology();

      // Get exodusII sideset metadata

      // Get the names (may not exist) of all sidesets and see if they are actually
      // side "blocks" (perhaps written by IO system for a restart).  In that case,
      // they were split by a previous run and we need to reconstruct the side "set"
      // that may contain one or more of them.
      Ioex::SideSetMap fs_map;
      Ioex::SideSetSet fs_set;

      {
        int error;

        for (int i = 0; i < m_groupCount[EX_SIDE_SET]; i++) {
          int64_t id = decomp->side_sets[i].id();
          std::vector<char> ss_name(maximumNameLength+1);
          error = ex_get_name(get_file_pointer(), EX_SIDE_SET, id, TOPTR(ss_name));
          if (error < 0) {
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }
          if (ss_name[0] != '\0') {
            Ioss::Utils::fixup_name(TOPTR(ss_name));
	    Ioex::decode_surface_name(fs_map, fs_set, TOPTR(ss_name));
          }
        }
      }

      // Create sidesets for each entry in the fs_set... These are the
      // sidesets which were probably written by a previous run of the
      // IO system and are already split into homogenous pieces...
      {
        Ioex::SideSetSet::iterator I = fs_set.begin();
        while (I != fs_set.end()) {
          const std::string fs_name = *I;
          Ioss::SideSet *side_set = new Ioss::SideSet(this, fs_name);
          get_region()->add(side_set);
          int64_t id = Ioex::extract_id(fs_name);
          if (id > 0)
            side_set->property_add(Ioss::Property("id", id));
          ++I;
        }
      }

      for (int iss = 0; iss < m_groupCount[EX_SIDE_SET]; iss++) {
        int64_t id = decomp->side_sets[iss].id();
        std::string sid = "";
	Ioex::TopologyMap topo_map;
	Ioex::TopologyMap side_map; // Used to determine side consistency

        Ioss::SurfaceSplitType split_type = splitType;
        std::string side_set_name;
        Ioss::SideSet *side_set = NULL;

        {
	  bool db_has_name = false;
          side_set_name = Ioex::get_entity_name(get_file_pointer(), EX_SIDE_SET, id, "surface",
						maximumNameLength, db_has_name);

	  std::string alias = Ioss::Utils::encode_entity_name("surface", id);

          if (side_set_name == "universal_sideset") {
            split_type = Ioss::SPLIT_BY_DONT_SPLIT;
          }

          bool in_fs_map = false;
	  Ioex::SideSetMap::iterator FSM = fs_map.find(side_set_name);
          if (FSM != fs_map.end()) {
            in_fs_map = true;
            std::string efs_name = (*FSM).second;
            side_set = get_region()->get_sideset(efs_name);
            Ioex::check_non_null(side_set, "sideset", efs_name);
          } else {
	    if (get_use_generic_canonical_name()) {
	      std::string temp = side_set_name;
	      side_set_name = alias;
	      alias = temp;
	    }
            side_set = new Ioss::SideSet(this, side_set_name);
            side_set->property_add(Ioss::Property("id", id));
	    if (db_has_name) {
	      std::string *db_name = &side_set_name;
	      if (get_use_generic_canonical_name()) {
		db_name = &alias;
	      }
	      side_set->property_add(Ioss::Property("db_name", *db_name));
	    }
            get_region()->add((Ioss::SideSet*)side_set);

            get_region()->add_alias(side_set_name, alias);
            get_region()->add_alias(side_set_name, Ioss::Utils::encode_entity_name("sideset", id));
          }

          //	  split_type = SPLIT_BY_ELEMENT_BLOCK;
          //	  split_type = SPLIT_BY_TOPOLOGIES;
          //	  split_type = SPLIT_BY_DONT_SPLIT;

          // Determine how many side blocks compose this side set.

          int64_t number_sides = decomp->side_sets[iss].ioss_count();
	  // FIXME: Support-  number_distribution_factors = decomp->side_sets[iss].df_count();

          Ioss::Int64Vector element(number_sides);
          Ioss::Int64Vector sides(number_sides);

          // Easier below here if the element and sides are a known 64-bit size...
          // Kluge here to do that...
          if (int_byte_size_api() == 4) {
            Ioss::Field side_field("sides",   Ioss::Field::INTEGER, SCALAR(), Ioss::Field::MESH, number_sides);
            Ioss::Field elem_field("ids_raw", Ioss::Field::INTEGER, SCALAR(), Ioss::Field::MESH, number_sides);

            Ioss::IntVector e32(number_sides);
            decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, TOPTR(e32));
            std::copy(e32.begin(), e32.end(), sides.begin());
            decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, TOPTR(e32));
            std::copy(e32.begin(), e32.end(), element.begin());
          } else {
            Ioss::Field side_field("sides",   Ioss::Field::INT64, SCALAR(), Ioss::Field::MESH, number_sides);
            Ioss::Field elem_field("ids_raw", Ioss::Field::INT64, SCALAR(), Ioss::Field::MESH, number_sides);
            decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, TOPTR(sides));
            decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, TOPTR(element));
          }

          if (!blockOmissions.empty()) {
	    Ioex::filter_element_list(get_region(), element, sides, true);
            number_sides = element.size();
            assert(element.size() == sides.size());
          }

          if (split_type == Ioss::SPLIT_BY_TOPOLOGIES && sideTopology.size() == 1) {
            // There is only one side type for all elements in the model
            topo_map[std::make_pair(sideTopology[0].first->name(),
                                    sideTopology[0].second)] = number_sides;

          } else if (split_type == Ioss::SPLIT_BY_DONT_SPLIT) {
            const Ioss::ElementTopology *mixed_topo = Ioss::ElementTopology::factory("unknown");
            topo_map[std::make_pair(std::string("unknown"),mixed_topo)] = number_sides;

          } else if (in_fs_map) {
            std::vector<std::string> tokens;
            Ioss::tokenize(side_set_name, "_", tokens);
            assert(tokens.size() >= 4);
            // The sideset should have only a single topology which is
            // given by the sideset name...
            const Ioss::ElementTopology *side_topo    = Ioss::ElementTopology::factory(tokens[tokens.size()-2]);
            assert(side_topo    != NULL);
            const Ioss::ElementTopology *element_topo = Ioss::ElementTopology::factory(tokens[tokens.size()-3], true);
            std::string name;
            if (element_topo != NULL) {
              name = element_topo->name();
            } else {
              //                           -4   -3   -2     -1
              // Name is of the form name_block_id_sidetopo_id
              name = tokens[tokens.size()-4] + "_" + tokens[tokens.size()-3];
            }

            topo_map[std::make_pair(name, side_topo)] = number_sides;

            // We want the id to match the id on the sideset in this
            // case so that the generated name will match the current
            // name.  Instead of converting from string to int back to
            // string, we just set a variable to query later.
            sid = tokens[tokens.size()-1];

          } else if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
            // There are multiple side types in the model.
            // Iterate through the elements in the sideset, determine
            // their parent element block using the blocks element
            // topology and the side number, determine the side
            // type.

            for (Ioss::TopoContainer::size_type i=0; i < sideTopology.size(); i++) {
              topo_map[std::make_pair(sideTopology[i].first->name(), sideTopology[i].second)] = 0;
              side_map[std::make_pair(sideTopology[i].first->name(), sideTopology[i].second)] = 0;
            }

	    Ioex::separate_surface_element_sides(element, sides, get_region(), topo_map, side_map, split_type);
          }
          else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
            // There are multiple side types in the model.  Iterate
            // through the elements in the sideset, determine their parent
            // element block using blocks element topology and the side
            // number, determine the side type.

            // Seed the topo_map map with <block->name, side_topo>
            // pairs so we are sure that all processors have the same
            // starting topo_map (size and order).
            Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
            assert(check_block_order(element_blocks));

            for (int i=0; i < m_groupCount[EX_ELEM_BLOCK]; i++) {
              Ioss::ElementBlock *block = element_blocks[i];
              if (!Ioss::Utils::block_is_omitted(block)) {
                std::string name = block->name();
                const Ioss::ElementTopology *common_ftopo = block->topology()->boundary_type(0);
                if (common_ftopo != NULL) {
                  // All sides of this element block's topology have the same topology
                  topo_map[std::make_pair(name, common_ftopo)] = 0;
                  side_map[std::make_pair(name, common_ftopo)] = 0;
                } else {
                  // The sides have different topology, iterate over
                  // them and create an entry for the unique side
                  // topology types
                  int par_dim = block->topology()->parametric_dimension();
                  if (par_dim == 2 || par_dim == 3) {
                    int64_t my_side_count = block->topology()->number_boundaries();
                    for (int64_t ii = 0; ii < my_side_count; ii++) {
                      const Ioss::ElementTopology *topo = block->topology()->boundary_type(ii+1);
                      topo_map[std::make_pair(name, topo)] = 0;
                      side_map[std::make_pair(name, topo)] = 0;
                    }
                  }
                }
              }
            }
	    Ioex::separate_surface_element_sides(element, sides, get_region(), topo_map, side_map, split_type);
          }
        }

        // End of first step in splitting.  Check among all processors
        // to see which potential splits have sides in them...
        Ioss::Int64Vector global_side_counts(topo_map.size());
        {
          int64_t i = 0;
          {
	    Ioex::TopologyMap::const_iterator I = topo_map.begin();
            while (I != topo_map.end()) {
              global_side_counts[i++] = (*I).second;
              ++I;
            }
          }

          // If splitting by element block, also sync the side_map
          // information which specifies whether the sideset has
          // consistent sides for all elements. Only really used for
          // shells, but easier to just set the value on all surfaces
          // in the element block split case.
          if (side_map.size() == topo_map.size()) {
            global_side_counts.resize(topo_map.size() + side_map.size());

            {
	      Ioex::TopologyMap::const_iterator I = side_map.begin();
              while (I != side_map.end()) {
                global_side_counts[i++] = (*I).second;
                ++I;
              }
            }
          }

          // See if any processor has non-zero count for the topo_map counts
          // For the side_map, need the max value.
          util().global_array_minmax(global_side_counts, Ioss::ParallelUtils::DO_MAX);
        }


        // Create Side Blocks

        int64_t i = 0;
	Ioex::TopologyMap::const_iterator I = topo_map.begin();
        assert(I != topo_map.end());
        while (I != topo_map.end()) {
          if (global_side_counts[i++] > 0) {
            const std::string topo_or_block_name   = (*I).first.first;
            const Ioss::ElementTopology *side_topo = (*I).first.second;
            assert(side_topo != NULL);
#if 0
            if (side_topo->parametric_dimension() == topology_dimension-1 ||
                split_type == Ioss::SPLIT_BY_DONT_SPLIT ) {
#else
              if (true) {
#endif
                int64_t my_side_count = (*I).second;

                std::string side_block_name = "surface_" + topo_or_block_name + "_" + side_topo->name();
                if (side_set_name == "universal_sideset") {
                  side_block_name = side_set_name;
                } else {
                  if (sid == "")
                    side_block_name = Ioss::Utils::encode_entity_name(side_block_name, id);
                  else {
                    side_block_name += "_";
                    side_block_name += sid;
                  }
                }

                Ioss::ElementBlock* block = NULL;
                // Need to get elem_topo....
                const Ioss::ElementTopology *elem_topo = NULL;
                if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
                  elem_topo = Ioss::ElementTopology::factory(topo_or_block_name);
                }
                else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
                  block = get_region()->get_element_block(topo_or_block_name);
                  if (block == NULL || Ioss::Utils::block_is_omitted(block)) {
                    std::ostringstream errmsg;
                    errmsg << "INTERNAL ERROR: Could not find element block '" << topo_or_block_name
			   << "' Something is wrong in the Iopx::DatabaseIO class. Please report.\n";
                    IOSS_ERROR(errmsg);
                  }
                  elem_topo = block->topology();
                }
                if (split_type == Ioss::SPLIT_BY_DONT_SPLIT) {
                  // Most likely this is "unknown", but can be a true
                  // topology if there is only a single element block in
                  // the model.
                  elem_topo = Ioss::ElementTopology::factory(topo_or_block_name);
                }
                assert(elem_topo != NULL);

                Ioss::SideBlock *side_block = new Ioss::SideBlock(this, side_block_name,
                                                                  side_topo->name(),
                                                                  elem_topo->name(),
                                                                  my_side_count);
                assert(side_block != NULL);
                side_block->property_add(Ioss::Property("id", id));
                side_set->add(side_block);

                // Note that all sideblocks within a specific
                // sideset might have the same id.

                // If splitting by element block, need to set the
                // element block member on this side block.
                if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
                  side_block->set_parent_element_block(block);

                }

                // If we calculated whether the element side is
                // consistent for all sides in this block, then
                // tell the block which side it is, or that they are
                // inconsistent. If it wasn't calculated above, then it
                // will be calculated on the fly when/if requested.
                // This is to avoid reading the sideset bulk data in
                // cases where we don't need to read it, but if we are
                // already reading it (to split the sidesets), then use
                // the data when we have it.
                if (!side_map.empty()) {
                  // Set a property indicating which element side
                  // (1-based) all sides in this block are applied to.
                  // If they are not all assigned to the same element
                  // side, indicate this with a side equal to 0.
                  //
                  // (note: 'i' has already been incremented earlier in
                  // the loop.  We need previous value here...)
                  int side = global_side_counts[i-1+topo_map.size()];
                  if (side == 999) side = 0;
                  assert(side <= elem_topo->number_boundaries());
                  side_block->set_consistent_side_number(side);
                }

                // Add an alias...
                get_region()->add_alias(side_block);

                if (split_type != Ioss::SPLIT_BY_DONT_SPLIT
                    && side_set_name != "universal_sideset") {
                  std::string storage = "Real[";
                  storage += Ioss::Utils::to_string(side_topo->number_nodes());
                  storage += "]";
                  side_block->field_add(Ioss::Field("distribution_factors",
                                                    Ioss::Field::REAL, storage,
                                                    Ioss::Field::MESH, my_side_count));
                }

                if (side_set_name == "universal_sideset") {
                  side_block->field_add(Ioss::Field("side_ids",
                                                    side_block->field_int_type(), "scalar",
                                                    Ioss::Field::MESH, my_side_count));
                }

                int num_attr = 0;
                {
                  int ierr = ex_get_attr_param(get_file_pointer(), EX_SIDE_SET, 1, &num_attr);
                  if (ierr < 0)
                    Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
                }
                // Add additional fields
                add_attribute_fields(EX_SIDE_SET, side_block, num_attr, "");
                add_results_fields(EX_SIDE_SET, side_block, iss);
              }
            }
            ++I;
          }
        }
      }
    }

    void DatabaseIO::compute_block_membership(int64_t id, std::vector<std::string> &block_membership) const
    {
      Ioss::Int64Vector block_ids(m_groupCount[EX_ELEM_BLOCK]);
      if (m_groupCount[EX_ELEM_BLOCK] == 1) {
        block_ids[0] = 1;
      } else {
        {
          ex_set set_param[1];
          set_param[0].id = id;
          set_param[0].type = EX_SIDE_SET;
          set_param[0].entry_list = NULL;
          set_param[0].extra_list = NULL;
          set_param[0].distribution_factor_list = NULL;
          int ierr = ex_get_sets(get_file_pointer(), 1, set_param);
          if (ierr < 0) {
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }

          int64_t number_sides = set_param[0].num_entry;

          // Get the element and element side lists.
          if (int_byte_size_api() == 4) {
            Ioss::IntVector element(number_sides);
            Ioss::IntVector sides(number_sides);

            Ioss::Field side_field("sides",   Ioss::Field::INTEGER, SCALAR(), Ioss::Field::MESH, number_sides);
            Ioss::Field elem_field("ids_raw", Ioss::Field::INTEGER, SCALAR(), Ioss::Field::MESH, number_sides);
            decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, TOPTR(sides));
            decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, TOPTR(element));

            Ioss::ElementBlock *block = NULL;
            for (int64_t iel = 0; iel < number_sides; iel++) {
              int64_t elem_id = element[iel];
              if (block == NULL || !block->contains(elem_id)) {
                block = get_region()->get_element_block(elem_id);
                assert(block != NULL);
                int block_order = block->get_property("original_block_order").get_int();
                block_ids[block_order] = 1;
              }
            }
          } else {
            Ioss::Int64Vector element(number_sides);
            Ioss::Int64Vector sides(number_sides);
            Ioss::Field side_field("sides",   Ioss::Field::INT64, SCALAR(), Ioss::Field::MESH, number_sides);
            Ioss::Field elem_field("ids_raw", Ioss::Field::INT64, SCALAR(), Ioss::Field::MESH, number_sides);
            decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, TOPTR(sides));
            decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, TOPTR(element));

            Ioss::ElementBlock *block = NULL;
            for (int64_t iel = 0; iel < number_sides; iel++) {
              int64_t elem_id = element[iel];
              if (block == NULL || !block->contains(elem_id)) {
                block = get_region()->get_element_block(elem_id);
                assert(block != NULL);
                int block_order = block->get_property("original_block_order").get_int();
                block_ids[block_order] = 1;
              }
            }
          }
        }

        // Synchronize among all processors....
        if (isParallel) {
          util().global_array_minmax(block_ids, Ioss::ParallelUtils::DO_MAX);
        }
      }
      Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
      assert(check_block_order(element_blocks));

      for (int64_t i=0; i < m_groupCount[EX_ELEM_BLOCK]; i++) {
        if (block_ids[i] == 1) {
          Ioss::ElementBlock *block = element_blocks[i];
          if (!Ioss::Utils::block_is_omitted(block)) {
            block_membership.push_back(block->name());
          }
        }
      }
    }

    template <typename T>
      void DatabaseIO::get_sets(ex_entity_type type, int64_t count, const std::string &base,
				const T* /* set_type */)
      {
	// Attributes of a Xset are:
	// -- id
	// -- name
	// -- number of nodes
	// -- number of distribution factors (see next comment)
	// ----the #distribution factors should equal #Xs or 0, any
	//     other value does not make sense. If it is 0, then a substitute
	//     list will be created returning 1.0 for the factor

	// In a parallel execution, it is possible that a Xset will have
	// no Xs or distribution factors on a particular processor...

	// Get exodusII Xset metadata
	for (int64_t ins=0; ins < count; ins++) {
	  int64_t id = decomp->node_sets[ins].id();

	  int num_attr = 0;
	  int ierr = ex_get_attr_param(get_file_pointer(), type, id, &num_attr);
	  if (ierr < 0)
	    Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

	  bool db_has_name = false;
	  std::string Xset_name = Ioex::get_entity_name(get_file_pointer(), type, id, base+"list",
							maximumNameLength, db_has_name);

	  std::string alias = Ioss::Utils::encode_entity_name(base+"list", id);
	  if (get_use_generic_canonical_name()) {
	    std::string temp = Xset_name;
	    Xset_name = alias;
	    alias = temp;
	  }

	  T* Xset = new T(this, Xset_name, decomp->node_sets[ins].ioss_count());
	  Xset->property_add(Ioss::Property("id", id));
	  if (db_has_name) {
	    std::string *db_name = &Xset_name;
	    if (get_use_generic_canonical_name()) {
	      db_name = &alias;
	    }
	    Xset->property_add(Ioss::Property("db_name", *db_name));
	  }
	  get_region()->add(Xset);
	  get_region()->add_alias(Xset_name, alias);
	  get_region()->add_alias(Xset_name, Ioss::Utils::encode_entity_name(base+"set",  id));
	  add_attribute_fields(type, Xset, num_attr, "");
	  add_results_fields(type, Xset, ins);
	}
      }

    void DatabaseIO::get_nodesets()
    {
      get_sets(EX_NODE_SET, m_groupCount[EX_NODE_SET], "node", (Ioss::NodeSet*)NULL);
    }

    void DatabaseIO::get_edgesets()
    {
      get_sets(EX_EDGE_SET, m_groupCount[EX_EDGE_SET], "edge", (Ioss::EdgeSet*)NULL);
    }

    void DatabaseIO::get_facesets()
    {
      get_sets(EX_FACE_SET, m_groupCount[EX_FACE_SET], "face", (Ioss::FaceSet*)NULL);
    }

    void DatabaseIO::get_elemsets()
    {
      get_sets(EX_ELEM_SET, m_groupCount[EX_ELEM_SET], "element", (Ioss::ElementSet*)NULL);
    }

    void DatabaseIO::get_commsets()
    {
      // Attributes of a commset are:
      // -- id (property)
      // -- name (property)
      // -- number of node--CPU pairs (field)

      // In a parallel execution, it is possible that a commset will have
      // no nodes on a particular processor...

      // If this is a serial execution, there will be no communication
      // nodesets, just return an empty container.

      if (isParallel || isSerialParallel) {

        // Create a single node commset
        Ioss::CommSet *commset = new Ioss::CommSet(this, "commset_node", "node",
                                                   decomp->get_commset_node_size());
        commset->property_add(Ioss::Property("id", 1));
        get_region()->add(commset);
      }
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock* nb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);

#ifndef NDEBUG
      int64_t my_node_count = field.raw_count();
      assert(my_node_count == nodeCount);
#endif

      Ioss::Field::RoleType role = field.get_role();
      if (role == Ioss::Field::MESH) {
        if (field.get_name() == "mesh_model_coordinates_x" ||
            field.get_name() == "mesh_model_coordinates_y" ||
            field.get_name() == "mesh_model_coordinates_z" ||
            field.get_name() == "mesh_model_coordinates") {
          decomp->get_node_coordinates(get_file_pointer(), (double*)data, field);
        }

        else if (field.get_name() == "ids") {
          // Map the local ids in this node block
          // (1...node_count) to global node ids.
          get_map(EX_NODE_BLOCK).map_implicit_data(data, field, num_to_get, 0);
        }

	// The 1..global_node_count id.  In a parallel-decomposed run,
	// it maps the node back to its implicit position in the serial
	// undecomposed mesh file.  This is ONLY provided for backward-
	// compatibility and should not be used unless absolutely required.
	else if (field.get_name() == "implicit_ids") {
	  size_t offset = decomp->nodeOffset;
	  size_t count = decomp->nodeCount;
          if (int_byte_size_api() == 4) {
	    std::vector<int> file_ids; file_ids.reserve(count);
	    for (size_t i=0; i<count; i++) {
	      file_ids.push_back(offset+i+1);
	    }
            decomp->communicate_node_data(TOPTR(file_ids), (int*)data, 1);
          } else {
	    std::vector<int64_t> file_ids; file_ids.reserve(count);
	    for (size_t i=0; i<count; i++) {
	      file_ids.push_back(offset+i+1);
	    }
            decomp->communicate_node_data(TOPTR(file_ids), (int64_t*)data, 1);
          }
	}

        else if (field.get_name() == "connectivity") {
          // Do nothing, just handles an idiosyncracy of the GroupingEntity
        }
        else if (field.get_name() == "connectivity_raw") {
          // Do nothing, just handles an idiosyncracy of the GroupingEntity
        }
        else if (field.get_name() == "node_connectivity_status") {
          compute_node_status();
          char *status = static_cast<char*>(data);
          std::copy(nodeConnectivityStatus.begin(), nodeConnectivityStatus.end(), status);
        }

        else if (field.get_name() == "owning_processor") {
          // If parallel, then set the "locally_owned" property on the nodeblocks.
          Ioss::CommSet *css = get_region()->get_commset("commset_node");
          if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
            int64_t *idata = static_cast<int64_t*>(data);
            for (int64_t i=0; i < nodeCount; i++) {
              idata[i] = myProcessor;
            }

            std::vector<int64_t> ent_proc;
            css->get_field_data("entity_processor_raw", ent_proc);
            for (size_t i=0; i < ent_proc.size(); i+=2) {
              int64_t node = ent_proc[i+0];
              int64_t proc = ent_proc[i+1];
              if (proc < idata[node-1]) {
                idata[node-1] = proc;
              }
            }
          }
          else {
            int *idata = static_cast<int*>(data);
            for (int64_t i=0; i < nodeCount; i++) {
              idata[i] = myProcessor;
            }

            std::vector<int> ent_proc;
            css->get_field_data("entity_processor_raw", ent_proc);
            for (size_t i=0; i < ent_proc.size(); i+=2) {
              int node = ent_proc[i+0];
              int proc = ent_proc[i+1];
              if (proc < idata[node-1]) {
                idata[node-1] = proc;
              }
            }
          }
        }

        else {
          num_to_get = Ioss::Utils::field_warning(nb, field, "input");
        }

      } else if (role == Ioss::Field::TRANSIENT) {
        // Check if the specified field exists on this node block.
        // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
        // exist on the database as scalars with the appropriate
        // extensions.

        // Read in each component of the variable and transfer into
        // 'data'.  Need temporary storage area of size 'number of
        // nodes in this block.
        num_to_get = read_transient_field(EX_NODE_BLOCK, m_variables[EX_NODE_BLOCK], field, nb, data);

      } else if (role == Ioss::Field::ATTRIBUTE) {
        num_to_get = read_attribute_field(EX_NODE_BLOCK, field, nb, data);
      }
      return num_to_get;
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {

      size_t num_to_get = field.verify(data_size);

      int64_t id = Ioex::get_id(eb, EX_ELEM_BLOCK, &ids_);
      size_t my_element_count = eb->get_property("entity_count").get_int();
      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
        // Handle the MESH fields required for an ExodusII file model.
        // (The 'genesis' portion)

        if (field.get_name() == "connectivity") {
          int element_nodes = eb->get_property("topology_node_count").get_int();
          assert(field.raw_storage()->component_count() == element_nodes);

          int order = eb->get_property("original_block_order").get_int();
          // The connectivity is stored in a 1D array.
          // The element_node index varies fastet

          decomp->get_block_connectivity(get_file_pointer(), data, id, order, element_nodes);
          get_map(EX_NODE_BLOCK).map_data(data, field, num_to_get*element_nodes);
        }
#if 0
        else if (field.get_name() == "connectivity_face") {
          int face_count = field.raw_storage()->component_count();

          // The connectivity is stored in a 1D array.
          // The element_face index varies fastest
          if (my_element_count > 0) {
            get_connectivity_data(get_file_pointer(), data, EX_ELEM_BLOCK, id, 2);
            get_map(EX_FACE_BLOCK).map_data(data, field, num_to_get*face_count);
          }
        }
        else if (field.get_name() == "connectivity_edge") {
          int edge_count = field.raw_storage()->component_count();

          // The connectivity is stored in a 1D array.
          // The element_edge index varies fastest
          if (my_element_count > 0) {
            get_connectivity_data(get_file_pointer(), data, EX_ELEM_BLOCK, id, 1);
            get_map(EX_EDGE_BLOCK).map_data(data, field, num_to_get*edge_count);
          }
        }
#endif
        else if (field.get_name() == "connectivity_raw") {
          // "connectivity_raw" has nodes in local id space (1-based)
          assert(field.raw_storage()->component_count() == eb->get_property("topology_node_count").get_int());

          // The connectivity is stored in a 1D array.
          // The element_node index varies fastest
          int element_nodes = eb->get_property("topology_node_count").get_int();
          int order = eb->get_property("original_block_order").get_int();
          // FIXME SIZE
          decomp->get_block_connectivity(get_file_pointer(), data, id, order, element_nodes);
        }
        else if (field.get_name() == "ids" || field.get_name() == "implicit_ids") {
          // Map the local ids in this element block
          // (eb_offset+1...eb_offset+1+my_element_count) to global element ids.
          get_map(EX_ELEM_BLOCK).map_implicit_data(data, field, num_to_get, eb->get_offset());
        }
        else if (field.get_name() == "skin") {
          // This is (currently) for the skinned body. It maps the
          // side element on the skin to the original element/local
          // side number.  It is a two component field, the first
          // component is the global id of the underlying element in
          // the initial mesh and its local side number (1-based).

          if (field.is_type(Ioss::Field::INTEGER)) {
            Ioss::IntVector element(my_element_count);
            Ioss::IntVector side(my_element_count);
            int *el_side = (int *)data;

            // FIX: Hardwired map ids....
            size_t eb_offset = eb->get_offset();
            assert(1==0 && "Unimplemented FIXME");
            ex_get_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 1, eb_offset+1, my_element_count,
                                   TOPTR(element)); // FIXME
            ex_get_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 2, eb_offset+1, my_element_count,
                                   TOPTR(side)); // FIXME

            int index = 0;
            for (size_t i=0; i < my_element_count; i++) {
              el_side[index++] = element[i];
              el_side[index++] = side[i];
            }
          } else {
            Ioss::Int64Vector element(my_element_count);
            Ioss::Int64Vector side(my_element_count);
            int64_t *el_side = (int64_t *)data;

            // FIX: Hardwired map ids....
            size_t eb_offset = eb->get_offset();
            assert(1==0 && "Unimplemented FIXME");
            ex_get_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 1, eb_offset+1, my_element_count,
                                   TOPTR(element)); // FIXME
            ex_get_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 2, eb_offset+1, my_element_count,
                                   TOPTR(side)); // FIXME

            size_t index = 0;
            for (size_t i=0; i < my_element_count; i++) {
              el_side[index++] = element[i];
              el_side[index++] = side[i];
            }
          }
        }
        else {
          num_to_get = Ioss::Utils::field_warning(eb, field, "input");
        }

      } else if (role == Ioss::Field::ATTRIBUTE) {
        num_to_get = read_attribute_field(EX_ELEM_BLOCK, field, eb, data);

      } else if (role == Ioss::Field::TRANSIENT) {
        // Check if the specified field exists on this element block.
        // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
        // exist on the database as scalars with the appropriate
        // extensions.

        // Read in each component of the variable and transfer into
        // 'data'.  Need temporary storage area of size 'number of
        // elements in this block.
        num_to_get = read_transient_field(EX_ELEM_BLOCK, m_variables[EX_ELEM_BLOCK], field, eb, data);
      } else if (role == Ioss::Field::REDUCTION) {
        num_to_get = Ioss::Utils::field_warning(eb, field, "input reduction");
      }
      return num_to_get;
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);

      int64_t id = Ioex::get_id(eb, EX_FACE_BLOCK, &ids_);
      size_t my_face_count = eb->get_property("entity_count").get_int();
      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
	// Handle the MESH fields required for an ExodusII file model.
	// (The 'genesis' portion)

	if (field.get_name() == "connectivity") {
	  int face_nodes = eb->get_property("topology_node_count").get_int();
	  assert(field.raw_storage()->component_count() == face_nodes);

	  // The connectivity is stored in a 1D array.
	  // The face_node index varies fastet
	  if (my_face_count > 0) {
	    get_connectivity_data(get_file_pointer(), data, EX_FACE_BLOCK, id, 0);
	    get_map(EX_NODE_BLOCK).map_data(data, field, num_to_get*face_nodes);
	  }
	}
	else if (field.get_name() == "connectivity_edge") {
	  int edge_count = field.raw_storage()->component_count();

	  // The connectivity is stored in a 1D array.
	  // The face_edge index varies fastest
	  if (my_face_count > 0) {
	    get_connectivity_data(get_file_pointer(), data, EX_FACE_BLOCK, id, 1);
	    get_map(EX_EDGE_BLOCK).map_data(data, field, num_to_get*edge_count);
	  }
	}
	else if (field.get_name() == "connectivity_raw") {
	  // "connectivity_raw" has nodes in local id space (1-based)
	  assert(field.raw_storage()->component_count() == eb->get_property("topology_node_count").get_int());

	  // The connectivity is stored in a 1D array.
	  // The face_node index varies fastet
	  if (my_face_count > 0) {
	    get_connectivity_data(get_file_pointer(), data, EX_FACE_BLOCK, id, 0);
	  }
	}
	else if (field.get_name() == "ids") {
	  // Map the local ids in this face block
	  // (eb_offset+1...eb_offset+1+my_face_count) to global face ids.
	  get_map(EX_FACE_BLOCK).map_implicit_data(data, field, num_to_get, eb->get_offset());
	}
	else {
	  num_to_get = Ioss::Utils::field_warning(eb, field, "input");
	}

      } else if (role == Ioss::Field::ATTRIBUTE) {
	num_to_get = read_attribute_field(EX_FACE_BLOCK, field, eb, data);

      } else if (role == Ioss::Field::TRANSIENT) {
	// Check if the specified field exists on this element block.
	// Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	// exist on the database as scalars with the appropriate
	// extensions.

	// Read in each component of the variable and transfer into
	// 'data'.  Need temporary storage area of size 'number of
	// elements in this block.
	num_to_get = read_transient_field(EX_FACE_BLOCK, m_variables[EX_FACE_BLOCK], field, eb, data);
      } else if (role == Ioss::Field::REDUCTION) {
	num_to_get = Ioss::Utils::field_warning(eb, field, "input reduction");
      }
      return num_to_get;
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);

      int64_t id = Ioex::get_id(eb, EX_EDGE_BLOCK, &ids_);
      int64_t my_edge_count = eb->get_property("entity_count").get_int();
      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
	// Handle the MESH fields required for an ExodusII file model.
	// (The 'genesis' portion)

	if (field.get_name() == "connectivity") {
	  int edge_nodes = eb->get_property("topology_node_count").get_int();
	  assert(field.raw_storage()->component_count() == edge_nodes);

	  // The connectivity is stored in a 1D array.
	  // The edge_node index varies fastet
	  if (my_edge_count > 0) {
	    get_connectivity_data(get_file_pointer(), data, EX_EDGE_BLOCK, id, 0);
	    get_map(EX_NODE_BLOCK).map_data(data, field, num_to_get*edge_nodes);
	  }
	}
	else if (field.get_name() == "connectivity_raw") {
	  // "connectivity_raw" has nodes in local id space (1-based)
	  assert(field.raw_storage()->component_count() == eb->get_property("topology_node_count").get_int());

	  // The connectivity is stored in a 1D array.
	  // The edge_node index varies fastet
	  if (my_edge_count > 0) {
	    get_connectivity_data(get_file_pointer(), data, EX_EDGE_BLOCK, id, 0);
	  }
	}
	else if (field.get_name() == "ids") {
	  // Map the local ids in this edge block
	  // (eb_offset+1...eb_offset+1+my_edge_count) to global edge ids.
	  get_map(EX_EDGE_BLOCK).map_implicit_data(data, field, num_to_get, eb->get_offset());
	}
	else {
	  num_to_get = Ioss::Utils::field_warning(eb, field, "input");
	}

      } else if (role == Ioss::Field::ATTRIBUTE) {
	num_to_get = read_attribute_field(EX_EDGE_BLOCK, field, eb, data);

      } else if (role == Ioss::Field::TRANSIENT) {
	// Check if the specified field exists on this element block.
	// Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	// exist on the database as scalars with the appropriate
	// extensions.

	// Read in each component of the variable and transfer into
	// 'data'.  Need temporary storage area of size 'number of
	// elements in this block.
	num_to_get = read_transient_field(EX_EDGE_BLOCK, m_variables[EX_EDGE_BLOCK], field, eb, data);

      } else if (role == Ioss::Field::REDUCTION) {
	num_to_get = Ioss::Utils::field_warning(eb, field, "input reduction");
      }
      return num_to_get;
    }

    int64_t DatabaseIO::get_Xset_field_internal(ex_entity_type type,
                                                const Ioss::EntitySet *ns,
                                                const Ioss::Field& field,
                                                void *data, size_t data_size) const
    {
      int ierr;
      size_t num_to_get = field.verify(data_size);
      Ioss::Field::RoleType role = field.get_role();

      // Find corresponding set in file decomp class...
      if (role == Ioss::Field::MESH) {
        int64_t id = Ioex::get_id(ns, type, &ids_);

        if (field.get_name() == "ids" ||
            field.get_name() == "ids_raw") {
          if (field.get_type() == Ioss::Field::INTEGER) {
            ierr = decomp32->get_set_mesh_var(get_file_pointer(), EX_NODE_SET, id, field,
					      static_cast<int*>(data));
          } else {
            ierr = decomp64->get_set_mesh_var(get_file_pointer(), EX_NODE_SET, id, field,
					      static_cast<int64_t*>(data));
          }
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

          if (field.get_name() == "ids") {
            // Convert the local node ids to global ids
            get_map(EX_NODE_BLOCK).map_data(data, field, num_to_get);
          }
        } else if (field.get_name() == "orientation") {
          if (field.get_type() == Ioss::Field::INTEGER) {
            ierr = decomp32->get_set_mesh_var(get_file_pointer(), EX_NODE_SET, id, field, static_cast<int*>(data));
          } else {
            ierr = decomp64->get_set_mesh_var(get_file_pointer(), EX_NODE_SET, id, field, static_cast<int64_t*>(data));
          }
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        } else if (field.get_name() == "distribution_factors") {
          ierr = decomp->get_set_mesh_double(get_file_pointer(), EX_NODE_SET, id, field, static_cast<double*>(data));
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        }
      } else if (role == Ioss::Field::ATTRIBUTE) {
        num_to_get = read_attribute_field(type, field, ns, data);

      } else if (role == Ioss::Field::TRANSIENT) {
        // Check if the specified field exists on this node block.
        // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
        // exist on the database as scalars with the appropriate
        // extensions.

        // Read in each component of the variable and transfer into
        // 'data'.  Need temporary storage area of size 'number of
        // nodes in this block.
        num_to_get = read_transient_field(type, m_variables[type], field, ns, data);
      }
      return num_to_get;
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet* ns,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      return get_Xset_field_internal(EX_NODE_SET, ns, field, data, data_size);
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::EdgeSet* ns,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      return get_Xset_field_internal(EX_EDGE_SET, ns, field, data, data_size);
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::FaceSet* ns,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      return get_Xset_field_internal(EX_FACE_SET, ns, field, data, data_size);
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::ElementSet* ns,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      return get_Xset_field_internal(EX_ELEM_SET, ns, field, data, data_size);
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::SideSet* fs,
                                           const Ioss::Field& field,
                                           void */* data */, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);
      if (field.get_name() == "ids") {
        // Do nothing, just handles an idiosyncracy of the GroupingEntity
      } else {
        num_to_get = Ioss::Utils::field_warning(fs, field, "input");
      }
      return num_to_get;
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::CommSet* cs,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);

      // Return the <entity (node or side), processor> pair
      if (field.get_name() == "entity_processor" || field.get_name() == "entity_processor_raw") {

	// Check type -- node or side
	std::string type = cs->get_property("entity_type").get_string();

	if (type == "node") {

	  bool do_map = field.get_name() == "entity_processor";
	  // Convert local node id to global node id and store in 'data'
	  const Ioss::MapContainer &map = get_map(EX_NODE_BLOCK).map;
	  if (int_byte_size_api() == 4) {
	    decomp32->get_node_entity_proc_data(static_cast<int*>(data), map, do_map);
	  } else {
	    decomp64->get_node_entity_proc_data(static_cast<int64_t*>(data), map, do_map);
	  }
	} else {
	  std::ostringstream errmsg;
	  errmsg << "ERROR: Invalid commset type " << type;
	  IOSS_ERROR(errmsg);
	}

      } else if (field.get_name() == "ids") {
	// Do nothing, just handles an idiosyncracy of the GroupingEntity
      } else {
	num_to_get = Ioss::Utils::field_warning(cs, field, "input");
      }
      return num_to_get;
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock* fb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      ssize_t num_to_get = field.verify(data_size);
      int ierr = 0;

      int64_t id = Ioex::get_id(fb, EX_SIDE_SET, &ids_);
      int64_t entity_count = fb->get_property("entity_count").get_int();
      if (num_to_get != entity_count) {
        std::ostringstream errmsg;
        errmsg << "ERROR: Partial field input not yet implemented for side blocks";
        IOSS_ERROR(errmsg);
      }

      const SetDecompositionData &set = decomp->get_decomp_set(EX_SIDE_SET, id);

      int64_t number_sides = set.ioss_count();
      int64_t number_distribution_factors = set.df_count();

      Ioss::Field::RoleType role = field.get_role();
      if (role == Ioss::Field::MESH) {


        // In exodusII, we may have split the sideset into multiple
        // side blocks if there are multiple side topologies in the
        // sideset.  Because of this, the passed in 'data' may not be
        // large enough to hold the data residing in the sideset and we
        // may need to allocate a temporary array...  This can be checked
        // by comparing the size of the sideset with the 'side_count' of
        // the side block.

        // Get size of data stored on the file...
        // FIX 64: FIX THIS -- STORING INT IN DOUBLE WON'T WORK
        if (field.get_name() == "side_ids" &&
            fb->name() == "universal_sideset") {
          // The side ids are being stored as the distribution factor
          // field on the universal sideset.  There should be no other
          // side sets that request this field...  (Eventually,
          // create an id field to store this info.

          if (number_distribution_factors == num_to_get) {
            std::vector<double> real_ids(num_to_get);
            Ioss::Field df_field("distribution_factor", Ioss::Field::REAL, "scalar", Ioss::Field::MESH, num_to_get);
            decomp->get_set_mesh_double(get_file_pointer(), EX_SIDE_SET, id, df_field, TOPTR(real_ids));

            // Need to convert 'double' to 'int' for Sierra use...
            int* ids = static_cast<int*>(data);
            for (ssize_t i = 0; i < num_to_get; i++) {
              ids[i] = static_cast<int>(real_ids[i]);
            }
          }
        }

        else if (field.get_name() == "side_ids") {
        }

        else if (field.get_name() == "ids") {
          // In exodusII, the 'side set' is stored as a sideset.  A
          // sideset has a list of elements and a corresponding local
          // element side (1-based) The side id is: side_id =
          // 10*element_id + local_side_number This assumes that all
          // sides in a sideset are boundary sides.  Since we
          // only have a single array, we need to allocate an extra array
          // to store all of the data.  Note also that the element_id is
          // the global id but only the local id is stored so we need to
          // map from local_to_global prior to generating the side id...

          Ioss::Field el_side = fb->get_field("element_side");
          std::vector<char> element_side(2* number_sides * int_byte_size_api());
          get_field_internal(fb, el_side, TOPTR(element_side), element_side.size());

          // At this point, have the 'element_side' data containing
          // the global element ids and the sides...  Iterate
          // through to generate the ids...
          if (int_byte_size_api() == 4) {
            int *ids = static_cast<int*>(data);
            int *els = (int*)TOPTR(element_side);
            size_t idx = 0;
            for (ssize_t iel = 0; iel < 2*entity_count; iel+=2) {
              int64_t new_id = 10*els[iel] + els[iel+1];
              ids[idx++] = new_id;
            }
          } else {
            int64_t *ids = static_cast<int64_t*>(data);
            int64_t *els = (int64_t*)TOPTR(element_side);
            size_t idx = 0;
            for (ssize_t iel = 0; iel < 2*entity_count; iel+=2) {
              int64_t new_id = 10*els[iel] + els[iel+1];
              ids[idx++] = new_id;
            }
          }

        } else if (field.get_name() == "element_side" ||
		   field.get_name() == "element_side_raw") {
          // In exodusII, the 'side set' is stored as a sideset.  A sideset
          // has a list of elements and a corresponding local element side
          // (1-based)

          // Since we only have a single array, we need to allocate an extra
          // array to store all of the data.  Note also that the element_id
          // is the global id but only the local id is stored so we need to
          // map from local_to_global prior to generating the side  id...

          // Get the element number map (1-based)...
          const Ioss::MapContainer &map = get_map(EX_ELEM_BLOCK).map;

          // Allocate space for local side number and element numbers
          // numbers.

          // See if edges or faces...
          int64_t side_offset = Ioss::Utils::get_side_offset(fb);


	  if (fb->owner()->block_count() == 1 && number_sides == entity_count) {

            if (int_byte_size_api() == 4) {
              int     *element_side = static_cast<int*>(data);
              decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, field, element_side);
              for (ssize_t iel = 1; iel < 2*entity_count; iel+=2) {
                element_side[iel] = element_side[iel] - side_offset;
              }
            } else {
              int64_t *element_side = static_cast<int64_t*>(data);
              decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, field, element_side);
              for (ssize_t iel = 1; iel < 2*entity_count; iel+=2) {
                element_side[iel] = element_side[iel] - side_offset;
              }
            }
          } else {
            // Need a larger vector to get the entire sideset and then filter down to the correct size...
            std::vector<char> element(number_sides * int_byte_size_api());
            std::vector<char> sides(number_sides * int_byte_size_api());
            if (int_byte_size_api() == 4) {
              Ioss::Field elem_field("ids",   Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, number_sides);
              Ioss::Field side_field("sides", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, number_sides);
              decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, (int*)TOPTR(element));
              decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, (int*)TOPTR(sides));
            } else {
              Ioss::Field elem_field("ids",   Ioss::Field::INT64, "scalar", Ioss::Field::MESH, number_sides);
              Ioss::Field side_field("sides", Ioss::Field::INT64, "scalar", Ioss::Field::MESH, number_sides);
              decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, (int64_t*)TOPTR(element));
              decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, (int64_t*)TOPTR(sides));
            }

            Ioss::IntVector is_valid_side;
            Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, int_byte_size_api(),
                                                        TOPTR(element), TOPTR(sides),
                                                        number_sides, get_region());

            ssize_t index = 0;
            if (int_byte_size_api() == 4) {
              int     *element_side = static_cast<int*>(data);
              int     *element32 = (int*)TOPTR(element);
              int     *sides32 = (int*)TOPTR(sides);
              for (int64_t iel = 0; iel < number_sides; iel++) {
                if (is_valid_side[iel] == 1) {
                  // This side  belongs in the side block
                  element_side[index++] = element32[iel];
                  element_side[index++] = sides32[iel] - side_offset;
                }
              }
            } else {
              int64_t *element_side = static_cast<int64_t*>(data);
              int64_t *element64 = (int64_t*)TOPTR(element);
              int64_t *sides64 = (int64_t*)TOPTR(sides);
              for (int64_t iel = 0; iel < number_sides; iel++) {
                if (is_valid_side[iel] == 1) {
                  // This side  belongs in the side block
                  element_side[index++] = element64[iel];
                  element_side[index++] = sides64[iel] - side_offset;
                }
              }
            }
            assert(index/2 == entity_count);
          }
          if (field.get_name() == "element_side") {
            if (int_byte_size_api() == 4) {
              int     *element_side = static_cast<int*>(data);
              for (ssize_t iel = 0; iel < 2*entity_count; iel+=2) {
                element_side[iel] = map[element_side[iel]];
              }
            } else {
              int64_t *element_side = static_cast<int64_t*>(data);
              for (ssize_t iel = 0; iel < 2*entity_count; iel+=2) {
                element_side[iel] = map[element_side[iel]];
              }
            }
          }

        } else if (field.get_name() == "connectivity") {
          // The side  connectivity needs to be generated 'on-the-fly' from
          // the element number and local side  of that element. A sideset
          // can span multiple element blocks, and contain multiple side
          // types; the side block contains side of similar topology.
          ierr = get_side_connectivity(fb, id, entity_count, data, true);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        } else if (field.get_name() == "connectivity_raw") {
          // The side  connectivity needs to be generated 'on-the-fly' from
          // the element number and local side  of that element. A sideset
          // can span multiple element blocks, and contain multiple side
          // types; the side block contains side of similar topology.
          ierr = get_side_connectivity(fb, id, entity_count, data, false);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        } else if (field.get_name() == "distribution_factors") {
          ierr = get_side_distributions(fb, id, entity_count,
					static_cast<double*>(data), data_size/sizeof(double));
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        } else {
          num_to_get = Ioss::Utils::field_warning(fb, field, "input");
        }
      } else if (role == Ioss::Field::TRANSIENT) {
        // Check if the specified field exists on this block.
        // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
        // exist on the database as scalars with the appropriate
        // extensions.

	if (fb->owner()->block_count() == 1 && number_sides == entity_count) {
          num_to_get = read_transient_field(EX_SIDE_SET, m_variables[EX_SIDE_SET], field, fb, data);
        } else {
          // Need to read all values for the specified field and then
          // filter down to the elements actually in this side block.

          Ioss::IntVector is_valid_side;
          // Need a larger vector to get the entire sideset and then filter down to the correct size...
          std::vector<char> element(number_sides * int_byte_size_api());
          std::vector<char> sides(number_sides * int_byte_size_api());
          if (int_byte_size_api() == 4) {
            Ioss::Field elem_field("ids",   Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, number_sides);
            Ioss::Field side_field("sides", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, number_sides);
            decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, (int*)TOPTR(element));
            decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, (int*)TOPTR(sides));
          } else {
            Ioss::Field elem_field("ids",   Ioss::Field::INT64, "scalar", Ioss::Field::MESH, number_sides);
            Ioss::Field side_field("sides", Ioss::Field::INT64, "scalar", Ioss::Field::MESH, number_sides);
            decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, (int64_t*)TOPTR(element));
            decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, (int64_t*)TOPTR(sides));
          }
          Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, int_byte_size_api(),
                                                      TOPTR(element), TOPTR(sides),
                                                      number_sides, get_region());

          num_to_get = read_ss_transient_field(field, id, data, is_valid_side);
        }
      }
      return num_to_get;
    }

    int64_t DatabaseIO::write_attribute_field(ex_entity_type type,
                                              const Ioss::Field& field,
                                              const Ioss::GroupingEntity *ge,
                                              void *data) const
    {
      std::string att_name = ge->name() + SEP() + field.get_name();
      ssize_t num_entity = ge->get_property("entity_count").get_int();
      ssize_t offset = field.get_index();

      int64_t id = Ioex::get_id(ge, type, &ids_);
      int attribute_count = ge->get_property("attribute_count").get_int();
      assert(offset > 0);
      assert(offset-1+field.raw_storage()->component_count() <= attribute_count);

      size_t proc_offset = 0;
      if (ge->property_exists("processor_offset"))
        proc_offset = ge->get_property("processor_offset").get_int();
      ssize_t file_count = num_entity;
      if (ge->property_exists("locally_owned_count"))
        file_count = ge->get_property("locally_owned_count").get_int();

      int comp_count = field.raw_storage()->component_count();
      double *rdata = static_cast<double*>(data);

      if (type == EX_NODAL) {
        for (int i=0; i < comp_count; i++) {
          std::vector<double> file_data; file_data.reserve(file_count);
          map_data(nodeOwningProcessor, myProcessor, rdata, file_data, i, comp_count);
          int ierr = ex_put_partial_one_attr(get_file_pointer(), type, id, proc_offset+1, file_count,
					     offset+i, TOPTR(file_data));
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        }
      } else if (type == EX_NODE_SET) {
        for (int i=0; i < comp_count; i++) {
          std::vector<double> file_data; file_data.reserve(file_count);
          map_nodeset_data(nodesetOwnedNodes[ge], rdata, file_data, i, comp_count);
          int ierr = ex_put_partial_one_attr(get_file_pointer(), type, id, proc_offset+1, file_count,
					     offset+i, TOPTR(file_data));
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        }
      } else {
        for (int i=0; i < comp_count; i++) {
          std::vector<double> file_data; file_data.reserve(file_count);
          assert(file_count == num_entity);
          size_t k = i;
          for (ssize_t j=0; j < num_entity; j++) {
            file_data.push_back(rdata[k]);
            k += comp_count;
          }

          int ierr = ex_put_partial_one_attr(get_file_pointer(), type, id, proc_offset+1, file_count,
					     offset+i, TOPTR(file_data));
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        }
      }
      return num_entity;
    }

    int64_t DatabaseIO::read_attribute_field(ex_entity_type type,
                                             const Ioss::Field& field,
                                             const Ioss::GroupingEntity *ge,
                                             void *data) const
    {
      int64_t num_entity = ge->get_property("entity_count").get_int();

      int attribute_count = ge->get_property("attribute_count").get_int();
      int64_t id = Ioex::get_id(ge, type, &ids_);

      std::string att_name = ge->name() + SEP() + field.get_name();
      ssize_t offset = field.get_index();
      assert(offset-1+field.raw_storage()->component_count() <= attribute_count);
      if (offset == 1 && field.raw_storage()->component_count() == attribute_count) {
        // Read all attributes in one big chunk...
        int ierr = decomp->get_attr(get_file_pointer(), type, id, attribute_count, static_cast<double*>(data));
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
      else {
        // Read a subset of the attributes.  If scalar, read one;
        // if higher-order (vector3d, ..) read each component and
        // put into correct location...
        if (field.raw_storage()->component_count() == 1) {
          int ierr = decomp->get_one_attr(get_file_pointer(), type, id,
					  offset, static_cast<double*>(data));
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        } else {
          // Multi-component...
          // Need a local memory space to read data into and
          // then push that into the user-supplied data block...
          std::vector<double> local_data(num_entity);
          int comp_count = field.raw_storage()->component_count();
          double *rdata = static_cast<double*>(data);
          for (int i=0; i < comp_count; i++) {
            int ierr = decomp->get_one_attr(get_file_pointer(), type, id,
					    offset+i, TOPTR(local_data));
            if (ierr < 0)
              Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

            size_t k = i;
            for (ssize_t j=0; j < num_entity; j++) {
              rdata[k] = local_data[j];
              k += comp_count;
            }
          }
        }
      }
      return num_entity;
    }

    int64_t DatabaseIO::read_transient_field(ex_entity_type type,
                                             const Ioex::VariableNameMap &variables,
                                             const Ioss::Field& field,
                                             const Ioss::GroupingEntity *ge,
                                             void *data) const
    {
      const Ioss::VariableType *var_type = field.raw_storage();

      // Read into a double variable since that is all ExodusII can store...
      size_t num_entity = ge->get_property("entity_count").get_int();
      std::vector<double> temp(num_entity);

      size_t step = get_current_state();

      // get number of components, cycle through each component
      // and add suffix to base 'field_name'.  Look up index
      // of this name in 'nodeVariables' map
      size_t comp_count = var_type->component_count();
      size_t var_index=0;

      char field_suffix_separator = get_field_separator();
      for (size_t i=0; i < comp_count; i++) {
        std::string var_name = var_type->label_name(field.get_name(), i+1, field_suffix_separator);

        // Read the variable...
        int64_t id = Ioex::get_id(ge, type, &ids_);
        int ierr = 0;
        var_index = variables.find(var_name)->second;
        assert(var_index > 0);
        ierr = decomp->get_var(get_file_pointer(), step, type,
			       var_index, id, num_entity, temp);
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        // Transfer to 'data' array.
        size_t k = 0;
        if (field.get_type() == Ioss::Field::INTEGER) {
          int *ivar = static_cast<int*>(data);
          for (size_t j=i; j < num_entity*comp_count; j+=comp_count) {
            ivar[j] = static_cast<int>(temp[k++]);
          }
        } else if (field.get_type() == Ioss::Field::INT64) {  // FIX 64 UNSAFE
          int64_t *ivar = static_cast<int64_t*>(data);
          for (size_t j=i; j < num_entity*comp_count; j+=comp_count) {
            ivar[j] = static_cast<int64_t>(temp[k++]);
          }
        } else if (field.get_type() == Ioss::Field::REAL) {
          double *rvar = static_cast<double*>(data);
          for (size_t j=i; j < num_entity*comp_count; j+=comp_count) {
            rvar[j] = temp[k++];
          }
        } else {
          std::ostringstream errmsg;
          errmsg << "IOSS_ERROR: Field storage type must be either integer or double.\n"
		 << "       Field '" << field.get_name() << "' is invalid.\n";
          IOSS_ERROR(errmsg);
        }
        assert(k == num_entity);
      }
      return num_entity;
    }

    int64_t DatabaseIO::read_ss_transient_field(const Ioss::Field& field,
                                                int64_t id, void *variables,
                                                Ioss::IntVector &is_valid_side) const
    {
      size_t num_valid_sides = 0;
      const Ioss::VariableType *var_type = field.raw_storage();
      size_t my_side_count = is_valid_side.size();
      std::vector<double> temp(my_side_count);

      size_t step = get_current_state();

      // get number of components, cycle through each component
      // and add suffix to base 'field_name'.  Look up index
      // of this name in 'nodeVariables' map
      size_t comp_count = var_type->component_count();
      size_t var_index=0;

      char field_suffix_separator = get_field_separator();
      for (size_t i=0; i < comp_count; i++) {
        std::string var_name = var_type->label_name(field.get_name(), i+1, field_suffix_separator);

        // Read the variable...
        int ierr = 0;
        var_index = m_variables[EX_SIDE_SET].find(var_name)->second;
        assert(var_index > 0);
        ierr = decomp->get_var(get_file_pointer(), step, EX_SIDE_SET,
			       var_index, id, my_side_count, temp);
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        // Transfer to 'variables' array.
        size_t j = i;
        if (field.get_type() == Ioss::Field::INTEGER) {
          int *ivar = static_cast<int*>(variables);
          for (size_t k = 0; k < my_side_count; k++) {
            if (is_valid_side[k] == 1) {
              ivar[j] = static_cast<int>(temp[k]);
              j += comp_count;
            }
          }
        } else if (field.get_type() == Ioss::Field::INT64) { // FIX 64 UNSAFE
          int64_t *ivar = static_cast<int64_t*>(variables);
          for (size_t k = 0; k < my_side_count; k++) {
            if (is_valid_side[k] == 1) {
              ivar[j] = static_cast<int64_t>(temp[k]);
              j += comp_count;
            }
          }
        } else if (field.get_type() == Ioss::Field::REAL) {
          double *rvar = static_cast<double*>(variables);
          for (size_t k = 0; k < my_side_count; k++) {
            if (is_valid_side[k] == 1) {
              rvar[j] = temp[k];
              j += comp_count;
            }
          }
        } else {
          std::ostringstream errmsg;
          errmsg << "IOSS_ERROR: Field storage type must be either integer or double.\n"
		 << "       Field '" << field.get_name() << "' is invalid.\n";
          IOSS_ERROR(errmsg);
        }
        if (i+1 == comp_count)
          num_valid_sides =  j / comp_count;
      }
      return num_valid_sides;
    }

    int64_t DatabaseIO::get_side_connectivity(const Ioss::SideBlock* fb,
                                              int64_t id, int64_t,
                                              void *fconnect,
                                              bool map_ids) const
    {
      // Get size of data stored on the file...
      ex_set set_param[1];
      set_param[0].id = id;
      set_param[0].type = EX_SIDE_SET;
      set_param[0].entry_list = NULL;
      set_param[0].extra_list = NULL;
      set_param[0].distribution_factor_list = NULL;
      int ierr = ex_get_sets(get_file_pointer(), 1, set_param);
      if (ierr < 0)
        Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

      int64_t number_sides = set_param[0].num_entry;

      // Allocate space for element and local side number
      assert(number_sides > 0);

      // Need a larger vector to get the entire sideset and then filter down to the correct size...
      std::vector<char> element(number_sides * int_byte_size_api());
      std::vector<char> side(number_sides * int_byte_size_api());
      if (int_byte_size_api() == 4) {
        Ioss::Field elem_field("ids",   Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, number_sides);
        Ioss::Field side_field("sides", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, number_sides);
        decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, (int*)TOPTR(element));
        decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, (int*)TOPTR(side));
      } else {
        Ioss::Field elem_field("ids",   Ioss::Field::INT64, "scalar", Ioss::Field::MESH, number_sides);
        Ioss::Field side_field("sides", Ioss::Field::INT64, "scalar", Ioss::Field::MESH, number_sides);
        decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, (int64_t*)TOPTR(element));
        decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, (int64_t*)TOPTR(side));
      }

      Ioss::IntVector is_valid_side;
      Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, int_byte_size_api(),
                                                  (void*)TOPTR(element), (void*)TOPTR(side),
                                                  number_sides, get_region());

      std::vector<char> elconnect;
      int64_t elconsize = 0; // Size of currently allocated connectivity block
      Ioss::ElementBlock *conn_block = NULL; // Block that we currently
      // have connectivity for

      Ioss::ElementBlock *block = NULL;

      int     *element32 = NULL;
      int64_t *element64 = NULL;
      int     *side32 = NULL;
      int64_t *side64 = NULL;

      int     *elconn32 = NULL;
      int64_t *elconn64 = NULL;
      int     *fconn32 = NULL;
      int64_t *fconn64 = NULL;

      if (int_byte_size_api() == 4) {
        element32 = (int*)TOPTR(element);
        side32 = (int*)TOPTR(side);
        fconn32 = (int*)fconnect;
      } else {
        element64 = (int64_t*)TOPTR(element);
        side64 = (int64_t*)TOPTR(side);
        fconn64 = (int64_t*)fconnect;
      }

      Ioss::IntVector side_elem_map; // Maps the side into the elements
      // connectivity array
      int64_t current_side = -1;
      int nelnode = 0;
      int nfnodes = 0;
      int ieb = 0;
      size_t offset = 0;
      for (ssize_t iel = 0; iel < number_sides; iel++) {
        if (is_valid_side[iel] == 1) {

          int64_t elem_id = 0;
          if (int_byte_size_api() == 4) {
            elem_id = element32[iel];
          } else {
            elem_id = element64[iel];
          }

          // ensure we have correct connectivity
          block = get_region()->get_element_block(elem_id);
          if (conn_block != block) {
            ssize_t nelem = block->get_property("entity_count").get_int();
            nelnode = block->topology()->number_nodes();
            // Used to map element number into position in connectivity array.
            // E.g., element 97 is the (97-offset)th element in this block and
            // is stored in array index (97-offset-1).
            offset = block->get_offset() + 1;
            if (elconsize < nelem * nelnode) {
              elconsize = nelem * nelnode;
              elconnect.resize(elconsize*int_byte_size_api());
              if (int_byte_size_api() == 4) {
                elconn32 = (int*)TOPTR(elconnect);
              } else {
                elconn64 = (int64_t*)TOPTR(elconnect);
              }
            }
            if (map_ids) {
              get_field_internal(block, block->get_field("connectivity"),
				 TOPTR(elconnect), nelem*nelnode*int_byte_size_api());
            } else {
              get_field_internal(block, block->get_field("connectivity_raw"),
				 TOPTR(elconnect), nelem*nelnode*int_byte_size_api());
            }
            conn_block = block;
            current_side = -1;
          }

          // NOTE: Element connectivity is returned with nodes in global id space if "map_ids" false,
          //       otherwise it is in local space.
          int64_t side_id = 0;
          if (int_byte_size_api() == 4) {
            side_id = side32[iel];
          } else {
            side_id = side64[iel];
          }

          if (current_side != side_id) {
            side_elem_map = block->topology()->boundary_connectivity(side_id);
            current_side = side_id;
            nfnodes = block->topology()->boundary_type(side_id)->number_nodes();
          }
          for (int inode = 0; inode < nfnodes; inode++) {
            size_t index = (elem_id-offset)*nelnode + side_elem_map[inode];
            if (int_byte_size_api() == 4) {
              fconn32[ieb++] = elconn32[index];
            } else {
              fconn64[ieb++] = elconn64[index];
            }
          }
        }
      }
      return ierr;
    }

    // Get distribution factors for the specified side block
    int64_t DatabaseIO::get_side_distributions(const Ioss::SideBlock* fb,
                                               int64_t id, int64_t my_side_count,
                                               double *dist_fact,
                                               size_t /* data_size */) const
    {
      // Allocate space for elements and local side numbers
      // Get size of data stored on the file...

      const SetDecompositionData &set = decomp->get_decomp_set(EX_SIDE_SET, id);
      int64_t number_sides = set.ioss_count();
      int64_t number_distribution_factors = set.df_count();

      const Ioss::ElementTopology *ftopo = fb->topology();
      int nfnodes = ftopo->number_nodes();

      if (set.distributionFactorConstant) {
        // Fill in the array with the constant value...
        for (int64_t i=0; i < nfnodes * my_side_count; i++)
          dist_fact[i] = set.distributionFactorValue;
        return 0;
      }

      // Take care of the easy situation -- If 'side_count' ==
      // 'number_sides' then the sideset is stored in a single sideblock
      // and all distribution factors on the database are transferred
      // 1-to-1 into 'dist_fact' array.
      int64_t entity_count = fb->get_property("entity_count").get_int();
      if (fb->owner()->block_count() == 1 && number_sides == entity_count) {
        assert(number_sides == 0 || number_distribution_factors % number_sides == 0);
        assert(number_sides == 0 || number_distribution_factors / number_sides == nfnodes);
        std::string storage = "Real["+Ioss::Utils::to_string(nfnodes)+"]";
        Ioss::Field dist("distribution_factors", Ioss::Field::REAL, storage, Ioss::Field::MESH, number_sides);
        decomp->get_set_mesh_double(get_file_pointer(), EX_SIDE_SET, id, dist, dist_fact);
        return 0;
      }

      std::string storage = "Real["+Ioss::Utils::to_string(nfnodes)+"]";
      Ioss::Field field("distribution_factors", Ioss::Field::REAL, storage, Ioss::Field::MESH, number_distribution_factors/nfnodes);
      std::vector<double> dist(number_distribution_factors);
      decomp->get_set_mesh_double(get_file_pointer(), EX_SIDE_SET, id, field, TOPTR(dist));

      // Another easy situation (and common for exodusII) is if the input
      // distribution factors are all the same value (typically 1).  In
      // that case, we only have to fill in the output array with that
      // value.
      {
        double value = number_distribution_factors > 0 ? dist[0] : 0.0;
        bool constant = true;
        for (int64_t i=1; i < number_distribution_factors; i++) {
          if (dist[i] != value) {
            constant = false;
            break;
          }
        }

        constant = util().global_minmax(constant ? 1 : 0, Ioss::ParallelUtils::DO_MIN);

        if (constant) {
          if (value == 0.0)
            value = 1.0;  // Take care of some buggy mesh generators
          for (ssize_t j=0; j < my_side_count * nfnodes; j++)
            dist_fact[j] = value;
          return 0;
        }
      }

      // If we get to here, the underlying sideset contains multiple side
      // topologies and the distribution factors are non-constant. Need to
      // allocate space to store all distribution factors and then pull
      // out those that are applied to sides with the correct topology.

      // Allocate space for element and local side number (this is bulk
      // data...)
      //----
      std::vector<char> element(number_sides * int_byte_size_api());
      std::vector<char> sides(number_sides * int_byte_size_api());
      if (int_byte_size_api() == 4) {
        Ioss::Field elem_field("ids",   Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, number_sides);
        Ioss::Field side_field("sides", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, number_sides);
        decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, (int*)TOPTR(element));
        decomp32->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, (int*)TOPTR(sides));
      } else {
        Ioss::Field elem_field("ids",   Ioss::Field::INT64, "scalar", Ioss::Field::MESH, number_sides);
        Ioss::Field side_field("sides", Ioss::Field::INT64, "scalar", Ioss::Field::MESH, number_sides);
        decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, elem_field, (int64_t*)TOPTR(element));
        decomp64->get_set_mesh_var(get_file_pointer(), EX_SIDE_SET, id, side_field, (int64_t*)TOPTR(sides));
      }
      //----

      Ioss::IntVector is_valid_side;
      Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, int_byte_size_api(),
                                                  TOPTR(element), TOPTR(sides),
                                                  number_sides, get_region());

      int64_t ieb = 0; // counter for distribution factors in this sideblock
      int64_t idb = 0; // counter for distribution factors read from database
      Ioss::ElementBlock *block = NULL;

      int     *element32 = NULL;
      int64_t *element64 = NULL;
      int     *side32 = NULL;
      int64_t *side64 = NULL;

      if (int_byte_size_api() == 4) {
        element32 = (int*)TOPTR(element);
        side32 = (int*)TOPTR(sides);
      } else {
        element64 = (int64_t*)TOPTR(element);
        side64 = (int64_t*)TOPTR(sides);
      }

      for (int64_t iel = 0; iel < number_sides; iel++) {
        int64_t elem_id = 0;
        int64_t side_id = 0;
        if (int_byte_size_api() == 4) {
          elem_id = element32[iel];
          side_id = side32[iel];
        } else {
          elem_id = element64[iel];
          side_id = side64[iel];
        }

        if (block == NULL || !block->contains(elem_id)) {
          block = get_region()->get_element_block(elem_id);
        }

        if (block == NULL) {
          std::ostringstream errmsg;
          errmsg << "INTERNAL ERROR: Could not find element block containing element with id " << elem_id
		 << "Something is wrong in the Iopx::DatabaseIO class. Please report.\n";
          IOSS_ERROR(errmsg);
        }

        const Ioss::ElementTopology *topo = block->topology()->boundary_type(side_id);

        if (topo == NULL) {
          std::ostringstream errmsg;
          errmsg << "INTERNAL ERROR: Could not find topology of element block boundary. "
		 << "Something is wrong in the Iopx::DatabaseIO class. Please report.\n";
          IOSS_ERROR(errmsg);
        }

        int nside_nodes = topo->number_nodes();

        if (is_valid_side[iel] == 1) {
          // This side belongs in the sideblock
          for (int64_t i=0; i < nside_nodes; i++) {
            dist_fact[ieb++] = dist[idb++];
          }
        } else {
          // Skip over unused 'dist' factors
          idb += topo->number_nodes();
        }
      }

      assert(ieb == my_side_count * nfnodes);
      // If the following assert fails, it may be due to bug in Patran
      // which writes too many distribution factors to the database in a
      // mixed element case. Note that this is checked earlier also with a
      // better error message.
      assert(idb == number_distribution_factors);
      return 0;
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock* nb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);

      size_t proc_offset = 0;
      if (nb->property_exists("processor_offset"))
	proc_offset = nb->get_property("processor_offset").get_int();
      size_t file_count = num_to_get;
      if (nb->property_exists("locally_owned_count"))
	file_count = nb->get_property("locally_owned_count").get_int();

      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
	if (field.get_name() == "owning_processor") {
	  // Set the nodeOwningProcessor vector for all nodes on this processor.
	  // Value is the processor that owns the node.
	  nodeOwningProcessor.reserve(num_to_get);
	  if (int_byte_size_api() == 4) {
	    int *owned = (int *)data;
	    for (size_t i=0; i < num_to_get; i++) {
	      nodeOwningProcessor.push_back(owned[i]);
	    }
	  } else {
	    int64_t *owned = (int64_t *)data;
	    for (size_t i=0; i < num_to_get; i++) {
	      nodeOwningProcessor.push_back(owned[i]);
	    }
	  }

	  // Now create the "implicit local" to "implicit global"
	  // map which maps data from its local implicit position
	  // to its implicit (1..num_global_node) position in the
	  // global file.  This is needed for the global-to-local
	  // mapping of element connectivity and nodeset nodelists.
	  create_implicit_global_map();
	}

	else if (field.get_name() == "mesh_model_coordinates_x") {
	  double *rdata = static_cast<double*>(data);
	  std::vector<double> file_data; file_data.reserve(file_count);
	  map_data(nodeOwningProcessor, myProcessor, rdata, file_data);

	  int ierr = ex_put_partial_coord(get_file_pointer(), proc_offset+1, file_count, rdata, NULL, NULL);
	  if (ierr < 0)
	    Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}

	else if (field.get_name() == "mesh_model_coordinates_y") {
	  double *rdata = static_cast<double*>(data);
	  std::vector<double> file_data; file_data.reserve(file_count);
	  map_data(nodeOwningProcessor, myProcessor, rdata, file_data);
	  int ierr = ex_put_partial_coord(get_file_pointer(), proc_offset+1, file_count, NULL, rdata, NULL);
	  if (ierr < 0)
	    Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}

	else if (field.get_name() == "mesh_model_coordinates_z") {
	  double *rdata = static_cast<double*>(data);
	  std::vector<double> file_data; file_data.reserve(file_count);
	  map_data(nodeOwningProcessor, myProcessor, rdata, file_data);
	  int ierr = ex_put_partial_coord(get_file_pointer(), proc_offset+1, file_count, NULL, NULL, rdata);
	  if (ierr < 0)
	    Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}

	else if (field.get_name() == "mesh_model_coordinates") {
	  // Data required by upper classes store x0, y0, z0, ... xn, yn, zn
	  // Data stored in exodusII file is x0, ..., xn, y0, ..., yn, z0, ..., zn
	  // so we have to allocate some scratch memory to read in the data
	  // and then map into supplied 'data'
	  std::vector<double> x;
	  std::vector<double> y;
	  std::vector<double> z;

	  x.reserve(num_to_get);
	  if (spatialDimension > 1)
	    y.reserve(num_to_get);
	  if (spatialDimension == 3)
	    z.reserve(num_to_get);

	  // Cast 'data' to correct size -- double
	  double *rdata = static_cast<double*>(data);
	  map_data(nodeOwningProcessor, myProcessor, rdata, x, 0, spatialDimension);
	  if (spatialDimension > 1)
	    map_data(nodeOwningProcessor, myProcessor, rdata, y, 1, spatialDimension);
	  if (spatialDimension == 3)
	    map_data(nodeOwningProcessor, myProcessor, rdata, z, 2, spatialDimension);

	  int ierr = ex_put_partial_coord(get_file_pointer(), proc_offset+1, file_count,
					  TOPTR(x), TOPTR(y), TOPTR(z));
	  if (ierr < 0)
	    Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

	} else if (field.get_name() == "ids") {
	  // The ids coming in are the global ids; their position is the
	  // local id -1 (That is, data[0] contains the global id of local
	  // node 1)

	  // Another 'const-cast' since we are modifying the database just
	  // for efficiency; which the client does not see...
	  handle_node_ids(data, num_to_get, proc_offset, file_count);

	} else if (field.get_name() == "connectivity") {
	  // Do nothing, just handles an idiosyncracy of the GroupingEntity
	} else if (field.get_name() == "connectivity_raw") {
	  // Do nothing, just handles an idiosyncracy of the GroupingEntity
	} else if (field.get_name() == "node_connectivity_status") {
	  // Do nothing, input only field.
	} else if (field.get_name() == "implicit_ids") {
	  // Do nothing, input only field.
	} else {
	  return Ioss::Utils::field_warning(nb, field, "mesh output");
	}

      } else if (role == Ioss::Field::TRANSIENT) {
	// Check if the specified field exists on this node block.
	// Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	// exist on the database as scalars with the appropriate
	// extensions.

	// Transfer each component of the variable into 'data' and then
	// output.  Need temporary storage area of size 'number of
	// nodes in this block.
	write_nodal_transient_field(EX_NODE_BLOCK, field, nb, num_to_get, data);

      } else if (role == Ioss::Field::REDUCTION) {
	store_reduction_field(EX_NODE_BLOCK, field, nb, data);
      }
      return num_to_get;
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);

      int ierr = 0;

      // Get the element block id and element count
      int64_t id = Ioex::get_id(eb, EX_ELEM_BLOCK, &ids_);
      int64_t my_element_count = eb->get_property("entity_count").get_int();
      Ioss::Field::RoleType role = field.get_role();

      size_t proc_offset = 0;
      if (eb->property_exists("processor_offset"))
        proc_offset = eb->get_property("processor_offset").get_int();
      size_t file_count = num_to_get;
      if (eb->property_exists("locally_owned_count"))
        file_count = eb->get_property("locally_owned_count").get_int();

      if (role == Ioss::Field::MESH) {
        // Handle the MESH fields required for an ExodusII file model.
        // (The 'genesis' portion)
        if (field.get_name() == "connectivity") {
          // Map element connectivity from global node id to local node id.
          int element_nodes = eb->get_property("topology_node_count").get_int();

          // Maps global to local
          nodeMap.reverse_map_data(data, field, num_to_get*element_nodes);

          // Maps local to "global_implicit"
          if (int_byte_size_api() == 4) {
            map_local_to_global_implicit((int*)data, num_to_get*element_nodes, nodeGlobalImplicitMap);
          } else {
            map_local_to_global_implicit((int64_t*)data, num_to_get*element_nodes, nodeGlobalImplicitMap);
          }

          ierr = ex_put_partial_elem_conn(get_file_pointer(), id, proc_offset+1, file_count, data);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        } else if (field.get_name() == "connectivity_edge") {
          // Map element connectivity from global edge id to local edge id.
          int element_edges = field.transformed_storage()->component_count();
          edgeMap.reverse_map_data(data, field, num_to_get*element_edges);
          ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, NULL, data, NULL);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        } else if (field.get_name() == "connectivity_face") {
          // Map element connectivity from global face id to local face id.
          int element_faces = field.transformed_storage()->component_count();
          faceMap.reverse_map_data(data, field, num_to_get*element_faces);
          ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, NULL, NULL, data);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        } else if (field.get_name() == "connectivity_raw") {
          // Element connectivity is already in local node id, map local to "global_implicit"
          int element_nodes = eb->get_property("topology_node_count").get_int();
          if (int_byte_size_api() == 4) {
            map_local_to_global_implicit((int*)data, num_to_get*element_nodes, nodeGlobalImplicitMap);
          } else {
            map_local_to_global_implicit((int64_t*)data, num_to_get*element_nodes, nodeGlobalImplicitMap);
          }

          ierr = ex_put_partial_elem_conn(get_file_pointer(), id, proc_offset+1, file_count, data);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        } else if (field.get_name() == "ids") {
          size_t glob_map_offset = eb->get_property("global_map_offset").get_int();
          handle_element_ids(eb, data, num_to_get, glob_map_offset+proc_offset, file_count);
        } else if (field.get_name() == "implicit_ids") {
          // Do nothing, input only field.
        } else if (field.get_name() == "skin") {
          // This is (currently) for the skinned body. It maps the
          // side element on the skin to the original element/local
          // side number.  It is a two component field, the first
          // component is the global id of the underlying element in
          // the initial mesh and its local side number (1-based).

          // FIX: Hardwired map ids....
          int map_count = ex_inquire_int(get_file_pointer(), EX_INQ_ELEM_MAP);
          if (map_count == 0) {
            // This needs to be fixed... Currently hardwired....
            ierr = ex_put_map_param(get_file_pointer(), 0, 2);
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }

          std::vector<char> element(my_element_count * int_byte_size_api());
          std::vector<char> side(my_element_count * int_byte_size_api());

          if (int_byte_size_api() == 4) {
            int *el_side = (int *)data;
            int *element32 = (int*)TOPTR(element);
            int *side32 = (int*)TOPTR(side);

            int index = 0;
            for (int i=0; i < my_element_count; i++) {
              element32[i] = el_side[index++];
              side32[i]    = el_side[index++];
            }
          } else {
            int64_t *el_side = (int64_t *)data;
            int64_t *element64 = (int64_t*)TOPTR(element);
            int64_t *side64 = (int64_t*)TOPTR(side);

            int64_t index = 0;
            for (int64_t i=0; i < my_element_count; i++) {
              element64[i] = el_side[index++];
              side64[i]    = el_side[index++];
            }
          }

          size_t eb_offset = eb->get_offset() + proc_offset;
          ierr = ex_put_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 1, eb_offset+1, file_count,
					TOPTR(element));
          if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

          ierr = ex_put_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 2, eb_offset+1, file_count,
					TOPTR(side));
          if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

          if (map_count == 0) {
            // NOTE: ex_put_*num_map must be called prior to defining the name...
            ierr = ex_put_name(get_file_pointer(), EX_ELEM_MAP, 1, "skin:parent_element_id");
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

            ierr = ex_put_name(get_file_pointer(), EX_ELEM_MAP, 2, "skin:parent_element_side_number");
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }

        } else {
          num_to_get = Ioss::Utils::field_warning(eb, field, "mesh output");
        }

      } else if (role == Ioss::Field::ATTRIBUTE) {
        num_to_get = write_attribute_field(EX_ELEM_BLOCK, field, eb, data);

      } else if (role == Ioss::Field::TRANSIENT) {
        // Check if the specified field exists on this element block.
        // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
        // exist on the database as scalars with the appropriate
        // extensions.

        // Transfer each component of the variable into 'data' and then
        // output.  Need temporary storage area of size 'number of
        // elements in this block.
        write_entity_transient_field(EX_ELEM_BLOCK, field, eb, my_element_count, data);

      } else if (role == Ioss::Field::REDUCTION) {
        store_reduction_field(EX_ELEM_BLOCK, field, eb, data);
      }
      return num_to_get;
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);

      int ierr = 0;

      // Get the face block id and face count
      int64_t id = Ioex::get_id(eb, EX_FACE_BLOCK, &ids_);
      int64_t my_face_count = eb->get_property("entity_count").get_int();
      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
	// Handle the MESH fields required for an ExodusII file model.
	// (The 'genesis' portion)
	if (field.get_name() == "connectivity") {
	  if (my_face_count > 0) {
	    // Map face connectivity from global node id to local node id.
	    int face_nodes = eb->get_property("topology_node_count").get_int();
	    nodeMap.reverse_map_data(data, field, num_to_get*face_nodes);
	    ierr = ex_put_conn(get_file_pointer(), EX_FACE_BLOCK, id, data, NULL, NULL);
	    if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }
	} else if (field.get_name() == "connectivity_edge") {
	  if (my_face_count > 0) {
	    // Map face connectivity from global edge id to local edge id.
	    // Do it in 'data' ...
	    int face_edges = field.transformed_storage()->component_count();
	    edgeMap.reverse_map_data(data, field, num_to_get*face_edges);
	    ierr = ex_put_conn(get_file_pointer(), EX_FACE_BLOCK, id, NULL, data, NULL);
	    if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }
	} else if (field.get_name() == "connectivity_raw") {
	  // Do nothing, input only field.
	} else if (field.get_name() == "ids") {
	  handle_face_ids(eb, data, num_to_get);

	} else {
	  num_to_get = Ioss::Utils::field_warning(eb, field, "mesh output");
	}

      } else if (role == Ioss::Field::ATTRIBUTE) {
	num_to_get = write_attribute_field(EX_FACE_BLOCK, field, eb, data);

      } else if (role == Ioss::Field::TRANSIENT) {
	// Check if the specified field exists on this face block.
	// Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	// exist on the database as scalars with the appropriate
	// extensions.

	// Transfer each component of the variable into 'data' and then
	// output.  Need temporary storage area of size 'number of
	// faces in this block.
	write_entity_transient_field(EX_FACE_BLOCK, field, eb, my_face_count, data);

      } else if (role == Ioss::Field::REDUCTION) {
	store_reduction_field(EX_FACE_BLOCK, field, eb, data);
      }
      return num_to_get;
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);

      int ierr = 0;

      // Get the edge block id and edge count
      int64_t id = Ioex::get_id(eb, EX_EDGE_BLOCK, &ids_);
      int64_t my_edge_count = eb->get_property("entity_count").get_int();
      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
	// Handle the MESH fields required for an ExodusII file model. (The 'genesis' portion)
	if (field.get_name() == "connectivity") {
	  if (my_edge_count > 0) {
	    // Map edge connectivity from global node id to local node id.
	    int edge_nodes = eb->get_property("topology_node_count").get_int();
	    nodeMap.reverse_map_data(data, field, num_to_get*edge_nodes);
	    ierr = ex_put_conn(get_file_pointer(), EX_EDGE_BLOCK, id, data, NULL, NULL);
	    if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }
	} else if (field.get_name() == "connectivity_raw") {
	  // Do nothing, input only field.
	} else if (field.get_name() == "ids") {
	  handle_edge_ids(eb, data, num_to_get);

	} else {
	  num_to_get = Ioss::Utils::field_warning(eb, field, "mesh output");
	}

      } else if (role == Ioss::Field::ATTRIBUTE) {
	num_to_get = write_attribute_field(EX_EDGE_BLOCK, field, eb, data);

      } else if (role == Ioss::Field::TRANSIENT) {
	// Check if the specified field exists on this edge block.
	// Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	// exist on the database as scalars with the appropriate
	// extensions.

	// Transfer each component of the variable into 'data' and then
	// output.  Need temporary storage area of size 'number of
	// edges in this block.
	write_entity_transient_field(EX_EDGE_BLOCK, field, eb, my_edge_count, data);

      } else if (role == Ioss::Field::REDUCTION) {
	store_reduction_field(EX_EDGE_BLOCK, field, eb, data);
      }
      return num_to_get;
    }

    int64_t DatabaseIO::handle_node_ids(void* ids, int64_t num_to_get, size_t offset, size_t count) const
    {
      /*!
       * There are two modes we need to support in this routine:
       * 1. Initial definition of node map (local->global) and
       * nodeMap.reverse (global->local).
       * 2. Redefinition of node map via 'reordering' of the original
       * map when the nodes on this processor are the same, but their
       * order is changed (or count because of ghosting)
       *
       * So, there will be two maps the 'nodeMap.map' map is a 'direct lookup'
       * map which maps current local position to global id and the
       * 'nodeMap.reverse' is an associative lookup which maps the
       * global id to 'original local'.  There is also a
       * 'nodeMap.reorder' which is direct lookup and maps current local
       * position to original local.

       * The ids coming in are the global ids; their position is the
       * "local id-1" (That is, data[0] contains the global id of local
       * node 1 in this node block).
       *
       * int local_position = nodeMap.reverse[NodeMap[i+1]]
       * (the nodeMap.map and nodeMap.reverse are 1-based)
       *
       * To determine which map to update on a call to this function, we
       * use the following hueristics:
       * -- If the database state is 'STATE_MODEL:', then update the
       *    'nodeMap.reverse' and 'nodeMap.map'
       *
       * -- If the database state is not STATE_MODEL, then leave the
       *    'nodeMap.reverse' and 'nodeMap.map' alone since they correspond to the
       *    information already written to the database. [May want to add a
       *    STATE_REDEFINE_MODEL]
       *
       * -- In both cases, update the nodeMap.reorder
       *
       * NOTE: The mapping is done on TRANSIENT fields only; MODEL fields
       *       should be in the orginal order...
       */
      if (dbState == Ioss::STATE_MODEL || dbState == Ioss::STATE_DEFINE_MODEL) {
        if (!nodeMap.defined) {
          if (nodeMap.map.empty()) {
            nodeMap.map.resize(num_to_get+1);
            nodeMap.map[0] = -1;
          }

          if (nodeMap.map[0] == -1) {
            if (int_byte_size_api() == 4) {
              nodeMap.set_map(static_cast<int*>(ids), num_to_get, 0);
            } else {
              nodeMap.set_map(static_cast<int64_t*>(ids), num_to_get, 0);
            }
          }

          nodeMap.build_reverse_map(myProcessor);
          nodeMap.build_reorder_map(0, num_to_get);

          // Only a single nodeblock and all set
          assert(nodeMap.map[0] == -1 || nodeMap.reverse.size() == (size_t)num_to_get);
          assert(get_region()->get_property("node_block_count").get_int() == 1);
          nodeMap.defined = true;
        }
      }
      return num_to_get;
    }

    namespace {
      size_t handle_block_ids(const Ioss::EntityBlock *eb,
                              ex_entity_type map_type,
                              Ioss::State db_state,
                              Ioss::Map &entity_map,
                              void* ids, size_t int_byte_size, size_t num_to_get,
                              size_t offset, size_t count, int file_pointer, int my_processor)
      {
        /*!
         * NOTE: "element" is generic for "element", "face", or "edge"
         *
         * There are two modes we need to support in this routine:
         * 1. Initial definition of element map (local->global) and
         * elemMap.reverse (global->local).
         * 2. Redefinition of element map via 'reordering' of the original
         * map when the elements on this processor are the same, but their
         * order is changed.
         *
         * So, there will be two maps the 'elemMap.map' map is a 'direct lookup'
         * map which maps current local position to global id and the
         * 'elemMap.reverse' is an associative lookup which maps the
         * global id to 'original local'.  There is also a
         * 'elemMap.reorder' which is direct lookup and maps current local
         * position to original local.

         * The ids coming in are the global ids; their position is the
         * local id -1 (That is, data[0] contains the global id of local
         * element 1 in this element block).  The 'model-local' id is
         * given by eb_offset + 1 + position:
         *
         * int local_position = elemMap.reverse[ElementMap[i+1]]
         * (the elemMap.map and elemMap.reverse are 1-based)
         *
         * But, this assumes 1..numel elements are being output at the same
         * time; we are actually outputting a blocks worth of elements at a
         * time, so we need to consider the block offsets.
         * So... local-in-block position 'i' is index 'eb_offset+i' in
         * 'elemMap.map' and the 'local_position' within the element
         * blocks data arrays is 'local_position-eb_offset'.  With this, the
         * position within the data array of this element block is:
         *
         * int eb_position =
         * elemMap.reverse[elemMap.map[eb_offset+i+1]]-eb_offset-1
         *
         * To determine which map to update on a call to this function, we
         * use the following hueristics:
         * -- If the database state is 'Ioss::STATE_MODEL:', then update the
         *    'elemMap.reverse'.
         * -- If the database state is not Ioss::STATE_MODEL, then leave
         *    the 'elemMap.reverse' alone since it corresponds to the
         *    information already written to the database. [May want to add
         *    a Ioss::STATE_REDEFINE_MODEL]
         * -- Always update elemMap.map to match the passed in 'ids'
         *    array.
         *
         * NOTE: the maps are built an element block at a time...
         * NOTE: The mapping is done on TRANSIENT fields only; MODEL fields
         *       should be in the orginal order...
         */

        // Overwrite this portion of the 'elemMap.map', but keep other
        // parts as they were.  We are adding elements starting at position
        // 'eb_offset+offset' and ending at
        // 'eb_offset+offset+num_to_get'. If the entire block is being
        // processed, this reduces to the range 'eb_offset..eb_offset+my_element_count'

        int64_t eb_offset = eb->get_offset();

        if (int_byte_size == 4) {
          entity_map.set_map(static_cast<int*>(ids), num_to_get, eb_offset);
        } else {
          entity_map.set_map(static_cast<int64_t*>(ids), num_to_get, eb_offset);
        }

        // Now, if the state is Ioss::STATE_MODEL, update the reverseEntityMap
        if (db_state == Ioss::STATE_MODEL) {
          entity_map.build_reverse_map(num_to_get, eb_offset, my_processor);

          // Output this portion of the entity number map
          int ierr = ex_put_partial_id_map(file_pointer, map_type, offset+1, count, ids);
          if (ierr < 0)
            Ioex::exodus_error(file_pointer, __LINE__, my_processor);
        }
        // Build the reorderEntityMap which does a direct mapping from
        // the current topologies local order to the local order
        // stored in the database...  This is 0-based and used for
        // remapping output and input TRANSIENT fields.
        entity_map.build_reorder_map(eb_offset, num_to_get);
        return num_to_get;
      }
    }

    int64_t DatabaseIO::handle_element_ids(const Ioss::ElementBlock *eb, void* ids, size_t num_to_get,
                                           size_t offset, size_t count) const
    {
      if (dbState == Ioss::STATE_MODEL) {
        if (elemGlobalImplicitMap.empty()) {
          elemGlobalImplicitMap.resize(elementCount);
        }
        // Build the implicit_global map used to map an elements
        // local-implicit position to the global-implicit
        // position. Primarily used for sideset elements.  'count'
        // Elements starting at 'eb_offset' map to the global implicit
        // position of 'offset'
        int64_t eb_offset = eb->get_offset();
        for (size_t i=0; i < count; i++) {
          elemGlobalImplicitMap[eb_offset+i] = offset+i+1;
        }
      }

      if (elemMap.map.empty()) {
        elemMap.map.resize(elementCount+1);
        elemMap.map[0] = -1;
      }
      return handle_block_ids(eb, EX_ELEM_MAP, dbState, elemMap,
			      ids, int_byte_size_api(), num_to_get, offset, count,
			      get_file_pointer(), myProcessor);
    }

    int64_t DatabaseIO::handle_face_ids(const Ioss::FaceBlock *eb, void* ids, size_t num_to_get) const
    {
      if (faceMap.map.empty()) {
        faceMap.map.resize(faceCount+1);
        faceMap.map[0] = -1;
      }
      return handle_block_ids(eb, EX_FACE_MAP, dbState, faceMap,
			      ids, int_byte_size_api(), num_to_get, 0, 0, get_file_pointer(), myProcessor);
    }

    int64_t DatabaseIO::handle_edge_ids(const Ioss::EdgeBlock *eb, void* ids, size_t num_to_get) const
    {
      if (edgeMap.map.empty()) {
        edgeMap.map.resize(edgeCount+1);
        edgeMap.map[0] = -1;
      }
      return handle_block_ids(eb, EX_EDGE_MAP, dbState, edgeMap,
			      ids, int_byte_size_api(), num_to_get, 0, 0, get_file_pointer(), myProcessor);
    }

    void DatabaseIO::write_nodal_transient_field(ex_entity_type /* type */,
                                                 const Ioss::Field &field,
                                                 const Ioss::NodeBlock *nb,
                                                 int64_t count, void *variables) const
    {
      Ioss::Field::BasicType ioss_type = field.get_type();
      assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER ||
             ioss_type == Ioss::Field::INT64 ||ioss_type == Ioss::Field::COMPLEX);

      // Note that if the field's basic type is COMPLEX, then each component of
      // the VariableType is a complex variable consisting of a real and
      // imaginary part.  Since exodus cannot handle complex variables,
      // we have to output a (real and imaginary) X (number of
      // components) fields. For example, if V is a 3d vector of complex
      // data, the data in the 'variables' array are v_x, v.im_x, v_y,
      // v.im_y, v_z, v.im_z which need to be output in six separate
      // exodus fields.  These fields were already defined in
      // "write_results_metadata".

      const Ioss::VariableType *var_type = field.transformed_storage();
      std::vector<double> temp(count);

      int step = get_current_state();
      step = get_database_step(step);

      // get number of components, cycle through each component
      // and add suffix to base 'field_name'.  Look up index
      // of this name in 'm_variables[EX_NODE_BLOCK]' map
      int comp_count = var_type->component_count();
      int var_index=0;

      int re_im = 1;
      if (ioss_type == Ioss::Field::COMPLEX)
        re_im = 2;
      for (int complex_comp = 0; complex_comp < re_im; complex_comp++) {
        std::string field_name = field.get_name();
        if (re_im == 2) {
          field_name += complex_suffix[complex_comp];
        }

        char field_suffix_separator = get_field_separator();
        for (int i=0; i < comp_count; i++) {
          std::string var_name = var_type->label_name(field_name, i+1, field_suffix_separator);

          if (m_variables[EX_NODE_BLOCK].find(var_name) == m_variables[EX_NODE_BLOCK].end()) {
            std::ostringstream errmsg;
            errmsg << "ERROR: Could not find nodal variable '" << var_name << "'\n";
            IOSS_ERROR(errmsg);
          }

          var_index = m_variables[EX_NODE_BLOCK].find(var_name)->second;

          size_t begin_offset = (re_im*i)+complex_comp;
          size_t stride = re_im*comp_count;
          size_t num_out = 0;

          if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX)
            num_out = nodeMap.map_field_to_db_scalar_order(static_cast<double*>(variables),
							   temp, begin_offset, count, stride, 0);
          else if (ioss_type == Ioss::Field::INTEGER)
            num_out = nodeMap.map_field_to_db_scalar_order(static_cast<int*>(variables),
							   temp, begin_offset, count, stride, 0);
          else if (ioss_type == Ioss::Field::INT64)
            num_out = nodeMap.map_field_to_db_scalar_order(static_cast<int64_t*>(variables),
							   temp, begin_offset, count, stride, 0);

          if (num_out != (size_t)nodeCount) {
            std::ostringstream errmsg;
            errmsg << "ERROR: Problem outputting nodal variable '" << var_name
		   << "' with index = " << var_index << " to file "
		   << util().decode_filename(get_filename(), isParallel) << "\n"
		   << "Should have output " << nodeCount << " values, but instead only output " << num_out << " values.\n";
            IOSS_ERROR(errmsg);
          }

          // Write the variable...
          size_t proc_offset = 0;
          if (nb->property_exists("processor_offset"))
            proc_offset = nb->get_property("processor_offset").get_int();
          size_t file_count = num_out;
          if (nb->property_exists("locally_owned_count"))
            file_count = nb->get_property("locally_owned_count").get_int();

          std::vector<double> file_temp; file_temp.reserve(file_count);
          map_data(nodeOwningProcessor, myProcessor, TOPTR(temp), file_temp);
          int ierr = ex_put_partial_var(get_file_pointer(), step, EX_NODE_BLOCK, var_index, 0,
					proc_offset+1, file_count, TOPTR(file_temp));
          if (ierr < 0) {
            std::ostringstream errmsg;
            errmsg << "ERROR: Problem outputting nodal variable '" << var_name
		   << "' with index = " << var_index << " to file "
		   << util().decode_filename(get_filename(), isParallel) << "\n";
            IOSS_ERROR(errmsg);
          }
        }
      }
    }

    void DatabaseIO::write_entity_transient_field(ex_entity_type type,
                                                  const Ioss::Field& field,
                                                  const Ioss::GroupingEntity *ge,
                                                  int64_t count,
                                                  void *variables) const
    {
      static Ioss::Map non_element_map; // Used as an empty map for ge->type() != element block.
      const Ioss::VariableType *var_type = field.transformed_storage();
      std::vector<double> temp(count);

      int step = get_current_state();
      step = get_database_step(step);

      Ioss::Map *map = NULL;
      ssize_t eb_offset = 0;
      if (ge->type() == Ioss::ELEMENTBLOCK) {
        const Ioss::ElementBlock *elb = dynamic_cast<const Ioss::ElementBlock*>(ge);
        assert(elb != NULL);
        eb_offset = elb->get_offset();
        map = &elemMap;
      } else {
        map = &non_element_map;
      }

      Ioss::Field::BasicType ioss_type = field.get_type();
      assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER ||
             ioss_type == Ioss::Field::INT64 || ioss_type == Ioss::Field::COMPLEX);

      // Note that if the field's basic type is COMPLEX, then each component of
      // the VariableType is a complex variable consisting of a real and
      // imaginary part.  Since exodus cannot handle complex variables,
      // we have to output a (real and imaginary) X (number of
      // components) fields. For example, if V is a 3d vector of complex
      // data, the data in the 'variables' array are v_x, v.im_x, v_y,
      // v.im_y, v_z, v.im_z which need to be output in six separate
      // exodus fields.  These fields were already defined in
      // "write_results_metadata".


      // get number of components, cycle through each component
      // and add suffix to base 'field_name'.  Look up index
      // of this name in 'm_variables[type]' map
      int comp_count = var_type->component_count();
      int var_index=0;

      int re_im = 1;
      if (ioss_type == Ioss::Field::COMPLEX)
        re_im = 2;
      for (int complex_comp = 0; complex_comp < re_im; complex_comp++) {
        std::string field_name = field.get_name();
        if (re_im == 2) {
          field_name += complex_suffix[complex_comp];
        }

        char field_suffix_separator = get_field_separator();
        for (int i=0; i < comp_count; i++) {
          std::string var_name = var_type->label_name(field_name, i+1, field_suffix_separator);

          var_index = m_variables[type].find(var_name)->second;
          assert(var_index > 0);

          // var is a [count,comp,re_im] array;  re_im = 1(real) or 2(complex)
          // beg_offset = (re_im*i)+complex_comp
          // number_values = count
          // stride = re_im*comp_count
          ssize_t begin_offset = (re_im*i)+complex_comp;
          ssize_t stride = re_im*comp_count;

          if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX)
            map->map_field_to_db_scalar_order(static_cast<double*>(variables),
					      temp, begin_offset, count, stride, eb_offset);
          else if (ioss_type == Ioss::Field::INTEGER)
            map->map_field_to_db_scalar_order(static_cast<int*>(variables),
					      temp, begin_offset, count, stride, eb_offset);
          else if (ioss_type == Ioss::Field::INT64)
            map->map_field_to_db_scalar_order(static_cast<int64_t*>(variables),
					      temp, begin_offset, count, stride, eb_offset);

          // Write the variable...
          size_t proc_offset = 0;
          if (ge->property_exists("processor_offset"))
            proc_offset = ge->get_property("processor_offset").get_int();
          size_t file_count = count;
          if (ge->property_exists("locally_owned_count"))
            file_count = ge->get_property("locally_owned_count").get_int();

          int64_t id = Ioex::get_id(ge, type, &ids_);
          int ierr;
          if (type == EX_SIDE_SET) {
            size_t offset = ge->get_property("set_offset").get_int();
            ierr = ex_put_partial_var(get_file_pointer(), step, type, var_index, id, proc_offset+offset+1,
				      count, TOPTR(temp));
          } else {
            // Write the variable...
            if (type == EX_NODE_SET) {
              std::vector<double> file_data; file_data.reserve(file_count);
              map_nodeset_data(nodesetOwnedNodes[ge], TOPTR(temp), file_data);
              ierr = ex_put_partial_var(get_file_pointer(), step, type, var_index, id,
					proc_offset+1, file_count, TOPTR(file_data));

            } else {
              ierr = ex_put_partial_var(get_file_pointer(), step, type, var_index, id,
					proc_offset+1, file_count, TOPTR(temp));
            }
          }

          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        }
      }
    }

    int64_t DatabaseIO::put_Xset_field_internal(ex_entity_type type,
                                                const Ioss::EntitySet* ns,
                                                const Ioss::Field& field,
                                                void *data, size_t data_size) const
    {
      size_t entity_count = ns->get_property("entity_count").get_int();
      size_t num_to_get = field.verify(data_size);

      int64_t id = Ioex::get_id(ns, type, &ids_);
      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {

        void *out_data = data;
        std::vector<int>     i32data;
        std::vector<int64_t> i64data;
        std::vector<double>  dbldata;

        size_t proc_offset = 0;
        if (ns->property_exists("processor_offset"))
          proc_offset = ns->get_property("processor_offset").get_int();
        size_t file_count = num_to_get;
        if (ns->property_exists("locally_owned_count"))
          file_count = ns->get_property("locally_owned_count").get_int();

        if (field.get_name() == "ids" ||
            field.get_name() == "ids_raw") {
          // Map node id from global node id to local node id.
          // Do it in 'data' ...

          if (field.get_name() == "ids") {
            nodeMap.reverse_map_data(data, field, num_to_get);
          }

          if (type == EX_NODE_SET) {
            nodesetOwnedNodes[ns].reserve(file_count);
            if (int_byte_size_api() == 4) {
              i32data.reserve(file_count);
              map_nodeset_id_data(nodeOwningProcessor, nodesetOwnedNodes[ns], myProcessor,
                                  (int*)data, num_to_get, i32data);
              assert(i32data.size() == file_count);
              // Maps local to "global_implicit"
              map_local_to_global_implicit(TOPTR(i32data), file_count, nodeGlobalImplicitMap);
              out_data = TOPTR(i32data);
            } else {
              i64data.reserve(file_count);
              assert(i64data.size() == file_count);
              map_nodeset_id_data(nodeOwningProcessor, nodesetOwnedNodes[ns], myProcessor,
                                  (int64_t*)data, num_to_get, i64data);
              map_local_to_global_implicit(TOPTR(i64data), file_count, nodeGlobalImplicitMap);
              out_data = TOPTR(i64data);
            }
          }
          int ierr = ex_put_partial_set(get_file_pointer(), type, id, proc_offset+1, file_count, out_data, NULL);
          if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        } else if (field.get_name() == "orientation") {
          int ierr = ex_put_partial_set(get_file_pointer(), type, id, proc_offset+1, file_count, NULL, out_data);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        } else if (field.get_name() == "distribution_factors") {
          int ierr = 0;
          if (type == EX_NODE_SET) {
            map_nodeset_data(nodesetOwnedNodes[ns], (double*)data, dbldata);
            ierr = ex_put_partial_set_dist_fact(get_file_pointer(), type, id, proc_offset+1, file_count,
						TOPTR(dbldata));
          } else {
            ierr = ex_put_partial_set_dist_fact(get_file_pointer(), type, id, proc_offset+1, num_to_get,
						static_cast<double*>(out_data));
          }
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        } else {
          num_to_get = Ioss::Utils::field_warning(ns, field, "output");
        }
      } else if (role == Ioss::Field::TRANSIENT) {
        // Check if the specified field exists on this element block.
        // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
        // exist on the database as scalars with the appropriate
        // extensions.

        // Transfer each component of the variable into 'data' and then
        // output.  Need temporary storage area of size 'number of
        // elements in this block.
        write_entity_transient_field(type, field, ns, entity_count, data);

      } else if (role == Ioss::Field::ATTRIBUTE) {
        num_to_get = write_attribute_field(type, field, ns, data);

      } else if (role == Ioss::Field::REDUCTION) {
        store_reduction_field(type, field, ns, data);
      }
      return num_to_get;
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet* ns,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      return put_Xset_field_internal(EX_NODE_SET, ns, field, data, data_size);
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet* ns,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      return put_Xset_field_internal(EX_EDGE_SET, ns, field, data, data_size);
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet* ns,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      return put_Xset_field_internal(EX_FACE_SET, ns, field, data, data_size);
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet* ns,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      return put_Xset_field_internal(EX_ELEM_SET, ns, field, data, data_size);
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::SideSet* fs,
                                           const Ioss::Field& field,
                                           void */* data */, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);
      if (field.get_name() == "ids") {
        // Do nothing, just handles an idiosyncracy of the GroupingEntity
      } else {
        num_to_get = Ioss::Utils::field_warning(fs, field, "output");
      }
      return num_to_get;
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::CommSet* cs,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);
      return num_to_get;
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock* fb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);
      int64_t id = Ioex::get_id(fb, EX_SIDE_SET, &ids_);

      size_t entity_count = fb->get_property("entity_count").get_int();
      size_t offset = fb->get_property("set_offset").get_int();

      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
        if (field.get_name() == "side_ids" &&
	    fb->name() == "universal_sideset") {
          // The side ids are being stored as the distribution factor
          // field on the universal sideset.  There should be no other
          // side sets that request this field...  (Eventually,
          // create an id field to store this info.

          // Need to convert 'ints' to 'double' for storage on mesh...
          // FIX 64
          int* ids = static_cast<int*>(data);
          std::vector<double> real_ids(num_to_get);
          for (size_t i = 0; i < num_to_get; i++) {
            real_ids[i] = static_cast<double>(ids[i]);
          }
          int ierr = ex_put_partial_set_dist_fact(get_file_pointer(),  EX_SIDE_SET, id,
						  offset+1, entity_count, TOPTR(real_ids));
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        }

        else if (field.get_name() == "side_ids") {
        }

        else if (field.get_name() == "ids") {
          // =============================================================
          // NOTE: Code is currently commented out since we have
          // redundant ways of getting the data (element/side) out to
          // the database.  The 'ids' field method relies on a numbering
          // kluge, so for now trying the 'element_side' field...
          // =============================================================
        }

        else if (field.get_name() == "distribution_factors") {
          int ierr;
          size_t df_offset = fb->get_property("set_df_offset").get_int();
          size_t proc_df_offset = fb->get_property("processor_df_offset").get_int();
          size_t df_count  = fb->get_property("distribution_factor_count").get_int();
          ierr = ex_put_partial_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id, proc_df_offset+df_offset+1,
					      df_count, static_cast<double*>(data));
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        } else if (field.get_name() == "element_side") {
          // In exodusII, the 'side block' is stored as a sideset.  A
          // sideset has a list of elements and a corresponding local
          // element side (1-based)

          // The 'data' passed into the function is stored as a
          // 2D vector e0,f0,e1,f1,... (e=element, f=side)

          // To avoid overwriting the passed in data, we allocate
          // two arrays to store the data for this sideset.

          // The element_id passed in is the global id; we need to
          // output the local id.

          // Allocate space for local side number and element numbers
          // numbers.
          // See if edges or faces...
          size_t side_offset = Ioss::Utils::get_side_offset(fb);

          size_t index = 0;

          size_t proc_offset = 0;
          if (fb->property_exists("processor_offset"))
            proc_offset = fb->get_property("processor_offset").get_int();

          if (field.get_type() == Ioss::Field::INTEGER) {
            Ioss::IntVector element(num_to_get);
            Ioss::IntVector side(num_to_get);
            int *el_side = (int *)data;

            for (size_t i=0; i < num_to_get; i++) {
              element[i] = elemMap.global_to_local(el_side[index++]);
              side[i]    = el_side[index++]+side_offset;
            }

            map_local_to_global_implicit(TOPTR(element), num_to_get, elemGlobalImplicitMap);
            int ierr = ex_put_partial_set(get_file_pointer(), EX_SIDE_SET, id, proc_offset+offset+1, num_to_get,
					  TOPTR(element), TOPTR(side));
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          } else {
            Ioss::Int64Vector element(num_to_get);
            Ioss::Int64Vector side(num_to_get);
            int64_t *el_side = (int64_t *)data;

            for (size_t i=0; i < num_to_get; i++) {
              element[i] = elemMap.global_to_local(el_side[index++]);
              side[i]    = el_side[index++]+side_offset;
            }

            map_local_to_global_implicit(TOPTR(element), num_to_get, elemGlobalImplicitMap);
            int ierr = ex_put_partial_set(get_file_pointer(), EX_SIDE_SET, id, proc_offset+offset+1, num_to_get,
					  TOPTR(element), TOPTR(side));
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }

        } else if (field.get_name() == "element_side_raw") {
          // In exodusII, the 'side block' is stored as a sideset.  A
          // sideset has a list of elements and a corresponding local
          // element side (1-based)

          // The 'data' passed into the function is stored as a
          // 2D vector e0,f0,e1,f1,... (e=element, f=side)

          // To avoid overwriting the passed in data, we allocate
          // two arrays to store the data for this sideset.

          // The element_id passed in is the local id.

          // See if edges or faces...
          size_t side_offset = Ioss::Utils::get_side_offset(fb);

          size_t index = 0;
          if (field.get_type() == Ioss::Field::INTEGER) {
            Ioss::IntVector element(num_to_get);
            Ioss::IntVector side(num_to_get);
            int *el_side = (int *)data;

            for (size_t i=0; i < num_to_get; i++) {
              element[i] = el_side[index++];
              side[i]    = el_side[index++]+side_offset;
            }

            int ierr = ex_put_partial_set(get_file_pointer(), EX_SIDE_SET, id, offset+1, entity_count, TOPTR(element), TOPTR(side));
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          } else {
            Ioss::Int64Vector element(num_to_get);
            Ioss::Int64Vector side(num_to_get);
            int64_t *el_side = (int64_t *)data;

            for (size_t i=0; i < num_to_get; i++) {
              element[i] = el_side[index++];
              side[i]    = el_side[index++]+side_offset;
            }

            int ierr = ex_put_partial_set(get_file_pointer(), EX_SIDE_SET, id, offset+1, entity_count, TOPTR(element), TOPTR(side));
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }

        } else if (field.get_name() == "connectivity") {
          // Do nothing, just handles an idiosyncracy of the GroupingEntity
        } else if (field.get_name() == "connectivity_raw") {
          // Do nothing, just handles an idiosyncracy of the GroupingEntity
        } else {
          num_to_get = Ioss::Utils::field_warning(fb, field, "output");
        }
      } else if (role == Ioss::Field::TRANSIENT) {
        // Check if the specified field exists on this block.
        // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
        // exist on the database as scalars with the appropriate
        // extensions.

        // Transfer each component of the variable into 'data' and then
        // output.  Need temporary storage area of size 'number of
        // entities in this block.
        write_entity_transient_field(EX_SIDE_SET, field, fb, entity_count, data);

      } else if (role == Ioss::Field::ATTRIBUTE) {
        num_to_get = write_attribute_field(EX_SIDE_SET, field, fb, data);

      } else if (role == Ioss::Field::REDUCTION) {
        store_reduction_field(EX_SIDE_SET, field, fb, data);
      }
      return num_to_get;
    }

    void DatabaseIO::write_meta_data()
    {
      Ioss::Region *region = get_region();

      Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();
      assert(node_blocks.size() == 1);
      nodeCount =        node_blocks[0]->get_property("entity_count").get_int();
      spatialDimension = node_blocks[0]->get_property("component_degree").get_int();

      // See if the nodeOwningProcessor vector has been populated...
      size_t locally_owned = nodeCount;
      if (node_blocks[0]->property_exists("locally_owned_count")) {
        locally_owned = node_blocks[0]->get_property("locally_owned_count").get_int();
      } else if (!nodeOwningProcessor.empty()) {
        locally_owned = std::count(nodeOwningProcessor.begin(), nodeOwningProcessor.end(), myProcessor);
      } 

      char the_title[max_line_length+1];

      // Title...
      if (region->property_exists("title")) {
        std::string title_str = region->get_property("title").get_string();
        std::strncpy(the_title, title_str.c_str(), max_line_length);
      } else {
        std::strncpy(the_title, "Sierra Output Default Title", max_line_length);
      }
      the_title[max_line_length] = '\0';

      Iopx::Mesh mesh(spatialDimension, the_title);
      Iopx::NodeBlock N(*node_blocks[0]);
      mesh.nodeblocks.push_back(N);

      // Edge Blocks --
      {
        Ioss::EdgeBlockContainer edge_blocks = region->get_edge_blocks();
        assert(check_block_order(edge_blocks));
        Ioss::EdgeBlockContainer::const_iterator I;
        // Set ids of all entities that have "id" property...
        for (I=edge_blocks.begin(); I != edge_blocks.end(); ++I) {
          Ioex::set_id(*I, EX_EDGE_BLOCK, &ids_);
        }

        edgeCount = 0;
        for (I=edge_blocks.begin(); I != edge_blocks.end(); ++I) {
          edgeCount += (*I)->get_property("entity_count").get_int();
          // Set ids of all entities that do not have "id" property...
          Ioex::get_id(*I, EX_EDGE_BLOCK, &ids_);
          Iopx::EdgeBlock T(*(*I));
          if (std::find(mesh.edgeblocks.begin(), mesh.edgeblocks.end(), T) == mesh.edgeblocks.end()) {
            mesh.edgeblocks.push_back(T);
          }
        }
        m_groupCount[EX_EDGE_BLOCK] = mesh.edgeblocks.size();
      }

      // Face Blocks --
      {
        Ioss::FaceBlockContainer face_blocks = region->get_face_blocks();
        assert(check_block_order(face_blocks));
        Ioss::FaceBlockContainer::const_iterator I;
        // Set ids of all entities that have "id" property...
        for (I=face_blocks.begin(); I != face_blocks.end(); ++I) {
          Ioex::set_id(*I, EX_FACE_BLOCK, &ids_);
        }

        faceCount = 0;
        for (I=face_blocks.begin(); I != face_blocks.end(); ++I) {
          faceCount += (*I)->get_property("entity_count").get_int();
          // Set ids of all entities that do not have "id" property...
          Ioex::get_id(*I, EX_FACE_BLOCK, &ids_);
          Iopx::FaceBlock T(*(*I));
          if (std::find(mesh.faceblocks.begin(), mesh.faceblocks.end(), T) == mesh.faceblocks.end()) {
            mesh.faceblocks.push_back(T);
          }
        }
        m_groupCount[EX_FACE_BLOCK] = mesh.faceblocks.size();
      }

      // Element Blocks --
      {
        Ioss::ElementBlockContainer element_blocks = region->get_element_blocks();
        assert(check_block_order(element_blocks));
        Ioss::ElementBlockContainer::const_iterator I;
        // Set ids of all entities that have "id" property...
        for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
          Ioex::set_id(*I, EX_ELEM_BLOCK, &ids_);
        }

        elementCount = 0;
        for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
          elementCount += (*I)->get_property("entity_count").get_int();
          // Set ids of all entities that do not have "id" property...
          Ioex::get_id(*I, EX_ELEM_BLOCK, &ids_);
          Iopx::ElemBlock T(*(*I));
          if (std::find(mesh.elemblocks.begin(), mesh.elemblocks.end(), T) == mesh.elemblocks.end()) {
            mesh.elemblocks.push_back(T);
          }
        }
        m_groupCount[EX_ELEM_BLOCK] = mesh.elemblocks.size();
      }

      // Nodesets ...
      {
        Ioss::NodeSetContainer nodesets = region->get_nodesets();
        Ioss::NodeSetContainer::const_iterator I;
        for (I=nodesets.begin(); I != nodesets.end(); ++I) {
          Ioex::set_id(*I, EX_NODE_SET, &ids_);
        }

        for (I=nodesets.begin(); I != nodesets.end(); ++I) {
          Ioex::get_id(*I, EX_NODE_SET, &ids_);
          const Iopx::NodeSet T(*(*I));
          if (std::find(mesh.nodesets.begin(), mesh.nodesets.end(), T) == mesh.nodesets.end()) {
            mesh.nodesets.push_back(T);
          }
        }
        m_groupCount[EX_NODE_SET] = mesh.nodesets.size();
      }

      // Edgesets ...
      {
        Ioss::EdgeSetContainer edgesets = region->get_edgesets();
        Ioss::EdgeSetContainer::const_iterator I;
        for (I=edgesets.begin(); I != edgesets.end(); ++I) {
          Ioex::set_id(*I, EX_EDGE_SET, &ids_);
        }

        for (I=edgesets.begin(); I != edgesets.end(); ++I) {
          Ioex::get_id(*I, EX_EDGE_SET, &ids_);
          const Iopx::EdgeSet T(*(*I));
          if (std::find(mesh.edgesets.begin(), mesh.edgesets.end(), T) == mesh.edgesets.end()) {
            mesh.edgesets.push_back(T);
          }
        }
        m_groupCount[EX_EDGE_SET] = mesh.edgesets.size();
      }

      // Facesets ...
      {
        Ioss::FaceSetContainer facesets = region->get_facesets();
        Ioss::FaceSetContainer::const_iterator I;
        for (I=facesets.begin(); I != facesets.end(); ++I) {
          Ioex::set_id(*I, EX_FACE_SET, &ids_);
        }

        for (I=facesets.begin(); I != facesets.end(); ++I) {
          Ioex::get_id(*I, EX_FACE_SET, &ids_);
          const Iopx::FaceSet T(*(*I));
          if (std::find(mesh.facesets.begin(), mesh.facesets.end(), T) == mesh.facesets.end()) {
            mesh.facesets.push_back(T);
          }
        }
        m_groupCount[EX_FACE_SET] = mesh.facesets.size();
      }

      // Elementsets ...
      {
        Ioss::ElementSetContainer elementsets = region->get_elementsets();
        Ioss::ElementSetContainer::const_iterator I;
        for (I=elementsets.begin(); I != elementsets.end(); ++I) {
          Ioex::set_id(*I, EX_ELEM_SET, &ids_);
        }

        for (I=elementsets.begin(); I != elementsets.end(); ++I) {
          Ioex::get_id(*I, EX_ELEM_SET, &ids_);
          const Iopx::ElemSet T(*(*I));
          if (std::find(mesh.elemsets.begin(), mesh.elemsets.end(), T) == mesh.elemsets.end()) {
            mesh.elemsets.push_back(T);
          }
        }
        m_groupCount[EX_ELEM_SET] = mesh.elemsets.size();
      }

      // SideSets ...
      {
        Ioss::SideSetContainer ssets = region->get_sidesets();
        Ioss::SideSetContainer::const_iterator I;

        for (I=ssets.begin(); I != ssets.end(); ++I) {
          Ioex::set_id(*I, EX_SIDE_SET, &ids_);
        }

        // Get entity counts for all face sets... Create SideSets.
        for (I=ssets.begin(); I != ssets.end(); ++I) {

          Ioex::get_id(*I, EX_SIDE_SET, &ids_);
          int64_t id = (*I)->get_property("id").get_int();
          int64_t entity_count = 0;
          int64_t df_count = 0;

          Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
          Ioss::SideBlockContainer::const_iterator J;
          for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
            // Add  "*_offset" properties to specify at what offset
            // the data for this block appears in the containing set.
            Ioss::SideBlock *new_block = const_cast<Ioss::SideBlock*>(*J);
            new_block->property_add(Ioss::Property("set_offset", entity_count));
            new_block->property_add(Ioss::Property("set_df_offset", df_count));

            // If combining sideblocks into sidesets on output, then
            // the id of the sideblock must be the same as the sideset
            // id.
            if (new_block->property_exists("id")) {
              new_block->property_erase("id");
            }
            new_block->property_add(Ioss::Property("id", id));

            entity_count += (*J)->get_property("entity_count").get_int();
            df_count     += (*J)->get_property("distribution_factor_count").get_int();
          }
          Ioss::SideSet *new_entity = const_cast<Ioss::SideSet*>(*I);
          new_entity->property_add(Ioss::Property("entity_count", entity_count));
          new_entity->property_add(Ioss::Property("distribution_factor_count", df_count));
        }

        for (I=ssets.begin(); I != ssets.end(); ++I) {
          // Add a SideSet corresponding to this SideSet/SideBlock
          Iopx::SideSet T(*(*I));
          if (std::find(mesh.sidesets.begin(), mesh.sidesets.end(), T) == mesh.sidesets.end()) {
            mesh.sidesets.push_back(T);
          }
        }
        m_groupCount[EX_SIDE_SET] = mesh.sidesets.size();
      }

      {
        // Write the metadata to the exodusII file...
        Iopx::Internals data(get_file_pointer(), maximumNameLength, util());
        int ierr = data.write_meta_data(mesh);

        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }

      metaDataWritten = true;

      // Output node map...
      output_node_map();

      output_other_meta_data();

      // Set the processor offset property. Specifies where in the global list, the data from this
      // processor begins...

      node_blocks[0]->property_add(Ioss::Property("processor_offset", mesh.nodeblocks[0].procOffset));
      Ioss::EdgeBlockContainer edge_blocks = region->get_edge_blocks();
      for (size_t i=0; i < edge_blocks.size(); i++) {
        edge_blocks[i]->property_add(Ioss::Property("processor_offset", mesh.edgeblocks[i].procOffset));
      }
      Ioss::FaceBlockContainer face_blocks = region->get_face_blocks();
      for (size_t i=0; i < face_blocks.size(); i++) {
        face_blocks[i]->property_add(Ioss::Property("processor_offset", mesh.faceblocks[i].procOffset));
      }

      int64_t offset = 0; // Offset into global element map...
      Ioss::ElementBlockContainer element_blocks = region->get_element_blocks();
      for (size_t i=0; i < element_blocks.size(); i++) {
        element_blocks[i]->property_add(Ioss::Property("global_map_offset", offset));
        offset += mesh.elemblocks[i].entityCount;
        element_blocks[i]->property_add(Ioss::Property("processor_offset", mesh.elemblocks[i].procOffset));
      }

      Ioss::NodeSetContainer nodesets = region->get_nodesets();
      for (size_t i=0; i < nodesets.size(); i++) {
        nodesets[i]->property_add(Ioss::Property("processor_offset", mesh.nodesets[i].procOffset));
      }
      Ioss::EdgeSetContainer edgesets = region->get_edgesets();
      for (size_t i=0; i < edgesets.size(); i++) {
        edgesets[i]->property_add(Ioss::Property("processor_offset", mesh.edgesets[i].procOffset));
      }
      Ioss::FaceSetContainer facesets = region->get_facesets();
      for (size_t i=0; i < facesets.size(); i++) {
        facesets[i]->property_add(Ioss::Property("processor_offset", mesh.facesets[i].procOffset));
      }
      Ioss::ElementSetContainer elementsets = region->get_elementsets();
      for (size_t i=0; i < facesets.size(); i++) {
        elementsets[i]->property_add(Ioss::Property("processor_offset", mesh.elemsets[i].procOffset));
      }

      Ioss::SideSetContainer ssets = region->get_sidesets();
      for (size_t i=0; i < ssets.size(); i++) {
        ssets[i]->property_add(Ioss::Property("processor_offset", mesh.sidesets[i].procOffset));
        ssets[i]->property_add(Ioss::Property("processor_df_offset", mesh.sidesets[i].dfProcOffset));

        // Propogate down to owned sideblocks...
        Ioss::SideBlockContainer side_blocks = ssets[i]->get_side_blocks();
        Ioss::SideBlockContainer::const_iterator J;
        for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
          (*J)->property_add(Ioss::Property("processor_offset", mesh.sidesets[i].procOffset));
          (*J)->property_add(Ioss::Property("processor_df_offset", mesh.sidesets[i].dfProcOffset));
        }
      }

    }

    void DatabaseIO::create_implicit_global_map() const
    {
      // If the node is locally owned, then its position is basically
      // determined by removing all shared nodes from the list and
      // then compressing the list. This location plus the proc_offset
      // gives its location in the global-implicit file.
      //
      // Do this over in the DecompositionData class since it has
      // several utilities in place for MPI communication.

      DecompositionData<int64_t> compose(Ioss::PropertyManager(), util().communicator());
      int64_t locally_owned_count = 0;
      int64_t processor_offset = 0;
      compose.create_implicit_global_map(nodeOwningProcessor, nodeGlobalImplicitMap,
					 nodeMap, &locally_owned_count, &processor_offset);

      Ioss::NodeBlockContainer node_blocks = get_region()->get_node_blocks();
      if (!node_blocks[0]->property_exists("locally_owned_count")) 
        node_blocks[0]->property_add(Ioss::Property("locally_owned_count", locally_owned_count));
      if (!node_blocks[0]->property_exists("processor_offset")) 
        node_blocks[0]->property_add(Ioss::Property("processor_offset", processor_offset));

      output_node_map();
    }      

    void DatabaseIO::output_node_map() const
    {
      // Write the partial nodemap to the database...  This is called
      // two times -- once from create_implicit_global_map() and once
      // from write_meta_data().  It will only output the map if
      // the metadata has been written to the output database AND if
      // the nodeMap.map and nodeGlobalImplicitMap are defined.

      if (metaDataWritten && !nodeMap.map.empty() && !nodeGlobalImplicitMap.empty()) {
        Ioss::NodeBlockContainer node_blocks = get_region()->get_node_blocks();
        assert(node_blocks[0]->property_exists("processor_offset"));
        assert(node_blocks[0]->property_exists("locally_owned_count"));
        size_t processor_offset = node_blocks[0]->get_property("processor_offset").get_int();
        size_t locally_owned_count =  node_blocks[0]->get_property("locally_owned_count").get_int();

        int ierr;
        if (int_byte_size_api() == 4) {
          std::vector<int> file_ids; file_ids.reserve(locally_owned_count);
          map_data(nodeOwningProcessor, myProcessor, &nodeMap.map[1], file_ids);
          ierr = ex_put_partial_id_map(get_file_pointer(), EX_NODE_MAP,
				       processor_offset+1, locally_owned_count, TOPTR(file_ids));
        } else {
          std::vector<int64_t> file_ids; file_ids.reserve(locally_owned_count);
          map_data(nodeOwningProcessor, myProcessor, &nodeMap.map[1], file_ids);
          ierr = ex_put_partial_id_map(get_file_pointer(), EX_NODE_MAP,
				       processor_offset+1, locally_owned_count, TOPTR(file_ids));
        }
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
    }

  } // End of namespace
