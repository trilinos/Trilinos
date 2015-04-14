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
#include <Ioss_ElementTopology.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_SerializeIO.h>
#include <Ioss_SurfaceSplit.h>
#include <Ioss_FileInfo.h>
#include <Ioss_Utils.h>
#include <assert.h>
#include <exodusII.h>
#include <exodus/Ioex_Utils.h>
#include <fpp_exo/Iofx_DatabaseIO.h>
#include <fpp_exo/Iofx_Internals.h>
#include <float.h>
#include <stddef.h>
#include <sys/select.h>
#include <time.h>
#include <tokenize.h>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
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

#ifdef HAVE_MPI    
#include "Ioss_FileInfo.h"
#endif

// ========================================================================
// Static internal helper functions
// ========================================================================
namespace {
  const size_t max_line_length   = MAX_LINE_LENGTH;

  const std::string SEP() {return std::string("@");} // Separator for attribute offset storage
  const std::string SCALAR()     {return std::string("scalar");}
  const std::string VECTOR3D()   {return std::string("vector_3d");}
  const std::string SYM_TENSOR() {return std::string("sym_tensor_33");}

  const char *complex_suffix[] = {".re", ".im"};

  void get_connectivity_data(int exoid, void *data, ex_entity_type type, ex_entity_id id, int position)
  {
    int ierr = 0;
    if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
      int64_t* conn[3];
      conn[0] = NULL;
      conn[1] = NULL;
      conn[2] = NULL;
      conn[position] = static_cast<int64_t*>(data);
      ierr = ex_get_conn(exoid, type, id, conn[0], conn[1], conn[2]);
    } else {
      int* conn[3];
      conn[0] = NULL;
      conn[1] = NULL;
      conn[2] = NULL;
      conn[position] = static_cast<int*>(data);
      ierr = ex_get_conn(exoid, type, id, conn[0], conn[1], conn[2]);
    }
    if (ierr < 0)
      Ioex::exodus_error(exoid, __LINE__, -1);
  }

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
}

namespace Iofx {
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
                         Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                         const Ioss::PropertyManager &props) :
    Ioex::DatabaseIO(region, filename, db_usage, communicator, props)
  {}

  bool DatabaseIO::ok(bool write_message, std::string *error_msg, int *bad_count) const
  {

    // Returns the number of processors on which this file is *NOT* ok.
    // If the file exists everywhere and is ok on all processors, it will return 0.
    // If the file does not exist anywhere or is bad on all processors, it returns parallel_size.
    // If the file is ok on some and bad on some processors, it will return the count of where it is bad.
    
    if (fileExists) {
      // File has already been opened at least once...
      return dbState != Ioss::STATE_INVALID;
    }

    int app_opt_val = ex_opts(EX_DEFAULT); 

    // File has not yet been opened, so verify that it is readable or
    // writable.  HOWEVER, *DO NOT* overwrite the file here if it exists
    // and the mode is WRITE.
    int cpu_word_size = sizeof(double);
    int io_word_size  = 0;
    float version;
    std::string decoded_filename = util().decode_filename(get_filename(), isParallel);

    int exodus_file_ptr = ex_open(decoded_filename.c_str(), EX_READ,
                                  &cpu_word_size, &io_word_size, &version);

    if (!is_input() && exodus_file_ptr < 0) {
      // File didn't exist above, but this OK if is an output
      // file. See if we can create it...
      int mode = 0;
      if (int_byte_size_api() == 8)
        mode |= EX_ALL_INT64_DB;

      exodus_file_ptr = ex_create(decoded_filename.c_str(), exodusMode|mode,
                                  &cpu_word_size, &dbRealWordSize);
    }

    ex_opts(app_opt_val); // Reset back to what it was.

    // Check for valid exodus_file_ptr (valid >= 0; invalid < 0)
    int global_file_ptr = exodus_file_ptr;
    if (isParallel) {
      global_file_ptr = util().global_minmax(exodus_file_ptr, Ioss::ParallelUtils::DO_MIN);
    }
    if (global_file_ptr < 0 && (write_message || error_msg != NULL || bad_count != NULL)) {
      Ioss::IntVector status;
      if (isParallel) {
        util().all_gather(exodus_file_ptr, status);
      }
      else {
        status.push_back(exodus_file_ptr);
      }

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

    if (Ioss::SerializeIO::isEnabled()) {
      if (!Ioss::SerializeIO::inBarrier()) {
        std::ostringstream errmsg;
        errmsg << "ERROR: Process " << Ioss::SerializeIO::getRank()
               << " is attempting to do I/O without serialized I/O";
        IOSS_ERROR(errmsg);
      }

      if (!Ioss::SerializeIO::inMyGroup()) {
        std::ostringstream errmsg;
        errmsg << "ERROR: Process " << Ioss::SerializeIO::getRank()
               << " is attempting to do I/O while " << Ioss::SerializeIO::getOwner()
               << " owns the token";
        IOSS_ERROR(errmsg);
      }
    }

    if (exodusFilePtr < 0) {
      int cpu_word_size = sizeof(double);
      int io_word_size  = 0;
      float version;
      std::string decoded_filename = util().decode_filename(get_filename(),
                                                            isParallel);

      int mode = exodusMode;
      if (int_byte_size_api() == 8)
        mode |= EX_ALL_INT64_API;

      if (is_input()) {
        exodusFilePtr = ex_open(decoded_filename.c_str(), EX_READ|mode,
                                &cpu_word_size, &io_word_size, &version);
      } else {
        if (fileExists) {
          exodusFilePtr = ex_open(decoded_filename.c_str(), EX_WRITE|mode,
                                  &cpu_word_size, &io_word_size, &version);
        } else {
          // If the first write for this file, create it...
          if (int_byte_size_api() == 8)
            mode |= EX_ALL_INT64_DB;
          exodusFilePtr = ex_create(decoded_filename.c_str(), mode,
                                    &cpu_word_size, &dbRealWordSize);
          if (exodusFilePtr < 0) {
            dbState = Ioss::STATE_INVALID;
            // NOTE: Code will not continue past this call...
            std::ostringstream errmsg;
            errmsg << "ERROR: Cannot create specified file '" << decoded_filename << "'";
            IOSS_ERROR(errmsg);
          }
        }
      }

      if (exodusFilePtr < 0) {
        dbState = Ioss::STATE_INVALID;
        fileExists = false;
        // NOTE: Code will not continue past this call...
        std::ostringstream errmsg;
        errmsg << "ERROR: Problem opening specified file '" << decoded_filename << "'";
        IOSS_ERROR(errmsg);
      }

      if (is_input()) {
        // Check for maximum name length used on the input file.
        int max_name_length = ex_inquire_int(exodusFilePtr, EX_INQ_DB_MAX_USED_NAME_LENGTH);
        if (max_name_length > maximumNameLength) {
          maximumNameLength = max_name_length;
        }
      }

      ex_set_max_name_length(exodusFilePtr, maximumNameLength);

      // Check properties handled post-create/open...
      if (properties.exists("COMPRESSION_LEVEL")) {
        int comp_level = properties.get("COMPRESSION_LEVEL").get_int();
        ex_set_option(exodusFilePtr, EX_OPT_COMPRESSION_LEVEL, comp_level);
      }
      if (properties.exists("COMPRESSION_SHUFFLE")) {
        int shuffle = properties.get("COMPRESSION_SHUFFLE").get_int();
        ex_set_option(exodusFilePtr, EX_OPT_COMPRESSION_SHUFFLE, shuffle);
      }

      if (!m_groupName.empty()) {
       ex_get_group_id(exodusFilePtr, m_groupName.c_str(), &exodusFilePtr);
      }
    }
    assert(exodusFilePtr >= 0);
    fileExists = true;
    return exodusFilePtr;
  }

  void DatabaseIO::read_meta_data()
  {
    // If this is a HISTORY file, there isn't really any metadata
    // Other than a single node and single element.  Just hardwire
    // it here (needed when appending to existing history file)
    if (dbUsage == Ioss::WRITE_HISTORY) {
      if (myProcessor == 0) {
        nodeCount = 1;
        elementCount = 1;
        Ioss::NodeBlock *nb = new Ioss::NodeBlock(this, "nodeblock_1", 1, 3);
        get_region()->add(nb);

        // Element Block
        Ioss::ElementBlock *eb = new Ioss::ElementBlock(this, "e1", "sphere", 1);
        eb->property_add(Ioss::Property("id", 1));
        get_region()->add(eb);
        get_step_times();
        add_region_fields();
      }
      return;
    }

    {
      Ioss::SerializeIO serializeIO__(this);

      Ioex::check_processor_info(get_file_pointer(), util().parallel_size(), myProcessor);

      read_region();
      read_communication_metadata();
    }

    get_step_times();

    get_nodeblocks();
    get_edgeblocks();
    get_faceblocks();
    get_elemblocks();

    check_side_topology();

    get_sidesets();
    get_nodesets();
    get_edgesets();
    get_facesets();
    get_elemsets();

    get_commsets();

    handle_groups();

    add_region_fields();

    if (!is_input() && open_create_behavior() == Ioss::DB_APPEND) {
      get_map(EX_NODE_BLOCK);
      get_map(EX_EDGE_BLOCK);
      get_map(EX_FACE_BLOCK);
      get_map(EX_ELEM_BLOCK);
    }

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

    spatialDimension = info.num_dim;
    nodeCount = info.num_nodes;
    edgeCount = info.num_edge;
    faceCount = info.num_face;
    elementCount = info.num_elem;

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
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      IOSS_WARNING << "No nodes were found in the model, file '" << decoded_filename << "'";
    } else if (nodeCount < 0) {
      // NOTE: Code will not continue past this call...
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      std::ostringstream errmsg;
      errmsg << "ERROR: Negative node count was found in the model\n"
             << "       File: '" << decoded_filename << "'.\n";
      IOSS_ERROR(errmsg);
    }

    if (elementCount == 0) {
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      IOSS_WARNING << "No elements were found in the model, file: '" << decoded_filename << "'";
    }

    if (elementCount < 0) {
      // NOTE: Code will not continue past this call...
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      std::ostringstream errmsg;
      errmsg << "ERROR: Negative element count was found in the model, file: '"
             << decoded_filename << "'";
      IOSS_ERROR(errmsg);
    }

    if (elementCount > 0 && m_groupCount[EX_ELEM_BLOCK] <= 0) {
      // NOTE: Code will not continue past this call...
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      std::ostringstream errmsg;
      errmsg << "ERROR: No element blocks were found in the model, file: '" << decoded_filename << "'";
      IOSS_ERROR(errmsg);
    }

    Ioss::Region *this_region = get_region();

    // See if any coordinate frames exist on mesh.  If so, define them on region.
    Ioex::add_coordinate_frames(get_file_pointer(), this_region);
    
    this_region->property_add(Ioss::Property(std::string("title"), info.title));
    this_region->property_add(Ioss::Property(std::string("spatial_dimension"),
                                             spatialDimension));

    // Get QA records from database and add to qaRecords...
    int num_qa = ex_inquire_int(get_file_pointer(), EX_INQ_QA);    
    if (num_qa > 0) {
      struct qa_element {
        char *qa_record[1][4];
      };
    
      qa_element *qa = new qa_element[num_qa];
      for (int i=0; i < num_qa; i++) {
        for (int j=0; j < 4; j++) {
          qa[i].qa_record[0][j] = new char[MAX_STR_LENGTH+1];
        }
      }

      ex_get_qa(get_file_pointer(), qa[0].qa_record);
      for (int i=0; i < num_qa; i++) {
        add_qa_record(qa[i].qa_record[0][0], qa[i].qa_record[0][1], qa[i].qa_record[0][2], qa[i].qa_record[0][3]);
      }
      for (int i=0; i < num_qa; i++) {
        for (int j=0; j < 4; j++) {
          delete [] qa[i].qa_record[0][j];
        }
      }
      delete [] qa;

    }
    
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
    bool exists = false;
    double last_time = DBL_MAX;
    int timestep_count = 0;
    std::vector<double> tsteps(0);

    if (dbUsage == Ioss::WRITE_HISTORY) {
      if (myProcessor == 0) {
        timestep_count = ex_inquire_int(get_file_pointer(), EX_INQ_TIME);
        if (timestep_count <= 0)
          return;
        
        // For an exodusII file, timesteps are global and are stored in the region.
        // A history file only stores that last time / step
        // Read the timesteps and add them to the region.
        // Since we can't access the Region's stateCount directly, we just add
        // all of the steps and assume the Region is dealing with them directly...
        tsteps.resize(timestep_count);
        int error = ex_get_all_times(get_file_pointer(), TOPTR(tsteps));
        if (error < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        Ioss::Region *this_region = get_region();
        for (int i=0; i < timestep_count; i++) {
          this_region->add_state(tsteps[i]*timeScaleFactor);
        }
      }
    } else {
      {
        Ioss::SerializeIO       serializeIO__(this);
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
	exists = Ioex::read_last_time_attribute(get_file_pointer(), &last_time);
      }
      if (exists && isParallel) {
        // Assume that if it exists on 1 processor, it exists on
        // all... Sync value among processors since could have a
        // corrupt step on only a single database.
        last_time = util().global_minmax(last_time, Ioss::ParallelUtils::DO_MIN);
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
                         << ".\n\tThe data for that step is possibly corrupt since the last time written successfully was "
                         << last_time << ".\n";
          }
        }
      }
    }
  }

  void DatabaseIO::read_communication_metadata()
  {
    // Sierra does not use/need the element communication map data,
    // but the Nemesis api call to get the nodal communication map
    // data also gets the element communication data. So,...,we need to
    // do some redundant reads to get back the element comm data
    // and also some extra memory allocations to give them a place
    // to go.

    // Check that file is nemesis.
    int num_proc;          // Number of processors file was decomposed for
    int num_proc_in_file;  // Number of processors this file has info for
    char file_type[2];         // "s" for scalar, "p" for parallel

    // Get global data (over all processors)
    int64_t global_nodes       = nodeCount;
    int64_t global_elements    = elementCount;
    int global_eblocks     = 0; // unused
    int global_nsets       = 0; // unused
    int global_ssets       = 0; // unused

    int64_t num_external_nodes; // unused
    int64_t num_elem_cmaps     = 0;
    int64_t num_node_cmaps     = 0;
    int64_t num_internal_nodes = nodeCount;
    int64_t num_border_nodes   = 0;
    int64_t num_internal_elems = elementCount;
    int64_t num_border_elems   = 0;

    bool nemesis_file = true;
    // If someone changed to EX_VERBOSE, temporarily change to default
    // so this call does not report an error in the serial case.
    // (See Trac 10774)
    int old_val = ex_opts(EX_DEFAULT); 
    int error = ex_get_init_info(get_file_pointer(),
        &num_proc, &num_proc_in_file, &file_type[0]);
    ex_opts(old_val); // Reset back to what it was.

    if (error < 0) {
      // Not a nemesis file
      nemesis_file = false;
      if (isParallel && util().parallel_size() > 1) {
        std::ostringstream errmsg;
        errmsg << "ERROR: Exodus file does not contain nemesis information.\n";
        IOSS_ERROR(errmsg);
      }
      num_proc = 1;
      num_proc_in_file = 1;
      file_type[0] = 'p';
    } else {
      if (!isParallel) {
        // The file contains nemesis parallel information.
        // Even though we are running in serial, make the information
        // available to the application.
        isSerialParallel = true;
        get_region()->property_add(Ioss::Property("processor_count", num_proc));
      }
    }

    if (isParallel && num_proc != util().parallel_size() && util().parallel_size() > 1) {
      std::ostringstream errmsg;
      errmsg <<  "ERROR: Exodus file was decomposed for " << num_proc
             << " processors; application is currently being run on "
             << util().parallel_size() << " processors";
      IOSS_ERROR(errmsg);

    }
    if (num_proc_in_file != 1) {
      std::ostringstream errmsg;
      errmsg <<"ERROR: Exodus file contains data for " << num_proc_in_file
             << " processors; application requires 1 processor per file.";
      IOSS_ERROR(errmsg);

    }
    if (file_type[0] != 'p') {
      std::ostringstream errmsg;
      errmsg << "ERROR: Exodus file contains scalar nemesis data; application requires parallel nemesis data.";
      IOSS_ERROR(errmsg);
    }

    if (nemesis_file) {
      if (int_byte_size_api() == 4) {
        int nin, nbn, nen, nie, nbe, nnc, nec;
        error = ex_get_loadbal_param(get_file_pointer(),
                                     &nin, &nbn, &nen, &nie, &nbe, &nnc, &nec, myProcessor);
        num_external_nodes = nen;
        num_elem_cmaps     = nec;
        num_node_cmaps     = nnc;
        num_internal_nodes = nin;
        num_border_nodes   = nbn;
        num_internal_elems = nie;
        num_border_elems   = nbe;
      } else {
        error = ex_get_loadbal_param(get_file_pointer(),
                                     &num_internal_nodes,
                                     &num_border_nodes,
                                     &num_external_nodes,
                                     &num_internal_elems,
                                     &num_border_elems,
                                     &num_node_cmaps,
                                     &num_elem_cmaps,
                                     myProcessor);
      }
      if (error < 0)
        Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // A nemesis file typically separates nodes into multiple
      // communication sets by processor.  (each set specifies
      // nodes/elements that communicate with only a single processor).
      // For Sierra, we want a single node commun. map and a single
      // element commun. map specifying all communications so we combine
      // all sets into a single set.

      if (int_byte_size_api() == 4) {
        int gn, ge, geb, gns, gss;
        error = ex_get_init_global(get_file_pointer(), &gn, &ge, &geb, &gns, &gss);
        global_nodes = gn;
        global_elements = ge;
        global_eblocks = geb;
        global_nsets = gns;
        global_ssets = gss;
      } else {
        error = ex_get_init_global(get_file_pointer(), &global_nodes, &global_elements,
                                   &global_eblocks, &global_nsets, &global_ssets);
      }
      if (error < 0)
        Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
    }

    commsetNodeCount = num_node_cmaps;
    commsetElemCount = num_elem_cmaps;

    Ioss::Region *region = get_region();
    region->property_add(Ioss::Property("internal_node_count",
                                        num_internal_nodes));
    region->property_add(Ioss::Property("border_node_count",
                                        num_border_nodes));
    region->property_add(Ioss::Property("internal_element_count",
                                        num_internal_elems));
    region->property_add(Ioss::Property("border_element_count",
                                        num_border_elems));
    region->property_add(Ioss::Property("global_node_count",    global_nodes));
    region->property_add(Ioss::Property("global_element_count", global_elements));
    region->property_add(Ioss::Property("global_element_block_count", global_eblocks));
    region->property_add(Ioss::Property("global_node_set_count", global_nsets));
    region->property_add(Ioss::Property("global_side_set_count", global_ssets));

    // Possibly, the following 4 fields should be nodesets and element
    // sets instead of fields on the region...
    region->field_add(Ioss::Field("internal_nodes", region->field_int_type(),
                                  SCALAR(),
                                  Ioss::Field::COMMUNICATION,
                                  num_internal_nodes));
    region->field_add(Ioss::Field("border_nodes", region->field_int_type(),
                                  SCALAR(),
                                  Ioss::Field::COMMUNICATION,
                                  num_border_nodes));
    region->field_add(Ioss::Field("internal_elements", region->field_int_type(),
                                  SCALAR(),
                                  Ioss::Field::COMMUNICATION,
                                  num_internal_elems));
    region->field_add(Ioss::Field("border_elements", region->field_int_type(),
                                  SCALAR(),
                                  Ioss::Field::COMMUNICATION,
                                  num_border_elems));

    assert(nodeCount    == num_internal_nodes + num_border_nodes);
    assert(elementCount == num_internal_elems + num_border_elems);
  }

  const Ioss::Map& DatabaseIO::get_map(ex_entity_type type) const
  {
    switch (type) {
    case EX_NODE_BLOCK:
    case EX_NODE_SET:
      return get_map(nodeMap, nodeCount, EX_NODE_MAP, EX_INQ_NODE_MAP);

    case EX_ELEM_BLOCK:
    case EX_ELEM_SET:
      return get_map(elemMap, elementCount, EX_ELEM_MAP, EX_INQ_ELEM_MAP);

    case EX_FACE_BLOCK:
    case EX_FACE_SET:
      return get_map(faceMap, faceCount, EX_FACE_MAP, EX_INQ_FACE_MAP);

    case EX_EDGE_BLOCK:
    case EX_EDGE_SET:
      return get_map(edgeMap, edgeCount, EX_EDGE_MAP, EX_INQ_EDGE_MAP);

    default:
      std::ostringstream errmsg;
      errmsg << "INTERNAL ERROR: Invalid map type. "
             << "Something is wrong in the Iofx::DatabaseIO::get_map() function. "
             << "Please report.\n";
      IOSS_ERROR(errmsg);
    }      
  }

  const Ioss::Map& DatabaseIO::get_map(Ioss::Map &entity_map,
                                       int64_t entityCount,
                                       ex_entity_type entity_type,
                                       ex_inquiry inquiry_type) const
  {
    // Allocate space for node number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (entity_map.map.empty()) {
      entity_map.map.resize(entityCount+1);

      if (is_input() || open_create_behavior() == Ioss::DB_APPEND) {

        Ioss::SerializeIO       serializeIO__(this);
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
            int error = 0;
            if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
              error = ex_get_num_map(get_file_pointer(), entity_type, 1, &entity_map.map[1]);
            } else {
              // Ioss stores as 64-bit, read as 32-bit and copy over...
              Ioss::IntVector tmp_map(entity_map.map.size());
              error = ex_get_num_map(get_file_pointer(), entity_type, 1, &tmp_map[1]);
              if (error >= 0)
                std::copy(tmp_map.begin(), tmp_map.end(), entity_map.map.begin());
            }
            if (error >= 0) {
              map_read = true;
            } else {
              // Clear out the vector...
              Ioss::MapContainer().swap(entity_map.map);
              Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              map_read = false;
            }
          }
	  Ioex::delete_exodus_names(names, map_count);
        }

        if (!map_read) {
          int error = 0;
          if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
            error = ex_get_id_map(get_file_pointer(), entity_type, &entity_map.map[1]);
          } else {
            // Ioss stores as 64-bit, read as 32-bit and copy over...
            Ioss::IntVector tmp_map(entity_map.map.size());
            error = ex_get_id_map(get_file_pointer(), entity_type, &tmp_map[1]);
            if (error >= 0)
              std::copy(tmp_map.begin(), tmp_map.end(), entity_map.map.begin());
          }
          if (error < 0) {
            // Clear out the vector...
            Ioss::MapContainer().swap(entity_map.map);
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }
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
        for (int64_t i=1; i < nodeCount+1; i++) {
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
    get_blocks(EX_FACE_BLOCK, 1, "faceblock");
  }

  void DatabaseIO::get_edgeblocks()
  {
    get_blocks(EX_EDGE_BLOCK, 2, "edgeblock");
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

    Ioss::Int64Vector X_block_ids(m_groupCount[entity_type]);

    int error;
    {
      Ioss::SerializeIO serializeIO__(this);

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
    }

    size_t all_X_type_length = m_groupCount[entity_type] * (MAX_STR_LENGTH+1);
    std::vector<char> all_X_type(all_X_type_length);

    Ioss::Int64Vector counts(m_groupCount[entity_type] * 4);
    Ioss::Int64Vector local_X_count(m_groupCount[entity_type]);
    Ioss::Int64Vector global_X_count(m_groupCount[entity_type]);
    int iblk ;

    {
      Ioss::SerializeIO serializeIO__(this);

      for ( iblk = 0; iblk < m_groupCount[entity_type]; iblk++) {
        int index = 4*iblk;
        int64_t id = X_block_ids[iblk];

        char * const X_type = TOPTR(all_X_type) + iblk * (MAX_STR_LENGTH+1);

        ex_block block;
        block.id = id;
        block.type = entity_type;
        error = ex_get_block_param(get_file_pointer(), &block);
        if (error < 0) {Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);}

        local_X_count[iblk] = block.num_entry;

        counts[index+0] = block.num_nodes_per_entry;
        counts[index+1] = block.num_edges_per_entry;
        counts[index+2] = block.num_faces_per_entry;
        counts[index+3] = block.num_attribute;

        if (block.num_entry == 0) {
          std::memset(X_type, 0, MAX_STR_LENGTH + 1);
        } else {
          strncpy(X_type, block.topology, MAX_STR_LENGTH+1);
        }
      }
    }

    // This is a collective call...
    util().attribute_reduction(all_X_type_length, TOPTR(all_X_type));

    // This is a collective call...
    util().global_array_minmax(counts, Ioss::ParallelUtils::DO_MAX);

    // Determine global X count for each X block....
    // Can also get this from an nemesis call, but the data may not always be there in all cases.
    util().global_count(local_X_count, global_X_count);

    // The 'offset' is used to map an X location within an X
    // block to the X 'file descriptor'.  For example, the file
    // descriptor of the 37th X in the 4th block is calculated by:
    // file_descriptor = offset of block 4 + 37 This can also be used to
    // determine which X block an X with a file_descriptor
    // maps into. An particular X block contains all Xs in
    // the range:
    //     offset < file_descriptor <= offset+number_Xs_per_block
    int64_t offset = 0;
    int used_blocks = 0;

    for ( iblk = 0; iblk < m_groupCount[entity_type]; iblk++) {
      int index = 4*iblk;
      int64_t nodes_per_X = counts[index+0];
      int64_t edges_per_X = counts[index+1];
      int64_t faces_per_X = counts[index+2];
      int64_t attributes  = counts[index+3];

      int64_t id = X_block_ids[iblk];
      std::string alias = Ioss::Utils::encode_entity_name(basename, id);
      char * const X_type = TOPTR(all_X_type) + iblk * (MAX_STR_LENGTH+1);

      bool db_has_name = false;
      std::string block_name;
      {
        Ioss::SerializeIO       serializeIO__(this);
        block_name = Ioex::get_entity_name(get_file_pointer(), entity_type, id, basename,
					   maximumNameLength, db_has_name);

      }
      if (get_use_generic_canonical_name()) {
        std::string temp = block_name;
        block_name = alias;
        alias = temp;
      }

      std::string save_type = X_type;
      std::string type = Ioss::Utils::fixup_type(X_type, nodes_per_X, spatialDimension-rank_offset);
      if (local_X_count[iblk] == 0 && type == "") {
        // For an empty block, exodusII does not store the X
        // type information and returns "NULL" If there are no
        // Xs on any processors for this block, it will have
        // an empty type which is invalid and will throw an
        // exception in the XBlock constructor. Try to discern
        // the correct X type based on the block_name.
        std::vector<std::string> tokens;
        Ioss::tokenize(block_name, "_", tokens);
        if (tokens.size() >= 2) {
          // Check whether last token names an X topology type...
          const Ioss::ElementTopology *topology = Ioss::ElementTopology::factory(tokens[tokens.size()-1], true);
          if (topology != NULL) {
            type = topology->name();
          }
        }
      }

      if (type == "null" || type == "") {
        // If we have no idea what the topology type for an empty
        // X block is, call it "unknown"
        type = "unknown";

        // If there are no Xs on any processor for this block and
        // we have no idea what the topology type is, skip it...
        if (global_X_count[iblk] == 0) {
          continue;
        }
      }

      Ioss::EntityBlock *block = NULL;
      if (entity_type == EX_ELEM_BLOCK) {
        Ioss::ElementBlock *eblock = new Ioss::ElementBlock(this, block_name, type, local_X_count[iblk]);
        block = eblock;
        get_region()->add(eblock);
      } else if (entity_type == EX_FACE_BLOCK) {
        Ioss::FaceBlock *fblock = new Ioss::FaceBlock(this, block_name, type, local_X_count[iblk]);
        block = fblock;
        get_region()->add(fblock);
      } else if (entity_type == EX_EDGE_BLOCK) {
        Ioss::EdgeBlock *eblock = new Ioss::EdgeBlock(this, block_name, type, local_X_count[iblk]);
        block = eblock;
        get_region()->add(eblock);
      } else {
        std::ostringstream errmsg;
        errmsg << "ERROR: Invalid type in get_blocks()";
        IOSS_ERROR(errmsg);
      }

      // See which connectivity options were defined for this block.
      // X -> Node is always defined.
      // X -> Face?
      if (faces_per_X > 0 && rank_offset < 1) {
        std::string storage = "Real["+Ioss::Utils::to_string(faces_per_X)+"]";
        block->field_add(Ioss::Field("connectivity_face",
                                     block->field_int_type(), storage, Ioss::Field::MESH,
                                     local_X_count[iblk]));
      }
      // X -> Edge?
      if (edges_per_X > 0 && rank_offset < 2) {
        std::string storage = "Real["+Ioss::Utils::to_string(edges_per_X)+"]";
        block->field_add(Ioss::Field("connectivity_edge",
                                     block->field_int_type(), storage, Ioss::Field::MESH,
                                     local_X_count[iblk]));
      }

      block->property_add(Ioss::Property("id", id)); // Do before adding for better error messages.
      if (db_has_name) {
        std::string *db_name = &block_name;
        if (get_use_generic_canonical_name()) {
          db_name = &alias;
        }
        if (alias != block_name)
          block->property_add(Ioss::Property("db_name", *db_name));
      }

      // Maintain block order on output database...
      block->property_add(Ioss::Property("original_block_order", used_blocks++));

      if (save_type != "null" && save_type != "") {
        if (block->property_exists("original_topology_type")) {
          if (block->get_property("original_topology_type").get_string() != save_type) {
            block->property_erase("original_topology_type");
            block->property_add(Ioss::Property("original_topology_type", save_type));
          }
        } else {
          if (block->get_property("topology_type").get_string() != save_type) {
            // Maintain original X type on output database if possible.
            block->property_add(Ioss::Property("original_topology_type", save_type));
          }
        }
      }

      block->property_add(Ioss::Property("global_entity_count", global_X_count[iblk]));


      offset += local_X_count[iblk];

      // See if this block is "omitted" by the calling code.
      // This only affects the generation of surfaces...
      if (!blockOmissions.empty()) {
        std::vector<std::string>::const_iterator I = std::find(blockOmissions.begin(),
                                                               blockOmissions.end(),
                                                               block_name);
        if (I != blockOmissions.end()) {
          block->property_add(Ioss::Property(std::string("omitted"), 1));
        }
      }

      get_region()->add_alias(block_name, alias);

      // Check for additional variables.
      add_attribute_fields(entity_type, block, attributes, type);
      add_results_fields(entity_type, block, iblk);

      if (entity_type == EX_ELEM_BLOCK) {
        Ioss::SerializeIO       serializeIO__(this);
	Ioex::add_map_fields(get_file_pointer(), (Ioss::ElementBlock*)block,
			     local_X_count[iblk], maximumNameLength);
      }
    }
    m_groupCount[entity_type] = used_blocks;
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
      Ioss::SerializeIO serializeIO__(this);
      for (int iblk = 0; iblk < m_groupCount[EX_ELEM_BLOCK]; iblk++) {
        Ioss::ElementBlock *eb = element_blocks[iblk];
        int blk_position =  eb->get_property("original_block_order").get_int();
        int64_t id =        eb->get_property("id").get_int();
        int element_nodes = eb->get_property("topology_node_count").get_int();
        int64_t my_element_count = eb->get_property("entity_count").get_int();
        if (my_element_count > 0) {
          if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
            std::vector<int64_t> conn(my_element_count * element_nodes);
            ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, TOPTR(conn), NULL, NULL);

            for (int64_t i=0; i < my_element_count * element_nodes; i++) {
              node_used[conn[i]-1] = blk_position+1;
            }
          } else {
            std::vector<int> conn(my_element_count * element_nodes);
            ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, TOPTR(conn), NULL, NULL);

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

#ifdef HAVE_MPI    
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
               << " in Iofx::DatabaseIO::compute_block_adjacencies";
        std::cerr << errmsg.str();
      }

      int local_error = (MPI_SUCCESS == result) ? 0 : 1 ;
      int global_error = util().global_minmax(local_error, Ioss::ParallelUtils::DO_MAX);

      if (global_error != 0) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Irecv error on some processor "
               << "in Iofx::DatabaseIO::compute_block_adjacencies";
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
               << " in Iofx::DatabaseIO::compute_block_adjacencies";
        std::cerr << errmsg.str();
      }

      local_error = (MPI_SUCCESS == result) ? 0 : 1 ;
      global_error = util().global_minmax(local_error, Ioss::ParallelUtils::DO_MAX);

      if (global_error != 0) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Rsend error on some processor "
               << "in Iofx::DatabaseIO::compute_block_adjacencies";
        IOSS_ERROR(errmsg);
      }

      result = MPI_Waitall(req_cnt, TOPTR(request), TOPTR(status));

      if (result != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Waitall error on processor " << util().parallel_rank()
               << " in Iofx::DatabaseIO::compute_block_adjacencies";
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
#endif    

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

#ifdef HAVE_MPI
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
#endif

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
      if (my_element_count > 0) {
        if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
          std::vector<int64_t> conn(my_element_count * element_nodes);
          ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, TOPTR(conn), NULL, NULL);
          for (int64_t j=0; j < my_element_count * element_nodes; j++) {
            nodeConnectivityStatus[conn[j]-1] |= status;
          }
        } else {
          std::vector<int> conn(my_element_count * element_nodes);
          ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, TOPTR(conn), NULL, NULL);
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

      Ioss::Int64Vector side_set_ids(m_groupCount[EX_SIDE_SET]);
      {
        Ioss::SerializeIO       serializeIO__(this);
        int error;
        if (ex_int64_status(get_file_pointer()) & EX_IDS_INT64_API) {
          error = ex_get_ids(get_file_pointer(), EX_SIDE_SET, TOPTR(side_set_ids));
        } else {
          Ioss::IntVector tmp_set_ids(side_set_ids.size());
          error = ex_get_ids(get_file_pointer(), EX_SIDE_SET, TOPTR(tmp_set_ids));
          if (error >= 0) std::copy(tmp_set_ids.begin(), tmp_set_ids.end(), side_set_ids.begin());
        }
        if (error < 0) {
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        }

        for (int i = 0; i < m_groupCount[EX_SIDE_SET]; i++) {
          std::vector<char> ss_name(maximumNameLength+1);
          error = ex_get_name(get_file_pointer(), EX_SIDE_SET, side_set_ids[i],
                              TOPTR(ss_name));
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
        int64_t id = side_set_ids[iss];
        std::string sid = "";
	Ioex::TopologyMap topo_map;
	Ioex::TopologyMap side_map; // Used to determine side consistency

        Ioss::SurfaceSplitType split_type = splitType;
        std::string side_set_name;
        Ioss::SideSet *side_set = NULL;

        {
          Ioss::SerializeIO     serializeIO__(this);

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
              if (alias != side_set_name)
                side_set->property_add(Ioss::Property("db_name", *db_name));
            }

            get_region()->add((Ioss::SideSet*)side_set);

            get_region()->add_alias(side_set_name, alias);
            get_region()->add_alias(side_set_name, Ioss::Utils::encode_entity_name("sideset", id));
          }

          //      split_type = SPLIT_BY_ELEMENT_BLOCK;
          //      split_type = SPLIT_BY_TOPOLOGIES;
          //      split_type = SPLIT_BY_DONT_SPLIT;

          // Determine how many side blocks compose this side set.
          ex_set set_param[1];
          set_param[0].id = id;
          set_param[0].type = EX_SIDE_SET;
          set_param[0].entry_list = NULL;
          set_param[0].extra_list = NULL;
          set_param[0].distribution_factor_list = NULL;

          int error = ex_get_sets(get_file_pointer(), 1, set_param);
          if (error < 0) {Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);}

          int64_t number_sides = set_param[0].num_entry;

          Ioss::Int64Vector element(number_sides);
          Ioss::Int64Vector sides(number_sides);

          // Easier below here is the element and sides are a know 64-bit size...
          // Kluge here to do that...
          if (int_byte_size_api() == 4) {
            Ioss::IntVector e32(number_sides);
            Ioss::IntVector s32(number_sides);
            int ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(e32), TOPTR(s32));
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            std::copy(e32.begin(), e32.end(), element.begin());
            std::copy(s32.begin(), s32.end(), sides.begin());
          } else {
            int ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(sides));
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
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
                           << "' Something is wrong in the Iofx::DatabaseIO class. Please report.\n";
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
                side_set->add(side_block);

                // Note that all sideblocks within a specific
                // sideset might have the same id.
                assert(side_block != NULL);
                side_block->property_add(Ioss::Property("id", id));

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
                  Ioss::SerializeIO serializeIO__(this);
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
          Ioss::SerializeIO     serializeIO__(this);

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

          if (number_sides > 0) {
            // Get the element and element side lists.
            if (int_byte_size_api() == 4) {
              Ioss::IntVector element(number_sides);
              Ioss::IntVector sides(number_sides);

              set_param[0].entry_list = TOPTR(element);
              set_param[0].extra_list = TOPTR(sides);
              ierr = ex_get_sets(get_file_pointer(), 1, set_param);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

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

              set_param[0].entry_list = TOPTR(element);
              set_param[0].extra_list = TOPTR(sides);
              ierr = ex_get_sets(get_file_pointer(), 1, set_param);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

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
        if (count > 0) {
          Ioss::Int64Vector Xset_ids(count);
          Ioss::IntVector attributes(count);
          std::vector<T*> Xsets(count);
          {
            Ioss::SerializeIO   serializeIO__(this);
            if (ex_int64_status(get_file_pointer()) & EX_IDS_INT64_API) {
              int error = ex_get_ids(get_file_pointer(), type, TOPTR(Xset_ids));
              if (error < 0) {Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);}
            } else {
              Ioss::IntVector tmp_set_ids(count);
              int error = ex_get_ids(get_file_pointer(), type, TOPTR(tmp_set_ids));
              if (error < 0) {Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);}
              std::copy(tmp_set_ids.begin(), tmp_set_ids.end(), Xset_ids.begin());
            }

            std::vector<ex_set> set_params(count);
            for (int ins = 0; ins < count; ins++) {
              set_params[ins].type = type;
              set_params[ins].id   = Xset_ids[ins];
              set_params[ins].entry_list = NULL;
              set_params[ins].extra_list = NULL;
              set_params[ins].distribution_factor_list = NULL;
            }

            int error = ex_get_sets(get_file_pointer(), count, TOPTR(set_params));
            if (error < 0) {
              Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }

          for (int ins = 0; ins < count; ins++) {
            int64_t id = set_params[ins].id;
            int num_attr = 0;
            int ierr = ex_get_attr_param(get_file_pointer(), type, id, &num_attr);
            if (ierr < 0)
              Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            attributes[ins] = num_attr;

            bool db_has_name = false;
            std::string Xset_name = Ioex::get_entity_name(get_file_pointer(), type, id, base+"list",
							  maximumNameLength, db_has_name);

            std::string alias = Ioss::Utils::encode_entity_name(base+"list", id);

            if (get_use_generic_canonical_name()) {
              std::string temp = Xset_name;
              Xset_name = alias;
              alias = temp;
            }

            bool filtered = false;
            int64_t original_set_size = set_params[ins].num_entry;
            Ioss::Int64Vector active_node_index;
            if (!blockOmissions.empty() && type == EX_NODE_SET) {
              active_node_index.resize(set_params[ins].num_entry);
              set_params[ins].entry_list = TOPTR(active_node_index);

              int old_status = ex_int64_status(get_file_pointer());
              ex_set_int64_status(get_file_pointer(), EX_BULK_INT64_API);             
              int error = ex_get_sets(get_file_pointer(), 1, &set_params[ins]);
	      if (error < 0) {Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);}
              ex_set_int64_status(get_file_pointer(), old_status);            

              compute_node_status();
              filtered = Ioex::filter_node_list(active_node_index, nodeConnectivityStatus);
              set_params[ins].num_entry = active_node_index.size();
            }
            T* Xset = new T(this, Xset_name, set_params[ins].num_entry);
            Xsets[ins] = Xset;
            Xset->property_add(Ioss::Property("id", id));
            if (db_has_name) {
              std::string *db_name = &Xset_name;
              if (get_use_generic_canonical_name()) {
                db_name = &alias;
              }
              if (alias != Xset_name)
                Xset->property_add(Ioss::Property("db_name", *db_name));
            }
            if (filtered && type == EX_NODE_SET) {
              Xset->property_add(Ioss::Property("filtered_db_set_size", original_set_size));
              activeNodesetNodesIndex[Xset_name].swap(active_node_index);
            }
            get_region()->add(Xset);
            get_region()->add_alias(Xset_name, alias);
            get_region()->add_alias(Xset_name, Ioss::Utils::encode_entity_name(base+"set",  id));
          }
        }

        // The attribute count will either be 0 if there are no
        // entities in the grouping entity on this processor, or it will be
        // the number of attributes (> 0). Therefore, if we take the 'max'
        // over all processors, each processor will then have the correct
        // attribute count...
        // This is a collective call...
        util().global_array_minmax(attributes, Ioss::ParallelUtils::DO_MAX);

        for (int ins = 0; ins < count; ins++) {
          add_attribute_fields(type, Xsets[ins], attributes[ins], "");
          add_results_fields(type, Xsets[ins], ins);
        }
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
        Ioss::SerializeIO       serializeIO__(this);
        // This is a parallel run. There should be communications data
        // Get nemesis commset metadata
        int64_t my_node_count = 0;
        int64_t elem_count = 0;

        // NOTE: It is possible for a parallel run to have no
        // communications maps if the decomposition occurs along contact
        // surfaces.  In this case, we create empty node and element
        // communication maps.
        if (commsetNodeCount > 0 || commsetElemCount > 0) {
          if (commsetNodeCount > 0) {
            nodeCmapIds.resize(commsetNodeCount);
            nodeCmapNodeCnts.resize(commsetNodeCount);
          }
          if (commsetElemCount > 0) {
            elemCmapIds.resize(commsetElemCount);
            elemCmapElemCnts.resize(commsetElemCount);
          }

          int error;
          if (int_byte_size_api() == 4) {
	    Ioss::IntVector nci(nodeCmapIds.size());
	    Ioss::IntVector ncnc(nodeCmapNodeCnts.size());
	    Ioss::IntVector eci(elemCmapIds.size());
	    Ioss::IntVector ecec(elemCmapElemCnts.size());
            error = ex_get_cmap_params(get_file_pointer(),
                                       TOPTR(nci), TOPTR(ncnc), TOPTR(eci), TOPTR(ecec),
                                       myProcessor);
            if (error >= 0) {
              std::copy(nci.begin(),  nci.end(),  nodeCmapIds.begin());
              std::copy(ncnc.begin(), ncnc.end(), nodeCmapNodeCnts.begin());
              std::copy(eci.begin(),  eci.end(),  elemCmapIds.begin());
              std::copy(ecec.begin(), ecec.end(), elemCmapElemCnts.begin());
            }
          } else {
            error = ex_get_cmap_params(get_file_pointer(),
                                       TOPTR(nodeCmapIds), TOPTR(nodeCmapNodeCnts),
                                       TOPTR(elemCmapIds), TOPTR(elemCmapElemCnts),
                                       myProcessor);
          }
          if (error < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

          // Count nodes, elements, and convert counts to offsets.
          //
          for (int64_t ics = 0; ics < commsetNodeCount; ics++) {
            my_node_count += nodeCmapNodeCnts[ics];
          }
          for (int64_t ecs = 0; ecs < commsetElemCount; ecs++) {
            elem_count += elemCmapElemCnts[ecs];
          }
        }
        // Create a single node commset and a single element commset
        Ioss::CommSet *commset = new Ioss::CommSet(this, "commset_node", "node",
                                                   my_node_count);
        commset->property_add(Ioss::Property("id", 1));
        get_region()->add(commset);

        commset = new Ioss::CommSet(this, "commset_side", "side",
                                    elem_count);
        commset->property_add(Ioss::Property("id", 1));
        get_region()->add(commset);
      }
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock* nb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      {
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);
        if (num_to_get > 0) {

#ifndef NDEBUG
          int64_t my_node_count = field.raw_count();
          assert(my_node_count == nodeCount);
#endif

          Ioss::Field::RoleType role = field.get_role();
          if (role == Ioss::Field::MESH) {
            if (field.get_name() == "mesh_model_coordinates_x") {
              double *rdata = static_cast<double*>(data);
              int ierr = ex_get_coord(get_file_pointer(), rdata, NULL, NULL);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }

            else if (field.get_name() == "mesh_model_coordinates_y") {
              double *rdata = static_cast<double*>(data);
              int ierr = ex_get_coord(get_file_pointer(), NULL, rdata, NULL);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }

            else if (field.get_name() == "mesh_model_coordinates_z") {
              double *rdata = static_cast<double*>(data);
              int ierr = ex_get_coord(get_file_pointer(), NULL, NULL, rdata);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }

            else if (field.get_name() == "mesh_model_coordinates") {
              // Data required by upper classes store x0, y0, z0, ... xn,
              // yn, zn. Data stored in exodusII file is x0, ..., xn, y0,
              // ..., yn, z0, ..., zn so we have to allocate some scratch
              // memory to read in the data and then map into supplied
              // 'data'
              std::vector<double> x(num_to_get);
              std::vector<double> y;
              if (spatialDimension > 1)
                y.resize(num_to_get);
              std::vector<double> z;
              if (spatialDimension == 3)
                z.resize(num_to_get);

              // Cast 'data' to correct size -- double
              double *rdata = static_cast<double*>(data);

              int ierr = ex_get_coord(get_file_pointer(), TOPTR(x), TOPTR(y), TOPTR(z));
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

              size_t index = 0;
              for (size_t i=0; i < num_to_get; i++) {
                rdata[index++] = x[i];
                if (spatialDimension > 1)
                  rdata[index++] = y[i];
                if (spatialDimension == 3)
                  rdata[index++] = z[i];
              }
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
              // If not parallel, then this is just 1..node_count
              // If parallel, then it is the data in the ex_get_id_map created by nem_spread.
              if (isParallel) {
                int error = ex_get_id_map(get_file_pointer(), EX_NODE_MAP, data);
                if (error < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              } else {
                if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
                  int64_t *idata = static_cast<int64_t*>(data);
                  for (int64_t i=0; i < nodeCount; i++) {
                    idata[i] = i+1;
                  }
                } else {
                  int *idata = static_cast<int*>(data);
                  for (int64_t i=0; i < nodeCount; i++) {
                    idata[i] = i+1;
                  }
                }
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
              if (isParallel) {
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
                    if (proc < myProcessor) {
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
                    if (proc < myProcessor) {
                      idata[node-1] = proc;
                    }
                  }
                }
              }
              else {
                // Serial case...
                if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
                  int64_t *idata = static_cast<int64_t*>(data);
                  for (int64_t i=0; i < nodeCount; i++) {
                    idata[i] = 0;
                  }
                } else {
                  int *idata = static_cast<int*>(data);
                  for (int64_t i=0; i < nodeCount; i++) {
                    idata[i] = 0;
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
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      {
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);
        if (num_to_get > 0) {

          int64_t id = Ioex::get_id(eb, EX_ELEM_BLOCK, &ids_);
          size_t my_element_count = eb->get_property("entity_count").get_int();
          Ioss::Field::RoleType role = field.get_role();

          if (role == Ioss::Field::MESH) {
            // Handle the MESH fields required for an ExodusII file model.
            // (The 'genesis' portion)

            if (field.get_name() == "connectivity") {
              int element_nodes = eb->get_property("topology_node_count").get_int();
              assert(field.raw_storage()->component_count() == element_nodes);

              // The connectivity is stored in a 1D array.
              // The element_node index varies fastet
              if (my_element_count > 0) {
                get_connectivity_data(get_file_pointer(), data, EX_ELEM_BLOCK, id, 0);
                get_map(EX_NODE_BLOCK).map_data(data, field, num_to_get*element_nodes);
              }
            }
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
            else if (field.get_name() == "connectivity_raw") {
              // "connectivity_raw" has nodes in local id space (1-based)
              assert(field.raw_storage()->component_count() == eb->get_property("topology_node_count").get_int());

              // The connectivity is stored in a 1D array.
              // The element_node index varies fastest
              if (my_element_count > 0) {
                get_connectivity_data(get_file_pointer(), data, EX_ELEM_BLOCK, id, 0);
              }
            }
            else if (field.get_name() == "ids") {
              // Map the local ids in this element block
              // (eb_offset+1...eb_offset+1+my_element_count) to global element ids.
              get_map(EX_ELEM_BLOCK).map_implicit_data(data, field, num_to_get, eb->get_offset());
            }
            else if (field.get_name() == "implicit_ids") {
              // If not parallel, then this is just one..element_count
              // If parallel, then it is the data in the ex_get_id_map created by nem_spread.
              size_t eb_offset_plus_one = eb->get_offset() + 1;
              if (isParallel) {
                int error = ex_get_partial_id_map(get_file_pointer(), EX_ELEM_MAP, eb_offset_plus_one, my_element_count, data);
                if (error < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              } else {
                if (ex_int64_status(get_file_pointer()) & EX_BULK_INT64_API) {
                  int64_t *idata = static_cast<int64_t*>(data);
                  for (size_t i=0; i < my_element_count; i++) {
                    idata[i] = eb_offset_plus_one + i;
                  }
                } else {
                  int *idata = static_cast<int*>(data);
                  for (size_t i=0; i < my_element_count; i++) {
                    idata[i] = eb_offset_plus_one + i;
                  }
                }
              }
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
                ex_get_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 1, eb_offset+1, my_element_count,
                                       TOPTR(element));
                ex_get_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 2, eb_offset+1, my_element_count,
                                       TOPTR(side));

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
                ex_get_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 1, eb_offset+1, my_element_count,
                                       TOPTR(element));
                ex_get_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 2, eb_offset+1, my_element_count,
                                       TOPTR(side));

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
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      {
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);
        if (num_to_get > 0) {

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
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      {
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);
        if (num_to_get > 0) {

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
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::get_Xset_field_internal(ex_entity_type type,
                                                const Ioss::EntitySet *ns,
                                                const Ioss::Field& field,
                                                void *data, size_t data_size) const
    {
      {
        int ierr;
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);
        if (num_to_get > 0) {

          int64_t id = Ioex::get_id(ns, type, &ids_);
          Ioss::Field::RoleType role = field.get_role();
          if (role == Ioss::Field::MESH) {

            if (field.get_name() == "ids" ||
                field.get_name() == "ids_raw") {
              if (field.get_type() == Ioss::Field::INTEGER) {
                ierr = ex_get_set(get_file_pointer(), type, id, static_cast<int*>(data), NULL);
              } else {
                ierr = ex_get_set(get_file_pointer(), type, id, static_cast<int64_t*>(data), NULL);
              }
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

              if (field.get_name() == "ids") {
                // Convert the local node ids to global ids
                get_map(type).map_data(data, field, num_to_get);
              }

            } else if (field.get_name() == "orientation") {
              if (field.get_type() == Ioss::Field::INTEGER) {
                ierr = ex_get_set(get_file_pointer(), type, id, NULL, static_cast<int*>(data));
              } else {
                ierr = ex_get_set(get_file_pointer(), type, id, NULL, static_cast<int64_t*>(data));
              }
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

            } else if (field.get_name() == "distribution_factors") {
              ex_set set_param[1];
              set_param[0].id = id;
              set_param[0].type = type;
              set_param[0].entry_list = NULL;
              set_param[0].extra_list = NULL;
              set_param[0].distribution_factor_list = NULL;
              ierr = ex_get_sets(get_file_pointer(), 1, set_param);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

              if (set_param[0].num_distribution_factor == 0) {
                double *rdata = static_cast<double*>(data);
                for (size_t i=0; i < num_to_get; i++)
                  rdata[i] = 1.0;
              } else {
                set_param[0].distribution_factor_list = static_cast<double*>(data);
                ierr = ex_get_sets(get_file_pointer(), 1, set_param);
                if (ierr < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              }
            } else {
              num_to_get = Ioss::Utils::field_warning(ns, field, "input");
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
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet* ns,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      if (!ns->property_exists("filtered_db_set_size")) {
        return get_Xset_field_internal(EX_NODE_SET, ns, field, data, data_size);
      }
      else {
        size_t db_size = ns->get_property("filtered_db_set_size").get_int();
        
        int ierr;
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);
        if (num_to_get > 0) {

          int64_t id = Ioex::get_id(ns, EX_NODE_SET, &ids_);
          Ioss::Field::RoleType role = field.get_role();
          if (role == Ioss::Field::MESH) {

            if (field.get_name() == "ids" ||
                field.get_name() == "ids_raw") {
              if (field.get_type() == Ioss::Field::INTEGER) {
		Ioss::IntVector dbvals(db_size);
                ierr = ex_get_set(get_file_pointer(), EX_NODE_SET, id, TOPTR(dbvals), NULL);
                if (ierr >= 0)
		  Ioex::filter_node_list(static_cast<int*>(data), dbvals, activeNodesetNodesIndex[ns->name()]);
              } else {
		Ioss::Int64Vector dbvals(db_size);
                ierr = ex_get_set(get_file_pointer(), EX_NODE_SET, id, TOPTR(dbvals), NULL);
                if (ierr >= 0)
		  Ioex::filter_node_list(static_cast<int64_t*>(data), dbvals, activeNodesetNodesIndex[ns->name()]);
              }
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

              if (field.get_name() == "ids") {
                // Convert the local node ids to global ids
                get_map(EX_NODE_SET).map_data(data, field, num_to_get);
              }

            } else if (field.get_name() == "distribution_factors") {
              ex_set set_param[1];
              set_param[0].id = id;
              set_param[0].type = EX_NODE_SET;
              set_param[0].entry_list = NULL;
              set_param[0].extra_list = NULL;
              set_param[0].distribution_factor_list = NULL;
              ierr = ex_get_sets(get_file_pointer(), 1, set_param);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

              if (set_param[0].num_distribution_factor == 0) {
                double *rdata = static_cast<double*>(data);
                for (size_t i=0; i < num_to_get; i++)
                  rdata[i] = 1.0;
              } else {
                std::vector<double>dbvals(db_size);
                set_param[0].distribution_factor_list = TOPTR(dbvals);
                ierr = ex_get_sets(get_file_pointer(), 1, set_param);
                if (ierr < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
		Ioex::filter_node_list(static_cast<double*>(data), dbvals, activeNodesetNodesIndex[ns->name()]);
              }
            } else {
              num_to_get = Ioss::Utils::field_warning(ns, field, "input");
            }
          } else if (role == Ioss::Field::ATTRIBUTE) {
            num_to_get = Ioss::Utils::field_warning(ns, field, "input");

          } else if (role == Ioss::Field::TRANSIENT) {
            // Filtered not currently implemented for transient or attributes....
            num_to_get = Ioss::Utils::field_warning(ns, field, "input");
          }
        }
        return num_to_get;
      }
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
      {
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);

        if (num_to_get > 0) {
          int64_t entity_count = cs->get_property("entity_count").get_int();

          // Return the <entity (node or side), processor> pair
          if (field.get_name() == "entity_processor" || field.get_name() == "entity_processor_raw") {

            // Check type -- node or side
            std::string type = cs->get_property("entity_type").get_string();

            // Allocate temporary storage space
            std::vector<char> entities(num_to_get * int_byte_size_api());
            std::vector<char> procs(num_to_get    * int_byte_size_api());

            if (type == "node") {
              int64_t cm_offset = 0;

              for (int64_t i=0; i < commsetNodeCount; i++) {
                int ierr = ex_get_node_cmap(get_file_pointer(), nodeCmapIds[i],
                                            &entities[cm_offset], &procs[cm_offset],
                                            myProcessor);
                if (ierr < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
                cm_offset += (nodeCmapNodeCnts[i] * int_byte_size_api());
              }
              assert(cm_offset == entity_count * int_byte_size_api());

              // Convert local node id to global node id and store in 'data'
              if (int_byte_size_api() == 4) {
                int* entity_proc = static_cast<int*>(data);
                int* ents = reinterpret_cast<int*>(&entities[0]);
                int* pros = reinterpret_cast<int*>(&procs[0]);

                size_t j=0;
                if (field.get_name() == "entity_processor") {
                  const Ioss::MapContainer &map = get_map(EX_NODE_BLOCK).map;

                  for (int64_t i=0; i < entity_count; i++) {
                    int local_id = ents[i];
                    entity_proc[j++] = map[local_id];
                    entity_proc[j++] = pros[i];
                  }
                } else {
                  for (int64_t i=0; i < entity_count; i++) {
                    entity_proc[j++] = ents[i];
                    entity_proc[j++] = pros[i];
                  }
                }
              } else {
                int64_t* entity_proc = static_cast<int64_t*>(data);
                int64_t* ents = reinterpret_cast<int64_t*>(&entities[0]);
                int64_t* pros = reinterpret_cast<int64_t*>(&procs[0]);

                size_t j=0;
                if (field.get_name() == "entity_processor") {
                  const Ioss::MapContainer &map = get_map(EX_NODE_BLOCK).map;

                  for (int64_t i=0; i < entity_count; i++) {
                    int64_t local_id = ents[i];
                    entity_proc[j++] = map[local_id];
                    entity_proc[j++] = pros[i];
                  }
                } else {
                  for (int64_t i=0; i < entity_count; i++) {
                    entity_proc[j++] = ents[i];
                    entity_proc[j++] = pros[i];
                  }
                }
              }
            }
            else if (type == "side") {
              std::vector<char> sides(entity_count * int_byte_size_api());
              int64_t cm_offset = 0;
              for (int64_t i=0; i < commsetElemCount; i++) {
                int ierr = ex_get_elem_cmap(get_file_pointer(), elemCmapIds[i],
                                            &entities[cm_offset], &sides[cm_offset],
                                            &procs[cm_offset], myProcessor);
                if (ierr < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
                cm_offset += (elemCmapElemCnts[i] * int_byte_size_api());
              }
              assert(cm_offset == entity_count * int_byte_size_api());

              if (int_byte_size_api() == 4) {
                int* entity_proc = static_cast<int*>(data);
                int* ents = reinterpret_cast<int*>(&entities[0]);
                int* pros = reinterpret_cast<int*>(&procs[0]);
                int* sids = reinterpret_cast<int*>(&sides[0]);

                size_t j=0;
                if (field.get_name() == "entity_processor") {
                  const Ioss::MapContainer &map = get_map(EX_ELEM_BLOCK).map;

                  for (ssize_t i=0; i < entity_count; i++) {
                    int local_id = ents[i];
                    entity_proc[j++] = 10*map[local_id]+sids[i];
                    entity_proc[j++] = pros[i];
                  }
                } else {
                  for (ssize_t i=0; i < entity_count; i++) {
                    entity_proc[j++] = 10*ents[i]+sids[i];
                    entity_proc[j++] = pros[i];
                  }
                }
              } else {
                int64_t* entity_proc = static_cast<int64_t*>(data);
                int64_t* ents = reinterpret_cast<int64_t*>(&entities[0]);
                int64_t* pros = reinterpret_cast<int64_t*>(&procs[0]);
                int64_t* sids = reinterpret_cast<int64_t*>(&sides[0]);

                size_t j=0;
                if (field.get_name() == "entity_processor") {
                  const Ioss::MapContainer &map = get_map(EX_ELEM_BLOCK).map;

                  for (ssize_t i=0; i < entity_count; i++) {
                    int64_t local_id = ents[i];
                    entity_proc[j++] = 10*map[local_id]+sids[i];
                    entity_proc[j++] = pros[i];
                  }
                } else {
                  for (ssize_t i=0; i < entity_count; i++) {
                    entity_proc[j++] = 10*ents[i]+sids[i];
                    entity_proc[j++] = pros[i];
                  }
                }
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
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock* fb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      Ioss::SerializeIO serializeIO__(this);
      ssize_t num_to_get = field.verify(data_size);
      if (num_to_get > 0) {

        int64_t id = Ioex::get_id(fb, EX_SIDE_SET, &ids_);
        int64_t entity_count = fb->get_property("entity_count").get_int();
        if (num_to_get != entity_count) {
          std::ostringstream errmsg;
          errmsg << "ERROR: Partial field input not yet implemented for side blocks";
          IOSS_ERROR(errmsg);
        }

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
        int64_t number_distribution_factors = set_param[0].num_distribution_factor;

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
              set_param[0].distribution_factor_list = TOPTR(real_ids);
              ierr = ex_get_sets(get_file_pointer(), 1, set_param);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

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

          } else if (field.get_name() == "element_side") {
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

            std::vector<char> element(number_sides * int_byte_size_api());
            std::vector<char> sides(number_sides * int_byte_size_api());

            ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(sides));
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

            if (number_sides == entity_count) {
              ssize_t index = 0;
              if (int_byte_size_api() == 4) {
                int     *element_side = static_cast<int*>(data);
                int     *element32 = (int*)TOPTR(element);
                int     *sides32 = (int*)TOPTR(sides);
                for (ssize_t iel = 0; iel < entity_count; iel++) {
                  element_side[index++] = map[element32[iel]];
                  element_side[index++] = sides32[iel] - side_offset;
                }
              } else {
                int64_t *element_side = static_cast<int64_t*>(data);
                int64_t *element64 = (int64_t*)TOPTR(element);
                int64_t *sides64 = (int64_t*)TOPTR(sides);
                for (ssize_t iel = 0; iel < entity_count; iel++) {
                  element_side[index++] = map[element64[iel]];
                  element_side[index++] = sides64[iel] - side_offset;
                }
              }
              assert(index/2 == entity_count);
            } else {
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
                    element_side[index++] = map[element32[iel]];
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
                    element_side[index++] = map[element64[iel]];
                    element_side[index++] = sides64[iel] - side_offset;
                  }
                }
              }
              assert(index/2 == entity_count);
            }

          } else if (field.get_name() == "element_side_raw") {
            // In exodusII, the 'side set' is stored as a sideset.  A sideset
            // has a list of elements and a corresponding local element side
            // (1-based)

            // Since we only have a single array, we need to allocate an extra
            // array to store all of the data.  Note also that the
            // element_id for the "raw" field is the local id, not the
            // global id.

            // See if edges or faces...
            int64_t side_offset = Ioss::Utils::get_side_offset(fb);

            std::vector<char> element(number_sides * int_byte_size_api());
            std::vector<char> sides(number_sides * int_byte_size_api());

            ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(sides));
            if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

            if (number_sides == entity_count) {
              ssize_t index = 0;
              if (int_byte_size_api() == 4) {
                int     *element_side = static_cast<int*>(data);
                int     *element32 = (int*)TOPTR(element);
                int     *sides32 = (int*)TOPTR(sides);
                for (ssize_t iel = 0; iel < entity_count; iel++) {
                  element_side[index++] = element32[iel];
                  element_side[index++] = sides32[iel] - side_offset;
                }
              } else {
                int64_t *element_side = static_cast<int64_t*>(data);
                int64_t *element64 = (int64_t*)TOPTR(element);
                int64_t *sides64 = (int64_t*)TOPTR(sides);
                for (ssize_t iel = 0; iel < entity_count; iel++) {
                  element_side[index++] = element64[iel];
                  element_side[index++] = sides64[iel] - side_offset;
                }
              }
              assert(index/2 == entity_count);
            } else {
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

          if (number_sides == entity_count) {
            num_to_get = read_transient_field(EX_SIDE_SET, m_variables[EX_SIDE_SET], field, fb, data);
          } else {
            // Need to read all values for the specified field and then
            // filter down to the elements actualy in this side block.

            // Determine which sides are member of this block
            Ioss::IntVector is_valid_side;
            {
              //----
              std::vector<char> element(number_sides * int_byte_size_api());
              std::vector<char> sides(number_sides * int_byte_size_api());
              ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(sides));
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              //----
              Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, int_byte_size_api(),
                                                          TOPTR(element), TOPTR(sides),
                                                          number_sides, get_region());
            }

            num_to_get = read_ss_transient_field(field, id, data, is_valid_side);
          }
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

      if (offset == 1 && field.raw_storage()->component_count() == attribute_count) {
        // Write all attributes in one big chunk...
        int ierr = ex_put_attr(get_file_pointer(), type, id, static_cast<double*>(data));
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
      else {
        // Write a subset of the attributes.  If scalar, write one;
        // if higher-order (vector3d, ..) write each component.
        if (field.raw_storage()->component_count() == 1) {
          int ierr = ex_put_one_attr(get_file_pointer(), type, id,
                                     offset, static_cast<double*>(data));
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        } else {
          // Multi-component...  Need a local memory space to push
          // data into and then write that out to the file...
          std::vector<double> local_data(num_entity);
          int comp_count = field.raw_storage()->component_count();
          double *rdata = static_cast<double*>(data);
          for (int i=0; i < comp_count; i++) {
            size_t k = i;
            for (ssize_t j=0; j < num_entity; j++) {
              local_data[j] = rdata[k];
              k += comp_count;
            }

            int ierr = ex_put_one_attr(get_file_pointer(), type, id,
                                       offset+i, TOPTR(local_data));
            if (ierr < 0)
              Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }
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
      if (num_entity == 0)
        return 0;

      int attribute_count = ge->get_property("attribute_count").get_int();
      int64_t id = Ioex::get_id(ge, type, &ids_);

      std::string att_name = ge->name() + SEP() + field.get_name();
      ssize_t offset = field.get_index();
      assert(offset-1+field.raw_storage()->component_count() <= attribute_count);
      if (offset == 1 && field.raw_storage()->component_count() == attribute_count) {
        // Read all attributes in one big chunk...
        int ierr = ex_get_attr(get_file_pointer(), type, id, static_cast<double*>(data));
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
      else {
        // Read a subset of the attributes.  If scalar, read one;
        // if higher-order (vector3d, ..) read each component and
        // put into correct location...
        if (field.raw_storage()->component_count() == 1) {
          int ierr = ex_get_one_attr(get_file_pointer(), type, id,
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
            int ierr = ex_get_one_attr(get_file_pointer(), type, id,
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
        ierr = ex_get_var(get_file_pointer(), step, type,
                          var_index, id, num_entity, TOPTR(temp));
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
        ierr = ex_get_var(get_file_pointer(), step, EX_SIDE_SET,
                          var_index, id, my_side_count, TOPTR(temp));
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

    template <typename INT>
      int64_t DatabaseIO::get_side_connectivity_internal(const Ioss::SideBlock* fb,
                                                         int64_t id, int64_t,
                                                         INT *fconnect,
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
        //----
        std::vector<INT> element(number_sides);
        std::vector<INT> side(number_sides);

        set_param[0].entry_list = TOPTR(element);
        set_param[0].extra_list = TOPTR(side);
        ierr = ex_get_sets(get_file_pointer(), 1, set_param);
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
        //----

        Ioss::IntVector is_valid_side;
        Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, int_byte_size_api(),
                                                    (void*)TOPTR(element), (void*)TOPTR(side),
                                                    number_sides, get_region());

        std::vector<INT> elconnect;
        int64_t elconsize = 0; // Size of currently allocated connectivity block
        Ioss::ElementBlock *conn_block = NULL; // Block that we currently
        // have connectivity for

        Ioss::ElementBlock *block = NULL;
        Ioss::IntVector side_elem_map; // Maps the side into the elements

        // connectivity array
        int64_t current_side = -1;
        int nelnode = 0;
        int nfnodes = 0;
        int ieb = 0;
        size_t offset = 0;
        for (ssize_t iel = 0; iel < number_sides; iel++) {
          if (is_valid_side[iel] == 1) {

            int64_t elem_id = element[iel];

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
                elconnect.resize(elconsize);
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
            int64_t side_id = side[iel];

            if (current_side != side_id) {
              side_elem_map = block->topology()->boundary_connectivity(side_id);
              current_side = side_id;
              nfnodes = block->topology()->boundary_type(side_id)->number_nodes();
            }
            for (int inode = 0; inode < nfnodes; inode++) {
              size_t index = (elem_id-offset)*nelnode + side_elem_map[inode];
              fconnect[ieb++] = elconnect[index];
            }
          }
        }
        return ierr;
      }

    int64_t DatabaseIO::get_side_connectivity(const Ioss::SideBlock* fb,
                                              int64_t id, int64_t my_side_count,
                                              void *fconnect,
                                              bool map_ids) const
    {
      if (int_byte_size_api() == 4) {
        return get_side_connectivity_internal(fb, id, my_side_count, (int*)fconnect, map_ids);
      } else {
        return get_side_connectivity_internal(fb, id, my_side_count, (int64_t*)fconnect, map_ids);
      }
    }

    // Get distribution factors for the specified side block
    int64_t DatabaseIO::get_side_distributions(const Ioss::SideBlock* fb,
                                               int64_t id, int64_t my_side_count,
                                               double *dist_fact,
                                               size_t /* data_size */) const
    {
      // Allocate space for elements and local side numbers
      // Get size of data stored on the file...
      ex_set set_param[1];
      set_param[0].id = id;
      set_param[0].type = EX_SIDE_SET;
      set_param[0].entry_list = NULL;
      set_param[0].extra_list = NULL;
      set_param[0].distribution_factor_list = NULL;

      int error = ex_get_sets(get_file_pointer(), 1, set_param);
      if (error < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      int64_t number_sides = set_param[0].num_entry;
      int64_t number_distribution_factors = set_param[0].num_distribution_factor;

      const Ioss::ElementTopology *ftopo = fb->topology();
      int nfnodes = ftopo->number_nodes();

      if (number_distribution_factors == 0) {
        // Fill in the array with '1.0'...
        for (int64_t i=0; i < nfnodes * my_side_count; i++)
          dist_fact[i] = 1.0;
        return 0;
      }

      // Take care of the easy situation -- If 'side_count' ==
      // 'number_sides' then the sideset is stored in a single sideblock
      // and all distribution factors on the database are transferred
      // 1-to-1 into 'dist_fact' array.
      if (my_side_count == number_sides) {
        return ex_get_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id, dist_fact);
      }

      // Allocate space for distribution factors.
      std::vector<double> dist(number_distribution_factors);
      int ierr = ex_get_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id, TOPTR(dist));
      if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // Another easy situation (and common for exodusII) is if the input
      // distribution factors are all the same value (typically 1).  In
      // that case, we only have to fill in the output array with that
      // value.
      {
        double value = dist[0];
        bool constant = true;
        for (int64_t i=1; i < number_distribution_factors; i++) {
          if (dist[i] != value) {
            constant = false;
            break;
          }
        }
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
      std::vector<char> side(number_sides * int_byte_size_api());

      ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(side));
      if (ierr < 0)
        Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      //----

      Ioss::IntVector is_valid_side;
      Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, int_byte_size_api(),
                                                  TOPTR(element), TOPTR(side),
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
        side32 = (int*)TOPTR(side);
      } else {
        element64 = (int64_t*)TOPTR(element);
        side64 = (int64_t*)TOPTR(side);
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
                 << "Something is wrong in the Iofx::DatabaseIO class. Please report.\n";
          IOSS_ERROR(errmsg);
        }

        const Ioss::ElementTopology *topo = block->topology()->boundary_type(side_id);

        if (topo == NULL) {
          std::ostringstream errmsg;
          errmsg << "INTERNAL ERROR: Could not find topology of element block boundary. "
                 << "Something is wrong in the Iofx::DatabaseIO class. Please report.\n";
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
      {
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);
        if (num_to_get > 0) {

          Ioss::Field::RoleType role = field.get_role();

          if (role == Ioss::Field::MESH) {
            if (field.get_name() == "mesh_model_coordinates_x") {
              double *rdata = static_cast<double*>(data);
              int ierr = ex_put_coord(get_file_pointer(), rdata, NULL, NULL);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }

            else if (field.get_name() == "mesh_model_coordinates_y") {
              double *rdata = static_cast<double*>(data);
              int ierr = ex_put_coord(get_file_pointer(), NULL, rdata, NULL);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }

            else if (field.get_name() == "mesh_model_coordinates_z") {
              double *rdata = static_cast<double*>(data);
              int ierr = ex_put_coord(get_file_pointer(), NULL, NULL, rdata);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }

            else if (field.get_name() == "mesh_model_coordinates") {
              // Data required by upper classes store x0, y0, z0, ... xn, yn, zn
              // Data stored in exodusII file is x0, ..., xn, y0, ..., yn, z0, ..., zn
              // so we have to allocate some scratch memory to read in the data
              // and then map into supplied 'data'
              std::vector<double> x; x.reserve(num_to_get);
              std::vector<double> y;
              if (spatialDimension > 1)
                y.reserve(num_to_get);
              std::vector<double> z;
              if (spatialDimension == 3)
                z.reserve(num_to_get);

              // Cast 'data' to correct size -- double
              double *rdata = static_cast<double*>(data);

              size_t index = 0;
              for (size_t i=0; i < num_to_get; i++) {
                x.push_back(rdata[index++]);
                if (spatialDimension > 1)
                  y.push_back(rdata[index++]);
                if (spatialDimension == 3)
                  z.push_back(rdata[index++]);
              }
              int ierr = ex_put_coord(get_file_pointer(), TOPTR(x), TOPTR(y), TOPTR(z));
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

            } else if (field.get_name() == "ids") {
              // The ids coming in are the global ids; their position is the
              // local id -1 (That is, data[0] contains the global id of local
              // node 1)
              handle_node_ids(data, num_to_get);
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

          } else if (role == Ioss::Field::ATTRIBUTE) {
            num_to_get = write_attribute_field(EX_NODE_BLOCK, field, nb, data);
          }
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      {
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);

        if (num_to_get > 0) {
          int ierr = 0;

          // Get the element block id and element count
          int64_t id = Ioex::get_id(eb, EX_ELEM_BLOCK, &ids_);
          size_t my_element_count = eb->get_property("entity_count").get_int();
          Ioss::Field::RoleType role = field.get_role();

          if (role == Ioss::Field::MESH) {
            // Handle the MESH fields required for an ExodusII file model.
            // (The 'genesis' portion)
            if (field.get_name() == "connectivity") {
              if (my_element_count > 0) {
                // Map element connectivity from global node id to local node id.
                int element_nodes = eb->get_property("topology_node_count").get_int();
                nodeMap.reverse_map_data(data, field, num_to_get*element_nodes);
                ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, data, NULL, NULL);
                if (ierr < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              }
            } else if (field.get_name() == "connectivity_edge") {
              if (my_element_count > 0) {
                // Map element connectivity from global edge id to local edge id.
                int element_edges = field.transformed_storage()->component_count();
                edgeMap.reverse_map_data(data, field, num_to_get*element_edges);
                ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, NULL, data, NULL);
                if (ierr < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              }
            } else if (field.get_name() == "connectivity_face") {
              if (my_element_count > 0) {
                // Map element connectivity from global face id to local face id.
                int element_faces = field.transformed_storage()->component_count();
                faceMap.reverse_map_data(data, field, num_to_get*element_faces);
                ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, NULL, NULL, data);
                if (ierr < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              }
            } else if (field.get_name() == "connectivity_raw") {
              if (my_element_count > 0) {
                // Element connectivity is already in local node id.
                ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, data, NULL, NULL);
                if (ierr < 0)
                  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              }
            } else if (field.get_name() == "ids") {
              handle_element_ids(eb, data, num_to_get);
            }
            else if (field.get_name() == "implicit_ids") {
              // Do nothing, input only field.
            }
            else if (field.get_name() == "skin") {
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
                for (size_t i=0; i < my_element_count; i++) {
                  element32[i] = el_side[index++];
                  side32[i]    = el_side[index++];
                }
              } else {
                int64_t *el_side = (int64_t *)data;
                int64_t *element64 = (int64_t*)TOPTR(element);
                int64_t *side64 = (int64_t*)TOPTR(side);

                int64_t index = 0;
                for (size_t i=0; i < my_element_count; i++) {
                  element64[i] = el_side[index++];
                  side64[i]    = el_side[index++];
                }
              }

              size_t eb_offset = eb->get_offset();
              ierr = ex_put_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 1, eb_offset+1, my_element_count,
                                            TOPTR(element));
              if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

              ierr = ex_put_partial_num_map(get_file_pointer(), EX_ELEM_MAP, 2, eb_offset+1, my_element_count,
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
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      {
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);

        if (num_to_get > 0) {
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
                // Do it in 'data' ...
                int face_nodes = eb->get_property("topology_node_count").get_int();
                nodeMap.reverse_map_data(data, field, num_to_get*face_nodes);
                ierr = ex_put_conn(get_file_pointer(), EX_FACE_BLOCK, id, data, NULL, NULL);
                if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
              }
            } else if (field.get_name() == "connectivity_edge") {
              if (my_face_count > 0) {
                // Map face connectivity from global edge id to local edge id.
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
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock* eb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      {
        Ioss::SerializeIO       serializeIO__(this);

        size_t num_to_get = field.verify(data_size);

        if (num_to_get > 0) {
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
                // Do it in 'data' ...
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
        }
        return num_to_get;
      }
    }

    int64_t DatabaseIO::handle_node_ids(void* ids, int64_t num_to_get) const
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
      assert(num_to_get == nodeCount);

      if (dbState == Ioss::STATE_MODEL) {
        if (nodeMap.map.empty()) {
          nodeMap.map.resize(nodeCount+1);
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

        // Only a single nodeblock and all set
        if (num_to_get == nodeCount) {
          assert(nodeMap.map[0] == -1 || nodeMap.reverse.size() == (size_t)nodeCount);
        }
        assert(get_region()->get_property("node_block_count").get_int() == 1);

        // Write to the database...
        int ierr = ex_put_id_map(get_file_pointer(), EX_NODE_MAP, ids);
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }

      nodeMap.build_reorder_map(0, num_to_get);
      return num_to_get;
    }

    namespace {
      size_t handle_block_ids(const Ioss::EntityBlock *eb,
                              ex_entity_type map_type,
                              Ioss::State db_state,
                              Ioss::Map &entity_map,
                              void* ids, size_t int_byte_size, size_t num_to_get, int file_pointer, int my_processor)
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
          int ierr = ex_put_partial_id_map(file_pointer, map_type, eb_offset+1, num_to_get, ids);
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

    int64_t DatabaseIO::handle_element_ids(const Ioss::ElementBlock *eb, void* ids, size_t num_to_get) const
    {
      if (elemMap.map.empty()) {
        elemMap.map.resize(elementCount+1);
        elemMap.map[0] = -1;
      }
      return handle_block_ids(eb, EX_ELEM_MAP, dbState, elemMap,
                              ids, int_byte_size_api(), num_to_get, get_file_pointer(), myProcessor);
    }

    int64_t DatabaseIO::handle_face_ids(const Ioss::FaceBlock *eb, void* ids, size_t num_to_get) const
    {
      if (faceMap.map.empty()) {
        faceMap.map.resize(faceCount+1);
        faceMap.map[0] = -1;
      }
      return handle_block_ids(eb, EX_FACE_MAP, dbState, faceMap,
                              ids, int_byte_size_api(), num_to_get, get_file_pointer(), myProcessor);
    }

    int64_t DatabaseIO::handle_edge_ids(const Ioss::EdgeBlock *eb, void* ids, size_t num_to_get) const
    {
      if (edgeMap.map.empty()) {
        edgeMap.map.resize(edgeCount+1);
        edgeMap.map[0] = -1;
      }
      return handle_block_ids(eb, EX_EDGE_MAP, dbState, edgeMap,
                              ids, int_byte_size_api(), num_to_get, get_file_pointer(), myProcessor);
    }

    void DatabaseIO::write_nodal_transient_field(ex_entity_type /* type */,
                                                 const Ioss::Field &field,
                                                 const Ioss::NodeBlock */* ge */,
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
          int64_t num_out = 0;

          if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX)
            num_out = nodeMap.map_field_to_db_scalar_order(static_cast<double*>(variables),
                                                           temp, begin_offset, count, stride, 0);
          else if (ioss_type == Ioss::Field::INTEGER)
            num_out = nodeMap.map_field_to_db_scalar_order(static_cast<int*>(variables),
                                                           temp, begin_offset, count, stride, 0);
          else if (ioss_type == Ioss::Field::INT64)
            num_out = nodeMap.map_field_to_db_scalar_order(static_cast<int64_t*>(variables),
                                                           temp, begin_offset, count, stride, 0);

          if (num_out != nodeCount) {
            std::ostringstream errmsg;
            errmsg << "ERROR: Problem outputting nodal variable '" << var_name
                   << "' with index = " << var_index << " to file "
                   << util().decode_filename(get_filename(), isParallel) << "\n"
                   << "Should have output " << nodeCount << " values, but instead only output " << num_out << " values.\n";
            IOSS_ERROR(errmsg);
          }

          // Write the variable...
          int ierr = ex_put_var(get_file_pointer(), step, EX_NODE_BLOCK, var_index, 0, num_out, TOPTR(temp));
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
          int64_t id = Ioex::get_id(ge, type, &ids_);
          int ierr;
          if (type == EX_SIDE_SET) {
            size_t offset = ge->get_property("set_offset").get_int();
            ierr = ex_put_partial_var(get_file_pointer(), step, type, var_index, id, offset+1, count, TOPTR(temp));
          } else {
            ierr = ex_put_var(get_file_pointer(), step, type, var_index, id, count, TOPTR(temp));
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
      {
        Ioss::SerializeIO       serializeIO__(this);
        ex_update(get_file_pointer());

        size_t entity_count = ns->get_property("entity_count").get_int();
        size_t num_to_get = field.verify(data_size);
        if (num_to_get > 0) {

          int64_t id = Ioex::get_id(ns, type, &ids_);
          Ioss::Field::RoleType role = field.get_role();

          if (role == Ioss::Field::MESH) {

            if (field.get_name() == "ids" ||
                field.get_name() == "ids_raw") {
              // Map node id from global node id to local node id.
              // Do it in 'data' ...

              if (field.get_name() == "ids") {
                nodeMap.reverse_map_data(data, field, num_to_get);
              }
              int ierr = ex_put_set(get_file_pointer(), type, id, data, NULL);
              if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

            } else if (field.get_name() == "orientation") {
              int ierr = ex_put_set(get_file_pointer(), type, id, NULL, data);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

            } else if (field.get_name() == "distribution_factors") {
              int ierr = ex_put_set_dist_fact(get_file_pointer(), type, id, static_cast<double*>(data));
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
        }
        return num_to_get;
      }
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

    int64_t DatabaseIO::put_field_internal(const Ioss::CommSet* cs,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      size_t num_to_get = field.verify(data_size);
      size_t entity_count = cs->get_property("entity_count").get_int();

      assert(num_to_get == entity_count);
      if (num_to_get == 0)
        return 0;

      // Return the <entity (node or side), processor> pair
      if (field.get_name() == "entity_processor") {

        // Check type -- node or side
        std::string type = cs->get_property("entity_type").get_string();

        // Allocate temporary storage space
        std::vector<char> entities(entity_count * int_byte_size_api());
        std::vector<char> procs(entity_count * int_byte_size_api());

        if (type == "node") {
          Ioss::SerializeIO     serializeIO__(this);
          // Convert global node id to local node id and store in 'entities'
          if (int_byte_size_api() == 4) {
            int* entity_proc = static_cast<int*>(data);
            int* ent = (int*)&entities[0];
            int* pro = (int*)&procs[0];
            int j=0;
            for (size_t i=0; i < entity_count; i++) {
              int global_id = entity_proc[j++];
              ent[i] = nodeMap.global_to_local(global_id, true);
              pro[i] = entity_proc[j++];
            }
          } else {
            int64_t* entity_proc = static_cast<int64_t*>(data);
            int64_t* ent = (int64_t*)&entities[0];
            int64_t* pro = (int64_t*)&procs[0];
            int64_t j=0;
            for (size_t i=0; i < entity_count; i++) {
              int64_t global_id = entity_proc[j++];
              ent[i] = nodeMap.global_to_local(global_id, true);
              pro[i] = entity_proc[j++];
            }
          }

          if (commsetNodeCount > 0) {
            int ierr = ex_put_node_cmap(get_file_pointer(), Ioex::get_id(cs, (ex_entity_type)0, &ids_),
                                        TOPTR(entities), TOPTR(procs), myProcessor);
            if (ierr < 0)
              Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }

          if (commsetNodeCount == 1) {
            // NOTE: The internal and border node maps must be output in one call.
            //       In this routine, we only have one commset at a time and can't
            //       construct the entire map at one time.  This is not really needed,
            //       so for now we just skip if there is more than one commset.  If
            //       this information is really needed, need to cache the information
            //       until all commsets have been processed.  Also need to change
            //       write_communication_metada() [Maybe, unless client sets correct
            //       properties.]

            // Construct the node map (internal vs. border).
            // Border nodes are those in the communication map (use entities array)
            // Internal nodes are the rest.  Allocate array to hold all nodes,
            // initialize all to '1', then zero out the nodes in 'entities'.
            // Iterate through array again and consolidate all '1's

            std::vector<char> internal(nodeCount * int_byte_size_api());
            if (int_byte_size_api() == 4) {
              compute_internal_border_maps((int*)&entities[0], (int*)&internal[0], nodeCount, entity_count);
            } else {
              compute_internal_border_maps((int64_t*)&entities[0], (int64_t*)&internal[0], nodeCount, entity_count);
            }

            int ierr = ex_put_processor_node_maps(get_file_pointer(), TOPTR(internal), TOPTR(entities), NULL,
                myProcessor);
            if (ierr < 0)
              Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }

        } else if (type == "side") {
          Ioss::SerializeIO     serializeIO__(this);
          std::vector<char> sides(entity_count * int_byte_size_api());
          if (int_byte_size_api() == 4) {
            int* entity_proc = static_cast<int*>(data);
            int* ent = (int*)&entities[0];
            int* sid = (int*)&sides[0];
            int* pro = (int*)&procs[0];
            int j=0;
            for (size_t i=0; i < entity_count; i++) {
              // Assume klugy side id generation.
              int global_id = entity_proc[j] / 10;
              ent[i] = elemMap.global_to_local(global_id);
              sid[i] = entity_proc[j++] % 10;
              pro[i] = entity_proc[j++];
            }
          } else {
            int64_t* entity_proc = static_cast<int64_t*>(data);
            int64_t* ent = (int64_t*)&entities[0];
            int64_t* sid = (int64_t*)&sides[0];
            int64_t* pro = (int64_t*)&procs[0];
            int64_t j=0;
            for (size_t i=0; i < entity_count; i++) {
              // Assume klugy side id generation.
              int64_t global_id = entity_proc[j] / 10;
              ent[i] = elemMap.global_to_local(global_id);
              sid[i] = entity_proc[j++] % 10;
              pro[i] = entity_proc[j++];
            }
          }

          int ierr = ex_put_elem_cmap(get_file_pointer(), Ioex::get_id(cs, (ex_entity_type)0, &ids_),
                                      TOPTR(entities), TOPTR(sides), TOPTR(procs), myProcessor);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

          // Construct the element map (internal vs. border).
          // Border elements are those in the communication map (use entities array)
          // Internal elements are the rest.  Allocate array to hold all elements,
          // initialize all to '1', then zero out the elements in 'entities'.
          // Iterate through array again and consolidate all '1's
          std::vector<char> internal(elementCount * int_byte_size_api());
          if (int_byte_size_api() == 4) {
            compute_internal_border_maps((int*)&entities[0], (int*)&internal[0], elementCount, entity_count);
          } else {
            compute_internal_border_maps((int64_t*)&entities[0], (int64_t*)&internal[0], elementCount, entity_count);
          }

          ierr = ex_put_processor_elem_maps(get_file_pointer(), TOPTR(internal),
                                            TOPTR(entities), myProcessor);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        } else {
          std::ostringstream errmsg;
          errmsg << "ERROR: Invalid commset type " << type;
          IOSS_ERROR(errmsg);
        }
      } else if (field.get_name() == "ids") {
        // Do nothing, just handles an idiosyncracy of the GroupingEntity
      } else {
        num_to_get = Ioss::Utils::field_warning(cs, field, "output");
      }
      return num_to_get;
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

    int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock* fb,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      Ioss::SerializeIO serializeIO__(this);
      size_t num_to_get = field.verify(data_size);
      if (num_to_get > 0) {

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
            int ierr = ex_put_partial_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id,
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
            size_t df_count  = fb->get_property("distribution_factor_count").get_int();
            ierr = ex_put_partial_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id,
                                                df_offset+1, df_count, static_cast<double*>(data));
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

            if (field.get_type() == Ioss::Field::INTEGER) {
              Ioss::IntVector element(num_to_get);
              Ioss::IntVector side(num_to_get);
              int *el_side = (int *)data;

              for (size_t i=0; i < num_to_get; i++) {
                element[i] = elemMap.global_to_local(el_side[index++]);
                side[i]    = el_side[index++]+side_offset;
              }

              int ierr = ex_put_partial_set(get_file_pointer(), EX_SIDE_SET, id,
                                            offset+1, entity_count, TOPTR(element), TOPTR(side));
              if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            } else {
              Ioss::Int64Vector element(num_to_get);
              Ioss::Int64Vector side(num_to_get);
              int64_t *el_side = (int64_t *)data;

              for (size_t i=0; i < num_to_get; i++) {
                element[i] = elemMap.global_to_local(el_side[index++]);
                side[i]    = el_side[index++]+side_offset;
              }

              int ierr = ex_put_partial_set(get_file_pointer(), EX_SIDE_SET, id,
                                            offset+1, entity_count, TOPTR(element), TOPTR(side));
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

              int ierr = ex_put_partial_set(get_file_pointer(), EX_SIDE_SET, id,
                                            offset+1, entity_count, TOPTR(element), TOPTR(side));
              if (ierr < 0) Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            } else {
              Ioss::Int64Vector element(num_to_get);
              Ioss::Int64Vector side(num_to_get);
              int64_t *el_side = (int64_t *)data;

              for (size_t i=0; i < num_to_get; i++) {
                element[i] = el_side[index++];
                side[i]    = el_side[index++]+side_offset;
              }

              int ierr = ex_put_partial_set(get_file_pointer(), EX_SIDE_SET, id,
                                            offset+1, entity_count, TOPTR(element), TOPTR(side));
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

      char the_title[max_line_length+1];

      // Title...
      if (region->property_exists("title")) {
        std::string title_str = region->get_property("title").get_string();
        std::strncpy(the_title, title_str.c_str(), max_line_length);
      } else {
        std::strncpy(the_title, "Sierra Output Default Title", max_line_length);
      }
      the_title[max_line_length] = '\0';


      Iofx::Mesh mesh(spatialDimension, the_title);

      Ioex::get_id(node_blocks[0], EX_NODE_BLOCK, &ids_);
      Iofx::NodeBlock N(*node_blocks[0]);
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
          Iofx::EdgeBlock T(*(*I));
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
          Iofx::FaceBlock T(*(*I));
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
          Iofx::ElemBlock T(*(*I));
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
          const Iofx::NodeSet T(*(*I));
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
          const Iofx::EdgeSet T(*(*I));
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
          const Iofx::FaceSet T(*(*I));
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
          const Iofx::ElemSet T(*(*I));
          if (std::find(mesh.elemsets.begin(), mesh.elemsets.end(), T) == mesh.elemsets.end()) {
            mesh.elemsets.push_back(T);
          }
        }
        m_groupCount[EX_ELEM_SET] = mesh.elemsets.size();
      }

      // SideSets ...
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
        Iofx::SideSet T(*(*I));
        if (std::find(mesh.sidesets.begin(), mesh.sidesets.end(), T) == mesh.sidesets.end()) {
          mesh.sidesets.push_back(T);
        }
      }
      m_groupCount[EX_SIDE_SET] = mesh.sidesets.size();

      gather_communication_metadata(&mesh.comm);

      // Write the metadata to the exodusII file...
      Iofx::Internals data(get_file_pointer(), maximumNameLength, util());
      int ierr = data.write_meta_data(mesh);
      
      if (ierr < 0)
	Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
    
      output_other_meta_data();
    }

    void DatabaseIO::gather_communication_metadata(CommunicationMetaData *meta)
    {
      // It's possible that we are a serial program outputting information
      // for later use by a parallel program.

      meta->processorCount = 0;
      meta->processorId = 0;
      meta->outputNemesis = false;

      if (isParallel) {
        meta->processorCount = util().parallel_size();
        meta->processorId = myProcessor;
        meta->outputNemesis = true;
      } else {
        if (get_region()->property_exists("processor_count")) {
          meta->processorCount = get_region()->get_property("processor_count").get_int();
        }
        if (get_region()->property_exists("my_processor")) {
          meta->processorId = get_region()->get_property("my_processor").get_int();
        }
        if (!get_region()->get_commsets().empty()) {
          isSerialParallel = true;
          meta->outputNemesis = true;
        }
      }

      if (isSerialParallel || meta->processorCount > 0) {
        meta->globalNodes    = 1; // Just need a nonzero value.
        meta->globalElements = 1; // Just need a nonzero value.

        if (get_region()->property_exists("global_node_count")) {
          meta->globalNodes = get_region()->get_property("global_node_count").get_int();
        }

        if (get_region()->property_exists("global_element_count")) {
          meta->globalElements = get_region()->get_property("global_element_count").get_int();
        }

        if (get_region()->property_exists("global_element_block_count")) {
          meta->globalElementBlocks = get_region()->get_property("global_element_block_count").get_int();
        } else {
          meta->globalElementBlocks = get_region()->get_element_blocks().size();
        }

        if (get_region()->property_exists("global_node_set_count")) {
          meta->globalNodeSets = get_region()->get_property("global_node_set_count").get_int();
        } else {
          meta->globalNodeSets = get_region()->get_nodesets().size();
        }

        if (get_region()->property_exists("global_side_set_count")) {
          meta->globalSideSets = get_region()->get_property("global_side_set_count").get_int();
        } else {
          meta->globalSideSets = get_region()->get_sidesets().size();
        }

        // ========================================================================
        // Load balance parameters (NEMESIS, p15)
        meta->nodesInternal = nodeCount;
        meta->nodesBorder   = 0;
        meta->nodesExternal = 0; // Shadow nodes == 0 for now
        meta->elementsInternal = elementCount;
        meta->elementsBorder   = 0;

        // Now, see if any of the above are redefined by a property...
        if (get_region()->property_exists("internal_node_count")) {
          meta->nodesInternal = get_region()->
            get_property("internal_node_count").get_int();
        }

        if (get_region()->property_exists("border_node_count")) {
          meta->nodesBorder = get_region()->
            get_property("border_node_count").get_int();
        }

        if (get_region()->property_exists("internal_element_count")) {
          meta->elementsInternal = get_region()->
            get_property("internal_element_count").get_int();
        }

        if (get_region()->property_exists("border_element_count")) {
          meta->elementsBorder = get_region()->
            get_property("border_element_count").get_int();
        }

        Ioss::CommSetContainer comm_sets = get_region()->get_commsets();
        Ioss::CommSetContainer::const_iterator I = comm_sets.begin();
        while (I != comm_sets.end()) {

          Ioss::CommSet *cs = *I;

          std::string type = cs->get_property("entity_type").get_string();
          size_t count = cs->get_property("entity_count").get_int();
          int64_t id = Ioex::get_id(cs, (ex_entity_type)0, &ids_);

          if (type == "node") {
            meta->nodeMap.push_back(Iofx::CommunicationMap(id, count, 'n'));
          } else if (type == "side") {
            meta->elementMap.push_back(Iofx::CommunicationMap(id, count, 'e'));
          } else {
            std::ostringstream errmsg;
            errmsg << "Internal Program Error...";
            IOSS_ERROR(errmsg);
          }
          ++I;
        }
      }
      commsetNodeCount = meta->nodeMap.size();
      commsetElemCount = meta->elementMap.size();
    }
  } // End of namespace

