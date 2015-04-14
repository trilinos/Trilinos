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
#include <exodus/Ioex_DatabaseIO.h>
#include <exodus/Ioex_Utils.h>
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

  void check_variable_consistency(const ex_var_params &exo_params,
                                  int my_processor, const std::string &filename,
                                  const Ioss::ParallelUtils &util);

  void check_attribute_index_order(Ioss::GroupingEntity *block);

  template <typename T>
  void write_attribute_names(int exoid, ex_entity_type type, const std::vector<T*>& entities,
                             const char suffix_separator);

  template <typename T>
  void generate_block_truth_table(Ioex::VariableNameMap &variables,
                                  Ioss::IntVector &truth_table,
                                  std::vector<T*> &blocks,
                                  char field_suffix_separator);

}

namespace Ioex {
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
                         Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                         const Ioss::PropertyManager &props) :
    Ioss::DatabaseIO(region, filename, db_usage, communicator, props),
    exodusFilePtr(-1), databaseTitle(""), exodusMode(EX_CLOBBER), dbRealWordSize(8),
    maximumNameLength(32), spatialDimension(0),
    nodeCount(0), edgeCount(0), faceCount(0), elementCount(0),
    commsetNodeCount(0), commsetElemCount(0),
    nodeMap("node"), edgeMap("edge"), faceMap("face"), elemMap("element"),
    timeLastFlush(0), fileExists(false),
    minimizeOpenFiles(false), blockAdjacenciesCalculated(false),
    nodeConnectivityStatusCalculated(false)
  {
    m_groupCount[EX_GLOBAL]     = 1; // To make some common code work more cleanly.
    m_groupCount[EX_NODE_BLOCK] = 1; // To make some common code work more cleanly.

    // A history file is only written on processor 0...
    if (db_usage == Ioss::WRITE_HISTORY)
      isParallel = false;

    dbState = Ioss::STATE_UNKNOWN;

    // Set exodusII warning level.
    if (util().get_environment("EX_DEBUG", isParallel)) {
      std::cerr << "IOEX: Setting EX_VERBOSE|EX_DEBUG because EX_DEBUG environment variable is set.\n";
      ex_opts(EX_VERBOSE|EX_DEBUG);
    }

    if (!is_input()) {
      // Check whether appending to existing file...
      if (open_create_behavior() == Ioss::DB_APPEND || open_create_behavior() == Ioss::DB_APPEND_GROUP) {
        // Append to file if it already exists -- See if the file exists.
        std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
        Ioss::FileInfo file = Ioss::FileInfo(decoded_filename);
        fileExists = file.exists();
      }

      if (util().get_environment("EX_MINIMIZE_OPEN_FILES", isParallel)) {
        std::cerr << "IOEX: Minimizing open files because EX_MINIMIZE_OPEN_FILES environment variable is set.\n";
        minimizeOpenFiles = true;
      }

      if (util().get_environment("EX_MODE", exodusMode, isParallel)) {
        std::cerr << "IOEX: Exodus create mode set to " << exodusMode
                  << " from value of EX_MODE environment variable.\n";
      }
    }

    // See if there are any properties that need to (or can) be
    // handled prior to opening/creating database...
    bool compress = ((properties.exists("COMPRESSION_LEVEL") &&
                      properties.get("COMPRESSION_LEVEL").get_int() > 0) ||
                     (properties.exists("COMPRESSION_SHUFFLE") &&
                      properties.get("COMPRESSION_SHUFFLE").get_int() > 0));

    if (compress) {
      exodusMode |= EX_NETCDF4;
    }

    if (properties.exists("FILE_TYPE")) {
      std::string type = properties.get("FILE_TYPE").get_string();
      if (type == "netcdf4" || type == "netcdf-4" || type == "hdf5") {
        exodusMode |= EX_NETCDF4;
      }
    }

    if (properties.exists("ENABLE_FILE_GROUPS")) {
        exodusMode |= EX_NETCDF4;
        exodusMode |= EX_NOCLASSIC;
    }

    if (properties.exists("MAXIMUM_NAME_LENGTH")) {
      maximumNameLength = properties.get("MAXIMUM_NAME_LENGTH").get_int();
    }

    if (properties.exists("REAL_SIZE_DB")) {
      int rsize = properties.get("REAL_SIZE_DB").get_int();
      if (rsize == 4) {
        dbRealWordSize = 4; // Only used for file create...
      }
    }

    if (properties.exists("INTEGER_SIZE_DB")) {
      int isize = properties.get("INTEGER_SIZE_DB").get_int();
      if (isize == 8) {
        exodusMode |= EX_ALL_INT64_DB;
      }
    }

    if (properties.exists("INTEGER_SIZE_API")) {
      int isize = properties.get("INTEGER_SIZE_API").get_int();
      if (isize == 8) {
        exodusMode |= EX_ALL_INT64_API;
        set_int_byte_size_api(Ioss::USE_INT64_API);
      }
    }

    // Don't open output files until they are actually going to be
    // written to.  This is needed for proper support of the topology
    // files and auto restart so we don't overwrite a file with data we
    // need to save...
  }

  // common
  DatabaseIO::~DatabaseIO()
  {
    try {
      free_file_pointer();
    } catch (...) {
    }
  }

  // common
  void DatabaseIO::release_memory()
  {
    nodeMap.release_memory();
    edgeMap.release_memory();
    faceMap.release_memory();
    elemMap.release_memory();
  }

  // common
  unsigned DatabaseIO::entity_field_support() const
  {
    return Ioss::NODEBLOCK | Ioss::EDGEBLOCK | Ioss::FACEBLOCK | Ioss::ELEMENTBLOCK |
      Ioss::NODESET   | Ioss::EDGESET   | Ioss::FACESET   | Ioss::ELEMENTSET   |
      Ioss::SIDESET   | Ioss::SIDEBLOCK | Ioss::REGION    | Ioss::SUPERELEMENT;
  }

  // common
  int DatabaseIO::free_file_pointer() const
  {
    if (exodusFilePtr != -1)
      ex_close(exodusFilePtr);
    exodusFilePtr = -1;

    return exodusFilePtr;
  }

  bool DatabaseIO::open_group(const std::string &group_name)
  {
    // Get existing file pointer...
    bool success = false;

    int exoid = get_file_pointer();

    m_groupName = group_name;
    ex_get_group_id(exoid, m_groupName.c_str(), &exodusFilePtr);

    if (exodusFilePtr < 0) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Could not open group named '" << m_groupName
            << "' in file '" << get_filename() << "'.\n";
      IOSS_ERROR(errmsg);
    } else {
      success = true;
    }
    return success;
  }

  bool DatabaseIO::create_subgroup(const std::string &group_name)
  {
    bool success = false;
    if (!is_input()) {
      // Get existing file pointer...
      int exoid = get_file_pointer();

      // Check name for '/' which is not allowed since it is the
      // separator character in a full group path
      if (group_name.find('/') != std::string::npos) {
       std::ostringstream errmsg;
       errmsg << "ERROR: Invalid group name '" << m_groupName
              << "' contains a '/' which is not allowed.\n";
       IOSS_ERROR(errmsg);
      }

      m_groupName = group_name;
      exoid = ex_create_group(exoid, m_groupName.c_str());
      if (exoid < 0) {
       std::ostringstream errmsg;
       errmsg << "ERROR: Could not create group named '" << m_groupName
              << "' in file '" << get_filename() << "'.\n";
       IOSS_ERROR(errmsg);
      }
      else {
       exodusFilePtr = exoid;
       success = true;
      }
    }
    return success;
  }

  // common
  void DatabaseIO::put_qa()
  {
    struct qa_element {
      char *qa_record[1][4];
    };
    
    size_t num_qa_records = qaRecords.size()/4;

    qa_element *qa = new qa_element[num_qa_records+1];
    for (size_t i=0; i < num_qa_records+1; i++) {
      for (int j=0; j < 4; j++) {
        qa[i].qa_record[0][j] = new char[MAX_STR_LENGTH+1];
      }
    }

    {
      int j = 0;
      for (size_t i=0; i < num_qa_records; i++) {
        std::strncpy(qa[i].qa_record[0][0], qaRecords[j++].c_str(), MAX_STR_LENGTH);
        std::strncpy(qa[i].qa_record[0][1], qaRecords[j++].c_str(), MAX_STR_LENGTH);
        std::strncpy(qa[i].qa_record[0][2], qaRecords[j++].c_str(), MAX_STR_LENGTH);
        std::strncpy(qa[i].qa_record[0][3], qaRecords[j++].c_str(), MAX_STR_LENGTH);
      }
    }

    Ioss::Utils::time_and_date(qa[num_qa_records].qa_record[0][3],
                               qa[num_qa_records].qa_record[0][2], MAX_STR_LENGTH);

    std::string codename = "unknown";
    std::string version  = "unknown";

    if (get_region()->property_exists("code_name")) {
      codename = get_region()->get_property("code_name").get_string();
    }
    if (get_region()->property_exists("code_version")) {
      version = get_region()->get_property("code_version").get_string();
    }

    char buffer[MAX_STR_LENGTH+1];
    std::strncpy(buffer, codename.c_str(), MAX_STR_LENGTH);
    buffer[MAX_STR_LENGTH] = '\0';
    std::strcpy(qa[num_qa_records].qa_record[0][0], buffer);

    std::strncpy(buffer, version.c_str(), MAX_STR_LENGTH);
    buffer[MAX_STR_LENGTH] = '\0';
    std::strcpy(qa[num_qa_records].qa_record[0][1], buffer);
    
    int ierr = ex_put_qa(get_file_pointer(), num_qa_records+1,
			 myProcessor == 0 ? qa[0].qa_record : NULL);
    if (ierr < 0)
      Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

    for (size_t i=0; i < num_qa_records+1; i++) {
      for (int j=0; j < 4; j++) {
        delete [] qa[i].qa_record[0][j];
      }
    }
    delete [] qa;
  }

  // common
  void DatabaseIO::put_info()
  {
    // dump info records, include the product_registry
    // See if the input file was specified as a property on the database...
    std::string filename;
    std::vector<std::string> input_lines;
    if (get_region()->property_exists("input_file_name")) {
      filename = get_region()->get_property("input_file_name").get_string();
      // Determine size of input file so can embed it in info records...
      Ioss::Utils::input_file(filename, &input_lines, max_line_length);
    }

    // See if the client added any "information_records" 
    size_t info_rec_size = informationRecords.size();
    size_t in_lines = input_lines.size();
    size_t qa_lines = 2; // Platform info and Version info...

    size_t total_lines = in_lines + qa_lines + info_rec_size;

    char** info = get_exodus_names(total_lines, max_line_length); // 'total_lines' pointers to char buffers

    int i = 0;
    std::strncpy(info[i++], Ioss::Utils::platform_information().c_str(),
                 max_line_length);

    std::strncpy(info[i++], Ioex::Version(), max_line_length);

    // Copy input file lines into 'info' array...
    for (size_t j=0; j < input_lines.size(); j++, i++) {
      std::strncpy(info[i], input_lines[j].c_str(), max_line_length);
      info[i][max_line_length] = '\0'; // Once more for good luck...
    }

    // Copy "information_records" property data ...
    for (size_t j=0; j < informationRecords.size(); j++, i++) {
      std::strncpy(info[i], informationRecords[j].c_str(), max_line_length);
      info[i][max_line_length] = '\0'; // Once more for good luck...
    }

    int ierr = ex_put_info(get_file_pointer(), total_lines, myProcessor == 0 ? info : NULL);
    if (ierr < 0)
      Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

    delete_exodus_names(info, total_lines);
  }

  // common
  int DatabaseIO::get_current_state() const
  {
    int step = get_region()->get_property("current_state").get_int();

    if (step <= 0) {
      std::ostringstream errmsg;
      errmsg << "ERROR: No currently active state.  The calling code must call Ioss::Region::begin_state(int step)\n"
             << "       to set the database timestep from which to read the transient data.\n"
             << "       [" << get_filename() << "]\n";
      IOSS_ERROR(errmsg);
    }
    return step;
  }

  // common
  void DatabaseIO::get_nodeblocks() 
  {
    // For exodusII, there is only a single node block which contains
    // all of the nodes.
    // The default id assigned is '1' and the name is 'nodeblock_1'

    std::string block_name = "nodeblock_1";
    Ioss::NodeBlock *block = new Ioss::NodeBlock(this, block_name,
                                                 nodeCount, spatialDimension);
    block->property_add(Ioss::Property("id", 1));
    // Check for results variables.

    int num_attr = 0;
    {
      Ioss::SerializeIO serializeIO__(this);
      int ierr = ex_get_attr_param(get_file_pointer(), EX_NODE_BLOCK, 1, &num_attr);
      if (ierr < 0)
        Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
    }

    add_attribute_fields(EX_NODE_BLOCK, block, num_attr, "");
    add_results_fields(EX_NODE_BLOCK, block);

    bool added = get_region()->add(block);
    if(!added) {
      delete block;
    }
  }

  // common
  void DatabaseIO::get_block_adjacencies(const Ioss::ElementBlock *eb,
                                         std::vector<std::string> &block_adjacency) const
  {
    if (!blockAdjacenciesCalculated) {
      compute_block_adjacencies();
    }

    // Extract the computed block adjacency information for this
    // element block:
    // Debug print...
    int blk_position  = eb->get_property("original_block_order").get_int();

    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    assert(check_block_order(element_blocks));
    for (int iblk = 0; iblk < m_groupCount[EX_ELEM_BLOCK]; iblk++) {

      Ioss::ElementBlock *leb = element_blocks[iblk];
      int lblk_position  = leb->get_property("original_block_order").get_int();

      if (blk_position != lblk_position &&
          blockAdjacency[blk_position][lblk_position] == 1) {
        block_adjacency.push_back(leb->name());
      }
    }
  }

  // common
  void DatabaseIO::compute_block_membership(Ioss::SideBlock *efblock,
                                              std::vector<std::string> &block_membership) const
    {
      Ioss::Int64Vector block_ids(m_groupCount[EX_ELEM_BLOCK]);
      if (m_groupCount[EX_ELEM_BLOCK] == 1) {
        block_ids[0] = 1;
      } else {
        Ioss::Int64Vector element_side;
        if (int_byte_size_api() == 4) {
          Ioss::IntVector es32;
          efblock->get_field_data("element_side", es32);
          element_side.resize(es32.size());
          std::copy(es32.begin(), es32.end(), element_side.begin());
        } else {
          efblock->get_field_data("element_side", element_side);
        }

        size_t number_sides = element_side.size() / 2;
        Ioss::ElementBlock *block = NULL;
        for (size_t iel = 0; iel < number_sides; iel++) {
          int64_t elem_id = element_side[2*iel];  // Vector contains both element and side.
          elem_id = elemMap.global_to_local(elem_id);
          if (block == NULL || !block->contains(elem_id)) {
            block = get_region()->get_element_block(elem_id);
            assert(block != NULL);
            int block_order = block->get_property("original_block_order").get_int();
            block_ids[block_order] = 1;
          }
        }
      }

      // Synchronize among all processors....
      if (isParallel) {
        util().global_array_minmax(block_ids, Ioss::ParallelUtils::DO_MAX);
      }

      Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
      assert(check_block_order(element_blocks));

      for (int i=0; i < m_groupCount[EX_ELEM_BLOCK]; i++) {
        if (block_ids[i] == 1) {
          Ioss::ElementBlock *block = element_blocks[i];
          if (!Ioss::Utils::block_is_omitted(block)) {
            block_membership.push_back(block->name());
          }
        }
      }
    }

    // common
    int64_t DatabaseIO::get_field_internal(const Ioss::Region* /* region */,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      // For now, assume that all TRANSIENT fields on a region
      // are REDUCTION fields (1 value).  We need to gather these
      // and output them all at one time.  The storage location is a
      // 'globalVariables' array
      {
        size_t num_to_get = field.verify(data_size);
        Ioss::SerializeIO       serializeIO__(this);

        Ioss::Field::RoleType role = field.get_role();

        if (role == Ioss::Field::TRANSIENT || role == Ioss::Field::REDUCTION) {
          get_reduction_field(EX_GLOBAL, field, get_region(), data);
        } else {
          std::ostringstream errmsg;
          errmsg << "ERROR: Can not handle non-TRANSIENT or non-REDUCTION fields on regions";
          IOSS_ERROR(errmsg);
        }
        return num_to_get;
      }
    }

  // common
    int64_t DatabaseIO::put_field_internal(const Ioss::Region* /* region */,
                                           const Ioss::Field& field,
                                           void *data, size_t data_size) const
    {
      // For now, assume that all TRANSIENT fields on a region
      // are REDUCTION fields (1 value).  We need to gather these
      // and output them all at one time.  The storage location is a
      // 'globalVariables' array
      {
        Ioss::SerializeIO       serializeIO__(this);

        Ioss::Field::RoleType role = field.get_role();
        size_t num_to_get = field.verify(data_size);

        if ((role == Ioss::Field::TRANSIENT || role == Ioss::Field::REDUCTION) &&
            num_to_get == 1) {
          store_reduction_field(EX_GLOBAL, field, get_region(), data);
        }
        else if (num_to_get != 1) {
          // There should have been a warning/error message printed to the
          // log file earlier for this, so we won't print anything else
          // here since it would be printed for each and every timestep....
          ;
        } else {
          std::ostringstream errmsg;
          errmsg << "ERROR: The variable named '" << field.get_name()
                 << "' is of the wrong type. A region variable must be of type"
                 << " TRANSIENT or REDUCTION.\n"
                 << "This is probably an internal error; please notify gdsjaar@sandia.gov";
          IOSS_ERROR(errmsg);
        }
        return num_to_get;
      }
    }

    namespace {
      // common
      template <typename T>
      void generate_block_truth_table(VariableNameMap &variables,
                                      Ioss::IntVector &truth_table,
                                      std::vector<T*> &blocks,
                                      char field_suffix_separator)
      {
        size_t block_count = blocks.size();
        size_t var_count = variables.size();

        if (var_count == 0 || block_count == 0)
          return;

        truth_table.resize(block_count*var_count);

        // Fill in the truth table.  It is conceptually a two-dimensional array of
        // the form 'array[num_blocks][num_element_var]'.  In C++,
        // the values for the first block are first, followed by
        // next block, ...
        size_t offset = 0;
        typename std::vector<T*>::const_iterator I;
        for (I=blocks.begin(); I != blocks.end(); ++I) {
          // Get names of all transient and reduction fields...
          Ioss::NameList results_fields;
          (*I)->field_describe(Ioss::Field::TRANSIENT, &results_fields);
          (*I)->field_describe(Ioss::Field::REDUCTION, &results_fields);

          Ioss::NameList::const_iterator IF;
          for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
            std::string field_name = *IF;

            Ioss::Field field = (*I)->get_field(field_name);
            const Ioss::VariableType *var_type = field.transformed_storage();
            Ioss::Field::BasicType ioss_type = field.get_type();

            int re_im = 1;
            if (ioss_type == Ioss::Field::COMPLEX)
              re_im = 2;
            for (int complex_comp = 0; complex_comp < re_im; complex_comp++) {
              field_name = field.get_name();
              if (re_im == 2) {
                field_name += complex_suffix[complex_comp];
              }

              for (int i=1; i <= var_type->component_count(); i++) {
                std::string var_string = var_type->label_name(field_name, i, field_suffix_separator);
                // Find position of 'var_string' in 'variables'
                VariableNameMap::iterator VN = variables.find(var_string);
                if (VN != variables.end()) {
                  // Index '(*VN).second' is 1-based...
                  truth_table[offset + (*VN).second-1] = 1;
                }
              }
            }
          }
          offset += var_count;
        }
        assert(offset == var_count * block_count);
      }
    }
  // common
    void DatabaseIO::store_reduction_field(ex_entity_type type,
                                           const Ioss::Field& field,
                                           const Ioss::GroupingEntity *ge,
                                           void *variables) const
    {
      const Ioss::VariableType *var_type = field.transformed_storage();

      Ioss::Field::BasicType ioss_type = field.get_type();
      assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER ||
             ioss_type == Ioss::Field::INT64 || ioss_type == Ioss::Field::COMPLEX);
      double  *rvar = static_cast<double*>(variables);
      int     *ivar = static_cast<int*>(variables);
      int64_t *ivar64 = static_cast<int64_t*>(variables);

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
      // of this name in 'm_variables[EX_GLOBAL]' map
      int comp_count = var_type->component_count();
      int var_index=0;

      int re_im = 1;
      if (field.get_type() == Ioss::Field::COMPLEX)
        re_im = 2;
      for (int complex_comp = 0; complex_comp < re_im; complex_comp++) {
        std::string field_name = field.get_name();
        if (re_im == 2) {
          field_name += complex_suffix[complex_comp];
        }

        char field_suffix_separator = get_field_separator();
        for (int i=0; i < comp_count; i++) {
          std::string var_name = var_type->label_name(field_name, i+1, field_suffix_separator);

          // If this is not a global variable, prepend the name to avoid
          // name collisions...
          if (type != EX_GLOBAL) {
            var_name = ge->name() + ":" + var_name;
          }
          assert(m_variables[EX_GLOBAL].find(var_name) != m_variables[EX_GLOBAL].end());
          var_index = m_variables[EX_GLOBAL].find(var_name)->second;

          assert(static_cast<int>(globalValues.size()) >= var_index);

          // Transfer from 'variables' array.
          if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX)
            globalValues[var_index-1] = rvar[i];
          else if (ioss_type == Ioss::Field::INTEGER)
            globalValues[var_index-1] = ivar[i];
          else if (ioss_type == Ioss::Field::INT64)
            globalValues[var_index-1] = ivar64[i]; // FIX 64 UNSAFE
        }
      }
    }

      // common
    void DatabaseIO::get_reduction_field(ex_entity_type,
                                         const Ioss::Field& field,
                                         const Ioss::GroupingEntity */* ge */,
                                         void *variables) const
    {
      const Ioss::VariableType *var_type = field.raw_storage();

      Ioss::Field::BasicType ioss_type = field.get_type();
      assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER || ioss_type == Ioss::Field::INT64);
      double *rvar = static_cast<double*>(variables);
      int  *ivar = static_cast<int*>(variables);
      int64_t *i64var = static_cast<int64_t*>(variables);

      // get number of components, cycle through each component
      // and add suffix to base 'field_name'.  Look up index
      // of this name in 'm_variables[EX_GLOBAL]' map
      char field_suffix_separator = get_field_separator();

      int comp_count = var_type->component_count();
      for (int i=0; i < comp_count; i++) {
        std::string field_name = field.get_name();
        std::string var_name = var_type->label_name(field_name, i+1, field_suffix_separator);

        assert(m_variables[EX_GLOBAL].find(var_name) != m_variables[EX_GLOBAL].end());
        int var_index = m_variables[EX_GLOBAL].find(var_name)->second;

        assert(static_cast<int>(globalValues.size()) >= var_index);

        // Transfer to 'variables' array.
        if (ioss_type == Ioss::Field::REAL)
          rvar[i] = globalValues[var_index-1];
        else if (ioss_type == Ioss::Field::INT64)
          i64var[i] = static_cast<int64_t>(globalValues[var_index-1]);
        else if (ioss_type == Ioss::Field::INTEGER)
          ivar[i] = static_cast<int>(globalValues[var_index-1]);
      }
    }

      // common
    void DatabaseIO::write_reduction_fields() const
    {
      int step = get_current_state();
      step = get_database_step(step);
      size_t count = globalValues.size();
      if (count > 0) {
        int ierr = ex_put_var(get_file_pointer(), step, EX_GLOBAL, 1, 0,
                              count, (double*)TOPTR(globalValues));
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
    }

      // common
    void DatabaseIO::read_reduction_fields() const
    {
      int step = get_current_state();
      size_t count = globalValues.size();
      if (count > 0) {
        int ierr = ex_get_var(get_file_pointer(), step, EX_GLOBAL, 1, 0,
                              count, TOPTR(globalValues));
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
    }

      // common
    bool DatabaseIO::begin(Ioss::State state)
    {
      dbState = state;
      return true;
    }

      // common
    bool DatabaseIO::end(Ioss::State state)
    {
      // Transitioning out of state 'state'
      assert(state == dbState);
      switch (state) {
      case Ioss::STATE_DEFINE_MODEL:
        if (!is_input() && open_create_behavior() != Ioss::DB_APPEND) {
          write_meta_data();
        }
        break;
      case Ioss::STATE_DEFINE_TRANSIENT:
        if (!is_input() && open_create_behavior() != Ioss::DB_APPEND) {
          write_results_metadata();
        }
        break;
      default: // ignore everything else...
        break;
      }

      {
        Ioss::SerializeIO       serializeIO__(this);

        if (!is_input()) {
          ex_update(get_file_pointer());
          if (minimizeOpenFiles)
            free_file_pointer();
        }
        dbState = Ioss::STATE_UNKNOWN;
      }

      return true;
    }

    // Default versions do nothing at this time...
    // Will be used for global variables...
      // common
    bool DatabaseIO::begin_state(Ioss::Region */* region */, int state, double time)
    {
      Ioss::SerializeIO serializeIO__(this);

      time /= timeScaleFactor;
      
      state = get_database_step(state);
      if (!is_input()) {
        int ierr = ex_put_time(get_file_pointer(), state, &time);
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        // Zero global variable array...
        std::fill(globalValues.begin(), globalValues.end(), 0.0);

      } else {
        // Store reduction variables
        read_reduction_fields();
      }
      return true;
    }
      
      // common
    bool DatabaseIO::end_state(Ioss::Region*, int, double time)
    {
      Ioss::SerializeIO serializeIO__(this);

      if (!is_input()) {
        write_reduction_fields();
        time /= timeScaleFactor;
        finalize_write(time);
        if (minimizeOpenFiles)
          free_file_pointer();
      }
      return true;
    }

      // common
    void DatabaseIO::add_region_fields()
    {
      int field_count = add_results_fields(EX_GLOBAL, get_region());
      globalValues.resize(field_count);
    }

      // common
    int64_t DatabaseIO::add_results_fields(ex_entity_type type, Ioss::GroupingEntity *entity, int64_t position)
    {
      return internal_add_results_fields(type, entity, position, m_groupCount[type], m_truthTable[type], m_variables[type]);
    }

      // common
    int64_t DatabaseIO::internal_add_results_fields(ex_entity_type type, Ioss::GroupingEntity *entity, int64_t position,
                                                    int64_t block_count, Ioss::IntVector &truth_table, Ioex::VariableNameMap &variables)
    {
      int nvar = 0;
      {
        Ioss::SerializeIO       serializeIO__(this);

        int ierr = ex_get_variable_param(get_file_pointer(), type, &nvar);
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }

      if (nvar > 0) {
        if (truth_table.empty()) {
          truth_table.resize(block_count * nvar);
        }

        // Read and store the truth table (Should be there since we only
        // get to this routine if there are variables...)
        if (!truth_table.empty()) {
          {
            Ioss::SerializeIO   serializeIO__(this);

            if (type == EX_NODE_BLOCK || type == EX_GLOBAL) {
              // These types don't have a truth table in the exodus api...
              // They do in Ioss just for some consistency...
              std::fill(truth_table.begin(), truth_table.end(), 1);
            }
            else {
              int ierr = ex_get_truth_table(get_file_pointer(), type, block_count, nvar, &truth_table[0]);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }
          }

          // If parallel, then synchronize the truth table among all processors...
          // Need to know that block_X has variable_Y even if block_X is
          // empty on a specific processor...  The truth table contains 0
          // if the variable doesn't exist and 1 if it does, so we just
          // take the maximum at each location...
          // This is a collective call...
          if (isParallel) {
            util().global_array_minmax(truth_table, Ioss::ParallelUtils::DO_MAX);
          }
        }

        // Get the variable names and add as fields. Need to decode these
        // into vector/tensor/... eventually, for now store all as
        // scalars.
        char **names = get_exodus_names(nvar, maximumNameLength);

        // Read the names...
        // (Currently, names are read for every block.  We could save them...)
        {
          Ioss::SerializeIO     serializeIO__(this);

          int ierr = ex_get_variable_names(get_file_pointer(), type, nvar, names);
          if (ierr < 0)
            Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

          // Add to VariableNameMap so can determine exodusII index given a
          // Sierra field name.  exodusII index is just 'i+1'
          for (int i=0; i < nvar; i++) {
            if (lowerCaseVariableNames)
              Ioss::Utils::fixup_name(names[i]);
            variables.insert(VNMValuePair(std::string(names[i]),   i+1));
          }

          int offset = position*nvar;
          int *local_truth = NULL;
          if (!truth_table.empty())
            local_truth = &truth_table[offset];

          std::vector<Ioss::Field> fields;
          int64_t count = entity->get_property("entity_count").get_int();
	  Ioex::get_fields(count, names, nvar, Ioss::Field::TRANSIENT,
			   get_field_separator(), local_truth, fields);

          std::vector<Ioss::Field>::const_iterator IF;
          for (IF = fields.begin(); IF != fields.end(); ++IF) {
            entity->field_add(*IF);
          }

          for (int i=0; i < nvar; i++) {
            // Verify that all names were used for a field...
            assert(names[i][0] == '\0' || local_truth[i] == 0);
            delete [] names[i];
          }
          delete [] names;
        }
      }
      return nvar;
    }

      // common
    void DatabaseIO::write_results_metadata()
    {
      int glob_index = 0;
      glob_index = gather_names(EX_GLOBAL, m_variables[EX_GLOBAL], get_region(),
                                glob_index, true);
      assert(glob_index == static_cast<int>(m_variables[EX_GLOBAL].size()));

      Ioss::NodeBlockContainer node_blocks = get_region()->get_node_blocks();
      assert(node_blocks.size() == 1);
      internal_write_results_metadata(EX_NODE_BLOCK, node_blocks, glob_index);

      Ioss::EdgeBlockContainer edge_blocks = get_region()->get_edge_blocks();
      internal_write_results_metadata(EX_EDGE_BLOCK, edge_blocks, glob_index);

      Ioss::FaceBlockContainer face_blocks = get_region()->get_face_blocks();
      internal_write_results_metadata(EX_FACE_BLOCK, face_blocks, glob_index);

      Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
      internal_write_results_metadata(EX_ELEM_BLOCK, element_blocks, glob_index);

      Ioss::NodeSetContainer nodesets = get_region()->get_nodesets();
      internal_write_results_metadata(EX_NODE_SET, nodesets, glob_index);

      Ioss::EdgeSetContainer edgesets = get_region()->get_edgesets();
      internal_write_results_metadata(EX_EDGE_SET, edgesets, glob_index);

      Ioss::FaceSetContainer facesets = get_region()->get_facesets();
      internal_write_results_metadata(EX_FACE_SET, facesets, glob_index);

      Ioss::ElementSetContainer elementsets = get_region()->get_elementsets();
      internal_write_results_metadata(EX_ELEM_SET, elementsets, glob_index);

      {
        int index = 0;
        Ioss::SideSetContainer sidesets = get_region()->get_sidesets();
        Ioss::SideSetContainer::const_iterator I;
        for (I=sidesets.begin(); I != sidesets.end(); ++I) {
          Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
          Ioss::SideBlockContainer::const_iterator J;

          for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
            glob_index = gather_names(EX_SIDE_SET, m_variables[EX_SIDE_SET], *J, glob_index, true);
            index = gather_names(EX_SIDE_SET, m_variables[EX_SIDE_SET], *J, index, false);
          }
        }
        assert(index == static_cast<int>(m_variables[EX_SIDE_SET].size()));
        generate_sideset_truth_table();
      }

      ex_var_params exo_params;
      exo_params.num_glob = m_variables[EX_GLOBAL].size();
      exo_params.num_node = m_variables[EX_NODE_BLOCK].size();
      exo_params.num_edge = m_variables[EX_EDGE_BLOCK].size();
      exo_params.num_face = m_variables[EX_FACE_BLOCK].size();
      exo_params.num_elem = m_variables[EX_ELEM_BLOCK].size();
      exo_params.num_nset = m_variables[EX_NODE_SET].size();
      exo_params.num_eset = m_variables[EX_EDGE_SET].size();
      exo_params.num_fset = m_variables[EX_FACE_SET].size();
      exo_params.num_sset = m_variables[EX_SIDE_SET].size();
      exo_params.num_elset= m_variables[EX_ELEM_SET].size();

      exo_params.edge_var_tab  = TOPTR(m_truthTable[EX_EDGE_BLOCK]);
      exo_params.face_var_tab  = TOPTR(m_truthTable[EX_FACE_BLOCK]);
      exo_params.elem_var_tab  = TOPTR(m_truthTable[EX_ELEM_BLOCK]);
      exo_params.nset_var_tab  = TOPTR(m_truthTable[EX_NODE_SET]);
      exo_params.eset_var_tab  = TOPTR(m_truthTable[EX_EDGE_SET]);
      exo_params.fset_var_tab  = TOPTR(m_truthTable[EX_FACE_SET]);
      exo_params.sset_var_tab  = TOPTR(m_truthTable[EX_SIDE_SET]);
      exo_params.elset_var_tab = TOPTR(m_truthTable[EX_ELEM_SET]);

      if (isParallel) {
        // Check consistency among all processors.  They should all
        // have the same number of each variable type...
        // The called function will throw an exception if the counts differ.
        check_variable_consistency(exo_params, myProcessor, get_filename(), util());
      }

      {
        Ioss::SerializeIO       serializeIO__(this);

        int ierr = ex_put_all_var_param_ext(get_file_pointer(), &exo_params);
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

        globalValues.resize(m_variables[EX_GLOBAL].size());
        output_results_names(EX_GLOBAL, m_variables[EX_GLOBAL]);
	if (nodeCount > 0)
	  output_results_names(EX_NODE_BLOCK, m_variables[EX_NODE_BLOCK]);
        output_results_names(EX_EDGE_BLOCK, m_variables[EX_EDGE_BLOCK]);
        output_results_names(EX_FACE_BLOCK, m_variables[EX_FACE_BLOCK]);
        output_results_names(EX_ELEM_BLOCK, m_variables[EX_ELEM_BLOCK]);
        output_results_names(EX_NODE_SET,   m_variables[EX_NODE_SET]);
        output_results_names(EX_EDGE_SET,   m_variables[EX_EDGE_SET]);
        output_results_names(EX_FACE_SET,   m_variables[EX_FACE_SET]);
        output_results_names(EX_ELEM_SET,   m_variables[EX_ELEM_SET]);
        output_results_names(EX_SIDE_SET,   m_variables[EX_SIDE_SET]);
      }
    }

      // common
    template <typename T>
      void DatabaseIO::internal_write_results_metadata(ex_entity_type type,
                                                       std::vector<T*> entities,
                                                       int &glob_index)
      {
        typename std::vector<T*>::const_iterator I;

        int index = 0;
        for (I=entities.begin(); I != entities.end(); ++I) {
          glob_index = gather_names(type, m_variables[type], *I, glob_index, true);
          index = gather_names(type, m_variables[type], *I, index, false);
        }
        assert(index == static_cast<int>(m_variables[type].size()));
        generate_block_truth_table(m_variables[type], m_truthTable[type], entities,
                                   get_field_separator());
      }

      // common
    int DatabaseIO::gather_names(ex_entity_type type,
                                 VariableNameMap &variables,
                                 const Ioss::GroupingEntity *ge,
                                 int index, bool reduction)
    {
      int new_index = index;

      bool nblock = (type == EX_NODE_BLOCK);

      // Get names of all transient and reduction fields...
      Ioss::NameList results_fields;
      if (reduction)
        ge->field_describe(Ioss::Field::REDUCTION, &results_fields);

      if (!reduction || type == EX_GLOBAL)
        ge->field_describe(Ioss::Field::TRANSIENT, &results_fields);

      // NOTE: For exodusII, the convention is that the displacement
      //       fields are the first 'ndim' fields in the file.
      //       Try to find a likely displacement field
      std::string disp_name;
      bool has_disp = false;
      if (!reduction && nblock && new_index == 0) {
        has_disp = find_displacement_field(results_fields, ge, spatialDimension, &disp_name);
        if (has_disp)
          new_index += spatialDimension;
      }

      int save_index = 0;
      Ioss::NameList::const_iterator IF;

      for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
        std::string name = *IF;

        if (has_disp && name == disp_name && new_index != 0) {
          save_index = new_index;
          new_index = 0;
        }

        Ioss::Field field = ge->get_field(name);
        const Ioss::VariableType *var_type = field.transformed_storage();

        int re_im = 1;
        if (field.get_type() == Ioss::Field::COMPLEX)
          re_im = 2;
        for (int complex_comp = 0; complex_comp < re_im; complex_comp++) {
          std::string field_name = field.get_name();
          if (re_im == 2) {
            field_name += complex_suffix[complex_comp];
          }

          char field_suffix_separator = get_field_separator();
          for (int i=1; i <= var_type->component_count(); i++) {
            std::string var_string = var_type->label_name(field_name, i, field_suffix_separator);

            // Add to 'VariableNameMap variables' so can determine
            // exodusII index given a Sierra field name.  exodusII index
            // is just 'i+1'
            if (reduction || type == EX_GLOBAL) {
              // If this is not a global (region) variable, need to prepend the block name
              // to avoid name collisions...
              if (type != EX_GLOBAL) {
                var_string = ge->name() + ":" + var_string;
              }
            }

            if (variables.find(var_string) == variables.end()) {
              variables.insert(VNMValuePair(var_string, ++new_index));
            }
          }
        }
        if (has_disp && name == disp_name) {
          new_index = save_index;
        }
      }
      return new_index;
    }

      // common
    void DatabaseIO::generate_sideset_truth_table()
    {
      size_t var_count = m_variables[EX_SIDE_SET].size();

      if (var_count == 0 || m_groupCount[EX_SIDE_SET] == 0)
        return;

      // Member variable.  Will be deleted in destructor...
      m_truthTable[EX_SIDE_SET].resize(m_groupCount[EX_SIDE_SET]*var_count);

      // Fill in the truth table.  It is conceptually a two-dimensional array of
      // the form 'array[num_blocks][num_var]'.  In C++,
      // the values for the first block are first, followed by
      // next block, ...
      size_t offset = 0;

      Ioss::SideSetContainer sidesets = get_region()->get_sidesets();
      Ioss::SideSetContainer::const_iterator I;

      char field_suffix_separator = get_field_separator();
      for (I=sidesets.begin(); I != sidesets.end(); ++I) {
        Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
        Ioss::SideBlockContainer::const_iterator J;

        for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
          // See if this sideblock has a corresponding entry in the sideset list.
          if ((*J)->property_exists("invalid"))
            continue;

          // Get names of all transient and reduction fields...
          Ioss::NameList results_fields;
          (*J)->field_describe(Ioss::Field::TRANSIENT, &results_fields);
          (*J)->field_describe(Ioss::Field::REDUCTION, &results_fields);

          Ioss::NameList::const_iterator IF;
          for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
            std::string field_name = *IF;

            Ioss::Field field = (*J)->get_field(field_name);
            const Ioss::VariableType *var_type = field.transformed_storage();
            Ioss::Field::BasicType ioss_type = field.get_type();

            int re_im = 1;
            if (ioss_type == Ioss::Field::COMPLEX)
              re_im = 2;
            for (int complex_comp = 0; complex_comp < re_im; complex_comp++) {
              field_name = field.get_name();
              if (re_im == 2) {
                field_name += complex_suffix[complex_comp];
              }

              for (int i=1; i <= var_type->component_count(); i++) {
                std::string var_string = var_type->label_name(field_name, i, field_suffix_separator);
                // Find position of 'var_string' in 'm_variables[]'
                VariableNameMap::iterator VN = m_variables[EX_SIDE_SET].find(var_string);
                if (VN != m_variables[EX_SIDE_SET].end()) {
                  // Index '(*VN).second' is 1-based...
                  m_truthTable[EX_SIDE_SET][offset + (*VN).second-1] = 1;
                }
              }
            }
          }
        }
        offset += var_count;
      }
      assert(offset == var_count * m_groupCount[EX_SIDE_SET]);
    }

      // common
    void DatabaseIO::output_results_names(ex_entity_type type,
					  VariableNameMap &variables) const
    {
      bool lowercase_names = (properties.exists("VARIABLE_NAME_CASE") &&
                              Ioss::Utils::lowercase(properties.get("VARIABLE_NAME_CASE").get_string()) == "lower");
      bool uppercase_names = (properties.exists("VARIABLE_NAME_CASE") &&
                              Ioss::Utils::lowercase(properties.get("VARIABLE_NAME_CASE").get_string()) == "upper");

      size_t var_count = variables.size();

      if (var_count > 0) {
        // Push into a char** array...
        std::vector<char*> var_names(var_count);
        std::vector<std::string> variable_names(var_count);
        VariableNameMap::const_iterator J  = variables.begin();
        VariableNameMap::const_iterator JE = variables.end();
        while (J != JE) {
          size_t index = (*J).second;
          assert(index > 0 && index <= var_count);
          variable_names[index-1] = (*J).first;
          if (uppercase_names)
            variable_names[index-1] = Ioss::Utils::uppercase(variable_names[index-1]);
          else if (lowercase_names)
            variable_names[index-1] = Ioss::Utils::lowercase(variable_names[index-1]);
          var_names[index-1] = (char*)variable_names[index-1].c_str();
          ++J;
        }

        int ierr = ex_put_variable_names(get_file_pointer(), type, var_count, TOPTR(var_names));
        if (ierr < 0)
          Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
    }

      // common
    // Handle special output time requests -- primarily restart (cycle, overwrite)
    // Given the global region step, return the step on the database...
    int DatabaseIO::get_database_step(int global_step) const
    {
      assert(overlayCount >= 0 && cycleCount >= 0);
      if (overlayCount == 0 && cycleCount == 0)
        return global_step;

      int local_step = global_step - 1;
      local_step /= (overlayCount + 1);
      if (cycleCount > 0)
        local_step %= cycleCount;
      return local_step + 1;
    }

      // common
    void DatabaseIO::finalize_write(double sim_time)
    {
      // Attempt to ensure that all data written up to this point has
      // actually made it out to disk.  We also write a special attribute
      // to the file to indicate that the current timestep should be
      // complete on the disk.
      // The attribute is a GLOBAL attribute named "last_written_time"
      // which is a double value which can be compared to the values in
      // the time array to make sure they match.  If they don't, then
      // hopefully the "last_written_time" is smaller than the time
      // array value and indicates that the last step is corrupt.

      // Update the attribute.
      Ioex::update_last_time_attribute(get_file_pointer(), sim_time);

      // Flush the files buffer to disk...
      // If a history file, then only flush if there is more
      // than 10 seconds since the last flush to avoid
      // the flush eating up cpu time for small fast jobs...
      // NOTE: If decide to do this on all files, need to sync across
      // processors to make sure they all flush at same time.

      // GDS: 2011/03/30 -- Use for all non-parallel files, but shorten
      // time for non history files.  Assume that can afford to lose ~10
      // seconds worth of data...  (Flush was taking long time on some
      // /scratch filesystems at SNL for short regression tests with
      // lots of steps)
      // GDS: 2011/07/27 -- shorten from 90 to 10.  Developers running
      // small jobs were not able to view output until job
      // finished. Hopefully the netcdf no-fsync fix along with this fix
      // results in negligible impact on runtime with more syncs.

      bool do_flush = true;
      if (dbUsage == Ioss::WRITE_HISTORY || !isParallel) {
        assert (myProcessor == 0);
        time_t cur_time = time(NULL);
        if (cur_time - timeLastFlush >= 10) {
          timeLastFlush = cur_time;
          do_flush = true;
        } else {
          do_flush = false;
        }
      }
      if (do_flush)
        ex_update(get_file_pointer());
    }

      // common
    void Ioex::DatabaseIO::add_attribute_fields(ex_entity_type entity_type,
                                                Ioss::GroupingEntity *block,
                                                int attribute_count, const std::string& type)
    {
      /// \todo REFACTOR Some of the attribute knowledge should be at
      /// the Ioss::ElementTopology level instead of here. That would
      /// make it easier for an application to register a new element
      /// type and its attributes.

      // Attribute "Conventions" to be used if there are no attribute names on the database:
      // from Table 1 in ExodusII manual
      //
      // Circle     1     Radius [Volume]
      // Sphere     1     Radius [Volume]
      // Truss      1     Area
      // 2D Beam    3     Area, I, J
      // 3D Beam    7     Area, I1, I2, J, V1, V2, V3 (V will be a 3D vector named "reference_axis")
      // Shell      1     Thickness
      //
      // Additional conventions not defined in ExodusII manual:
      // * If a "beam" has 1 attribute, call it "area"
      // * Treat "bar" and "rod" as aliases for "truss"
      // * Treat "trishell" as alias for "shell"
      // * All "shell" or "trishell" elements -- If #attributes == #node/element, the
      //                                         attribute is "nodal_thickness"
      //
      // If there are attribute names on the database, use those names.
      // Always create a variable "attribute" which contains a single
      // field for all attributes...

      assert(block != NULL);
      if (attribute_count > 0) {
        std::string block_name = block->name();
        size_t my_element_count = block->get_property("entity_count").get_int();

        // Get the attribute names. May not exist or may be blank...
        char **names = get_exodus_names(attribute_count, maximumNameLength);
        int64_t id = block->get_property("id").get_int();
        
        // Some older applications do not want to used named
        // attributes; in this case, just create a field for each
        // attribute named "attribute_1", "attribute_2", ..., "attribute_#"
        // This is controlled by the database property
        // "IGNORE_ATTRIBUTE_NAMES"
        char field_suffix_separator = get_field_separator();
        bool attributes_named = true; // Possibly reset below; note that even if ignoring
                                      // attribute names, they are still 'named'

        if (properties.exists("IGNORE_ATTRIBUTE_NAMES")) {
          field_suffix_separator = ' '; // Do not combine into a
          // higher-order storage type.
          
          for (int i=0; i < attribute_count; i++) {
            int writ = ::snprintf(names[i], maximumNameLength+1, "attribute_%d", i+1);
            if (writ > maximumNameLength)
              names[i][maximumNameLength] = '\0';
          }
        }
        else {
          // Use attribute names if they exist.
          {
            Ioss::SerializeIO   serializeIO__(this);
            if (block->get_property("entity_count").get_int() != 0) {
              int ierr = ex_get_attr_names(get_file_pointer(), entity_type, id, &names[0]);
              if (ierr < 0)
                Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }
          }

          // Sync names across processors...
          if (isParallel) {
            std::vector<char> cname(attribute_count * (maximumNameLength+1));
            if (block->get_property("entity_count").get_int() != 0) {
              for (int i=0; i < attribute_count; i++) {
                std::memcpy(&cname[i*(maximumNameLength+1)], names[i], maximumNameLength+1);
              }
            }
            util().attribute_reduction(attribute_count * (maximumNameLength+1), TOPTR(cname));
            for (int i=0; i < attribute_count; i++) {
              std::memcpy(names[i], &cname[i*(maximumNameLength+1)], maximumNameLength+1);
            }
          }
          
          // Convert to lowercase.
          attributes_named = true;
          for (int i=0; i < attribute_count; i++) {
            fix_bad_name(names[i]);
            Ioss::Utils::fixup_name(names[i]);
            if (names[i][0] == '\0' || names[i][0] == ' ' || !std::isalnum(names[i][0])) {
              attributes_named = false;
            }
          }
        }

        if (attributes_named) {
          std::vector<Ioss::Field> attributes;
	  Ioex::get_fields(my_element_count, names, attribute_count,
			   Ioss::Field::ATTRIBUTE, field_suffix_separator, NULL,
			   attributes);
          int offset = 1;
          std::vector<Ioss::Field>::const_iterator IF;
          for (IF = attributes.begin(); IF != attributes.end(); ++IF) {
            Ioss::Field field = *IF;
            if (block->field_exists(field.get_name())) {
              std::ostringstream errmsg;
              errmsg << "ERROR: In block '" << block->name() << "', attribute '" << field.get_name()
                     << "' is defined multiple times which is not allowed.\n";
              IOSS_ERROR(errmsg);
            }
            block->field_add(field);
            const Ioss::Field &tmp_field = block->get_fieldref(field.get_name());
            tmp_field.set_index(offset);
            offset += field.raw_storage()->component_count();
          }
        }
        else {
          // Attributes are not named....
          // Try to assign some meaningful names based on conventions...
          std::string att_name = "attribute";  // Default
          int unknown_attributes = 0;

          if (type_match(type, "shell") || type_match(type, "trishell")) {
            if (attribute_count == block->get_property("topology_node_count").get_int()) {
              att_name = "nodal_thickness";

              std::string storage = "Real[";
              storage += Ioss::Utils::to_string(attribute_count);
              storage += "]";

              block->field_add(Ioss::Field(att_name, Ioss::Field::REAL, storage,
                                           Ioss::Field::ATTRIBUTE, my_element_count, 1));
            } else {
              att_name = "thickness";
              block->field_add(Ioss::Field(att_name, Ioss::Field::REAL, SCALAR(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, 1));
              unknown_attributes = attribute_count - 1;
            }

          }

          // NOTE: This must appear before the "sphere" check since
          // sphere is substring of "sphere-mass"
          // Want an exact match here, not substring match...
          else if (Ioss::Utils::case_strcmp(type, "sphere-mass") == 0) {
            if (attribute_count != 10) {
              if (myProcessor == 0) {
                IOSS_WARNING << "For element block '" << block->name()
                             << "' of type '" << type << "' there were "
                             << attribute_count
                             << " attributes instead of the expected 10 attributes "
                             << "known to the IO Subsystem. "
                             << " The attributes can be accessed as the field named 'attribute'";
              }
            } else {
              // First attribute is concentrated mass...
              size_t offset = 1;
              block->field_add(Ioss::Field("mass", Ioss::Field::REAL, SCALAR(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, offset));
              offset += 1;

              // Next six attributes are moment of inertia -- symmetric tensor
              block->field_add(Ioss::Field("inertia", Ioss::Field::REAL, SYM_TENSOR(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, offset));
              offset += 6;

              // Next three attributes are offset from node to CG
              block->field_add(Ioss::Field("offset", Ioss::Field::REAL, VECTOR3D(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, offset));
            }
          }

          else if (type_match(type, "circle") || type_match(type, "sphere")) {
            att_name = "radius";
            size_t offset = 1;
            block->field_add(Ioss::Field(att_name, Ioss::Field::REAL, SCALAR(),
                                         Ioss::Field::ATTRIBUTE, my_element_count, offset++));
            if (attribute_count > 1) {
              // Default second attribute (from sphgen3d) is "volume"
              // which is the volume of the cube which would surround a
              // sphere of the given radius.
              att_name = "volume";
              block->field_add(Ioss::Field(att_name, Ioss::Field::REAL, SCALAR(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, offset++));
            }
            unknown_attributes = attribute_count - 2;
          }

          else if (type_match(type, "truss") ||
                   type_match(type, "bar")   ||
                   type_match(type, "beam")   ||
                   type_match(type, "rod")) {
            // Technically, truss, bar, rod should all only have 1 attribute; however,
            // there are some mesh generation codes that treat all of these types the
            // same and put "beam-type" attributes on bars...
            int index = 1;
            att_name = "area";
            block->field_add(Ioss::Field(att_name, Ioss::Field::REAL, SCALAR(),
                                         Ioss::Field::ATTRIBUTE, my_element_count, index++));

            if (spatialDimension == 2 && attribute_count >= 3) {
              block->field_add(Ioss::Field("i", Ioss::Field::REAL, SCALAR(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, index++));
              block->field_add(Ioss::Field("j", Ioss::Field::REAL, SCALAR(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, index++));
            }
            else if (spatialDimension == 3 && attribute_count >= 7) {
              block->field_add(Ioss::Field("i1", Ioss::Field::REAL, SCALAR(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, index++));
              block->field_add(Ioss::Field("i2", Ioss::Field::REAL, SCALAR(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, index++));
              block->field_add(Ioss::Field("j", Ioss::Field::REAL, SCALAR(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, index++));
              block->field_add(Ioss::Field("reference_axis", Ioss::Field::REAL, VECTOR3D(),
                                           Ioss::Field::ATTRIBUTE, my_element_count, index));
              index += 3;
              if (attribute_count >= 10) {
                // Next three attributes would (hopefully) be offset vector...
                // This is typically from a NASGEN model.
                block->field_add(Ioss::Field("offset", Ioss::Field::REAL, VECTOR3D(),
                                             Ioss::Field::ATTRIBUTE, my_element_count, index));
                index += 3;
              }
            }
            unknown_attributes = attribute_count - (index-1);
          }

          else {
            unknown_attributes = attribute_count;
          }

          if (unknown_attributes > 0) {
            att_name = "extra_attribute_";
            att_name += Ioss::Utils::to_string(unknown_attributes);
            std::string storage = "Real[";
            storage += Ioss::Utils::to_string(unknown_attributes);
            storage += "]";
            size_t index = attribute_count - unknown_attributes + 1;
            block->field_add(Ioss::Field(att_name, Ioss::Field::REAL, storage,
                                         Ioss::Field::ATTRIBUTE, my_element_count, index));
          }
        }

        // Always create a field called "attribute" containing data
        // for all attributes on the mesh
        std::string att_name = "attribute";  // Default
        std::string storage = "Real[";
        storage += Ioss::Utils::to_string(attribute_count);
        storage += "]";

        block->field_add(Ioss::Field(att_name, Ioss::Field::REAL, storage,
                                     Ioss::Field::ATTRIBUTE, my_element_count, 1));

        // Release memory...
        delete_exodus_names(names, attribute_count);
      }
    }

  void DatabaseIO::output_other_meta_data()
  {
    Ioss::SerializeIO       serializeIO__(this);
    
    if (!properties.exists("OMIT_QA_RECORDS")) {
      put_qa();
    }
    if (!properties.exists("OMIT_INFO_RECORDS")) {
      put_info();
    }

    // Write attribute names (if any)...
    char field_suffix_separator = get_field_separator();
    write_attribute_names(get_file_pointer(), EX_NODE_SET,   get_region()->get_nodesets(),       field_suffix_separator);
    write_attribute_names(get_file_pointer(), EX_EDGE_SET,   get_region()->get_edgesets(),       field_suffix_separator);
    write_attribute_names(get_file_pointer(), EX_FACE_SET,   get_region()->get_facesets(),       field_suffix_separator);
    write_attribute_names(get_file_pointer(), EX_ELEM_SET,   get_region()->get_elementsets(),    field_suffix_separator);
    write_attribute_names(get_file_pointer(), EX_NODE_BLOCK, get_region()->get_node_blocks(),    field_suffix_separator);
    write_attribute_names(get_file_pointer(), EX_EDGE_BLOCK, get_region()->get_edge_blocks(),    field_suffix_separator);
    write_attribute_names(get_file_pointer(), EX_FACE_BLOCK, get_region()->get_face_blocks(),    field_suffix_separator);
    write_attribute_names(get_file_pointer(), EX_ELEM_BLOCK, get_region()->get_element_blocks(), field_suffix_separator);
    
    // Write coordinate names...
    char const *labels[3];
    labels[0] = "x";
    labels[1] = "y";
    labels[2] = "z";
    int ierr = ex_put_coord_names(get_file_pointer(), (char**)labels);
    if (ierr < 0)
      Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
    
    // Write coordinate frame data...
    write_coordinate_frames(get_file_pointer(), get_region()->get_coordinate_frames());
  }
  } // End of namespace

  namespace {
    // common
    template <typename T>
    void write_attribute_names(int exoid, ex_entity_type type, const std::vector<T*>& entities,
                               const char suffix_separator)
    {
      // For the entity, determine the attribute fields and the correct
      // order. Write the names of these fields.  However, be aware that
      // the field "attribute" always exists to contain all attributes
      // and its name should not be used even if it is the only
      // attribute field.
      typename std::vector<T*>::const_iterator I;
      for (I=entities.begin(); I != entities.end(); ++I) {
        Ioss::GroupingEntity *ge = *I;

        std::string ge_name = ge->name();
        int attribute_count = ge->get_property("attribute_count").get_int();
        if (attribute_count > 0) {

          check_attribute_index_order(ge);

          std::vector<char*> names(attribute_count);
          std::vector<std::string> names_str(attribute_count);

          // Get the attribute fields...
          Ioss::NameList results_fields;
          ge->field_describe(Ioss::Field::ATTRIBUTE, &results_fields);

          Ioss::NameList::const_iterator IF;
          for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
            std::string field_name = *IF;
            const Ioss::Field &field = ge->get_fieldref(field_name);
            assert(field.get_index() != 0);

            if (field_name == "attribute") {
              field.set_index(1);
              continue;
            }

            const Ioss::VariableType *vtype = field.raw_storage();
            int comp_count = vtype->component_count();
            int field_offset = field.get_index();
            for (int i=0; i < comp_count; i++) {
              names_str[field_offset-1+i] = vtype->label_name(field_name, i+1, suffix_separator);
              names[field_offset-1+i] = (char*)names_str[field_offset-1+i].c_str();
            }
          }
          size_t ge_id = ge->get_property("id").get_int();
          int ierr = ex_put_attr_names(exoid, type, ge_id, TOPTR(names));
          if (ierr < 0)
            Ioex::exodus_error(exoid, __LINE__, -1);
        }
      }
    }

    // common
    void check_attribute_index_order(Ioss::GroupingEntity *block)
    {
      int attribute_count = block->get_property("attribute_count").get_int();
      if (attribute_count == 0)
        return;
      int component_sum = 0;

      std::vector<int> attributes(attribute_count+1);

      // Get the attribute fields...
      Ioss::NameList results_fields;
      block->field_describe(Ioss::Field::ATTRIBUTE, &results_fields);

      bool all_attributes_indexed = true;
      bool some_attributes_indexed = false;

      Ioss::NameList::const_iterator IF;
      for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
        std::string field_name = *IF;
        const Ioss::Field &field = block->get_fieldref(field_name);

        if (field_name == "attribute") {
          field.set_index(1);
          if (results_fields.size() == 1) {
            return;
          }
          continue;
        }

        int field_offset = field.get_index();
        if (field_offset == 0) {
          all_attributes_indexed = false;
        } else {
          some_attributes_indexed = true;
        }

        const Ioss::VariableType *type = field.raw_storage();
        int comp_count = type->component_count();
        component_sum += comp_count;

        if (field_offset == 0)
          continue;

        if (field_offset + comp_count - 1 > attribute_count) {
          std::ostringstream errmsg;
          errmsg << "INTERNAL ERROR: For block '" << block->name() << "', attribute '" << field_name
                 << "', the indexing is incorrect.\n"
                 << "Something is wrong in the Ioex::DatabaseIO class, function check_attribute_index_error. Please report.\n";
          IOSS_ERROR(errmsg);
        }

        for (int i=field_offset; i < field_offset+comp_count; i++) {
          if (attributes[i] != 0) {
            std::ostringstream errmsg;
            errmsg << "INTERNAL ERROR: For block '" << block->name() << "', attribute '" << field_name
                   << "', indexes into the same location as a previous attribute.\n"
                   << "Something is wrong in the Ioex::DatabaseIO class, function check_attribute_index_error. Please report.\n";
            IOSS_ERROR(errmsg);
          } else {
            attributes[i] = 1;
          }
        }
      }

      if (component_sum > attribute_count) {
        std::ostringstream errmsg;
        errmsg << "INTERNAL ERROR: Block '" << block->name() << "' is supposed to have " << attribute_count
               << " attributes, but " << component_sum << " attributes were counted.\n"
               << "Something is wrong in the Ioex::DatabaseIO class, function check_attribute_index_error. Please report.\n";
        IOSS_ERROR(errmsg);
      }

      // Take care of the easy cases first...
      if (all_attributes_indexed) {
        // Check that all attributes are defined.  This should have
        // caught above in the duplicate index check.
        for (int i=1; i <= attribute_count; i++) {
          if (attributes[i] == 0) {
            std::ostringstream errmsg;
            errmsg << "INTERNAL ERROR: Block '" << block->name() << "' has an incomplete set of attributes.\n"
                   << "Something is wrong in the Ioex::DatabaseIO class, function check_attribute_index_error. Please report.\n";
            IOSS_ERROR(errmsg);
          }
        }
        return;
      }

      if (!all_attributes_indexed && !some_attributes_indexed) {
        // Index was not set for any of the attributes; set them all...
        size_t offset = 1;
        for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
          std::string field_name = *IF;
          const Ioss::Field &field = block->get_fieldref(field_name);

          if (field_name == "attribute") {
            field.set_index(1);
            continue;
          }

          const Ioss::VariableType *type = field.raw_storage();
          int comp_count = type->component_count();

          assert(field.get_index() == 0);
          field.set_index(offset);
          offset += comp_count;
        }
        assert((int)offset == attribute_count+1);
        return;
      }

      // At this point, we have a partially indexed set of attributes.  Some have an index and some don't
      // The easy case is if the missing indices are at the end of the list...
      assert(!all_attributes_indexed && some_attributes_indexed);
      int last_defined = 0;
      for (int i=1; i < attribute_count+1; i++) {
        if (attributes[i] != 0)
          last_defined = i;
      }
      int first_undefined = attribute_count;
      for (int i=attribute_count; i > 0; i--) {
        if (attributes[i] == 0)
          first_undefined = i;
      }
      if (last_defined < first_undefined) {
        for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
          std::string field_name = *IF;
          const Ioss::Field &field = block->get_fieldref(field_name);

          if (field_name == "attribute") {
            field.set_index(1);
            continue;
          }

          if (field.get_index() == 0) {
            field.set_index(first_undefined);
            const Ioss::VariableType *type = field.raw_storage();
            int comp_count = type->component_count();
            first_undefined += comp_count;
          }
        }
        assert(first_undefined == attribute_count+1);
        return;
      }

      // Take the easy way out... Just reindex all attributes.
      size_t offset = 1;
      for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
        std::string field_name = *IF;
        const Ioss::Field &field = block->get_fieldref(field_name);

        if (field_name == "attribute") {
          field.set_index(1);
          continue;
        }

        const Ioss::VariableType *type = field.raw_storage();
        int comp_count = type->component_count();

        assert(field.get_index() == 0);
        field.set_index(offset);
        offset += comp_count;
      }
      assert((int)offset == attribute_count+1);
      return;
    }

#ifndef NDEBUG
    template <typename T>
    bool check_block_order(const std::vector<T*> &blocks)
    {
      // Verify that element blocks are defined in sorted offset order...
      typename std::vector<T*>::const_iterator I;

      int64_t eb_offset = -1;
      for (I=blocks.begin(); I != blocks.end(); ++I) {
        int64_t this_off = (*I)->get_offset();
        if (this_off < eb_offset)
          return false;
        eb_offset = this_off;
      }
      return true;
    }
#endif

    void check_variable_consistency(const ex_var_params &exo_params,
                                    int my_processor, const std::string& filename,
                                    const Ioss::ParallelUtils &util)
    {
#ifdef HAVE_MPI    
      const int num_types=10;
      std::vector<int> var_counts(num_types);
      var_counts[0] = exo_params.num_glob;
      var_counts[1] = exo_params.num_node;
      var_counts[2] = exo_params.num_edge;
      var_counts[3] = exo_params.num_face;
      var_counts[4] = exo_params.num_elem;
      var_counts[5] = exo_params.num_nset;
      var_counts[6] = exo_params.num_eset;
      var_counts[7] = exo_params.num_fset;
      var_counts[8] = exo_params.num_sset;
      var_counts[9] = exo_params.num_elset;

      Ioss::IntVector all_counts;
      util.gather(var_counts, all_counts);

      bool any_diff = false;
      std::ostringstream errmsg;
      if (my_processor == 0) {
        bool diff[num_types];
        // See if any differ...
        for (int iv = 0; iv < 10; iv++) {
          diff[iv] = false;
          std::string type;
          switch (iv) {
          case 0:  type = "global";     break;
          case 1:  type = "nodal";      break;
          case 2:  type = "edge";             break;
          case 3:  type = "face";             break;
          case 4:  type = "element";    break;
          case 5:  type = "nodeset";    break;
          case 6:  type = "edgeset";    break;
          case 7:  type = "faceset";    break;
          case 8:  type = "sideset";    break;
          case 9:  type = "elementset"; break;
          }

          for (int ip = 1; ip < util.parallel_size(); ip++) {
            if (var_counts[iv] != all_counts[ip*num_types+iv]) {
              any_diff = true;
              if (!diff[iv]) {
                Ioss::FileInfo db(filename);
                diff[iv] = true;
                errmsg << "\nERROR: Number of " << type
                       << " variables is not consistent on all processors.\n"
                       << "       Database: " << db.tailname() << "\n"
                       << "\tProcessor 0 count = " << var_counts[iv] << "\n";
              }
              errmsg << "\tProcessor " << ip << " count = " << all_counts[ip*num_types+iv] << "\n";
            }
          }
        }
      } else {
        // Give the other processors something to say...
        errmsg << "ERROR: Variable type counts are inconsistent. See processor 0 output for more details.\n";
      }
      int idiff = any_diff ? 1 : 0;
      MPI_Bcast(&idiff, 1, MPI_INT, 0, util.communicator());
      any_diff = idiff == 1;

      if (any_diff) {
        std::runtime_error x(errmsg.str());
        throw x;
      }
#endif
    }
  }

