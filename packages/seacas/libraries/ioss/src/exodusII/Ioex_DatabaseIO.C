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
#include <exodusII/Ioex_DatabaseIO.h>
#include <exodusII/Ioex_Internals.h>
#include <tokenize.h>

#include <string>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <iterator>
#include <time.h>
#include <float.h>

#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_SerializeIO.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_FileInfo.h>
#include <Ioss_SurfaceSplit.h>

#include <exodusII.h>
#include <ne_nemesisI.h>
#include <assert.h>

namespace Ioex {
  // ========================================================================
  // Internal typedefs/structs/classes used in some algorithms.
  // ========================================================================
  typedef std::set<std::string> SideSetSet;
  typedef std::map<std::string, const std::string, std::less<const std::string> > SideSetMap;

  struct TopologyMapCompare {
    bool operator() (const std::pair<std::string, const Ioss::ElementTopology*> &lhs,
                     const std::pair<std::string, const Ioss::ElementTopology*> &rhs) const
    {
      assert(lhs.second != NULL);
      assert(rhs.second != NULL);
      return lhs.first < rhs.first || (!(rhs.first < lhs.first) && lhs.second->name() < rhs.second->name());
    }
  };

  typedef std::map<std::pair<std::string, const Ioss::ElementTopology*>, int, TopologyMapCompare > TopologyMap;
  typedef TopologyMap::value_type TopoMapPair;
}


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

  const char *Version() {return "Ioex_DatabaseIO.C 2011/04/14 gdsjaar";}

  bool type_match(const std::string& type, const char *substring);
  int extract_id(const std::string &name_id);
  bool set_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Ioex::EntityIdSet *idset);
  int  get_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Ioex::EntityIdSet *idset);
  void decode_surface_name(Ioex::SideSetMap &fs_map, Ioex::SideSetSet &fs_set, const std::string &name);
  void fix_bad_name(char* name);

  void exodus_error(int exoid, int lineno, int /* processor */) {
    std::ostringstream errmsg;
    // Create errmsg here so that the exerrval doesn't get cleared by
    // the ex_close call.
    // Try to interpret exodus error messages...
    std::string error_type;
    switch (exerrval) {
    case -31:
      error_type = "System Error -- Usually disk full or filesystem issue"; break;
    case -33:
      error_type = "Not a netcdf id"; break;
    case -34:
      error_type = "Too many files open"; break;
    case -41:
    case -44:
    case -48:
    case -53:
    case -62:
      error_type = "Internal netcdf/exodusII dimension exceeded"; break;
    case -51:
      error_type = "Not an exodusII/netcdf file"; break;
    case -59:
      error_type = "Attribute of variable name contains illegal characters"; break;
    case -60:
      error_type = "Memory allocation (malloc) failure"; break;
    case -64:
      error_type = "Filesystem issue; File likely truncated or possibly corrupted"; break;
    default:
      ;
    }
    errmsg << "ExodusII error (" << exerrval << ")" << error_type << " at line " << lineno
	   << " in file '" << Version()
	   << "' Please report to gdsjaar@sandia.gov if you need help.";

    ex_err(NULL, NULL, EX_PRTLASTMSG);
    if (exoid > 0)
      ex_close(exoid);
    IOSS_ERROR(errmsg);
  }

  void check_non_null(void *ptr, const char *type, const std::string &name)
  {
    if (ptr == NULL) {
      std::ostringstream errmsg;
      errmsg << "INTERNAL ERROR: Could not find " << type << " '" << name << "'."
		<< " Something is wrong in the Ioex::DatabaseIO class. Please report.\n";
      IOSS_ERROR(errmsg);
    }
  }

  bool block_is_omitted(Ioss::ElementBlock *block) {
    bool omitted = false;
    if (block->property_exists("omitted"))
      omitted = (block->get_property("omitted").get_int() == 1);
    return omitted;
  }

  void filter_element_list(Ioss::Region *region,
			   Ioss::IntVector &elements, Ioss::IntVector &sides,
			   bool remove_omitted_elements);
  
  void set_sideblock_ids(Ioss::SideSetContainer &ssets,
			 std::vector<Ioex::SideSet> &sidesets,
			 const Ioss::ParallelUtils &util,
			 Ioex::EntityIdSet *idset);

  void separate_surface_element_sides(Ioss::IntVector &element,
				      Ioss::IntVector &sides,
				      Ioss::Region *region,
				      Ioex::TopologyMap &topo_map,
				      Ioex::TopologyMap &side_map,
				      Ioss::SurfaceSplitType split_type);
    
  bool find_displacement_field(Ioss::NameList &fields,
			       const Ioss::GroupingEntity *block,
			       int ndim, std::string *disp_name);

  void get_fields(int entity_count, char** names, size_t num_names,
		  Ioss::Field::RoleType fld_role,
		  const char suffix_separator, int *local_truth,
		  std::vector<Ioss::Field> &fields);

  void add_map_fields(int exoid, Ioss::ElementBlock *block, int my_element_count);

  template <typename T>
  bool check_block_order(const std::vector<T*> &blocks);

  void check_variable_consistency(const ex_var_params &exo_params,
				  int my_processor, const std::string &filename,
				  const Ioss::ParallelUtils &util);

  char ** get_exodus_names(size_t count, int size)
  {
    char **names = new char* [count];
    for (size_t i=0; i < count; i++) {
      names[i] = new char [size+1];
      std::memset(names[i], '\0', size+1);
    }
    return names;
  }

  void delete_exodus_names(char **names, int count)
  {
    for (int i=0; i < count; i++) {delete [] names[i];}
    delete [] names;
  }

  std::string get_entity_name(int exoid, ex_entity_type type, int id,
			      const std::string &basename, int length)
  {
    std::vector<char> buffer(length+1);
    buffer[0] = '\0';
    int error = ex_get_name(exoid, type, id, TOPTR(buffer));
    if (error < 0)
      exodus_error(exoid, __LINE__, -1);
    if (buffer[0] != '\0') {
      Ioss::Utils::fixup_name(TOPTR(buffer));
      return (std::string(TOPTR(buffer)));
    } else {
      return Ioss::Utils::encode_entity_name(basename, id);
    }
  }
  
  void check_attribute_index_order(Ioss::GroupingEntity *block);

  template <typename T>
  void write_attribute_names(int exoid, ex_entity_type type, const std::vector<T*>& entities,
			     const char suffix_separator);

  template <typename T>
  void generate_block_truth_table(Ioex::VariableNameMap &variables,
				  Ioex::IntVector &truth_table,
				  std::vector<T*> &blocks,
				  char field_suffix_separator);
}

namespace Ioex {
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
			 Ioss::DatabaseUsage db_usage, MPI_Comm communicator) :
    Ioss::DatabaseIO(region, filename, db_usage, communicator),
    exodusFilePtr(-1), databaseTitle(""), exodusMode(EX_CLOBBER),
    maximumNameLength(32), spatialDimension(0),
    nodeCount(0), edgeCount(0), faceCount(0), elementCount(0),
    commsetNodeCount(0), commsetElemCount(0),
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
#if 0
    ex_opts(EX_VERBOSE|EX_DEBUG);
#endif
    
    if (!is_input()) {
      if (util().get_environment("EX_MINIMIZE_OPEN_FILES", isParallel)) {
	std::cerr << "IOEX: Minimizing open files because EX_MINIMIZE_OPEN_FILES environment variable is set.\n";
	minimizeOpenFiles = true;
      }

      if (util().get_environment("EX_MODE", exodusMode, isParallel)) {
	std::cerr << "IOEX: Exodus create mode set to " << exodusMode
		  << " from value of EX_MODE environment variable.\n";
      }
    }

    // Don't open output files until they are actually going to be
    // written to.  This is needed for proper support of the topology
    // files and auto restart so we don't overwrite a file with data we
    // need to save...
  }

  DatabaseIO::~DatabaseIO()
  {
    try {
      free_file_pointer();
    } catch (...) {
    }
  }

  unsigned DatabaseIO::entity_field_support() const
  {
    return Ioss::NODEBLOCK | Ioss::EDGEBLOCK | Ioss::FACEBLOCK | Ioss::ELEMENTBLOCK |
      Ioss::NODESET   | Ioss::EDGESET   | Ioss::FACESET   | Ioss::ELEMENTSET   |
      Ioss::SIDESET   | Ioss::SIDEBLOCK | Ioss::REGION    | Ioss::SUPERELEMENT;
  }

  bool DatabaseIO::ok(bool write_message) const
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
    std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
    
    int exodus_file_ptr = ex_open(decoded_filename.c_str(), EX_READ,
				  &cpu_word_size, &io_word_size, &version);
    
    if (!is_input() && exodus_file_ptr < 0) {
      // File didn't exist above, but this OK if is an output
      // file. See if we can create it...
      io_word_size = cpu_word_size;
      exodus_file_ptr = ex_create(decoded_filename.c_str(), exodusMode,
				  &cpu_word_size, &io_word_size);
    }
    
    // Check for valid exodus_file_ptr (valid >= 0; invalid < 0)
    int global_file_ptr = util().global_minmax(exodus_file_ptr, Ioss::ParallelUtils::DO_MIN);
    if (write_message && global_file_ptr < 0) {
      // See which processors could not open/create the file...
      Ioss::IntVector status;
      util().gather(exodus_file_ptr, status);
      if (myProcessor == 0) {
	std::string open_create = is_input() ? "open input" : "create output";
	bool first = true;
	std::ostringstream errmsg;
	errmsg << "ERROR: Unable to " << open_create << " database '" << get_filename() << "' of type 'exodusII'";
	if (isParallel) {
	  errmsg << "\n       on processor(s): ";
	  for (int i=0; i < util().parallel_size(); i++) {
	    if (status[i] < 0) {
	      if (!first) errmsg << ", ";
	      errmsg << i;
	      first = false;
	    }
	  }
	}
	errmsg << "\n";
	std::cerr << errmsg.str();
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

    if (Ioss::SerializeIO::isEnabled() && !Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Process " << Ioss::SerializeIO::getRank()
	     << " is attempting to do I/O without serialized I/O";
      IOSS_ERROR(errmsg);
    }

    if (Ioss::SerializeIO::isEnabled() && !Ioss::SerializeIO::inMyGroup()) {
      std::ostringstream errmsg;
      errmsg << "Process " << Ioss::SerializeIO::getRank()
	     << " is attempting to do I/O while " << Ioss::SerializeIO::getOwner()
	     << " owns the token";
      IOSS_ERROR(errmsg);
    }


    if (exodusFilePtr < 0) {
      int cpu_word_size = sizeof(double);
      int io_word_size  = 0;
      float version;
      std::string decoded_filename = util().decode_filename(get_filename(),
							    isParallel);

      if (is_input()) {
	exodusFilePtr = ex_open(decoded_filename.c_str(), EX_READ,
				&cpu_word_size, &io_word_size, &version);
      } else {
	if (fileExists) {
	  exodusFilePtr = ex_open(decoded_filename.c_str(), EX_WRITE,
				  &cpu_word_size, &io_word_size, &version);
	} else {
	  // If the first write for this file, create it...
	  io_word_size = cpu_word_size;
	  exodusFilePtr = ex_create(decoded_filename.c_str(), exodusMode,
				    &cpu_word_size, &io_word_size);
	  if (exodusFilePtr < 0) {
	    dbState = Ioss::STATE_INVALID;
	    // NOTE: Code will not continue past this call...
	    std::ostringstream errmsg;
	    errmsg << "Cannot create specified file '" << decoded_filename << "'";
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
	errmsg << "Problem opening specified file '" << decoded_filename << "'";
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
    }
    assert(exodusFilePtr >= 0);
    fileExists = true;
    return exodusFilePtr;
  }

  int DatabaseIO::free_file_pointer() const
  {
    if (exodusFilePtr != -1)
      ex_close(exodusFilePtr);
    exodusFilePtr = -1;

    return exodusFilePtr;
  }

  void DatabaseIO::put_qa()
  {
    static char qa_temp[4][MAX_STR_LENGTH+1];
    static char *qa[1][4] =
      {{qa_temp[0],qa_temp[1],qa_temp[2],qa_temp[3]}};

    Ioss::Utils::time_and_date(qa[0][3], qa[0][2], MAX_STR_LENGTH);

    std::string codename = "unknown";
    std::string version  = "unknown";

    if (get_region()->property_exists("code_name")) {
      codename = get_region()->get_property("code_name").get_string();
    }
    if (get_region()->property_exists("code_version")) {
      version = get_region()->get_property("code_version").get_string();
    }

    std::strncpy(qa[0][0], codename.c_str(), MAX_STR_LENGTH);
    std::strncpy(qa[0][1], version.c_str(),  MAX_STR_LENGTH);
    qa[0][0][MAX_STR_LENGTH] = '\0';
    qa[0][1][MAX_STR_LENGTH] = '\0';

    int ierr = ex_put_qa(get_file_pointer(), 1, qa);
    if (ierr < 0)
      exodus_error(get_file_pointer(), __LINE__, myProcessor);
  }

  void DatabaseIO::put_info()
  {
    // dump info records, include the product_registry
    assert(myProcessor == 0);

    // See if the input file was specified as a property on the database...
    std::string filename;
    if (get_region()->property_exists("input_file_name")) {
      filename = get_region()->get_property("input_file_name").get_string();
    }
    
    // Determine size of input file so can embed it in info records...
    std::vector<std::string> input_lines;
    Ioss::Utils::input_file(filename, &input_lines, max_line_length);

    size_t in_lines = input_lines.size();
    size_t qa_lines = 2; // Platform info and Version info...

    size_t total_lines = in_lines + qa_lines;

    char** info = get_exodus_names(total_lines, max_line_length); // 'total_lines' pointers to char buffers

    int i = 0;
    std::strncpy(info[i++], Ioss::Utils::platform_information().c_str(),
		 max_line_length);

    std::strncpy(info[i++], Version(), max_line_length);

    // Copy input file lines into 'info' array...
    for (size_t j=0; j < input_lines.size(); j++, i++) {
      std::strncpy(info[i], input_lines[j].c_str(), max_line_length);
      info[i][max_line_length] = '\0'; // Once more for good luck...
    }

    int ierr = ex_put_info(get_file_pointer(), total_lines, info);
    if (ierr < 0)
      exodus_error(get_file_pointer(), __LINE__, myProcessor);

    delete_exodus_names(info, total_lines);
  }

  void DatabaseIO::read_meta_data()
  {
    {
      Ioss::SerializeIO	serializeIO__(this);

      get_file_pointer();

      Ioex::Internals data(get_file_pointer(), maximumNameLength);
      data.check_processor_info(util().parallel_size(), myProcessor);

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
      exodus_error(get_file_pointer(), __LINE__, myProcessor);

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
      errmsg << "Negative node count was found in the model\n"
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
      errmsg << "Negative element count was found in the model, file: '"
	     << decoded_filename << "'";
      IOSS_ERROR(errmsg);
    }

    if (elementCount > 0 && m_groupCount[EX_ELEM_BLOCK] <= 0) {
      // NOTE: Code will not continue past this call...
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      std::ostringstream errmsg;
      errmsg << "No element blocks were found in the model, file: '" << decoded_filename << "'";
      IOSS_ERROR(errmsg);
    }

    Ioss::Region *this_region = get_region();
    this_region->property_add(Ioss::Property(std::string("title"), info.title));
    this_region->property_add(Ioss::Property(std::string("spatial_dimension"),
					     spatialDimension));
  }

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

  void DatabaseIO::get_step_times()
  {
    bool exists = false;
    double last_time = DBL_MAX;
    int timestep_count = 0;
    std::vector<double> tsteps(0);

    {
      Ioss::SerializeIO	serializeIO__(this);
      timestep_count = ex_inquire_int(get_file_pointer(), EX_INQ_TIME);
      if (timestep_count <= 0)
	return;
      
      // For an exodusII file, timesteps are global and are stored in the region.
      // Read the timesteps and add to the region
      tsteps.resize(timestep_count);
      int error = ex_get_all_times(get_file_pointer(), TOPTR(tsteps));
      if (error < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // See if the "last_written_time" attribute exists and if it
      // does, check that it matches the largest time in 'tsteps'.
      Ioex::Internals data(get_file_pointer(), maximumNameLength);
      exists = data.read_last_time_attribute(&last_time);
    }

    if (exists) {
      // Assume that if it exists on 1 processor, it exists on
      // all... Sync value among processors since could have a
      // corrupt step on only a single database.
      last_time = util().global_minmax(last_time, Ioss::ParallelUtils::DO_MIN);

      // Increase value slightly in case there is some small
      // roundoff and to allow for floating point equality issues.
      last_time *= 1.0001;
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
	this_region->add_state(tsteps[i]);
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
    int global_nodes       = nodeCount;
    int global_elements    = elementCount;
    int global_eblocks     = 0; // unused
    int global_nsets       = 0; // unused
    int global_ssets       = 0; // unused
    
    int num_external_nodes; // unused
    int num_elem_cmaps     = 0;
    int num_node_cmaps     = 0;
    int num_internal_nodes = nodeCount;
    int num_border_nodes   = 0;
    int num_internal_elems = elementCount;
    int num_border_elems   = 0;
    
    if (isParallel) {
      int error = ne_get_init_info(get_file_pointer(),
				   &num_proc, &num_proc_in_file, &file_type[0]);
      if (error < 0) {
	// Not a nemesis file
	if (util().parallel_size() > 1) {
	  std::ostringstream errmsg;
	  errmsg << "Exodus file does not contain nemesis information.\n";
	  IOSS_ERROR(errmsg);
	}
      } else if (num_proc != util().parallel_size() && util().parallel_size() > 1) {
	std::ostringstream errmsg;
	errmsg <<  "Exodus file was decomposed for " << num_proc
	       << " processors; application is currently being run on "
	       << util().parallel_size() << " processors";
	IOSS_ERROR(errmsg);

      } else if (num_proc_in_file != 1) {
	std::ostringstream errmsg;
	errmsg <<"Exodus file contains data for " << num_proc_in_file
	       << " processors; application requires 1 processor per file.";
	IOSS_ERROR(errmsg);

      } else if (file_type[0] != 'p') {
	std::ostringstream errmsg;
	errmsg << "Exodus file contains scalar nemesis data; application requires parallel nemesis data.";
	IOSS_ERROR(errmsg);
      }

      error = ne_get_loadbal_param(get_file_pointer(),
				   &num_internal_nodes,
				   &num_border_nodes,
				   &num_external_nodes,
				   &num_internal_elems,
				   &num_border_elems,
				   &num_node_cmaps,
				   &num_elem_cmaps,
				   myProcessor);
      if (error < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      commsetNodeCount = num_node_cmaps;
      commsetElemCount = num_elem_cmaps;

      // A nemesis file typically separates nodes into multiple
      // communication sets by processor.  (each set specifies
      // nodes/elements that communicate with only a single processor).
      // For Sierra, we want a single node commun. map and a single
      // element commun. map specifying all communications so we combine
      // all sets into a single set.
      error = ne_get_init_global(get_file_pointer(), &global_nodes, &global_elements,
				 &global_eblocks, &global_nsets, &global_ssets);
      if (error < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);
    }

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

    // Possibly, the following 4 fields should be nodesets and element
    // sets instead of fields on the region...
    region->field_add(Ioss::Field("internal_nodes", Ioss::Field::INTEGER,
				  SCALAR(),
				  Ioss::Field::COMMUNICATION,
				  num_internal_nodes));
    region->field_add(Ioss::Field("border_nodes", Ioss::Field::INTEGER,
				  SCALAR(),
				  Ioss::Field::COMMUNICATION,
				  num_border_nodes));
    region->field_add(Ioss::Field("internal_elements", Ioss::Field::INTEGER,
				  SCALAR(),
				  Ioss::Field::COMMUNICATION,
				  num_internal_elems));
    region->field_add(Ioss::Field("border_elements", Ioss::Field::INTEGER,
				  SCALAR(),
				  Ioss::Field::COMMUNICATION,
				  num_border_elems));

    assert(nodeCount    == num_internal_nodes + num_border_nodes);
    assert(elementCount == num_internal_elems + num_border_elems);
  }

  const Ioss::MapContainer& DatabaseIO::get_map(ex_entity_type type) const
  {
    if (type == EX_NODE_BLOCK || type == EX_NODE_SET) {
      return get_node_map();
    } else if (type == EX_ELEM_BLOCK || type == EX_ELEM_SET) {
      return get_element_map();
    } else if (type == EX_FACE_BLOCK || type == EX_FACE_SET) {
      return get_face_map();
    } else if (type == EX_EDGE_BLOCK || type == EX_EDGE_SET) {
      return get_edge_map();
    }
    std::ostringstream errmsg;
    errmsg << "INTERNAL ERROR: Invalid map type. "
           << "Something is wrong in the Ioex::DatabaseIO::get_map() function. "
           << "Please report.\n";
    IOSS_ERROR(errmsg);
  }

  const Ioss::MapContainer& DatabaseIO::get_node_map() const
  {
    // Allocate space for node number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (nodeMap.empty()) {
      nodeMap.resize(nodeCount+1);

      if (is_input()) {
	bool backward_compat = get_node_global_id_backward_compatibility() &&
	  !isParallel && dbUsage == Ioss::READ_MODEL;
	if (!backward_compat) {
	  Ioss::SerializeIO	serializeIO__(this);
	  // Check whether there is a "original_global_id_map" map on
	  // the database. If so, use it instead of the "node_num_map".
	  bool map_read = false;
	  int map_count = ex_inquire_int(get_file_pointer(), EX_INQ_NODE_MAP);
	  if (map_count > 0 && !get_node_global_id_backward_compatibility()) {
	    char **names = get_exodus_names(map_count, maximumNameLength);
	    int ierr = ex_get_names(get_file_pointer(), EX_NODE_MAP, names);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, -1);

	    if (map_count == 1 && Ioss::Utils::case_strcmp(names[0], "original_global_id_map") == 0) {
	      int error = ex_get_num_map(get_file_pointer(), EX_NODE_MAP, 1, &nodeMap[1]);
	      if (error >= 0) {
		map_read = true;
	      } else {
		// Clear out the vector...
		Ioss::MapContainer().swap(nodeMap);
		exodus_error(get_file_pointer(), __LINE__, myProcessor);
		map_read = false;
	      }
	    }
	    delete_exodus_names(names, map_count);
	  }

	  if (!map_read) {
	    int error = ex_get_node_num_map(get_file_pointer(), &nodeMap[1]);
	    if (error < 0) {
	      // Clear out the vector...
	      Ioss::MapContainer().swap(nodeMap);
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	    }
	  }
	} else {
	  // The node map is ignored in the serial case (Unless
	  // reading a restart file). This is to provide a consistent
	  // id space for a serial model and a model that has been
	  // decomposed using nem_slice/nem_spread.
	  for (int i=1; i < nodeCount+1; i++) {
	    nodeMap[i] = i;
	  }
	}

	// Check for sequential node map.
	// If not, build the reverse G2L node map...
	nodeMap[0] = -1;
	for (int i=1; i < nodeCount+1; i++) {
	  if (i != nodeMap[i]) {
	    nodeMap[0] = 1;
	    break;
	  }
	}

	if (nodeMap[0] == 1) {
	  Ioss::Map::build_reverse_map(&reverseNodeMap, &nodeMap[1], nodeCount, 0,
				       "NodeBlock", myProcessor);
	}

      } else {
	// Output database; nodeMap not set yet... Build a default map.
	for (int i=1; i < nodeCount+1; i++) {
	  nodeMap[i] = i;
	}
	// Sequential map
	nodeMap[0] = -1;
      }
    }
    return nodeMap;
  }

  const Ioss::MapContainer& DatabaseIO::get_element_map() const
  {
    // Allocate space for elemente number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (elementMap.empty()) {
      elementMap.resize(elementCount+1);

      if (is_input()) {
	bool backward_compat = get_node_global_id_backward_compatibility() &&
	  !isParallel && dbUsage == Ioss::READ_MODEL;
	if (!backward_compat) {
	  Ioss::SerializeIO	serializeIO__(this);
	  
	  // Check whether there is a "original_global_id_map" map on
	  // the database. If so, use it instead of the "elem_num_map".
	  bool map_read = false;
	  int map_count = ex_inquire_int(get_file_pointer(), EX_INQ_ELEM_MAP);
	  if (map_count > 0) {
	    int max_name_length = ex_inquire_int(get_file_pointer(), EX_INQ_DB_MAX_USED_NAME_LENGTH);
	    char **names = get_exodus_names(map_count, max_name_length);
	    int ierr = ex_get_names(get_file_pointer(), EX_ELEM_MAP, names);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, -1);

	    if (map_count == 1 && Ioss::Utils::case_strcmp(names[0], "original_global_id_map") == 0) {
	      int error = ex_get_num_map(get_file_pointer(), EX_ELEM_MAP, 1, &elementMap[1]);
	      if (error >= 0) {
		map_read = true;
	      } else {
		// Clear out the vector...
		Ioss::MapContainer().swap(elementMap);
		exodus_error(get_file_pointer(), __LINE__, myProcessor);
		map_read = false;
	      }
	    }
	    delete_exodus_names(names, map_count);
	  }

	  if (!map_read) {
	    int error = ex_get_id_map(get_file_pointer(), EX_ELEM_MAP, &elementMap[1]);
	    if (error < 0) {
	      // Clear out the vector...
	      Ioss::MapContainer().swap(elementMap);
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	    }
	  }
	} else {
	  // The element map is ignored in the serial case. This is to
	  // provide a consistent id space for a serial model and a
	  // model that has been decomposed using nem_slice/nem_spread.
	  for (int i=1; i < elementCount+1; i++) {
	    elementMap[i] = i;
	  }
	}

	// Check for sequential element map.
	// If not, build the reverse G2L element map...
	elementMap[0] = -1;
	for (int i=1; i < elementCount+1; i++) {
	  if (i != elementMap[i]) {
	    elementMap[0] = 1;
	    break;
	  }
	}

	if (elementMap[0] == 1) {
	  Ioss::Map::build_reverse_map(&reverseElementMap,
				       &elementMap[1], elementCount, 0,
				       "element", myProcessor);
	}

      } else {
	// Output database; elementMap not set yet... Build a default map.
	for (int i=1; i < elementCount+1; i++) {
	  elementMap[i] = i;
	}
	// Sequential map
	elementMap[0] = -1;
      }
    }
    return elementMap;
  }

  const Ioss::MapContainer& DatabaseIO::get_face_map() const
  {
    // Allocate space for face number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (faceMap.empty()) {
      faceMap.resize(faceCount+1);

      if (is_input()) {
	Ioss::SerializeIO	serializeIO__(this);
	
	int error = ex_get_id_map(get_file_pointer(), EX_FACE_MAP, &faceMap[1]);
	if (error < 0) {
	  // Clear out the vector...
	  Ioss::MapContainer().swap(faceMap);
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}
	
	// Check for sequential face map.
	// If not, build the reverse G2L face map...
	faceMap[0] = -1;
	for (int i=1; i < faceCount+1; i++) {
	  if (i != faceMap[i]) {
	    faceMap[0] = 1;
	    break;
	  }
	}
	
	if (faceMap[0] == 1) {
	  Ioss::Map::build_reverse_map(&reverseFaceMap,
				       &faceMap[1], faceCount, 0,
				       "face", myProcessor);
	}
	
      } else {
	// Output database; faceMap not set yet... Build a default map.
	for (int i=1; i < faceCount+1; i++) {
	  faceMap[i] = i;
	}
	// Sequential map
	faceMap[0] = -1;
      }
    }
    return faceMap;
  }

  const Ioss::MapContainer& DatabaseIO::get_edge_map() const
  {
    // Allocate space for edge number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (edgeMap.empty()) {
      edgeMap.resize(edgeCount+1);

      if (is_input()) {
	Ioss::SerializeIO	serializeIO__(this);
	
	int error = ex_get_id_map(get_file_pointer(), EX_EDGE_MAP, &edgeMap[1]);
	if (error < 0) {
	  // Clear out the vector...
	  Ioss::MapContainer().swap(edgeMap);
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}
	
	// Check for sequential edge map.
	// If not, build the reverse G2L edge map...
	edgeMap[0] = -1;
	for (int i=1; i < edgeCount+1; i++) {
	  if (i != edgeMap[i]) {
	    edgeMap[0] = 1;
	    break;
	  }
	}
	
	if (edgeMap[0] == 1) {
	  Ioss::Map::build_reverse_map(&reverseEdgeMap,
				       &edgeMap[1], edgeCount, 0,
				       "edge", myProcessor);
	}
	
      } else {
	// Output database; edgeMap not set yet... Build a default map.
	for (int i=1; i < edgeCount+1; i++) {
	  edgeMap[i] = i;
	}
	// Sequential map
	edgeMap[0] = -1;
      }
    }
    return edgeMap;
  }

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
    int ierr = ex_get_attr_param(get_file_pointer(), EX_NODE_BLOCK, 1, &num_attr);
    if (ierr < 0)
      exodus_error(get_file_pointer(), __LINE__, myProcessor);

    add_attribute_fields(EX_NODE_BLOCK, block, num_attr, "");
    add_results_fields(EX_NODE_BLOCK, block);

    get_region()->add(block);
  }

  void DatabaseIO::get_elemblocks()
  {
    get_blocks(EX_ELEM_BLOCK, 0);
  }
  
  void DatabaseIO::get_faceblocks()
  {
    get_blocks(EX_FACE_BLOCK, 1);
  }
  
  void DatabaseIO::get_edgeblocks()
  {
    get_blocks(EX_EDGE_BLOCK, 2);
  }
  
  void DatabaseIO::get_blocks(ex_entity_type entity_type, int rank_offset)
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

    Ioss::IntVector X_block_ids(m_groupCount[entity_type]);

    int error;
    {
      Ioss::SerializeIO	serializeIO__(this);

      error = ex_get_ids(get_file_pointer(), entity_type, TOPTR(X_block_ids));
      if (error < 0) {
	exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
    }

    size_t all_X_type_length = m_groupCount[entity_type] * (MAX_STR_LENGTH+1);
    std::vector<char> all_X_type(all_X_type_length);

    Ioss::IntVector counts(m_groupCount[entity_type] * 4);
    Ioss::IntVector local_X_count(m_groupCount[entity_type]);
    Ioss::IntVector global_X_count(m_groupCount[entity_type]);
    int iblk ;

    {
      Ioss::SerializeIO	serializeIO__(this);

      for ( iblk = 0; iblk < m_groupCount[entity_type]; iblk++) {
	int index = 4*iblk;
	int id = X_block_ids[iblk];
	int number_Xs   = 0;
	int nodes_per_X = 0;
	int edges_per_X = 0;
	int faces_per_X = 0;
	int attributes_per_X = 0;

	char * const X_type = TOPTR(all_X_type) + iblk * (MAX_STR_LENGTH+1);

	error = ex_get_block(get_file_pointer(), entity_type, id, X_type,
			     &number_Xs, &nodes_per_X, &edges_per_X, &faces_per_X,
			     &attributes_per_X);
	if (error < 0) {
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}

	local_X_count[iblk] = number_Xs;

	counts[index+0] = nodes_per_X;
	counts[index+1] = edges_per_X;
	counts[index+2] = faces_per_X;
	counts[index+3] = attributes_per_X;

	if (number_Xs == 0) {
	  counts[index+0] = 0;
	  counts[index+3] = 0;
	  std::memset(X_type, 0, MAX_STR_LENGTH + 1);
	} else {
	  counts[index+0] = nodes_per_X;
	  counts[index+3] = attributes_per_X;
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
    int offset = 0;
    int used_blocks = 0;
    
    for ( iblk = 0; iblk < m_groupCount[entity_type]; iblk++) {
      int index = 4*iblk;
      int nodes_per_X = counts[index+0];
      int edges_per_X = counts[index+1];
      int faces_per_X = counts[index+2];
      int attributes        = counts[index+3];

      int id = X_block_ids[iblk];
      std::string alias = Ioss::Utils::encode_entity_name("block", id);
      char * const X_type = TOPTR(all_X_type) + iblk * (MAX_STR_LENGTH+1);

      std::string block_name;
      {
	Ioss::SerializeIO	serializeIO__(this);
	block_name = get_entity_name(get_file_pointer(), entity_type, id, "block", maximumNameLength);
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
	get_region()->add(fblock);
	block = fblock;
      } else if (entity_type == EX_EDGE_BLOCK) {
	Ioss::EdgeBlock *eblock = new Ioss::EdgeBlock(this, block_name, type, local_X_count[iblk]);
	get_region()->add(eblock);
	block = eblock;
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
				     Ioss::Field::INTEGER, storage, Ioss::Field::MESH,
				     local_X_count[iblk]));
      }
      // X -> Edge?
      if (edges_per_X > 0 && rank_offset < 2) {
	std::string storage = "Real["+Ioss::Utils::to_string(edges_per_X)+"]";
	block->field_add(Ioss::Field("connectivity_edge",
				     Ioss::Field::INTEGER, storage, Ioss::Field::MESH,
				     local_X_count[iblk]));
      }

      block->property_add(Ioss::Property("id", id));
      
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
      if (block_name != alias) {
	get_region()->add_alias(block_name, alias);
      }
    
      // Check for additional variables.
      add_attribute_fields(entity_type, block, attributes, type);
      add_results_fields(entity_type, block, iblk);
      
      if (entity_type == EX_ELEM_BLOCK) {
	Ioss::SerializeIO	serializeIO__(this);
	add_map_fields(get_file_pointer(), (Ioss::ElementBlock*)block, local_X_count[iblk]);
      }
    }
    m_groupCount[entity_type] = used_blocks;
  }

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

    std::vector<int> node_used(nodeCount);
    std::vector<std::vector<int> > inv_con(nodeCount);
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    assert(check_block_order(element_blocks));

    {
      Ioss::SerializeIO	serializeIO__(this);
      for (int iblk = 0; iblk < m_groupCount[EX_ELEM_BLOCK]; iblk++) {
	Ioss::ElementBlock *eb = element_blocks[iblk];
	int blk_position =  eb->get_property("original_block_order").get_int();
	int id =            eb->get_property("id").get_int();
	int element_nodes = eb->get_property("topology_node_count").get_int();
	int my_element_count = eb->get_property("entity_count").get_int();
	if (my_element_count > 0) {
	  std::vector<int> conn(my_element_count * element_nodes);
	  ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, TOPTR(conn), NULL, NULL);
	  
	  for (int i=0; i < my_element_count * element_nodes; i++) {
	    node_used[conn[i]-1] = blk_position+1;
	  }

	  for (int i=0; i < nodeCount; i++) {
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
      check_non_null(css, "communication map", "commset_node");
      std::vector<std::pair<int,int> > proc_node;
      {
	std::vector<int> entity_processor;
	css->get_field_data("entity_processor", entity_processor);
	proc_node.resize(entity_processor.size()/2);
	size_t j=0;
	for (size_t i=0; i < entity_processor.size(); j++, i+=2) {
	  proc_node[j].second = entity_processor[i];
	  proc_node[j].first  = entity_processor[i+1];
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
	int glob_id = proc_node[i].second;
	int proc    = proc_node[i].first;
	procs[proc]++;
	send[offset++] = glob_id;
	int loc_id = node_global_to_local(glob_id, true) -1;
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
	       << " in Ioex::DatabaseIO::compute_block_adjacencies";
	std::cerr << errmsg.str();
      }

      int local_error = (MPI_SUCCESS == result) ? 0 : 1 ;
      int global_error = util().global_minmax(local_error, Ioss::ParallelUtils::DO_MAX);

      if (global_error != 0) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Irecv error on some processor "
	       << "in Ioex::DatabaseIO::compute_block_adjacencies";
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
	       << " in Ioex::DatabaseIO::compute_block_adjacencies";
	std::cerr << errmsg.str();
      }

      local_error = (MPI_SUCCESS == result) ? 0 : 1 ;
      global_error = util().global_minmax(local_error, Ioss::ParallelUtils::DO_MAX);

      if (global_error != 0) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Rsend error on some processor "
	       << "in Ioex::DatabaseIO::compute_block_adjacencies";
	IOSS_ERROR(errmsg);
      }
      
      result = MPI_Waitall(req_cnt, TOPTR(request), TOPTR(status));
      
      if (result != MPI_SUCCESS) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Waitall error on processor " << util().parallel_rank()
	       << " in Ioex::DatabaseIO::compute_block_adjacencies";
	std::cerr << errmsg.str();
      }

      // Unpack the data and update the inv_con arrays for boundary
      // nodes...
      offset = 0;
      for (size_t i=0; i < proc_node.size(); i++) {
	int glob_id = recv[offset++];
	int loc_id = node_global_to_local(glob_id, true) -1;
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

    for (int i = 0; i < nodeCount; i++) {
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
      int offset = 0;
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
      if (block_is_omitted(block)) {
	status = 1;
      }

      int id               = block->get_property("id").get_int();
      int element_nodes    = block->get_property("topology_node_count").get_int();
      int my_element_count = block->get_property("entity_count").get_int();
      if (my_element_count > 0) {
	std::vector<int> conn(my_element_count * element_nodes);
	ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, TOPTR(conn), NULL, NULL);
	
	for (int j=0; j < my_element_count * element_nodes; j++) {
	  nodeConnectivityStatus[conn[j]-1] |= status;
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
      SideSetMap fs_map;
      SideSetSet fs_set;

      Ioss::IntVector side_set_ids(m_groupCount[EX_SIDE_SET]);
      {
	Ioss::SerializeIO	serializeIO__(this);
	int error = ex_get_ids(get_file_pointer(), EX_SIDE_SET, TOPTR(side_set_ids));
	if (error < 0) {
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}

	for (int i = 0; i < m_groupCount[EX_SIDE_SET]; i++) {
          std::vector<char> ss_name(maximumNameLength+1);
	  error = ex_get_name(get_file_pointer(), EX_SIDE_SET, side_set_ids[i],
                              TOPTR(ss_name));
	  if (error < 0) {
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }
	  if (ss_name[0] != '\0') {
	    Ioss::Utils::fixup_name(TOPTR(ss_name));
	    decode_surface_name(fs_map, fs_set, TOPTR(ss_name));
	  }
	}
      }

      // Create sidesets for each entry in the fs_set... These are the
      // sidesets which were probably written by a previous run of the
      // IO system and are already split into homogenous pieces...
      {
	SideSetSet::iterator I = fs_set.begin();
	while (I != fs_set.end()) {
	  const std::string fs_name = *I;
	  Ioss::SideSet *side_set = new Ioss::SideSet(this, fs_name);
	  get_region()->add(side_set);
	  int id = extract_id(fs_name);
	  if (id > 0)
	    side_set->property_add(Ioss::Property("id", id));
	  ++I;
	}
      }

      for (int iss = 0; iss < m_groupCount[EX_SIDE_SET]; iss++) {
	int id = side_set_ids[iss];
	std::string sid = "";
	int number_sides;
	int number_distribution_factors;
	TopologyMap topo_map;
	TopologyMap side_map; // Used to determine side consistency

	Ioss::SurfaceSplitType split_type = splitType;
	std::string side_set_name;
	Ioss::SideSet *side_set = NULL;

	{
	  Ioss::SerializeIO	serializeIO__(this);

	  side_set_name = get_entity_name(get_file_pointer(), EX_SIDE_SET, id, "surface",
					  maximumNameLength);

	  if (side_set_name == "universal_sideset") {
	    split_type = Ioss::SPLIT_BY_DONT_SPLIT;
	  }

	  bool in_fs_map = false;
	  SideSetMap::iterator FSM = fs_map.find(side_set_name);
	  if (FSM != fs_map.end()) {
	    in_fs_map = true;
	    std::string efs_name = (*FSM).second;
	    side_set = get_region()->get_sideset(efs_name);
	    check_non_null(side_set, "sideset", efs_name);
	  } else {
	    side_set = new Ioss::SideSet(this, side_set_name);
	    get_region()->add((Ioss::SideSet*)side_set);
	    side_set->property_add(Ioss::Property("id", id));

	    get_region()->add_alias(side_set_name, Ioss::Utils::encode_entity_name("surface", id));
	    get_region()->add_alias(side_set_name, Ioss::Utils::encode_entity_name("sideset", id));
	  }

	  //	  split_type = SPLIT_BY_ELEMENT_BLOCK;
	  //	  split_type = SPLIT_BY_TOPOLOGIES;
	  //	  split_type = SPLIT_BY_DONT_SPLIT;

	  // Determine how many side blocks compose this side set.
	  int error = ex_get_set_param(get_file_pointer(), EX_SIDE_SET, id,
				       &number_sides, &number_distribution_factors);
	  if (error < 0) {
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }

	  // Get the element list and filter out any sides connected to an omitted element block...
	  Ioss::IntVector element(number_sides);
	  Ioss::IntVector sides(number_sides);
	  
	  int ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(sides));
	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  
	  if (!blockOmissions.empty()) {
	    // Potentially changes element and sides...
	    filter_element_list(get_region(), element, sides, true);
	    number_sides = (int)element.size();
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

	    separate_surface_element_sides(element, sides, get_region(), topo_map, side_map, split_type);
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
	      if (!block_is_omitted(block)) {
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
		    int my_side_count = block->topology()->number_boundaries();
		    for (int ii = 0; ii < my_side_count; ii++) {
		      const Ioss::ElementTopology *topo = block->topology()->boundary_type(ii+1);
		      topo_map[std::make_pair(name, topo)] = 0;
		      side_map[std::make_pair(name, topo)] = 0;
		    }
		  }
		}
	      }
	    }
	    separate_surface_element_sides(element, sides, get_region(), topo_map, side_map, split_type);
	  }
	}

	// End of first step in splitting.  Check among all processors
	// to see which potential splits have sides in them...
	Ioss::IntVector global_side_counts(topo_map.size());
	{
	  int i = 0;
	  {
	    TopologyMap::const_iterator I = topo_map.begin();
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
	      TopologyMap::const_iterator I = side_map.begin();
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

	int i = 0;
	TopologyMap::const_iterator I = topo_map.begin();
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
		int my_side_count = (*I).second;

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
		  if (block == NULL || block_is_omitted(block)) {
		    std::ostringstream errmsg;
		    errmsg << "INTERNAL ERROR: Could not find element block '" << topo_or_block_name 
			   << "' Something is wrong in the Ioex::DatabaseIO class. Please report.\n";
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
		if (side_map.size() > 0) {
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
		    && (number_distribution_factors > 0 || isParallel)
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
						    Ioss::Field::INTEGER, "scalar",
						    Ioss::Field::MESH, my_side_count));
		}

		int num_attr = 0;
		int ierr = ex_get_attr_param(get_file_pointer(), EX_SIDE_SET, 1, &num_attr);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      
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

    void DatabaseIO::compute_block_membership(int id, std::vector<std::string> &block_membership) const
    {
      Ioss::IntVector block_ids(m_groupCount[EX_ELEM_BLOCK]);
      if (m_groupCount[EX_ELEM_BLOCK] == 1) {
	block_ids[0] = 1;
      } else {
	{
	  Ioss::SerializeIO	serializeIO__(this);
	
	  int number_sides;
	  int number_distribution_factors;
	  int error = ex_get_set_param(get_file_pointer(), EX_SIDE_SET, id,
				       &number_sides, &number_distribution_factors);
	  if (error < 0) {
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }
      
	  if (number_sides > 0) {
	    // Get the element and element side lists.
	    Ioss::IntVector element(number_sides);
	    Ioss::IntVector sides(number_sides);
	
	    int ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(sides));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	
	    Ioss::ElementBlock *block = NULL;
	    for (int iel = 0; iel < number_sides; iel++) {
	      int elem_id = element[iel];
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

      for (int i=0; i < m_groupCount[EX_ELEM_BLOCK]; i++) {
	if (block_ids[i] == 1) {
	  Ioss::ElementBlock *block = element_blocks[i];
	  if (!block_is_omitted(block)) {
	    block_membership.push_back(block->name());
	  }
	}
      }
    }
  
    void DatabaseIO::compute_block_membership(Ioss::SideBlock *efblock,
					      std::vector<std::string> &block_membership) const
    {
      Ioss::IntVector block_ids(m_groupCount[EX_ELEM_BLOCK]);
      if (m_groupCount[EX_ELEM_BLOCK] == 1) {
	block_ids[0] = 1;
      } else {
	Ioss::IntVector element_side;
	efblock->get_field_data("element_side", element_side);
      
	size_t number_sides = element_side.size() / 2;
	Ioss::ElementBlock *block = NULL;
	for (size_t iel = 0; iel < number_sides; iel++) {
	  int elem_id = element_side[2*iel];  // Vector contains both element and side.
	  elem_id = element_global_to_local(elem_id);
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
	  if (!block_is_omitted(block)) {
	    block_membership.push_back(block->name());
	  }
	}
      }
    }
  
    template <typename T>
      void DatabaseIO::get_sets(ex_entity_type type, int count, const std::string &base, const T*)
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
	  Ioss::IntVector Xset_ids(count);
	  Ioss::IntVector attributes(count);
	  std::vector<T*> Xsets(count);

	  {
	    Ioss::SerializeIO	serializeIO__(this);
	    int error = ex_get_ids(get_file_pointer(), type, TOPTR(Xset_ids));
	    if (error < 0) {
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	    }

	    for (int ins = 0; ins < count; ins++) {
	      int id = Xset_ids[ins];
	      int number_Xs;
	      int number_distribution_factors;
	  
	      error = ex_get_set_param(get_file_pointer(), type, id,
				       &number_Xs, &number_distribution_factors);
	      if (error < 0) {
		exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }

	      int num_attr = 0;
	      int ierr = ex_get_attr_param(get_file_pointer(), type, id,
					   &num_attr);
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      attributes[ins] = num_attr;

	      std::string Xset_name = get_entity_name(get_file_pointer(), type, id, base+"list",
						      maximumNameLength);

	      T* Xset = new T(this, Xset_name, number_Xs);
	      Xsets[ins] = Xset;
	      Xset->property_add(Ioss::Property("id", id));
	      get_region()->add(Xset);
	      get_region()->add_alias(Xset_name, Ioss::Utils::encode_entity_name(base+"list", id));
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

      if (isParallel) {
	Ioss::SerializeIO	serializeIO__(this);
	// This is a parallel run. There should be communications data
	// Get nemesis commset metadata
	int my_node_count = 0;
	int elem_count = 0;

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

	  int error = ne_get_cmap_params(get_file_pointer(),
					 TOPTR(nodeCmapIds), TOPTR(nodeCmapNodeCnts),
					 TOPTR(elemCmapIds), TOPTR(elemCmapElemCnts),
					 myProcessor);
	  if (error < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);

	  // Count nodes, elements, and convert counts to offsets.
	  //
	  for (int ics = 0; ics < commsetNodeCount; ics++) {
	    my_node_count += nodeCmapNodeCnts[ics];
	  }
	  for (int ecs = 0; ecs < commsetElemCount; ecs++) {
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

    int DatabaseIO::get_field_internal(const Ioss::Region* /* region */,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      // For now, assume that all TRANSIENT fields on a region
      // are REDUCTION fields (1 value).  We need to gather these
      // and output them all at one time.  The storage location is a
      // 'globalVariables' array
      {
	size_t num_to_get = field.verify(data_size);
	Ioss::SerializeIO	serializeIO__(this);

	Ioss::Field::RoleType role = field.get_role();

	if (role == Ioss::Field::TRANSIENT || role == Ioss::Field::REDUCTION) {
	  get_reduction_field(EX_GLOBAL, field, get_region(), data);
	} else {
	  std::ostringstream errmsg;
	  errmsg << "Can not handle non-TRANSIENT or non-REDUCTION fields on regions";
	  IOSS_ERROR(errmsg);
	}
	return num_to_get;
      }
    }

    int DatabaseIO::get_field_internal(const Ioss::NodeBlock* nb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);
	if (num_to_get > 0) {

#ifndef NDEBUG
	  int my_node_count = field.raw_count();
	  assert(my_node_count == nodeCount);
#endif

	  Ioss::Field::RoleType role = field.get_role();
	  if (role == Ioss::Field::MESH) {
	    if (field.get_name() == "mesh_model_coordinates") {
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
		exodus_error(get_file_pointer(), __LINE__, myProcessor);

	      int index = 0;
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
	      const Ioss::MapContainer &map = get_node_map();
	      int *ids = static_cast<int*>(data);

	      for (size_t i=0; i < num_to_get; i++) {
		ids[i] = map[1 + i];
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

    int DatabaseIO::get_field_internal(const Ioss::ElementBlock* eb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);
	if (num_to_get > 0) {

	  int ierr = 0;
	  int id = get_id(eb, EX_ELEM_BLOCK, &ids_);
	  int my_element_count = eb->get_property("entity_count").get_int();
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
		ierr = ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, static_cast<int*>(data), NULL, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);

		// Now, map the nodes in the connectivity from local to global ids
		const Ioss::MapContainer &map = get_node_map();
		if (!Ioss::Map::is_sequential(map)) {
		  int *connect = static_cast<int*>(data);
		  for (size_t i=0; i < num_to_get*element_nodes; i++)
		    connect[i] = map[connect[i]];
		}
	      }
	    }
	    else if (field.get_name() == "connectivity_face") {
	      int face_count = field.raw_storage()->component_count();

	      // The connectivity is stored in a 1D array.
	      // The element_face index varies fastest
	      if (my_element_count > 0) {
		ierr = ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, NULL, NULL, static_cast<int*>(data));
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);

		// Now, map the faces in the connectivity from local to global ids
		const Ioss::MapContainer &map = get_face_map();
		if (!Ioss::Map::is_sequential(map)) {
		  int *connect = static_cast<int*>(data);
		  for (size_t i=0; i < num_to_get*face_count; i++)
		    connect[i] = map[connect[i]];
		}
	      }
	    }
	    else if (field.get_name() == "connectivity_edge") {
	      int edge_count = field.raw_storage()->component_count();

	      // The connectivity is stored in a 1D array.
	      // The element_edge index varies fastest
	      if (my_element_count > 0) {
		ierr = ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, NULL, static_cast<int*>(data), NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);

		// Now, map the edges in the connectivity from local to global ids
		const Ioss::MapContainer &map = get_edge_map();
		if (!Ioss::Map::is_sequential(map)) {
		  int *connect = static_cast<int*>(data);
		  for (size_t i=0; i < num_to_get*edge_count; i++)
		    connect[i] = map[connect[i]];
		}
	      }
	    }
	    else if (field.get_name() == "connectivity_raw") {
	      // "connectivity_raw" has nodes in local id space (1-based)
	      assert(field.raw_storage()->component_count() == eb->get_property("topology_node_count").get_int());

	      // The connectivity is stored in a 1D array.
	      // The element_node index varies fastet
	      if (my_element_count > 0) {
		ierr = ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, static_cast<int*>(data), NULL, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }
	    }
	    else if (field.get_name() == "ids") {
	      // Map the local ids in this element block
	      // (eb_offset+1...eb_offset+1+my_element_count) to global element ids.
	      const Ioss::MapContainer &map = get_element_map();
	      int eb_offset = eb->get_offset();
	      int *ids = static_cast<int*>(data);

	      for (size_t i=0; i < num_to_get; i++) {
		ids[i] = map[eb_offset + 1 + i];
	      }
	    }
	    else if (field.get_name() == "skin") {
	      // This is (currently) for the skinned body. It maps the
	      // side element on the skin to the original element/local
	      // side number.  It is a two component field, the first
	      // component is the global id of the underlying element in
	      // the initial mesh and its local side number (1-based).

	      Ioss::IntVector element(my_element_count);
	      Ioss::IntVector side(my_element_count);
	      int *el_side = (int *)data;

	      // FIX: Hardwired map ids....
	      int eb_offset = eb->get_offset();
	      ex_get_partial_elem_map(get_file_pointer(), 1, eb_offset+1, my_element_count,
				      TOPTR(element));
	      ex_get_partial_elem_map(get_file_pointer(), 2, eb_offset+1, my_element_count,
				      TOPTR(side));

	      int index = 0;
	      for (int i=0; i < my_element_count; i++) {
		el_side[index++] = element[i];
		el_side[index++] = side[i];
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

    int DatabaseIO::get_field_internal(const Ioss::FaceBlock* eb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);
	if (num_to_get > 0) {

	  int ierr = 0;
	  int id = get_id(eb, EX_FACE_BLOCK, &ids_);
	  int my_face_count = eb->get_property("entity_count").get_int();
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
		ierr = ex_get_conn(get_file_pointer(), EX_FACE_BLOCK, id, static_cast<int*>(data), NULL, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);

		// Now, map the nodes in the connectivity from local to global ids
		const Ioss::MapContainer &map = get_node_map();
		if (!Ioss::Map::is_sequential(map)) {
		  int *connect = static_cast<int*>(data);
		  for (size_t i=0; i < num_to_get*face_nodes; i++)
		    connect[i] = map[connect[i]];
		}
	      }
	    }
	    else if (field.get_name() == "connectivity_edge") {
	      int edge_count = field.raw_storage()->component_count();

	      // The connectivity is stored in a 1D array.
	      // The face_edge index varies fastest
	      if (my_face_count > 0) {
		ierr = ex_get_conn(get_file_pointer(), EX_FACE_BLOCK, id, NULL, static_cast<int*>(data), NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);

		// Now, map the edges in the connectivity from local to global ids
		const Ioss::MapContainer &map = get_edge_map();
		if (!Ioss::Map::is_sequential(map)) {
		  int *connect = static_cast<int*>(data);
		  for (size_t i=0; i < num_to_get*edge_count; i++)
		    connect[i] = map[connect[i]];
		}
	      }
	    }
	    else if (field.get_name() == "connectivity_raw") {
	      // "connectivity_raw" has nodes in local id space (1-based)
	      assert(field.raw_storage()->component_count() == eb->get_property("topology_node_count").get_int());

	      // The connectivity is stored in a 1D array.
	      // The face_node index varies fastet
	      if (my_face_count > 0) {
		ierr = ex_get_conn(get_file_pointer(), EX_FACE_BLOCK, id, static_cast<int*>(data), NULL, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }
	    }
	    else if (field.get_name() == "ids") {
	      // Map the local ids in this face block
	      // (eb_offset+1...eb_offset+1+my_face_count) to global face ids.
	      const Ioss::MapContainer &map = get_map(EX_FACE_BLOCK);
	      int eb_offset = eb->get_offset();
	      int *ids = static_cast<int*>(data);

	      for (size_t i=0; i < num_to_get; i++) {
		ids[i] = map[eb_offset + 1 + i];
	      }
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

    int DatabaseIO::get_field_internal(const Ioss::EdgeBlock* eb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);
	if (num_to_get > 0) {

	  int ierr = 0;
	  int id = get_id(eb, EX_EDGE_BLOCK, &ids_);
	  int my_edge_count = eb->get_property("entity_count").get_int();
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
		ierr = ex_get_conn(get_file_pointer(), EX_EDGE_BLOCK, id, static_cast<int*>(data), NULL, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);

		// Now, map the nodes in the connectivity from local to global ids
		const Ioss::MapContainer &map = get_node_map();
		if (!Ioss::Map::is_sequential(map)) {
		  int *connect = static_cast<int*>(data);
		  for (size_t i=0; i < num_to_get*edge_nodes; i++)
		    connect[i] = map[connect[i]];
		}
	      }
	    }
	    else if (field.get_name() == "connectivity_raw") {
	      // "connectivity_raw" has nodes in local id space (1-based)
	      assert(field.raw_storage()->component_count() == eb->get_property("topology_node_count").get_int());

	      // The connectivity is stored in a 1D array.
	      // The edge_node index varies fastet
	      if (my_edge_count > 0) {
		ierr = ex_get_conn(get_file_pointer(), EX_EDGE_BLOCK, id, static_cast<int*>(data), NULL, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }
	    }
	    else if (field.get_name() == "ids") {
	      // Map the local ids in this edge block
	      // (eb_offset+1...eb_offset+1+my_edge_count) to global edge ids.
	      const Ioss::MapContainer &map = get_map(EX_EDGE_BLOCK);
	      int eb_offset = eb->get_offset();
	      int *ids = static_cast<int*>(data);

	      for (size_t i=0; i < num_to_get; i++) {
		ids[i] = map[eb_offset + 1 + i];
	      }
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

    int DatabaseIO::get_Xset_field_internal(ex_entity_type type,
					    const Ioss::EntitySet *ns,
					    const Ioss::Field& field,
					    void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);
	if (num_to_get > 0) {

	  int id = get_id(ns, type, &ids_);
	  Ioss::Field::RoleType role = field.get_role();
	  if (role == Ioss::Field::MESH) {

	    if (field.get_name() == "ids") {
	      int ierr = ex_get_set(get_file_pointer(), type, id, static_cast<int*>(data), NULL);
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);

	      // Convert the local node ids to global ids
	      const Ioss::MapContainer &map = get_map(type);
	      if (!Ioss::Map::is_sequential(map)) {
		int *ids = static_cast<int*>(data);
		for (size_t i=0; i < num_to_get; i++) {
		  ids[i] = map[ids[i]];
		}
	      }
            } else if (field.get_name() == "orientation") {
	      int ierr = ex_get_set(get_file_pointer(), type, id, NULL, static_cast<int*>(data));
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);

	    } else if (field.get_name() == "distribution_factors") {
	      int number_nodes;
	      int number_distribution_factors;

	      int ierr = ex_get_set_param(get_file_pointer(), type, id,
					  &number_nodes, &number_distribution_factors);
	      if (number_distribution_factors == 0) {
		double *rdata = static_cast<double*>(data);
		for (size_t i=0; i < num_to_get; i++)
		  rdata[i] = 1.0;
	      } else {
		ierr = ex_get_set_dist_fact(get_file_pointer(), type, id, static_cast<double*>(data));
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
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
					    
    int DatabaseIO::get_field_internal(const Ioss::NodeSet* ns,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      return get_Xset_field_internal(EX_NODE_SET, ns, field, data, data_size);
    }

    int DatabaseIO::get_field_internal(const Ioss::EdgeSet* ns,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      return get_Xset_field_internal(EX_EDGE_SET, ns, field, data, data_size);
    }

    int DatabaseIO::get_field_internal(const Ioss::FaceSet* ns,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      return get_Xset_field_internal(EX_FACE_SET, ns, field, data, data_size);
    }

    int DatabaseIO::get_field_internal(const Ioss::ElementSet* ns,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      return get_Xset_field_internal(EX_ELEM_SET, ns, field, data, data_size);
    }

    int DatabaseIO::get_field_internal(const Ioss::SideSet* fs,
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

    int DatabaseIO::get_field_internal(const Ioss::CommSet* cs,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);

	if (num_to_get > 0) {
	  int entity_count = cs->get_property("entity_count").get_int();

	  // Return the <entity (node or side), processor> pair
	  if (field.get_name() == "entity_processor") {

	    // Check type -- node or side
	    std::string type = cs->get_property("entity_type").get_string();

	    // Allocate temporary storage space
	    Ioss::IntVector entities(num_to_get);
	    Ioss::IntVector procs(num_to_get);

	    if (type == "node") {
	      int cm_offset = 0;
	      int i;
	      for (i=0; i < commsetNodeCount; i++) {
		int ierr = ne_get_node_cmap(get_file_pointer(), nodeCmapIds[i],
					    &entities[cm_offset], &procs[cm_offset],
					    myProcessor);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
		cm_offset += nodeCmapNodeCnts[i];
	      }
	      assert(cm_offset == entity_count);

	      // Convert local node id to global node id and store in 'data'
	      int* entity_proc = static_cast<int*>(data);
	      const Ioss::MapContainer &map = get_node_map();

	      int j=0;
	      for (i=0; i < entity_count; i++) {
		int local_id = entities[i];
		entity_proc[j++] = map[local_id];
		entity_proc[j++] = procs[i];
	      }
	    } else if (type == "side") {
	      Ioss::IntVector sides(entity_count);
	      int cm_offset = 0;
	      int i;
	      for (i=0; i < commsetElemCount; i++) {
		int ierr = ne_get_elem_cmap(get_file_pointer(), elemCmapIds[i],
					    &entities[cm_offset], &sides[cm_offset],
					    &procs[cm_offset], myProcessor);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
		cm_offset += elemCmapElemCnts[i];
	      }
	      assert(cm_offset == entity_count);

	      int* entity_proc = static_cast<int*>(data);
	      const Ioss::MapContainer &map = get_element_map();

	      int j=0;
	      for (i=0; i < entity_count; i++) {
		int local_id = entities[i];
		entity_proc[j++] = 10*map[local_id]+sides[i];
		entity_proc[j++] = procs[i];
	      }
	    } else {
	      std::ostringstream errmsg;
	      errmsg << "Invalid commset type " << type;
	      IOSS_ERROR(errmsg);
	    }

	  } else {
	    num_to_get = Ioss::Utils::field_warning(cs, field, "input");
	  }
	}
	return num_to_get;
      }
    }

    int DatabaseIO::get_field_internal(const Ioss::SideBlock* fb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      Ioss::SerializeIO	serializeIO__(this);
      size_t num_to_get = field.verify(data_size);
      if (num_to_get > 0) {

	int id = get_id(fb, EX_SIDE_SET, &ids_);
	size_t entity_count = fb->get_property("entity_count").get_int();
	if (num_to_get != entity_count) {
	  std::ostringstream errmsg;
	  errmsg << "Partial field input not yet implemented for side blocks";
	  IOSS_ERROR(errmsg);
	}

	int number_sides;
	int number_distribution_factors;
	int ierr = ex_get_set_param(get_file_pointer(), EX_SIDE_SET, id,
				    &number_sides, &number_distribution_factors);
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
      
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
	  if (field.get_name() == "side_ids" &&
	      fb->name() == "universal_sideset") {
	    // The side ids are being stored as the distribution factor
	    // field on the universal sideset.  There should be no other
	    // side sets that request this field...  (Eventually,
	    // create an id field to store this info.

	    if (number_distribution_factors == static_cast<int>(num_to_get)) {
	      std::vector<double> real_ids(num_to_get);
	      ierr = ex_get_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id, TOPTR(real_ids));
	    
	      // Need to convert 'double' to 'int' for Sierra use...
	      int* ids = static_cast<int*>(data);
	      for (size_t i = 0; i < num_to_get; i++) {
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

	    // Get the element number map (1-based)...
	    const Ioss::MapContainer &map = get_element_map();

	    // Allocate space for local side number, use 'data' as temporary
	    // storage for element numbers and overwrite with the side
	    // numbers.
	    Ioss::IntVector sides;
	    int *element = NULL;
	    int *ids = static_cast<int*>(data);
	    if (number_sides == static_cast<int>(entity_count)) {
	      // Only 1 side block in this sideset
	      sides.resize(entity_count);
	      element = static_cast<int*>(data);
	    } else {
	      // Multiple side blocks in sideset. Need to allocate memory to
	      // store sideset
	      sides.resize(number_sides);
	      element = new int[number_sides];
	    }

	    ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, element, TOPTR(sides));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

	    if (number_sides == static_cast<int>(entity_count)) {
	      for (size_t iel = 0; iel < entity_count; iel++) {
		int new_id = 10*map[element[iel]] + sides[iel];
		ids[iel] = new_id;
	      }
	    } else {
	      Ioss::IntVector is_valid_side(number_sides);
	      Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, element, TOPTR(sides),
							  number_sides, get_region());
	      size_t ieb = 0;
	      for (int iel = 0; iel < number_sides; iel++) {
		if (is_valid_side[iel] == 1) {
		  int new_id = 10*map[element[iel]] + sides[iel];
		  ids[ieb++] = new_id;
		}
	      }
	      assert(ieb == entity_count);
	    }

	    if (number_sides != static_cast<int>(entity_count))
	      delete [] element;

	  } else if (field.get_name() == "element_side") {
	    // In exodusII, the 'side set' is stored as a sideset.  A sideset
	    // has a list of elements and a corresponding local element side
	    // (1-based)

	    // Since we only have a single array, we need to allocate an extra
	    // array to store all of the data.  Note also that the element_id
	    // is the global id but only the local id is stored so we need to
	    // map from local_to_global prior to generating the side  id...

	    // Get the element number map (1-based)...
	    const Ioss::MapContainer &map = get_element_map();

	    // Allocate space for local side number and element numbers
	    // numbers.
	    int *element_side = static_cast<int*>(data);

	    Ioss::IntVector element(number_sides);
	    Ioss::IntVector sides(number_sides);

	    ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(sides));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

	    if (number_sides == static_cast<int>(entity_count)) {
	      size_t index = 0;
	      for (size_t iel = 0; iel < entity_count; iel++) {
		element_side[index++] = map[element[iel]];
		element_side[index++] = sides[iel];
	      }
	      assert(index/2 == entity_count);
	    } else {
	      Ioss::IntVector is_valid_side(number_sides);
	      Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, TOPTR(element), TOPTR(sides),
							  number_sides, get_region());

	      size_t index = 0;
	      for (int iel = 0; iel < number_sides; iel++) {
		if (is_valid_side[iel] == 1) {
		  // This side  belongs in the side block
		  element_side[index++] = map[element[iel]];
		  element_side[index++] = sides[iel];
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

	    // Allocate space for local side number and element numbers
	    // numbers.
	    int *element_side = static_cast<int*>(data);

	    Ioss::IntVector element(number_sides);
	    Ioss::IntVector sides(number_sides);

	    ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(sides));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

	    if (number_sides == static_cast<int>(entity_count)) {
	      size_t index = 0;
	      for (size_t iel = 0; iel < entity_count; iel++) {
		element_side[index++] = element[iel];
		element_side[index++] = sides[iel];
	      }
	      assert(index/2 == entity_count);
	    } else {
	      Ioss::IntVector is_valid_side(number_sides);
	      Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, TOPTR(element), TOPTR(sides),
							  number_sides, get_region());

	      size_t index = 0;
	      for (int iel = 0; iel < number_sides; iel++) {
		if (is_valid_side[iel] == 1) {
		  // This side  belongs in the side block
		  element_side[index++] = element[iel];
		  element_side[index++] = sides[iel];
		}
	      }
	      assert(index/2 == entity_count);
	    }

	  } else if (field.get_name() == "connectivity") {
	    // The side  connectivity needs to be generated 'on-the-fly' from
	    // the element number and local side  of that element. A sideset
	    // can span multiple element blocks, and contain multiple side
	    // types; the side block contains side of similar topology.
	    ierr = get_side_connectivity(fb, id, entity_count, 
					 static_cast<int*>(data), data_size/sizeof(int), true);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

	  } else if (field.get_name() == "connectivity_raw") {
	    // The side  connectivity needs to be generated 'on-the-fly' from
	    // the element number and local side  of that element. A sideset
	    // can span multiple element blocks, and contain multiple side
	    // types; the side block contains side of similar topology.
	    ierr = get_side_connectivity(fb, id, entity_count, 
					 static_cast<int*>(data), data_size/sizeof(int), false);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

	  } else if (field.get_name() == "distribution_factors") {
	    ierr = get_side_distributions(fb, id, entity_count, 
					  static_cast<double*>(data), data_size/sizeof(double));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  } else {
	    num_to_get = Ioss::Utils::field_warning(fb, field, "input");
	  }
	} else if (role == Ioss::Field::TRANSIENT) {
	  // Check if the specified field exists on this block.
	  // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	  // exist on the database as scalars with the appropriate
	  // extensions.

	  if (number_sides == static_cast<int>(entity_count)) {
	    num_to_get = read_transient_field(EX_SIDE_SET, m_variables[EX_SIDE_SET], field, fb, data);
	  } else {
	    // Need to read all values for the specified field and then
	    // filter down to the elements actualy in this side block.

	    // Determine which sides are member of this block
	    Ioss::IntVector is_valid_side(number_sides);
	    {
	      //----
	      Ioss::IntVector element(number_sides);
	      Ioss::IntVector sides(number_sides);
	      ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(sides));
	      //----
	      Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb,
							  TOPTR(element), TOPTR(sides),
							  number_sides, get_region());
	    }
	  
	    num_to_get = read_ss_transient_field(field, id, data, is_valid_side);
	  }
	}
      }
      return num_to_get;
    }

    int DatabaseIO::write_attribute_field(ex_entity_type type,
					  const Ioss::Field& field,
					  const Ioss::GroupingEntity *ge,
					  void *data) const
    {
      std::string att_name = ge->name() + SEP() + field.get_name();
      int num_entity = ge->get_property("entity_count").get_int();
      int offset = (int)field.get_index();
      assert(offset > 0);
    
      int id = get_id(ge, type, &ids_);
      int attribute_count = ge->get_property("attribute_count").get_int();
      assert(offset > 0);
      assert(offset-1+field.raw_storage()->component_count() <= attribute_count);
    
      if (offset == 1 && field.raw_storage()->component_count() == attribute_count) {
	// Write all attributes in one big chunk...
	int ierr = ex_put_attr(get_file_pointer(), type, id, static_cast<double*>(data));
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
      else {
	// Write a subset of the attributes.  If scalar, write one;
	// if higher-order (vector3d, ..) write each component.
	if (field.raw_storage()->component_count() == 1) {
	  int ierr = ex_put_one_attr(get_file_pointer(), type, id,
				     offset, static_cast<double*>(data));
	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);
	} else {
	  // Multi-component...  Need a local memory space to push
	  // data into and then write that out to the file...
	  std::vector<double> local_data(num_entity);
	  int comp_count = field.raw_storage()->component_count();
	  double *rdata = static_cast<double*>(data);
	  for (int i=0; i < comp_count; i++) {
	    int k = i;
	    for (int j=0; j < num_entity; j++) {
	      local_data[j] = rdata[k];
	      k += comp_count;
	    }
	  
	    int ierr = ex_put_one_attr(get_file_pointer(), type, id,
				       offset+i, TOPTR(local_data));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }
	}
      }
      return num_entity;
    }
  
    int DatabaseIO::read_attribute_field(ex_entity_type type,
					 const Ioss::Field& field,
					 const Ioss::GroupingEntity *ge,
					 void *data) const
    {
      int num_entity = ge->get_property("entity_count").get_int();
      if (num_entity == 0)
	return 0;
    
      int attribute_count = ge->get_property("attribute_count").get_int();
      int id = get_id(ge, type, &ids_);
    
      std::string att_name = ge->name() + SEP() + field.get_name();
      int offset = (int)field.get_index();
      assert(offset > 0);
      assert(offset-1+field.raw_storage()->component_count() <= attribute_count);
      if (offset == 1 && field.raw_storage()->component_count() == attribute_count) {
	// Read all attributes in one big chunk...
	int ierr = ex_get_attr(get_file_pointer(), type, id, static_cast<double*>(data));
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
      else {
	// Read a subset of the attributes.  If scalar, read one;
	// if higher-order (vector3d, ..) read each component and
	// put into correct location...
	if (field.raw_storage()->component_count() == 1) {
	  int ierr = ex_get_one_attr(get_file_pointer(), type, id,
				     offset, static_cast<double*>(data));
	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);
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
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  
	    int k = i;
	    for (int j=0; j < num_entity; j++) {
	      rdata[k] = local_data[j];
	      k += comp_count;
	    }
	  }
	}
      }
      return num_entity;
    }

    int DatabaseIO::read_transient_field(ex_entity_type type,
					 const VariableNameMap &variables,
					 const Ioss::Field& field,
					 const Ioss::GroupingEntity *ge,
					 void *data) const
    {
      const Ioss::VariableType *var_type = field.raw_storage();

      // Read into a double variable since that is all ExodusII can store...
      int num_entity = ge->get_property("entity_count").get_int();
      std::vector<double> temp(num_entity);

      int step = get_current_state();

      // get number of components, cycle through each component
      // and add suffix to base 'field_name'.  Look up index
      // of this name in 'nodeVariables' map
      int comp_count = var_type->component_count();
      int var_index=0;

      for (int i=0; i < comp_count; i++) {
	std::string var_name = var_type->label_name(field.get_name(), i+1, fieldSuffixSeparator);

	// Read the variable...
	int id = get_id(ge, type, &ids_);
	int ierr = 0;
	var_index = variables.find(var_name)->second;
	assert(var_index > 0);
	ierr = ex_get_var(get_file_pointer(), step, type,
			  var_index, id, num_entity, TOPTR(temp));
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);

	// Transfer to 'data' array.
	int k = 0;
	if (field.get_type() == Ioss::Field::INTEGER) {
	  int *ivar = static_cast<int*>(data);
	  for (int j=i; j < num_entity*comp_count; j+=comp_count) {
	    ivar[j] = static_cast<int>(temp[k++]);
	  }
	} else if (field.get_type() == Ioss::Field::REAL) {
	  double *rvar = static_cast<double*>(data);
	  for (int j=i; j < num_entity*comp_count; j+=comp_count) {
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

    int DatabaseIO::read_ss_transient_field(const Ioss::Field& field,
					    int id, void *variables,
					    Ioss::IntVector &is_valid_side) const
    {
      int num_valid_sides = 0;
      const Ioss::VariableType *var_type = field.raw_storage();
      size_t my_side_count = is_valid_side.size();
      std::vector<double> temp(my_side_count);

      int step = get_current_state();

      // get number of components, cycle through each component
      // and add suffix to base 'field_name'.  Look up index
      // of this name in 'nodeVariables' map
      int comp_count = var_type->component_count();
      int var_index=0;

      for (int i=0; i < comp_count; i++) {
	std::string var_name = var_type->label_name(field.get_name(), i+1, fieldSuffixSeparator);

	// Read the variable...
	int ierr = 0;
	var_index = m_variables[EX_SIDE_SET].find(var_name)->second;
	assert(var_index > 0);
	ierr = ex_get_var(get_file_pointer(), step, EX_SIDE_SET,
			  var_index, id, my_side_count, TOPTR(temp));
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);

	// Transfer to 'variables' array.
	int j = i;
	if (field.get_type() == Ioss::Field::INTEGER) {
	  int *ivar = static_cast<int*>(variables);
	  for (size_t k = 0; k < my_side_count; k++) {
	    if (is_valid_side[k] == 1) {
	      ivar[j] = static_cast<int>(temp[k]);
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

    int DatabaseIO::get_side_connectivity(const Ioss::SideBlock* fb,
					  int id, int,
					  int *fconnect,
					  size_t /* data_size */,
					  bool map_ids) const
    {
      // Get size of data stored on the file...
      int ierr;
      int number_sides;
      int number_distribution_factors;
      ierr = ex_get_set_param(get_file_pointer(), EX_SIDE_SET, id,
			      &number_sides, &number_distribution_factors);
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // Allocate space for element and local side number
      assert(number_sides > 0);
      //----
      Ioss::IntVector element(number_sides);
      Ioss::IntVector side(number_sides);

      ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(side));
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);
      //----

      Ioss::IntVector is_valid_side(number_sides);
      Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, TOPTR(element), TOPTR(side),
						  number_sides, get_region());

      Ioss::IntVector elconnect;
      int elconsize = 0; // Size of currently allocated connectivity block
      Ioss::ElementBlock *conn_block = NULL; // Block that we currently
      // have connectivity for

      Ioss::ElementBlock *block = NULL;

      Ioss::MapContainer side_elem_map; // Maps the side into the elements
      // connectivity array
      int current_side = -1;
      int nelnode = 0;
      int nfnodes = 0;
      int ieb = 0;
      int offset = 0;
      for (int iel = 0; iel < number_sides; iel++) {
	if (is_valid_side[iel] == 1) {

	  int elem_id = element[iel];

	  // ensure we have correct connectivity
	  block = get_region()->get_element_block(elem_id);
	  if (conn_block != block) {
	    int nelem = block->get_property("entity_count").get_int();
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
				 TOPTR(elconnect), nelem*nelnode*sizeof(int));
	    } else {
	      get_field_internal(block, block->get_field("connectivity_raw"),
				 TOPTR(elconnect), nelem*nelnode*sizeof(int));
	    }
	    conn_block = block;
	    current_side = -1;
	  }

	  // NOTE: Element connectivity is returned with nodes in global id space if "map_ids" false,
	  //       otherwise it is in local space.
	  if (current_side != side[iel]) {
	    side_elem_map = block->topology()->boundary_connectivity(side[iel]);
	    current_side = side[iel];
	    nfnodes = block->topology()->boundary_type(side[iel])->number_nodes();
	  }
	  for (int inode = 0; inode < nfnodes; inode++) {
	    int global_node = elconnect[(elem_id-offset)*nelnode +
					side_elem_map[inode]];
	    fconnect[ieb++] = global_node;
	  }
	}
      }
      return ierr;
    }

    // Get distribution factors for the specified side block
    int DatabaseIO::get_side_distributions(const Ioss::SideBlock* fb,
					   int id, int my_side_count,
					   double *dist_fact,
					   size_t /* data_size */) const
    {
      // Allocate space for elements and local side numbers
      // Get size of data stored on the file...
      int number_sides;
      int number_distribution_factors;
      int ierr = ex_get_set_param(get_file_pointer(), EX_SIDE_SET, id,
				  &number_sides, &number_distribution_factors);
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      const Ioss::ElementTopology *ftopo = fb->topology();
      int nfnodes = ftopo->number_nodes();

      if (number_distribution_factors == 0) {
	// Fill in the array with '1.0'...
	for (int i=0; i < nfnodes * my_side_count; i++)
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
      ierr = ex_get_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id, TOPTR(dist));
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // Another easy situation (and common for exodusII) is if the input
      // distribution factors are all the same value (typically 1).  In
      // that case, we only have to fill in the output array with that
      // value.
      {
	double value = dist[0];
	bool constant = true;
	for (int i=1; i < number_distribution_factors; i++) {
	  if (dist[i] != value) {
	    constant = false;
	    break;
	  }
	}
	if (constant) {
	  if (value == 0.0)
	    value = 1.0;  // Take care of some buggy mesh generators
	  for (int j=0; j < my_side_count * nfnodes; j++)
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
      Ioss::IntVector element(number_sides);
      Ioss::IntVector side(number_sides);

      ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, TOPTR(element), TOPTR(side));
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);
      //----

      Ioss::IntVector is_valid_side(number_sides);
      Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, TOPTR(element), TOPTR(side),
						  number_sides, get_region());

      int ieb = 0; // counter for distribution factors in this sideblock
      int idb = 0; // counter for distribution factors read from database
      Ioss::ElementBlock *block = NULL;

      for (int iel = 0; iel < number_sides; iel++) {
	int elem_id = element[iel];

	if (block == NULL || !block->contains(elem_id)) {
	  block = get_region()->get_element_block(elem_id);
	}

	if (block == NULL) {
	  std::ostringstream errmsg;
	  errmsg << "INTERNAL ERROR: Could not find element block containing element with id " << elem_id
		 << "Something is wrong in the Ioex::DatabaseIO class. Please report.\n";
	  IOSS_ERROR(errmsg);
	}

	const Ioss::ElementTopology *topo = block->topology()->boundary_type(side[iel]);

	if (topo == NULL) {
	  std::ostringstream errmsg;
	  errmsg << "INTERNAL ERROR: Could not find topology of element block boundary. "
		 << "Something is wrong in the Ioex::DatabaseIO class. Please report.\n";
	  IOSS_ERROR(errmsg);
	}

	int nside_nodes = topo->number_nodes();

	if (is_valid_side[iel] == 1) {
	  // This side belongs in the sideblock
	  for (int i=0; i < nside_nodes; i++) {
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

    //------------------------------------------------------------------------
    int DatabaseIO::put_field_internal(const Ioss::Region* /* region */,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      // For now, assume that all TRANSIENT fields on a region
      // are REDUCTION fields (1 value).  We need to gather these
      // and output them all at one time.  The storage location is a
      // 'globalVariables' array
      {
	Ioss::SerializeIO	serializeIO__(this);

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
	  errmsg << "The variable named '" << field.get_name()
		 << "' is of the wrong type. A region variable must be of type"
		 << " TRANSIENT or REDUCTION.\n"
		 << "This is probably an internal error; please notify gdsjaar@sandia.gov";
	  IOSS_ERROR(errmsg);
	}
	return num_to_get;
      }
    }

    int DatabaseIO::put_field_internal(const Ioss::NodeBlock* nb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);
	if (num_to_get > 0) {

	  Ioss::Field::RoleType role = field.get_role();

	  if (role == Ioss::Field::MESH) {
	    if (field.get_name() == "mesh_model_coordinates") {
	      // Data required by upper classes store x0, y0, z0, ... xn, yn, zn
	      // Data stored in exodusII file is x0, ..., xn, y0, ..., yn, z0, ..., zn
	      // so we have to allocate some scratch memory to read in the data
	      // and then map into supplied 'data'
	      std::vector<double> x(num_to_get);
	      std::vector<double> y;
	      if (spatialDimension > 1)
		y.resize(num_to_get);
	      std::vector<double> z;
	      if (spatialDimension == 3)
		z.resize(num_to_get);

	      // Cast 'data' to correct size -- double
	      double *rdata = static_cast<double*>(data);

	      int index = 0;
	      for (size_t i=0; i < num_to_get; i++) {
		x[i] = rdata[index++];
		if (spatialDimension > 1)
		  y[i] = rdata[index++];
		if (spatialDimension == 3)
		  z[i] = rdata[index++];
	      }
	      int ierr = ex_put_coord(get_file_pointer(), TOPTR(x), TOPTR(y), TOPTR(z));
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);

	    } else if (field.get_name() == "ids") {
	      // The ids coming in are the global ids; their position is the
	      // local id -1 (That is, data[0] contains the global id of local
	      // node 1)

	      // Another 'const-cast' since we are modifying the database just
	      // for efficiency; which the client does not see...
	      DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
	      new_this->handle_node_ids(static_cast<int*>(data), num_to_get);
	    } else if (field.get_name() == "connectivity") {
	      // Do nothing, just handles an idiosyncracy of the GroupingEntity
	    } else if (field.get_name() == "connectivity_raw") {
	      // Do nothing, just handles an idiosyncracy of the GroupingEntity
	    } else if (field.get_name() == "node_connectivity_status") {
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
	}
	return num_to_get;
      }
    }

    int DatabaseIO::put_field_internal(const Ioss::ElementBlock* eb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);

	if (num_to_get > 0) {
	  int ierr = 0;

	  // Get the element block id and element count
	  int id = get_id(eb, EX_ELEM_BLOCK, &ids_);
	  int my_element_count = eb->get_property("entity_count").get_int();
	  Ioss::Field::RoleType role = field.get_role();

	  if (role == Ioss::Field::MESH) {
	    // Handle the MESH fields required for an ExodusII file model.
	    // (The 'genesis' portion)
	    if (field.get_name() == "connectivity") {
	      if (my_element_count > 0) {
		// Map element connectivity from global node id to local node id.
		// Do it in 'data' ...
		int* connect = static_cast<int*>(data);
		
		assert(!nodeMap.empty());
		if (nodeMap[0] != -1) {
		  int element_nodes =
		    eb->get_property("topology_node_count").get_int();
		  assert(field.transformed_storage()->component_count() == element_nodes);
		  
		  for (size_t i=0; i < num_to_get * element_nodes; i++) {
		    int global_id = connect[i];
		    connect[i] = node_global_to_local(global_id, true);
		  }
		}
		ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, connect, NULL, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }
	    } else if (field.get_name() == "connectivity_edge") {
		if (my_element_count > 0) {
		// Map element connectivity from global edge id to local edge id.
		// Do it in 'data' ...
		int* connect = static_cast<int*>(data);

		if (edgeMap[0] != -1) {
		  int element_edges = field.transformed_storage()->component_count();

		  for (size_t i=0; i < num_to_get * element_edges; i++) {
		    int global_id = connect[i];
		    connect[i] = global_to_local(EX_EDGE_BLOCK, global_id);
		  }
		}
		ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, NULL, connect, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
		}
	    } else if (field.get_name() == "connectivity_face") {
		if (my_element_count > 0) {
		// Map element connectivity from global face id to local face id.
		// Do it in 'data' ...
		int* connect = static_cast<int*>(data);

		if (faceMap[0] != -1) {
		  int element_faces = field.transformed_storage()->component_count();

		  for (size_t i=0; i < num_to_get * element_faces; i++) {
		    int global_id = connect[i];
		    connect[i] = global_to_local(EX_FACE_BLOCK, global_id);
		  }
		}
		ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, NULL, NULL, connect);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }
	    } else if (field.get_name() == "connectivity_raw") {
	      // Do nothing, input only field.
	    } else if (field.get_name() == "ids") {
	      // Another 'const-cast' since we are modifying the database just
	      // for efficiency; which the client does not see...
	      DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
	      new_this->handle_element_ids(eb, static_cast<int*>(data), num_to_get);

	    } else if (field.get_name() == "skin") {
	      // This is (currently) for the skinned body. It maps the
	      // side element on the skin to the original element/local
	      // side number.  It is a two component field, the first
	      // component is the global id of the underlying element in
	      // the initial mesh and its local side number (1-based).

	      Ioss::IntVector element(my_element_count);
	      Ioss::IntVector side(my_element_count);
	      int *el_side = (int *)data;

	      int index = 0;
	      for (int i=0; i < my_element_count; i++) {
		element[i] = el_side[index++];
		side[i]    = el_side[index++];
	      }

	      // FIX: Hardwired map ids....
	      int map_count = ex_inquire_int(get_file_pointer(), EX_INQ_ELEM_MAP);
	      if (map_count == 0) {
		// This needs to be fixed... Currently hardwired....
		ierr = ex_put_map_param(get_file_pointer(), 0, 2);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }

	      int eb_offset = eb->get_offset();
	      ierr = ex_put_partial_elem_map(get_file_pointer(), 1, eb_offset+1, my_element_count,
					     TOPTR(element));
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      ierr = ex_put_partial_elem_map(get_file_pointer(), 2, eb_offset+1, my_element_count,
					     TOPTR(side));
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);

	      if (map_count == 0) {
		ierr = ex_put_name(get_file_pointer(), EX_ELEM_MAP, 1, "skin:parent_element_id");
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
		ierr = ex_put_name(get_file_pointer(), EX_ELEM_MAP, 2, "skin:parent_element_side_number");
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
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

    int DatabaseIO::put_field_internal(const Ioss::FaceBlock* eb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);

	if (num_to_get > 0) {
	  int ierr = 0;

	  // Get the face block id and face count
	  int id = get_id(eb, EX_FACE_BLOCK, &ids_);
	  int my_face_count = eb->get_property("entity_count").get_int();
	  Ioss::Field::RoleType role = field.get_role();

	  if (role == Ioss::Field::MESH) {
	    // Handle the MESH fields required for an ExodusII file model.
	    // (The 'genesis' portion)
	    if (field.get_name() == "connectivity") {
	      if (my_face_count > 0) {
		// Map face connectivity from global node id to local node id.
		// Do it in 'data' ...
		int* connect = static_cast<int*>(data);

		assert(!nodeMap.empty());
		if (nodeMap[0] != -1) {
		  int face_nodes =
		    eb->get_property("topology_node_count").get_int();
		  assert(field.transformed_storage()->component_count() == face_nodes);

		  for (size_t i=0; i < num_to_get * face_nodes; i++) {
		    int global_id = connect[i];
		    connect[i] = node_global_to_local(global_id, true);
		  }
		}
		ierr = ex_put_conn(get_file_pointer(), EX_FACE_BLOCK, id, connect, NULL, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }
	    } else if (field.get_name() == "connectivity_edge") {
		if (my_face_count > 0) {
		// Map face connectivity from global edge id to local edge id.
		// Do it in 'data' ...
		int* connect = static_cast<int*>(data);

		if (edgeMap[0] != -1) {
		  int face_edges = field.transformed_storage()->component_count();

		  for (size_t i=0; i < num_to_get * face_edges; i++) {
		    int global_id = connect[i];
		    connect[i] = global_to_local(EX_EDGE_BLOCK, global_id);
		  }
		}
		ierr = ex_put_conn(get_file_pointer(), EX_FACE_BLOCK, id, NULL, connect, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
		}
	    } else if (field.get_name() == "connectivity_raw") {
	      // Do nothing, input only field.
	    } else if (field.get_name() == "ids") {
	      // Another 'const-cast' since we are modifying the database just
	      // for efficiency; which the client does not see...
	      DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
	      new_this->handle_face_ids(eb, static_cast<int*>(data), num_to_get);

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

    int DatabaseIO::put_field_internal(const Ioss::EdgeBlock* eb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);

	size_t num_to_get = field.verify(data_size);

	if (num_to_get > 0) {
	  int ierr = 0;

	  // Get the edge block id and edge count
	  int id = get_id(eb, EX_EDGE_BLOCK, &ids_);
	  int my_edge_count = eb->get_property("entity_count").get_int();
	  Ioss::Field::RoleType role = field.get_role();

	  if (role == Ioss::Field::MESH) {
	    // Handle the MESH fields required for an ExodusII file model.
	    // (The 'genesis' portion)
	    if (field.get_name() == "connectivity") {
	      if (my_edge_count > 0) {
		// Map edge connectivity from global node id to local node id.
		// Do it in 'data' ...
		int* connect = static_cast<int*>(data);

		assert(!nodeMap.empty());
		if (nodeMap[0] != -1) {
		  int edge_nodes =
		    eb->get_property("topology_node_count").get_int();
		  assert(field.transformed_storage()->component_count() == edge_nodes);

		  for (size_t i=0; i < num_to_get * edge_nodes; i++) {
		    int global_id = connect[i];
		    connect[i] = node_global_to_local(global_id, true);
		  }
		}
		ierr = ex_put_conn(get_file_pointer(), EX_EDGE_BLOCK, id, connect, NULL, NULL);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }
	    } else if (field.get_name() == "connectivity_raw") {
	      // Do nothing, input only field.
	    } else if (field.get_name() == "ids") {
	      // Another 'const-cast' since we are modifying the database just
	      // for efficiency; which the client does not see...
	      DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
	      new_this->handle_edge_ids(eb, static_cast<int*>(data), num_to_get);

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

    int DatabaseIO::handle_node_ids(int* ids, size_t num_to_get)
    {
      /*!
       * There are two modes we need to support in this routine:
       * 1. Initial definition of node map (local->global) and
       * reverseNodeMap (global->local).
       * 2. Redefinition of node map via 'reordering' of the original
       * map when the nodes on this processor are the same, but their
       * order is changed (or count because of ghosting)
       *
       * So, there will be two maps the 'nodeMap' map is a 'direct lookup'
       * map which maps current local position to global id and the
       * 'reverseNodeMap' is an associative lookup which maps the
       * global id to 'original local'.  There is also a
       * 'reorderNodeMap' which is direct lookup and maps current local
       * position to original local.

       * The ids coming in are the global ids; their position is the
       * "local id-1" (That is, data[0] contains the global id of local
       * node 1 in this node block).
       *
       * int local_position = reverseNodeMap[NodeMap[i+1]]
       * (the nodeMap and reverseNodeMap are 1-based)
       *
       * To determine which map to update on a call to this function, we
       * use the following hueristics:
       * -- If the database state is 'STATE_MODEL:', then update the
       *    'reverseNodeMap' and 'nodeMap'
       *
       * -- If the database state is not STATE_MODEL, then leave the
       *    'reverseNodeMap' and 'nodeMap' alone since they correspond to the
       *    information already written to the database. [May want to add a
       *    STATE_REDEFINE_MODEL]
       *
       * -- In both cases, update the reorderNodeMap
       *
       * NOTE: The mapping is done on TRANSIENT fields only; MODEL fields
       *       should be in the orginal order...
       */
      assert(static_cast<int>(num_to_get) == nodeCount);
    
      if (dbState == Ioss::STATE_MODEL) {
	if (nodeMap.empty()) {
	  nodeMap.resize(nodeCount+1);
	  nodeMap[0] = -1;
	}

	if (nodeMap[0] == -1) {
	  for (size_t i=0; i < num_to_get; i++) {
	    nodeMap[i+1] = ids[i];
	    if (static_cast<int>(i+1) != ids[i]) {
	      assert(ids[i] != 0);
	      nodeMap[0] = 1;
	    }
	  }
	}

	if (nodeMap[0] != -1 || static_cast<int>(num_to_get) != nodeCount) {
	  Ioss::Map::build_reverse_map(&reverseNodeMap, ids, num_to_get, 0,
				       "node", myProcessor);
	}

	// Only a single nodeblock and all set
	if (static_cast<int>(num_to_get) == nodeCount) {
	  assert(nodeMap[0] == -1 || static_cast<int>(reverseNodeMap.size()) == nodeCount);
	}
	assert(get_region()->get_property("node_block_count").get_int() == 1);

	// Write to the database...
	int ierr = ex_put_node_num_map(get_file_pointer(), ids);
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }

      build_node_reorder_map(ids, num_to_get);
      return num_to_get;
    }

    namespace {
      template <typename T>
      void generate_block_truth_table(VariableNameMap &variables,
				      IntVector &truth_table,
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
    
      int internal_global_to_local(Ioss::ReverseMapContainer &reverseEntityMap, bool sequential,
				   size_t entity_count, int global)
      {
	int local = global;
	if (!sequential) {
	  std::pair<RMapI, RMapI> iter = std::equal_range(reverseEntityMap.begin(),
							  reverseEntityMap.end(),
							  global,
							  Ioss::IdPairCompare());
	  if (iter.first == iter.second) {
	    std::ostringstream errmsg;
	    errmsg << "Element with global id equal to " << global
		   << " does not exist in this mesh on this processor\n";
	    IOSS_ERROR(errmsg);
	  }
	  local = (iter.first)->second;
	}
	if (local > (int)entity_count || local <= 0) {
	  std::ostringstream errmsg;
	  errmsg << "Entity (element, face, edge, node) with global id equal to " << global
		 << " returns a local id of " << local
		 << " which is invalid. This should not happen, please report.\n";
	  IOSS_ERROR(errmsg);
	}
	return local;
      }

      void build_entity_reorder_map(Ioss::MapContainer &entityMap,
				    Ioss::ReverseMapContainer &reverseEntityMap,
				    Ioss::MapContainer &reorderEntityMap,
				    int start, int count)
      {
	// Note: To further add confusion, the reorderEntityaMap is 0-based
	// and the reverseEntityMap and entityMap are 1-baed. This is
	// just a consequence of how they are intended to be used...
	//
	// start is based on a 0-based array -- start of the reorderMap to build.
      
	if (reorderEntityMap.empty())
	  reorderEntityMap.resize(entityMap.size()-1);
      
	int my_end = start+count;
	for (int i=start; i < my_end; i++) {
	  int global_id = entityMap[i+1];
	  int orig_local_id = internal_global_to_local(reverseEntityMap, entityMap[0] == -1, entityMap.size()-1, global_id) - 1;
	
	  // If we assume that partial output is not being used (it
	  // currently isn't in Sierra), then the reordering should only be
	  // a permutation of the original ordering within this entity block...
	  assert(orig_local_id >= start && orig_local_id <= my_end);
	  reorderEntityMap[i] = orig_local_id;
	}
      }
    
      int handle_block_ids(const Ioss::EntityBlock *eb,
			   ex_entity_type map_type,
			   Ioss::State db_state,
			   Ioss::MapContainer &entityMap,
			   Ioss::ReverseMapContainer &reverseEntityMap,
			   Ioss::MapContainer &reorderEntityMap,
			   int* ids, size_t num_to_get, int file_pointer, int my_processor)
      {
	/*!
	 * NOTE: "element" is generic for "element", "face", or "edge"
	 *
	 * There are two modes we need to support in this routine:
	 * 1. Initial definition of element map (local->global) and
	 * reverseElementMap (global->local).
	 * 2. Redefinition of element map via 'reordering' of the original
	 * map when the elements on this processor are the same, but their
	 * order is changed.
	 *
	 * So, there will be two maps the 'elementMap' map is a 'direct lookup'
	 * map which maps current local position to global id and the
	 * 'reverseElementMap' is an associative lookup which maps the
	 * global id to 'original local'.  There is also a
	 * 'reorderElementMap' which is direct lookup and maps current local
	 * position to original local.

	 * The ids coming in are the global ids; their position is the
	 * local id -1 (That is, data[0] contains the global id of local
	 * element 1 in this element block).  The 'model-local' id is
	 * given by eb_offset + 1 + position:
	 *
	 * int local_position = reverseElementMap[ElementMap[i+1]]
	 * (the elementMap and reverseElementMap are 1-based)
	 *
	 * But, this assumes 1..numel elements are being output at the same
	 * time; we are actually outputting a blocks worth of elements at a
	 * time, so we need to consider the block offsets.
	 * So... local-in-block position 'i' is index 'eb_offset+i' in
	 * 'elementMap' and the 'local_position' within the element
	 * blocks data arrays is 'local_position-eb_offset'.  With this, the
	 * position within the data array of this element block is:
	 *
	 * int eb_position =
	 * reverseElementMap[elementMap[eb_offset+i+1]]-eb_offset-1
	 *
	 * To determine which map to update on a call to this function, we
	 * use the following hueristics:
	 * -- If the database state is 'Ioss::STATE_MODEL:', then update the
	 *    'reverseElementMap'.
	 * -- If the database state is not Ioss::STATE_MODEL, then leave
	 *    the 'reverseElementMap' alone since it corresponds to the
	 *    information already written to the database. [May want to add
	 *    a Ioss::STATE_REDEFINE_MODEL]
	 * -- Always update elementMap to match the passed in 'ids'
	 *    array.
	 *
	 * NOTE: the maps are built an element block at a time...
	 * NOTE: The mapping is done on TRANSIENT fields only; MODEL fields
	 *       should be in the orginal order...
	 */

	// Overwrite this portion of the 'elementMap', but keep other
	// parts as they were.  We are adding elements starting at position
	// 'eb_offset+offset' and ending at
	// 'eb_offset+offset+num_to_get'. If the entire block is being
	// processed, this reduces to the range 'eb_offset..eb_offset+my_element_count'

	int eb_offset = eb->get_offset();

	for (size_t i=0; i < num_to_get; i++) {
	  int local_id = eb_offset + i + 1;
	  entityMap[local_id] = ids[i];
	  if (local_id != ids[i]) {
	    entityMap[0] = 1;
	    assert(ids[i] != 0);
	  }
	}

	// Now, if the state is Ioss::STATE_MODEL, update the reverseEntityMap
	if (db_state == Ioss::STATE_MODEL) {
	  Ioss::Map::build_reverse_map(&reverseEntityMap, ids, num_to_get,
				       eb_offset, eb->type_string(), my_processor);

	  // Output this portion of the entity number map
	  int ierr = ex_put_partial_id_map(file_pointer, map_type, eb_offset+1, num_to_get, ids);
	  if (ierr < 0)
	    exodus_error(file_pointer, __LINE__, my_processor);
	}
	// Build the reorderEntityMap which does a direct mapping from
	// the current topologies local order to the local order
	// stored in the database...  This is 0-based and used for
	// remapping output and input TRANSIENT fields.
	build_entity_reorder_map(entityMap, reverseEntityMap, reorderEntityMap,
                                 eb_offset, num_to_get);
	return num_to_get;
      }
    }
  
    int  DatabaseIO::global_to_local(ex_entity_type type, int global) const
    {
      if (type == EX_NODE_BLOCK || type == EX_NODE_SET) {
	assert(!nodeMap.empty());
	return internal_global_to_local(reverseNodeMap, nodeMap[0] == -1, nodeCount, global);
      } else if (type == EX_ELEM_BLOCK || type == EX_ELEM_SET) {
	assert(!elementMap.empty());
	return internal_global_to_local(reverseElementMap, elementMap[0] == -1, elementCount, global);
      } else if (type == EX_FACE_BLOCK || type == EX_FACE_SET) {
	assert(!faceMap.empty());
	return internal_global_to_local(reverseFaceMap, faceMap[0] == -1, faceCount, global);
      } else if (type == EX_EDGE_BLOCK || type == EX_EDGE_SET) {
	assert(!edgeMap.empty());
	return internal_global_to_local(reverseEdgeMap, edgeMap[0] == -1, edgeCount, global);
      }
      return 0;
    }

    int DatabaseIO::handle_element_ids(const Ioss::ElementBlock *eb,
                                       int* ids, size_t num_to_get)
    {
      if (elementMap.empty()) {
	elementMap.resize(elementCount+1);
	elementMap[0] = -1;
      }
      return handle_block_ids(eb, EX_ELEM_MAP, dbState,
			      elementMap, reverseElementMap, reorderElementMap,
			      ids, num_to_get, get_file_pointer(), myProcessor);
    }

    int DatabaseIO::handle_face_ids(const Ioss::FaceBlock *eb, int* ids, size_t num_to_get)
    {
      if (faceMap.empty()) {
	faceMap.resize(faceCount+1);
	faceMap[0] = -1;
      }
      return handle_block_ids(eb, EX_FACE_MAP, dbState,
			      faceMap, reverseFaceMap, reorderFaceMap,
			      ids, num_to_get, get_file_pointer(), myProcessor);
    }

    int DatabaseIO::handle_edge_ids(const Ioss::EdgeBlock *eb, int* ids, size_t num_to_get)
    {
      if (edgeMap.empty()) {
	edgeMap.resize(edgeCount+1);
	edgeMap[0] = -1;
      }
      return handle_block_ids(eb, EX_EDGE_MAP, dbState,
			      edgeMap, reverseEdgeMap, reorderEdgeMap,
			      ids, num_to_get, get_file_pointer(), myProcessor);
    }

    void DatabaseIO::write_nodal_transient_field(ex_entity_type /* type */,
						 const Ioss::Field &field,
						 const Ioss::NodeBlock */* ge */,
						 int count, void *variables) const
    {
      Ioss::Field::BasicType ioss_type = field.get_type();
      assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER ||
	     ioss_type == Ioss::Field::COMPLEX);
      double *rvar = static_cast<double*>(variables);
      int  *ivar = static_cast<int*>(variables);

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

	for (int i=0; i < comp_count; i++) {
	  std::string var_name = var_type->label_name(field_name, i+1, fieldSuffixSeparator);

	  if (m_variables[EX_NODE_BLOCK].find(var_name) == m_variables[EX_NODE_BLOCK].end()) {
	    std::ostringstream errmsg;
	    errmsg << "ERROR: Could not find nodal variable '" << var_name << "'\n";
	    IOSS_ERROR(errmsg);
	  }

	  var_index = m_variables[EX_NODE_BLOCK].find(var_name)->second;

	  // Transfer from 'variables' array.
	  int k = 0;
	  int num_out = 0;
	  for (int j=(re_im*i)+complex_comp; j < re_im*count*comp_count; j+=(re_im*comp_count)) {
	    int where = reorderNodeMap[k++];
	    if (where >= 0) {
	      assert(where < count);
	      if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX)
		temp[where] = rvar[j];
	      else
		temp[where] = ivar[j];
	      num_out++;
	    }
	  }
	  assert(k == count);
	  assert(num_out == nodeCount);

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
						  int count,
						  void *variables) const
    {
      const Ioss::VariableType *var_type = field.transformed_storage();
      std::vector<double> temp(count);

      int step = get_current_state();
      step = get_database_step(step);

      int eb_offset = 0;
      if (ge->type() == Ioss::ELEMENTBLOCK) {
	const Ioss::ElementBlock *elb = dynamic_cast<const Ioss::ElementBlock*>(ge);
	assert(elb != NULL);
	eb_offset = elb->get_offset();
      }

      Ioss::Field::BasicType ioss_type = field.get_type();
      assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER ||
	     ioss_type == Ioss::Field::COMPLEX);
      double *rvar = static_cast<double*>(variables);
      int  *ivar = static_cast<int*>(variables);

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

	for (int i=0; i < comp_count; i++) {
	  std::string var_name = var_type->label_name(field_name, i+1, fieldSuffixSeparator);

	  var_index = m_variables[type].find(var_name)->second;
	  assert(var_index > 0);

	  // Transfer from 'variables' array.  Note that the
	  // 'reorderElementMap has '1..numel' ids in it, but the 'temp'
	  // array is stored in 'element block local id space', so we need
	  // to add/subtract the element block offset to keep things
	  // straight.
	  int k = eb_offset;
	  int where = 0;
	  for (int j=(re_im*i)+complex_comp; j < re_im*count*comp_count; j+=(re_im*comp_count)) {
	    // Map to storage location.

	    if (type == EX_ELEM_BLOCK)
	      where = reorderElementMap[k++] - eb_offset;
	    else
	      where = k++;

	    assert(where >= 0 && where < count);
	    if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX)
	      temp[where] = rvar[j];
	    else
	      temp[where] = ivar[j];
	  }
	  assert(k-eb_offset == count);

	  // Write the variable...
	  int id = get_id(ge, type, &ids_);
	  int ierr;
	  if (type == EX_SIDE_SET) {
	    int offset = ge->get_property("set_offset").get_int();
	    ierr = ex_put_n_var(get_file_pointer(), step, type, var_index, id, offset+1, count, TOPTR(temp));
	  } else {
	    ierr = ex_put_var(get_file_pointer(), step, type, var_index, id, count, TOPTR(temp));
	  }

	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}
      }
    }

    void DatabaseIO::store_reduction_field(ex_entity_type type, 
					   const Ioss::Field& field,
					   const Ioss::GroupingEntity *ge,
					   void *variables) const
    {
      const Ioss::VariableType *var_type = field.transformed_storage();

      Ioss::Field::BasicType ioss_type = field.get_type();
      assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER ||
	     ioss_type == Ioss::Field::COMPLEX);
      double *rvar = static_cast<double*>(variables);
      int  *ivar = static_cast<int*>(variables);

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

	for (int i=0; i < comp_count; i++) {
	  std::string var_name = var_type->label_name(field_name, i+1, fieldSuffixSeparator);

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
	  else
	    globalValues[var_index-1] = ivar[i];
	}
      }
    }

    void DatabaseIO::get_reduction_field(ex_entity_type,
					 const Ioss::Field& field,
					 const Ioss::GroupingEntity */* ge */,
					 void *variables) const
    {
      const Ioss::VariableType *var_type = field.raw_storage();

      Ioss::Field::BasicType ioss_type = field.get_type();
      assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER);
      double *rvar = static_cast<double*>(variables);
      int  *ivar = static_cast<int*>(variables);

      // get number of components, cycle through each component
      // and add suffix to base 'field_name'.  Look up index
      // of this name in 'm_variables[EX_GLOBAL]' map
      int comp_count = var_type->component_count();
      for (int i=0; i < comp_count; i++) {
	std::string field_name = field.get_name();
	std::string var_name = var_type->label_name(field_name, i+1, fieldSuffixSeparator);

	assert(m_variables[EX_GLOBAL].find(var_name) != m_variables[EX_GLOBAL].end());
	int var_index = m_variables[EX_GLOBAL].find(var_name)->second;

	assert(static_cast<int>(globalValues.size()) >= var_index);

	// Transfer to 'variables' array.
	if (ioss_type == Ioss::Field::REAL)
	  rvar[i] = globalValues[var_index-1];
	else
	  ivar[i] = static_cast<int>(globalValues[var_index-1]);
      }
    }

    void DatabaseIO::write_reduction_fields() const
    {
      int step = get_current_state();
      step = get_database_step(step);
      size_t count = globalValues.size();
      if (count > 0) {
	int ierr = ex_put_glob_vars(get_file_pointer(), step, count,
				    (double*)TOPTR(globalValues));
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
    }

    void DatabaseIO::read_reduction_fields() const
    {
      int step = get_current_state();
      size_t count = globalValues.size();
      if (count > 0) {
	int ierr = ex_get_glob_vars(get_file_pointer(), step, count,
				    TOPTR(globalValues));
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
    }

    int DatabaseIO::put_Xset_field_internal(ex_entity_type type,
					    const Ioss::EntitySet* ns,
					    const Ioss::Field& field,
					    void *data, size_t data_size) const
    {
      {
	Ioss::SerializeIO	serializeIO__(this);
	ex_update(get_file_pointer());

	int entity_count = ns->get_property("entity_count").get_int();
	size_t num_to_get = field.verify(data_size);
	if (num_to_get > 0) {

	  int id = get_id(ns, type, &ids_);
	  Ioss::Field::RoleType role = field.get_role();

	  if (role == Ioss::Field::MESH) {

	    if (field.get_name() == "ids") {
	      // Map node id from global node id to local node id.
	      // Do it in 'data' ...
	      int* ids = static_cast<int*>(data);

	      assert(!nodeMap.empty());
	      if (nodeMap[0] != -1) {
		for (size_t i=0; i < num_to_get; i++) {
		  int global_id = ids[i];
		  ids[i] = node_global_to_local(global_id, true);
		}
	      }
	      int ierr = ex_put_set(get_file_pointer(), type, id, static_cast<int*>(data), NULL);
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);

           } else if (field.get_name() == "orientation") {
	      int ierr = ex_put_set(get_file_pointer(), type, id, NULL, static_cast<int*>(data));
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);

	    } else if (field.get_name() == "distribution_factors") {
	      int ierr = ex_put_set_dist_fact(get_file_pointer(), type, id, static_cast<double*>(data));
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);

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

    int DatabaseIO::put_field_internal(const Ioss::NodeSet* ns,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      return put_Xset_field_internal(EX_NODE_SET, ns, field, data, data_size);
    }

    int DatabaseIO::put_field_internal(const Ioss::EdgeSet* ns,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      return put_Xset_field_internal(EX_EDGE_SET, ns, field, data, data_size);
    }

    int DatabaseIO::put_field_internal(const Ioss::FaceSet* ns,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      return put_Xset_field_internal(EX_FACE_SET, ns, field, data, data_size);
    }

    int DatabaseIO::put_field_internal(const Ioss::ElementSet* ns,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      return put_Xset_field_internal(EX_ELEM_SET, ns, field, data, data_size);
    }

    int DatabaseIO::put_field_internal(const Ioss::SideSet* fs,
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

    int DatabaseIO::put_field_internal(const Ioss::CommSet* cs,
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
	Ioss::IntVector entities(entity_count);
	Ioss::IntVector procs(entity_count);

	if (type == "node") {
	  // Convert global node id to local node id and store in 'entities'
	  int* entity_proc = static_cast<int*>(data);
	  int j=0;
	  for (size_t i=0; i < entity_count; i++) {
	    int global_id = entity_proc[j++];
	    entities[i] = node_global_to_local(global_id, true);
	    procs[i] = entity_proc[j++];
	  }

	  if (commsetNodeCount > 0) {
	    int ierr = ne_put_node_cmap(get_file_pointer(), get_id(cs, (ex_entity_type)0, &ids_),
					TOPTR(entities), TOPTR(procs), myProcessor);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
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
	    Ioss::IntVector internal(nodeCount);
	    for (int ij=0; ij < nodeCount; ij++)
	      internal[ij] = 1;
	    for (size_t J=0; J < entity_count; J++)
	      internal[entities[J]-1] = 0;

	    int b = 0;
	    for (int ij=0; ij < nodeCount; ij++) {
	      if (internal[ij] == 0) {
		entities[b++] = ij+1;
	      }
	    }

	    int k = 0;
	    for (int ij=0; ij < nodeCount; ij++) {
	      if (internal[ij] == 1) {
		internal[k++] = ij+1;
	      }
	    }

#ifndef NDEBUG
	    int border_nodes   = 0;
	    int internal_nodes = 0;
	    if (get_region()->property_exists("border_node_count")) {
	      border_nodes = get_region()->get_property("border_node_count").get_int();
	      assert(border_nodes   == b);
	    }

	    if (get_region()->property_exists("internal_node_count")) {
	      internal_nodes = get_region()->
		get_property("internal_node_count").get_int();
	      assert(internal_nodes == k);
	    }
#endif

	    int ierr = ne_put_node_map(get_file_pointer(), TOPTR(internal), TOPTR(entities), NULL,
				       myProcessor);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }

	} else if (type == "side") {
	  Ioss::IntVector sides(entity_count);
	  int* entity_proc = static_cast<int*>(data);
	  int j=0;
	  for (size_t i=0; i < entity_count; i++) {
	    // Assume klugy side id generation.
	    int global_id = entity_proc[j] / 10;
	    entities[i] = element_global_to_local(global_id);
	    sides[i] = entity_proc[j++] % 10;
	    procs[i] = entity_proc[j++];
	  }

#if 0
	  int offset = 0;
	  for (int ics=0; ics < commsetElemCount; ics++) {
	    int ierr = ne_put_elem_cmap(get_file_pointer(), elemCmapIds[ics],
					&entities[offset], &sides[offset],
					&procs[offset], myProcessor);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	    offset += elemCmapElemCnts[ics];
	  }
#else
	  int ierr = ne_put_elem_cmap(get_file_pointer(), get_id(cs, (ex_entity_type)0, &ids_),
				      TOPTR(entities), TOPTR(sides), TOPTR(procs), myProcessor);
	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);
#endif

	  // Construct the element map (internal vs. border).
	  // Border elements are those in the communication map (use entities array)
	  // Internal elements are the rest.  Allocate array to hold all elements,
	  // initialize all to '1', then zero out the elements in 'entities'.
	  // Iterate through array again and consolidate all '1's
	  Ioss::IntVector internal(elementCount);
	  for (j=0; j < elementCount; j++)
	    internal[j] = 1;
	  for (size_t ij=0; ij < entity_count; ij++)
	    internal[entities[ij]-1] = 0;

	  int b = 0;
	  for (j=0; j < elementCount; j++) {
	    if (internal[j] == 0) {
	      entities[b++] = j+1;
	    }
	  }

	  int k = 0;
	  for (j=0; j < elementCount; j++) {
	    if (internal[j] == 1) {
	      internal[k++] = j+1;
	    }
	  }

#ifndef NDEBUG
	  int border_elems   = 0;
	  int internal_elems = 0;
	  if (get_region()->property_exists("border_elem_count")) {
	    border_elems =
	      get_region()->get_property("border_elem_count").get_int();
	    assert(border_elems   == b);
	  }

	  if (get_region()->property_exists("internal_elem_count")) {
	    internal_elems =
	      get_region()->get_property("internal_elem_count").get_int();
	    assert(internal_elems == k);
	  }
#endif

	  ierr = ne_put_elem_map(get_file_pointer(), TOPTR(internal),
				 TOPTR(entities), myProcessor);
	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);

	} else {
	  std::ostringstream errmsg;
	  errmsg << "Invalid commset type " << type;
	  IOSS_ERROR(errmsg);
	}
      } else {
	num_to_get = Ioss::Utils::field_warning(cs, field, "output");
      }
      return num_to_get;
    }

    int DatabaseIO::put_field_internal(const Ioss::SideBlock* fb,
				       const Ioss::Field& field,
				       void *data, size_t data_size) const
    {
      Ioss::SerializeIO	serializeIO__(this);
      size_t num_to_get = field.verify(data_size);
      if (num_to_get > 0) {

	int id = get_id(fb, EX_SIDE_SET, &ids_);

	int entity_count = fb->get_property("entity_count").get_int();
	int offset = 0;
	offset = fb->get_property("set_offset").get_int();
      
	Ioss::Field::RoleType role = field.get_role();

	if (role == Ioss::Field::MESH) {
	  if (field.get_name() == "side_ids" &&
	      fb->name() == "universal_sideset") {
	    // The side ids are being stored as the distribution factor
	    // field on the universal sideset.  There should be no other
	    // side sets that request this field...  (Eventually,
	    // create an id field to store this info.

	    // Need to convert 'ints' to 'double' for storage on mesh...
	    int* ids = static_cast<int*>(data);
	    std::vector<double> real_ids(num_to_get);
	    for (size_t i = 0; i < num_to_get; i++) {
	      real_ids[i] = static_cast<double>(ids[i]);
	    }
	    //	  int ierr = ex_put_partial_set_dist_fact(get_file_pointer(),  EX_SIDE_SET, id, TOPTR(real_ids));
	    int ierr = ne_put_n_side_set_df(get_file_pointer(), id, offset+1, entity_count, TOPTR(real_ids));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
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
	    int df_offset = fb->get_property("set_df_offset").get_int();
	    int df_count  = fb->get_property("distribution_factor_count").get_int();
	    ierr = ne_put_n_side_set_df(get_file_pointer(), id, df_offset+1, df_count, static_cast<double*>(data));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

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
	    Ioss::IntVector element(num_to_get);
	    Ioss::IntVector side(num_to_get);
	    int *el_side = (int *)data;

	    int index = 0;
	    for (size_t i=0; i < num_to_get; i++) {
	      element[i] = element_global_to_local(el_side[index++]);
	      side[i]    = el_side[index++];
	    }

	    int ierr;
	    ierr = ne_put_n_side_set(get_file_pointer(), id, offset+1, entity_count, TOPTR(element), TOPTR(side));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

	  } else if (field.get_name() == "element_side_raw") {
	    // In exodusII, the 'side block' is stored as a sideset.  A
	    // sideset has a list of elements and a corresponding local
	    // element side (1-based)

	    // The 'data' passed into the function is stored as a
	    // 2D vector e0,f0,e1,f1,... (e=element, f=side)

	    // To avoid overwriting the passed in data, we allocate
	    // two arrays to store the data for this sideset.

	    // The element_id passed in is the local id.

	    // Allocate space for local side number and element numbers numbers.
	    Ioss::IntVector element(num_to_get);
	    Ioss::IntVector side(num_to_get);
	    int *el_side = (int *)data;

	    int index = 0;
	    for (size_t i=0; i < num_to_get; i++) {
	      element[i] = el_side[index++];
	      side[i]    = el_side[index++];
	    }

	    int ierr;
	    ierr = ne_put_n_side_set(get_file_pointer(), id, offset+1, entity_count, TOPTR(element), TOPTR(side));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

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

    // ------------------------------------------------------------------------
    // Node and Element mapping functions.  The ExodusII database
    // stores ids in a local-id system (1..NUMNP), (1..NUMEL) but
    // Sierra wants entities in a global system. These routines
    // take care of the mapping from local <-> global

    int DatabaseIO::local_to_global(ex_entity_type type, int local)  const
    {
      const Ioss::MapContainer &map = get_map(type);
      int global = map[local];
      return global;
    }

    int DatabaseIO::node_local_to_global(int local)  const
    {
      assert(local <= nodeCount && local > 0);
      const Ioss::MapContainer &node_map = get_node_map();
      int global = node_map[local];
      return global;
    }

    int DatabaseIO::element_local_to_global(int local)  const
    {
      assert(local <= elementCount && local > 0);
      const Ioss::MapContainer &element_map = get_element_map();
      int global = element_map[local];
      return global;
    }

    bool DatabaseIO::begin(Ioss::State state)
    {
      dbState = state;
      return true;
    }

    bool DatabaseIO::end(Ioss::State state)
    {
      // Transitioning out of state 'state'
      assert(state == dbState);
      switch (state) {
      case Ioss::STATE_DEFINE_MODEL:
	if (!is_input())
	  write_meta_data();
	break;
      case Ioss::STATE_DEFINE_TRANSIENT:
	if (!is_input())
	  write_results_metadata();
	break;
      default: // ignore everything else...
	break;
      }

      {
	Ioss::SerializeIO	serializeIO__(this);

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
    bool DatabaseIO::begin_state(Ioss::Region */* region */, int state, double time)
    {
      Ioss::SerializeIO	serializeIO__(this);

      state = get_database_step(state);
      if (!is_input()) {
	int ierr = ex_put_time(get_file_pointer(), state, &time);
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);

#ifdef REDS
	// This flush is a workaround for an issue on redstorm-s.  In a
	// 1600 processor run, 1 or 2 of the restart files had a time
	// value of 0.0 on the file even though it was verified that the
	// time value in the ex_put_time call above was the correct
	// time. Adding this ex_update call "fixed" this problem and
	// gives valid time values on all of the restart files.
	if (dbUsage != Ioss::WRITE_HISTORY)
	  ex_update(get_file_pointer());
#endif

	// Zero global variable array...
	std::fill(globalValues.begin(), globalValues.end(), 0.0);

      } else {
	// Store reduction variables
	read_reduction_fields();
      }
      return true;
    }

    bool DatabaseIO::end_state(Ioss::Region*, int, double time)
    {
      Ioss::SerializeIO	serializeIO__(this);

      if (!is_input()) {
	write_reduction_fields();
	finalize_write(time);
	if (minimizeOpenFiles)
	  free_file_pointer(); 
      }
      return true;
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


      Ioex::Mesh mesh(spatialDimension, the_title);
      Ioex::NodeBlock N(*node_blocks[0]);
      mesh.nodeblocks.push_back(N);

      // Edge Blocks --
      {
	Ioss::EdgeBlockContainer edge_blocks = region->get_edge_blocks();
	assert(check_block_order(edge_blocks));
	Ioss::EdgeBlockContainer::const_iterator I;
	// Set ids of all entities that have "id" property...
	for (I=edge_blocks.begin(); I != edge_blocks.end(); ++I) {
	  set_id(*I, EX_EDGE_BLOCK, &ids_);
	}
      
	edgeCount = 0;
	for (I=edge_blocks.begin(); I != edge_blocks.end(); ++I) {
	  edgeCount += (*I)->get_property("entity_count").get_int();
	  // Set ids of all entities that do not have "id" property...
	  get_id(*I, EX_EDGE_BLOCK, &ids_);
	  Ioex::EdgeBlock T(*(*I));
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
	  set_id(*I, EX_FACE_BLOCK, &ids_);
	}
      
	faceCount = 0;
	for (I=face_blocks.begin(); I != face_blocks.end(); ++I) {
	  faceCount += (*I)->get_property("entity_count").get_int();
	  // Set ids of all entities that do not have "id" property...
	  get_id(*I, EX_FACE_BLOCK, &ids_);
	  Ioex::FaceBlock T(*(*I));
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
	  set_id(*I, EX_ELEM_BLOCK, &ids_);
	}
      
	elementCount = 0;
	for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
	  elementCount += (*I)->get_property("entity_count").get_int();
	  // Set ids of all entities that do not have "id" property...
	  get_id(*I, EX_ELEM_BLOCK, &ids_);
	  Ioex::ElemBlock T(*(*I));
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
	  set_id(*I, EX_NODE_SET, &ids_);
	}

	for (I=nodesets.begin(); I != nodesets.end(); ++I) {
	  get_id(*I, EX_NODE_SET, &ids_);
	  const Ioex::NodeSet T(*(*I));
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
	  set_id(*I, EX_EDGE_SET, &ids_);
	}

	for (I=edgesets.begin(); I != edgesets.end(); ++I) {
	  get_id(*I, EX_EDGE_SET, &ids_);
	  const Ioex::EdgeSet T(*(*I));
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
	  set_id(*I, EX_FACE_SET, &ids_);
	}

	for (I=facesets.begin(); I != facesets.end(); ++I) {
	  get_id(*I, EX_FACE_SET, &ids_);
	  const Ioex::FaceSet T(*(*I));
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
	  set_id(*I, EX_ELEM_SET, &ids_);
	}

	for (I=elementsets.begin(); I != elementsets.end(); ++I) {
	  get_id(*I, EX_ELEM_SET, &ids_);
	  const Ioex::ElemSet T(*(*I));
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
	set_id(*I, EX_SIDE_SET, &ids_);
      }
      
      // Get entity counts for all face sets... Create SideSets.
      for (I=ssets.begin(); I != ssets.end(); ++I) {
	
	get_id(*I, EX_SIDE_SET, &ids_);
	int id = (*I)->get_property("id").get_int();
	int entity_count = 0;
	int df_count = 0;
	
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
	Ioex::SideSet T(*(*I));
	if (std::find(mesh.sidesets.begin(), mesh.sidesets.end(), T) == mesh.sidesets.end()) {
	  mesh.sidesets.push_back(T);
	}
      }
      m_groupCount[EX_SIDE_SET] = mesh.sidesets.size();

      gather_communication_metadata(&mesh.comm);

      {
	Ioss::SerializeIO	serializeIO__(this);

	if (myProcessor == 0) {
	  put_qa();
	  put_info();
	}
      
	// Write the metadata to the exodusII file...
	Ioex::Internals data(get_file_pointer(), maximumNameLength);
	int ierr = data.write_meta_data(mesh);

	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);

	// Write attribute names (if any)...
	write_attribute_names(get_file_pointer(), EX_NODE_SET,   get_region()->get_nodesets(),       fieldSuffixSeparator);
	write_attribute_names(get_file_pointer(), EX_EDGE_SET,   get_region()->get_edgesets(),       fieldSuffixSeparator);
	write_attribute_names(get_file_pointer(), EX_FACE_SET,   get_region()->get_facesets(),       fieldSuffixSeparator);
	write_attribute_names(get_file_pointer(), EX_ELEM_SET,   get_region()->get_elementsets(),    fieldSuffixSeparator);
	write_attribute_names(get_file_pointer(), EX_NODE_BLOCK, get_region()->get_node_blocks(),    fieldSuffixSeparator);
	write_attribute_names(get_file_pointer(), EX_EDGE_BLOCK, get_region()->get_edge_blocks(),    fieldSuffixSeparator);
	write_attribute_names(get_file_pointer(), EX_FACE_BLOCK, get_region()->get_face_blocks(),    fieldSuffixSeparator);
	write_attribute_names(get_file_pointer(), EX_ELEM_BLOCK, get_region()->get_element_blocks(), fieldSuffixSeparator);
      
	// Write coordinate names...
	char const *labels[3];
	labels[0] = "x";
	labels[1] = "y";
	labels[2] = "z";
	ierr = ex_put_coord_names(get_file_pointer(), (char**)labels);
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);

      }
    }

    void DatabaseIO::gather_communication_metadata(Ioex::CommunicationMetaData *meta)
    {
      if (isParallel) {

	meta->processorCount = util().parallel_size();
	meta->processorId = myProcessor;
	meta->globalNodes    = 1; // Just need a nonzero value.
	meta->globalElements = 1; // Just need a nonzero value.

	if (get_region()->property_exists("global_node_count")) {
	  meta->globalNodes = get_region()->get_property("global_node_count").get_int();
	}

	if (get_region()->property_exists("global_element_count")) {
	  meta->globalElements = get_region()->get_property("global_element_count").get_int();
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
	  int id = get_id(cs, (ex_entity_type)0, &ids_);

	  if (type == "node") {
	    meta->nodeMap.push_back(Ioex::CommunicationMap(id, count, 'n'));
	  } else if (type == "side") {
	    meta->elementMap.push_back(Ioex::CommunicationMap(id, count, 'e'));
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

    void DatabaseIO::add_region_fields()
    {
      int field_count = add_results_fields(EX_GLOBAL, get_region());
      globalValues.resize(field_count);
    }

    int DatabaseIO::add_results_fields(ex_entity_type type, Ioss::GroupingEntity *entity, int position)
    {
      return internal_add_results_fields(type, entity, position, m_groupCount[type], m_truthTable[type], m_variables[type]);
    }

    int DatabaseIO::internal_add_results_fields(ex_entity_type type, Ioss::GroupingEntity *entity, int position,
						int block_count, IntVector &truth_table, Ioex::VariableNameMap &variables)
    {
      int nvar = 0;
      {
	Ioss::SerializeIO	serializeIO__(this);

	int ierr = ex_get_variable_param(get_file_pointer(), type, &nvar);
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }

      if (nvar > 0) {
	if (truth_table.empty()) {
	  truth_table.resize(block_count * nvar);
	}
      
	// Read and store the truth table (Should be there since we only
	// get to this routine if there are variables...)
	if (!truth_table.empty()) {
	  {
	    Ioss::SerializeIO	serializeIO__(this);

	    if (type == EX_NODE_BLOCK || type == EX_GLOBAL) {
	      // These types don't have a truth table in the exodus api...
	      // They do in Ioss just for some consistency...
	      std::fill(truth_table.begin(), truth_table.end(), 1);
	    }
	    else {
	      int ierr = ex_get_truth_table(get_file_pointer(), type, block_count, nvar, &truth_table[0]);
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);
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
	  Ioss::SerializeIO	serializeIO__(this);

	  int ierr = ex_get_variable_names(get_file_pointer(), type, nvar, names);
	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);

	  // Add to VariableNameMap so can determine exodusII index given a
	  // Sierra field name.  exodusII index is just 'i+1'
	  for (int i=0; i < nvar; i++) {
	    Ioss::Utils::fixup_name(names[i]);
	    variables.insert(VNMValuePair(std::string(names[i]),   i+1));
	  }

	  int offset = position*nvar;
	  int *local_truth = NULL;
	  if (!truth_table.empty())
	    local_truth = &truth_table[offset];

	  std::vector<Ioss::Field> fields;
	  int count = entity->get_property("entity_count").get_int();
	  get_fields(count, names, nvar, Ioss::Field::TRANSIENT,
		     fieldSuffixSeparator, local_truth, fields);

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
    
      exo_params.edge_var_tab  = &m_truthTable[EX_EDGE_BLOCK][0];
      exo_params.face_var_tab  = &m_truthTable[EX_FACE_BLOCK][0];
      exo_params.elem_var_tab  = &m_truthTable[EX_ELEM_BLOCK][0];
      exo_params.nset_var_tab  = &m_truthTable[EX_NODE_SET][0];
      exo_params.eset_var_tab  = &m_truthTable[EX_EDGE_SET][0];
      exo_params.fset_var_tab  = &m_truthTable[EX_FACE_SET][0];
      exo_params.sset_var_tab  = &m_truthTable[EX_SIDE_SET][0];
      exo_params.elset_var_tab = &m_truthTable[EX_ELEM_SET][0];
    
      if (isParallel) {
	// Check consistency among all processors.  They should all
	// have the same number of each variable type...
	// The called function will throw an exception if the counts differ.
	check_variable_consistency(exo_params, myProcessor, get_filename(), util());
      }

      {
	Ioss::SerializeIO	serializeIO__(this);

	int ierr = ex_put_all_var_param_ext(get_file_pointer(), &exo_params);
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);

	globalValues.resize(m_variables[EX_GLOBAL].size());
	output_results_names(EX_GLOBAL, m_variables[EX_GLOBAL]);
	if (nodeCount > 0)
	  output_results_names(EX_NODE_BLOCK,    m_variables[EX_NODE_BLOCK]);
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
                                 fieldSuffixSeparator);
    }

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

	  for (int i=1; i <= var_type->component_count(); i++) {
	    std::string var_string = var_type->label_name(field_name, i, fieldSuffixSeparator);

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

      for (I=sidesets.begin(); I != sidesets.end(); ++I) {
	Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
	Ioss::SideBlockContainer::const_iterator J;

	bool found_one = false;
	for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	  // See if this sideblock has a corresponding entry in the sideset list.
	  if ((*J)->property_exists("invalid"))
	    continue;

	  found_one = true;
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
		std::string var_string = var_type->label_name(field_name, i, fieldSuffixSeparator);
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

    void
      DatabaseIO::output_results_names(ex_entity_type type,
				       VariableNameMap &variables) const
    {
      size_t var_count = variables.size();

      if (var_count > 0) {
	// Push into a char** array...
	std::vector<char*> var_names(var_count);
	VariableNameMap::const_iterator J  = variables.begin();
	VariableNameMap::const_iterator JE = variables.end();
	while (J != JE) {
	  size_t index = (*J).second;
	  assert(index > 0 && index <= var_count);
	  var_names[index-1] = (char*)(*J).first.c_str();
	  ++J;
	}

	int ierr = ex_put_variable_names(get_file_pointer(), type, var_count, TOPTR(var_names));
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
    }

    void DatabaseIO::build_node_reorder_map(int *new_ids, int count)
    {
      // This routine builds a map that relates the current node id order
      // to the original node ordering in affect at the time the file was
      // created. That is, the node map used to define the topology of the
      // model.  Now, if there are changes in node ordering at the
      // application level, we build the node reorder map to map the
      // current order into the original order.  An added complication is
      // that this is more than just a reordering... It may be that the
      // application has 'ghosted' nodes that it doesnt want put out on
      // the database, so the reorder map must handle a node that is not
      // in the original mesh and map that to an invalid value (currently
      // using -1 as invalid value...)


      // Note: To further add confusion,
      // the reorderNodeMap and new_ids are 0-based
      // the reverseNodeMap and nodeMap are 1-based. This is
      // just a consequence of how they are intended to be used...

      reorderNodeMap.resize(count);

      for (int i=0; i < count; i++) {
	int global_id = new_ids[i];

	// This will return 0 if node is not found in list.
	int orig_local_id = node_global_to_local(global_id, false) - 1;

	reorderNodeMap[i] = orig_local_id;
      }
    }

    // Handle special output time requests -- primarily restart (cycle, keep, overwrite)
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
      Ioex::Internals data(get_file_pointer(), maximumNameLength);
      data.update_last_time_attribute(sim_time);

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
	int my_element_count = block->get_property("entity_count").get_int();
      

	// Get the attribute names. May not exist or may be blank...
	char **names = get_exodus_names(attribute_count, maximumNameLength);
	int id = block->get_property("id").get_int();
	{
	  Ioss::SerializeIO	serializeIO__(this);
	  if (block->get_property("entity_count").get_int() != 0) {
	    int ierr = ex_get_attr_names(get_file_pointer(), entity_type, id, &names[0]);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
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
	for (int i=0; i < attribute_count; i++) {
	  fix_bad_name(names[i]);
	  Ioss::Utils::fixup_name(names[i]);
	}
      
	if (names[0][0] != '\0' && names[0][0] != ' ' && std::isalnum(names[0][0])) {
	  std::vector<Ioss::Field> attributes;
	  get_fields(my_element_count, names, attribute_count,
		     Ioss::Field::ATTRIBUTE, fieldSuffixSeparator, NULL,
		     attributes);
	  int offset = 1;
	  std::vector<Ioss::Field>::const_iterator IF;
	  for (IF = attributes.begin(); IF != attributes.end(); ++IF) {
	    Ioss::Field field = *IF;
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
	      offset += 3;
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
	    unknown_attributes = attribute_count - index;
	  }

	  else {
	    unknown_attributes = attribute_count;
	  }

	  if (unknown_attributes > 0) {
	    att_name = "extra_attribute";
	    std::string storage = "Real[";
	    storage += Ioss::Utils::to_string(unknown_attributes);
	    storage += "]";
	    size_t index = attribute_count - unknown_attributes + 1;
	    block->field_add(Ioss::Field(att_name, Ioss::Field::REAL, storage,
					 Ioss::Field::ATTRIBUTE, my_element_count, index));
	  
	    if (myProcessor == 0) {
	      IOSS_WARNING << "For element block '" << block->name()
			   << "' of type '" << type << "'\n\tthere were "
			   << unknown_attributes << " attributes that are not known to the IO Subsystem\n\t"
			   << "in addition to the " << attribute_count - unknown_attributes
			   << " known. The extra attributes can be accessed\n\tas the field '"
			   << att_name << "' with " << unknown_attributes << " components.\n\n";
	    }
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
} // End of namespace

namespace {
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
	int ge_id = ge->get_property("id").get_int();
	int ierr = ex_put_attr_names(exoid, type, ge_id, TOPTR(names));
	if (ierr < 0)
	  exodus_error(exoid, __LINE__, -1);
      }
    }
  }
  
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

  size_t match(const char *name1, const char *name2)
  {
    size_t l1 = std::strlen(name1);
    size_t l2 = std::strlen(name2);
    size_t len = l1 < l2 ? l1 : l2;
    for (size_t i=0; i < len; i++) {
      if (name1[i] != name2[i]) {
	while (i > 0 && isdigit(name1[i-1]) && isdigit(name2[i-1])) {
	  i--;
	  // Back up to first non-digit so to handle "evar0000, evar0001, ..., evar 1123"
	}
	return i;
      }
    }
    return len;
  }
  
  const Ioss::VariableType *match_composite_field(char** names, Ioss::IntVector &which_names,
						  const char suffix_separator)
  {
    // ASSUME: Fields are in order...
    // The field we are trying to match will be a composite field of
    // type base_x_1, base_y_1, base_z_1, ...., base_y_3, base_z_3.
    // The composite field type currently always has a numeric Real[N]
    // type field as the last suffix and the other field as the first
    // suffix.
    // If we take the last suffix of the last name, it should give us
    // the 'N' in the Real[N] field.  Dividing 'which_names.size()' by
    // 'N' will give the number of components in the inner field.
    
    char suffix[2];
    suffix[0] = suffix_separator;
    suffix[1] = 0;

    std::vector<std::string> tokens;
    Ioss::tokenize(names[which_names[which_names.size()-1]] ,suffix, tokens);

    if (tokens.size() <= 2)
      return NULL;
    
    assert(tokens.size() > 2);
    size_t N = std::strtol(tokens[tokens.size()-1].c_str(), NULL, 10);
    if (N == 0)
      return NULL;

    assert(N > 0 && which_names.size() % N == 0);

    size_t inner_token = tokens.size() - 2;
    size_t inner_comp = which_names.size() / N;
    
    // Gather the first 'inner_ccomp' inner field suffices...
    std::vector<Ioss::Suffix> suffices;
    for (size_t i=0; i < inner_comp; i++) {
      std::vector<std::string> ltokens;
      Ioss::tokenize(names[which_names[i]], suffix, ltokens);
      // The second-last token is the suffix for this component...
      Ioss::Suffix tmp(ltokens[inner_token]);
      suffices.push_back(tmp);
    }

    // check that the suffices on the next copies of the inner field
    // match the first copy...
    size_t j = inner_comp;
    for (size_t copy = 1; copy < N; copy++) {
      for (size_t i=0; i < inner_comp; i++) {
	std::vector<std::string> ltokens;
	Ioss::tokenize(names[which_names[j++]], suffix, ltokens);
	// The second-last token is the suffix for this component...
	if (suffices[i] != ltokens[inner_token]) {
	  return NULL;
	}
      }
    }
    
    // All 'N' copies of the inner field match, now see the
    // suffices actually defines a field...
    const Ioss::VariableType *type = Ioss::VariableType::factory(suffices);
    if (type != NULL) {
      type = Ioss::VariableType::factory(type->name(), N);
    }
    return type;
  }

  const Ioss::VariableType *match_single_field(char** names, Ioss::IntVector &which_names,
					       const char suffix_separator)
  {
    // Strip off the suffix from each name indexed in 'which_names'
    // and see if it defines a valid type...
    std::vector<Ioss::Suffix> suffices;
    
    char suffix[2];
    suffix[0] = suffix_separator;
    suffix[1] = 0;

    for (size_t i=0; i < which_names.size(); i++) {
      std::vector<std::string> tokens;
      Ioss::tokenize(names[which_names[i]], suffix, tokens);
      size_t num_tokens = tokens.size();
      // The last token is the suffix for this component...
      Ioss::Suffix tmp(tokens[num_tokens-1]);
      suffices.push_back(tmp);
    }
    const Ioss::VariableType *type = Ioss::VariableType::factory(suffices);
    return type;
  }

  Ioss::Field get_next_field(char** names, int num_names, int count,
			     Ioss::Field::RoleType fld_role,
			     const char suffix_separator, int *truth_table)
  {
    // NOTE: 'names' are all lowercase at this point.
    
    // Assumption 1: To convert to a non-SCALAR type, the variable name
    // must have an field_suffix_sep in the name separating the suffixes from
    // the main name.

    // Find first unused name (used names have '\0' as first character...
    int index = 0;
    bool found_valid = false;
    for (index = 0; index < num_names; index++) {
      assert(truth_table == NULL || truth_table[index] == 1 || truth_table[index] == 0);
      if ((truth_table == NULL || truth_table[index] == 1) && names[index][0] != '\0') {
	found_valid = true;
	break;
      }
    }

    if (!found_valid) {
      // Return an invalid field...
      return Ioss::Field("", Ioss::Field::INVALID, SCALAR(), fld_role, 1);
    }

    // At this point, name[index] should be a valid potential field
    // name and all names[i] with i < index are either already used or
    // not valid for this grouping entity (truth_table entry == 0).
    assert (index < num_names && names[index][0] != '\0' && (truth_table == NULL || truth_table[index] == 1));
    char *name = names[index];

    // Split the name up into tokens separated by the
    // 'suffix_separator'.  Note that the basename itself could
    // contain a suffix_separator (back_stress_xx or
    // back_stress_xx_01). Need to ignore embedded separators
    // (back_stress) and also recognize composite variable types
    // (back_stress_xx_01). At the current time, a composite variable
    // type can only contain two non-composite variable types, so we
    // only need to look to be concerned with the last 1 or 2 tokens...
    char suffix[2];
    suffix[0] = suffix_separator;
    suffix[1] = 0;
    std::vector<std::string> tokens;
    Ioss::tokenize(name, suffix, tokens);
    size_t num_tokens = tokens.size();
    
    // Check that separator is not first or last character of the name...
    bool invalid = tokens[0].empty() || tokens[num_tokens-1].empty();
    if (num_tokens == 1 || invalid) {
      // It is not a (Sierra-generated) name for a non-SCALAR variable
      // Return a SCALAR field
      Ioss::Field field(name, Ioss::Field::REAL, SCALAR(), fld_role, count);
      names[index][0] = '\0';
      return field;
    }

    // KNOW: The field_suffix_sep is not in first or last position.
    // KNOW: num_tokens > 1 at this point.  Possible that we still
    // just have a scalar with an embedded separator character...
    int suffix_size = 1;
    if (num_tokens > 2)
      suffix_size = 2;
    
    // If num_tokens > 2, then we can potentially have a composite
    // variable type which would have a double suffix (_xx_01).

    // Gather all names which match in the first
    // (num_tokens-suffix_size) tokens and see if their suffices form
    // a valid variable type...
    while (suffix_size > 0) {
      Ioss::IntVector which_names; // Contains index of names that
				    // potentially match as components
				    // of a higher-order type.

      std::string base_name = tokens[0];
      for (size_t i=1; i < num_tokens-suffix_size; i++) {
	base_name += suffix_separator;
	base_name += tokens[i];
      }
      base_name += suffix_separator;
      size_t bn_len = base_name.length(); // Length of basename portion only
      size_t length = std::strlen(name); // Length of total name (with suffix)

      // Add the current name...
      which_names.push_back(index);

      // Gather all other names that are valid for this entity, and
      // have the same overall length and match in the first 'bn_len'
      // characters. 
      //
      // Check that they have the same number of tokens,
      // It is possible that the first name(s) that match with two
      // suffices have a basename that match other names with only a
      // single suffix lc_cam_x, lc_cam_y, lc_sfarea.
      for (int i = index+1; i < num_names; i++) {
	char *tst_name = names[i];
	std::vector<std::string> subtokens;
	Ioss::tokenize(tst_name,suffix,subtokens);
	if ((truth_table == NULL || truth_table[i] == 1) &&  // Defined on this entity
	    std::strlen(tst_name) == length &&              // names must be same length
	    std::strncmp(name, tst_name, bn_len) == 0 &&   // base portion must match
	    subtokens.size() == num_tokens) {
	  which_names.push_back(i);
	}
      }

      const Ioss::VariableType *type = NULL;
      if (suffix_size == 2) {
	if (which_names.size() > 1) 
	  type = match_composite_field(names, which_names, suffix_separator);
      } else {
	assert(suffix_size == 1);
	type = match_single_field(names, which_names, suffix_separator);
      }

      if (type != NULL) {
	// A valid variable type was recognized.
	// Mark the names which were used so they aren't used for another field on this entity.
	// Create a field of that variable type.
	assert(type->component_count() == static_cast<int>(which_names.size()));
	Ioss::Field field(base_name.substr(0,bn_len-1), Ioss::Field::REAL, type, fld_role, count);
	for (size_t i=0; i < which_names.size(); i++) {
	  names[which_names[i]][0] = '\0';
	}
	return field;
      } else {
	if (suffix_size == 1) {
	  Ioss::Field field(name, Ioss::Field::REAL, SCALAR(), fld_role, count);
	  names[index][0] = '\0';
	  return field;
	}
      }
      suffix_size--;
    }
    return Ioss::Field("", Ioss::Field::INVALID, SCALAR(), fld_role, 1);
  }

  bool define_field(size_t nmatch, size_t match_length,
		    char **names, std::vector<Ioss::Suffix> &suffices,
		    int entity_count, Ioss::Field::RoleType fld_role,
		    std::vector<Ioss::Field> &fields)
  {
    // Try to define a field of size 'nmatch' with the suffices in 'suffices'.
    // If this doesn't define a known field, then assume it is a scalar instead
    // and return false.
    if (nmatch > 1) {
      const Ioss::VariableType *type = Ioss::VariableType::factory(suffices);
      if (type == NULL) {
	nmatch = 1;
      } else {
	char *name = names[0];
	name[match_length] = '\0';
	Ioss::Field field(name, Ioss::Field::REAL, type, fld_role, entity_count);
	if (field.is_valid()) {
	  fields.push_back(field);
	}
	for (size_t j = 0; j < nmatch; j++)
	  names[j][0] = '\0';
	return true;
      }
    }
    
    // NOTE: nmatch could be reset inside previous if block.
    // This is not an 'else' block, it is a new if block.
    if (nmatch == 1) {
      Ioss::Field field(names[0], Ioss::Field::REAL, SCALAR(), fld_role, entity_count);
      if (field.is_valid()) {
	fields.push_back(field);
      }
      names[0][0] = '\0';
      return false;
    }
    return false; // Can't get here...  Quiet the compiler
  }

  // Read scalar fields off an input database and determine whether
  // they are components of a higher order type (vector, tensor, ...).
  // This routine is used if there is no field component separator.  E.g.,
  // fieldx, fieldy, fieldz instead of field_x field_y field_z 

  void get_fields(int entity_count, // The number of objects in this entity.
		  char** names,     // Raw list of field names from exodus
		  size_t num_names, // Number of names in list
		  Ioss::Field::RoleType fld_role, // Role of field
		  const char suffix_separator,
		  int *local_truth, // Truth table for this entity;
				    // null if not applicable.
		  std::vector<Ioss::Field> &fields) // The fields that were found.
  {
    if (suffix_separator != 0) {
      while (1) {
	// NOTE: 'get_next_field' determines storage type (vector, tensor,...)
	Ioss::Field field = get_next_field(names, num_names, entity_count, fld_role,
					   suffix_separator, local_truth);
	if (field.is_valid()) {
	  fields.push_back(field);
	} else {
	  break;
	}
      }
    } else {
      size_t nmatch = 1;
      size_t ibeg   = 0;
      size_t pmat   = 0;
      std::vector<Ioss::Suffix> suffices;
    top:
      
      while (ibeg+nmatch < num_names) {
	if (local_truth != NULL) {
	  while (ibeg < num_names && local_truth[ibeg] == 0)
	    ibeg++;
	}
	for (size_t i=ibeg+1; i < num_names; i++) {
	  size_t mat = match(names[ibeg], names[i]);
	  if (local_truth != NULL && local_truth[i] == 0)
	    mat = 0;
	  
	  // For all fields, the total length of the name is the same
	  // for all components of that field.  The 'basename' of the
	  // field will also be the same for all cases.
	  //
	  // It is possible that the length of the match won't be the
	  // same for all components of a field since the match may
	  // include a portion of the suffix; (sigxx, sigxy, sigyy
	  // should match only 3 characters of the basename (sig), but
	  // sigxx and sigxy will match 4 characters) so consider a
	  // valid match if the match length is >= previous match length.
	  if ((std::strlen(names[ibeg]) == std::strlen(names[i])) &&
	      mat > 0 && (pmat == 0 || mat >= pmat)) {
	    nmatch++;
	    if (nmatch == 2) {
	      // Get suffix for first field in the match
	      pmat = mat;
	      Ioss::Suffix tmp(&names[ibeg][pmat]);
	      suffices.push_back(tmp);
	    }
	    // Get suffix for next fields in the match
	    Ioss::Suffix tmp(&names[i][pmat]);
	    suffices.push_back(tmp);
	  } else {
	    
	    bool multi_component = define_field(nmatch, pmat, &names[ibeg], suffices,
						entity_count, fld_role, fields);
	    if (!multi_component) {
	      // Although we matched multiple suffices, it wasn't a
	      // higher-order field, so we only used 1 name instead of
	      // the 'nmatch' we thought we might use.
	      i = ibeg + 1; 
	    }
	    
	    // Cleanout the suffices vector.
	    std::vector<Ioss::Suffix>().swap(suffices);
	    
	    // Reset for the next time through the while loop...
	    nmatch=1;
	    pmat = 0;
	    ibeg=i;
	    break;
	  }
	}
      }
      // We've gone through the entire list of names; see if what we
      // have forms a multi-component field; if not, then define a
      // scalar field and jump up to the loop again to handle the others
      // that had been gathered.
      if (ibeg < num_names) {
	if (local_truth == NULL || local_truth[ibeg] == 1) {
	  bool multi_component = define_field(nmatch, pmat, &names[ibeg], suffices,
					      entity_count, fld_role, fields);
	  std::vector<Ioss::Suffix>().swap(suffices);
	  if (nmatch > 1 && !multi_component) {
	    ibeg++;
	    goto top;
	  }
	} else {
	  ibeg++;
	  goto top;
	}
      }
    }
  }

  bool type_match(const std::string& type, const char *substring)
  {
    // Returns true if 'substring' is a sub-string of 'type'.
    // The comparisons are case-insensitive
    // 'substring' is required to be in all lowercase.
    const char *s = substring;
    const char *t = type.c_str();

    assert(s != NULL && t != NULL);
    while (*s != '\0' && *t != '\0') {
      if (*s++ != tolower(*t++)) {
	return false;
      }
    }
    return true;
  }

  void decode_surface_name(Ioex::SideSetMap &fs_map, Ioex::SideSetSet &fs_set, const std::string &name)
  {
    std::vector<std::string> tokens;
    Ioss::tokenize(name, "_", tokens);
    if (tokens.size() >= 4) {
      // Name of form: "name_eltopo_sidetopo_id" or
      // "name_block_id_sidetopo_id" "name" is typically "surface".
      // The sideset containing this should then be called "name_id"
      
      // Check whether the second-last token is a side topology and
      // the third-last token is an element topology.
      const Ioss::ElementTopology *side_topo = Ioss::ElementTopology::factory(tokens[tokens.size()-2], true);
      if (side_topo != NULL) {
	const Ioss::ElementTopology *element_topo = Ioss::ElementTopology::factory(tokens[tokens.size()-3], true);
	if (element_topo != NULL || tokens[tokens.size()-4] == "block") {
	  // The remainder of the tokens will be used to create
	  // a side set name and then this sideset will be
	  // a side block in that set.
	  std::string fs_name;
	  size_t last_token = tokens.size()-3;
	  if (element_topo == NULL)
	    last_token--;
	  for (size_t tok=0; tok < last_token; tok++) {
	    fs_name += tokens[tok];
	  }
	  fs_name += "_";
	  fs_name += tokens[tokens.size()-1]; // Add on the id.
	  
	  fs_set.insert(fs_name);
	  fs_map.insert(Ioex::SideSetMap::value_type(name,fs_name));
	}
      }
    }
  }

bool set_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Ioex::EntityIdSet *idset)
{
  // See description of 'get_id' function.  This function just primes
  // the idset with existing ids so that when we start generating ids,
  // we don't overwrite an existing one.

  // Avoid a few string constructors/destructors
  static std::string prop_name("name");
  static std::string id_prop("id");

  bool succeed = false;
  if (entity->property_exists(id_prop)) {
    int id = entity->get_property(id_prop).get_int();

    // See whether it already exists...
    succeed = idset->insert(std::make_pair((int)type,id)).second;
    if (!succeed) {
      // Need to remove the property so it doesn't cause problems
      // later...
      Ioss::GroupingEntity *new_entity = const_cast<Ioss::GroupingEntity*>(entity);
      new_entity->property_erase(id_prop);
      assert(!entity->property_exists(id_prop));
    }
  }
  return succeed;
}

// Potentially extract the id from a name possibly of the form name_id.
// If not of this form, return 0;
int extract_id(const std::string &name_id)
{
  std::vector<std::string> tokens;
  Ioss::tokenize(name_id,"_",tokens);
  
  if (tokens.size() == 1)
    return 0;

  // Check whether last token is an integer...
  std::string str_id = tokens[tokens.size()-1];
  size_t len = str_id.length();
  bool is_int = true;
  for (size_t i=0; i < len; i++) {
    if (str_id[i] < '0' || str_id[i] > '9') {
      is_int = false;
      break;
    }
  }      
  if (is_int)
    return std::atoi(str_id.c_str());

  return 0;
}

int get_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Ioex::EntityIdSet *idset)
{
  // Sierra uses names to refer to grouping entities; however,
  // exodusII requires integer ids.  When reading an exodusII file,
  // the DatabaseIO creates a name by concatenating the entity
  // type (e.g., 'block') and the id separated by an underscore.  For
  // example, an exodusII element block with an id of 100 would be
  // encoded into "block_100"

  // This routine tries to determine the id of the entity using 3
  // approaches:
  //
  // 1. If the entity contains a property named 'id', this is used.
  // The DatabaseIO actually stores the id in the "id" property;
  // however, other grouping entity creators are not required to do
  // this so the property is not guaranteed to exist.
  //
  // 2.If property does not exist, it tries to decode the entity name
  // based on the above encoding.  Again, it is not required that the
  // name follow this convention so success is not guaranteed.
  //
  // 3. If all other schemes fail, the routine picks an id for the entity
  // and returns it.  It also stores this id in the "id" property so an
  // entity will always return the same id for multiple calls.
  // Note that this violates the 'const'ness of the entity so we use
  // a const-cast.

  // Avoid a few string constructors/destructors
  static std::string prop_name("name");
  static std::string id_prop("id");

  int id = 1;

  if (entity->property_exists(id_prop)) {
    id = entity->get_property(id_prop).get_int();
    return id;

  } else {
    // Try to decode an id from the name.
    std::string name_string = entity->get_property(prop_name).get_string();
    id = extract_id(name_string);
    if (id <= 0) id = 1;
  }

  // At this point, we either have an id equal to '1' or we have an id
  // extracted from the entities name. Increment it until it is
  // unique...
  while (idset->find(std::make_pair(int(type), id)) != idset->end()) {
    ++id;
  }

  // 'id' is a unique id for this entity type...
  idset->insert(std::make_pair((int)type,id));
  Ioss::GroupingEntity *new_entity = const_cast<Ioss::GroupingEntity*>(entity);
  new_entity->property_add(Ioss::Property(id_prop, id));
  return id;
}

bool find_displacement_field(Ioss::NameList &fields,
			     const Ioss::GroupingEntity *block,
			     int ndim,
			     std::string *disp_name)
{
  // This is a kluge to work with many of the SEACAS codes.  The
  // convention used (in Blot and others) is that the first 'ndim'
  // nodal variables are assumed to be displacements *if* the first
  // character of the names is 'D' and the last characters match the
  // coordinate labels (typically 'X', 'Y', and 'Z').  This routine
  // looks for the field that has the longest match with the string
  // "displacement" and is of the correct storage type (VECTOR_2D or
  // VECTOR_3D).  If found, it returns the name.
  //

  static char displace[] = "displacement";
  
  Ioss::NameList::const_iterator IF;
  Ioss::NameList::const_iterator IFend = fields.end();
  size_t max_span = 0;
  
  for (IF = fields.begin(); IF != IFend; ++IF) {
    const char *name = (*IF).c_str();
    std::string lc_name(name);
    
    Ioss::Utils::fixup_name(lc_name);
    size_t span = match(lc_name.c_str(), displace);
    if (span > max_span) {
      const Ioss::VariableType *var_type =
	block->get_field((*IF)).transformed_storage();
      int comp_count = var_type->component_count();
      if (comp_count == ndim) {
	max_span  = span;
	*disp_name = *IF;
      }
    }
  }
  return max_span > 0 ? true : false;
}

void fix_bad_name(char* name)
{
  assert(name != NULL);

  size_t len = std::strlen(name);
  for (size_t i=0; i < len; i++) {
    if (name[i] < 32 || name[i] > 126) {
      // Zero out entire name if a bad character found anywhere in the name.
      for (size_t j=0; j < len; j++) {
	name[j] = '\0';
      }
      return;
    }
  }
}

void filter_element_list(Ioss::Region *region,
			 Ioss::IntVector &elements, Ioss::IntVector &sides,
			 bool remove_omitted_elements)
{
  // Iterate through 'elements' and remove the elements which are in an omitted block.
  // Precondition is that there is at least one omitted element block.
  // The 'elements' list contains local element ids, not global.
  // Since there are typically a small number of omitted blocks, do
  // the following:
  // For each omitted block, determine the min and max element id in
  // that block.  Iterate 'elements' vector and set the id to zero if
  // min <= id <= max.  Once all omitted blocks have been processed,
  // then iterate the vector and compress out all zeros.  Keep 'sides'
  // array consistent.

  // Get all element blocks in region...
  bool omitted = false;
  Ioss::ElementBlockContainer element_blocks = region->get_element_blocks();
  for (size_t blk=0; blk < element_blocks.size(); blk++) {
    Ioss::ElementBlock *block = element_blocks[blk];
    if (block_is_omitted(block)) {
      int min_id = block->get_offset() + 1;
      int max_id = min_id + block->get_property("entity_count").get_int() - 1;
      for (size_t i=0; i < elements.size(); i++) {
	if (min_id <= elements[i] && elements[i] <= max_id) {
	  omitted = true;
	  elements[i] = 0;
	  sides[i]    = 0;
	}
      }
    }
  }
  if (remove_omitted_elements && omitted) {
    elements.erase(std::remove(elements.begin(), elements.end(), 0), elements.end());
    sides.erase(   std::remove(sides.begin(),    sides.end(),    0), sides.end());
  }
}

void separate_surface_element_sides(Ioss::IntVector &element,
				    Ioss::IntVector &sides,
				    Ioss::Region *region,
				    Ioex::TopologyMap &topo_map,
				    Ioex::TopologyMap &side_map,
				    Ioss::SurfaceSplitType split_type)
{
  if (element.size() > 0) {
    Ioss::ElementBlock *block = NULL;
    // Topology of sides in current element block
    const Ioss::ElementTopology *common_ftopo = NULL;
    const Ioss::ElementTopology *topo = NULL; // Topology of current side
    int current_side = -1;

    for (size_t iel = 0; iel < element.size(); iel++) {
      int elem_id = element[iel];
      if (block == NULL || !block->contains(elem_id)) {
	block = region->get_element_block(elem_id);
	assert(block != NULL);
	assert(!block_is_omitted(block)); // Filtered out above.

	// NULL if hetero sides on element
	common_ftopo = block->topology()->boundary_type(0);
	if (common_ftopo != NULL)
	  topo = common_ftopo;
	current_side = -1;
      }

      if (common_ftopo == NULL && sides[iel] != current_side) {
	current_side = sides[iel];
	assert(current_side > 0 && current_side <= block->topology()->number_boundaries());
	topo = block->topology()->boundary_type(sides[iel]);
	assert(topo != NULL);
      }
      std::pair<std::string, const Ioss::ElementTopology*> name_topo;
      if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
	name_topo = std::make_pair(block->topology()->name(), topo);
      } else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
	name_topo = std::make_pair(block->name(), topo);
      }
      topo_map[name_topo]++;
      if (side_map[name_topo] == 0)
	side_map[name_topo] = sides[iel];
      else if (side_map[name_topo] != sides[iel]) {
	// Not a consistent side for all sides in this
	// sideset. Set to large number. Note that maximum
	// sides/element is 6, so don't have to worry about
	// a valid element having 999 sides (unless go to
	// arbitrary polyhedra some time...) Using a large
	// number instead of -1 makes it easier to check the
	// parallel consistency...
	side_map[name_topo] = 999;
      }
    }
  }
}

void set_sideblock_ids(Ioss::SideSetContainer &ssets,
		       std::vector<Ioex::SideSet> &sidesets,
		       const Ioss::ParallelUtils &util,
		       Ioex::EntityIdSet *ids)
{
  bool is_parallel  = util.parallel_size() > 1;
  Ioss::SideSetContainer::const_iterator I;
  // Set ids of all entities that have an existing "id" property

  // Get entity counts for all side blocks...
  Ioss::IntVector sset_entity_count;
  for (I=ssets.begin(); I != ssets.end(); ++I) {
    Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
    Ioss::SideBlockContainer::const_iterator J;
    
    for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
      size_t count = (*J)->get_property("entity_count").get_int();
      sset_entity_count.push_back(count);
    }
  }

  // Resolve count among all side blocks on all processors...
  // NOTE: This is a collective call.
  if (is_parallel && !sset_entity_count.empty()) {
    util.global_array_minmax(sset_entity_count, Ioss::ParallelUtils::DO_MAX);
  }

  // If count > 0 on any processor, set the id if they have an id property.
  // If the containing faceset has an id and the sideblock doesn't have an id,
  // then set the sideblock id to the faceset id...
  size_t fb_index = 0;
  for (I=ssets.begin(); I != ssets.end(); ++I) {
    int fs_id = 0;
    if ((*I)->property_exists("id")) {
      fs_id = (*I)->get_property("id").get_int();
    }
    
    Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
    Ioss::SideBlockContainer::const_iterator J;
    
    if (fs_id > 0) {
      for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	if (!(*J)->property_exists("id")) {
	  (*J)->property_add(Ioss::Property("id", fs_id));
	}
	size_t count = sset_entity_count[fb_index++];
	if (count > 0) {
	  set_id((*J), EX_SIDE_SET, ids);
	}
      }
    } else {
      for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	size_t count = sset_entity_count[fb_index++];
	if (count > 0) {
	  set_id((*J), EX_SIDE_SET, ids);
	}
      }
    }
  }
  assert(fb_index == sset_entity_count.size());
  
  // The id has been set on all side blocks that had the id property.
  // Now, go through again and set the id on all side blocks.
  fb_index = 0;
  for (I=ssets.begin(); I != ssets.end(); ++I) {
    Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
    Ioss::SideBlockContainer::const_iterator J;
    
    for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
      size_t count = sset_entity_count[fb_index++];
      if (count > 0) {
	get_id((*J), EX_SIDE_SET, ids);
	Ioex::SideSet T(*(*J));
	if (std::find(sidesets.begin(),
		      sidesets.end(), T) == sidesets.end()) {
	  sidesets.push_back(T);
	}
      } else {
	// Set the "invalid" property.
	Ioss::SideBlock *new_entity = const_cast<Ioss::SideBlock*>(*J);
	new_entity->property_add(Ioss::Property("invalid", 1));
      }
    }
  }
  assert(fb_index == sset_entity_count.size());
}

void add_map_fields(int exoid, Ioss::ElementBlock *block, int my_element_count)
{
  // Check for optional element maps...
  int map_count = ex_inquire_int(exoid, EX_INQ_ELEM_MAP);
  if (map_count <= 0)
    return;

  // Get the names of the maps...
  int max_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  char **names = get_exodus_names(map_count, max_length);
  int ierr = ex_get_names(exoid, EX_ELEM_MAP, names);
  if (ierr < 0)
    exodus_error(exoid, __LINE__, -1);

  // Convert to lowercase.
  for (int i=0; i < map_count; i++) {
    Ioss::Utils::fixup_name(names[i]);
  }

  if (map_count == 2 && std::strncmp(names[0], "skin:", 5) == 0 && std::strncmp(names[1], "skin:", 5) == 0) {
    // Currently, only support the "skin" map -- It will be a 2
    // component field consisting of "parent_element":"local_side"
    // pairs.  The parent_element is an element in the original mesh,
    // not this mesh.
    block->field_add(Ioss::Field("skin", Ioss::Field::INTEGER, "Real[2]",
				 Ioss::Field::MESH, my_element_count));
  }
  delete_exodus_names(names, map_count);
}

#ifndef NDEBUG
template <typename T>
bool check_block_order(const std::vector<T*> &blocks)
{
  // Verify that element blocks are defined in sorted offset order...
  typename std::vector<T*>::const_iterator I;

  int eb_offset = -1;
  for (I=blocks.begin(); I != blocks.end(); ++I) {
    int this_off = (*I)->get_offset();
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
	case 2:  type = "edge";	      break;
	case 3:  type = "face";	      break;
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

