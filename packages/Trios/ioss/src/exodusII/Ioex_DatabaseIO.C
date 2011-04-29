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
  const size_t max_string_length = MAX_STR_LENGTH;
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

#ifndef NDEBUG
  bool check_block_order(const Ioss::ElementBlockContainer &blocks);
#endif
  void check_variable_consistency(size_t gv_count, size_t nv_count, size_t ev_count,
				  size_t nsv_count, size_t ssv_count,
				  int my_processor, const std::string &filename,
				  const Ioss::ParallelUtils &util);

  char ** get_exodus_names(size_t count, int size=max_string_length)
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

  std::string get_entity_name(int exoid, ex_entity_type type, int id, const std::string &basename)
  {
    char buffer[max_string_length+1];
    buffer[0] = '\0';
    int error = ex_get_name(exoid, type, id, buffer);
    if (error < 0)
      exodus_error(exoid, __LINE__, -1);
    if (buffer[0] != '\0') {
      Ioss::Utils::fixup_name(buffer);
      return (std::string(buffer));
    } else {
      return Ioss::Utils::encode_entity_name(basename, id);
    }
  }
  
  void check_attribute_index_order(Ioss::GroupingEntity *block);
}

namespace Ioex {
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
			 Ioss::DatabaseUsage db_usage, MPI_Comm communicator) :
    Ioss::DatabaseIO(region, filename, db_usage, communicator),
    exodusFilePtr(-1), databaseTitle(""), exodusMode(EX_CLOBBER),
    spatialDimension(0), nodeCount(0), elementCount(0),
    nodeBlockCount(0), elementBlockCount(0), nodesetCount(0), sidesetCount(0),
    nodeCmapIds(NULL), nodeCmapNodeCnts(NULL),
    elemCmapIds(NULL), elemCmapElemCnts(NULL), commsetNodeCount(0),
    commsetElemCount(0),
    elementTruthTable(NULL), nodesetTruthTable(NULL), sidesetTruthTable(NULL),
    timeLastFlush(0),
    sequentialNG2L(true), sequentialEG2L(true), fileExists(false),
    minimizeOpenFiles(false), blockAdjacenciesCalculated(false), nodeConnectivityStatusCalculated(false)
  {

    // A history file is only written on processor 0...
    if (db_usage == Ioss::WRITE_HISTORY)
      isParallel = false;

    dbState = Ioss::STATE_UNKNOWN;

    // Set exodusII warning level.
#if 0
    ex_opts(EX_VERBOSE|EX_DEBUG);
#endif
    
    if (!is_input() && util().get_environment("EX_MINIMIZE_OPEN_FILES", isParallel)) {
      std::cerr << "IOEX: Minimizing open files because EX_MINIMIZE_OPEN_FILES environment variable is set.\n";
      minimizeOpenFiles = true;
    }

    if (!is_input() && util().get_environment("EX_SURFACE_SPLIT_BACKWARD_COMPATIBILITY", isParallel)) {
      std::cerr << "IOEX: Splitting surfaces because EX_SURFACE_SPLIT_BACKWARD_COMPATIBILITY environment variable is set.\n";
      surfaceSplitBackwardCompatibility = true;
    }

    if (!is_input() && util().get_environment("EX_MODE", exodusMode, isParallel)) {
      std::cerr << "IOEX: Exodus create mode set to " << exodusMode
		<< " from value of EX_MODE environment variable.\n";
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
    
      delete [] elementTruthTable;
      delete [] nodesetTruthTable;
      delete [] sidesetTruthTable;
      delete [] nodeCmapIds;
      delete [] nodeCmapNodeCnts;
      delete [] elemCmapIds;
      delete [] elemCmapElemCnts;
    } catch (...) {
    }
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
      }

      if (exodusFilePtr < 0) {
	dbState = Ioss::STATE_INVALID;
	fileExists = false;
	// NOTE: Code will not continue past this call...
	std::ostringstream errmsg;
	errmsg << "Problem opening specified file '" << decoded_filename << "'";
	IOSS_ERROR(errmsg);
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
    static char qa_temp[4][max_string_length+1];
    static char *qa[1][4] =
      {{qa_temp[0],qa_temp[1],qa_temp[2],qa_temp[3]}};

    Ioss::Utils::time_and_date(qa[0][3], qa[0][2], max_string_length);

    std::string codename = "unknown";
    std::string version  = "unknown";

    if (get_region()->property_exists("code_name")) {
      codename = get_region()->get_property("code_name").get_string();
    }
    if (get_region()->property_exists("code_version")) {
      version = get_region()->get_property("code_version").get_string();
    }

    std::strncpy(qa[0][0], codename.c_str(), max_string_length);
    std::strncpy(qa[0][1], version.c_str(),  max_string_length);
    qa[0][0][max_string_length] = '\0';
    qa[0][1][max_string_length] = '\0';

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

      Ioex::Internals data(get_file_pointer());
      data.check_processor_info(util().parallel_size(), myProcessor);

      read_region();
      read_communication_metadata();
    }

    get_step_times();
    get_elemblocks();
    check_side_topology();

    {
      Ioss::SerializeIO	serializeIO__(this);
      get_nodeblocks();
    }
    get_sidesets();
    {
      Ioss::SerializeIO	serializeIO__(this);
      get_nodesets();
      get_commsets();

      add_region_fields();
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

    char the_title[max_line_length+1];
    std::memset(the_title, 0, max_line_length+1);

    int error = ex_get_init(get_file_pointer(), the_title, &spatialDimension,
			    &nodeCount, &elementCount, &elementBlockCount,
			    &nodesetCount, &sidesetCount);
    if (error < 0)
      exodus_error(get_file_pointer(), __LINE__, myProcessor);

    nodeBlockCount = 1;

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

    if (elementCount > 0 && elementBlockCount <= 0) {
      // NOTE: Code will not continue past this call...
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      std::ostringstream errmsg;
      errmsg << "No element blocks were found in the model, file: '" << decoded_filename << "'";
      IOSS_ERROR(errmsg);
    }

    Ioss::Region *this_region = get_region();
    this_region->property_add(Ioss::Property(std::string("title"), the_title));
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
      int error = ex_get_all_times(get_file_pointer(), &tsteps[0]);
      if (error < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // See if the "last_written_time" attribute exists and if it
      // does, check that it matches the largest time in 'tsteps'.
      Ioex::Internals data(get_file_pointer());
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

  const Ioss::MapContainer& DatabaseIO::get_node_map() const
  {
    // Allocate space for node number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (nodeMap.empty()) {
      nodeMap.resize(nodeCount+1);

      if (is_input()) {
	bool backward_compat = get_node_global_id_backward_compatibility() && !isParallel && dbUsage == Ioss::READ_MODEL;
	if (!backward_compat) {
	  Ioss::SerializeIO	serializeIO__(this);
	  // Check whether there is a "original_global_id_map" map on
	  // the database. If so, use it instead of the "node_num_map".
	  bool map_read = false;
	  int map_count = ex_inquire_int(get_file_pointer(), EX_INQ_NODE_MAP);
	  if (map_count > 0 && !get_node_global_id_backward_compatibility()) {
	    char **names = get_exodus_names(map_count);
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
	sequentialNG2L = true;
	for (int i=1; i < nodeCount+1; i++) {
	  if (i != nodeMap[i]) {
	    sequentialNG2L = false;
	    nodeMap[0] = 1;
	    break;
	  }
	}

	if (!sequentialNG2L) {
	  Ioss::Map::build_reverse_map(&reverseNodeMap, &nodeMap[1], nodeCount, 0,
				       "node", myProcessor);
	} else {
	  // Sequential map
	  nodeMap[0] = -1;
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
	bool backward_compat = get_node_global_id_backward_compatibility() && !isParallel && dbUsage == Ioss::READ_MODEL;
	if (!backward_compat) {
	  Ioss::SerializeIO	serializeIO__(this);
	  
	  // Check whether there is a "original_global_id_map" map on
	  // the database. If so, use it instead of the "elem_num_map".
	  bool map_read = false;
	  int map_count = ex_inquire_int(get_file_pointer(), EX_INQ_ELEM_MAP);
	  if (map_count > 0) {
	    char **names = get_exodus_names(map_count);
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
	    int error = ex_get_elem_num_map(get_file_pointer(), &elementMap[1]);
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
	sequentialEG2L = true;
	for (int i=1; i < elementCount+1; i++) {
	  if (i != elementMap[i]) {
	    sequentialEG2L = false;
	    elementMap[0] = 1;
	    break;
	  }
	}

	if (!sequentialEG2L) {
	  Ioss::Map::build_reverse_map(&reverseElementMap,
				       &elementMap[1], elementCount, 0,
				       "element", myProcessor);
	} else {
	  // Sequential map
	  sequentialEG2L = true;
	  elementMap[0] = -1;
	}

      } else {
	// Output database; elementMap not set yet... Build a default map.
	for (int i=1; i < elementCount+1; i++) {
	  elementMap[i] = i;
	}
	// Sequential map
	sequentialEG2L = true;
	elementMap[0] = -1;
      }
    }
    return elementMap;
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
    add_results_fields(EX_NODAL, block, nodeCount);

    get_region()->add(block);
  }

  void DatabaseIO::get_elemblocks() 
  {
    // Attributes of an element block are:
    // -- id
    // -- name
    // -- element type
    // -- number of elements
    // -- number of attributes per element
    // -- number of nodes per element (derivable from type)
    // -- number of faces per element (derivable from type)
    // -- number of edges per element (derivable from type)

    // In a parallel execution, it is possible that an element block will have
    // no elements on a particular processor...

    // NOTE: This routine may be called multiple times on a single database.
    //       make sure it is not dependent on being called one time only...

    // Get exodusII element block metadata
    if (elementBlockCount == 0)
      return;

    Ioss::IntVector element_block_ids(elementBlockCount);

    int error;
    {
      Ioss::SerializeIO	serializeIO__(this);

      error = ex_get_ids(get_file_pointer(), EX_ELEM_BLOCK, &element_block_ids[0]);
      if (error < 0) {
	exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
    }


    size_t all_element_type_length = elementBlockCount * (max_string_length+1);
    std::vector<char> all_element_type(all_element_type_length);

    Ioss::IntVector attributes(elementBlockCount);
    Ioss::IntVector my_node_count(elementBlockCount);
    Ioss::IntVector local_element_count(elementBlockCount);
    Ioss::IntVector global_element_count(elementBlockCount);
    int iblk ;

    {
      Ioss::SerializeIO	serializeIO__(this);

      for ( iblk = 0; iblk < elementBlockCount; iblk++) {
	int id = element_block_ids[iblk];
	int number_elements;
	int nodes_per_element;
	int attributes_per_element;

	char * const element_type = &all_element_type[0] + iblk * (max_string_length+1);

	error = ex_get_elem_block(get_file_pointer(), id,
				  element_type,
				  &number_elements, &nodes_per_element,
				  &attributes_per_element);
	if (error < 0) {
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}

	local_element_count[iblk] = number_elements;

	if (number_elements == 0) {
	  attributes[iblk] = 0;
	  my_node_count[iblk] = 0;
	  std::memset( element_type , 0 , max_string_length + 1 );
	} else {
	  attributes[iblk] = attributes_per_element;
	  my_node_count[iblk] = nodes_per_element;
	}
      }
    }

    // This is a collective call...
    util().attribute_reduction(all_element_type_length, &all_element_type[0]);

    // The element attribute count will either be 0 if there are no
    // elements in the element block on this processor, or it will be
    // the number of attributes (> 0). Therefore, if we take the 'max'
    // over all processors, each processor will then have the correct
    // attribute count...
    // This is a collective call...
    util().global_array_minmax(&attributes[0], attributes.size(), Ioss::ParallelUtils::DO_MAX);

    // Similar for number of nodes per element...
    // This is a collective call... Could combine with previous call to eliminate one call...
    util().global_array_minmax(&my_node_count[0], my_node_count.size(), Ioss::ParallelUtils::DO_MAX);

    // Determine global element count for each element block....
    // Can also get this from an nemesis call, but the data may not always be there in all cases.
    util().global_count(local_element_count, global_element_count);
    
    // The 'offset' is used to map an element location within an element
    // block to the element 'file descriptor'.  For example, the file
    // descriptor of the 37th element in the 4th block is calculated by:
    // file_descriptor = offset of block 4 + 37 This can also be used to
    // determine which element block an element with a file_descriptor
    // maps into. An particular element block contains all elements in
    // the range:
    //     offset < file_descriptor <= offset+number_elements_per_block
    int offset = 0;
    int used_blocks = 0;
    
    for ( iblk = 0; iblk < elementBlockCount; iblk++) {
      int id = element_block_ids[iblk];
      std::string alias = Ioss::Utils::encode_entity_name("block", id);
      char * const element_type = &all_element_type[0] + iblk * (max_string_length+1);

      Ioss::ElementBlock *block = NULL;
      std::string block_name = get_entity_name(get_file_pointer(), EX_ELEM_BLOCK, id, "block");

      std::string save_type = element_type;
      std::string type = Ioss::Utils::fixup_element_type(element_type, my_node_count[iblk],
							 spatialDimension);
      if (local_element_count[iblk] == 0 && type == "") {
	// For an empty block, exodusII does not store the element
	// type information and returns "NULL" If there are no
	// elements on any processors for this block, it will have
	// an empty type which is invalid and will throw an
	// exception in the ElementBlock constructor. Try to discern
	// the correct element type based on the block_name.
	std::vector<std::string> tokens;
	Ioss::tokenize(block_name, "_", tokens);
	if (tokens.size() >= 2) {
	  // Check whether last token names an element topology type...
	  const Ioss::ElementTopology *topology = Ioss::ElementTopology::factory(tokens[tokens.size()-1], true);
	  if (topology != NULL) {
	    type = topology->name();
	  }
	}	    
      }
      
      if (type == "null" || type == "") {
	// If we have no idea what the topology type for an empty
	// element block is, call it "unknown"
	type = "unknown";
	
	// If there are no elements on any processor for this block and
	// we have no idea what the topology type is, skip it...
	if (global_element_count[iblk] == 0) {
	  continue;
	}
      }
      
      block = new Ioss::ElementBlock(this, block_name, type,
				     local_element_count[iblk],
				     attributes[iblk]);

      block->property_add(Ioss::Property("id", id));
      
      // Maintain block order on output database...
      block->property_add(Ioss::Property("original_block_order", used_blocks++));
      
      if (save_type != "null" && save_type != "") {
	if (block->property_exists("original_element_type")) {
	  if (block->get_property("original_element_type").get_string() != save_type) {
	    block->property_erase("original_element_type");
	    block->property_add(Ioss::Property("original_element_type", save_type));
	  }
	} else {
	  if (block->get_property("topology_type").get_string() != save_type) {
	    // Maintain original element type on output database if possible.
	    block->property_add(Ioss::Property("original_element_type", save_type));
	  }
	}
      }
      
      block->property_add(Ioss::Property("global_entity_count", global_element_count[iblk]));
      
      add_attribute_fields(block, attributes[iblk], type);
      
      offset += local_element_count[iblk];
      
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
      get_region()->add(block);
      if (block_name != alias) {
	get_region()->add_alias(block_name, alias);
      }
    
      // Check for results variables.
      add_results_fields(EX_ELEM_BLOCK, block, local_element_count[iblk], iblk);
      
      {
	Ioss::SerializeIO	serializeIO__(this);
	add_map_fields(get_file_pointer(), block, local_element_count[iblk]);
      }
    }
    elementBlockCount = used_blocks;
    assert(elementCount == offset);
    assert(elementBlockCount == (int)get_region()->get_element_blocks().size());
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
    for (int iblk = 0; iblk < elementBlockCount; iblk++) {

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
    
    if (elementBlockCount == 1) {
      blockAdjacency.resize(elementBlockCount);
      blockAdjacency[0].resize(elementBlockCount);
      blockAdjacency[0][0] = 0;
      return;
    }

    // Each processor first calculates the adjacencies on their own
    // processor...

    std::vector<int> node_used(nodeCount);
    std::vector<std::vector<int> > inv_con(nodeCount);
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    assert(check_block_order(element_blocks));

    for (int iblk = 0; iblk < elementBlockCount; iblk++) {
      Ioss::ElementBlock *eb = element_blocks[iblk];
      int blk_position =  eb->get_property("original_block_order").get_int();
      int id =            eb->get_property("id").get_int();
      int element_nodes = eb->get_property("topology_node_count").get_int();
      int my_element_count = eb->get_property("entity_count").get_int();
      if (my_element_count > 0) {
	std::vector<int> conn(my_element_count * element_nodes);
	ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, &conn[0], NULL, NULL);
	
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
      size_t bits_size = (elementBlockCount+word_size-1)/word_size;
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
      
      result = MPI_Waitall(req_cnt, &request[0], &status[0]);
      
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
	for (int iblk = 0; iblk < elementBlockCount; iblk++) {
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
    blockAdjacency.resize(elementBlockCount);
    for (int iblk = 0; iblk < elementBlockCount; iblk++) {
      blockAdjacency[iblk].resize(elementBlockCount);
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
      size_t bits_size = (elementBlockCount+word_size-1)/word_size;

      std::vector<unsigned> data(elementBlockCount * bits_size);
      int offset = 0;
      for (int jblk = 0; jblk < elementBlockCount; jblk++) {
	for (int iblk = 0; iblk < elementBlockCount; iblk++) {
	  if (blockAdjacency[jblk][iblk] == 1) {
	    size_t wrd_off = iblk / word_size;
	    size_t bit     = iblk % word_size;
	    data[offset+wrd_off] |= (1 << bit);
	  }
	}
	offset += bits_size;
      }

      std::vector<unsigned> out_data(elementBlockCount * bits_size);
      MPI_Allreduce((void*)&data[0], &out_data[0], static_cast<int>(data.size()),
		    MPI_UNSIGNED, MPI_BOR, util().communicator());

      offset = 0;
      for (int jblk = 0; jblk < elementBlockCount; jblk++) {
	for (int iblk = 0; iblk < elementBlockCount; iblk++) {
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
    for (int iblk = 0; iblk < elementBlockCount; iblk++) {
      for (int jblk = iblk; jblk < elementBlockCount; jblk++) {
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

    for (int i=0; i < elementBlockCount; i++) {
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
	ex_get_conn(get_file_pointer(), EX_ELEM_BLOCK, id, &conn[0], NULL, NULL);
	
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

    if (sidesetCount > 0) {
      check_side_topology();

      // Get exodusII sideset metadata

      // Get the names (may not exist) of all sidesets and see if they are actually
      // side "blocks" (perhaps written by IO system for a restart).  In that case,
      // they were split by a previous run and we need to reconstruct the side "set"
      // that may contain one or more of them.
      SideSetMap fs_map;
      SideSetSet fs_set;

      Ioss::IntVector side_set_ids(sidesetCount);
      {
	Ioss::SerializeIO	serializeIO__(this);
	int error = ex_get_ids(get_file_pointer(), EX_SIDE_SET, &side_set_ids[0]);
	if (error < 0) {
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}

	for (int i = 0; i < sidesetCount; i++) {
	  char ss_name[max_string_length+1];
	  error = ex_get_name(get_file_pointer(), EX_SIDE_SET, side_set_ids[i], ss_name);
	  if (error < 0) {
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }
	  if (ss_name[0] != '\0') {
	    Ioss::Utils::fixup_name(ss_name);
	    decode_surface_name(fs_map, fs_set, ss_name);
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

      for (int iss = 0; iss < sidesetCount; iss++) {
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

	  side_set_name = get_entity_name(get_file_pointer(), EX_SIDE_SET, id, "surface");

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
	  
	  int ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &sides[0]);
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

	    for (int i=0; i < elementBlockCount; i++) {
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
	  util().global_array_minmax(&global_side_counts[0], global_side_counts.size(),
				     Ioss::ParallelUtils::DO_MAX);
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
	      int side_count = (*I).second;

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
								side_count);
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
						Ioss::Field::MESH, side_count));
	      }

	      if (side_set_name == "universal_sideset") {
		side_block->field_add(Ioss::Field("side_ids",
						Ioss::Field::INTEGER, "scalar",
						Ioss::Field::MESH, side_count));
	      }

	      // Add results fields 
	      Ioss::SerializeIO	serializeIO__(this);
	      add_results_fields(EX_SIDE_SET, side_block, side_count, iss);
	    }
	  }
	  ++I;
	}
      }
    }
  }

  void DatabaseIO::compute_block_membership(int id, std::vector<std::string> &block_membership) const
  {
    Ioss::IntVector block_ids(elementBlockCount);
    if (elementBlockCount == 1) {
      block_ids[0] = 1;
    } else {
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
	
	int ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &sides[0]);
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

      // Synchronize among all processors....
      if (isParallel) {
	util().global_array_minmax(&block_ids[0],
					 block_ids.size(),
					 Ioss::ParallelUtils::DO_MAX);
      }
    }
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    assert(check_block_order(element_blocks));

    for (int i=0; i < elementBlockCount; i++) {
      if (block_ids[i] == 1) {
	Ioss::ElementBlock *block = element_blocks[i];
	if (!block_is_omitted(block)) {
	  block_membership.push_back(block->name());
	}
      }
    }
  }
  
  void DatabaseIO::compute_block_membership(Ioss::EntityBlock *efblock,
					    std::vector<std::string> &block_membership) const
  {
    Ioss::IntVector block_ids(elementBlockCount);
    if (elementBlockCount == 1) {
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
      util().global_array_minmax(&block_ids[0],
				 block_ids.size(),
				 Ioss::ParallelUtils::DO_MAX);
    }
    
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    assert(check_block_order(element_blocks));
    
    for (int i=0; i < elementBlockCount; i++) {
      if (block_ids[i] == 1) {
	Ioss::ElementBlock *block = element_blocks[i];
	if (!block_is_omitted(block)) {
	  block_membership.push_back(block->name());
	}
      }
    }
  }
  
  void DatabaseIO::get_nodesets()
  {
    // Attributes of a nodeset are:
    // -- id
    // -- name
    // -- number of nodes
    // -- number of distribution factors (see next comment)
    // ----the #distribution factors should equal #nodes or 0, any
    //     other value does not make sense. If it is 0, then a substitute
    //     list will be created returning 1.0 for the factor

    // In a parallel execution, it is possible that a nodeset will have
    // no nodes or distribution factors on a particular processor...

    // Get exodusII nodeset metadata
    if (nodesetCount > 0) {
      Ioss::IntVector nodeset_ids(nodesetCount);
      int error = ex_get_ids(get_file_pointer(), EX_NODE_SET, &nodeset_ids[0]);
      if (error < 0) {
	exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }

      for (int ins = 0; ins < nodesetCount; ins++) {
	int id = nodeset_ids[ins];
	int number_nodes;
	int number_distribution_factors;

	error = ex_get_set_param(get_file_pointer(), EX_NODE_SET, id,
				 &number_nodes, &number_distribution_factors);
	if (error < 0) {
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	}

	std::string nodeset_name = get_entity_name(get_file_pointer(), EX_NODE_SET, id, "nodelist");
	Ioss::NodeSet *nodeset = new Ioss::NodeSet(this, nodeset_name, number_nodes);
	nodeset->property_add(Ioss::Property("id", id));
	get_region()->add(nodeset);

	add_results_fields(EX_NODE_SET, nodeset, number_nodes, ins);

	get_region()->add_alias(nodeset_name, Ioss::Utils::encode_entity_name("nodelist", id));
	get_region()->add_alias(nodeset_name, Ioss::Utils::encode_entity_name("nodeset", id));
      }

    }
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
	  nodeCmapIds      = new int[commsetNodeCount];
	  nodeCmapNodeCnts = new int[commsetNodeCount];
	}
	if (commsetElemCount > 0) {
	  elemCmapIds      = new int[commsetElemCount];
	  elemCmapElemCnts = new int[commsetElemCount];
	}

	int error = ne_get_cmap_params(get_file_pointer(),
				       nodeCmapIds, nodeCmapNodeCnts,
				       elemCmapIds, elemCmapElemCnts,
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
	    std::vector<double> y(num_to_get);
	    std::vector<double> z;
	    if (spatialDimension == 3)
	      z.resize(num_to_get);

	    // Cast 'data' to correct size -- double
	    double *rdata = static_cast<double*>(data);

	    int ierr = ex_get_coord(get_file_pointer(), &x[0], &y[0], &z[0]);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

	    int index = 0;
	    for (size_t i=0; i < num_to_get; i++) {
	      rdata[index++] = x[i];
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
	  else if (field.get_name() == "node_connectivity_status") {
	    compute_node_status();
	    char *status = static_cast<char*>(data);
	    std::copy(nodeConnectivityStatus.begin(), nodeConnectivityStatus.end(), &status[0]);
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
	  num_to_get = read_transient_field(EX_NODAL, field, nb, data);
	}
      }
      return num_to_get;
    }
  }

  int DatabaseIO::read_transient_field(ex_entity_type type,
				       const Ioss::Field& field,
				       const Ioss::GroupingEntity *ge,
				       void *variables) const
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
      if (type == EX_NODAL) {
	var_index = nodalVariables.find(var_name)->second;
      } else if (type == EX_ELEM_BLOCK) {
	var_index = elementVariables.find(var_name)->second;
      } else if (type == EX_NODE_SET) {
	var_index = nodesetVariables.find(var_name)->second;
      } else if (type == EX_SIDE_SET) {
	var_index = sidesetVariables.find(var_name)->second;
      }
      assert(var_index > 0);
      ierr = ex_get_var(get_file_pointer(), step, type,
			var_index, id, num_entity, &temp[0]);
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // Transfer to 'variables' array.
      int k = 0;
      if (field.get_type() == Ioss::Field::INTEGER) {
	int *ivar = static_cast<int*>(variables);
	for (int j=i; j < num_entity*comp_count; j+=comp_count) {
	  ivar[j] = static_cast<int>(temp[k++]);
	}
      } else if (field.get_type() == Ioss::Field::REAL) {
	double *rvar = static_cast<double*>(variables);
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
    size_t side_count = is_valid_side.size();
    std::vector<double> temp(side_count);

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
      var_index = sidesetVariables.find(var_name)->second;
      assert(var_index > 0);
      ierr = ex_get_var(get_file_pointer(), step, EX_SIDE_SET,
			var_index, id, side_count, &temp[0]);
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // Transfer to 'variables' array.
      int j = i;
      if (field.get_type() == Ioss::Field::INTEGER) {
	int *ivar = static_cast<int*>(variables);
	for (size_t k = 0; k < side_count; k++) {
	  if (is_valid_side[k] == 1) {
	    ivar[j] = static_cast<int>(temp[k]);
	    j += comp_count;
	  }
	}
      } else if (field.get_type() == Ioss::Field::REAL) {
	double *rvar = static_cast<double*>(variables);
	for (size_t k = 0; k < side_count; k++) {
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
				    &element[0]);
	    ex_get_partial_elem_map(get_file_pointer(), 2, eb_offset+1, my_element_count,
				    &side[0]);

	    int index = 0;
	    for (int i=0; i < my_element_count; i++) {
	      el_side[index++] = element[i];
	      el_side[index++] = side[i];
	    }

	  }
	  else {
	    num_to_get = Ioss::Utils::field_warning(eb, field, "input");
	  }
	}

	else if (role == Ioss::Field::ATTRIBUTE) {
	  if (my_element_count > 0) {
	    int attribute_count = eb->get_property("attribute_count").get_int();

	    std::string att_name = eb->name() + SEP() + field.get_name();
	    int offset = (int)field.get_index();
	    assert(offset > 0);
	    assert(offset-1+field.raw_storage()->component_count() <= attribute_count);
	    if (offset == 1 && field.raw_storage()->component_count() == attribute_count) {
	      // Read all attributes in one big chunk...
	      ierr = ex_get_attr(get_file_pointer(), EX_ELEM_BLOCK, id, static_cast<double*>(data));
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);
	    }
	    else {
	      // Read a subset of the attributes.  If scalar, read one;
	      // if higher-order (vector3d, ..) read each component and
	      // put into correct location...
	      if (field.raw_storage()->component_count() == 1) {
		ierr = ex_get_one_attr(get_file_pointer(), EX_ELEM_BLOCK, id,
				       offset, static_cast<double*>(data));
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      } else {
		// Multi-component...
		// Need a local memory space to read data into and
		// then push that into the user-supplied data block...
		std::vector<double> local_data(my_element_count);
		int comp_count = field.raw_storage()->component_count();
		double *rdata = static_cast<double*>(data);
		for (int i=0; i < comp_count; i++) {
		  ierr = ex_get_one_attr(get_file_pointer(), EX_ELEM_BLOCK, id,
					 offset+i, &local_data[0]);
		  if (ierr < 0)
		    exodus_error(get_file_pointer(), __LINE__, myProcessor);

		  int k = i;
		  for (int j=0; j < my_element_count; j++) {
		    rdata[k] = local_data[j];
		    k += comp_count;
		  }
		}
	      }
	    }
	  }
	} else if (role == Ioss::Field::TRANSIENT) {
	  // Check if the specified field exists on this element block.
	  // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	  // exist on the database as scalars with the appropriate
	  // extensions.

	  // Read in each component of the variable and transfer into
	  // 'data'.  Need temporary storage area of size 'number of
	  // elements in this block.
	  num_to_get = read_transient_field(EX_ELEM_BLOCK, field, eb, data);
	} else if (role == Ioss::Field::REDUCTION) {
	  num_to_get = Ioss::Utils::field_warning(eb, field, "input reduction");
	}
      }
      return num_to_get;
    }
  }

  int DatabaseIO::get_field_internal(const Ioss::NodeSet* ns,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    {
      Ioss::SerializeIO	serializeIO__(this);

      size_t num_to_get = field.verify(data_size);
      if (num_to_get > 0) {

	int id = get_id(ns, EX_NODE_SET, &ids_);
	Ioss::Field::RoleType role = field.get_role();
	if (role == Ioss::Field::MESH) {

	  if (field.get_name() == "ids") {
	    int ierr = ex_get_set(get_file_pointer(), EX_NODE_SET, id, static_cast<int*>(data), NULL);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

	    // Convert the local node ids to global ids
	    const Ioss::MapContainer &map = get_node_map();
	    if (!Ioss::Map::is_sequential(map)) {
	      int *ids = static_cast<int*>(data);
	      for (size_t i=0; i < num_to_get; i++) {
		ids[i] = map[ids[i]];
	      }
	    }
	  } else if (field.get_name() == "distribution_factors") {
	    int number_nodes;
	    int number_distribution_factors;

	    int ierr = ex_get_set_param(get_file_pointer(), EX_NODE_SET, id,
					&number_nodes, &number_distribution_factors);
	    if (number_distribution_factors == 0) {
	      double *rdata = static_cast<double*>(data);
	      for (size_t i=0; i < num_to_get; i++)
		rdata[i] = 1.0;
	    } else {
	      ierr = ex_get_set_dist_fact(get_file_pointer(), EX_NODE_SET, id, static_cast<double*>(data));
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);
	    }
	  } else {
	    num_to_get = Ioss::Utils::field_warning(ns, field, "input");
	  }
	} else if (role == Ioss::Field::TRANSIENT) {
	  // Check if the specified field exists on this node block.
	  // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	  // exist on the database as scalars with the appropriate
	  // extensions.

	  // Read in each component of the variable and transfer into
	  // 'data'.  Need temporary storage area of size 'number of
	  // nodes in this block.
	  num_to_get = read_transient_field(EX_NODE_SET, field, ns, data);
	}
      }
      return num_to_get;
    }
  }

  int DatabaseIO::get_field_internal(const Ioss::SideSet* fs,
				     const Ioss::Field& field,
				     void */* data */, size_t /* data_size */) const
  {
    size_t num_to_get = Ioss::Utils::field_warning(fs, field, "input");
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
	    ierr = ex_get_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id, &real_ids[0]);
	    
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

	  ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, element, &sides[0]);
	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);

	  if (number_sides == static_cast<int>(entity_count)) {
	    for (size_t iel = 0; iel < entity_count; iel++) {
	      int new_id = 10*map[element[iel]] + sides[iel];
	      ids[iel] = new_id;
	    }
	  } else {
	    Ioss::IntVector is_valid_side(number_sides);
	    Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, element, &sides[0],
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

	  ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &sides[0]);
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
	    Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, &element[0], &sides[0],
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

	  ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &sides[0]);
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
	    Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, &element[0], &sides[0],
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
				       static_cast<int*>(data), data_size/sizeof(int));
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
	  num_to_get = read_transient_field(EX_SIDE_SET, field, fb, data);
	} else {
	  // Need to read all values for the specified field and then
	  // filter down to the elements actualy in this side block.

	  // Determine which sides are member of this block
	  Ioss::IntVector is_valid_side(number_sides);
	  {
//----
	    Ioss::IntVector element(number_sides);
	    Ioss::IntVector sides(number_sides);
	    ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &sides[0]);
//----
	    Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb,
							&element[0], &sides[0],
							number_sides, get_region());
	  }
	  
	  num_to_get = read_ss_transient_field(field, id, data, is_valid_side);
	}
      }
    }
    return num_to_get;
  }

  int DatabaseIO::get_side_connectivity(const Ioss::EntityBlock* fb,
					int id, int,
					int *fconnect,
					size_t /* data_size */) const
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

    ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &side[0]);
    if (ierr < 0)
      exodus_error(get_file_pointer(), __LINE__, myProcessor);
//----

    Ioss::IntVector is_valid_side(number_sides);
    Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, &element[0], &side[0],
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
	  get_field_internal(block, block->get_field("connectivity"),
			     &elconnect[0], nelem*nelnode*sizeof(int));
	  conn_block = block;
	  current_side = -1;
	}

	// NOTE: Element connectivity is returned with nodes in global id space
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
  int DatabaseIO::get_side_distributions(const Ioss::EntityBlock* fb,
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
    ierr = ex_get_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id, &dist[0]);
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

    ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &side[0]);
    if (ierr < 0)
      exodus_error(get_file_pointer(), __LINE__, myProcessor);
//----

    Ioss::IntVector is_valid_side(number_sides);
    Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, &element[0], &side[0],
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
	    std::vector<double> y(num_to_get);
	    std::vector<double> z;
	    if (spatialDimension == 3)
	      z.resize(num_to_get);

	    // Cast 'data' to correct size -- double
	    double *rdata = static_cast<double*>(data);

	    int index = 0;
	    for (size_t i=0; i < num_to_get; i++) {
	      x[i] = rdata[index++];
	      y[i] = rdata[index++];
	      if (spatialDimension == 3)
		z[i] = rdata[index++];
	    }
	    int ierr = ex_put_coord(get_file_pointer(), &x[0], &y[0], &z[0]);
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
	  write_nodal_transient_field(EX_NODAL, field, nb, num_to_get, data);

	} else if (role == Ioss::Field::REDUCTION) {
	  store_reduction_field(EX_NODAL, field, nb, data);
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

	      if (!sequentialNG2L) {
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
					   &element[0]);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	    ierr = ex_put_partial_elem_map(get_file_pointer(), 2, eb_offset+1, my_element_count,
					   &side[0]);
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
	}

	else if (role == Ioss::Field::ATTRIBUTE) {
	  std::string att_name = eb->name() + SEP() + field.get_name();
	  int offset = (int)field.get_index();
	  assert(offset > 0);

	  int attribute_count = eb->get_property("attribute_count").get_int();
	  assert(offset > 0);
	  assert(offset-1+field.raw_storage()->component_count() <= attribute_count);
	  
	  if (offset == 1 && field.raw_storage()->component_count() == attribute_count) {
	    // Write all attributes in one big chunk...
	    ierr = ex_put_attr(get_file_pointer(), EX_ELEM_BLOCK, id, static_cast<double*>(data));
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);
	  }
	  else {
	    // Write a subset of the attributes.  If scalar, write one;
	    // if higher-order (vector3d, ..) write each component.
	    if (field.raw_storage()->component_count() == 1) {
	      ierr = ex_put_one_attr(get_file_pointer(), EX_ELEM_BLOCK, id,
				     offset, static_cast<double*>(data));
	      if (ierr < 0)
		exodus_error(get_file_pointer(), __LINE__, myProcessor);
	    } else {
	      // Multi-component...  Need a local memory space to push
	      // data into and then write that out to the file...
	      std::vector<double> local_data(my_element_count);
	      int comp_count = field.raw_storage()->component_count();
	      double *rdata = static_cast<double*>(data);
	      for (int i=0; i < comp_count; i++) {
		int k = i;
		for (int j=0; j < my_element_count; j++) {
		  local_data[j] = rdata[k];
		  k += comp_count;
		}
		
		ierr = ex_put_one_attr(get_file_pointer(), EX_ELEM_BLOCK, id,
				       offset+i, &local_data[0]);
		if (ierr < 0)
		  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	      }
	    }
	  }
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

    // The ids coming in are the global ids; their position is the
    // local id -1 (That is, data[0] contains the global id of local
    // node 1)

    assert(static_cast<int>(num_to_get) == nodeCount);
    
    if (dbState == Ioss::STATE_MODEL) {
      if (nodeMap.empty()) {
	nodeMap.resize(nodeCount+1);
	sequentialNG2L = true;
	nodeMap[0] = -1;
      }

      if (sequentialNG2L) {
	for (size_t i=0; i < num_to_get; i++) {
	  nodeMap[i+1] = ids[i];
	  if (static_cast<int>(i+1) != ids[i]) {
	    assert(ids[i] != 0);
	    sequentialNG2L = false;
	    nodeMap[0] = 1;
	  }
	}
      }

      if (!sequentialNG2L || static_cast<int>(num_to_get) != nodeCount) {
	sequentialNG2L = false;
	Ioss::Map::build_reverse_map(&reverseNodeMap, ids, num_to_get, 0,
				     "node", myProcessor);
      }

      // Only a single nodeblock and all set
      if (static_cast<int>(num_to_get) == nodeCount) {
	assert(sequentialNG2L || static_cast<int>(reverseNodeMap.size()) == nodeCount);
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

  int DatabaseIO::handle_element_ids(const Ioss::ElementBlock *eb,
				     int* ids, size_t num_to_get)
  {
    /*!
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
    if (elementMap.empty()) {
      elementMap.resize(elementCount+1);
      sequentialEG2L = true;
    }

    assert(static_cast<int>(elementMap.size()) == elementCount+1);
    elementMap[0] = -1;

    int eb_offset = eb->get_offset();

    for (size_t i=0; i < num_to_get; i++) {
      int local_id = eb_offset + i + 1;
      elementMap[local_id] = ids[i];
      if (local_id != ids[i]) {
	sequentialEG2L = false;
	elementMap[0] = 1;
	assert(ids[i] != 0);
      }
    }

    // Now, if the state is Ioss::STATE_MODEL, update the reverseElementMap
    if (dbState == Ioss::STATE_MODEL) {
      Ioss::Map::build_reverse_map(&reverseElementMap, ids, num_to_get,
				   eb_offset, "element", myProcessor);

      // Strongest assertion we can make is that size of map <=
      // elementCount
      assert(static_cast<int>(reverseElementMap.size()) <= elementCount);

      // Output this portion of the element number map
      int ierr = ne_put_n_elem_num_map(get_file_pointer(), eb_offset+1, num_to_get, ids);
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);
    }
    // Build the reorderElementMap which does a direct mapping from
    // the current topologies local order to the local order stored in
    // the database...  This is 0-based and used for remapping output
    // TRANSIENT fields. (Will also need one on input once we read fields)
    build_element_reorder_map(eb_offset, num_to_get);
    return num_to_get;
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
    // of this name in 'nodalVariables' map
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

	if (nodalVariables.find(var_name) == nodalVariables.end()) {
	  std::ostringstream errmsg;
	  errmsg << "ERROR: Could not find nodal variable '" << var_name << "'\n";
	  IOSS_ERROR(errmsg);
	}

	var_index = nodalVariables.find(var_name)->second;

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
	int ierr = ex_put_var(get_file_pointer(), step, EX_NODAL, var_index, 0, num_out, &temp[0]);
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
    assert(type == EX_ELEM_BLOCK || type == EX_NODE_SET || type == EX_SIDE_SET);

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
    // of this name in 'elementVariables' map
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

	if (type == EX_ELEM_BLOCK)
	  var_index = elementVariables.find(var_name)->second;
	else if (type == EX_NODE_SET)
	  var_index = nodesetVariables.find(var_name)->second;
	else if (type == EX_SIDE_SET)
	  var_index = sidesetVariables.find(var_name)->second;

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
	if (!surfaceSplitBackwardCompatibility && type == EX_SIDE_SET) {
	  int offset = ge->get_property("set_offset").get_int();
	  ierr = ex_put_n_var(get_file_pointer(), step, type, var_index, id, offset+1, count, &temp[0]);
	} else {
	  ierr = ex_put_var(get_file_pointer(), step, type, var_index, id, count, &temp[0]);
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
    // of this name in 'globalVariables' map
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
	assert(globalVariables.find(var_name) != globalVariables.end());
	var_index = globalVariables.find(var_name)->second;

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
    // of this name in 'globalVariables' map
    int comp_count = var_type->component_count();
    for (int i=0; i < comp_count; i++) {
      std::string field_name = field.get_name();
      std::string var_name = var_type->label_name(field_name, i+1, fieldSuffixSeparator);

      assert(globalVariables.find(var_name) != globalVariables.end());
      int var_index = globalVariables.find(var_name)->second;

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
				  (double*)&globalValues[0]);
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
				  &globalValues[0]);
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);
    }
  }

  int DatabaseIO::put_field_internal(const Ioss::NodeSet* ns,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    {
      ex_update(get_file_pointer());
      Ioss::SerializeIO	serializeIO__(this);

      int entity_count = ns->get_property("entity_count").get_int();
      size_t num_to_get = field.verify(data_size);
      if (num_to_get > 0) {

	int id = get_id(ns, EX_NODE_SET, &ids_);
	Ioss::Field::RoleType role = field.get_role();

	if (role == Ioss::Field::MESH) {

	  if (field.get_name() == "ids") {
	    // Map node id from global node id to local node id.
	    // Do it in 'data' ...
	    int* ids = static_cast<int*>(data);

	    if (!sequentialNG2L) {
	      for (size_t i=0; i < num_to_get; i++) {
		int global_id = ids[i];
		ids[i] = node_global_to_local(global_id, true);
	      }
	    }
	    int ierr = ex_put_set(get_file_pointer(), EX_NODE_SET, id, static_cast<int*>(data), NULL);
	    if (ierr < 0)
	      exodus_error(get_file_pointer(), __LINE__, myProcessor);

	  } else if (field.get_name() == "distribution_factors") {
	    int ierr = ex_put_set_dist_fact(get_file_pointer(), EX_NODE_SET, id, static_cast<double*>(data));
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
	  write_entity_transient_field(EX_NODE_SET, field, ns, entity_count, data);

	} else if (role == Ioss::Field::REDUCTION) {
	  store_reduction_field(EX_NODE_SET, field, ns, data);
	}
      }
      return num_to_get;
    }
  }

  int DatabaseIO::put_field_internal(const Ioss::SideSet* fs,
				     const Ioss::Field& field,
				     void */* data */, size_t /* data_size */) const
  {
    size_t num_to_get = Ioss::Utils::field_warning(fs, field, "output");
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
				      &entities[0], &procs[0], myProcessor);
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

	  int ierr = ne_put_node_map(get_file_pointer(), &internal[0], &entities[0], NULL,
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
				    &entities[0], &sides[0], &procs[0], myProcessor);
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

	ierr = ne_put_elem_map(get_file_pointer(), &internal[0],
			       &entities[0], myProcessor);
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
      if (!surfaceSplitBackwardCompatibility)
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
	  //	  int ierr = ex_put_partial_set_dist_fact(get_file_pointer(),  EX_SIDE_SET, id, &real_ids[0]);
	  int ierr = ne_put_n_side_set_df(get_file_pointer(), id, offset+1, entity_count, &real_ids[0]);
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
	  if (!surfaceSplitBackwardCompatibility) {
	    int df_offset = fb->get_property("set_df_offset").get_int();
	    int df_count  = fb->get_property("distribution_factor_count").get_int();
	    ierr = ne_put_n_side_set_df(get_file_pointer(), id, df_offset+1, df_count, static_cast<double*>(data));
	  } else {
	    ierr = ex_put_set_dist_fact(get_file_pointer(), EX_SIDE_SET, id, static_cast<double*>(data));
	  }
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
	  if (!surfaceSplitBackwardCompatibility) {
	    ierr = ne_put_n_side_set(get_file_pointer(), id, offset+1, entity_count, &element[0], &side[0]);
	  } else {
	    ierr = ex_put_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &side[0]);
	  }
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
	  if (!surfaceSplitBackwardCompatibility) {
	    ierr = ne_put_n_side_set(get_file_pointer(), id, offset+1, entity_count, &element[0], &side[0]);
	  } else {
	    ierr = ex_put_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &side[0]);
	  }
	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);

	} else if (field.get_name() == "connectivity") {
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

    // Get count of nodes, element blocks, nodesets, sidesets, ...
    std::vector<Ioex::Block>   blocks;
    std::vector<Ioex::NodeSet> nsets;
    std::vector<Ioex::SideSet> sidesets;

    char the_title[max_line_length+1];

    // Title...
    if (region->property_exists("title")) {
      std::string title_str = region->get_property("title").get_string();
      std::strncpy(the_title, title_str.c_str(), max_line_length);
    } else {
      std::strncpy(the_title, "Sierra Output Default Title", max_line_length);
    }
    the_title[max_line_length] = '\0';

    {
      Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();
      assert(node_blocks.size() == 1);
      spatialDimension = node_blocks[0]->
	get_property("component_degree").get_int();
      nodeCount =        node_blocks[0]->
	get_property("entity_count").get_int();
    }
    
    // Element Blocks --
    {
      Ioss::ElementBlockContainer element_blocks =
	region->get_element_blocks();
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
	Ioex::Block T(*(*I));
	if (std::find(blocks.begin(), blocks.end(), T) == blocks.end()) {
	  blocks.push_back(T);
	}
      }
      elementBlockCount = blocks.size();
    }
    
    // SideSets ...
    if (surfaceSplitBackwardCompatibility) {
      Ioss::SideSetContainer ssets = region->get_sidesets();
	set_sideblock_ids(ssets, sidesets, util(), &ids_);
    } else {
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
	if (std::find(sidesets.begin(), sidesets.end(), T) == sidesets.end()) {
	  sidesets.push_back(T);
	}
      }
    }
    sidesetCount = sidesets.size();

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
	if (std::find(nsets.begin(), nsets.end(), T) == nsets.end()) {
	  nsets.push_back(T);
	}
      }
      nodesetCount = nsets.size();
    }

    Ioex::Mesh mesh(spatialDimension, nodeCount, elementCount,
		    elementBlockCount, nodesetCount, sidesetCount,
		    the_title);

    Ioex::CommunicationMetaData comm_meta;
    gather_communication_metadata(&comm_meta);

    {
      Ioss::SerializeIO	serializeIO__(this);

      if (myProcessor == 0) {
	put_qa();
	put_info();
      }
      
      // Write the metadata to the exodusII file...
      Ioex::Internals data(get_file_pointer());
      int ierr = data.write_meta_data(mesh, blocks, nsets, sidesets, comm_meta);

      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      // Write element block attribute names (if any)...
      write_attribute_names();

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

  void DatabaseIO::write_attribute_names()
  {
    // For each block, determine the attribute fields and the correct
    // order. Write the names of these fields.  However, be aware that
    // the field "attribute" always exists to contain all attributes
    // and its name should not be used even if it is the only
    // attribute field.
    Ioss::ElementBlockContainer::const_iterator I;
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
      Ioss::ElementBlock *block = *I;
      std::string block_name = block->name();
      int attribute_count = block->get_property("attribute_count").get_int();
      if (attribute_count > 0) {
	
	check_attribute_index_order(block);
	
	char **names = get_exodus_names(attribute_count);

	// Get the attribute fields...
	Ioss::NameList results_fields;
	block->field_describe(Ioss::Field::ATTRIBUTE, &results_fields);

	Ioss::NameList::const_iterator IF;
	for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
	  std::string field_name = *IF;
	  const Ioss::Field &field = block->get_fieldref(field_name);
	  assert(field.get_index() != 0);
	  
	  if (field_name == "attribute") {
	    field.set_index(1);
	    continue;
	  }

	  const Ioss::VariableType *type = field.raw_storage();
	  int comp_count = type->component_count();
	  int field_offset = field.get_index();
	  for (int i=0; i < comp_count; i++) {
	    std::string var_name = type->label_name(field_name, i+1, fieldSuffixSeparator);
	    std::strncpy(names[field_offset-1+i], var_name.c_str(), max_string_length);
	    names[field_offset-1+i][max_string_length] = '\0';
	  }
	}
	int block_id = block->get_property("id").get_int();
	int ierr = ex_put_attr_names(get_file_pointer(), EX_ELEM_BLOCK, block_id, names);
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
	delete_exodus_names(names, attribute_count);
      }
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
    int field_count = add_results_fields(EX_GLOBAL, get_region(), 1);
    globalValues.resize(field_count);
  }

  int DatabaseIO::add_results_fields(ex_entity_type type,
				     Ioss::GroupingEntity *entity,
				     int count, int position)
  {
    int nvar = 0;
    {
      Ioss::SerializeIO	serializeIO__(this);

      int ierr = ex_get_variable_param(get_file_pointer(), type, &nvar);
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);
    }

    if (nvar > 0) {
      int *truth_table = NULL;
      int block_count = 0;
      // If this is for an element block, we need to make sure we have
      // the element truth table first so we can tell whether the
      // variable exists for this element block...
      if (type == EX_ELEM_BLOCK && elementTruthTable == NULL) {
	elementTruthTable = new int[elementBlockCount*nvar];
	truth_table = elementTruthTable;
	block_count = elementBlockCount;
      }
      else if (type == EX_NODE_SET && nodesetTruthTable == NULL) {
	nodesetTruthTable = new int[nodesetCount*nvar];
	truth_table = nodesetTruthTable;
	block_count = nodesetCount;
      }
      else if (type == EX_SIDE_SET && sidesetTruthTable == NULL) {
	sidesetTruthTable = new int[sidesetCount*nvar];
	truth_table = sidesetTruthTable;
	block_count = sidesetCount;
      }

      // Read and store the truth table (Should be there since we only
      // get to this routine if there are variables...)
      if (truth_table != NULL) {
	  Ioss::SerializeIO	serializeIO__(this);

	  int ierr = ex_get_truth_table(get_file_pointer(), type, block_count, nvar, truth_table);
	  if (ierr < 0)
	    exodus_error(get_file_pointer(), __LINE__, myProcessor);

	  // If parallel, then synchronize the truth table among all processors...
	  // Need to know that block_X has variable_Y even if block_X is
	  // empty on a specific processor...  The truth table contains 0
	  // if the variable doesn't exist and 1 if it does, so we just
	  // take the maximum at each location...
	  // This is a collective call...
	  if (isParallel) {
	    util().global_array_minmax(truth_table, block_count * nvar,
					     Ioss::ParallelUtils::DO_MAX);
	  }
      }

      // Get the variable names and add as fields. Need to decode these
      // into vector/tensor/... eventually, for now store all as
      // scalars.

      
      char **names = get_exodus_names(nvar);

      // Read the names...
      // (Currently, names are read for every block.  We could save them...)
      {
	Ioss::SerializeIO	serializeIO__(this);

	int ierr = ex_get_variable_names(get_file_pointer(), type, nvar, names);
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);

	// Convert to lowercase.
	for (int i=0; i < nvar; i++) {
	  Ioss::Utils::fixup_name(names[i]);
	}

	// Add to VariableNameMap elementVariables, nodalVariables,
	// or globalVariables so can determine exodusII index given a
	// Sierra field name.  exodusII index is just 'i+1'
	for (int i=0; i < nvar; i++) {
	  if (type == EX_ELEM_BLOCK) {
	    elementVariables.insert(VNMValuePair(std::string(names[i]), i+1));
	  } else if (type == EX_NODAL) {
	    nodalVariables.insert(VNMValuePair(std::string(names[i]),   i+1));
	  } else if (type == EX_NODE_SET) {
	    nodesetVariables.insert(VNMValuePair(std::string(names[i]), i+1));
	  } else if (type == EX_SIDE_SET) {
	    sidesetVariables.insert(VNMValuePair(std::string(names[i]), i+1));
	  } else if (type == EX_GLOBAL) {
	    globalVariables.insert(VNMValuePair(std::string(names[i]),  i+1));
	  }
	}

	// Add as fields...
	if (type == EX_ELEM_BLOCK) truth_table = elementTruthTable;
	if (type == EX_NODE_SET)   truth_table = nodesetTruthTable;
	if (type == EX_SIDE_SET)   truth_table = sidesetTruthTable;

	int offset = position*nvar;
	int *local_truth = NULL;
	if (truth_table != NULL)
	  local_truth = &truth_table[offset];

	std::vector<Ioss::Field> fields;
	get_fields(count, names, nvar, Ioss::Field::TRANSIENT,
		   fieldSuffixSeparator, local_truth, fields);

	std::vector<Ioss::Field>::const_iterator IF;
	for (IF = fields.begin(); IF != fields.end(); ++IF) {
	  Ioss::Field field = *IF;
	  entity->field_add(field);
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
    // Does not yet support results on sideblocks or any sets
    int glob_index = 0;
    glob_index = gather_names(EX_GLOBAL, get_region(), glob_index, true);
    assert(glob_index == static_cast<int>(globalVariables.size()));

    {
      Ioss::NodeBlockContainer node_blocks = get_region()->get_node_blocks();
      assert(node_blocks.size() == 1);
      Ioss::NodeBlockContainer::const_iterator I;

      int index = 0;
      for (I=node_blocks.begin(); I != node_blocks.end(); ++I) {
	glob_index = gather_names(EX_NODAL, *I, glob_index, true);
	index = gather_names(EX_NODAL, *I, index, false);
      }
      assert(index == static_cast<int>(nodalVariables.size()));
    }

    {
      Ioss::ElementBlockContainer element_blocks =
	get_region()->get_element_blocks();
      assert(check_block_order(element_blocks));
      Ioss::ElementBlockContainer::const_iterator I;

      int index = 0;
      for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
	glob_index = gather_names(EX_ELEM_BLOCK, *I, glob_index, true);
	index = gather_names(EX_ELEM_BLOCK, *I, index, false);
      }
      assert(index == static_cast<int>(elementVariables.size()));
      generate_truth_table(EX_ELEM_BLOCK);
    }

    {
      Ioss::NodeSetContainer nodesets = get_region()->get_nodesets();
      Ioss::NodeSetContainer::const_iterator I;

      int index = 0;
      for (I=nodesets.begin(); I != nodesets.end(); ++I) {
	glob_index = gather_names(EX_NODE_SET, *I, glob_index, true);
	index = gather_names(EX_NODE_SET, *I, index, false);
      }
      assert(index == static_cast<int>(nodesetVariables.size()));
      generate_truth_table(EX_NODE_SET);
    }

    {
      int index = 0;
	Ioss::SideSetContainer sidesets = get_region()->get_sidesets();
	Ioss::SideSetContainer::const_iterator I;
	for (I=sidesets.begin(); I != sidesets.end(); ++I) {
	  Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
	  Ioss::SideBlockContainer::const_iterator J;

	  for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	    glob_index = gather_names(EX_SIDE_SET, *J, glob_index, true);
	    index = gather_names(EX_SIDE_SET, *J, index, false);
	  }
	}
      assert(index == static_cast<int>(sidesetVariables.size()));
      generate_truth_table(EX_SIDE_SET);
    }

    {
      if (isParallel) {
	// Check consistency among all processors.  They should all
	// have the same number of each variable type...
	// The called function will thrown an exception if the counts differ.
	check_variable_consistency(globalVariables.size(),
				   nodalVariables.size(),
				   elementVariables.size(),
				   nodesetVariables.size(),
				   sidesetVariables.size(),
				   myProcessor, get_filename(), util());
      }

    }

    {
      Ioss::SerializeIO	serializeIO__(this);

      // Defined the results variable metadata...
      // Need special handling in the case where there are no nodes...
      // This should be handled in exodusII API, but currently is not
      size_t node_var_size = nodalVariables.size();
      if (nodeCount == 0)
	node_var_size = 0;
      int ierr = ex_put_all_var_param(get_file_pointer(),
				      globalVariables.size(),
				      node_var_size,
				      elementVariables.size(), elementTruthTable,
				      nodesetVariables.size(), nodesetTruthTable,
				      sidesetVariables.size(), sidesetTruthTable);
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      globalValues.resize(globalVariables.size());
      output_results_names(EX_GLOBAL, globalVariables);
      if (nodeCount > 0)
	output_results_names(EX_NODAL, nodalVariables);
      output_results_names(EX_ELEM_BLOCK, elementVariables);
      output_results_names(EX_NODE_SET,   nodesetVariables);
      output_results_names(EX_SIDE_SET,   sidesetVariables);
    }
  }

  int DatabaseIO::gather_names(ex_entity_type type,
			       const Ioss::GroupingEntity *ge,
			       int index, bool reduction)
  {
    int new_index = index;

    bool nblock = (type == EX_NODAL);

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

	  // Add to 'VariableNameMap elementVariables' or
	  // 'VariableNameMap nodalVariables' so can determine exodusII
	  // index given a Sierra field name.  exodusII index is just 'i+1'
	  if (reduction || type == EX_GLOBAL) {
	    // If this is not a global (region) variable, need to prepend the block name
	    // to avoid name collisions...
	    if (type != EX_GLOBAL) {
	      var_string = ge->name() + ":" + var_string;
	    }
	    if (globalVariables.find(var_string) == globalVariables.end()) {
	      globalVariables.insert(VNMValuePair(var_string, ++new_index));
	    }
	  } else if (type == EX_NODAL) {
	    if (nodalVariables.find(var_string) == nodalVariables.end()) {
	      nodalVariables.insert(VNMValuePair(var_string, ++new_index));
	    }
	  } else if (type == EX_ELEM_BLOCK) {
	    if (elementVariables.find(var_string) == elementVariables.end()) {
	      elementVariables.insert(VNMValuePair(var_string, ++new_index));
	    }
	  } else if (type == EX_NODE_SET) {
	    if (nodesetVariables.find(var_string) == nodesetVariables.end()) {
	      nodesetVariables.insert(VNMValuePair(var_string, ++new_index));
	    }
	  } else if (type == EX_SIDE_SET) {
	    if (sidesetVariables.find(var_string) == sidesetVariables.end()) {
	      sidesetVariables.insert(VNMValuePair(var_string, ++new_index));
	    }
	  }
	}
      }
      if (has_disp && name == disp_name) {
	new_index = save_index;
      }
    }
    return new_index;
  }

  void DatabaseIO::generate_truth_table(ex_entity_type type)
  {
    if (type == EX_ELEM_BLOCK)
      generate_element_truth_table();
    else if (type == EX_NODE_SET)
      generate_nodeset_truth_table();
    else if (type == EX_SIDE_SET)
      generate_sideset_truth_table();
  }

  void DatabaseIO::generate_element_truth_table()
  {
    size_t var_count = elementVariables.size();

    if (var_count == 0 || elementBlockCount == 0)
      return;

    // Member variable.  Will be deleted in destructor...
    elementTruthTable = new int[elementBlockCount*var_count];

    // Zero-fill truth table; set to '1' when variable found...
    std::fill_n(&elementTruthTable[0], elementBlockCount*var_count, 0);

    // Fill in the truth table.  It is conceptually a two-dimensional array of
    // the form 'array[num_element_blocks][num_element_var]'.  In C++,
    // the values for the first element block are first, followed by
    // next element block, ...
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator I;

    assert(check_block_order(element_blocks));

    size_t offset = 0;
    for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
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
	    std::string var_string = var_type->label_name(field_name, i, fieldSuffixSeparator);
	    // Find position of 'var_string' in 'elementVariables'
	    VariableNameMap::iterator VN = elementVariables.find(var_string);
	    if (VN != elementVariables.end()) {
	      // Index '(*VN).second' is 1-based...
	      elementTruthTable[offset + (*VN).second-1] = 1;
	    }
	  }
	}
      }
      offset += var_count;
    }
    assert(offset == var_count * elementBlockCount);
  }

  void DatabaseIO::generate_sideset_truth_table()
  {
    size_t var_count = sidesetVariables.size();

    if (var_count == 0 || sidesetCount == 0)
      return;

    // Member variable.  Will be deleted in destructor...
    sidesetTruthTable = new int[sidesetCount*var_count];

    // Zero-fill truth table; set to '1' when variable found...
    std::fill_n(&sidesetTruthTable[0], sidesetCount*var_count, 0);

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
	      // Find position of 'var_string' in 'elementVariables'
	      VariableNameMap::iterator VN = sidesetVariables.find(var_string);
	      if (VN != sidesetVariables.end()) {
		// Index '(*VN).second' is 1-based...
		sidesetTruthTable[offset + (*VN).second-1] = 1;
	      }
	    }
	  }
	}
	if (found_one && surfaceSplitBackwardCompatibility)
	  offset += var_count;
      }
      if (!surfaceSplitBackwardCompatibility)
	offset += var_count;
    }
    assert(offset == var_count * sidesetCount);
  }

  void DatabaseIO::generate_nodeset_truth_table()
  {
    size_t var_count = nodesetVariables.size();

    if (var_count == 0 || nodesetCount == 0)
      return;

    // Member variable.  Will be deleted in destructor...
    nodesetTruthTable = new int[nodesetCount*var_count];

    // Zero-fill truth table; set to '1' when variable found...
    std::fill_n(&nodesetTruthTable[0], nodesetCount*var_count, 0);

    // Fill in the truth table.  It is conceptually a two-dimensional array of
    // the form 'array[num_blocks][num_var]'.  In C++,
    // the values for the first block are first, followed by
    // next block, ...
    Ioss::NodeSetContainer nodesets = get_region()->get_nodesets();
    Ioss::NodeSetContainer::const_iterator I;

    size_t offset = 0;
    for (I=nodesets.begin(); I != nodesets.end(); ++I) {
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
	    std::string var_string = var_type->label_name(field_name, i, fieldSuffixSeparator);
	    // Find position of 'var_string' in 'elementVariables'
	    VariableNameMap::iterator VN = nodesetVariables.find(var_string);
	    if (VN != nodesetVariables.end()) {
	      // Index '(*VN).second' is 1-based...
	      nodesetTruthTable[offset + (*VN).second-1] = 1;
	    }
	  }
	}
      }
      offset += var_count;
    }
    assert(offset == var_count * nodesetCount);
  }

  void
  DatabaseIO::output_results_names(ex_entity_type type,
				   VariableNameMap &variables) const
  {
    size_t var_count = variables.size();

    if (var_count > 0) {
      // Push into a char** array...
      char **var_names = get_exodus_names(var_count);
      VariableNameMap::const_iterator J  = variables.begin();
      VariableNameMap::const_iterator JE = variables.end();
      while (J != JE) {
	size_t index = (*J).second;
	assert(index > 0 && index <= var_count);
	std::strncpy(var_names[index-1], (*J).first.c_str(), max_string_length);
	var_names[index-1][max_string_length] = '\0';

	if ((*J).first.length() > 32) {
	  if (dbUsage == Ioss::WRITE_RESTART) {
	    std::ostringstream errmsg;
	    errmsg << "INTERNAL ERROR: The variable name for field '" << (*J).first
		      << "' is longer than 32 characters which is not allowed and will cause errors.\n"
		      << "Something is wrong in the Ioex::DatabaseIO class. Please report.\n";
	    IOSS_ERROR(errmsg);
	  } else {
	    if (myProcessor == 0) {
	      IOSS_WARNING << "WARNING: The variable name for field '" << (*J).first
			<< "' is longer than 32 characters which is not allowed and may cause errors.\n"
			<< "The name will be truncated to '" << var_names[index-1]
			<< "' which may result in a non-unique name.\n";
	    }
	  }
	}
	
	++J;
      }

      int ierr = ex_put_variable_names(get_file_pointer(), type, var_count, var_names);
      if (ierr < 0)
	exodus_error(get_file_pointer(), __LINE__, myProcessor);

      delete_exodus_names(var_names, var_count);
    }
  }

  void DatabaseIO::build_element_reorder_map(int start, int count)
  {
    // Note: To further add confusion, the reorderElementMap is 0-based
    // and the reverseElementMap and elementMap are 1-based. This is
    // just a consequence of how they are intended to be used...
    //
    // start is based on a 0-based array -- start of the reorderMap to build.

    if (reorderElementMap.empty())
      reorderElementMap.resize(elementCount);
    assert(static_cast<int>(elementMap.size()) == elementCount+1);
    assert(static_cast<int>(reverseElementMap.size()) <= elementCount); // Built in pieces

    int my_end = start+count;
    for (int i=start; i < my_end; i++) {
      int global_id = elementMap[i+1];
      int orig_local_id = element_global_to_local(global_id) - 1;

      // If we assume that partial output is not being used (it
      // currently isn't in Sierra), then the reordering should only be
      // a permutation of the original ordering within this element block...
      assert(orig_local_id >= start && orig_local_id <= my_end);
      reorderElementMap[i] = orig_local_id;
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
    Ioex::Internals data(get_file_pointer());
    data.update_last_time_attribute(sim_time);

    // Flush the files buffer to disk...
    // If a history file, then only flush if there is more
    // than 90 seconds since the last flush to avoid
    // the flush eating up cpu time for small fast jobs...
    // NOTE: If decide to do this on all files, need to sync across
    // processors to make sure they all flush at same time.

    // GDS: 2011/03/30 -- Use for all non-parallel files, but shorten
    // time for non history files.  Assume that can afford to lose ~90
    // seconds worth of data...  (Flush was taking long time on some
    // /scratch filesystems at SNL for short regression tests with
    // lots of steps)
    bool do_flush = true;
    if (dbUsage == Ioss::WRITE_HISTORY || !isParallel) {
      assert (myProcessor == 0);
      time_t cur_time = time(NULL);
      if (cur_time - timeLastFlush >= 90) {
	timeLastFlush = cur_time;
	do_flush = true;
      } else {
	do_flush = false;
      }
    }
    if (do_flush)
      ex_update(get_file_pointer());
  }

  void Ioex::DatabaseIO::add_attribute_fields(Ioss::ElementBlock *block,
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
      char **names = get_exodus_names(attribute_count);
      int id = block->get_property("id").get_int();
      if (block->get_property("entity_count").get_int() != 0) {
	int ierr = ex_get_attr_names(get_file_pointer(), EX_ELEM_BLOCK, id, &names[0]);
	if (ierr < 0)
	  exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }

      // Sync names across processors...
      std::vector<char> cname(attribute_count * (max_string_length+1));
      if (block->get_property("entity_count").get_int() != 0) {
	for (int i=0; i < attribute_count; i++) {
	  std::memcpy(&cname[i*(max_string_length+1)], names[i], max_string_length+1);
	}
      }
      util().attribute_reduction(attribute_count * (max_string_length+1), &cname[0]);
      for (int i=0; i < attribute_count; i++) {
	std::memcpy(names[i], &cname[i*(max_string_length+1)], max_string_length+1);
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
  static char lc_name[max_string_length+1];
  
  Ioss::NameList::const_iterator IF;
  Ioss::NameList::const_iterator IFend = fields.end();
  size_t max_span = 0;
  
  for (IF = fields.begin(); IF != IFend; ++IF) {
    const char *name = (*IF).c_str();
    std::strncpy(lc_name, name, max_string_length);
    lc_name[max_string_length] = '\0';
    
    Ioss::Utils::fixup_name(lc_name);
    size_t span = match(lc_name, displace);
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
    util.global_array_minmax(&sset_entity_count[0], sset_entity_count.size(),
			     Ioss::ParallelUtils::DO_MAX);
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
  char **names = get_exodus_names(map_count);
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
bool check_block_order(const Ioss::ElementBlockContainer &blocks)
{
  // Verify that element blocks are defined in sorted offset order...
  Ioss::ElementBlockContainer::const_iterator I;

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

  void check_variable_consistency(size_t gv_count, size_t nv_count, size_t ev_count,
				  size_t nsv_count, size_t ssv_count,
				  int my_processor, const std::string& filename,
				  const Ioss::ParallelUtils &util)
  {
#ifdef HAVE_MPI    
    std::vector<int> var_counts(5);
    var_counts[0] = gv_count;
    var_counts[1] = nv_count;
    var_counts[2] = ev_count;
    var_counts[3] = nsv_count;
    var_counts[4] = ssv_count;
    
    Ioss::IntVector all_counts;
    util.gather(var_counts, all_counts);

    bool any_diff = false;
    std::ostringstream errmsg;
    if (my_processor == 0) {
      bool diff[5];
      // See if any differ...
      for (int iv = 0; iv < 5; iv++) {
	diff[iv] = false;
	std::string type;
	switch (iv) {
	case 0:
	  type = "global";
	  break;
	case 1:
	  type = "nodal";
	  break;
	case 2:
	  type = "element";
	  break;
	case 3:
	  type = "nodeset";
	  break;
	case 4:
	  type = "sideset";
	  break;
	}
	
	for (int ip = 1; ip < util.parallel_size(); ip++) {
	  if (var_counts[iv] != all_counts[ip*5+iv]) {
	    any_diff = true;
	    if (!diff[iv]) {
	      Ioss::FileInfo db(filename);
	      diff[iv] = true;
	      errmsg << "\nERROR: Number of " << type 
		     << " variables is not consistent on all processors.\n"
		     << "       Database: " << db.tailname() << "\n"
		     << "\tProcessor 0 count = " << var_counts[iv] << "\n";
	    }
	    errmsg << "\tProcessor " << ip << " count = " << all_counts[ip*5+iv] << "\n";
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

