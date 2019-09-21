// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//
//     * Neither the name of NTESS nor the names of its
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

#ifndef IOSS_Ioss_DatabaseIO_h
#define IOSS_Ioss_DatabaseIO_h

#include <Ioss_BoundingBox.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>    // for DatabaseUsage, etc
#include <Ioss_DataSize.h>   // for DataSize
#include <Ioss_EntityType.h> // for EntityType
#include <Ioss_Map.h>
#include <Ioss_ParallelUtils.h>   // for ParallelUtils
#include <Ioss_PropertyManager.h> // for PropertyManager
#include <Ioss_State.h>           // for State, State::STATE_INVALID
#include <Ioss_SurfaceSplit.h>    // for SurfaceSplitType
#include <chrono>
#include <cstddef> // for size_t, nullptr
#include <cstdint> // for int64_t
#include <map>     // for map
#include <string>  // for string
#include <utility> // for pair
#include <vector>  // for vector

namespace Ioss {
  class CommSet;
  class EdgeBlock;
  class EdgeSet;
  class ElementBlock;
  class ElementSet;
  class ElementTopology;
  class FaceBlock;
  class FaceSet;
  class Field;
  class GroupingEntity;
  class NodeBlock;
  class NodeSet;
  class Region;
  class SideBlock;
  class SideSet;
  class StructuredBlock;
} // namespace Ioss

namespace Ioss {
  class EntityBlock;

  // Contains (parent_element, side) topology pairs
  using TopoContainer = std::vector<std::pair<const ElementTopology *, const ElementTopology *>>;

  /** \brief An input or output Database.
   *
   */
  class DatabaseIO
  {
  public:
    friend class SerializeIO;

    /** \brief Check to see if database state is OK.
     *
     *  \param[in] write_message If true, then output a warning message indicating the problem.
     *  \param[in,out] error_message If non-null on input, then a warning message on output.
     *  \param[in,out] bad_count If non-null on input, then count of the number of processors
     *                 where the file does not exist on output. If ok returns false, but
     * *bad_count==0,
     *                 then the routine does not support this argument.
     *  \returns True if database state is OK. False if not.
     */
    bool ok(bool write_message = false, std::string *error_message = nullptr,
            int *bad_count = nullptr) const
    {
      IOSS_FUNC_ENTER(m_);
      return ok__(write_message, error_message, bad_count);
    }

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    virtual unsigned entity_field_support() const = 0;

    bool using_parallel_io() const { return usingParallelIO; }

    /** \brief Get the local (process-specific) node number corresponding to a global node number.
     *
     *  \param[in] global The global node number
     *  \param[in] must_exist If true, error will occur if the global node number is not
     *             mapped to any local node number on the current process.
     *  \returns The local node number
     */
    int64_t node_global_to_local(int64_t global, bool must_exist) const
    {
      IOSS_FUNC_ENTER(m_);
      return node_global_to_local__(global, must_exist);
    }

    int64_t element_global_to_local(int64_t global) const
    {
      IOSS_FUNC_ENTER(m_);
      return element_global_to_local__(global);
    }

    /** If there is a single block of nodes in the model, then it is
     *  considered a node_major() database.  If instead the nodes are
     * local to each element block or structured block, then it is
     *   not a node_major database.  Exodus is node major, CGNS is not.
     */
    virtual bool node_major() const { return true; }

    virtual ~DatabaseIO();

    // Eliminate as much memory as possible, but still retain meta data information
    // Typically, eliminate the maps...
    void release_memory()
    {
      IOSS_FUNC_ENTER(m_);
      release_memory__();
    }

    // Do anything that might be needed to the database prior to it
    // being closed and destructed.
    virtual void finalize_database() {}

    // Let's save the name on disk after Filename gets modified, e.g: decoded_filename
    void set_pfsname(const std::string &name) const { pfsName = name; }

    std::string get_pfsname() const { return pfsName; }

    /** \brief this will be the name in BB namespace
     */
    void set_dwname(const std::string &name) const { bbName = name; }

    std::string get_dwname() const
    {
      if (!bbName.empty() && !is_input() && using_dw()) {
        return bbName;
      }
      else {
        return get_filename();
      }
    }

    /** \brief We call this ONLY after we assure that using_dw() is TRUE
     *  \ returns mount point of Datawarp namespace, e.g: `/opt/cray/....<jobid>`
     */
    std::string get_dwPath() const { return dwPath; }

    /** Determine whether Cray Datawarp module is loaded and we have BB capacity allocated for this
     * job ( i.e: DW_JOB_STRIPED is set) && IOSS property to use DATAWARP is set to Y/YES (i.e
     * environmental variable ENABLE_DATAWARP) . If we are using DW then set some pathnames for
     * writing directly to BB instead of PFS(i.e Lustre)
     */
    void check_setDW() const;

    /** Set if  Cray Datawarp exists and allocated space is found , i.e PATH to DW name space
     */
    bool using_dw() const { return usingDataWarp; }

    /** \brief Get the file name associated with the database.
     *
     *  \returns The database file name.
     */
    std::string get_filename() const { return DBFilename; }

    /** For the database types that support it, return an integer `handle`
     * through which a client can directly access the underlying file.
     * Please use sparingly and with discretion. Basically, a kluge
     */
    virtual int get_file_pointer() const { return 0; }

    /** \brief Get a file-per-processor filename associated with the database.
     *
     * \ returns The file-per-processor name for a file on this processor.
     */
    const std::string &decoded_filename() const;

    /** Return a string specifying underlying format of database (exodus, cgns, ...) */
    virtual const std::string get_format() const = 0;

    /** \brief Determine whether the database is an input database.
     *
     *  \returns True if the database is an input database. False otherwise.
     */
    bool is_input() const { return isInput; }

    /** \brief Get the Ioss::DatabaseUsage type of the database.
     *
     *  \returns The Ioss::DatabaseUsage type of the database.
     */
    Ioss::DatabaseUsage usage() const { return dbUsage; }

    /** \brief Determine whether the database needs information about process ownership of nodes.
     *
     *  \returns True if database needs information about process ownership of nodes.
     */
    virtual bool needs_shared_node_information() const { return false; }

    Ioss::IfDatabaseExistsBehavior open_create_behavior() const;

    //! This function is used to create the path to an output directory (or history, restart, etc.)
    //  if it does not exist.  Called by all processors. Will throw exception if path does not
    //  specify a valid directory or if the path cannot be created.
    void create_path(const std::string &filename) const;

    void set_region(Region *region) { region_ = region; }

    /** \brief If we are planning to use BB(aka Burst Buffer) service, we will call
     *   simple C API provided by Cray DataWarp module.
     *
     *   Note:  We would only like to use BB during write to avoid cache coherency issues during
     * read. `dwPath` is the DW namespace which gets set during runtime from environment variable
     * DW_JOB_STRIPED. Afterword, `using_dw()` function is used extensively for filename redirection
     * for all related subsequent functions(e.g: get_filename, get_file_ptr etc) once burst buffer
     * is found and set to be used.
     */
    void openDW(const std::string &filename) const;

    /** \brief  Function which invokes stageout  from BB to Disk, prior to completion of final close
     */
    void closeDW() const;

    void openDatabase() const
    {
      IOSS_FUNC_ENTER(m_);
      openDatabase__();
    }

    void closeDatabase() const
    {
      IOSS_FUNC_ENTER(m_);
      closeDatabase__();
    }

    void flush_database() const
    {
      IOSS_FUNC_ENTER(m_);
      flush_database__();
    }

    /** \brief If a database type supports groups and if the database
     *         contains groups, open the specified group.
     *
     *  If the group_name begins with '/', it specifies the absolute path
     *  name from the root with '/' separating groups.  Otherwise, the
     *  group_name specifies a child group of the currently active group.
     *  If group_name == "/" then the root group is opened.
     *
     *  \param[in] group_name The name of the group to open.
     *  \returns True if successful.
     */
    bool open_group(const std::string &group_name)
    {
      IOSS_FUNC_ENTER(m_);
      return open_group__(group_name);
    }

    /** \brief If a database type supports groups, create the specified
     *        group as a child of the current group.
     *
     *  The name of the group must not contain a '/' character.
     *  If the command is successful, then the group will be the
     *  active group for all subsequent writes to the database.
     *
     *  \param[in] group_name The name of the subgroup to create.
     *  \returns True if successful.
     */
    bool create_subgroup(const std::string &group_name)
    {
      IOSS_FUNC_ENTER(m_);
      return create_subgroup__(group_name);
    }

    /** \brief Set the database to the given State.
     *
     *  All transitions must begin from the 'STATE_CLOSED' state or be to
     *  the 'STATE_CLOSED' state (There are no nested begin/end pairs at
     *  this time.)
     *
     *  The database state is automatically set when Region::begin_mode is called
     *  for its associated region, so it may not be necessary to call this method
     *  directly.
     *
     *  \param[in] state The new State to which the database should be set.
     *  \returns True if successful.
     *
     */
    bool begin(Ioss::State state)
    {
      IOSS_FUNC_ENTER(m_);
      return begin__(state);
    }

    /** \brief Return the database to STATE_CLOSED.
     *
     *  The database is automatically set to STATE_CLOSED when Region::end_mode
     *  is called for its associated region, so it may not be necessary to call this
     *  method directly.
     *
     *  \param[in] state The State to end, i.e. the current state.
     *  \returns True if successful.
     *
     */
    bool end(Ioss::State state)
    {
      IOSS_FUNC_ENTER(m_);
      return end__(state);
    }

    bool begin_state(int state, double time);
    bool end_state(int state, double time);

    // Metadata-related functions.
    void read_meta_data()
    {
      IOSS_FUNC_ENTER(m_);
      return read_meta_data__();
    }

    void get_step_times()
    {
      IOSS_FUNC_ENTER(m_);
      return get_step_times__();
    }

    virtual bool internal_edges_available() const { return false; }
    virtual bool internal_faces_available() const { return false; }

    // Information Records:

    /** \brief Get all information records (informative strings) for the database.
     *
     *  \returns The informative strings.
     */
    const std::vector<std::string> &get_information_records() const { return informationRecords; }
    void                            add_information_records(const std::vector<std::string> &info);
    void                            add_information_record(const std::string &info);

    // QA Records:

    /** \brief Get all QA records, each of which consists of 4 strings, from the database
     *
     *  The 4 strings that make up a database QA record are:
     *
     *  1. A descriptive code name, such as the application that modified the database.
     *
     *  2. A descriptive string, such as the version of the application that modified the database.
     *
     *  3. A relevant date, such as the date the database was modified.
     *
     *  4. A relevant time, such as the time the database was modified.
     *
     *  \returns All QA records in a single vector. Every 4 consecutive elements of the
     *           vector make up a single QA record.
     */
    const std::vector<std::string> &get_qa_records() const { return qaRecords; }
    void add_qa_record(const std::string &code, const std::string &code_qa, const std::string &date,
                       const std::string &time);

    bool get_logging() const { return doLogging && !singleProcOnly; }
    void set_logging(bool on_off) { doLogging = on_off; }

    // The get_field and put_field functions are just a wrapper around the
    // pure virtual get_field_internal and put_field_internal functions,
    // but this lets me add some debug/checking/common code to the
    // functions without having to do it in the calling code or in the
    // derived classes code.  This also fulfills the heuristic that a
    // public interface should not contain pure virtual functions.
    template <typename T>
    int64_t get_field(const T *reg, const Field &field, void *data, size_t data_size) const
    {
      IOSS_FUNC_ENTER(m_);
      verify_and_log(reg, field, 1);
      int64_t retval = get_field_internal(reg, field, data, data_size);
      verify_and_log(nullptr, field, 1);
      return retval;
    }

    template <typename T>
    int64_t put_field(const T *reg, const Field &field, void *data, size_t data_size) const
    {
      IOSS_FUNC_ENTER(m_);
      verify_and_log(reg, field, 0);
      int64_t retval = put_field_internal(reg, field, data, data_size);
      verify_and_log(nullptr, field, 0);
      return retval;
    }

    /** Determine whether application will make field data get/put calls parallel consistently.
     *
     *  True is default and required for parallel-io databases.
     *  Even if false, metadata operations must be called by all processors.
     *
     *  \returns True if application will make field data get/put calls parallel consistently.
     *
     */
    bool is_parallel_consistent() const { return isParallelConsistent; }
    void set_parallel_consistency(bool on_off) { isParallelConsistent = on_off; }

    bool get_use_generic_canonical_name() const { return useGenericCanonicalName; }
    void set_use_generic_canonical_name(bool yes_no) { useGenericCanonicalName = yes_no; }

    bool ignore_database_names() const { return ignoreDatabaseNames; }
    void ignore_database_names(bool yes_no) { ignoreDatabaseNames = yes_no; }

    /** \brief Get the length of the longest name in the database file.
     *
     *  \returns The length, or 0 for unlimited.
     */
    virtual int  maximum_symbol_length() const { return 0; } // Default is unlimited...
    virtual void set_maximum_symbol_length(int /* requested_symbol_size */) {
    } // Default does nothing...

    char get_field_separator() const { return fieldSeparator; }
    bool get_field_recognition() const { return enableFieldRecognition; }
    void set_field_separator(char separator);
    void set_field_recognition(bool yes_no) { enableFieldRecognition = yes_no; }

    void set_lower_case_variable_names(bool true_false) const
    {
      lowerCaseVariableNames = true_false;
    }

    /* \brief Set the method used to split sidesets into homogeneous blocks.
     *
     *  \param[in] split_type The desired method.
     *
     */
    void set_surface_split_type(Ioss::SurfaceSplitType split_type) { splitType = split_type; }
    Ioss::SurfaceSplitType get_surface_split_type() const { return splitType; }

    void set_block_omissions(const std::vector<std::string> &omissions,
                             const std::vector<std::string> &inclusions = {});

    void get_block_adjacencies(const Ioss::ElementBlock *eb,
                               std::vector<std::string> &block_adjacency) const
    {
      return get_block_adjacencies__(eb, block_adjacency);
    }
    void compute_block_membership(Ioss::SideBlock *         efblock,
                                  std::vector<std::string> &block_membership) const
    {
      return compute_block_membership__(efblock, block_membership);
    }

    AxisAlignedBoundingBox get_bounding_box(const Ioss::ElementBlock *eb) const;
    AxisAlignedBoundingBox get_bounding_box(const Ioss::StructuredBlock *sb) const;

    virtual int  int_byte_size_db() const = 0; //! Returns 4 or 8
    int          int_byte_size_api() const;    //! Returns 4 or 8
    virtual void set_int_byte_size_api(Ioss::DataSize size) const;

    /*!
     * The owning region of this database.
     */
    Region *get_region() const { return region_; }

    /*!
     *     The overlay_count specifies the number of restart outputs
     *     which will be overlaid on top of the currently written
     *     step before advancing to the next step on the restart
     *     database.
     *
     *     For example, if restarts are being output every 0.1
     *     seconds and the overlay count is specified as 2, then
     *     restart will write time 0.1 to step 1 of the database.
     *     It will then write 0.2 and 0.3 also to step 1.  It will
     *     then increment the database step and write 0.4 to step 2;
     *     overlay 0.5 and 0.6 on step 2...  At the end of the
     *     analysis, assuming it runs to completion, the database
     *     would have times 0.3, 0.6, 0.9, ... However, if there
     *     were a problem during the analysis, the last step on the
     *     database would contain an intermediate step.
     *
     *     The cycle_count specifies the number of restart steps
     *     which will be written to the restart database before
     *     previously written steps are overwritten.
     *
     *     For example, if the cycle count is 5 and restart is written every 0.1
     *     seconds, the restart system will write data at times 0.1, 0.2, 0.3,
     *     0.4, 0.5 to the database.
     *
     *     It will then overwrite the first step with data from time 0.6, the
     *     second with time 0.7.  At time 0.8, the database would contain data at
     *     times 0.6, 0.7, 0.8, 0.4, 0.5.  Note that time will not necessarily be
     *     monotonically increasing on a database that specifies the cycle
     *     count.
     *
     *     The cycle count and overlay count can both be used at the same time
     *     also.  The basic formula is:
     *
     *            db_step = (((output_step - 1) / overlay) % cycle) + 1
     *
     *     where "output_step" is the step that this would have been on the
     *     database in a normal write (1,2,3,....) and "db_step" is the step
     *     number that this will be written to.
     *
     *     If you only want the last step available on the database,
     *     use set_cycle_count(1)
     */
    void set_cycle_count(int count) const { cycleCount = count; }
    int  get_cycle_count() const { return cycleCount; }
    void set_overlay_count(int count) const { overlayCount = count; }
    int  get_overlay_count() const { return overlayCount; }
    void set_file_per_state(bool yes_no) const { filePerState = yes_no; }
    bool get_file_per_state() const { return filePerState; }

    void set_time_scale_factor(double factor) { timeScaleFactor = factor; }

    const Ioss::ParallelUtils &  util() const { return util_; }
    const Ioss::PropertyManager &get_property_manager() const { return properties; }
    /** \brief Get the processor that this mesh database is on.
     *
     *  \returns The processor that this mesh database is on.
     */
    int  parallel_rank() const { return myProcessor; }
    int  parallel_size() const { return util().parallel_size(); }
    bool is_parallel() const { return isParallel; }

    void progress(const std::string &output) const
    {
      if (m_enableTracing) {
        util().progress(output);
      }
    }

  protected:
    DatabaseIO(Region *region, std::string filename, Ioss::DatabaseUsage db_usage,
               MPI_Comm communicator, const Ioss::PropertyManager &props);

    /*!
     * The properties member data contains properties that can be
     * used to set database-specific options.  Examples include
     * compression, name lengths, integer sizes, floating point
     * sizes. By convention, the property name is all
     * uppercase. Some existing properties recognized by the Exodus
     * DatabaseIO class are:
     *
     * | Property              | Value
     * |-----------------------|-------------------
     * | COMPRESSION_LEVEL     | In the range [0..9]. A value of 0 indicates no compression
     * | COMPRESSION_SHUFFLE   | (true/false) to enable/disable hdf5's shuffle compression
     * algorithm.
     * | FILE_TYPE             | netcdf4
     * | MAXIMUM_NAME_LENGTH   | Maximum length of names that will be returned/passed via api call.
     * | INTEGER_SIZE_DB       | 4 or 8 indicating byte size of integers stored on the database.
     * | INTEGER_SIZE_API      | 4 or 8 indicating byte size of integers used in api functions.
     * | LOGGING               | (true/false) to enable/disable logging of field input/output
     */

    Ioss::PropertyManager properties;

    /*!
     * Utility function that may be used by derived classes.
     * Determines whether all elements in the model have the same side
     * topology.  This can be used to speed-up certain algorithms since
     * they don't have to check each side (or group of sides)
     * individually.
     */
    void             set_common_side_topology() const;
    ElementTopology *commonSideTopology{nullptr};

    template <typename T>
    void create_groups(const std::string &property_name, EntityType type,
                       const std::string &type_name, const T *set_type);
    template <typename T>
    void create_group(EntityType type, const std::string &type_name,
                      const std::vector<std::string> &group_spec, const T *set_type);

    // Create new sets as groups of existing exodus sets...
    void handle_groups();

    /*!
     * Filename that this Database is connected with.  Derived
     * DatabaseIO classes may need to change this if the passed  in
     * filename is not the same as the filename actually used  E.g.,
     * the Ioex_DatabaseIO (exodusII) changes if this is a parallel
     * run since the passed in filename is just the basename, not the
     * processor-specific filename.
     */
    std::string         originalDBFilename;
    std::string         DBFilename;
    mutable std::string decodedFilename;

    /*!
     * `bbName` is a temporary swizzled name which resides inside Burst Buffer namespace.
     * This is a private trivial mapped name vs original `DBFilename` (which resides in
     * permament storage backed by parallel filesystem.
     * `dwPath` is global BB mountpoint for current job with requested capacity via SLURM \c \#DW
     * directive. `usingDataWarp`  -- a boolean, for convenience of use so that we don't have to do
     * getenv() calls to see if BB present.
     */
    mutable std::string bbName{};
    mutable std::string pfsName{};
    mutable std::string dwPath{};
    mutable bool        usingDataWarp{false};

    mutable Ioss::State dbState{STATE_INVALID};

    bool isParallel{false}; //!< true if running in parallel
    int  myProcessor{0};    //!< number of processor this database is for

    int64_t nodeCount{0};
    int64_t elementCount{0};

    /*!
     * Check the topology of all face/element pairs in the model and
     * fill the "TopoContainer faceTopology" variable with the
     * unique pairs.  This information is used for the
     * faceblock/facesets and edgeblock/edgesets. If the
     * 'topo_dimension' is 2, then face/element pairs are generated; if
     * 'topo_dimension' is 1, then edge/element pairs are generated.
     */
    void check_side_topology() const;

    /// Used to speed up faceblock/edgeblock calculations.
    TopoContainer sideTopology;

    /*! Typically used for restart output, but can be used for all output...
     * Maximum number of states on the output file.  Overwrite the existing
     * steps in a cyclic manner once exceed this count.  Note that this breaks
     * the convention that times be monotonically increasing on an exodusII file.
     * Used by derived classes if they support this capability...
     */
    mutable int cycleCount{0};

    mutable int overlayCount{0};

    /*! EXPERIMENTAL If this is true, then each state (timestep)
     *  output will be directed to a separate file.  Currently this is
     *  only implemented for the exodus (parallel or serial, single
     *  file or fpp) database type.
     */
    mutable bool filePerState{false};

    /*! Scale the time read/written from/to the file by the specified
      scaleFactor.  If the datbase times are 0.1, 0.2, 0.3 and the
      scaleFactor is 20, then the application will think that the
      times read are 20, 40, 60.

      If specified for an output database, then the analysis time
      is divided by the scaleFactor time prior to output.
    */
    double timeScaleFactor{1.0};

    Ioss::SurfaceSplitType splitType{SPLIT_BY_TOPOLOGIES};
    Ioss::DatabaseUsage    dbUsage;
    mutable Ioss::DataSize dbIntSizeAPI{USE_INT32_API};
    mutable bool           lowerCaseVariableNames{true};
    bool                   usingParallelIO{false};

    // List of element blocks that should be omitted or included from
    // this model.  Surfaces will take this into account while
    // splitting; however, node and nodesets will not be filtered
    // (perhaps this will be done at a later time...)  NOTE: All local
    // element ids and offsets are still calculated assuming that the
    // blocks exist in the model...
    // Only one of these can have values and the other must be empty.
    std::vector<std::string> blockOmissions;
    std::vector<std::string> blockInclusions;

    std::vector<std::string> informationRecords;
    std::vector<std::string> qaRecords;

    //---Node Map -- Maps internal (1..NUMNP) ids to global ids used on the
    //               application side.   global = nodeMap[local]
    mutable Ioss::Map nodeMap{"node", DBFilename, myProcessor};
    mutable Ioss::Map edgeMap{"edge", DBFilename, myProcessor};
    mutable Ioss::Map faceMap{"face", DBFilename, myProcessor};
    mutable Ioss::Map elemMap{"element", DBFilename, myProcessor};

    mutable std::vector<std::vector<bool>> blockAdjacency;

    virtual void openDatabase__() const;
    virtual void closeDatabase__() const;
    virtual void flush_database__() const {}

  private:
    virtual bool ok__(bool /* write_message */, std::string * /* error_message */,
                      int *bad_count) const
    {
      if (bad_count != nullptr) {
        *bad_count = 0;
      }
      return dbState != Ioss::STATE_INVALID;
    }

    virtual int64_t node_global_to_local__(int64_t global, bool must_exist) const
    {
      return nodeMap.global_to_local(global, must_exist);
    }

    virtual int64_t element_global_to_local__(int64_t global) const
    {
      return elemMap.global_to_local(global);
    }

    virtual void release_memory__()
    {
      nodeMap.release_memory();
      edgeMap.release_memory();
      faceMap.release_memory();
      elemMap.release_memory();
    }

    virtual bool open_group__(const std::string & /* group_name */) { return false; }
    virtual bool create_subgroup__(const std::string & /* group_name */) { return false; }

    virtual bool begin__(Ioss::State state) = 0;
    virtual bool end__(Ioss::State state)   = 0;

    virtual void read_meta_data__() = 0;
    virtual void get_step_times__() {}

    virtual bool begin_state__(int state, double time);
    virtual bool end_state__(int state, double time);

    void get_block_adjacencies__(const Ioss::ElementBlock *eb,
                                 std::vector<std::string> &block_adjacency) const;

    virtual void compute_block_membership__(Ioss::SideBlock * /* efblock */,
                                            std::vector<std::string> & /* block_membership */) const
    {
    }

    void compute_block_adjacencies() const;

    void verify_and_log(const GroupingEntity *ge, const Field &field, int in_out) const;

    virtual int64_t get_field_internal(const Region *reg, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const NodeBlock *nb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const EdgeBlock *nb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const FaceBlock *nb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const ElementBlock *eb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const SideBlock *fb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const NodeSet *ns, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const EdgeSet *ns, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const FaceSet *ns, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const ElementSet *ns, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const SideSet *fs, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const CommSet *cs, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t get_field_internal(const StructuredBlock * /*sb*/, const Field & /*field*/,
                                       void * /*data*/, size_t /*data_size*/) const
    {
      return 0;
    }

    virtual int64_t put_field_internal(const Region *reg, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const NodeBlock *nb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const EdgeBlock *nb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const FaceBlock *nb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const ElementBlock *eb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const SideBlock *fb, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const NodeSet *ns, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const EdgeSet *ns, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const FaceSet *ns, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const ElementSet *ns, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const SideSet *fs, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const CommSet *cs, const Field &field, void *data,
                                       size_t data_size) const = 0;
    virtual int64_t put_field_internal(const StructuredBlock * /*sb*/, const Field & /*field*/,
                                       void * /*data*/, size_t /*data_size*/) const
    {
      return 0;
    }

    DatabaseIO()                   = delete;
    DatabaseIO(const DatabaseIO &) = delete;
    DatabaseIO &operator=(const DatabaseIO &) = delete;

    mutable std::map<std::string, AxisAlignedBoundingBox> elementBlockBoundingBoxes;

    Ioss::ParallelUtils util_; // Encapsulate parallel and other utility functions.
#if defined(IOSS_THREADSAFE)
  protected:
    mutable std::mutex m_;

  private:
#endif
    Region *region_{nullptr};
    char    fieldSeparator{'_'};
    bool    enableFieldRecognition{true};
    bool    isInput;
    bool    isParallelConsistent{
        true}; // True if application will make field data get/put calls parallel
                  // consistently.
                  // True is default and required for parallel-io databases.
    // Even if false, metadata operations must be called by all processors

    bool singleProcOnly;   // True if history or heartbeat which is only written from proc 0...
    bool doLogging{false}; // True if logging field input/output
    bool useGenericCanonicalName{
        false}; // True if "block_id" is used as canonical name instead of the name
    // given on the mesh file e.g. "fireset".  Both names are still aliases.
    bool ignoreDatabaseNames{false}; // True if "block_{id}" used as canonical name; ignore any
                                     // names on database.
    mutable bool blockAdjacenciesCalculated{false}; // True if the lazy creation of
    // block adjacencies has been calculated.

    bool m_timeStateInOut{false};
    bool m_enableTracing{false};
    std::chrono::time_point<std::chrono::high_resolution_clock>
        m_stateStart; // Used for optional output step timing.
  };
} // namespace Ioss
#endif
