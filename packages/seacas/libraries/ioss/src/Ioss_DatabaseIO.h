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

#ifndef IOSS_Ioss_DatabaseIO_h
#define IOSS_Ioss_DatabaseIO_h

#include <Ioss_BoundingBox.h> // for AxisAlignedBoundingBox
#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>         // for DatabaseUsage, etc
#include <Ioss_DataSize.h>        // for DataSize
#include <Ioss_EntityType.h>      // for EntityType
#include <Ioss_ParallelUtils.h>   // for ParallelUtils
#include <Ioss_PropertyManager.h> // for PropertyManager
#include <Ioss_State.h>           // for State, State::STATE_INVALID
#include <Ioss_SurfaceSplit.h>    // for SurfaceSplitType
#include <map>                    // for map
#include <stddef.h>               // for size_t, nullptr
#include <stdint.h>               // for int64_t
#include <string>                 // for string
#include <utility>                // for pair
#include <vector>                 // for vector
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
}

namespace Ioss {
  class EntityBlock;

  // Contains (parent_element, side) topology pairs
  typedef std::vector<std::pair<const ElementTopology *, const ElementTopology *>> TopoContainer;

  /** \brief An input or output Database.
   *
   */
  class DatabaseIO
  {
  public:
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
    virtual bool ok(bool write_message = false, std::string *error_message = nullptr,
                    int *bad_count = nullptr) const
    {
      if (bad_count != nullptr) {
        *bad_count = 0;
      }
      return dbState != Ioss::STATE_INVALID;
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
    virtual int64_t node_global_to_local(int64_t global, bool must_exist) const = 0;
    virtual int64_t element_global_to_local(int64_t global) const = 0;

    virtual ~DatabaseIO();

    // Eliminate as much memory as possible, but still retain meta data information
    // Typically, eliminate the maps...
    virtual void release_memory() {}

    /** \brief Get the file name associated with the database.
     *
     *  \returns The database file name.
     */
    std::string get_filename() const { return DBFilename; }

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

    /** \brief Determine whether the database needs information about proces ownership of nodes.
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

    virtual void openDatabase() const {}
    virtual void closeDatabase() const {}
    virtual void flush_database() const {}

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
    virtual bool open_group(const std::string &group_name) { return false; }

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
    virtual bool create_subgroup(const std::string &group_name) { return false; }

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
    virtual bool begin(Ioss::State state) = 0;

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
    virtual bool end(Ioss::State state) = 0;

    virtual bool begin_state(Region *region, int state, double time);
    virtual bool end_state(Region *region, int state, double time);

    // Metadata-related functions.
    virtual void read_meta_data() = 0;
    virtual void get_step_times() {}

    virtual bool internal_edges_available() const { return false; }
    virtual bool internal_faces_available() const { return false; }

    // Information Records:

    /** \brief Get all information records (informative strings) for the database.
     *
     *  \returns The informative strings.
     */
    const std::vector<std::string> &get_information_records() const { return informationRecords; }
    void add_information_records(const std::vector<std::string> &info);
    void add_information_record(const std::string &info);

    // QA Records:

    /** \brief Get all QA records, each of which consists of 4 strings, from the database
     *
     *  The 4 strings that make up a databse QA record are:
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
    // derived classes code.  This also fulfills the hueristic that a
    // public interface should not contain pure virtual functions.
    template <typename T>
    int64_t get_field(const T *reg, const Field &field, void *data, size_t data_size) const
    {
      verify_and_log(reg, field, 1);
      int64_t retval = get_field_internal(reg, field, data, data_size);
      verify_and_log(nullptr, field, 1);
      return retval;
    }

    template <typename T>
    int64_t put_field(const T *reg, const Field &field, void *data, size_t data_size) const
    {
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
    virtual int maximum_symbol_length() const { return 0; } // Default is unlimited...
    virtual void
    set_maximum_symbol_length(int /* requested_symbol_size */){}; // Default does nothing...

    char get_field_separator() const;
    void set_field_separator(const char separator);
    void set_lower_case_variable_names(bool true_false) const
    {
      lowerCaseVariableNames = true_false;
    }

    /* \brief Set the method used to split sidesets into homogenous blocks.
     *
     *  \param[in] split_type The desired method.
     *
     */
    void set_surface_split_type(Ioss::SurfaceSplitType split_type) { splitType = split_type; }
    Ioss::SurfaceSplitType get_surface_split_type() const { return splitType; }

    void set_block_omissions(const std::vector<std::string> &omissions);

    virtual void get_block_adjacencies(const Ioss::ElementBlock *eb,
                                       std::vector<std::string> &block_adjacency) const
    {
    }
    virtual void compute_block_membership(Ioss::SideBlock *         efblock,
                                          std::vector<std::string> &block_membership) const
    {
    }

    AxisAlignedBoundingBox get_bounding_box(const Ioss::ElementBlock *eb) const;
    AxisAlignedBoundingBox get_bounding_box(const Ioss::StructuredBlock *sb) const;

    int          int_byte_size_api() const; //! Returns 4 or 8
    virtual void set_int_byte_size_api(Ioss::DataSize size) const;

    /*!
     * The owning region of this database.
     */
    Region *get_region() const { return region_; }

    /*!
     *     The overlay_count specifies the number of restart outputs
     *     which will be overlayed on top of the currently written
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
    void set_overlay_count(int count) const { overlayCount = count; }

    void set_time_scale_factor(double factor) { timeScaleFactor = factor; }

    const Ioss::ParallelUtils &util() const { return util_; }

    /** \brief Get the processor that this mesh database is on.
     *
     *  \returns The processor that this mesh database is on.
     */
    int  parallel_rank() const { return myProcessor; }
    int  parallel_size() const { return util().parallel_size(); }
    bool is_parallel() const { return isParallel; }

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
    ElementTopology *commonSideTopology;

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
    std::string DBFilename;

    mutable Ioss::State dbState;

    bool isParallel;  //!< true if running in parallel
    int  myProcessor; //!< number of processor this database is for

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
    mutable int cycleCount;

    mutable int overlayCount;

    /*! Scale the time read/written from/to the file by the specified
      scaleFactor.  If the datbase times are 0.1, 0.2, 0.3 and the
      scaleFactor is 20, then the application will think that the
      times read are 20, 40, 60.

      If specified for an output database, then the analysis time
      is divided by the scaleFactor time prior to output.
    */
    double timeScaleFactor;

    Ioss::SurfaceSplitType splitType;
    Ioss::DatabaseUsage    dbUsage;
    mutable Ioss::DataSize dbIntSizeAPI;
    mutable bool           lowerCaseVariableNames;
    bool usingParallelIO;

    // List of element blocks that should be omitted from this model.
    // Surfaces will take this into account while splitting;
    // however, node and nodesets will not be filtered
    // (perhaps this will be done at a later time...)
    // NOTE: All local element ids and offsets are still calculated
    //       assuming that the blocks exist in the model...
    std::vector<std::string> blockOmissions;

    std::vector<std::string> informationRecords;
    std::vector<std::string> qaRecords;

  private:
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
    virtual int64_t get_field_internal(const StructuredBlock *sb, const Field &field, void *data,
                                       size_t data_size) const
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
    virtual int64_t put_field_internal(const StructuredBlock *sb, const Field &field, void *data,
                                       size_t data_size) const
    {
      return 0;
    }

    DatabaseIO()                   = delete;
    DatabaseIO(const DatabaseIO &) = delete;
    DatabaseIO &operator=(const DatabaseIO &) = delete;

    mutable std::map<std::string, AxisAlignedBoundingBox> elementBlockBoundingBoxes;

    Ioss::ParallelUtils util_; // Encapsulate parallel and other utility functions.
    Region *            region_;
    bool                isInput;
    bool isParallelConsistent; // True if application will make field data get/put calls parallel
                               // consistently.
                               // True is default and required for parallel-io databases.
    // Even if false, metadata operations must be called by all processors

    bool singleProcOnly; // True if history or heartbeat which is only written from proc 0...
    bool doLogging;      // True if logging field input/output
    bool
        useGenericCanonicalName; // True if "block_id" is used as canonical name instead of the name
    // given on the mesh file e.g. "fireset".  Both names are still aliases.
    bool ignoreDatabaseNames; // True if "block_{id}" used as canonical name; ignore any names on
                              // database.
  };
}
#endif
