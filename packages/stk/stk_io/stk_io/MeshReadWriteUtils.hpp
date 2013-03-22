/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_IO_MESHREADWRITEUTILS_HPP
#define STK_IO_MESHREADWRITEUTILS_HPP
#include <string>
#include <vector>
#include <set>
#include <stk_util/parallel/Parallel.hpp>
#include <Ioss_PropertyManager.h>

namespace Ioss {
  class Region;
}

namespace stk {
  namespace mesh {
    class Part;
    class BulkData;
    class Selector;
    namespace fem {
      class FEMMetaData;
    }
  }
  namespace io {
    class MeshData {
      // Used to maintain state between the meta data and bulk data
      // portions of the mesh generation process for use cases.
    public:
      MeshData() : m_input_region(NULL), m_output_region(NULL),
		   m_anded_selector(NULL)
      {}

      ~MeshData();

      /*! 
       * The `m_property_manager` member data contains properties that
       * can be used to set database-specific options in the
       * Ioss::DatabaseIO class.  Examples include compression, name
       * lengths, integer sizes, floating point sizes. By convention,
       * the property name is all uppercase. Some existing properties
       * recognized by the Exodus Ioex::DatabaseIO class are:
       *
       * | Property              | Value
       * |-----------------------|-------------------
       * | COMPRESSION_LEVEL     | In the range [0..9]. A value of 0 indicates no compression
       * | COMPRESSION_SHUFFLE   | (true/false) to enable/disable hdf5's shuffle compression algorithm.
       * | FILE_TYPE             | netcdf4 
       * | MAXIMUM_NAME_LENGTH   | Maximum length of names that will be returned/passed via api call.
       * | INTEGER_SIZE_DB       | 4 or 8 indicating byte size of integers stored on the database.
       * | INTEGER_SIZE_API      | 4 or 8 indicating byte size of integers used in api functions.
       * | LOGGING               | (true/false) to enable/disable logging of field input/output
       * | FIELD_SUFFIX_SEPARATOR| character separating field basename from suffix in auto field recognition.
       * |                       | E.g., if '_' (default), then a_x, a_y, a_z would be a 3d vector field 'a'.
       */
      Ioss::PropertyManager m_property_manager;

      Ioss::Region *m_input_region;
      Ioss::Region *m_output_region;

      /*!
       * An optional selector used for filtering entities on the
       * output database. This can be used for specifying
       * active/inactive entities.  If present, then this selector is
       * *anded* with the normal selectors use for output
       */
      stk::mesh::Selector *m_anded_selector;

    private:
      MeshData(const MeshData&); // Do not implement
      MeshData& operator=(const MeshData&); // Do not implement

    };

    /** Output a help message showing the valid options for the mesh
	read and options for generated mesh.
    */
    void show_mesh_help();

    /**
     * Read/Generate the metadata for mesh of the specified type. By
     * default, all entities in the mesh (nodeblocks, element blocks,
     * nodesets, sidesets) will have an associated stk mesh part
     * created for it.
     *
     * If the mesh_data argument contains a non-null m_input_region
     * data member, then this is assumed to be a valid Ioss::Region*
     * that should be used instead of opening the file and creating a
     * new Ioss::Region.
     *
     * Following this call, the 'populate_bulk_data()' function should
     * be called to read the bulk data from the mesh and generate the
     * corresponding stk mesh entities (nodes, elements, faces, ...)
     *
     * Only the non-transient data stored in the mesh database will be
     * accessed in this function.  To access any transient field data
     * that may be on the mesh database, use the
     * 'define_input_fields()' function.
     *
     * \param[in] type   The format of the mesh that will be
     * "read".  Valid types are "exodus", "generated", "pamgen".
     *
     * \param[in] filename If the mesh type is file based ("exodus"),
     * then this contains the full pathname to the file containing the
     * mesh information.  If the mesh type is a generated type, then
     * this parameter contains data used by the generation routines.
     * See the output from the show_mesh_help() function for details.
     *
     * \param[in] comm  MPI Communicator to be used for all parallel
     * communication needed to generate the mesh.
     *
     * \param[in,out] meta_data The STK meta data object which will
     * be populated with parts and fields based on the mesh model
     * described by the mesh in filename.  The meta_data will not be
     * committed by this function, so the caller will need to call
     * meta_data.commit() after the function returns.
     *
     * \param[in,out] mesh_data A small class used for maintaining
     * some state used by the stk_io routines.
     *
     */
    void create_input_mesh(const std::string &type,
			   const std::string &filename,
			   MPI_Comm comm,
			   stk::mesh::fem::FEMMetaData &metadata,
			   MeshData &mesh_data);

    /**
     * Read/Generate the bulk data for the mesh.  The bulk_data must
     * have been constructed using the meta_data passed to the
     * create_input_mesh() function and the mesh_data must also be the
     * same. This function will create all stk mesh entities (nodes,
     * elements) with the correct nodeal coordinates, element
     * connectivity, element attribute data, and nodeset and sideset
     * membership.  Note that meta_data.commit() needs to be called
     * prior to calling this function.
     */
    void populate_bulk_data(stk::mesh::BulkData &bulk_data, stk::io::MeshData &mesh_data);

    /**
     * Read/Generate the bulk data for the mesh.  The bulk_data must
     * have been constructed using the meta_data passed to the
     * create_input_mesh() function and the mesh_data must also be the
     * same. This function will create all stk mesh entities (nodes,
     * elements) with the correct nodal coordinates, element
     * connectivity, element attribute data, and nodeset and sideset
     * membership.  Note that meta_data.commit() followed by 
     * bulk_data.modification_begin() needs to be called
     * prior to calling this function. Further, bulk_data.modification_end()
     * must be called upon return from this function. The above populate_bulk_data call
     * is a wrapper for this function.
     */
    void process_mesh_bulk_data(Ioss::Region *region, stk::mesh::BulkData &bulk_data);

    /**
     * Iterate over all Ioss entities in the input mesh database and
     * define a stk field for each transient field found.  The stk
     * field will have the same name as the field on the database.
     *
     * Note that all transient fields found on the mesh database will
     * have a corresponding stk field defined.  If you want just a
     * selected subset of the database fields defined in the stk mesh,
     * you need to define the fields manually.
     *
     * To populate the stk field with data from the database, call
     * process_input_request().
     *
     */
    void define_input_fields(MeshData &mesh_data, stk::mesh::fem::FEMMetaData &meta_data);

    /**
     * For all transient input fields defined either manually or via
     * the define_input_fields() function, read the data at the
     * specified database step 'step' (1-based) and populate the stk
     * data structures with those values.
     */
    void process_input_request(MeshData &mesh_data, stk::mesh::BulkData &bulk, int step);

    /**
     * For all transient input fields defined either manually or via
     * the define_input_fields() function, read the data at the
     * specified database time 'time' and populate the stk
     * data structures with those values.  The database time closest
     * to the specified time will be used with no interpolation (yet).
     */
    void process_input_request(MeshData &mesh_data, stk::mesh::BulkData &bulk, double time);

    /**
     * For all transient input fields defined either manually or via
     * the define_input_fields() function, read the data at the
     * specified database step 'step' (1-based) and populate the stk
     * data structures with those values. Note that 
     * bulk_data.modification_begin() needs to be called prior to 
     * calling this function. Further, bulk_data.modification_end()
     * must be called upon return from this function. Also note 
     * that the two above functions are wrappers for this one.
     */
    void input_mesh_fields(Ioss::Region *region, stk::mesh::BulkData &bulk_data, int step);

    /**
     * For all transient input fields defined either manually or via
     * the define_input_fields() function, read the data at the
     * specified database time 'time' and populate the stk
     * data structures with those values.  The database time closest
     * to the specified time will be used with no interpolation (yet).  Note that 
     * bulk_data.modification_begin() needs to be called prior to 
     * calling this function. Further, bulk_data.modification_end()
     * must be called upon return from this function.
     */
    void input_mesh_fields(Ioss::Region *region, stk::mesh::BulkData &bulk_data, double time);

    /**
     * Create an exodus mesh database with the specified
     * filename. This function creates the exodus metadata which
     * is the number and type of element blocks, nodesets, and
     * sidesets; and then outputs the mesh bulk data such as the
     * node coordinates, id maps, element connectivity.  When the
     * function returns, the non-transient portion of the mesh will
     * have been defined.
     *
     * A stk part will have a corresponding exodus entity (element
     * block, nodeset, sideset) defined if the "is_io_part()" function
     * returns true.  By default, all parts read from the mesh
     * database in the create_input_mesh() function will return true
     * as will all stk parts on which the function
     * stk::io::put_io_part_attribute() was called.  The function
     * stk::io::remove_io_part_attribute(part) can be called to omit a
     * part from being output.
     *
     * \param[in] filename The full pathname to the file which will be
     * created and the mesh data written to. If the file already
     * exists, it will be overwritten.
     *
     * \param[in] comm  MPI Communicator to be used for all parallel
     * communication needed by the mesh routines.
     *
     * \param[in] bulk The STK bulk data object defining the stk mesh.
     *
     * \param[in,out] mesh_data A small class used for maintaining
     * some state used by the stk_io routines.
     *
     */
    void create_output_mesh(const std::string &filename,
			    MPI_Comm comm,
			    stk::mesh::BulkData &bulk_data,
			    MeshData &mesh_data);

    /**
     * Iterate over all stk fields and for each transient field
     * defined on a part that is output to the mesh file, define a
     * corresponding database field. The database field will have the
     * same name as the stk field.  A transient field will be defined
     * if the stk::io::is_valid_part_field() returns true.  This can
     * be set via a call to stk::io::set_field_role().
     *
     * If the 'add_all_fields' param is true, then all transient
     * stk fields will have a corresponding database field defined.
     */
    void define_output_fields(const MeshData &mesh_data,
			      const stk::mesh::fem::FEMMetaData &fem_meta,
			      bool add_all_fields = false);

    /**
     * Add a transient step to the mesh database at time 'time' and
     * output the data for all defined fields to the database.
     */
    int process_output_request(MeshData &mesh_data,
			       stk::mesh::BulkData &bulk,
			       double time,
                               const std::set<const stk::mesh::Part*> &exclude=std::set<const stk::mesh::Part*>());
    /**
     * Method to query a MeshData for the number of element blocks and the
     * number of elements in each. MeshData is input, std:vector is output
     */
    template <typename INT>
    void get_element_block_sizes(MeshData &mesh_data,
                                 std::vector<INT>& el_blocks);
  }
}
#endif
