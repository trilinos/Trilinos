/*------------------------------------------------------------------------*/
/*  Copyright 2010, 2011, 2012, 2013 Sandia Corporation.                  */
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
#include <stk_util/environment/ReportHandler.hpp>
#include <Ioss_PropertyManager.h>
#include <init/Ionit_Initializer.h>
#include <Teuchos_RCP.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace Ioss {
  class DatabaseIO;
  class Region;
}

namespace stk {
  namespace mesh {
    class Part;
    class BulkData;
    class Selector;
  }
  namespace io {
    typedef stk::mesh::Field<double, stk::mesh::Cartesian> CoordinateFieldType;
    static const std::string CoordinateFieldName("coordinates");

    class MeshData {
        // Used to maintain state between the meta data and bulk data
        // portions of the mesh generation process for use cases.
      public:
        /**
         * \param[in] comm  MPI Communicator to be used for all parallel
         * communication needed to generate the mesh.
         */
        MeshData(MPI_Comm comm);
        MeshData();

        ~MeshData();

        /**
         * Set the output Ioss::Region directly instead of letting it be
         * created by MeshData during the create_output_mesh() call.
         */
        void set_output_io_region(Teuchos::RCP<Ioss::Region> ioss_output_region);

        /**
         * Set the input Ioss::Region directly instead of letting it be
         * created by MeshData during the create_input_mesh(type,
         * filename) call. After setting the input io region, you would
         * then either set the metadata manually using the
         * set_meta_data() call, or call the no-argument
         * create_input_mesh() function which will then create a meta
         * data corresponding to the data in the Ioss::Region.
         */
        void set_input_io_region(Teuchos::RCP<Ioss::Region> ioss_input_region);

        Teuchos::RCP<Ioss::DatabaseIO> input_io_database()  { return m_input_database;   }
        Teuchos::RCP<Ioss::Region> input_io_region()
        {
          if (Teuchos::is_null(m_input_region) && !Teuchos::is_null(m_input_database)) {
            create_input_mesh();
          }
          return m_input_region;
        }

        Teuchos::RCP<Ioss::Region> output_io_region() { return m_output_region;  }

        Teuchos::RCP<stk::mesh::Selector> selector()  { return m_anded_selector; }

        /**
         * The selector is an optional selector that can be provided
         * which will be used to filter the mesh entities (node,
         * element, face, edge) associated with an
         * Ioss Entity (element block, node block, node set, ...). The
         * optional selector will be 'anded' with the normal selector
         * (typically locally owned part) used to associate entities
         * when generating the output database.
         */
        void set_selector(Teuchos::RCP<stk::mesh::Selector> selector)
        {	m_anded_selector = selector; }
        void set_selector(stk::mesh::Selector &selector)
        { set_selector(Teuchos::rcpFromRef(selector));}

        /**
         * Set meta data directly with your own meta data. If this is
         * not called, a meta data will be created during the
         * create_input_mesh() call.
         */
        void set_meta_data(Teuchos::RCP<stk::mesh::MetaData> arg_meta_data);
        void set_meta_data(stk::mesh::MetaData &arg_meta_data)
        { set_meta_data(Teuchos::rcpFromRef(arg_meta_data));}

        /**
         * Set bulk data directly with your own bulk data. If this is
         * not called prior to the populate_bulk_data() call, it will be
         * created automatically using the communicator of the m_input_region.
         * If meta data is not already set, then set the meta data from the
         * bulk data's metadata
         */
        void set_bulk_data(Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data);
        void set_bulk_data(stk::mesh::BulkData &arg_bulk_data)
        { set_bulk_data(Teuchos::rcpFromRef(arg_bulk_data));}

        /**
         * Create the Ioss::DatabaseIO associated with the specified filename
         * and type (exodus by default). The routine checks that the
         * file exists and is readable and will return false if not.
         * No meta_data or bulk_data is created at this time, but the
         * input_io_database() will be valid.
         *
         * The Ioss::DatabaseIO can be accessed vi the input_io_database()
         * method if you need to set some options on the database prior
         * to it being read.
         *
         * \param[in] filename If the mesh type is file based ("exodus"),
         * then this contains the full pathname to the file containing the
         * mesh information.  If the mesh type is a generated type, then
         * this parameter contains data used by the generation routines.
         * See the output from the show_mesh_help() function for details.
         *
         *          * \param[in] type   The format of the mesh that will be
         * "read".  Valid types are "exodus", "generated", "pamgen".
         *
         */
        bool open_mesh_database(const std::string &filename,
                                const std::string &type);

        /**
         * Create Ioss::DatabaseIO associated with the specified filename using
         * the default filetype (typically "exodus"). If the filename is
         * prepended with a type followed by a colon (e.g., "generated:10x10x10"),
         * then use that type.  Valid input types are:
         *   exodus, dof, pamgen, and possibly others.
         *
         * @param filename If the mesh type is file based ("exodus"),
         * then this contains the full pathname to the file containing the
         * mesh information.  If the mesh type is a generated type, then
         * this parameter contains data used by the generation routines.
         * Optionally prepended by a filetype and a colon.
         */
        bool open_mesh_database(const std::string &filename);

        /**
         * Read/Generate the metadata for mesh of the specified type. By
         * default, all entities in the mesh (nodeblocks, element blocks,
         * nodesets, sidesets) will have an associated stk mesh part
         * created for it.
         *
         * If m_input_region is non-null then this is assumed to be a
         * valid Ioss::Region* that should be used instead of opening
         * the file and creating a new Ioss::Region. This can be set via
         * the set_input_io_region(Ioss::Region *region) call prior to
         * calling this function.
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
         * The meta_data will not be committed by this function, so the
         * caller will need to call meta_data.commit() after the
         * function returns.
         */
        void create_input_mesh();

        /**
         * Read/Generate the bulk data for the mesh.  This function will
         * create all stk mesh entities (nodes, elements) with the
         * correct nodeal coordinates, element connectivity, element
         * attribute data, and nodeset and sideset membership.  Note
         * that meta_data.commit() needs to be called prior to calling
         * this function.
         */
        void populate_bulk_data();

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
        void define_input_fields();

        /**
         * For all transient input fields defined either manually or via
         * the define_input_fields() function, read the data at the
         * specified database step 'step' (1-based) and populate the stk
         * data structures with those values.
         */
        void process_input_request(int step);

        /**
         * For all transient input fields defined either manually or via
         * the define_input_fields() function, read the data at the
         * specified database time 'time' and populate the stk
         * data structures with those values.  The database time closest
         * to the specified time will be used with no interpolation (yet).
         */
        void process_input_request(double time);

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
         */
        void create_output_mesh(const std::string &filename);
        void create_output_mesh();

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
        void define_output_fields(bool add_all_fields = false);

        /**
         * Add a transient step to the mesh database at time 'time' and
         * output the data for all defined fields to the database.
         */
        int process_output_request(double time);
        int process_output_request(double time, const std::set<const stk::mesh::Part*> &exclude);

        bool meta_data_is_set() const
        {
          return !Teuchos::is_null(m_meta_data);
        }
        bool bulk_data_is_set() const
        {
          return !Teuchos::is_null(m_bulk_data);
        }

        stk::mesh::MetaData &meta_data()
        {
          ThrowRequire( !Teuchos::is_null(m_meta_data)) ;
          return *m_meta_data;
        }

        stk::mesh::BulkData &bulk_data()
        {
          ThrowRequire( !Teuchos::is_null(m_bulk_data)) ;
          return *m_bulk_data;
        }

        const stk::mesh::MetaData &meta_data() const
        {
          ThrowRequire( !Teuchos::is_null(m_meta_data)) ;
          return *m_meta_data;
        }

        const stk::mesh::BulkData &bulk_data() const
        {
          ThrowRequire( !Teuchos::is_null(m_bulk_data)) ;
          return *m_bulk_data;
        }

        /*!
         * Return the coordinate field for this mesh.
         */
        CoordinateFieldType & get_coordinate_field()
        {
          CoordinateFieldType * coord_field = meta_data().get_field<CoordinateFieldType>(CoordinateFieldName);
          ThrowRequire( coord_field != NULL);
          return * coord_field;
        }

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
         */
        Ioss::PropertyManager m_property_manager;

      private:
        void internal_process_output_request(int step, const std::set<const stk::mesh::Part*> &exclude);

        MPI_Comm communicator_;

        Teuchos::RCP<Ioss::DatabaseIO> m_input_database;
        Teuchos::RCP<Ioss::Region> m_input_region;
        Teuchos::RCP<Ioss::Region> m_output_region;

        Teuchos::RCP<stk::mesh::MetaData>  m_meta_data;
        Teuchos::RCP<stk::mesh::BulkData>  m_bulk_data;

        /*!
         * An optional selector used for filtering entities on the
         * output database. This can be used for specifying
         * active/inactive entities.  If present, then this selector is
         * *anded* with the normal selectors use for output
         */
        Teuchos::RCP<stk::mesh::Selector> m_anded_selector;

        MeshData(const MeshData&); // Do not implement
        MeshData& operator=(const MeshData&); // Do not implement
    };

    /** Output a help message showing the valid options for the mesh
	read and options for generated mesh.
     */
    void show_mesh_help();
  }
}
#endif
