#ifndef STK_IO_Inputfile
#define STK_IO_Inputfile

#include <Teuchos_RCP.hpp>              // for is_null, RCP::operator->, etc
#include <stk_mesh/base/Types.hpp>
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose
#include <stk_io/MeshField.hpp>   

namespace Ioss {
  class PropertyManager;
  class GroupingEntity;
  class Region;
  class DatabaseIO;
}
  
namespace stk {
  namespace mesh {
    class MetaData;
    class BulkData;
    class Part;
  }
  
  namespace io {
    class InputFile
    {
    public:
      enum PeriodType {
	CYCLIC,   /*< Cycles 0 1 2 0 1 2 0 1 2 ... */
	REVERSING /*< Cycles 0 1 2 1 0 1 2 1 0 ... */
      };

      InputFile(std::string filename,
		MPI_Comm communicator,
		const std::string &type,
		DatabasePurpose purpose,
		Ioss::PropertyManager& property_manager);
      InputFile(Teuchos::RCP<Ioss::Region> ioss_input_region);

      ~InputFile()
      {}

      void create_ioss_region();
      //      void add_input_field(stk::mesh::FieldBase &field, const std::string &db_name);
      void add_input_field(const stk::io::MeshField &mesh_field);
      void add_all_mesh_fields_as_input_fields(stk::mesh::MetaData &meta, MeshField::TimeMatchOption tmo);
      double read_defined_input_fields(double time, std::vector<stk::io::MeshField> *missing,
				       stk::mesh::BulkData &bulk);
      double read_defined_input_fields(int step, std::vector<stk::io::MeshField> *missing,
				       stk::mesh::BulkData &bulk);
      void get_global_variable_names(std::vector<std::string> &names);

      void build_field_part_associations(stk::mesh::BulkData &bulk);

      Teuchos::RCP<Ioss::Region> get_input_io_region()
      {
	if (Teuchos::is_null(m_region) && !Teuchos::is_null(m_database)) {
	  create_ioss_region();
	}
	return m_region;
      }

      InputFile& set_offset_time(double offset_time);
      InputFile& set_periodic_time(double period_length, double startup_time = 0.0,
				   PeriodType ptype = CYCLIC);
      InputFile& set_start_time(double start_time);
      InputFile& set_stop_time(double stop_time);
  
    private:
      void build_field_part_associations(stk::io::MeshField &mesh_field,
					 const stk::mesh::Part &part,
					 const stk::mesh::EntityRank rank,
					 Ioss::GroupingEntity *io_entity);

      void report_missing_fields(std::vector<stk::io::MeshField> *missing) const;

      DatabasePurpose m_db_purpose;
      Teuchos::RCP<Ioss::DatabaseIO> m_database;
      Teuchos::RCP<Ioss::Region> m_region;
      std::vector<stk::io::MeshField> m_fields;

      /*@{*/

      /**
       * For input interpolation only: The 'startupTime' and
       * 'periodLength' are used to support input of periodic data.  The
       * 'startupTime' specifies the length of time prior to the start of
       * the periodic behavior; if the application time is less than
       * 'startupTime', then it is passed unchanged to the database.  Once
       * the 'startupTime' is exceeded, then the variables on the database
       * enter their periodic behavior the length of the period is
       * specified by 'periodLength'.
       *
       * If:  t_app < t_start_time || t_app > t_stop_time then
       *  don't do anything...
       *
       * If:  t_app < t_startup, then
       *  \code
       *      t_db = t_app + t_offset
       *  \endcode
       * Else:
       *   If: ptype == CYCLIC
       *    \code
       *      t_db = t_startup + mod(t_app-t_startup, period_length) + t_offset
       *    \endcode
       *   If: ptype == REVERSING
       *    \code
       *      tpm = mod(t_app-t_startup, 2*period_length)
       *      if (tpm <= period_length)
       *         t_db = t_startup + tpm + t_offset
       *      else 
       *         t_db = t_startup + (2*period_length - tpm) + t_offset
       *    \endcode
       *
       * Currently only used for Input interpolation.
       * Ignored for all other cases.
       */
      double m_startupTime;

      /** See InputFile::startupTime */
      double m_periodLength;
      /** See InputFile::startupTime */
      double m_offsetTime;
      /** See InputFile::startupTime */
      double m_startTime;
      /** See InputFile::startupTime */
      double m_stopTime;
      /** See InputFile::startupTime */
      PeriodType m_periodType;
      
      /*@}*/

    public:
      bool m_fieldsInitialized;
      
    private:
      InputFile(const InputFile &);
      const InputFile & operator=(const InputFile &);
    };
  }
}
#endif
