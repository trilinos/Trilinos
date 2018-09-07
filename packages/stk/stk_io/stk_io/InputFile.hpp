// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 

#ifndef STK_IO_Inputfile
#define STK_IO_Inputfile

#include <Teuchos_RCP.hpp>              // for is_null, RCP::operator->, etc
#include <stk_mesh/base/Types.hpp>
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose
#include <stk_io/MeshField.hpp>
#include <stk_io/IossBridge.hpp>
#include "Ioss_EntityType.h"

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
      {delete m_multiStateSuffixes;}

      void create_ioss_region();
      FieldNameToPartVector get_var_names(Ioss::EntityType type, stk::mesh::MetaData& meta);
      void add_input_field(const stk::io::MeshField &mesh_field);
      void add_all_mesh_fields_as_input_fields(stk::mesh::MetaData &meta, MeshField::TimeMatchOption tmo);
      bool read_input_field(stk::io::MeshField &mf, stk::mesh::BulkData &bulk);
      double read_defined_input_fields(double time, std::vector<stk::io::MeshField> *missingFields,
				       stk::mesh::BulkData &bulk);
      double read_defined_input_fields(int step, std::vector<stk::io::MeshField> *missingFields,
				       stk::mesh::BulkData &bulk);
      double read_defined_input_fields_at_step(int step, std::vector<stk::io::MeshField> *missingFields,
                                       stk::mesh::BulkData &bulk);
      void get_global_variable_names(std::vector<std::string> &names);

      void build_field_part_associations(stk::mesh::BulkData &bulk, std::vector<stk::io::MeshField> *missing);

      void build_field_part_associations_from_grouping_entity(stk::mesh::BulkData &bulk, std::vector<stk::io::MeshField> *missingFields);

      Teuchos::RCP<Ioss::Region> get_input_io_region()
      {
	if (Teuchos::is_null(m_region) && !Teuchos::is_null(m_database)) {
	  create_ioss_region();
	}
	return m_region;
      }

      InputFile& set_offset_time(double offset_time);
      InputFile& set_scale_time(double scale_time);
      InputFile& set_periodic_time(double period_length, double startup_time = 0.0,
				   PeriodType ptype = CYCLIC);
      InputFile& set_start_time(double start_time);
      InputFile& set_stop_time(double stop_time);
  
      // Only public so easier to test...
      double map_analysis_to_db_time(double time) const;

      void set_surface_split_type(Ioss::SurfaceSplitType split_type) {
          if(!Teuchos::is_null(m_database)) {
              m_database->set_surface_split_type(split_type);
          }
      }
      Ioss::SurfaceSplitType get_surface_split_type() const {
          if(!Teuchos::is_null(m_database)) {
              return m_database->get_surface_split_type();
          }

          return Ioss::SPLIT_INVALID;
      }

      Teuchos::RCP<Ioss::DatabaseIO> get_input_database()
      {
	return m_database;
      }

      bool set_multistate_suffixes(std::vector<std::string>& multiStateSuffixes)
      {
          if(nullptr != m_multiStateSuffixes) {
              delete m_multiStateSuffixes;
              m_multiStateSuffixes = nullptr;
          }

          m_multiStateSuffixes = new std::vector<std::string>(multiStateSuffixes);
          return true;
      }

    private:
      bool process_fields_for_grouping_entity(stk::io::MeshField &mesh_field,
                                              const stk::mesh::Part &part,
                                              Ioss::GroupingEntity *io_entity,
                                              std::map<stk::mesh::FieldBase *, const stk::io::MeshField *> *missing_fields_collector_ptr = nullptr);

      bool build_field_part_associations(stk::io::MeshField &mesh_field,
					 const stk::mesh::Part &part,
					 const stk::mesh::EntityRank rank,
					 Ioss::GroupingEntity *io_entity,
					 std::map<stk::mesh::FieldBase *, const stk::io::MeshField *> *missing_fields = nullptr);

      void build_field_part_associations_for_part(Ioss::Region *region,
                                                  const stk::mesh::FieldBase *f,
                                                  const stk::mesh::Part * part,
                                                  stk::io::MeshField &mf);

      DatabasePurpose m_db_purpose;
      Teuchos::RCP<Ioss::DatabaseIO> m_database;
      Teuchos::RCP<Ioss::Region> m_region;
      std::vector<stk::io::MeshField> m_fields;

      /*@{*/

      /**
       * The 'startupTime' and 'periodLength' are used to support input of periodic data.
       * The 'startupTime' specifies the length of time prior to the start of
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
       *      t_db = t_app
       *  \endcode
       * Else:
       *   If: ptype == CYCLIC
       *    \code
       *      t_db = t_startup + mod(t_app-t_startup, period_length)
       *    \endcode
       *   If: ptype == REVERSING
       *    \code
       *      tpm = mod(t_app-t_startup, 2*period_length)
       *      if (tpm <= period_length)
       *         t_db = t_startup + tpm
       *      else 
       *         t_db = t_startup + (2*period_length - tpm)
       *    \endcode
       *
       * \code
       *  t_db = t_db * t_scale + t_offset
       * \endcode
       *  
       */
      double m_startupTime;

      /** See InputFile::startupTime */
      double m_periodLength;
      /** See InputFile::startupTime */
      double m_scaleTime;
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
      std::vector<std::string>* m_multiStateSuffixes = nullptr;

    private:
      InputFile(const InputFile &);
      const InputFile & operator=(const InputFile &);
    };
  } // namespace io
} // namespace stk
#endif
