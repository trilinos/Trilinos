// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
      InputFile(std::shared_ptr<Ioss::Region> ioss_input_region);

      ~InputFile()
      {delete m_multiStateSuffixes;}

      void create_ioss_region();
      FieldNameToPartVector get_var_names(Ioss::EntityType type, const stk::mesh::MetaData& meta);
      void add_input_field(const stk::io::MeshField &mesh_field);
      void add_all_mesh_fields_as_input_fields(stk::mesh::MetaData &meta, MeshField::TimeMatchOption tmo);
      bool read_input_field(stk::io::MeshField &mf, stk::mesh::BulkData &bulk);
      double read_defined_input_fields(double time, std::vector<stk::io::MeshField> *missingFields,
				       stk::mesh::BulkData &bulk);
      double read_defined_input_fields(int step, std::vector<stk::io::MeshField> *missingFields,
				       stk::mesh::BulkData &bulk);
      double read_defined_input_fields_at_step(int step, std::vector<stk::io::MeshField> *missingFields,
                                               stk::mesh::BulkData &bulk, bool useEntityListCache = false);
      void get_global_variable_names(std::vector<std::string> &names);

      std::shared_ptr<Ioss::Region> get_input_ioss_region()
      {
	      if (m_region.get() == nullptr && m_database.get() != nullptr) {
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
          if(m_database.get() != nullptr) {
              m_database->set_surface_split_type(split_type);
          }
      }
      Ioss::SurfaceSplitType get_surface_split_type() const {
          if(m_database.get() != nullptr) {
              return m_database->get_surface_split_type();
          }

          return Ioss::SPLIT_INVALID;
      }

      std::shared_ptr<Ioss::DatabaseIO> get_ioss_input_database() const
      {
	      return m_database;
      }

      bool set_multistate_suffixes(const std::vector<std::string>& multiStateSuffixes)
      {
          if(nullptr != m_multiStateSuffixes) {
              delete m_multiStateSuffixes;
              m_multiStateSuffixes = nullptr;
          }

          m_multiStateSuffixes = new std::vector<std::string>(multiStateSuffixes);
          return true;
      }

      const std::vector<std::string>& get_multistate_suffixes() const
      {
        static std::vector<std::string> emptyVector;

        if(nullptr != m_multiStateSuffixes) {
          return *m_multiStateSuffixes;
        }

        return emptyVector;
      }

      DatabasePurpose get_database_purpose() const { return m_db_purpose; }

      void initialize_input_fields();

      // Query the database for the maximum time. If useCache is true, then it will
      // return the cached value. The only exception to this is if the database was
      // opened in APPEND mode where it is possible to open the file for both read
      // and write with possibly changing max times. In this case, ignore useCache
      // and always query the database
      double get_max_time(bool useCache = true);

    private:

      DatabasePurpose m_db_purpose;
      std::shared_ptr<Ioss::DatabaseIO> m_database;
      std::shared_ptr<Ioss::Region> m_region;
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
  double m_maxTime;
  PeriodType m_periodType;

  /*@}*/

public:
  bool m_fieldsInitialized;

private:
  bool m_haveCachedEntityList = false;
  std::vector<std::string>* m_multiStateSuffixes = nullptr;


private:
  InputFile(const InputFile &);
  const InputFile & operator=(const InputFile &);
};
} // namespace io
} // namespace stk
#endif
