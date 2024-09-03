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

#ifndef STK_IO_MeshField_h
#define STK_IO_MeshField_h

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <cstddef>                  // for size_t
#include <string>                   // for string
#include <vector>                   // for vector
#include "stk_mesh/base/Types.hpp"  // for EntityRank
#include "stk_mesh/base/Entity.hpp"
namespace Ioss { class GroupingEntity; }
namespace Ioss { class Region; }
namespace stk { namespace io { class DBStepTimeInterval; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
namespace stk { namespace io { class InputFile; } }

namespace stk {
namespace io {

class MeshFieldPart {
public:
  MeshFieldPart(stk::mesh::EntityRank rank, const stk::mesh::Part *part,
                Ioss::GroupingEntity *io_entity, const std::string db_field_name)
    : m_rank(rank), m_stkPart(part), m_ioEntity(io_entity), m_dbName(db_field_name),
      m_preStep(0), m_postStep(0)
  {}

  void get_interpolated_field_data(const DBStepTimeInterval &sti, std::vector<double> &values);
  void release_field_data();

  stk::mesh::EntityRank get_entity_rank() const {return m_rank;}
  Ioss::GroupingEntity* get_io_entity() const {return m_ioEntity;}
  const stk::mesh::Part* get_stk_part() const {return m_stkPart;}
  void set_cached_entity_list(std::vector<stk::mesh::Entity> && entityList) { m_cachedEntityList = std::move(entityList); }
  const std::vector<stk::mesh::Entity> & get_cached_entity_list() const { return m_cachedEntityList; }

private:
  void load_field_data(const DBStepTimeInterval &sti);

  stk::mesh::EntityRank m_rank;
  const stk::mesh::Part *m_stkPart;
  Ioss::GroupingEntity  *m_ioEntity;
  std::string m_dbName;
  std::vector<double> m_preData;
  std::vector<double> m_postData;
  std::vector<stk::mesh::Entity> m_cachedEntityList;
  size_t m_preStep;
  size_t m_postStep;
};

class MeshField
{
public:

  friend class InputFile;

  // Options:
  // * Frequency:
  //   -- one time only
  //   -- multiple times
  // * Matching of time
  //   -- Linear Interpolation
  //   -- Closest
  //   -- specified time

  enum TimeMatchOption {
    LINEAR_INTERPOLATION,
    CLOSEST,
    SPECIFIED }; // Use time specified on MeshField

  // Read 'db_name' field data into 'field' using 'tmo' (default CLOSEST) time on database.
  // Analysis time will be mapped to db time.
  MeshField(stk::mesh::FieldBase *field,
            const std::string &db_name="",
            TimeMatchOption tmo = CLOSEST);
  MeshField(stk::mesh::FieldBase &field,
            const std::string &db_name="",
            TimeMatchOption tmo = CLOSEST);

  ~MeshField();

  // MeshField(const MeshField&); Default version is good.
  // MeshField& operator=(const MeshField&); Default version is good

  MeshField& set_read_time(double time_to_read);
  MeshField& set_active();
  MeshField& set_inactive();
  MeshField& set_single_state(bool yesno);
  MeshField& set_read_once(bool yesno);
  MeshField& set_classic_restart();

  double get_read_time() const {return m_timeToRead;}

  // Limit the field to part(s) specified by this call.
  // Default is to restore field on all parts that it is defined on.
  MeshField &add_subset(const stk::mesh::Part &part);

  bool is_active() const {return m_isActive;}

  // Returns the time at which the field data was restored.
  // Either the closest time on the database, or the interopolated
  // time.
  double restore_field_data(stk::mesh::BulkData &bulk,
                            const DBStepTimeInterval &sti,
                            bool ignore_missing_fields = false,
                            std::vector<std::string>* multiStateSuffixes=nullptr);

  void fill_entity_list_cache(const stk::mesh::BulkData &bulk);

  double restore_field_data_at_step(Ioss::Region *region,
                                    stk::mesh::BulkData &bulk,
                                    int step,
                                    bool ignore_missing_fields = false,
                                    std::vector<std::string>* multiStateSuffixes = nullptr,
                                    bool useEntityListCache = false);

  const std::string &db_name() const {return m_dbName;}
  stk::mesh::FieldBase *field() const {return m_field;}


  void add_part(const stk::mesh::EntityRank rank,
                const stk::mesh::Part &part,
                Ioss::GroupingEntity *io_entity);

  bool operator==(const MeshField &other) const;

  bool field_restored() const {return m_fieldRestored;}
  double time_restored() const {return m_timeRestored;}

private:
  MeshField();

  std::vector<const stk::mesh::Part*> m_subsetParts;
  std::vector<MeshFieldPart> m_fieldParts;

  stk::mesh::FieldBase *m_field;
  std::string m_dbName; ///<  Name of the field on the input/output database.

  double m_timeToRead;

  TimeMatchOption m_timeMatch;
  bool m_oneTimeOnly;
  bool m_singleState;
  bool m_isActive;

  bool m_fieldRestored;
  double m_timeRestored;
};
}
} 
#endif
