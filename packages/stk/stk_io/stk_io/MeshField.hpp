/*--------------------------------------------------------------------*/
/*    Copyright (c) 2013, Sandia Corporation.
/*    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*    the U.S. Governement retains certain rights in this software.
/*    
/*    Redistribution and use in source and binary forms, with or without
/*    modification, are permitted provided that the following conditions are
/*    met:
/*    
/*        * Redistributions of source code must retain the above copyright
/*          notice, this list of conditions and the following disclaimer.
/*    
/*        * Redistributions in binary form must reproduce the above
/*          copyright notice, this list of conditions and the following
/*          disclaimer in the documentation and/or other materials provided
/*          with the distribution.
/*    
/*        * Neither the name of Sandia Corporation nor the names of its
/*          contributors may be used to endorse or promote products derived
/*          from this software without specific prior written permission.
/*    
/*    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*    
/*--------------------------------------------------------------------*/

#ifndef STK_IO_MeshField_h
#define STK_IO_MeshField_h

#include <stddef.h>                     // for size_t
#include <string>                       // for string, basic_string
#include <vector>                       // for vector
#include "stk_mesh/base/Types.hpp"      // for EntityRank
namespace Ioss { class GroupingEntity; }
namespace stk { namespace io { class DBStepTimeInterval; } }
namespace stk { namespace io { class InputFile; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class Part; } }

namespace Ioss {
  class Region;
}

namespace stk {
  namespace mesh {
  }
  
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

    private:
      void load_field_data(const DBStepTimeInterval &sti);

      stk::mesh::EntityRank m_rank;
      const stk::mesh::Part *m_stkPart;
      Ioss::GroupingEntity  *m_ioEntity;
      std::string m_dbName;
      std::vector<double> m_preData;
      std::vector<double> m_postData;
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

      MeshField();

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
      
      // Limit the field to part(s) specified by this call.
      // Default is to restore field on all parts that it is defined on.
      MeshField &add_subset(const stk::mesh::Part &part);

      bool is_active() const {return m_isActive;}
      
      void restore_field_data(stk::mesh::BulkData &bulk,
			      const DBStepTimeInterval &sti);
      
      const std::string &db_name() const {return m_dbName;}
      stk::mesh::FieldBase *field() const {return m_field;}
	

      void add_part(const stk::mesh::EntityRank rank,
		    const stk::mesh::Part &part,
		    Ioss::GroupingEntity *io_entity);
      
      bool operator==(const MeshField &other) const;

    private:
      std::vector<const stk::mesh::Part*> m_subsetParts;
      std::vector<MeshFieldPart> m_fieldParts;
      
      stk::mesh::FieldBase *m_field;
      std::string m_dbName; ///<  Name of the field on the input/output database.

      double m_timeToRead;
      
      TimeMatchOption m_timeMatch;
      bool m_oneTimeOnly;
      bool m_singleState;
      bool m_isActive;
    };
  } 
} 
#endif
