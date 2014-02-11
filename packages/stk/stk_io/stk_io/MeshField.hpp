/*--------------------------------------------------------------------*/
/*    Copyright 2014 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_IO_MeshField_h
#define STK_IO_MeshField_h

#include <string>
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/Part.hpp"  

namespace Ioss {
  class GroupingEntity;
  class Region;
}

namespace stk {
  namespace mesh {
    class FieldBase;
  }
  
  namespace io {
    class DBStepTimeInterval;
    
    class MeshFieldPart {
    public:
      MeshFieldPart(stk::mesh::EntityRank rank, Ioss::GroupingEntity *io_entity, const std::string db_field_name)
	: m_rank(rank), m_ioEntity(io_entity), m_dbName(db_field_name),
	  m_preStep(0), m_postStep(0)
      {}

      void get_interpolated_field_data(const DBStepTimeInterval &sti, std::vector<double> &values);
      void release_field_data();
      
    private:
      void load_field_data(const DBStepTimeInterval &sti);

    public:
      stk::mesh::EntityRank m_rank;
      Ioss::GroupingEntity *m_ioEntity;

    private:
      std::string m_dbName;

      std::vector<double> m_preData;
      std::vector<double> m_postData;
      size_t m_preStep;
      size_t m_postStep;
    };
    
    class MeshField
    {
    public:

      // Options:
      // * Frequency:
      //   -- one time only
      //   -- multiple times
      // * Database Time:
      //   -- non-transient variable
      //   -- analysis time (Offset/Period/... from analysis time)
      //   -- specified time
      // * Matching of time
      //   -- Linear Interpolation
      //   -- Closest
      //   -- Less Equal
      //   -- Greater Equal
      
      enum DatabaseTimeMapping {
	MODEL,      // Non-transient variable on db
	ANALYSIS,   // Use analysis time as db time (with offset,period,...)
	SPECIFIED,  // Use time specified on MeshField
      };

      enum TimeMatchOption {
	LINEAR_INTERPOLATION,
	CLOSEST, 
	LESS_EQUAL,
	GREATER_EQUAL,
      };

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
      void restore_field_data(stk::mesh::BulkData &bulk,
			      const DBStepTimeInterval &sti);
      
      const std::string &db_name() const {return m_dbName;}
      stk::mesh::FieldBase *field() const {return m_field;}
	
      void add_part(const stk::mesh::EntityRank rank,
		    Ioss::GroupingEntity *io_entity);
      
    private:
      std::vector<MeshFieldPart> m_fieldParts;
      
      stk::mesh::FieldBase *m_field;
      std::string m_dbName; ///<  Name of the field on the input/output database.

      DatabaseTimeMapping m_timeMap;
    public:
      TimeMatchOption m_timeMatch;
    private:
      bool m_oneTimeOnly;
    public:
      bool m_singleState;
      bool m_wasFound;
    };
  } 
} 
#endif
