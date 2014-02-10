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
    class Interpolator;
    
    struct MeshFieldPart {
      MeshFieldPart(const stk::mesh::Part &part, stk::mesh::EntityRank rank, Ioss::GroupingEntity *io_entity)
	: m_ordinal(part.mesh_meta_data_ordinal()), m_rank(rank), m_ioEntity(io_entity)
      {}

      unsigned m_ordinal;
      stk::mesh::EntityRank m_rank;
      Ioss::GroupingEntity *m_ioEntity;
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

      enum PeriodType {
	CYCLIC,   /*< Cycles 0 1 2 0 1 2 0 1 2 ... */
	REVERSING /*< Cycles 0 1 2 1 0 1 2 1 0 ... */
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

      /**
       * The MeshField class uses the "Named Parameter Idiom" which supports "method chaining".
       * For example, you can do the following:
       *
       * MeshField my_field = MeshField(stress_fld, "stress")
       *                               .set_read_time(1.0)
       *                               .set_offset_time(2.0)
       *                               .set_periodic_time(11.5)
       *                               .set_start_time(0.5)
       *                               .set_stop_time(12.3);
       */

      MeshField& set_read_time(double time_to_read);
      MeshField& set_offset_time(double offset_time);
      MeshField& set_periodic_time(double period_length, double startup_time = 0.0,
				   PeriodType ptype = CYCLIC);
      MeshField& set_start_time(double start_time);
      MeshField& set_stop_time(double stop_time);
  
      void restore_field_data(Ioss::Region *region,
			      stk::mesh::BulkData &bulk,
			      const stk::io::DBStepTimeInterval &sti);
      
      void get_data(double time, bool use_state_n = false);
      void initialize(); // Initialize for interpolation role.

      const std::string &db_name() const {return m_dbName;}
      stk::mesh::FieldBase *field() const {return m_field;}
	
      void add_part(const stk::mesh::Part &part,
		    const stk::mesh::EntityRank rank,
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
    private:
      double m_time;
      double m_timeRead;
      
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

      /** See MeshField::startupTime */
      double m_periodLength;
      /** See MeshField::startupTime */
      double m_offsetTime;
      /** See MeshField::startupTime */
      double m_startTime;
      /** See MeshField::startupTime */
      double m_stopTime;
      /** See MeshField::startupTime */
      PeriodType m_periodType;
      
      /*@}*/
      Interpolator *m_interpolator;
    };
  } 
} 
#endif
