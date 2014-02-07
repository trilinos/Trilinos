/*--------------------------------------------------------------------*/
/*    Copyright 2004, 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <stk_io/MeshField.hpp>
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, etc

#include <limits>
#include <math.h>
#include <assert.h>

namespace stk {
  namespace io {

MeshField::MeshField(stk::mesh::FieldBase *field,
		     const std::string &db_name,
		     TimeMatchOption tmo)
  : m_field(field),
    m_dbName(db_name),
    m_timeMap(ANALYSIS),
    m_timeMatch(tmo),
    m_oneTimeOnly(false),
    m_wasFound(false),
    m_time(0.0),
    m_timeRead(0.0),
    m_startupTime(0.0),
    m_periodLength(0.0),
    m_offsetTime(0.0),
    m_startTime(-std::numeric_limits<double>::max()),
    m_stopTime(std::numeric_limits<double>::max()),
    m_periodType(CYCLIC),
    m_interpolator(NULL)
{
  if (db_name == "") {
    m_dbName = field->name();
  }
}
MeshField::MeshField(stk::mesh::FieldBase &field,
		     const std::string &db_name,
		     TimeMatchOption tmo)
  : m_field(&field),
    m_dbName(db_name),
    m_timeMap(ANALYSIS),
    m_timeMatch(tmo),
    m_oneTimeOnly(false),
    m_wasFound(false),
    m_time(0.0),
    m_timeRead(0.0),
    m_startupTime(0.0),
    m_periodLength(0.0),
    m_offsetTime(0.0),
    m_startTime(-std::numeric_limits<double>::max()),
    m_stopTime(std::numeric_limits<double>::max()),
    m_periodType(CYCLIC),
    m_interpolator(NULL)
{
  if (db_name == "") {
    m_dbName = field.name();
  }
}

MeshField::~MeshField()
{
#if 0
  delete m_interpolator;
  m_interpolator = NULL;
#endif
}

#if 0
MeshField& MeshField::set_field(const std::string &entity_name,
				const std::string &field_name,
				const std::string &fmwk_field_name,
				UInt type)
{
  assert(entity_name     != "");
  assert(field_name      != "");
  assert(fmwk_field_name != "");
  entityName = entity_name;
  fieldName = field_name;
  fmwkFieldName = fmwk_field_name;
  type_ = type;
  return *this;
}
#endif

MeshField& MeshField::set_read_time(double time_to_read)
{
  m_time = time_to_read;
  return *this;
}

MeshField& MeshField::set_offset_time(double offset_time)
{
  m_offsetTime = offset_time;
  return *this;
}

MeshField& MeshField::set_periodic_time(double period_length, double startup_time, PeriodType ptype)
{
  m_periodLength = period_length;
  m_startupTime = startup_time;
  m_periodType = ptype;
  return *this;
}

MeshField& MeshField::set_start_time(double start_time)
{
  m_startTime = start_time;
  return *this;
}
  
MeshField& MeshField::set_stop_time(double stop_time)
{
  m_stopTime = stop_time;
  return *this;
}
  
void MeshField::get_data(double time, bool use_state_n)
{
#if 0
  if (role_ == Ioss::FemIO::Bridge::COPY_TRANSFER_PROXIMITY ||
      role_ == Ioss::FemIO::Bridge::COPY_TRANSFER_BY_ID ||
      role_ == Ioss::FemIO::Bridge::THREE_STATE)
    return;
  
  if (m_interpolator == NULL)
    m_interpolator = InterpolatorFactory::create(this);
  assert(m_interpolator != NULL);

  if (role_ == Ioss::FemIO::Bridge::TRANSIENT ||
      role_ == Ioss::FemIO::Bridge::CLOSEST_TIME_ONCE ||
      role_ == Ioss::FemIO::Bridge::MODEL) {
    time = m_time;
    m_interpolator->get_data(time, fieldState_);
    if (use_state_n && fieldState_ == sierra::STATE_NP1) {
      // If initializing (use_state_n == true), then also load data
      // into STATE_N if this is a two state field.
      sierra::State state = sierra::STATE_N;
      m_interpolator->get_data(time, state);
    }
  }

  else if (role_ == Ioss::FemIO::Bridge::FIXED_TIME_INTERPOLATION) {
    time = m_time;
    m_interpolator->get_data(time, fieldState_);

  }  else if (role_ == Ioss::FemIO::Bridge::TIME_INTERPOLATION ||
	      role_ == Ioss::FemIO::Bridge::CLOSEST_TIME) {
    // See if 'time' is outside active range 
    if (time < startTime || time > stopTime)
      return;
    
    // Modify input time if startupTime and periodLength specified
    if (time > startupTime && periodLength > 0.0) {

      if (periodType == CYCLIC) {
	time = startupTime + fmod((time-startupTime), periodLength);

      } else if (periodType == REVERSING) {
	double pmod = fmod((time-startupTime), 2.0*periodLength);
	if (pmod <= periodLength)
	  time = startupTime + pmod;
	else
	  time = startupTime + 2.0 * periodLength - pmod;
      }	  
    }
    time += offsetTime;
    
    m_interpolator->get_data(time, fieldState_);
    if (use_state_n && fieldState_ == sierra::STATE_NP1) {
      // If initializing (use_state_n == true), then also load data
      // into STATE_N if this is a two state field.
      sierra::State state = sierra::STATE_N;
      m_interpolator->get_data(time, state);
    }
  }
#endif
}

}}  // Close namespaces
