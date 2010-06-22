/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_FieldDofMap_hpp_
#define _fei_FieldDofMap_hpp_

#include <fei_macros.hpp>

namespace fei {

/** A simple class to map fieldIDs to dof_ids.
  fieldIDs are arbitrary numbers, while dof_ids are members of a
  zero-based contiguous set, one dof_id for each scalar component of each field.
*/
template<class LocalOrdinal>
class FieldDofMap {
 public:
  FieldDofMap()
    : m_dof_id_map(), m_need_to_compute_dof_ids(true)
  {}

  ~FieldDofMap(){}

  void add_field(LocalOrdinal fieldID, LocalOrdinal fieldSize);

  LocalOrdinal get_dof_id(LocalOrdinal fieldID, LocalOrdinal offset);

 private:
  void compute_dof_ids();

  //dof_id_map maps fields to pairs which each pair is <fieldSize, dof_id>.
  typedef std::map<LocalOrdinal,std::pair<LocalOrdinal,LocalOrdinal> > dof_id_map;

  dof_id_map m_dof_id_map;
  bool m_need_to_compute_dof_ids;
};//class FieldDofMap

template<class LocalOrdinal>
void FieldDofMap<LocalOrdinal>::add_field(LocalOrdinal fieldID, LocalOrdinal fieldSize)
{
  //initially store 0 for the dof_id, will compute dof_ids later.
  m_dof_id_map.insert(std::make_pair(fieldID,std::make_pair(fieldSize,0)));
  m_need_to_compute_dof_ids = true;
}

template<class LocalOrdinal>
LocalOrdinal FieldDofMap<LocalOrdinal>::get_dof_id(LocalOrdinal fieldID, LocalOrdinal offset)
{
  if (m_need_to_compute_dof_ids) {
    compute_dof_ids();
  }

  typename dof_id_map::const_iterator
    iter = m_dof_id_map.find(fieldID);
  if (iter == m_dof_id_map.end()) {
    throw std::runtime_error("fei::FieldDofMap ERROR, specified fieldID not found.");
  }

  if (offset >= iter->second.first) {
    throw std::runtime_error("FieldDofMap::get_dof_id ERROR, specified offset is greater than field-size.");
  }

  //dof_id_map maps fieldID to a pair<fieldSize,dof_id>.
  //We want to return the dof_id plus the input offset:
  return iter->second.second + offset;
}

template<class LocalOrdinal>
void FieldDofMap<LocalOrdinal>::compute_dof_ids()
{
  if (!m_need_to_compute_dof_ids) return;

  LocalOrdinal dof_id = 0;
  typename dof_id_map::iterator
    iter = m_dof_id_map.begin(), iter_end = m_dof_id_map.end();

  for(; iter!=iter_end; ++iter) {
    LocalOrdinal fieldSize = iter->second.first;
    iter->second.second = dof_id;
    dof_id += fieldSize;
  }

  m_need_to_compute_dof_ids = false;
}

}//namespace fei

#endif

