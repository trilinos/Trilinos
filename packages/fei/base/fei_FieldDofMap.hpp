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
#include <fei_constants.hpp>
#include <sstream>

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

  void add_field(LocalOrdinal fieldID, LocalOrdinal fieldSize, LocalOrdinal fieldType=fei::UNKNOWN);

  LocalOrdinal get_dof_id(LocalOrdinal fieldID, LocalOrdinal offset);

 private:
  void compute_dof_ids();

  //dof_id_map maps fields to pairs which each pair is <fieldSize, dof_id>.
  typedef std::map<LocalOrdinal,std::pair<LocalOrdinal,LocalOrdinal> > dof_id_map;

  dof_id_map m_dof_id_map;
  bool m_need_to_compute_dof_ids;
};//class FieldDofMap

template<class LocalOrdinal>
void FieldDofMap<LocalOrdinal>::add_field(LocalOrdinal fieldID, LocalOrdinal fieldSize, LocalOrdinal fieldType)
{
  //initially store 0 for the dof_id, will compute dof_ids later.
  m_dof_id_map.insert(std::make_pair(fieldID,std::make_pair(fieldSize,fieldType)));
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

  //first make sure dof_ids are unique. so make a second map which
  //maps dof_ids to iterators into the first map, and watch for
  //duplicates as we fill the second map. When we find duplicates,
  //we'll shift dof_ids up as necessary. (throw an exception if the
  //duplicate is not fei::UNKNOWN or bigger.)

  typedef std::map<LocalOrdinal, typename dof_id_map::iterator> dof_iter_map;
  dof_iter_map dof_2_iter;

  typename dof_id_map::iterator
    iter = m_dof_id_map.begin(), iter_end = m_dof_id_map.end();

  for(; iter!=iter_end; ++iter) {
    LocalOrdinal this_dof_id = iter->second.second;

    typename dof_iter_map::iterator di_iter = dof_2_iter.find(this_dof_id);

    if (di_iter != dof_2_iter.end()) {
      if (this_dof_id < fei::UNKNOWN) {
        std::ostringstream osstr;
        osstr << "fei::FieldDofMap::compute_dof_ids ERROR, duplicate field types found (";
        osstr << this_dof_id << " used more than once.)";
        std::string str = osstr.str();
        throw std::runtime_error(str);
      }

      //now we need to get fieldSize, and this is a little ugly:
      //di_iter->second is the iterator to the other map.
      //so di_iter->second->second is the value in the other map, which is a pair.
      //so di_iter->second->second.first is the fieldSize that we want.
      std::pair<LocalOrdinal,LocalOrdinal>& fsize_and_dof = di_iter->second->second;
      LocalOrdinal fieldSize = fsize_and_dof.first;
      LocalOrdinal last_dof_id = fsize_and_dof.second;
      di_iter->second = iter;
      iter->second.second = last_dof_id + fieldSize;
    }
    else dof_2_iter.insert(std::make_pair(this_dof_id, iter));
  }

  m_need_to_compute_dof_ids = false;
}

}//namespace fei

#endif

