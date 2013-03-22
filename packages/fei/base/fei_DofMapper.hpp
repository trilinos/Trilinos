/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_DofMapper_hpp_
#define _fei_DofMapper_hpp_

#include <fei_macros.hpp>

#include <set>
#include <map>
#include <sstream>
#include <stdexcept>

#include <fei_Dof.hpp>

namespace fei {

/** A mapping from mesh-degrees-of-freedom (dofs) to equation-numbers.
 *
 * Mesh-dofs are represented by fei::Dof, see fei_Dof.hpp for details. 
 *
 * Equation numbers are also called global indices. They are globally-unique
 * and zero-based.
 *
 * A 'dof' may correspond to multiple global-indices if the dof's field
 * has multiple scalar components, e.g. a vector field such as velocity in 3D.
 * - The index returned by getGlobalIndex is the eqn-number for the first
 *   component of the dof's field.
 * - The return-value of getDof is a pair<Dof,component> where component 
 *   indicates which component of the dof's field corresponds to the input
 *   global-index.
 *
 * Fields are assumed to be scalar fields (have 1 component) unless a field-size
 * is set using the setFieldSize method.
 */
template<class LocalOrdinal, class GlobalOrdinal, class DofOrder=less_rank_id_field<LocalOrdinal, GlobalOrdinal> >
class DofMapper {
 public:
  /** constructor */
  DofMapper()
  : m_dof_idx(), m_idx_dof(), m_maps_are_valid(false), m_field_sizes() {}

  /** destructor */
  ~DofMapper() {}

  void addDOF(LocalOrdinal rank, GlobalOrdinal id, LocalOrdinal field)
  {
	if (m_field_sizes.find(field) == m_field_sizes.end()) m_field_sizes.insert(std::make_pair(field,1));
    //m_maps_are_valid is false when a new Dof is inserted.
    m_maps_are_valid =
      m_dof_idx.insert(std::make_pair(Dof<LocalOrdinal,GlobalOrdinal>(rank, id, field), 0)).second;
  }

  /** Set the specified field to have the specified field_size.
   * 'field' is added to the internal field map if not already present.
   * If 'field' is already present, its field_size is reset to the new value.
   */
  void setFieldSize(LocalOrdinal field, LocalOrdinal field_size);

  LocalOrdinal getFieldSize(LocalOrdinal field) const;

  GlobalOrdinal getGlobalIndex(LocalOrdinal rank, GlobalOrdinal id, LocalOrdinal field) const;

  std::pair<const Dof<LocalOrdinal,GlobalOrdinal>*,LocalOrdinal> getDof(GlobalOrdinal global_index) const;

  bool maps_are_valid() const { return m_maps_are_valid; }
  void set_maps_are_valid(bool flag) { m_maps_are_valid = flag; }

  typedef typename std::map<Dof<LocalOrdinal,GlobalOrdinal>,GlobalOrdinal,DofOrder> DofMap;

  typename DofMap::const_iterator begin_dof() const
  { return m_dof_idx.begin(); }

  typename DofMap::const_iterator end_dof() const
  { return m_dof_idx.end(); }

  typename DofMap::iterator begin_dof()
  { return m_dof_idx.begin(); }

  typename DofMap::iterator end_dof()
  { return m_dof_idx.end(); }

  typedef typename std::map<GlobalOrdinal,const Dof<LocalOrdinal,GlobalOrdinal>*> IdxMap;

  typename IdxMap::const_iterator begin_idx() const
  { return m_idx_dof.begin(); }

  typename IdxMap::const_iterator end_idx() const
  { return m_idx_dof.end(); }

  typename IdxMap::iterator begin_idx()
  { return m_idx_dof.begin(); }

  typename IdxMap::iterator end_idx()
  { return m_idx_dof.end(); }

  const DofMap& get_dof_idx_map() const {return m_dof_idx;}
  DofMap& get_dof_idx_map() {return m_dof_idx;}

  const IdxMap& get_idx_dof_map() const {return m_idx_dof;}
  IdxMap& get_idx_dof_map() {return m_idx_dof;}

  typedef typename std::map<LocalOrdinal,LocalOrdinal> FieldSizeMap;
  const FieldSizeMap& getFieldSizeMap() const {return m_field_sizes;}

 private:
  std::map<Dof<LocalOrdinal, GlobalOrdinal>, GlobalOrdinal, DofOrder > m_dof_idx;

  std::map<GlobalOrdinal, const Dof<LocalOrdinal, GlobalOrdinal>*> m_idx_dof;
  bool m_maps_are_valid;

  std::map<LocalOrdinal,LocalOrdinal> m_field_sizes;

  DofMapper(const DofMapper<LocalOrdinal,GlobalOrdinal>& src);
  DofMapper& operator=(const DofMapper<LocalOrdinal,GlobalOrdinal>& src);
};//class DofMapper

template<class LocalOrdinal,class GlobalOrdinal,class DofOrder>
void DofMapper<LocalOrdinal,GlobalOrdinal,DofOrder>::setFieldSize(LocalOrdinal field, LocalOrdinal field_size)
{
	typename FieldSizeMap::iterator f_iter = m_field_sizes.find(field);
	if (f_iter == m_field_sizes.end()) {
	  m_field_sizes.insert(std::make_pair(field, field_size));
	}
	else {
	  //field already present, resetting field_size:
	  f_iter->second = field_size;
	}
}

template<class LocalOrdinal,class GlobalOrdinal,class DofOrder>
LocalOrdinal DofMapper<LocalOrdinal,GlobalOrdinal,DofOrder>::getFieldSize(LocalOrdinal field) const
{
  typename FieldSizeMap::const_iterator f_iter = m_field_sizes.find(field);
  if (f_iter == m_field_sizes.end()) {
	  std::ostringstream os;
	  os << "fei::DofMapper::getFieldSize ERROR, field=="
		  << field << " not found";
	  std::string str = os.str();
	  throw std::runtime_error(str);
  }
  return f_iter->second;
}

template<class LocalOrdinal,class GlobalOrdinal,class DofOrder>
GlobalOrdinal DofMapper<LocalOrdinal,GlobalOrdinal,DofOrder>::getGlobalIndex(LocalOrdinal rank, GlobalOrdinal id, LocalOrdinal field) const
{
  typename DofMap::const_iterator iter = m_dof_idx.find(Dof<LocalOrdinal,GlobalOrdinal>(rank,id,field));
  if (iter == m_dof_idx.end()) {
    std::ostringstream osstr;
    osstr << "fei::DofMapper::getGlobalIndex ERROR, dof("
        << rank << "," << id << "," << field << ") not found.";
    std::string str = osstr.str();
    throw std::runtime_error(str);
  }

  return iter->second;
}

template<class LocalOrdinal,class GlobalOrdinal,class DofOrder>
std::pair<const Dof<LocalOrdinal,GlobalOrdinal>*,LocalOrdinal>
DofMapper<LocalOrdinal,GlobalOrdinal,DofOrder>::getDof(GlobalOrdinal global_index) const
{
  typename IdxMap::const_iterator iter = m_idx_dof.lower_bound(global_index);
  if (iter == m_idx_dof.begin()) {
	if (iter->first == global_index) {
	  return std::make_pair(iter->second, global_index);
	}
	else {
      std::ostringstream osstr;
      osstr << "fei::DofMapper::getDof ERROR, dof not found for global_index=="
        << global_index;
      std::string str = osstr.str();
      throw std::runtime_error(str);
	}
  }
  else if (iter != m_idx_dof.end() && iter->first == global_index) {
    //return pair(dof,component-of-field)
	return std::make_pair(iter->second, 0);
  }

  bool last_dof = iter == m_idx_dof.end();
  --iter;
  //return pair(dof,component-of-field)
  LocalOrdinal component = global_index - iter->first;
  bool check_range_of_component = last_dof && !m_field_sizes.empty();
  if (check_range_of_component) {
	typename std::map<LocalOrdinal,LocalOrdinal>::const_iterator f_iter = m_field_sizes.find(iter->second->field());
	if (f_iter == m_field_sizes.end() || f_iter->second <= component) {
	  std::ostringstream os;
	  os << "fei::DofMapper::getDof ERROR2, dof not found for global_index=="
		  << global_index;
	  std::string str = os.str();
	  throw std::runtime_error(str);
	}
  }

  return std::make_pair(iter->second, component);
}

}//namespace fei

#endif

