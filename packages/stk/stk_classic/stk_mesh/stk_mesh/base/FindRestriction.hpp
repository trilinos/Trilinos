/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_base_FindRestriction_hpp
#define stk_mesh_base_FindRestriction_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldRestriction.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace stk {
namespace mesh {

namespace {
inline
const FieldBase::Restriction & empty_field_restriction()
{
  //NOT THREAD SAFE
  static const FieldBase::Restriction empty ;
  return empty ;
}
}//namespace anonymous

//Given a field and a range of parts, determine whether the field has a restriction
//for the part range. (Common usage is to provide the part-range from a bucket; i.e.,
//determine whether the field should be allocated in the bucket.)
template<typename PartIterator, class Compare>
const FieldBase::Restriction& find_restriction(const FieldBase& field,
                                               EntityRank erank,
                                               PartIterator pbegin,
                                               PartIterator pend,
                                               Compare comp)
{
  GetPartIterOrdinal<PartIterator> get_part_ordinal;

  const FieldBase::Restriction & empty = empty_field_restriction();
  const FieldBase::Restriction * restriction = & empty ;

  const std::vector<FieldBase::Restriction> & restr_vec = field.restrictions();
  const std::vector<FieldBase::Restriction>::const_iterator iend = restr_vec.end();
        std::vector<FieldBase::Restriction>::const_iterator ibeg = restr_vec.begin();

  //NOT THREAD SAFE
  static FieldRestriction restr;

  for ( PartIterator it = pbegin; it != pend && iend != ibeg ; ++it ) {

    restr.set_entity_rank(erank);
    restr.set_part_ordinal(get_part_ordinal(it));

    //lower_bound returns an iterator to either the insertion point for the
    //'restr' argument, or to a matching restriction.
    //It only returns the 'end' iterator if 'restr' is past the end of the
    //vector of restrictions being searched.
    //This depends on the input part ordinals being sorted, and on the restriction
    //vector being sorted by part ordinal.

    ibeg = std::lower_bound( ibeg , iend , restr );

    if ( (iend != ibeg) && (*ibeg == restr) ) {
      if ( restriction == & empty ) { restriction = & *ibeg ; }

      if ( ibeg->not_equal_stride(*restriction) ) {

        Part & p_old = MetaData::get(field).get_part( ibeg->part_ordinal() );
        Part & p_new = MetaData::get(field).get_part( restriction->part_ordinal() );

        std::ostringstream msg ;
        msg << " FAILED WITH INCOMPATIBLE DIMENSIONS FOR " ;
        msg << field ;
        msg << " Part[" << p_old.name() ;
        msg << "] and Part[" << p_new.name() ;
        msg << "]" ;

        ThrowErrorMsg( msg.str() );
      }
    }
  }

  const std::vector<FieldBase::Restriction> & sel_res = field.selector_restrictions();
  std::pair<PartIterator,PartIterator> part_range = std::make_pair(pbegin, pend);
  for(std::vector<FieldBase::Restriction>::const_iterator it=sel_res.begin(), it_end=sel_res.end(); it != it_end; ++it) {
    const Selector& selector = it->selector();
    if (it->entity_rank() == erank && selector.apply(part_range, comp)) {
      if (restriction == &empty) {
        restriction = &*it;
      }
      if (it->not_equal_stride(*restriction)) {
        ThrowErrorMsg("find_restriction calculation failed with different field-restriction selectors giving incompatible sizes.");
      }
    }
  }

  return *restriction ;
}

} // namespace mesh
} // namespace stk

#endif // stk_mesh_base_FindRestriction_hpp
