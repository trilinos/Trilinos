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

#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>
#include <algorithm>                    // for lower_bound, sort, unique
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/FieldBase.hpp>  // for FieldBase, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Part.hpp>       // for Part
#include <stk_mesh/base/Types.hpp>      // for ::MaximumFieldDimension
#include <stk_util/util/SimpleArrayOps.hpp>  // for Copy
#include <vector>                       // for vector, etc
#include "Shards_Array.hpp"             // for ArrayDimTag
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/FieldRestriction.hpp"
#include "stk_mesh/base/FieldState.hpp"  // for ::MaximumFieldStates, etc
#include "stk_mesh/base/NgpField.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowErrorMsgIf, etc
#include "stk_util/util/SortAndUnique.hpp"

namespace stk {
namespace mesh {
namespace impl {

FieldBaseImpl::FieldBaseImpl(MetaData                   * arg_mesh_meta_data,
                             stk::topology::rank_t        arg_entity_rank,
                             unsigned                     arg_ordinal,
                             const std::string          & arg_name,
                             const DataTraits           & arg_traits,
                             unsigned                     arg_rank,
                             const shards::ArrayDimTag  * const * arg_dim_tags,
                             unsigned                     arg_number_of_states,
                             FieldState                   arg_this_state)
  : m_entity_rank(arg_entity_rank),
    m_name( arg_name ),
    m_attribute(),
    m_data_traits( arg_traits ),
    m_meta_data( arg_mesh_meta_data ),
    m_ordinal( arg_ordinal ),
    m_num_states( arg_number_of_states ),
    m_this_state( arg_this_state ),
    m_field_rank( arg_rank ),
    m_restrictions(),
    m_initial_value(nullptr),
    m_initial_value_num_bytes(0),
    m_ngpField(nullptr),
    m_numSyncsToHost(0),
    m_numSyncsToDevice(0)
{
  FieldBase * const pzero = nullptr ;
  const shards::ArrayDimTag * const dzero = nullptr ;
  Copy<MaximumFieldStates>(    m_field_states , pzero );
  Copy<MaximumFieldDimension>( m_dim_tags ,     dzero );

  for ( unsigned i = 0 ; i < arg_rank ; ++i ) {
    m_dim_tags[i] = arg_dim_tags[i];
  }
}

//----------------------------------------------------------------------
FieldBaseImpl::~FieldBaseImpl()
{
  delete m_ngpField;

  if (state() == StateNone) {
    void*& init_val = m_initial_value;

    delete [] reinterpret_cast<char*>(init_val);
    init_val = nullptr;
  }
}

//----------------------------------------------------------------------
const FieldRestrictionVector & FieldBaseImpl::restrictions() const
{ return m_field_states[0]->m_impl.m_restrictions ; }

FieldRestrictionVector & FieldBaseImpl::restrictions()
{ return m_field_states[0]->m_impl.m_restrictions ; }

//----------------------------------------------------------------------

template<typename PARTVECTOR>
void get_parts_and_all_subsets(const PARTVECTOR& parts, OrdinalVector& parts_and_all_subsets)
{
  parts_and_all_subsets.clear();
  for(const Part* part : parts) {
    stk::util::insert_keep_sorted_and_unique(part->mesh_meta_data_ordinal(), parts_and_all_subsets);

    const PartVector& subsets = part->subsets();
    for(const Part* subset : subsets) {
      stk::util::insert_keep_sorted_and_unique(subset->mesh_meta_data_ordinal(), parts_and_all_subsets);
    }
  }
}

std::pair<bool,bool> check_for_existing_subsets_or_supersets(FieldRestriction& tmp,
                                                             FieldRestrictionVector::iterator i,
                                                             PartVector& selectorI_parts,
                                                             OrdinalVector& selectorI_parts_and_subsets,
                                                             PartVector& arg_selector_parts,
                                                             OrdinalVector& arg_selector_parts_and_subsets,
                                                             bool arg_selector_is_all_unions
                                                            )
{
  bool found_subset = false;
  bool found_superset = false;

  const Selector& selectorI = i->selector();
  selectorI_parts.clear();
  selectorI.get_parts(selectorI_parts);
  get_parts_and_all_subsets(selectorI_parts, selectorI_parts_and_subsets);
  const bool selectorI_is_all_unions = selectorI.is_all_unions();
  
  const bool both_selectors_are_unions = arg_selector_is_all_unions && selectorI_is_all_unions;
  found_superset = both_selectors_are_unions ? is_subset(arg_selector_parts_and_subsets, selectorI_parts_and_subsets) : false;
  
  found_subset = both_selectors_are_unions ? is_subset(selectorI_parts_and_subsets, arg_selector_parts_and_subsets) : false;
  if (found_subset) {
    *i = tmp;
  }

  return std::make_pair(found_superset, found_subset);
}

void FieldBaseImpl::insert_restriction(
  const char     * arg_method ,
  const Selector & arg_selector ,
  const unsigned   arg_num_scalars_per_entity ,
  const unsigned   arg_first_dimension ,
  const void*      arg_init_value )
{
  FieldRestriction tmp( arg_selector );

  tmp.set_num_scalars_per_entity(arg_num_scalars_per_entity);
  if ( m_field_rank ) {
    tmp.set_dimension(arg_first_dimension);
  }
  else { // Scalar field has m_field_rank==0
    tmp.set_dimension(1);
  }

  if (arg_init_value != NULL) {
    //insert_restriction can be called multiple times for the same field, giving
    //the field different lengths on different mesh-parts.
    //We will only store one initial-value array, we need to store the one with
    //maximum length for this field so that it can be used to initialize data
    //for all field-restrictions. For the parts on which the field is shorter,
    //a subset of the initial-value array will be used.
    //
    //We want to end up storing the longest arg_init_value array for this field.
    //
    //Thus, we call set_initial_value only if the current length is longer
    //than what's already been stored.

    //length in bytes is num-scalars X sizeof-scalar:

    size_t num_scalars = arg_num_scalars_per_entity;
    size_t sizeof_scalar = m_data_traits.size_of;
    size_t nbytes = sizeof_scalar * num_scalars;

    size_t old_nbytes = 0;
    if (get_initial_value() != NULL) {
      old_nbytes = get_initial_value_num_bytes();
    }
    if (nbytes > old_nbytes) {
      set_initial_value(arg_init_value, num_scalars, nbytes);
    }
  }

  {
    FieldRestrictionVector & restrs = restrictions();

    FieldRestrictionVector::iterator restr = restrs.begin();
    FieldRestrictionVector::iterator last_restriction = restrs.end();

    restr = std::lower_bound(restr,last_restriction,tmp);

    const bool new_restriction = ( ( restr == last_restriction ) || !(*restr == tmp) );

    if ( new_restriction ) {
      if (m_meta_data->is_commit()) {
        ThrowRequireMsg(m_meta_data->are_late_fields_enabled(),
                        "Attempting to register Field '" << m_name << "' after MetaData is" << std::endl <<
                        "committed. If you are willing to accept the performance implications, call" << std::endl <<
                        "MetaData::enable_late_fields() before adding these Fields.");
      }

      PartVector arg_selector_parts, selectorI_parts;
      OrdinalVector arg_selector_parts_and_subsets, selectorI_parts_and_subsets;
      arg_selector.get_parts(arg_selector_parts);
      get_parts_and_all_subsets(arg_selector_parts, arg_selector_parts_and_subsets);

      const bool arg_selector_is_all_unions = arg_selector.is_all_unions();
      bool found_superset = false;
      bool found_subset = false;

      for(FieldRestrictionVector::iterator i=restrs.begin(), iend=restrs.end(); i!=iend; ++i) {

        const unsigned i_num_scalars_per_entity = i->num_scalars_per_entity();

        bool shouldCheckForExistingSubsetsOrSupersets =
          arg_num_scalars_per_entity != i_num_scalars_per_entity;
#ifndef NDEBUG
        shouldCheckForExistingSubsetsOrSupersets = true;
#endif
        if (shouldCheckForExistingSubsetsOrSupersets) {
          std::pair<bool,bool> result =
             check_for_existing_subsets_or_supersets(tmp, i,
                                                     selectorI_parts, selectorI_parts_and_subsets,
                                                     arg_selector_parts, arg_selector_parts_and_subsets,
                                                     arg_selector_is_all_unions
                                                      );
          if (result.first) {
            found_superset = true;
          }
          if (result.second) {
            found_subset = true;
          }
        }

        if (found_superset) {
          ThrowErrorMsgIf( i_num_scalars_per_entity != arg_num_scalars_per_entity,
                       " FAILED to add new field-restriction " <<
                       print_restriction( *i, arg_selector, m_field_rank ) <<
                       " WITH INCOMPATIBLE REDECLARATION " <<
                       print_restriction( tmp, arg_selector, m_field_rank ));
          return;
        }
        if (found_subset) {
          ThrowErrorMsgIf( i_num_scalars_per_entity != arg_num_scalars_per_entity,
                           arg_method << " FAILED for " << *this << " " <<
                           print_restriction( *i, arg_selector, m_field_rank ) <<
                           " WITH INCOMPATIBLE REDECLARATION " <<
                           print_restriction( tmp, arg_selector, m_field_rank ));
        }
      }
      if (!found_subset) {
        restrs.insert( restr , tmp );
      }
      else {
        //if subsets were found, we replaced them with the new restriction. so now we need
        //to sort and unique the vector, and trim it to remove any duplicates:
        stk::util::sort_and_unique(restrs);
      }
    }
    else {
      ThrowErrorMsgIf( restr->num_scalars_per_entity() != tmp.num_scalars_per_entity(),
                       arg_method << " FAILED for " << *this << " " <<
                       print_restriction( *restr, arg_selector, m_field_rank ) <<
                       " WITH INCOMPATIBLE REDECLARATION " <<
                       print_restriction( tmp, arg_selector, m_field_rank ));
    }
  }
}

void FieldBaseImpl::verify_and_clean_restrictions(const Part& superset, const Part& subset)
{
  FieldRestrictionVector & restrs = restrictions();

  //Check whether restriction contains subset part, if so, it may now be redundant
  //with another restriction.
  //If they are, make sure they are compatible and remove the subset restrictions.
  PartVector scratch1, scratch2;
  for (size_t r = 0; r < restrs.size(); ++r) {
    FieldRestriction const& curr_restriction = restrs[r];

    if (curr_restriction.selector()(subset)) {
      bool delete_me = false;
      for (size_t i = 0, ie = restrs.size(); i < ie; ++i) {
        FieldRestriction const& check_restriction = restrs[i];
        if (i != r &&
            check_restriction.num_scalars_per_entity() != curr_restriction.num_scalars_per_entity() &&
            check_restriction.selector()(subset) &&
            is_subset(curr_restriction.selector(), check_restriction.selector(), scratch1, scratch2)) {
          ThrowErrorMsgIf( check_restriction.num_scalars_per_entity() != curr_restriction.num_scalars_per_entity(),
                           "Incompatible field restrictions for parts "<< superset.name() << " and "<< subset.name());
          delete_me = true;
          break;
        }
      }
      if (delete_me) {
        restrs.erase(restrs.begin() + r);
        --r;
      }
    }
  }
}

const void* FieldBaseImpl::get_initial_value() const
{
  return m_field_states[0]->m_impl.m_initial_value;
}

void* FieldBaseImpl::get_initial_value() {
  return m_field_states[0]->m_impl.m_initial_value;
}

unsigned FieldBaseImpl::get_initial_value_num_bytes() const {
  return m_field_states[0]->m_impl.m_initial_value_num_bytes;
}

void FieldBaseImpl::set_initial_value(const void* new_initial_value, unsigned num_scalars, unsigned num_bytes) {
  void*& init_val = m_field_states[0]->m_impl.m_initial_value;

  delete [] reinterpret_cast<char*>(init_val);
  init_val = new char[num_bytes];

  m_field_states[0]->m_impl.m_initial_value_num_bytes = num_bytes;

  m_data_traits.copy(init_val, new_initial_value, num_scalars);
}

//----------------------------------------------------------------------

unsigned FieldBaseImpl::max_size( unsigned ent_rank ) const
{
  unsigned max = 0 ;

  if(static_cast<unsigned>(entity_rank()) == ent_rank)
  {
      const FieldRestrictionVector & rMap = restrictions();
      const FieldRestrictionVector::const_iterator ie = rMap.end() ;
            FieldRestrictionVector::const_iterator i = rMap.begin();

      for ( ; i != ie ; ++i ) {
          const unsigned len = i->num_scalars_per_entity();
          if ( max < len ) { max = len ; }
      }
  }
  return max ;
}

void
FieldBaseImpl::modify_on_host() const
{
  if (m_ngpField != nullptr) {
    m_ngpField->modify_on_host();
  }
}

void
FieldBaseImpl::modify_on_device() const
{
  if (m_ngpField != nullptr) {
    m_ngpField->modify_on_device();
  }
}

void
FieldBaseImpl::sync_to_host() const
{
  if (m_ngpField != nullptr) {
    m_ngpField->sync_to_host();
  }
}

void
FieldBaseImpl::sync_to_device() const
{
  if (m_ngpField != nullptr) {
    m_ngpField->sync_to_device();
  }
}

void
FieldBaseImpl::clear_sync_state() const
{
  if (m_ngpField != nullptr) {
    m_ngpField->clear_sync_state();
  }
}

NgpFieldBase *
FieldBaseImpl::get_ngp_field() const
{
  return m_ngpField;
}

void
FieldBaseImpl::set_ngp_field(NgpFieldBase * ngpField) const
{
  ThrowRequireMsg(m_ngpField == nullptr || m_ngpField == ngpField,
                  "Error: Only one NgpField may be set on a StkField(" << m_name << ")");
  m_ngpField = ngpField;
}

size_t
FieldBaseImpl::num_syncs_to_host() const
{
  return m_numSyncsToHost;
}

size_t
FieldBaseImpl::num_syncs_to_device() const
{
  return m_numSyncsToDevice;
}

void
FieldBaseImpl::increment_num_syncs_to_host() const
{
  ++m_numSyncsToHost;
}

void
FieldBaseImpl::increment_num_syncs_to_device() const
{
  ++m_numSyncsToDevice;
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const FieldBaseImpl & field )
{
  s << "FieldBaseImpl<" ;
  s << field.data_traits().name ;
  for ( unsigned i = 0 ; i < field.field_array_rank() ; ++i ) {
    s << "," << field.dimension_tags()[i]->name();
  }
  s << ">" ;

  s << "[ name = \"" ;
  s << field.name() ;
  s << "\" , #states = " ;
  s << field.number_of_states();
  s << " ]" ;
  return s ;
}

std::ostream & print( std::ostream & s ,
                      const char * const b ,
                      const FieldBase & field )
{
  const std::vector<FieldBase::Restriction> & rMap = field.restrictions();
  s << field.name() ;
  s << " {" ;
  for ( FieldBase::RestrictionVector::const_iterator
        i = rMap.begin() ; i != rMap.end() ; ++i ) {
    s << std::endl << b << "  " ;
    i->print( s, i->selector(), field.field_array_rank() );
    s << std::endl;
  }
  s << std::endl << b << "}" ;
  return s ;
}

} // namespace impl
} // namespace mesh
} // namespace stk
