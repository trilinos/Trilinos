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

#include <stk_mesh/base/FieldBase.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <vector>                       // for vector, etc
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/MetaData.hpp"  // for FieldRestriction
#include "stk_mesh/base/FieldRestriction.hpp"  // for FieldRestriction
#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_mesh/base/NgpFieldBase.hpp>
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg


namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

template<typename PARTVECTOR>
void get_parts_and_all_subsets(const PARTVECTOR& parts, OrdinalVector& parts_and_all_subsets)
{
  parts_and_all_subsets.clear();
  for(const Part* part : parts) {
    parts_and_all_subsets.push_back(part->mesh_meta_data_ordinal());
  }
  stk::util::sort_and_unique(parts_and_all_subsets);

  OrdinalVector scratch;
  scratch.reserve(256);
  for(const Part* part : parts) {
    const PartVector& subsets = part->subsets();
    for(const Part* subset : subsets) {
      unsigned ord = subset->mesh_meta_data_ordinal();
      scratch.push_back(ord);
    }

    if (scratch.size() >= parts_and_all_subsets.size()) {
      stk::util::sort_and_unique(scratch);
      stk::util::insert_keep_sorted_and_unique(scratch, parts_and_all_subsets, std::less<unsigned>());
      scratch.clear();
    }
  }

  if (!scratch.empty()) {
    stk::util::sort_and_unique(scratch);
    stk::util::insert_keep_sorted_and_unique(scratch, parts_and_all_subsets, std::less<unsigned>());
  }
}

template<typename PARTVECTOR>
std::pair<bool,bool> check_for_existing_subsets_or_supersets(FieldRestriction& tmp,
                                                             FieldRestrictionVector::iterator i,
                                                             PARTVECTOR& selectorI_parts,
                                                             OrdinalVector& selectorI_parts_and_subsets,
                                                             PARTVECTOR& arg_selector_parts,
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

std::ostream & operator<<(std::ostream & s, const FieldBase & field)
{
  s << "Field<" << field.data_traits().name << ">[\"" << field.name() << "\", #states: "
    << field.number_of_states() << "]";
  return s ;
}

std::ostream & print(std::ostream & s, const char * const b, const FieldBase & field)
{
  s << b << field << std::endl;
  std::string indent = b;
  indent += "  ";
  print_restrictions(s, indent.c_str(), field);
  return s ;
}

std::ostream & print_restrictions(std::ostream & s, const char * const b, const FieldBase & field)
{
  const std::vector<FieldBase::Restriction> & rMap = field.restrictions();

  for ( const FieldBase::Restriction& r : rMap ) {
    s << b;
    r.print(s, r.selector());
    s << std::endl;
  }
  return s;
}

void FieldBase::set_initial_value(const void* new_initial_value, unsigned num_scalars, unsigned num_bytes) {
  void*& init_val = field_state(StateNone)->m_initial_value;

  delete [] reinterpret_cast<char*>(init_val);
  init_val = new char[num_bytes];

  m_field_states[0]->m_initial_value_num_bytes = num_bytes;

  data_traits().copy(init_val, new_initial_value, num_scalars);
}

void FieldBase::insert_restriction(const char     * arg_method,
                                   const Part     & arg_part,
                                   const unsigned   arg_num_scalars_per_entity,
                                   const unsigned   arg_first_dimension,
                                   const void*      arg_init_value)
{
  FieldRestriction tmp( arg_part );

  tmp.set_num_scalars_per_entity(arg_num_scalars_per_entity);
  tmp.set_dimension(arg_first_dimension);

  if (arg_init_value != nullptr) {
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
    size_t sizeof_scalar = data_traits().size_of;
    size_t nbytes = sizeof_scalar * num_scalars;

    size_t old_nbytes = 0;
    if (get_initial_value() != nullptr) {
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
      ConstPartVector arg_parts, selectorI_parts;
      OrdinalVector arg_parts_and_subsets, selectorI_parts_and_subsets;
      arg_parts = {&arg_part};
      get_parts_and_all_subsets(arg_parts, arg_parts_and_subsets);

      const bool arg_is_all_unions = true;
      bool found_superset = false;
      bool found_subset = false;


      for(FieldRestrictionVector::iterator i=restrs.begin(), iend=restrs.end(); i!=iend; ++i) {
        const unsigned i_num_scalars_per_entity = i->num_scalars_per_entity();

        if (i->selects(arg_part)) {
          STK_ThrowRequireMsg(i_num_scalars_per_entity == arg_num_scalars_per_entity && i->dimension() == tmp.dimension(),
            "Error, adding field-restriction to field "<<name()
            <<", existing field-restriction with selector "<<i->selector()
            <<" selects incoming part "<<arg_part.name()
            <<" but scalars-per-entity is "<<i_num_scalars_per_entity
            <<" inconsistent with incoming scalars-per-entity "<<arg_num_scalars_per_entity
            <<" or dimension "<<i->dimension()<<" inconsistent with incoming dimension "<<tmp.dimension());
          return;
        }

        if (i_num_scalars_per_entity == arg_num_scalars_per_entity && i->dimension() == tmp.dimension()) {
          if (mesh_meta_data().is_commit()) {
            STK_ThrowRequireMsg(mesh_meta_data().are_late_fields_enabled(),
                          "Attempting to register Field '" << m_name << "' after MetaData is" << std::endl <<
                          "committed. If you are willing to accept the performance implications, call" << std::endl <<
                          "MetaData::enable_late_fields() before adding these Fields.");
          }

          std::pair<bool,bool> result =
            check_for_existing_subsets_or_supersets(tmp, i,
                                                    selectorI_parts, selectorI_parts_and_subsets,
                                                    arg_parts, arg_parts_and_subsets,
                                                    arg_is_all_unions);
          found_superset = result.first;
          found_subset = result.second;

          //if found_superset, nothing to do because the incoming restriction is already part of the
          //existing restriction.
          //if found_subsset, the 'check_for_existing..' function already replaced the existing
          //restriction with the incoming restriction. (Confusing. Fix the check function to only check.)
          if (!found_superset && !found_subset) {
            i->add_union(arg_part);
          }
          return;
        }
      }

      if (mesh_meta_data().is_commit()) {
        STK_ThrowRequireMsg(mesh_meta_data().are_late_fields_enabled(),
                      "Attempting to register Field '" << m_name << "' after MetaData is" << std::endl <<
                      "committed. If you are willing to accept the performance implications, call" << std::endl <<
                      "MetaData::enable_late_fields() before adding these Fields.");
      }
      restrs.push_back(tmp);
      stk::util::sort_and_unique(restrs);
    }
    else {
      STK_ThrowErrorMsgIf(restr->num_scalars_per_entity() != tmp.num_scalars_per_entity(),
                      arg_method << " FAILED for " << *this << " " << print_restriction(*restr, restr->selector()) <<
                      " WITH INCOMPATIBLE REDECLARATION " << print_restriction(tmp, Selector(arg_part)));
    }
  }
}

void FieldBase::insert_restriction(const char     * arg_method,
                                   const Selector & arg_selector,
                                   const unsigned   arg_num_scalars_per_entity,
                                   const unsigned   arg_first_dimension,
                                   const void*      arg_init_value)
{
  FieldRestriction tmp( arg_selector );

  tmp.set_num_scalars_per_entity(arg_num_scalars_per_entity);
  tmp.set_dimension(arg_first_dimension);

  if (arg_init_value != nullptr) {
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
    size_t sizeof_scalar = data_traits().size_of;
    size_t nbytes = sizeof_scalar * num_scalars;

    size_t old_nbytes = 0;
    if (get_initial_value() != nullptr) {
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
      PartVector arg_selector_parts, selectorI_parts;
      OrdinalVector arg_selector_parts_and_subsets, selectorI_parts_and_subsets;
      arg_selector.get_parts(arg_selector_parts);
      get_parts_and_all_subsets(arg_selector_parts, arg_selector_parts_and_subsets);

      const bool arg_selector_is_all_unions = arg_selector.is_all_unions();
      bool found_superset = false;
      bool found_subset = false;

      for(FieldRestrictionVector::iterator i=restrs.begin(), iend=restrs.end(); i!=iend; ++i) {

        const unsigned i_num_scalars_per_entity = i->num_scalars_per_entity();

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

        if (found_superset || found_subset) {
          STK_ThrowErrorMsgIf(i_num_scalars_per_entity != arg_num_scalars_per_entity,
                          "FAILED to add new field-restriction " << print_restriction(*i, arg_selector) <<
                          " WITH INCOMPATIBLE REDECLARATION " << print_restriction(tmp, arg_selector));
        }
        if (found_superset) {
          return;
        }
      }

      if (!found_subset) {
        if (mesh_meta_data().is_commit()) {
          if (!mesh_meta_data().are_late_fields_enabled()) {
            std::cerr<<"field "<<m_name<<", arg sel: "<<arg_selector<<", old restrs: "<<std::endl;
            for(const FieldRestriction& r : restrs) {
              std::cerr<<" -  "<<r.selector()<<std::endl;
            }
          }
          STK_ThrowRequireMsg(mesh_meta_data().are_late_fields_enabled(),
                          "Attempting to register Field '" << m_name << "' after MetaData is" << std::endl <<
                          "committed. If you are willing to accept the performance implications, call" << std::endl <<
                          "MetaData::enable_late_fields() before adding these Fields.");
        }   

        bool addedToUnion = false;
        for(FieldRestriction& r : restrs) {
          if (r.num_scalars_per_entity() == tmp.num_scalars_per_entity() && r.dimension() == tmp.dimension()) {
            r.add_union(tmp.selector());
            addedToUnion = true;
            break;
          }
        }
        if (!addedToUnion) {
          restrs.insert( restr , tmp );
        }
      }
      else {
        //if subsets were found, we replaced them with the new restriction. so now we need
        //to sort and unique the vector, and trim it to remove any duplicates:
        stk::util::sort_and_unique(restrs);
      }
    }
    else {
      STK_ThrowErrorMsgIf(restr->num_scalars_per_entity() != tmp.num_scalars_per_entity(),
                      arg_method << " FAILED for " << *this << " " << print_restriction(*restr, restr->selector()) <<
                      " WITH INCOMPATIBLE REDECLARATION " << print_restriction(tmp, tmp.selector()));
    }
  }
}

void FieldBase::verify_and_clean_restrictions(const Part& superset, const Part& subset)
{
  FieldRestrictionVector & restrs = restrictions();
    
  //Check whether restriction contains subset part, if so, it may now be redundant
  //with another restriction.
  //If they are, make sure they are compatible and remove the subset restrictions.
  PartVector scratch1, scratch2;
  const FieldRestriction* restrData = restrs.data();
  std::vector<unsigned> scratch;
  for (size_t r = 0; r < restrs.size(); ++r) {
    FieldRestriction const& curr_restriction = restrData[r];
      
    if (curr_restriction.selector()(subset)) {
      scratch.push_back(r);
    }
  } 
  
  for (size_t r = 0; r < scratch.size(); ++r) {
    FieldRestriction const& curr_restriction = restrData[scratch[r]];
  
    bool delete_me = false;
    for (size_t i = 0, ie = scratch.size(); i < ie; ++i) {
      FieldRestriction const& check_restriction = restrData[scratch[i]];
      if (i != r &&
          check_restriction.num_scalars_per_entity() != curr_restriction.num_scalars_per_entity() &&
          is_subset(curr_restriction.selector(), check_restriction.selector(), scratch1, scratch2)) {
        STK_ThrowErrorMsgIf( check_restriction.num_scalars_per_entity() != curr_restriction.num_scalars_per_entity(),
                         "Incompatible field restrictions for parts "<< superset.name() << " and "<< subset.name());
        delete_me = true;
        break;
      }
    }
    if (delete_me) {
      restrs.erase(restrs.begin() + scratch[r]);
      for(unsigned j=r+1; j<scratch.size(); ++j) {
        --scratch[j];
      }
    }
  }
}

void FieldBase::set_mesh(stk::mesh::BulkData* bulk)
{
  if (m_mesh == nullptr || bulk == nullptr) {
    m_mesh = bulk;
  }
  else {
    STK_ThrowRequireMsg(bulk == m_mesh, "Internal Error: Trying to use field " << name() << " on more than one bulk data");
  }
}

bool FieldBase::defined_on(const stk::mesh::Part& part) const
{
  return (length(part) > 0);
}

unsigned FieldBase::length(const stk::mesh::Part& part) const
{
  const stk::mesh::FieldRestriction& restriction = stk::mesh::find_restriction(*this, entity_rank(), part);
  return restriction.num_scalars_per_entity();
}

unsigned FieldBase::max_size() const
{
  FieldRestriction::size_type max = 0 ; 

  for (const FieldRestriction& restriction : restrictions()) {
    max = std::max(max, restriction.num_scalars_per_entity());
  }   

  return max ;
}

unsigned FieldBase::max_extent(unsigned dimension) const
{
  if (dimension == 0) {
    FieldRestriction::size_type max = 0;
    for (const FieldRestriction& restriction : restrictions()) {
      max = std::max(max, restriction.dimension());
    }
    return max;
  }
  else if (dimension == 1) {
    FieldRestriction::size_type max = 0;
    for (const FieldRestriction& restriction : restrictions()) {
      if (restriction.dimension() != 0) {
        max = std::max(max, restriction.num_scalars_per_entity() / restriction.dimension());
      }
    }
    return max;

  }
  else {
    for (const FieldRestriction& restriction : restrictions()) {
      if (restriction.num_scalars_per_entity() > 0) {
        return 1;
      }
    }
    return 0;
  }
}

void FieldBase::rotate_multistate_data(bool rotateNgpFieldViews)
{
  const int numStates = number_of_states();
  if (numStates > 1 && StateNew == state()) {
    bool allStatesHaveNgpFields = true;
    for(int s = 0; s < numStates; ++s) {
      if (field_state(static_cast<FieldState>(s))->get_ngp_field() == nullptr) {
        allStatesHaveNgpFields = false;
      }
    }

    for (int s = 1; s < numStates; ++s) {
      FieldBase* sField = field_state(static_cast<FieldState>(s));
      m_field_meta_data.swap(sField->m_field_meta_data);

      std::swap(m_numSyncsToDevice, sField->m_numSyncsToDevice);
      std::swap(m_numSyncsToHost, sField->m_numSyncsToHost);
      std::swap(m_modifiedOnHost, sField->m_modifiedOnHost);
      std::swap(m_modifiedOnDevice, sField->m_modifiedOnDevice);
    }

    for(int s = 0; s < numStates; ++s) {
      NgpFieldBase* ngpField = field_state(static_cast<FieldState>(s))->get_ngp_field();
      if (ngpField != nullptr) {
        ngpField->update_bucket_pointer_view();
        ngpField->fence();
      }
    }

    if (rotateNgpFieldViews && allStatesHaveNgpFields) {
      for (int s = 1; s < numStates; ++s) {
        NgpFieldBase* ngpField_sminus1 = field_state(static_cast<FieldState>(s-1))->get_ngp_field();
        NgpFieldBase* ngpField_s = field_state(static_cast<FieldState>(s))->get_ngp_field();
        ngpField_s->swap_field_views(ngpField_sminus1);
      }
    }
  }
}

void
FieldBase::modify_on_host() const
{ 
  STK_ThrowRequireMsg(m_modifiedOnDevice == false,
                  "Modify on host called for Field: " << name() << " but it has an uncleared modified_on_device");

  m_modifiedOnHost = true;
}

void
FieldBase::modify_on_device() const
{ 
  STK_ThrowRequireMsg(m_modifiedOnHost == false,
                  "Modify on device called for Field: " << name() << " but it has an uncleared modified_on_host");

  m_modifiedOnDevice = true;
}

void
FieldBase::modify_on_host(const Selector& s) const
{ 
  modify_on_host();
}

void
FieldBase::modify_on_device(const Selector& s) const
{ 
  modify_on_device();
}

bool
FieldBase::need_sync_to_device() const
{
  return m_modifiedOnHost;
}

bool
FieldBase::need_sync_to_host() const
{
  return m_modifiedOnDevice;
}

void
FieldBase::sync_to_host() const
{ 
  if (m_ngpField != nullptr) {
    m_ngpField->sync_to_host();
  } else {
    clear_device_sync_state();
  }
}

void FieldBase::sync_to_host(const stk::ngp::ExecSpace& exec_space) const
{
  if (m_ngpField != nullptr) {
    m_ngpField->sync_to_host(exec_space);
  } else {
    clear_device_sync_state();
  }
}

void
FieldBase::sync_to_device() const
{
  if (m_ngpField != nullptr) {
    m_ngpField->sync_to_device();
  } else {
    clear_host_sync_state();
  }
}

void FieldBase::sync_to_device(const stk::ngp::ExecSpace& exec_space) const
{
  if (m_ngpField != nullptr) {
    m_ngpField->sync_to_device(exec_space);
  } else {
    clear_host_sync_state();
  }
}

void
FieldBase::clear_sync_state() const
{
  if(m_ngpField != nullptr) {
    m_ngpField->notify_sync_debugger_clear_sync_state();
  }
  m_modifiedOnHost = false;
  m_modifiedOnDevice = false;
}

void
FieldBase::clear_host_sync_state() const
{
  if(m_ngpField != nullptr) {
    m_ngpField->notify_sync_debugger_clear_host_sync_state();
  }
  m_modifiedOnHost = false;
}

void
FieldBase::clear_device_sync_state() const
{
  if(m_ngpField != nullptr) {
    m_ngpField->notify_sync_debugger_clear_device_sync_state();
  }
  m_modifiedOnDevice = false;
}

void
FieldBase::set_ngp_field(NgpFieldBase * ngpField) const
{
  STK_ThrowRequireMsg(m_ngpField == nullptr || m_ngpField == ngpField,
                  "Error: Only one NgpField may be set on a StkField(" /*<< m_name <<*/ ")");
  m_ngpField = ngpField;
}

void
FieldBase::fence() const
{
  if(m_ngpField != nullptr) {
    m_ngpField->fence();
  }
}

size_t
FieldBase::num_syncs_to_host() const
{
  return m_numSyncsToHost;
}

size_t
FieldBase::num_syncs_to_device() const
{
  return m_numSyncsToDevice;
}

void
FieldBase::increment_num_syncs_to_host() const
{
  ++m_numSyncsToHost;
}

void
FieldBase::increment_num_syncs_to_device() const
{
  ++m_numSyncsToDevice;
}

namespace impl {

stk::CSet & get_attributes(FieldBase & field) {
  return field.get_attributes();
}

} // namespace impl

} // namespace mesh
} // namespace stk

