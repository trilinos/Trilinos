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

#ifndef stk_mesh_baseImpl_FieldBaseImpl_hpp
#define stk_mesh_baseImpl_FieldBaseImpl_hpp

#include <stddef.h>                     // for NULL
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/DataTraits.hpp>  // for DataTraits
#include <stk_mesh/base/FieldRestriction.hpp>
#include <stk_mesh/base/FieldState.hpp>  // for FieldState, etc
#include <stk_util/util/CSet.hpp>       // for CSet
#include <string>                       // for string
#include <typeinfo>                     // for type_info
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for ::MaximumFieldDimension
#include "stk_topology/topology.hpp"    // for topology, topology::rank_t
namespace shards { class ArrayDimTag; }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class NgpFieldBase; } }

namespace stk {
namespace mesh {
namespace impl {

class FieldBaseImpl {
public:
  FieldBaseImpl(MetaData                   * arg_mesh_meta_data,
                stk::topology::rank_t        arg_entity_rank,
                unsigned                     arg_ordinal,
                const std::string          & arg_name,
                const DataTraits           & arg_traits,
                unsigned                     arg_rank,
                const shards::ArrayDimTag  * const * arg_dim_tags,
                unsigned                     arg_number_of_states,
                FieldState                   arg_this_state);

  ~FieldBaseImpl();

  FieldBaseImpl() = delete;
  FieldBaseImpl(const FieldBase &) = delete;
  FieldBaseImpl(const FieldBaseImpl &) = delete;
  FieldBaseImpl & operator=(const FieldBaseImpl &) = delete;

  MetaData & meta_data() const {
    return *m_meta_data ;
  }

  unsigned ordinal() const {
    return m_ordinal ;
  }

  const std::string & name() const {
    return m_name ;
  }

  template<class Type>
  bool type_is() const {
    return m_data_traits.type_info == typeid(Type) ;
  }

  const DataTraits & data_traits() const {
    return m_data_traits ;
  }

  unsigned number_of_states() const {
    return m_num_states ;
  }

  FieldState state() const {
    return m_this_state ;
  }

  unsigned field_array_rank() const {
    return m_field_rank ;
  }

  EntityRank entity_rank() const {
      return m_entity_rank;
  }

  const shards::ArrayDimTag * const * dimension_tags() const {
    return m_dim_tags ;
  }

  //not inline
  unsigned max_size( unsigned ent_rank) const ;

  template<class A>
  const A * attribute() const {
    return m_attribute.template get<A>();
  }

  template<class A>
    const A * declare_attribute_no_delete(const A * a) {
      return m_attribute.template insert_no_delete<A>(a);
    }

  template<class A>
    const A * declare_attribute_with_delete(const A * a) {
      return m_attribute.template insert_with_delete<A>(a);
    }

  template<class A>
    bool remove_attribute(const A * a) {
      return m_attribute.template remove<A>(a);
    }

  const FieldRestrictionVector &restrictions() const ;

  FieldBase * field_state(FieldState fstate) const {
    return m_field_states[fstate];
  }

  void insert_restriction( const char       * arg_method ,
                           const Selector   & arg_selector ,
                           const unsigned     arg_num_scalars_per_entity ,
                           const unsigned     arg_first_dimension ,
                           const void*        arg_init_value = NULL);

  void verify_and_clean_restrictions( const Part& superset, const Part& subset );

  const void* get_initial_value() const;

  void* get_initial_value();

  unsigned get_initial_value_num_bytes() const;

  template<typename FieldType>
  void set_field_states( FieldType ** field_states)
  {
    for (unsigned i = 0; i < m_num_states; ++i) {
      m_field_states[i] = field_states[i];
    }
  }

  void modify_on_host() const;
  void modify_on_device() const;
  void sync_to_host() const;
  void sync_to_device() const;
  void clear_sync_state() const;

  NgpFieldBase * get_ngp_field() const;
  void set_ngp_field(NgpFieldBase * ngpField) const;

  size_t num_syncs_to_host() const;
  size_t num_syncs_to_device() const;

  void increment_num_syncs_to_host() const;
  void increment_num_syncs_to_device() const;

private:

  void set_initial_value(const void* new_initial_value, unsigned num_scalars, unsigned num_bytes);

  FieldRestrictionVector & restrictions();

  EntityRank                   m_entity_rank ;             ///< entity-rank that this field can be allocated for
  const std::string            m_name ;                    ///< Name of the field
  CSet                         m_attribute ;               ///< User's attributes
  const DataTraits           & m_data_traits ;             ///< Data type traits
  MetaData             * const m_meta_data ;               ///< Owner of this field
  const unsigned               m_ordinal ;                 ///< Ordinal in the field set
  const unsigned               m_num_states ;              ///< Number of states
  const FieldState             m_this_state ;              ///< Field state of this field
  unsigned                     m_field_rank ;              ///< Number of dimensions
  FieldRestrictionVector       m_restrictions ;            ///< Only valid on StateNone
  FieldBase                  * m_field_states[ MaximumFieldStates ];
  const shards::ArrayDimTag  * m_dim_tags[ MaximumFieldDimension ];
  void*                        m_initial_value;
  unsigned                     m_initial_value_num_bytes;
  mutable NgpFieldBase       * m_ngpField;
  mutable size_t               m_numSyncsToHost;
  mutable size_t               m_numSyncsToDevice;
};


/** \brief  Print the field type, text name, and number of states. */
std::ostream & operator << ( std::ostream & , const FieldBaseImpl & );

/** \brief  Print field and field restrictions on new lines. */
std::ostream & print( std::ostream & ,
                      const char * const , const FieldBase & );


} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_FieldBaseIml_hpp
