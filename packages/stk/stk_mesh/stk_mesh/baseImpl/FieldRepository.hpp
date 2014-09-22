// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef stk_mesh_baseImpl_FieldRepository_hpp
#define stk_mesh_baseImpl_FieldRepository_hpp

#include <stddef.h>                     // for NULL
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>  // for FieldBaseImpl
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_topology/topology.hpp"    // for topology, topology::rank_t, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowAssert
namespace shards { class ArrayDimTag; }
namespace stk { namespace mesh { class DataTraits; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }


namespace stk {
namespace mesh {

typedef std::vector<FieldBase *> FieldVector;

namespace impl {

class FieldRepository {

  public:
    FieldRepository() {}
    ~FieldRepository();

    void add_field(FieldBase* new_field) {
      m_fields.push_back(new_field);
      m_rankedFields[new_field->entity_rank()].push_back(new_field);
    }

    FieldBase * get_field(
        stk::topology::rank_t               arg_entity_rank ,
        const std::string                 & arg_name ,
        const DataTraits                  & arg_traits ,
        unsigned                            arg_array_rank ,
        const shards::ArrayDimTag * const * arg_dim_tags ,
        unsigned                            arg_num_states
        ) const;

    FieldBase * declare_field(
        const std::string                 & arg_name ,
        stk::topology::rank_t               arg_entity_rank ,
        const DataTraits                  & arg_traits ,
        unsigned                            arg_rank ,
        const shards::ArrayDimTag * const * arg_dim_tags ,
        unsigned                            arg_num_states ,
        MetaData                          * arg_meta_data
        );

    void verify_and_clean_restrictions(const Part& superset, const Part& subset);

    const FieldVector & get_fields() const {
      return m_fields;
    }

    // return all fields of a given topological rank (node/face/elem, etc.)
    const FieldVector & get_fields(stk::topology::rank_t rank) const {
      ThrowAssert(rank <= stk::topology::NUM_RANKS);
      return m_rankedFields[rank];
    }


    template<class T>
      const T *
      declare_attribute_with_delete( FieldBase & f , const T * a )
      {
        return f.m_impl.declare_attribute_with_delete(a);
      }

    template<class T>
      const T *
      declare_attribute_no_delete( FieldBase & f , const T * a )
      {
        return f.m_impl.declare_attribute_no_delete(a);
      }

    template<class T>
      bool remove_attribute( FieldBase & f , const T * a )
      {
	return f.m_impl.remove_attribute<T>( a );
      }

    void declare_field_restriction(
        const char     * arg_method,
        FieldBase      & arg_field ,
        const Part     & arg_part ,
        const PartVector & arg_all_parts,
        const unsigned   arg_num_scalars_per_entity ,
        const unsigned   arg_first_dimension ,
        const void     * arg_init_value = NULL)
    {
      arg_field.m_impl.insert_restriction( arg_method, arg_part, arg_num_scalars_per_entity, arg_first_dimension, arg_init_value);
    }

    void declare_field_restriction(
        const char     * arg_method,
        FieldBase      & arg_field ,
        const Selector & arg_selector ,
        const PartVector & arg_all_parts,
        const unsigned   arg_num_scalars_per_entity ,
        const unsigned   arg_first_dimension ,
        const void     * arg_init_value = NULL)
    {
      arg_field.m_impl.insert_restriction( arg_method, arg_selector, arg_num_scalars_per_entity, arg_first_dimension, arg_init_value);
    }

  private:
    FieldVector m_fields;

    // Fields assocated with each topological rank  
    FieldVector m_rankedFields[stk::topology::NUM_RANKS];

    //disallow copy and assignment
    FieldRepository( const FieldRepository &);
    FieldRepository & operator = ( const FieldRepository &);
};

} // namespace impl
} // namespace mesh
} // namespace stk

#endif //stk_mesh_baseImpl_FieldRepository_hpp
