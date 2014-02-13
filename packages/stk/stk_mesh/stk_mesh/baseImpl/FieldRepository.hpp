/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
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
        const char                        * arg_method ,
        const std::string                 & arg_name ,
        const DataTraits                  & arg_traits ,
        unsigned                            arg_rank ,
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

    // Fields assocated with each topologyical rank  
    FieldVector m_rankedFields[stk::topology::NUM_RANKS];

    //disallow copy and assignment
    FieldRepository( const FieldRepository &);
    FieldRepository & operator = ( const FieldRepository &);
};

} // namespace impl
} // namespace mesh
} // namespace stk

#endif //stk_mesh_baseImpl_FieldRepository_hpp
