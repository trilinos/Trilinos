/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_PartRepository_hpp
#define stk_mesh_PartRepository_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/baseImpl/PartImpl.hpp>


namespace stk_classic {
namespace mesh {

class MetaData;

namespace impl {


class PartRepository {
public:
  explicit PartRepository(MetaData * meta);
  ~PartRepository();

  Part * universal_part() const;

  const PartVector & get_all_parts() const;

  Part * declare_part( const std::string & arg_name , EntityRank arg_rank );
  Part * declare_part( const PartVector & part_intersect );
  void declare_subset( Part & superset, Part & subset );
  void declare_part_relation( Part & root_part, PartRelation relation, Part & target_part );

  template<class T>
  const T * declare_attribute_with_delete( Part & , const T *);
  template<class T>
  const T * declare_attribute_no_delete( Part & , const T *);
  template<class T>
  bool remove_attribute( Part & , const T *);

private:
  PartRepository();
  PartRepository(const PartRepository & );
  PartRepository & operator = ( const PartRepository & );
  
  Part * declare_part_impl( const std::string & name, EntityRank rank);
  void declare_subset_impl( Part & superset, Part & subset );

  MetaData * m_meta_data;
  Part * m_universal_part;
  PartVector m_all_parts;
};

template<class T>
inline
const T *
PartRepository::declare_attribute_with_delete( Part & p, const T * a )
{
  return p.m_partImpl.declare_attribute_with_delete<T>( a );
}

template<class T>
inline
const T *
PartRepository::declare_attribute_no_delete( Part & p, const T * a )
{
  return p.m_partImpl.declare_attribute_no_delete<T>( a );
}

template<class T>
inline
bool
PartRepository::remove_attribute( Part & p, const T * a )
{
  return p.m_partImpl.remove_attribute<T>( a );
}


} // namespace impl 
} // namespace mesh
} // namespace stk_classic


#endif // stk_mesh_PartRepository_hpp
