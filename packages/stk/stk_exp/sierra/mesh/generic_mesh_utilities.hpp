#ifndef STK_SIERRA_MESH_GENERIC_MESH_UTILITIES_HPP
#define STK_SIERRA_MESH_GENERIC_MESH_UTILITIES_HPP

namespace sierra {

namespace private_details {

template<class Iter, class Mesh, class Functor>
inline
bool is_selected(Iter i, Iter j,
                 const typename mesh::mesh_traits<Mesh>::bucket_key& bucket,
                 Functor has_part,
                 const Mesh& m)
{
  bool result = i != j;
  while( result && i != j ) {
    if ( i->m_count ) { // Compound statement
      result = i->m_unary ^ is_selected(i+1, i+i->m_count, bucket, has_part, m);
      i += i->m_count;
    }
    else { // Test for containment of bucket in this part
      result = i->m_unary ^ has_part(bucket, i->m_part_id, m);
      ++i;
    }
  }
  return result;
}

}//private_details

namespace mesh {

template<class Mesh>
inline
bool bucket_has_part(const typename mesh_traits<Mesh>::bucket_key& bucket,
                     const typename mesh_traits<Mesh>::part_key& part_id,
                     const Mesh& m)
{
  typedef typename mesh_traits<Mesh>::bucket_part_range part_range;
  typedef typename mesh_traits<Mesh>::bucket_part_iterator part_iterator;
  typedef typename mesh_traits<Mesh>::part_key part_key;

  part_range parts = get_parts(bucket, m);
  for(part_iterator p_it=parts.first, p_end=parts.second; p_it!=p_end; ++p_it) {
    if (part_id == *p_it) {
      return true;
    }
  }
  return false;
}

  template<class MESH>
  struct has_part {
    bool operator()(const typename mesh_traits<MESH>::bucket_key& bucket,
                     const typename mesh_traits<MESH>::part_key& part_id,
                     const MESH& m) const
    {
      return bucket_has_part(bucket, part_id, m);
    }
  };

template<class Selector, class Mesh>
inline
bool is_selected(const typename mesh_traits<Mesh>::bucket_key& bucket,
                 const Selector& select,
                 const Mesh& m)
{
  const std::vector<typename Selector::OpType>& ops = select.get_ops();
  typename std::vector<typename Selector::OpType>::const_iterator
      i = ops.begin(), j = ops.end();

  has_part<Mesh> hp;
  return sierra::private_details::is_selected(i, j, bucket, hp, m);
}

}//namespace mesh
}//namespace sierra

#endif

