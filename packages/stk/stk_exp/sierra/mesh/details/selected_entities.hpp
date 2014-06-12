#ifndef SIERRA_SIERRA_MESH_SELECTED_ENTITIES_HPP
#define SIERRA_SIERRA_MESH_SELECTED_ENTITIES_HPP

#include <sierra/mesh/details/selected_buckets.hpp>

#include <boost/range.hpp>

namespace sierra {
namespace mesh {
namespace details {

//selected_entity_iterator is an iterator that uses a selector to decide
//whether the next bucket in a sequence should be used.
//Users shouldn't construct this class directly, they should call
//the get_selected_buckets function below.
template <class Mesh>
class selected_entity_iterator
  : public std::iterator< std::forward_iterator_tag,
                          entity_descriptor,
                          ptrdiff_t,
                          const entity_descriptor *,
                          entity_descriptor
                        >
{
  typedef mesh_traits<Mesh> MeshTraits;
 public:

  selected_entity_iterator( const Mesh& mesh,
                            typename MeshTraits::selected_bucket_range arg_rng
                          )
   : m_mesh(&mesh)
   , m_bucket_range(arg_rng)
   , m_entity_range()
   , m_current_entity()
  {
    if ( !boost::empty(m_bucket_range) ) {
      m_entity_range = get_entities(*m_bucket_range.first, *m_mesh);
      m_current_entity = *m_entity_range.first;
    }
  }

  selected_entity_iterator()
   : m_mesh(NULL)
   , m_bucket_range()
   , m_entity_range()
   , m_current_entity()
  {
  }

  bool operator!=(const selected_entity_iterator& itr) const {
    return m_current_entity != itr.m_current_entity;
  }

  bool operator==(const selected_entity_iterator& itr) const {
    return m_current_entity == itr.m_current_entity;
  }

  selected_entity_iterator operator++() {
    update();
    return *this;
  }

  selected_entity_iterator operator++(int) {
    selected_entity_iterator temp(*this);
    update();
    return temp;
  }

  entity_descriptor operator *() {
    return *m_entity_range.first;
  }

 private:

  void update() {
    ++m_entity_range.first;
    //test whether we're at the end of the current bucket:
    if (m_entity_range.first == m_entity_range.second) {
      ++m_bucket_range.first;
      //test whether we're at the end of all the buckets:
      if (m_bucket_range.first != m_bucket_range.second) {
        m_entity_range = get_entities(*m_bucket_range.first, *m_mesh);
        m_current_entity = *m_entity_range.first;
      }
      else {
        typename mesh_traits<Mesh>::entity_descriptor invalid;
        m_current_entity = invalid;
      }
    }
  }

  const Mesh           * m_mesh;
  typename mesh_traits<Mesh>::selected_bucket_range   m_bucket_range;
  typename mesh_traits<Mesh>::bucket_entity_range   m_entity_range;
  typename mesh_traits<Mesh>::entity_descriptor   m_current_entity;
};

template <class Mesh>
inline
std::pair<selected_entity_iterator<Mesh>, selected_entity_iterator<Mesh> >
get_entities(const selector& select, const Mesh& mesh)
{
  typename mesh_traits<Mesh>::selected_bucket_range selected_buckets =
     get_selected_buckets(select, mesh);

  selected_entity_iterator<Mesh> begin(mesh, selected_buckets);
  selected_entity_iterator<Mesh> end;

  return std::make_pair(begin, end);
}

} // namespace details
} // namespace mesh
} // namespace sierra

#endif //SIERRA_SIERRA_MESH_SELECTED_ENTITIES_HPP

