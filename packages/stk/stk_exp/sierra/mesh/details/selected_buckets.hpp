#ifndef SIERRA_SIERRA_MESH_SELECTEDBUCKETS_HPP
#define SIERRA_SIERRA_MESH_SELECTEDBUCKETS_HPP

#include <sierra/mesh/details/selector.hpp>
#include <sierra/mesh/mesh_traits.hpp>

#include <boost/range.hpp>

namespace sierra {
namespace mesh {
namespace details {

//selected_bucket_iterator is an iterator that uses a selector to decide
//whether the next bucket in a sequence should be used.
//Users shouldn't construct this class directly, they should call
//the get_selected_buckets function below.
template <class Mesh>
class selected_bucket_iterator
  : public std::iterator< std::forward_iterator_tag,
                          bucket_key,
                          ptrdiff_t,
                          const bucket_key *,
                          bucket_key
                        >
{
 public:
  selected_bucket_iterator( const Mesh& mesh,
                            const selector & select,
                            typename mesh_traits<Mesh>::bucket_iterator arg_itr,
                            typename mesh_traits<Mesh>::bucket_iterator arg_end
                          )
   : m_mesh(&mesh)
   , m_select(select)
   , m_itr(arg_itr)
   , m_end(arg_end)
  {
    if ( m_itr != m_end && !contains(m_select,*m_itr,*m_mesh) ) {
      update();
    }
  }

  selected_bucket_iterator()
   : m_mesh(NULL)
   , m_select()
   , m_itr()
   , m_end()
  {
  }

  bool operator!=(const selected_bucket_iterator& itr) const {
    return m_itr != itr.m_itr;
  }

  bool operator==(const selected_bucket_iterator& itr) const {
    return m_itr == itr.m_itr;
  }

  selected_bucket_iterator operator++() {
    update();
    return *this;
  }

  selected_bucket_iterator operator++(int) {
    selected_bucket_iterator temp(*this);
    update();
    return temp;
  }

  bucket_key operator *() {
    return *m_itr;
  }

 private:

  void update() {
    do {
      ++m_itr;
    }
    while (m_itr != m_end && !contains(m_select,*m_itr,*m_mesh));
  }

  const Mesh           * m_mesh;
  selector               m_select;
  typename mesh_traits<Mesh>::bucket_iterator   m_itr;
  typename mesh_traits<Mesh>::bucket_iterator   m_end;
};

template <class Mesh>
inline
std::pair<selected_bucket_iterator<Mesh>, selected_bucket_iterator<Mesh> >
get_selected_buckets(const selector& select, const Mesh& mesh)
{
  return std::make_pair(
      selected_bucket_iterator<Mesh>(mesh, select, boost::begin(mesh.get_buckets()), boost::end(mesh.get_buckets())),
      selected_bucket_iterator<Mesh>(mesh, select, boost::end(mesh.get_buckets()), boost::end(mesh.get_buckets())));
}

} // namespace details
} // namespace mesh
} // namespace sierra

#endif //SIERRA_SIERRA_MESH_SELECTEDBUCKETS_HPP
