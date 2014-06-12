#ifndef STKTOPOLOGY_DETAIL_FILL_CONTAINER_HPP
#define STKTOPOLOGY_DETAIL_FILL_CONTAINER_HPP

#include <vector>

// functors used with boost::mpl::for_each
// to extract values from boost::mpl::vectors
namespace stk { namespace topology_detail {

template <typename OrdinalOutputIterator>
struct fill_ordinal_container {

  template <typename Ordinal>
  BOOST_GPU_ENABLED
  void operator()(Ordinal i)
  { *m_itr = i; ++m_itr; }

  BOOST_GPU_ENABLED
  fill_ordinal_container( OrdinalOutputIterator itr)
    : m_itr(itr)
  {}

  OrdinalOutputIterator m_itr;
};

template <typename T, typename A>
struct fill_ordinal_container< std::vector<T,A> >
{
  template <typename Ordinal>
  BOOST_GPU_ENABLED
  void operator()(Ordinal i)
  { *m_itr = i; ++m_itr; }

  BOOST_GPU_ENABLED
  fill_ordinal_container( std::vector<T,A> & vec)
    : m_itr(vec.begin())
  {}

  typename std::vector<T,A>::iterator m_itr;

};

template <typename NodeArray, typename NodeOutputIterator>
struct fill_node_container {

  template <typename Ordinal>
  BOOST_GPU_ENABLED
  void operator()(Ordinal i)
  { *m_itr = m_nodes[i]; ++m_itr; }

  BOOST_GPU_ENABLED
  fill_node_container( const NodeArray & nodes, NodeOutputIterator itr)
    : m_nodes(nodes)
    , m_itr(itr)
  {}

  const NodeArray    & m_nodes;
  NodeOutputIterator   m_itr;
};

template <typename NodeArray, typename T, typename A>
struct fill_node_container<NodeArray, std::vector<T,A> > {

  template <typename Ordinal>
  BOOST_GPU_ENABLED
  void operator()(Ordinal i)
  { *m_itr = m_nodes[i]; ++m_itr; }

  BOOST_GPU_ENABLED
  fill_node_container( const NodeArray & nodes, std::vector<T,A> & vec)
    : m_nodes(nodes)
    , m_itr(vec.begin())
  {}

  const NodeArray    & m_nodes;
  typename std::vector<T,A>::iterator   m_itr;
};

}} // namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_FILL_UNSIGNED_CONTAINER_HPP

