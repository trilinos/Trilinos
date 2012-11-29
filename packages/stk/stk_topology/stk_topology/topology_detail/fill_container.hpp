#ifndef STKTOPOLOGY_DETAIL_FILL_CONTAINER_HPP
#define STKTOPOLOGY_DETAIL_FILL_CONTAINER_HPP


// functors used with boost::mpl::for_each
// to extract values from boost::mpl::vectors
namespace stk { namespace topology_detail {

template <typename OrdinalOutputIterator>
struct fill_ordinal_container {

  template <typename Ordinal>
  void operator()(Ordinal i)
  { *m_itr = i; ++m_itr; }

  fill_ordinal_container( OrdinalOutputIterator itr)
    : m_itr(itr)
  {}

  OrdinalOutputIterator m_itr;
};

template <typename NodeArray, typename NodeOutputIterator>
struct fill_node_container {

  template <typename Ordinal>
  void operator()(Ordinal i)
  { *m_itr = m_nodes[i]; ++m_itr; }

  fill_node_container( const NodeArray & nodes, NodeOutputIterator itr)
    : m_nodes(nodes)
    , m_itr(itr)
  {}

  const NodeArray    & m_nodes;
  NodeOutputIterator   m_itr;
};


}} // namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_FILL_UNSIGNED_CONTAINER_HPP

