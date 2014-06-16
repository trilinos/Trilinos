#ifndef SAMBA_SAMBA_MESH_PARTITION_ITERATOR_IMPL_TCC
#define SAMBA_SAMBA_MESH_PARTITION_ITERATOR_IMPL_TCC

#include <boost/iterator/iterator_facade.hpp>

namespace samba {  namespace detail {

class partition_iterator
  : public boost::iterator_facade< partition_iterator
                                  ,entity_proxy
                                  ,boost::random_access_traversal_tag
                                  ,entity_proxy
                                  ,size_t
                                 >
{
  public:

    partition_iterator(partition_impl const* arg_partition, partition_offset arg_offset)
      : m_partition(arg_partition)
      , m_offset(arg_offset)
    {}

  private:
    friend class boost::iterator_core_access;

    entity_proxy dereference() const
    { return entity_proxy(m_partition,m_offset); }

    bool equal(partition_iterator rhs) const
    {
      return    m_partition == rhs.m_partition
             && m_offset == rhs.m_offset;
    }

    void increment()
    { ++m_offset; }

    void decrement()
    { ++m_offset; }

    void advance(int n)
    { m_offset += n; }

    size_t distance_to(partition_iterator rhs) const
    { return rhs.m_offset() - m_offset(); }


  private:
    partition_impl const* m_partition;
    partition_offset         m_offset;
};

//*****************************************************************************
//defined here to break cycles
//*****************************************************************************
inline partition_iterator partition_impl::begin() const
{ return partition_iterator(this,partition_offset::create(0)); }

inline partition_iterator partition_impl::end() const
{ return partition_iterator(this,partition_offset::create(size())); }

} //namespace detail

inline detail::partition_iterator partition_proxy::begin() const
{ return detail::partition_iterator(this->m_partition,partition_offset::create(0)); }

inline detail::partition_iterator partition_proxy::end() const
{ return detail::partition_iterator(this->m_partition,partition_offset::create(size())); }

} //namespace samba

#endif //SAMBA_SAMBA_MESH_PARTITION_ITERATOR_IMPL_TCC
