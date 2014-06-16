#ifndef SIERRA_SIERRA_UTIL_DESCRIPTOR_MANAGER_HPP
#define SIERRA_SIERRA_UTIL_DESCRIPTOR_MANAGER_HPP

#include <sierra/util/descriptor.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace sierra {
namespace util {

template <class DescriptorType>
struct descriptor_manager_traits
{
  typedef descriptor<DescriptorType>                             descriptor_type;
  typedef typename descriptor_traits<DescriptorType>::base_type  base_type;

  static const base_type first_valid() { return DescriptorType::first_valid(); }
  static const base_type last_valid()  { return DescriptorType::last_valid(); }
};


template <
          class DescriptorType,
          class BaseType = typename descriptor_manager_traits<DescriptorType>::base_type
         >
class descriptor_manager
{
  private:
    typedef boost::icl::discrete_interval< BaseType > BaseInterval;
    typedef boost::icl::interval_set< BaseType > BaseSet;
    typedef typename BaseSet::element_const_iterator  base_iterator;

  public:
    typedef BaseType base_type;
    typedef typename descriptor_manager_traits<DescriptorType>::descriptor_type descriptor_type;

  private:
    struct to_descriptor
      : std::unary_function<base_type, descriptor_type>
    {
      descriptor_type operator() ( const base_type & base) const {
        return descriptor_type(base);
      }
    };

  public:
    static const base_type first_valid() { return descriptor_manager_traits<DescriptorType>::first_valid(); }
    static const base_type last_valid() { return descriptor_manager_traits<DescriptorType>::last_valid(); }

    typedef boost::transform_iterator<
                                       to_descriptor,
                                       base_iterator
                                     > iterator;

    typedef std::pair<iterator,iterator> iterator_range;

    descriptor_manager()
      : m_all_descriptors(BaseInterval(first_valid(),last_valid()))
      , m_used_descriptors()
      , m_unused_descriptors()
    {
      m_unused_descriptors += m_all_descriptors;
    }

    bool add( descriptor_type d) {
      if (!in_use(d)) {
        m_unused_descriptors -= base_type(d);
        m_used_descriptors += base_type(d);
        return true;
      }
      return false;
    }

    bool remove( descriptor_type d) {
      if (in_use(d)) {
        m_unused_descriptors += base_type(d);
        m_used_descriptors -= base_type(d);
        return true;
      }
      return false;
    }

    descriptor_type create() {
      base_type b = boost::icl::first(m_unused_descriptors);
      m_unused_descriptors -= b;
      m_used_descriptors += b;
      return descriptor_type(b);
    }

    bool in_use( descriptor_type d ) const {
      return (boost::icl::contains(m_used_descriptors, base_type(d) ));
    }

    size_t size() const {
      return m_used_descriptors.size();

    }

    iterator begin() const { return iterator(boost::icl::elements_begin(m_used_descriptors), to_descriptor()); }
    iterator end()   const { return iterator(boost::icl::elements_end(m_used_descriptors), to_descriptor()); }

    iterator_range range() const {
      return std::make_pair( begin(), end() );
    }

  private:
    const BaseSet  m_all_descriptors;
    BaseSet        m_used_descriptors;
    BaseSet        m_unused_descriptors;

    descriptor_manager(const descriptor_manager&);
    descriptor_manager & operator=(const descriptor_manager&);
};

} // namespace util
} // namespace sierra

#endif // SIERRA_SIERRA_UTIL_DESCRIPTOR_MANAGER_HPP
