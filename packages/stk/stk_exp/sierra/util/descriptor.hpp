#ifndef SIERRA_SIERRA_UTIL_DESCRIPTOR_HPP
#define SIERRA_SIERRA_UTIL_DESCRIPTOR_HPP

#include <cstddef>

namespace sierra {
namespace util {

struct parameters_by_value {};
struct parameters_by_ref {};

template <class DescriptorType>
struct descriptor_traits
{
  typedef typename DescriptorType::base_type          base_type;
  typedef typename DescriptorType::parameter_catagory parameter_catagory;

  static const base_type invalid() { return DescriptorType::invalid(); }
};


template <
          class DescriptorType,
          class BaseType = typename descriptor_traits<DescriptorType>::base_type,
          class ParameterType = typename descriptor_traits<DescriptorType>::parameter_catagory
         >
class descriptor;

template <class DescriptorType, class BaseType >
class descriptor<DescriptorType,BaseType, parameters_by_value>
{
  public:
    typedef BaseType base_type;
    typedef parameters_by_value parameter_catagory;
    static const base_type invalid() { return descriptor_traits<DescriptorType>::invalid(); }

    explicit
    descriptor(base_type base = invalid() )
      : m_base(base)
    {}

    bool operator==(descriptor rhs) const
    { return m_base == rhs.m_base; }

    bool operator!=(descriptor rhs) const
    { return m_base != rhs.m_base; }

    bool operator<(descriptor rhs) const
    { return m_base < rhs.m_base;  }

    bool operator<=(descriptor rhs) const
    { return m_base <= rhs.m_base; }

    bool operator>(descriptor rhs) const
    { return m_base > rhs.m_base;  }

    bool operator>=(descriptor rhs) const
    { return m_base >= rhs.m_base; }

    operator base_type() const
    { return m_base; }

    base_type base() const { return m_base; }

  private:
      base_type m_base;
};

template <class DescriptorType, class BaseType >
class descriptor<DescriptorType,BaseType, parameters_by_ref>
{
  public:
    typedef BaseType base_type;
    typedef parameters_by_ref parameter_catagory;
    static const base_type invalid() { return descriptor_traits<DescriptorType>::invalid(); }

    explicit
    descriptor(const base_type & base = invalid() )
      : m_base(base)
    {}

    bool operator==(const descriptor & rhs) const
    { return m_base == rhs.m_base; }

    bool operator!=(const descriptor & rhs) const
    { return m_base != rhs.m_base; }

    bool operator<(const descriptor & rhs) const
    { return m_base < rhs.m_base;  }

    bool operator<=(const descriptor & rhs) const
    { return m_base <= rhs.m_base; }

    bool operator>(const descriptor & rhs) const
    { return m_base > rhs.m_base;  }

    bool operator>=(const descriptor & rhs) const
    { return m_base >= rhs.m_base; }

    operator base_type() const
    { return m_base; }

    const base_type & base() const { return m_base; }

  private:
      base_type m_base;
};


} // namespace util
} // namespace sierra

#endif // SIERRA_SIERRA_UTIL_DESCRIPTOR_HPP
