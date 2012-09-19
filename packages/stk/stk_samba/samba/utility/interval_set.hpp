#ifndef SAMBA_SAMBA_UTILITY_INTERVAL_SET_HPP
#define SAMBA_SAMBA_UTILITY_INTERVAL_SET_HPP

#include <samba/utility/is_primitive.hpp>
#include <samba/utility/value_of.hpp>
#include <samba/utility/interval.hpp>

#ifdef SAMBA_ENABLE_PARALLEL
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#endif

#include <boost/icl/interval_set.hpp>

namespace samba {

/**
 * Need to wrap the boost inteval_set so that we can get an interval_set
 * implementation that works with samba primitives.
 */
template <typename Primitive>
struct interval_set
{
  BOOST_MPL_ASSERT((is_primitive<Primitive>));

  typedef Primitive value_type;
  typedef typename value_of<Primitive>::type integer_type;

  BOOST_MPL_ASSERT((boost::is_integral<integer_type>));

  typedef interval<value_type> interval_type;
  typedef size_t size_type;

  typedef interval_set<Primitive> self;

  private:
    typedef boost::icl::interval_set<integer_type> set_type;
    typedef typename set_type::interval_type set_interval_type;

    typedef typename set_type::element_const_iterator iterator_impl;

    struct to_interval
    {
      typedef interval_type result_type;

      result_type operator() ( set_interval_type i) const
      {
        return result_type(i.lower(), i.upper());
      }
    };

    struct to_value
    {
      typedef value_type result_type;
      value_type operator()(integer_type i) const
      {
        value_type v = value_type::create(i);
        return v;
      }
    };

  public:

  typedef boost::transform_iterator<
              to_value
            , iterator_impl
          > iterator;
  typedef iterator const_iterator;

  typedef boost::transform_iterator<
              to_interval
            , typename set_type::const_iterator
          > interval_iterator;

  typedef interval_iterator const_interval_iterator;

  interval_set()
    : m_set()
  {}

  interval_set( interval_type i)
    : m_set( set_interval_type(*(i.lower()),*(i.upper())))
  {}

  //empty

  value_type upper() const
  { value_type v = value_type::create(boost::icl::upper(m_set)); return v; }

  value_type lower() const
  { value_type v = value_type::create(boost::icl::lower(m_set)); return v; }

  bool empty() const
  { return m_set.empty(); }

  //number of elements
  size_type size() const
  { return m_set.size(); }

  //clear
  void clear() { m_set.clear(); }

  const_iterator begin() const
  {
    return const_iterator( boost::icl::elements_begin(m_set),
                          to_value()
                         );
  }

  const_iterator end() const
  {
    return const_iterator( boost::icl::elements_end(m_set),
                          to_value()
                         );
  }

  //number of intervals
  size_type interval_size() const { return m_set.iterative_size(); }

  const_interval_iterator interval_begin() const
  { return const_interval_iterator(m_set.begin(), to_interval()); }

  const_interval_iterator interval_end() const
  { return const_interval_iterator(m_set.end(), to_interval()); }

  //add element
  interval_set & operator += (value_type v)
  {
    m_set += v();
    return *this;
  }

  //add interval
  interval_set & operator += (interval_type i)
  {
    m_set += set_interval_type(i.lower()(),i.upper()());
    return *this;
  }

  //add set
  interval_set & operator += (self const& s)
  {
    m_set += s;
    return *this;
  }

  //remove element
  interval_set & operator -= (value_type v)
  {
    m_set -= v();
    return *this;
  }

  //remove interval
  interval_set & operator -= (interval_type i)
  {
    m_set -= set_interval_type(i.lower()(),i.upper()());
    return *this;
  }

  //remove set
  interval_set & operator -= (self const& s)
  {
    m_set -= s;
    return *this;
  }

  // comparsion operators
  // use lexigraphical ordering
  bool operator < (self const& rhs) const
  { return m_set < rhs.m_set; }
  bool operator <= (self const& rhs) const
  { return m_set <= rhs.m_set; }

  bool operator > (self const& rhs) const
  { return m_set > rhs.m_set; }
  bool operator >= (self const& rhs) const
  { return m_set >= rhs.m_set; }

  bool operator == (self const& rhs) const
  { return m_set == rhs.m_set; }
  bool operator != (self const& rhs) const
  { return m_set != rhs.m_set; }

  friend
  inline
  bool contains(self const& iset, value_type& v)
  {
    return boost::icl::contains(iset.m_set, v());
  }

  //swap sets
  void swap(self & s)
  { m_set.swap(s.m_set); }

  private:
    set_type m_set;

#ifdef SAMBA_ENABLE_PARALLEL
    friend class boost::serialization::access;

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const
    {
      std::vector<interval_type> vec;
      vec.reserve(interval_size());
      for (typename set_type::const_iterator i=m_set.begin(), e=m_set.end(); i != e; ++i)
        vec.push_back(interval_type(i->lower(),i->upper()));

      ar << vec;
    }

    template <typename Archive>
    void load(Archive &ar, const unsigned int version)
    {
      clear();

      std::vector<interval_type> vec;
      ar >> vec;

      for (size_t i=0; i<vec.size(); ++i)
        *this += vec[i];
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
};

template <typename Primitive>
inline
std::ostream& operator << (std::ostream& out, interval_set<Primitive> const& s)
{
  typedef interval_set<Primitive> set_type;
  typedef typename set_type::const_interval_iterator iterator_type;
  out << '{' << typename tag_of<Primitive>::type() << ":";
  if (!s.empty()) {
    for (iterator_type i=s.interval_begin(), e=s.interval_end(); i != e; ++i)
    {
      print_interval_range(out,*i);
      out << ",";
    }
    out << "\b}";
  }
  else {
    out << "empty}";
  }
  return out;
}

} // namespace samba

#endif //SAMBA_SAMBA_UTILITY_INTERVAL_SET_HPP
