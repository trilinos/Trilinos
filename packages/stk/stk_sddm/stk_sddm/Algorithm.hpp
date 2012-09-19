#ifndef STK_SDDM_ALGORITHM_HPP
#define STK_SDDM_ALGORITHM_HPP

#include <iterator>
#include <map>
#include <string>
#include <vector>

#include <stk_sddm/Property.hpp>

namespace stk {
namespace sddm {

struct find_pred {
  
  static std::string name(const std::string &s) {
    std::string::size_type i = s.find("[");
    
    if (i != std::string::npos) {
      return s.substr(0, i);      
    }
    else
      return s;
  }
  
  find_pred(const std::string &query)
    : m_name(name(query))
  {}

  bool operator()(const Property *property) {
    const size_t property_name_size = property->getName().size();
    const size_t name_size = m_name.size();
    const size_t len = std::min(property_name_size, name_size);
    
    int r = ignorecase_traits::compare(property->getName().data(), m_name.data(), len);
    if (!r)
      r =  property_name_size - name_size;
    return r == 0;
  }
    
  const std::string     m_name;
  const std::string     m_;
};
  
// Copy_if was dropped from the standard library by accident.
template<typename In, typename Out, typename Pred>
Out copy_if(In first, In last, Out res, Pred Pr)
{
  while (first != last) {
    if (Pr(*first))
      *res++ = *first;
    ++first;
  }

  return res;
}


class counter_iterator : public std::iterator<std::output_iterator_tag, void, void, void, void>
{
protected:
  size_t &      m_counter;

public:
  explicit counter_iterator(size_t &counter)
    : m_counter(counter)
  {}

  template <class T>
  counter_iterator &operator=(const T &t) {
    ++m_counter;
    return *this;
  }
  
  counter_iterator &operator*() {
    return *this;
  }

  counter_iterator &operator++() {
    return *this;
  }

 counter_iterator &operator++(int) {
   return *this;
 }
};

inline counter_iterator counter(size_t &counter) {
  return counter_iterator(counter);
}

template <class It, class Out>
void
query(
  It                            first,
  It                            last,
  Out                           out,
  const std::string &           query)
{
  find_pred pred(query);
  
  copy_if(first, last, out, pred);
}

void remove(PropertyVector &property_vector, const std::string &name);
void remove(ConstPropertyVector &property_vector, const std::string &name);

template <class It, class Key, class Value>
void
map_properties(
  It                            begin,
  It                            end,
  const std::string &           key_name,
  const std::string &           value_name,
  std::map<Key, Value> &        map)
{
  typedef std::map<Key, Value> Map;

  for (It it = begin; it != end; ++it) {
    const stk::sddm::Property &property = *(*it);
    map[property.value<Key>(key_name)] = property.value<Value>(value_name);
  }
}

template <class T>
std::pair<Property *, bool>
setDefaultProperty(
  Property &            parent_property,
  const std::string &   name,
  const T &             t,
  const char *          function) 
{
  Property *property = parent_property.find(name);
  if (property)
    return std::pair<Property *, bool>(property, false);
  else {
    property = &parent_property.create(name, t);
    property->addAttribute(new DefaultAttribute());
    property->addAttribute(new MethodSourceAttribute(function));
    return std::pair<Property *, bool>(property, true);
  }
}

bool isDefaultProperty(const Property &parent_property, const std::string &name);

} // namespace sddm
} // namespace stk

#endif // STK_SDDM_ALGORITHM_HPP
