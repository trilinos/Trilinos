#ifndef STK_SDDM_STATEVALUE_HPP
#define STK_SDDM_STATEVALUE_HPP

#include <stk_sddm/Value.hpp>
#include <stk_sddm/Property.hpp>

std::runtime_error state_error(const std::string &state);

/**
 * @brief template<b>StateValue</b> defines a resource which stores
 * it's data by value indexed by the state of another property.
 *
 */
template<typename T>
class StateValue : public AnyValue
{
public:
  typedef typename DataTypeTraits<T>::Type DataType;
  typedef std::map<std::string, DataType> StateMap;
  
  /**
   * Creates a new <b>Value</b> instance.
   *
   */
  StateValue(const Property &state_property)
    : AnyValue(),
      m_stateProperty(state_property),
      m_default()
  {}

  StateValue(const Property &state_property, const DataType &default_t)
    : AnyValue(),
      m_stateProperty(state_property),
      m_default(default_t)
  {}

  /**
   * Destroys a <b>StateValue</b> instance.
   *
   */
  virtual ~StateValue()
  {}

  virtual const std::type_info &type() const {
    return typeid(DataType);
  }
  
  virtual const char *name() const {
    return DataTypeTraits<T>::name();
  }
  
  virtual const DataType &value() const;

  virtual DataType &value();

  StateValue &setValue(const std::string &state, const DataType &t) {
    m_stateMap[state] = t;
    return *this;
  }
  
  const Property &getStateProperty() {
    return m_stateProperty;
  }
  
  DataType &value(const std::string &state) {
    typename StateMap::iterator it = m_stateMap.find(state);
    if (it != m_stateMap.end())
      return (*it).second;

    throw state_error(state);
  }
  
  virtual std::ostream &dump(std::ostream &os) const {
    return DataTypeTraits<T>::dump(os, m_default);
  }
  
private:
  const Property &      m_stateProperty;        ///< Property of current state
  StateMap              m_stateMap;             ///< Mapping from state to value
  DataType              m_default;              ///< Default value
};

template <class T>
const typename StateValue<T>::DataType &
StateValue<T>::value() const
{
  typename StateMap::const_iterator it = m_stateMap.find(m_stateProperty.value<std::string>());
  if (it != m_stateMap.end())
    return (*it).second;
    
  return m_default;
}

template <class T>
typename StateValue<T>::DataType &
StateValue<T>::value()
{
  typename StateMap::iterator it = m_stateMap.find(m_stateProperty.value<std::string>());
  if (it != m_stateMap.end())
    return (*it).second;
    
  return m_default;
}

#endif // #define STK_SDDM_STATEVALUE_HPP
