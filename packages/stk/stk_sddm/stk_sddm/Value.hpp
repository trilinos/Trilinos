#ifndef STK_SDDM_VALUE_HPP
#define STK_SDDM_VALUE_HPP

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <iomanip>

#include <stk_sddm/DataTypeTraits.hpp>

namespace stk {
namespace sddm {

/**
 * @brief Function <b>cast_error</b> creates a <b>bad_any_data_cast</b>, message
 * describing the invalid conversion.
 *
 * @param of_type		a <b>std::type_info</b> const ...
 *
 * @param to_type		a <b>std::type_info</b> const ...
 *
 * @return			a <b>bad_any_data_cast</b> ...
 */
std::runtime_error cast_error(const char *of_type_name, const char *to_type_name);

/**
 * @brief Interface class <b>Data&lt;void&gt;</b> is used as the base class of all values.
 *
 */
class AnyValue
{
public:
  /**
   * Creates a new <b>AnyValue</b> instance.
   *
   */
  AnyValue()
  {}

  /**
   * Destroys a <b>AnyValue</b> instance.
   *
   */
  virtual ~AnyValue()
  {}

// protected:
//   template<typename T>
//   inline void setValue(const typename DataTypeTraits<T>::Type &value);

//   template<typename T>
//   inline void setText(const typename DataTypeTraits<T>::TextType &text);


public:
  /** 
   * @brief <b>value</b> 
   * 
   * 
   * @return 
   */
  template<typename T>
  inline const typename DataTypeTraits<T>::Type &value() const;

  template<typename T>
  inline const typename DataTypeTraits<T>::TextType &text() const;

  /**
   * @brief Pure virtual member function <b>type</b> returns the type of the value
   * stored in the <b>data</b> object.
   *
   * @return		a <b>std::type_info</b> const reference to the type
   */
  virtual const std::type_info &type() const = 0;

  virtual const char *typeName() const = 0;

  virtual std::ostream &dump(std::ostream &os) const = 0;

  virtual std::ostream &xml(std::ostream &os, size_t indent) const = 0;
};


/**
 * @brief template<b>Value</b> defines a resource which stores it's data by
 * value.
 *
 */
template<typename T>
class Value : public AnyValue
{
public:
  typedef Value<T> ValueType;
  typedef typename DataTypeTraits<T>::Type DataType;
  typedef std::string TextType;
  
  Value(const DataType &v)
    : AnyValue(),
      m_value(v)
  {}

  Value(const DataType &v, const TextType &t)
    : AnyValue(),
      m_value(v),
      m_text(t)
  {}

  /**
   * Destroys a <b>Value</b> instance.
   *
   */
  virtual ~Value()
  {}

  virtual const std::type_info &type() const {
    return typeid(DataType);
  }
  
  virtual const char *typeName() const {
    return DataTypeTraits<T>::name();
  }
  
  virtual const DataType &value() const {
    return m_value;
  }

  virtual const TextType &text() const {
    return m_text;
  }

  virtual std::ostream &dump(std::ostream &os) const {
    if (m_text.empty())
      return DataTypeTraits<T>::dump(os, m_value);
    else
      return os << m_text;
  }
  
  virtual std::ostream &xml(std::ostream &os, size_t indent) const {
    os << std::setw(indent*2) << "" << "<" << DataTypeTraits<T>::name() << ">";
    if (m_text.empty())
      DataTypeTraits<T>::dump(os, m_value);
    else
      os << m_text;
    os << "</" << DataTypeTraits<T>::name() << ">" << std::endl;

    return os;
  }
  
private:
  DataType      m_value;                ///< Data value
  TextType      m_text;                 ///< Text value
};


/**
 * @brief template<b>Value</b> defines a resource which stores it's data by
 * value.
 *
 */
template<>
class Value<std::string> : public AnyValue
{
public:
  typedef Value<std::string> ValueType;
  typedef DataTypeTraits<std::string>::Type DataType;
  typedef std::string TextType;
  
  Value(const DataType &v)
    : AnyValue(),
      m_value(v)
  {}

  Value(const DataType &v, const TextType &t)
    : AnyValue(),
      m_value(v),
      m_text(t)
  {
    if (m_value.size() > 1)
      if (m_value[0] == '\"' && m_value[m_value.size() - 1] == '\"')
        m_value = m_value.substr(1, m_value.size() - 2);
  }

  /**
   * Destroys a <b>Value</b> instance.
   *
   */
  virtual ~Value()
  {}

  virtual const std::type_info &type() const {
    return typeid(DataType);
  }
  
  virtual const char *typeName() const {
    return DataTypeTraits<std::string>::name();
  }
  
  virtual const DataType &value() const {
    return m_value;
  }

  virtual const TextType &text() const {
    return m_text;
  }

  virtual std::ostream &dump(std::ostream &os) const {
    os << "\"";
    if (m_text.empty())
      DataTypeTraits<std::string>::dump(os, m_value);
    else
      os << m_text;
    os << "\"";
    
    return os;
  }
  
  virtual std::ostream &xml(std::ostream &os, size_t indent) const {
    os << std::setw(indent*2) << "" << "<" << DataTypeTraits<std::string>::name() << ">";
    if (m_text.empty())
      DataTypeTraits<std::string>::dump(os, m_value);
    else
      os << m_text;
    os << "</" << DataTypeTraits<std::string>::name() << ">" << std::endl;

    return os;
  }
  
private:
  DataType      m_value;                ///< Data value
  TextType      m_text;                 ///< Text value
};


template<typename T>
class Value<std::vector<T> > : public AnyValue
{
public:
  typedef Value<std::vector<T> > ValueType;
  typedef typename DataTypeTraits<std::vector<T> >::Type DataType;
  typedef std::vector<std::string> TextType;
  
  Value(const DataType &v)
    : AnyValue(),
      m_value(v)
  {}

  Value(const DataType &v, const TextType &t)
    : AnyValue(),
      m_value(v),
      m_text(t)
  {}

  /**
   * Destroys a <b>Value</b> instance.
   *
   */
  virtual ~Value()
  {}

  virtual const std::type_info &type() const {
    return typeid(DataType);
  }
  
  virtual const char *typeName() const {
    return DataTypeTraits<std::vector<T> >::name();
  }
  
  virtual const DataType &value() const {
    return m_value;
  }

  virtual const TextType &text() const {
    return m_text;
  }

  virtual std::ostream &dump(std::ostream &os) const {
    if (m_text.empty()) {
      for (typename DataType::const_iterator it = m_value.begin(); it != m_value.end(); ++it) {
        if (it != m_value.begin())
          os << " ";
        DataTypeTraits<T>::dump(os, (*it));
      }
    }
    
    else {
      for (TextType::const_iterator it = m_text.begin(); it != m_text.end(); ++it) {
        if (it != m_text.begin())
          os << " ";
        os << (*it);
      }
    }
    
    return os;
  }
  
  virtual std::ostream &xml(std::ostream &os, size_t indent) const {
    os << std::setw(indent*2) << "" << "<" << DataTypeTraits<std::vector<T> >::name() << ">" << std::endl;

    if (m_text.empty()) {
      for (typename DataType::const_iterator it = m_value.begin(); it != m_value.end(); ++it) {
        os << std::setw(indent*2 + 2) << "" << "<" << DataTypeTraits<T>::name() << ">";
        DataTypeTraits<T>::dump(os, (*it));
        os << "</" << DataTypeTraits<T>::name() << ">" << std::endl;
      }
    }
    
    else {
      for (TextType::const_iterator it = m_text.begin(); it != m_text.end(); ++it) {
        os << std::setw(indent*2 + 2) << "" << "<" << DataTypeTraits<T>::name() << ">";
        os << (*it);
        os << "</" << DataTypeTraits<T>::name() << ">";
      }
    }
    os << std::setw(indent*2) << "" << "</" << DataTypeTraits<std::vector<T> >::name() << ">" << std::endl;
      
    return os;
  }
  
private:
  DataType      m_value;                ///< Data value
  TextType      m_text;                 ///< Text value
};

template<>
class Value<std::vector<std::string> > : public AnyValue
{
public:
  typedef Value<std::vector<std::string> > ValueType;
  typedef DataTypeTraits<std::vector<std::string> >::Type DataType;
  typedef std::vector<std::string> TextType;
  
  Value(const DataType &v)
    : AnyValue(),
      m_value(v)
  {}

  Value(const DataType &v, const TextType &t)
    : AnyValue(),
      m_value(v),
      m_text(t)
  {}

  /**
   * Destroys a <b>Value</b> instance.
   *
   */
  virtual ~Value()
  {}

  virtual const std::type_info &type() const {
    return typeid(DataType);
  }
  
  virtual const char *typeName() const {
    return DataTypeTraits<std::vector<std::string> >::name();
  }
  
  virtual const DataType &value() const {
    return m_value;
  }

  virtual const TextType &text() const {
    return m_text;
  }

  virtual std::ostream &dump(std::ostream &os) const {
    if (m_text.empty()) {
      for (DataType::const_iterator it = m_value.begin(); it != m_value.end(); ++it) {
        if (it != m_value.begin())
          os << " ";
        os << "\"" << (*it) << "\"";
      }
    }
    
    else {
      for (TextType::const_iterator it = m_text.begin(); it != m_text.end(); ++it) {
        if (it != m_text.begin())
          os << " ";
        os << "\"" << (*it) << "\"";
      }  
    }
    
    return os;
  }
  
  virtual std::ostream &xml(std::ostream &os, size_t indent) const {
    os << std::setw(indent*2) << "" << "<" << DataTypeTraits<std::vector<std::string> >::name() << ">" << std::endl;

    if (m_text.empty())
      for (DataType::const_iterator it = m_value.begin(); it != m_value.end(); ++it) {
        os << std::setw(indent*2 + 2) << "" << "<" << DataTypeTraits<std::string>::name() << ">";
        DataTypeTraits<std::string>::dump(os, (*it));
        os << "</" << DataTypeTraits<std::string>::name() << ">" << std::endl;
      }
    else {
      for (TextType::const_iterator it = m_text.begin(); it != m_text.end(); ++it) {
        os << std::setw(indent*2 + 2) << "" << "<" << DataTypeTraits<std::string>::name() << ">";
        os << (*it);
        os << "</" << DataTypeTraits<std::string>::name() << ">" << std::endl;
      }
    }
    os << std::setw(indent*2) << "" << "</" << DataTypeTraits<std::vector<std::string> >::name() << ">" << std::endl;
    return os;
  }
  
private:
  DataType      m_value;                ///< Data value
  TextType      m_text;                 ///< Text value
};

/**
 * @brief Member function <b>value</b> attempts to cast the data to the
 * specified type.  If the data is not of the specified type, an exception is thrown.
 *
 * @return		a <b>T</b> reference to the data.
 *
 * @throws		a <b>std::runtime_error</b> exception with description.
 */
template<typename T>
inline const typename DataTypeTraits<T>::Type &
AnyValue::value() const
{
  if (typeid(typename DataTypeTraits<T>::Type) != type())
    throw cast_error(type().name(), DataTypeTraits<T>::name());
  
  const typename Value<T>::ValueType *t = static_cast<const typename Value<T>::ValueType *>(this);
  return t->value();
}


// /**
//  * @brief Member function <b>value</b> attempts to cast the data to the
//  * specified type.  If the data is not of the specified type, an exception is thrown.
//  *
//  * @return		a <b>T</b> reference to the data.
//  *
//  * @throws		a <b>std::runtime_error</b> exception with description.
//  */
// template<typename T>
// inline void
// AnyValue::setValue(const typename DataTypeTraits<T>::Type &t)
// {
//   if (typeid(DataTypeTraits<T>::Type) != type())
//     throw cast_error(typeName(), DataTypeTraits<T>::name());

//   typename Value<T>::Type *data = static_cast<typename Value<T>::Type *>(this);
//   data->value() = t;
// //  DataTypeTraits<T>::dump(t, data->value());
// }


// /**
//  * @brief Member function <b>text</b> attempts to cast the data to the
//  * specified type.  If the data is not of the specified type, an exception is thrown.
//  *
//  * @return		a <b>T</b> reference to the data.
//  *
//  * @throws		a <b>std::runtime_error</b> exception with description.
//  */
// template<typename T>
// inline void
// AnyValue::setText(const typename DataTypeTraits<T>::TextType &t)
// {
//   if (typeid(DataTypeTraits<T>::Type) != type())
//     throw cast_error(type(), DataTypeTraits<T>::name());
  
//   typename Value<T>::Type *data = static_cast<const typename Value<T>::Type *>(this);
//   data->text() = t;
// //  DataTypeTraits<T>::load(t, data->value());
// }

} // namespace sddm
} // namespace stk

#endif // STK_SDDM_VALUE_HPP
