#ifndef SDDM_TYPE_HPP
#define SDDM_TYPE_HPP

#include <iosfwd>
#include <vector>
#include <typeinfo>

namespace stk {
namespace sddm {

/**
 * @brief Interface class <b>Data&lt;void&gt;</b> is used as the base class of all type descriptors.
 *
 */
class AnyType
{
public:
  /**
   * Creates a new <b>AnyType</b> instance.
   *
   */
  AnyType()
  {}

  /**
   * Destroys a <b>AnyType</b> instance.
   *
   */
  virtual ~AnyType()
  {}

  /**
   * @brief Pure virtual member function <b>type</b> returns the type of the value
   * stored in the <b>data</b> object.
   *
   * @return		a <b>std::type_info</b> const reference to the type
   */
  virtual const std::type_info &type() const = 0;

  virtual const char *name() const = 0;

  virtual std::ostream &dump(std::ostream &os) const = 0;
};


/**
 * @brief template<b>Value</b> defines a resource which stores it's data by
 * value.
 *
 */
template<typename T>
class Type : public AnyType
{
public:
  typedef typename DataTypeTraits<T>::Type DataType;
  
  static AnyType *instance();

  /**
   * Creates a new <b>Type</b> instance.
   *
   */
  Type()
    : AnyType()
  {}

  Type(const Type &t)
    : AnyType()
  {}

  /**
   * Destroys a <b>Type</b> instance.
   *
   */
  virtual ~Type()
  {}

  virtual const std::type_info &type() const {
    return typeid(DataType);
  }
  
  virtual const char *name() const {
    return DataTypeTraits<T>::name();
  }
  
  virtual std::ostream &dump(std::ostream &os) const {
    return os << DataTypeTraits<T>::name();
  }
};

/**
 * @brief template<b>Value</b> defines a resource which stores it's data by
 * value.
 *
 */
template<typename T>
class Type<std::vector<T> > : public AnyType
{
public:
  typedef std::vector<typename DataTypeTraits<T>::Type> DataType;
  
  static AnyType *instance();

  /**
   * Creates a new <b>Type</b> instance.
   *
   */
  Type()
    : AnyType()
  {}

  Type(const Type &t)
    : AnyType()
  {}

  /**
   * Destroys a <b>Type</b> instance.
   *
   */
  virtual ~Type()
  {}

  virtual const std::type_info &type() const {
    return typeid(DataType);
  }
  
  virtual const char *name() const {
    return DataTypeTraits<std::vector<T> >::name();
  }
  
  virtual std::ostream &dump(std::ostream &os) const {
    return os << DataTypeTraits<std::vector<T> >::name();
  }
};

inline
std::ostream &operator<<(std::ostream &os, const AnyType &any_type) {
  return any_type.dump(os);
}

std::vector<const AnyType *> &defaultTypes(std::vector<const AnyType *> &default_types);

} // namespace sddm
} // namespace stk

#endif // SDDM_TYPE_HPP
