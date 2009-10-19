#ifndef STK_UTIL_UTIL_MARSHAL_HPP
#define STK_UTIL_UTIL_MARSHAL_HPP

#include <stdint.h>

#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <typeinfo>

namespace stk {

/**
 * @brief Struct <code>Marshal</code> is a data packer for sending and receiving parallel messages.
 * The data put-to (<<) is appended to the stream as a string of bytes, likewise data gotten-from
 * (>>) is extracted from the stream into the object as a string of bytes.
 *
 * The write() and read() functions perform the data movements to and from the packed stream.
 *
 * The common implementation is the create a << and >> operator for an object which properly appends
 * and extracts the object's members.
 *
 * The object can put-to and get-from it's typeid() to add type checking.  This operation ensures
 * that the data types being read was the data type written before the data is extracted.  This type
 * checking can be disabled since it may be desired to put-to an object of one type, but get-from
 * into an object of an extractable but different type.
 *
 * The TYPE_CHECK bit masks can be provided at put-to Marshal construction to activate the type
 * checking.  The Marshaller send the type check code as the first message to allow the get-from to
 * initialize properly.
 *
 * The put-to operator and get-from operators for plain old data, std::string, std::vector and
 * std::list have been implemented.  Additional ones could be added here, or left to the developer
 * using the marshaller.
 *
 * The stream and type_check members were left as public due to the extensive use.  If this proves
 * bothersome, getter/setter methods could be introduced.
 *
 */
struct Marshal
{
  /**
   * @brief Enumeration to activate type checking for std classes and plain old data.
   *
   */
  enum {
    TYPE_CHECK_NONE     = 0x00000000,
    TYPE_CHECK_POD      = 0x00000001,
    TYPE_CHECK_LIST     = 0x00000002,
    TYPE_CHECK_VECTOR   = 0x00000004,
    TYPE_CHECK_ALL      = 0xFFFFFFFF
  };
  
  /**
   * Creates a new <code>Marshal</code> instance for put-to operations.
   *
   */
  Marshal(unsigned type_check = TYPE_CHECK_NONE);

  /**
   * Creates a new <code>Marshal</code> instance for get-from operations.
   *
   * @param s			a <code>std::string</code> constant variable of packed bytes to
   *                            extract using the get-from operators.
   */
  explicit Marshal(const std::string &s);

  /**
   * @brief Member function <code>str</code> returns the string of packed bytes created by put-to
   * operations to the stream.
   *
   * @return			a <code>std::string</code> created from the packed byte stream. 
   */
  std::string str() const;

  /**
   * @brief Member function <code>size</code> returns the byte count of the string of packed bytes
   * creates by put-to operations to the stream.
   *
   * @return			a <code>size_t</code> in bytes of the packed byte stream.
   */
  size_t size() const;

  /**
   * @brief Member function <code>write</code> writer bytes to the packed byte stream.
   *
   * @param byte_count		a <code>size_t</code> value of the number of packed bytes to write.
   *
   * @param address		a <code>char</code> constant pointer to get the bytes from.
   *
   */
  void write(const char *address, size_t byte_count);
  
  /**
   * @brief Member function <code>read</code> reads bytes from the packed byte stream.
   *
   * @param byte_count		a <code>size_t</code> value of the number of packed bytes to read.
   *
   * @param address		a <code>char</code> constant pointer to put the bytes to.
   *
   */
  void read(char *address, size_t byte_count);

  /**
   * @brief Member function <code>operator void *</code> returns the state of the packed byte stream.
   *
   * @return			a <code>void</code> const pointer which is non-zero if status is
   *                            good. 
   */
  operator void * () const;
  
private:
  Marshal(const Marshal &marshal);                      ///< Not copyable
  Marshal &operator=(const Marshal &);                  ///< Not assignable

public:
  std::stringstream     stream;                         ///< Packed byte stream to put-to or get-from
  unsigned              m_typeCheck;                    ///< Type checking to activate
};


/**
 * @brief Function <code>operator<< </code> writes the object to the packed byte stream.  This is
 * the template class and has no implementation.  You must specialize this class to write an
 * object. 
 *
 * @param mout  		a <code>Marshal</code> reference to the marshaller.
 *
 * @param t       		a <code>T</code> const reference to the object to write.
 *
 * @return			a <code>Marshal</code> reference to the marhsaller.
 */
template <typename T>
Marshal &operator<<(Marshal &mout, const T &t);

/**
 * @brief Function <code>operator>> </code> reads the object from the packed byte stream.  This is
 * the template class and has no implementation.  You must specialize this class to read an object.
 *
 * @param min     		a <code>Marshal</code> reference to the marshaller.
 *
 * @param t       		a <code>T</code> const reference to the object to read.
 *
 * @return			a <code>Marshal</code> reference to the marhsaller.
 */
template <typename T>
Marshal &operator>>(Marshal &min, T &t);

/**
 * @brief Function <code>operator<< </code> write the crc32 encoding of the name from the type
 * information to the packed byte stream.  When the bytes are read, the crc32 encoding of the type
 * being read is varified.
 *
 * @param mout  		a <code>Marshal</code> reference to the marshaller.
 *
 * @param t       		a <code>std::type_info</code> const reference to the type
 *                              information to write for verification when read on extraction. 
 *
 * @return			a <code>Marshal</code> reference to the marhsaller.
 */
template<>
Marshal &operator<<(Marshal &mout, const std::type_info &t);

/**
 * @brief Function <code>operator<< </code> reads the crc32 encoding of the name from the type
 * information from the packed byte stream.  The read crc32 is compared to the crc32 encoding of the
 * name from the type information passed.  If the two are different and exception is thrown.
 *
 * @param min     		a <code>Marshal</code> reference to the marshaller.
 *
 * @param t       		a <code>std::type_info</code> const reference to the type
 *                              information to compare with the what was read from the packed byte
 *                              stream.
 *
 * @return			a <code>Marshal</code> reference to the marhsaller.
 */
template<>
Marshal &operator>>(Marshal &min, const std::type_info &t);

template<>
Marshal &operator<<(Marshal &mout, const signed char &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned char &t);
template<>
Marshal &operator<<(Marshal &mout, const char &t);
template<>
Marshal &operator<<(Marshal &mout, const short &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned short &t);
template<>
Marshal &operator<<(Marshal &mout, const int &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned int &t);
template<>
Marshal &operator<<(Marshal &mout, const long &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned long &t);
template<>
Marshal &operator<<(Marshal &mout, const long long &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned long long &t);
template<>
Marshal &operator<<(Marshal &mout, const float &t);
template<>
Marshal &operator<<(Marshal &mout, const double &t);
template<>
Marshal &operator<<(Marshal &mout, const std::string &s);

template<>
Marshal &operator>>(Marshal &min, signed char &t);
template<>
Marshal &operator>>(Marshal &min, unsigned char &t);
template<>
Marshal &operator>>(Marshal &min, char &t);
template<>
Marshal &operator>>(Marshal &min, short &t);
template<>
Marshal &operator>>(Marshal &min, unsigned short &t);
template<>
Marshal &operator>>(Marshal &min, int &t);
template<>
Marshal &operator>>(Marshal &min, unsigned int &t);
template<>
Marshal &operator>>(Marshal &min, long &t);
template<>
Marshal &operator>>(Marshal &min, unsigned long &t);
template<>
Marshal &operator>>(Marshal &min, long long &t);
template<>
Marshal &operator>>(Marshal &min, unsigned long long &t);
template<>
Marshal &operator>>(Marshal &min, float &t);
template<>
Marshal &operator>>(Marshal &min, double &t);
template<>
Marshal &operator>>(Marshal &min, std::string &s);


template <class T>
Marshal &operator<<(Marshal &mout, const std::vector<T> &v)  {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_VECTOR)
    mout << typeid(v);

  size_t size = v.size();
  mout << size;
  for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it)
    mout << (*it);

  return mout;
}

template <class T>
Marshal &operator>>(Marshal &min, std::vector<T> &v)  {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_VECTOR)
    min >> typeid(v);
   
  size_t size = 0;
  min >> size;
  v.reserve(size);
  for (size_t i = 0; i < size; ++i) {
    T t;
    min >> t;
    v.push_back(t);
  }

  return min;
}

template <class T>
Marshal &operator<<(Marshal &mout, const std::list<T> &l)  {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_LIST)
    mout << typeid(l);

  size_t size = l.size();
  mout << size;
  for (typename std::list<T>::const_iterator it = l.begin(); it != l.end(); ++it)
    mout << (*it);

  return mout;
}

template <class T>
Marshal &operator>>(Marshal &min, std::list<T> &l)  {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_LIST)
    min >> typeid(l);

  size_t size;
  min >> size;
  for (size_t i = 0; i < size; ++i) {
    T t;
    min >> t;
    l.push_back(t);
  }

  return min;
}

template <class T>
Marshal &write(Marshal &mout, const T &t) {
  mout.write((const char *) &t, sizeof(T));

  return mout;
}

template <typename T>
Marshal &read(Marshal &min, T &t) {
  t = T();
  
  min.read((char *) &t, sizeof(T));
  return min;
}

} // namespace stk

#endif // STK_UTIL_UTIL_MARSHAL_HPP
