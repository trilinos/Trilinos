#ifndef SAMBA_SAMBA_UTILITY_MACROS_HPP
#define SAMBA_SAMBA_UTILITY_MACROS_HPP

#include <samba/utility/is_primitive.hpp>
#include <samba/utility/primitive_arithmetic_operators.hpp>
#include <samba/utility/primitive_comparison.hpp>

#include <boost/type_traits/is_pod.hpp>
#include <boost/type_traits/has_nothrow_copy.hpp>
#include <boost/type_traits/has_nothrow_constructor.hpp>

#ifdef SAMBA_ENABLE_PARALLEL
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/mpi/datatype.hpp>
#endif

//****************************************************************************
//setup primitives type traits for broken compilers
//****************************************************************************
#ifdef SAMBA_ENABLE_PARALLEL

#define SAMBA_IS_PRIMITIVE(type)                   \
                                                   \
  namespace samba {                                \
  template<> struct is_primitive<type>             \
    : public boost::mpl::true_ {};                 \
  } /*namespace samba*/                            \
  namespace boost {                                \
  template <> struct is_pod<type>                  \
    : public mpl::true_ {};                        \
  } /*namespace boost*/                            \
  namespace boost { namespace serialization {      \
  template<> struct is_bitwise_serializable<type>  \
    : public mpl::true_ {};                        \
  }} /*namespace boost::serialization*/            \
  namespace boost { namespace mpi {                \
  template<> struct is_mpi_datatype<type>          \
    : public mpl::true_ {};                        \
  }} /*namespace boost::mpi*/

#else

#define SAMBA_IS_PRIMITIVE(type)                   \
                                                   \
  namespace samba {                                \
  template<> struct is_primitive<type>             \
    : public boost::mpl::true_ {};                 \
  } /*namespace samba*/                            \
  namespace boost {                                \
  template <> struct is_pod<type>                  \
    : public mpl::true_ {};                        \
  } /*namespace boost*/

#endif

//****************************************************************************
//setup type traits for broken compilers
//****************************************************************************
#define SAMBA_HAS_NO_THROW_COPY(type)      \
  namespace boost {                        \
  template<> struct has_nothrow_copy<type> \
    : public mpl::true_ {};                \
  }

#define SAMBA_HAS_NO_THROW_CONSTRUCTOR(type)      \
  namespace boost {                               \
  template<> struct has_nothrow_constructor<type> \
    : public mpl::true_ {};                       \
  }

#endif //SAMBA_SAMBA_UTILITY_MACROS_HPP
