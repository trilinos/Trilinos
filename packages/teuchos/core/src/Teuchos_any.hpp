/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#ifndef TEUCHOS_ANY_HPP
#define TEUCHOS_ANY_HPP

/*! \file Teuchos_any.hpp
   \brief Modified boost::any class for holding a templated value
*/

#include <utility>
#include <type_traits>
#include <exception>

#include "Teuchos_Assert.hpp"
#include "Teuchos_TypeNameTraits.hpp"

//
// This file was taken from the boost library which contained the
// following notice:
//
// *************************************************************
//
// what:  variant type boost::any
// who:   contributed by Kevlin Henney,
//        with features contributed and bugs found by
//        Ed Brey, Mark Rodgers, Peter Dimov, and James Curran
// when:  July 2001
// where: tested with BCC 5.5, MSVC 6.0, and g++ 2.95
//
// Copyright Kevlin Henney, 2000, 2001, 2002. All rights reserved.
//
// Permission to use, copy, modify, and distribute this software for any
// purpose is hereby granted without fee, provided that this copyright and
// permissions notice appear in all copies and derivatives.
//
// This software is provided "as is" without express or implied warranty.
//
// *************************************************************
//
// RAB modified the file for use in Teuchos.  I changed the nature of
// the any_cast<> to be easier to use.
//

namespace Teuchos {

template<class T>
struct is_comparable
{
    template<class X>
    static auto test(int) -> decltype(std::declval<X>() == std::declval<X>(),
                                      void(), std::true_type());
    template<class X>
    static auto test(...) -> std::false_type;
    using type = decltype(test<T>(0));
};

template<class T>
struct is_printable
{
    template<class X>
    static auto test(int) -> decltype(std::declval<std::ostream&>() << std::declval<X>(),
                                      void(), std::true_type());
    template<class X>
    static auto test(...) -> std::false_type;
    using type = decltype(test<T>(0));
};

template <class T, class ok = typename is_comparable<T>::type>
struct compare;

template <class T>
struct compare<T, std::false_type> {
  bool operator()(T const&, T const&) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Trying to compare type " << typeid(T).name() << " which is not comparable");
#ifndef __CUDACC__
    return false;
#endif
  }
};

template <class T>
struct compare<T, std::true_type> {
  bool operator()(T const& a, T const& b) const {
    return a == b;
  }
};

template <class T, class ok = typename is_printable<T>::type>
struct print;

template <class T>
struct print<T, std::false_type> {
  std::ostream& operator()(std::ostream& s, T const&) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Trying to print type " << typeid(T).name() << " which is not printable");
#ifndef __CUDACC__
    return s;
#endif
  }
};

template <class T>
struct print<T, std::true_type> {
  std::ostream& operator()(std::ostream& a, T const& b) const {
    return a << b;
  }
};

/** \brief Modified boost::any class, which is a container for a templated
 * value.
 */
class TEUCHOSCORE_LIB_DLL_EXPORT any
{
public:
  //! Empty constructor
  any()
    : content(0)
    {}

  //! Templated constructor
  template<typename ValueType>
  explicit any(const ValueType & value)
    : content(new holder<ValueType>(value))
    {}

  //! Copy constructor
  any(const any & other)
    : content(other.content ? other.content->clone() : 0)
    {}

  //! Destructor
  ~any()
    {
      delete content;
    }

  //! Method for swapping the contents of two any classes
  any & swap(any & rhs)
    {
      std::swap(content, rhs.content);
      return *this;
    }

  //! Copy the value <tt>rhs</tt>
  template<typename ValueType>
  any & operator=(const ValueType & rhs)
    {
      any(rhs).swap(*this);
      return *this;
    }

  //! Copy the value held in <tt>rhs</tt>
  any & operator=(const any & rhs)
    {
      any(rhs).swap(*this);
      return *this;
    }

  //! Return true if nothing is being stored
  bool empty() const
    {
      return !content;
    }

  //! Return the type of value being stored
  const std::type_info & type() const
    {
      return content ? content->type() : typeid(void);
    }

  //! Return the name of the type
  std::string typeName() const
    {
      return content ? content->typeName() : "NONE";
    }

  /*! \brief Return if two any objects are the same or not.
   *  \warning This function with throw an exception if
   *           operator== can't be applied to the held type!
   */
  bool same( const any &other ) const
    {
      if( this->empty() && other.empty() )
        return true;
      else if( this->empty() && !other.empty() )
        return false;
      else if( !this->empty() && other.empty() )
        return false;
      // !this->empty() && !other.empty()
      return content->same(*other.content);
    }

  /*! \brief Print this value to the output stream <tt>os</tt>
   *  \warning This function with throw an exception if
   *           the held type can't be printed via operator<< !
   */
  void print(std::ostream& os) const
    {
      if (content) content->print(os);
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  /** @name Private??? types */
  //@{

  /** \brief . */
  class placeholder
  {
  public:
    /** \brief . */
    virtual ~placeholder() {}
    /** \brief . */
    virtual const std::type_info & type() const = 0;
    /** \brief . */
    virtual std::string typeName() const = 0;
    /** \brief . */
    virtual placeholder * clone() const = 0;
    /** \brief . */
    virtual bool same( const placeholder &other ) const = 0;
    /** \brief . */
    virtual void print(std::ostream & os) const = 0;
  };

  /** \brief . */
  template<typename ValueType>
  class holder : public placeholder
  {
  public:
    /** \brief . */
    holder(const ValueType & value)
      : held(value)
      {}
    /** \brief . */
    const std::type_info & type() const
      { return typeid(ValueType); }
    /** \brief . */
    std::string typeName() const
      { return TypeNameTraits<ValueType>::name(); }
    /** \brief . */
    placeholder * clone() const
      { return new holder(held); }
    /** \brief . */
    bool same( const placeholder &other ) const
      {
        if( type() != other.type() ) {
          return false;
        }
        // type() == other.type()
        const ValueType
          &other_held = dynamic_cast<const holder<ValueType>&>(other).held;
        return ::Teuchos::compare<ValueType>{}(held, other_held);
      }
    /** \brief . */
    void print(std::ostream & os) const
      { ::Teuchos::print<ValueType>{}(os, held); }
    /** \brief . */
    ValueType held;
  };

  //@}

public:
  // Danger: This is made public to allow any_cast to be non-friend
  placeholder* access_content()
    { return content; }
  const placeholder* access_content() const
    { return content; }
#endif

private:

  // /////////////////////////
  // Private data members

  placeholder * content;

};

/*! \relates any
    \brief Thrown if any_cast is attempted between two incompatable types.
*/
class bad_any_cast : public std::runtime_error
{
public:
  bad_any_cast( const std::string msg ) : std::runtime_error(msg) {}
};

/*! \relates any
    \brief Used to extract the templated value held in Teuchos::any to a given value type.

    \note <ul>   <li> If the templated value type and templated type are not the same then a
    bad_any_cast is thrown.
    <li> If the dynamic cast fails, then a Teuchos::bad_any_cast std::exception is thrown.
    </ul>
*/
template<typename ValueType>
ValueType& any_cast(any &operand)
{
  const std::string ValueTypeName = TypeNameTraits<ValueType>::name();
  TEUCHOS_TEST_FOR_EXCEPTION(
    operand.type() != typeid(ValueType), bad_any_cast,
    "any_cast<"<<ValueTypeName<<">(operand): Error, cast to type "
    << "any::holder<"<<ValueTypeName<<"> failed since the actual underlying type is \'"
    << typeName(*operand.access_content()) << "!"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    !operand.access_content(), bad_any_cast
    ,"any_cast<"<<ValueTypeName<<">(operand): Error, cast to type "
    << "any::holder<"<<ValueTypeName<<"> failed because the content is NULL"
    );
  any::holder<ValueType>
    *dyn_cast_content = dynamic_cast<any::holder<ValueType>*>(operand.access_content());
  TEUCHOS_TEST_FOR_EXCEPTION(
    !dyn_cast_content, std::logic_error
    ,"any_cast<"<<ValueTypeName <<">(operand): Error, cast to type "
    << "any::holder<"<<ValueTypeName<<"> failed but should not have and the actual underlying type is \'"
    << typeName(*operand.access_content()) << "!"
    << "  The problem might be related to incompatible RTTI systems in static and shared libraries!"
    );
  return dyn_cast_content->held;
}

/*! \relates any
    \brief Used to extract the const templated value held in Teuchos::any to a given
  const value type.

    \note <ul>   <li> If the templated value type and templated type are not the same then a
    bad_any_cast is thrown.
    <li> If the dynamic cast fails, then a logic_error is thrown.
    </ul>
*/
template<typename ValueType>
const ValueType& any_cast(const any &operand)
{
  return any_cast<ValueType>(const_cast<any&>(operand));
}

/*! \relates any
    \brief Keep the convenient behavior of Teuchos::any_cast w.r.t. references, but don't confuse it
      with the behavior for C++17 std::any_cast.

    \note In C++17, one must use std::any_cast<T&> to get a reference. This function will ensure
    that uses of Teuchos::any_ref_cast<T> are consciously replaced with std::any_cast<T&> when C++17 is used.
*/
template<typename ValueType>
ValueType& any_ref_cast(any &operand)
{
  return Teuchos::any_cast<ValueType>(operand);
}

/*! \relates any
    \brief Converts the value in <tt>any</tt> to a std::string.
    \warning This function with throw an exception if
             the held type can't be printed via operator<< !
*/
inline std::string toString(const any &rhs)
{
  std::ostringstream oss;
  rhs.print(oss);
  return oss.str();
}

/*! \relates any
    \brief Returns true if two any objects have the same value.
    \warning This function with throw an exception if
             operator== can't be applied to the held type!
*/
inline bool operator==( const any &a, const any &b )
{
  return a.same(b);
}

/*! \relates any
    \brief Returns true if two any objects <b>do not</b> have the same value.
    \warning This function with throw an exception if
             operator== can't be applied to the held type!
*/
inline bool operator!=( const any &a, const any &b )
{
  return !a.same(b);
}

/*! \relates any
    \brief Writes "any" input <tt>rhs</tt> to the output stream <tt>os</tt>.
    \warning This function with throw an exception if
             the held type can't be printed via operator<< !
*/
inline std::ostream & operator<<(std::ostream & os, const any &rhs)
{
  rhs.print(os);
  return os;
}

/*! \relates any
    \brief Special swap for other code to find via Argument Dependent Lookup
*/
inline void swap(Teuchos::any& a, Teuchos::any& b) {
  a.swap(b);
}

/*! \relates any
    \brief Default constructs a new T value and returns a reference to it.
*/
template <typename T>
T & make_any_ref(any &rhs)
{
  any(T()).swap(rhs);
  return any_cast<T>(rhs);
}

} // namespace Teuchos

#endif // TEUCHOS_ANY_HPP
