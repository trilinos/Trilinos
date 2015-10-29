#ifndef PHALANX_ANY_HPP
#define PHALANX_ANY_HPP

#include <algorithm>
#include <typeinfo>
#include <type_traits>
#include "Teuchos_Assert.hpp"
#include "Sacado_mpl_if.hpp"

#  define BOOST_AUX_ANY_TYPE_ID_NAME
#include <cstring>

namespace PHX
{
  template <bool B, class T = void>
  struct disable_if_c {
    typedef T type;
  };
  
  template <class T>
  struct disable_if_c<true, T> {};
  
  template <class Cond, class T = void> 
  struct disable_if : public disable_if_c<Cond::value, T> {};

    class any
    {
    public: // structors

        any() noexcept
          : content(0)
        {
        }

        template<typename ValueType>
        any(const ValueType & value)
          : content(new holder<typename std::decay<const ValueType>::type>(value))
        {
        }

        any(const any & other)
          : content(other.content ? other.content->clone() : 0)
        {
        }

        // Move constructor
        any(any&& other) noexcept
          : content(other.content)
        {
            other.content = 0;
        }

        // Perfect forwarding of ValueType
        template<typename ValueType>
        any(ValueType&& value
            , typename PHX::disable_if<std::is_same<any&, ValueType> >::type* = 0 // disable if value has type `any&`
            , typename PHX::disable_if<std::is_const<ValueType> >::type* = 0) // disable if value has type `const ValueType&&`
          : content(new holder< typename std::decay<ValueType>::type >(static_cast<ValueType&&>(value)))
        {
        }

        ~any() noexcept
        {
            delete content;
        }

    public: // modifiers

        any & swap(any & rhs) noexcept
        {
            std::swap(content, rhs.content);
            return *this;
        }


        any & operator=(const any& rhs)
        {
            any(rhs).swap(*this);
            return *this;
        }

        // move assignement
        any & operator=(any&& rhs) noexcept
        {
            rhs.swap(*this);
            any().swap(rhs);
            return *this;
        }

        // Perfect forwarding of ValueType
        template <class ValueType>
        any & operator=(ValueType&& rhs)
        {
            any(static_cast<ValueType&&>(rhs)).swap(*this);
            return *this;
        }

    public: // queries

        bool empty() const noexcept
        {
            return !content;
        }

        void clear() noexcept
        {
            any().swap(*this);
        }

        const std::type_info & type() const noexcept
        {
            return content ? content->type() : typeid(void);
        }

#ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
    private: // types
#else
    public: // types (public so any_cast can be non-friend)
#endif

        class placeholder
        {
        public: // structors

            virtual ~placeholder()
            {
            }

        public: // queries

            virtual const std::type_info & type() const noexcept = 0;

            virtual placeholder * clone() const = 0;

        };

        template<typename ValueType>
        class holder : public placeholder
        {
        public: // structors

            holder(const ValueType & value)
              : held(value)
            {
            }

            holder(ValueType&& value)
              : held(static_cast< ValueType&& >(value))
            {
            }

        public: // queries

            virtual const std::type_info & type() const noexcept
            {
                return typeid(ValueType);
            }

            virtual placeholder * clone() const
            {
                return new holder(held);
            }

        public: // representation

            ValueType held;

        private: // intentionally left unimplemented
            holder & operator=(const holder &);
        };

#ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS

    private: // representation

        template<typename ValueType>
        friend ValueType * any_cast(any *) noexcept;

        template<typename ValueType>
        friend ValueType * unsafe_any_cast(any *) noexcept;

#else

    public: // representation (public so any_cast can be non-friend)

#endif

        placeholder * content;

    };
 
    inline void swap(any & lhs, any & rhs) noexcept
    {
        lhs.swap(rhs);
    }

    class bad_any_cast : public std::bad_cast
    {
    public:
        virtual const char * what() const noexcept
        {
            return "PHX::any::bad_any_cast: "
                   "failed conversion using PHX::any_cast";
        }
    };

    template<typename ValueType>
    ValueType * any_cast(any * operand) noexcept
    {
        return operand && 
            operand->type() == typeid(ValueType)
            ? &static_cast<any::holder<ValueType> *>(operand->content)->held
            : 0;
    }

    template<typename ValueType>
    inline const ValueType * any_cast(const any * operand) noexcept
    {
        return any_cast<ValueType>(const_cast<any *>(operand));
    }

    template<typename ValueType>
    ValueType any_cast(any & operand)
    {
      typedef typename std::remove_reference<ValueType>::type nonref;

#ifdef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
        // If 'nonref' is still reference type, it means the user has not
        // specialized 'remove_reference'.

        // Please use BOOST_BROKEN_COMPILER_TYPE_TRAITS_SPECIALIZATION macro
        // to generate specialization of remove_reference for your class
        // See type traits library documentation for details
      static_assert(!std::is_reference<nonref>::value);
#endif

        nonref * result = any_cast<nonref>(&operand);
	TEUCHOS_TEST_FOR_EXCEPTION(!result,
				   std::runtime_error,
				   "ERROR: failed to cast from any object to requested object type!");
	
        // Attempt to avoid construction of a temporary object in cases when 
        // `ValueType` is not a reference. Example:
        // `static_cast<std::string>(*result);` 
        // which is equal to `std::string(*result);`
        typedef typename Sacado::mpl::mpl_if<
	  std::is_reference<ValueType>,
            ValueType,
	  typename std::add_lvalue_reference<ValueType>::type
        >::type ref_type;
	

        return static_cast<ref_type>(*result);
    }

    template<typename ValueType>
    inline ValueType any_cast(const any & operand)
    {
      typedef typename std::remove_reference<ValueType>::type nonref;

#ifdef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
        // The comment in the above version of 'any_cast' explains when this
        // assert is fired and what to do.
        static_assert(!is_reference<nonref>::value);
#endif

        return any_cast<const nonref &>(const_cast<any &>(operand));
    }

    template<typename ValueType>
    inline ValueType&& any_cast(any&& operand)
    {
        BOOST_STATIC_ASSERT_MSG(
            std::is_rvalue_reference<ValueType&&>::value 
            || std::is_const< typename std::remove_reference<ValueType>::type >::value,
            "PHX::any_cast shall not be used for getting nonconst references to temporary objects" 
        );
        return any_cast<ValueType&&>(operand);
    }
}

// Copyright Kevlin Henney, 2000, 2001, 2002. All rights reserved.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#endif
