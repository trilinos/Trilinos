// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_util_util_random_access_iterator_wrapper_hpp
#define stk_util_util_random_access_iterator_wrapper_hpp

#include <boost/iterator/iterator_adaptor.hpp>

#ifndef BOOST_NO_SFINAE
# include <boost/type_traits/is_convertible.hpp>
# include <boost/utility/enable_if.hpp>
#endif

namespace stk_util {


template <class Value>
class random_access_iterator_wrapper
  : public boost::iterator_adaptor<
        random_access_iterator_wrapper<Value>   // Derived
      , Value*                                  // Base
      , boost::use_default                      // Value
      , boost::random_access_traversal_tag      // CategoryOrTraversal
    >
{
  private:

    typedef boost::iterator_adaptor<
                                     random_access_iterator_wrapper<Value>,
                                     Value*,
                                     boost::use_default,
                                     boost::random_access_traversal_tag
                                   > base_type;

    struct enabler {};  // used to enable coversion constructor (if SFINAE)

  public:
    random_access_iterator_wrapper()
      : base_type(0) {}

    explicit random_access_iterator_wrapper(Value* p)
      : base_type(p) {}

    Value& operator[] (ptrdiff_t index)
    { return *(*this + index); }

    template <class OtherValue>
      random_access_iterator_wrapper(
          random_access_iterator_wrapper<OtherValue> const& other
# ifndef BOOST_NO_SFINAE
          , typename boost::enable_if<
          boost::is_convertible<OtherValue*,Value*>
          , enabler
          >::type = enabler()
# endif
          )
      : base_type(other.base()) {}
};

}//namespace stk_util

#endif

