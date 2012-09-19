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

