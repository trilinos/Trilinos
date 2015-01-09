#ifndef stk_percept_Name_hpp
#define stk_percept_Name_hpp

#include <string>

namespace stk_classic
{
  namespace percept
  {

    /// this is to avoid a common bug where the name of the String Function is given instead of function_string,
    ///    ie., first two args are accidentally reversed - this is essentially a model of a "named argument",
    ///    as opposed to a positional one; of course, it is still only a hint (though a strong one) to the user
    
    /// Useful in other places where two strings are passed into a function or constructor
    class Name 
    {
      const std::string m_name;
    public:
      explicit 
      Name(const std::string name) : m_name(name) {}
      const std::string& getName() const { return m_name; }
    };

  }
}
#endif
