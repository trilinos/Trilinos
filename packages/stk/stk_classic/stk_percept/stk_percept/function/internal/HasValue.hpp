#ifndef stk_percept_function_HasValue_hpp
#define stk_percept_function_HasValue_hpp

namespace stk_classic
{
  namespace percept
  {
    template<typename ValueType>
    class HasValue
    {
    public:
      virtual ValueType& getValue() = 0;
      virtual void setValue(ValueType& ) = 0;
      virtual ~HasValue() {}
    };

    template<typename ValueType>
    class HasConstValue
    {
    public:
      virtual ValueType& getValue() = 0;
      virtual ~HasConstValue() {}
    };


  }
}
#endif
