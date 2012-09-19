#ifndef stk_encr_Dimensions_hpp
#define stk_encr_Dimensions_hpp

#include <vector>
#include <iostream>

namespace stk
{
  namespace percept
  {
    class Dimensions : public std::vector<int>
    {
    public:
      typedef std::vector<int> base_type;
      //Dimensions() : base_type(1, 0) {}
      Dimensions() : base_type( 0) {}
      Dimensions(int i0) : base_type(1) { (*this)[0] = i0; }
      Dimensions(int i0, int i1) : base_type(2) 
      { 
        (*this)[0] = i0; 
        (*this)[1] = i1; 
      }
      Dimensions(int i0, int i1, int i2) : base_type(3) 
      { 
        (*this)[0] = i0; 
        (*this)[1] = i1; 
        (*this)[2] = i2; 
      }
    };

#ifndef SWIG
    std::ostream &operator<<(std::ostream& out, const Dimensions& dim);
#endif

  }
}

#endif
