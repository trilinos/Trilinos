#ifndef STK_MIDDLE_MESH_IS_POINT
#define STK_MIDDLE_MESH_IS_POINT

namespace stk {
namespace middle_mesh {
namespace utils {

// SFINAE test for if a thing has the Point interface
template <typename T>
class IsPoint
{
  private:
    typedef char one;
    struct two { char x[2]; };

    template <typename C>
    static one test( decltype(&C::getX));

    template <typename C>
    static two test(...);

  public:
    enum
    {
      Value = sizeof(test<T>(0)) == sizeof(char)
    };
};


}
}
}

#endif