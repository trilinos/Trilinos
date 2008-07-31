/*
  This is a simple test to check the portability of the typedef declarations for ScalarT.  If portable, this allows us to add calculation types with minimal additions to the code.
*/

#include <iostream>

class Traits {
public:

  struct ResidualType {};
  struct JacobianType {};

  template <typename CalcT> class CalcTMap {};
};


template <> struct Traits::CalcTMap<Traits::ResidualType> { typedef double
ScalarT; };
template <> struct Traits::CalcTMap<Traits::JacobianType> { typedef int
ScalarT; };

template <typename CalcT, typename Traits>
class Density {
public:

  typedef typename Traits::template CalcTMap<CalcT>::ScalarT ScalarT;

  Density(ScalarT x);

  ScalarT x_;

};

template <typename CalcT, typename Traits>
Density<CalcT, Traits>::Density(ScalarT x) :
  x_(x)
{
}

int main() {
  Density<Traits::ResidualType, Traits> density(1.0);
  std::cout << "Hello World!" << std::endl;
  return 0;
}
