/*
  This is a simple test to check the portability of the typedef declarations for ScalarT.  If portable, this allows us to add calculation types with minimal additions to the code.
*/

#include <iostream>

class Traits {
public:

  struct Residual { typedef double ScalarT; };
  struct Jacobian { typedef int ScalarT; };

  //template <typename CalcT> class CalcTMap {};
};


//template <> struct Traits::CalcTMap<Traits::ResidualType> { typedef double
//ScalarT; };
//template <> struct Traits::CalcTMap<Traits::JacobianType> { typedef int
//ScalarT; };

template <typename CalcT, typename Traits>
class Density {
public:

  //typedef typename Traits::template CalcTMap<CalcT>::ScalarT ScalarT;
  typedef typename CalcT::ScalarT ScalarT;

  Density(ScalarT x);

  void evaluate();

  ScalarT x_;

};

template <typename CalcT, typename Traits>
Density<CalcT, Traits>::Density(ScalarT x) :
  x_(x)
{
}

template <typename CalcT, typename Traits>
void Density<CalcT, Traits>::evaluate()
{
  std::cout << "Start evaluate" << std::endl;
  ScalarT tmp_val;
  std::cout << "Finish evaluate" << std::endl;
}

int main() {
  Density<Traits::Residual, Traits> r_density(1.0);
  Density<Traits::Jacobian, Traits> j_density(1);

  r_density.evaluate();

  std::cout << "Hello World!" << std::endl;
  return 0;
}
