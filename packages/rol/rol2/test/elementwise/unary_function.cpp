#include "ROL2.hpp"
#include <algorithm>

template<typename Real,typename Function>
void index_function( Function f, std::vector<Real>& x ) {
  for( int i=0; i<x.size(); ++i ) x[i] = f(i);
}

template<typename Real>
void apply_unary( std::vector<Real>& x, 
                  const ROL2::Elementwise::UnaryFunction<Real>& uf ) {
  for( auto& e : x ) e = uf.apply(e);
}

template<typename Real> 
std::ostream& operator << ( std::ostream& os, const std::vector<Real>& x ) {
  for( const auto& e : x ) os << e <<  std::endl;
  return os;
}

//auto square = []( auto x ) { return x*x; };

//template<typename Real>
//struct Square : public ROL2::Elementwise::UnaryFromFunction<Real,decltype(square)>{
//  Square() : ROL2::Elementwise::UnaryFromFunction<Real,decltype(square)>::UnaryFromFunction(square) {};
//};



int main( int argc, char* argv[] ) {

  using RealT = double;

  auto os_ptr = ROL2::makeStreamPtr(std::cout, argc);
  auto& os = *os_ptr;
  int errorFlag = 0;


  try {

    int n = 10;
    std::vector<RealT> x(n);
    index_function( [](auto i){ return std::pow(-1,i)/std::sqrt(1+i*i); }, x );
    ROL2::Elementwise::Fill<RealT> fill(1.0);
    ROL2::Elementwise::AbsoluteValue<RealT> abs;
    ROL2::Elementwise::Heaviside<RealT> heaviside;

    os << x << std::endl;
    apply_unary(x,heaviside);
    os << x << std::endl;

  } catch( std::exception& e ) {
    errorFlag = -1000;
    os << e.what() << std::endl;
  };

  if( errorFlag ) std::cout << "End Result: TEST FAILED!\n";
  else            std::cout << "End Result: TEST PASSED!\n";

  return 0;
}
