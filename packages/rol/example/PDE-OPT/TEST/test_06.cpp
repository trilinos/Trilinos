// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file test_06.cpp
    \brief Unit test for Automatic Differentiation evalutators
*/

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


#include "../TOOLS/template_tools.hpp"
#include "ROL_Ptr.hpp"

// Example of ScalarFunction 
// f(x,y) = <x,x> + 2*<y,y> - 2 y[0]*(x[0]+x[1]) + 3*x[0]*(y[1]-y[0])
//
// coeff = [1,2,-2,3]
template<class Param, template<class> class Array>
struct ExampleQuadratic {

//  using Real = TemplateTools::ElementType<Param>;



  template<class Real, class ScalarX, class ScalarY>
  static AD::ResultType<Real,ScalarX,ScalarY> 
  eval( const Param &p,
        const Array<ScalarX> &x, 
        const Array<ScalarY> &y ) {

    using Index = TemplateTools::index_t<Array<ScalarX>>;    
    using TemplateTools::get;

    ScalarX xdotx(0); 
    ScalarY ydoty(0);

    auto nx = TemplateTools::size(x);
    auto ny = TemplateTools::size(y);     

    for( Index i=0; i<nx; ++i ) {
      xdotx += get(x, i);
    }
    for( Index j=0; j<ny; ++j ) {
      ydoty = get(y,j);
    }

    return get(p,0)*xdotx + 
           get(p,1)*ydoty + 
           get(p,2)*get(y,0)*(get(x,0)+get(x,1)) +
           get(p,3)*get(x,0)*(get(y,1)-get(y,0));

  }
}; // Example Quadratic


struct Foo {
  int size(){return 0;} 
  void bar(int i) {}
};

int main(int argc, char *argv[] ) {

  using RealT = double;

    
  ROL::nullstream bhs;

  ROL::Ptr<std::ostream> os;
  if(argc>1)   os = ROL::makePtrFromRef(std::cout);
  else         os = ROL::makePtrFromRef(bhs);
  
  using DFad = Sacado::Fad::DFad<RealT>;

  int errorFlag = 0;

  // Univariate tests
  {
    RealT x = 2.0;

    DFad x_fad(1,0,x);
    DFad f_fad = x_fad*x_fad;

    errorFlag += static_cast<int>(AD::has_dx<RealT>::value);
    errorFlag += static_cast<int>(!AD::has_dx<DFad>::value);

    *os << "Does type RealT = double have a method .dx(int)?... " 
        << AD::has_dx<RealT>::value << std::endl;

    *os << "Does type DFad<RealT> have a method .dx(int)?...... " 
        << AD::has_dx<DFad>::value << std::endl;

    auto xval  = AD::value(x);
    auto xfval = AD::value(x_fad);
    auto ffval = AD::value(f_fad);

    auto xder  = AD::derivative(x,0);
    auto xfder = AD::derivative(x_fad,0);
    auto ffder = AD::derivative(f_fad,0);

    *os << "RealT x(2.0);"             << std::endl;
    *os << "DFad x_fad(1,x); "         << std::endl;
    *os << "DFad f_fad = x_fad*x_fad " << std::endl;  

    *os << "AD::value(x)     = " << xval  << std::endl;
    *os << "AD::value(x_fad) = " << xfval << std::endl;
    *os << "AD::value(f_fad) = " << ffval << std::endl;

    errorFlag += (xval  != 2);
    errorFlag += (xfval != 2);
    errorFlag += (ffval != 4);

    *os << "AD::derivative(x,0)      = "  << xder  << std::endl;
    *os << "AD::derivative(x_fad,0)  = "  << xfder << std::endl;
    *os << "AD::derivative(f_fad,0)  = "  << ffder << std::endl;

    errorFlag += (xder  != 0);
    errorFlag += (xfder != 1);
    errorFlag += (ffder != 4);
  }

  // Multivariate tests
  {
    using size_type = typename std::vector<RealT>::size_type;

    size_type Nx = 5;
    size_type Ny = 3;

    std::vector<RealT> coeff({1.0,2.0,-2.0,3.0});

    std::vector<RealT> x(Nx);
    std::vector<RealT> y(Ny); 

//    std::cout << TemplateTools::size(coeff) << std::endl;
//    std::cout << TemplateTools::max( coeff ) << std::endl;
//    std::cout << TemplateTools::getElement( coeff, 1 ) << std::endl;

  }



  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}
