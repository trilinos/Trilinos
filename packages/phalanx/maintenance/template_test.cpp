// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

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
