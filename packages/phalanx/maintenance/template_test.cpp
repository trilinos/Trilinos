// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
