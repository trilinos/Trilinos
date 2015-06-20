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


#ifndef ROL_SANDIARULES2_HPP
#define ROL_SANDIARULES2_HPP

# include <fstream>
# include <string>
# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>
# include <ctime>

# include "ROL_SandiaRules.hpp"

namespace ROL {

class SandiaRules2 : public SandiaRules {
public: 

  SandiaRules2(void) {};

  void ccn_points ( int n, int dim, double x[] );
  void ccn_weights ( int n, int dim, double w[] );

  void clenshaw_curtis_points ( int n, int dim, double x[] );
  void clenshaw_curtis_weights ( int n, int dim, double w[] );

  void fejer2_points ( int n, int dim, double x[] );
  void fejer2_weights ( int n, int dim, double w[] );

  void gen_hermite_points ( int n, int dim, double x[] );
  void gen_hermite_weights ( int n, int dim, double w[] );

  void gen_laguerre_points ( int n, int dim, double x[] );
  void gen_laguerre_weights ( int n, int dim, double w[] );

  void hcc_points ( int n, int dim, double x[] );
  void hcc_weights ( int n, int dim, double w[] );

  void hce_points ( int n, int dim, double x[] );
  void hce_weights ( int n, int dim, double w[] );

  void hermite_genz_keister_points ( int n, int dim, double x[] );
  void hermite_genz_keister_weights ( int n, int dim, double w[] );

  void hermite_points ( int n, int dim, double x[] );
  void hermite_weights ( int n, int dim, double w[] );

  void jacobi_points ( int n, int dim, double x[] );
  void jacobi_weights ( int n, int dim, double w[] );

  void laguerre_points ( int n, int dim, double x[] );
  void laguerre_weights ( int n, int dim, double w[] );

  void legendre_points ( int n, int dim, double x[] );
  void legendre_weights ( int n, int dim, double w[] );

  void ncc_points ( int n, int dim, double x[] );
  void ncc_weights ( int n, int dim, double w[] );

  void nco_points ( int n, int dim, double x[] );
  void nco_weights ( int n, int dim, double w[] );

  double parameter ( int dim, int offset );

  void patterson_points ( int n, int dim, double x[] );
  void patterson_weights ( int n, int dim, double w[] );
}; // class SandiaRules2

} // namespace ROL

#include "ROL_SandiaRules2Def.hpp"
#endif
