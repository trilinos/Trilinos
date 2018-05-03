// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (width_14) Sandia Corporation
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


#pragma once
#ifndef ROL_FINITEDIFFERENCEDEF_HPP
#define ROL_FINITEDIFFERENCEDEF_HPP

#include <iomanip>

namespace ROL {

namespace details {

using namespace std;
using ::ROL::Finite_Difference_Arrays::shifts;
using ::ROL::Finite_Difference_Arrays::weights;

template<typename Real>
FiniteDifference<Real>::FiniteDifference( Teuchos::ParameterList& pl,
                                          ostream& os ) :
  order_(1), numSteps_(ROL_NUM_CHECKDERIV_STEPS), width_(width_), precision_(11),
  printToStream_(true), steps_(numSteps_), os_(cout) {

  auto &fdlist = pl.sublist("Finite Difference Check");
  
  order_         = fdlist.get("Finite Difference Order", 1);
  numSteps_      = fdlist.get("Number of Steps", ROL_NUM_CHECKDERIV_STEPS);
  width_         = fdlist.get("Output Display Field Width", 8);
  precision_     = fdlist.get("Output Display Precision", 16);
  printToStream_ = fdlist.get("Output Display", true);


  TEUCHOS_TEST_FOR_EXCEPTION( order_<1 || order_>4, std::invalid_argument, 
                              "Error: finite difference order must be 1,2,3, or 4" );

  for(int i=0;i<numSteps_;++i) steps_[i] = pow(10,-i);

}


template<typename Real>
FiniteDifference<Real>::FiniteDifference( const int order, 
                                          const int numSteps,
                                          const int width,
                                          const int precision,
                                          const bool printToStream,
                                          ostream& os ) :
  order_(order), numSteps_(numSteps), width_(width), precision_(precision), 
  printToStream_(printToStream), steps_(numSteps_), os_(os) {

  TEUCHOS_TEST_FOR_EXCEPTION( order_<1 || order_>4, std::invalid_argument, 
                              "Error: finite difference order must be 1,2,3, or 4" );

  for(int i=0;i<numSteps_;++i) steps_[i] = pow(10,-i);
} 
 
template<typename Real>

vector<vector<Real>>
FiniteDifference<Real>::scalar_check( f_scalar_t<Real> f_ref,
                                      f_vector_t<Real> f_test,
                                      f_update_t<Real> f_update,
                                      const Vector<Real>& result,
                                      const Vector<Real>& direction,
                                      const Vector<Real>& input,
                                      const string& label ) const {
  int numVals = 4;
  vector<Real> tmp(numVals);
  vector<vector<Real>> vCheck(numSteps_,tmp);
  
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os_);  

  // Evaluate reference scalar at input
  f_update(input);
  Real val = f_ref(input);

  // Evaluate vector quantity
  auto r = workspace_.clone(result);
  
  // Compute the derivative in the given direction
  f_test( *r, input );
  Real dr = r->dot(direction.dual());
  
  auto x = workspace_.clone(input);
  
  for (int i=0; i<numSteps_; i++) {

    Real eta = steps_[i];
    x->set(input);

    vCheck[i][0] = eta;
    vCheck[i][1] = dr;
    vCheck[i][2] = weights[order_-1][0] * val;

    for(int j=0; j<order_; ++j) {
      x->axpy(eta*shifts[order_-1][j], direction);

      // Only evaluate at shifts where the weight is nonzero  
      if( weights[order_-1][j+1] != 0 ) {
        f_update(*x);
        vCheck[i][2] += weights[order_-1][j+1] * f_ref(*x);
      }
    }

    vCheck[i][2] /= eta;

    vCheck[i][3] = abs(vCheck[i][2] - vCheck[i][1]);

    if (printToStream_) {
      if (i==0) {
        os_ << right
             << setw(width_) << "Step size"
             << setw(width_) << label
             << setw(width_) << "FD approx"
             << setw(width_) << "abs error"
             << "\n"
             << setw(width_) << "---------"
             << setw(width_) << "---------"
             << setw(width_) << "---------"
             << setw(width_) << "---------"
             << "\n";
      }
      os_ << scientific << setprecision(precision_) << right
           << setw(width_) << vCheck[i][0]
           << setw(width_) << vCheck[i][1]
           << setw(width_) << vCheck[i][2]
           << setw(width_) << vCheck[i][3]
           << "\n";
    }
  }
  // Reset format state of outStream.
  os_.copyfmt(oldFormatState);
  return vCheck;
}


template<typename Real>
vector<vector<Real>>
FiniteDifference<Real>::vector_check( f_vector_t<Real> f_ref,
                                      f_dderiv_t<Real> f_test,
                                      f_update_t<Real> f_update,
                                      const Vector<Real>& result,
                                      const Vector<Real>& direction,
                                      const Vector<Real>& input,
                                      const string& label ) const {


  int numVals = 4;
  vector<Real> tmp(numVals);
  vector<vector<Real>> dCheck(numSteps_,tmp);
  
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os_);  

  
  auto r = workspace_.clone(result);

  // Evaluate reference vector at input
  f_update(input);
  f_ref(*r, input);

//  auto rdir = workspace_.clone(r);
  auto rdir = result.clone(); 

  // Compute directional derivative
  f_test( *rdir, direction, input );

  Real norm_rdir = rdir->norm();  

  auto rdif = workspace_.clone(r);
  auto rnew = workspace_.clone(r);
  auto x    = workspace_.clone(input);
  
  for (int i=0; i<numSteps_; i++) {

    Real eta = steps_[i];

    x->set(input);

    rdif->set(*r);
    rdif->scale(weights[order_-1][0]);

    for(int j=0; j<order_; ++j) {

      x->axpy(eta*shifts[order_-1][j], direction);

      if( weights[order_-1][j+1] != 0 ) {
        f_update(*x);
        f_ref(*rnew, *x);  
        rdif->axpy(weights[order_-1][j+1],*rnew);
      }
    }

   rdif->scale(1.0/eta);

   dCheck[i][0] = eta;
   dCheck[i][1] = norm_rdir;
   dCheck[i][2] = rdif->norm();

   rdif->axpy(-1.0, *rdir);
   dCheck[i][3] = rdif->norm();

   if (printToStream_) {
      if (i==0) {
        os_ << right
            << setw(width_) << "Step size"
            << setw(width_) << label
            << setw(width_) << "FD approx"
            << setw(width_) << "abs error"
            << "\n"
            << setw(width_) << "---------"
            << setw(width_) << "---------"
            << setw(width_) << "---------"
            << setw(width_) << "---------"
            << "\n";
      }

      os_ << scientific << setprecision(precision_) << right
          << setw(width_) << dCheck[i][0]
          << setw(width_) << dCheck[i][1]
          << setw(width_) << dCheck[i][2]
          << setw(width_) << dCheck[i][3]
          << "\n";
    }
  }
  // Reset format state of outStream.
  os_.copyfmt(oldFormatState);
  return dCheck;
}

template<typename Real>
vector<Real>
FiniteDifference<Real>::symmetry_check( f_dderiv_t<Real> A, 
                                        f_update_t<Real> A_update,
                                        const Vector<Real>& u, 
                                        const Vector<Real>& v, 
                                        const Vector<Real>& x,
                                        const string& name, 
                                        const string& symbol ) const {

  auto Au = workspace_.clone(u.dual());
  auto Av = workspace_.clone(v.dual());

  A_update(x);
  A( *Au, u, x );
  A( *Av, v, x );

  Real vAu = v.dot(Au->dual());
  Real uAv = u.dot(Av->dual());

  vector<Real> symCheck(3,0);
  symCheck[0] = vAu;
  symCheck[1] = uAv;
  symCheck[2] = abs(vAu-uAv);

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os_);  

if (printToStream_) {
    os_ << "\nTest symmetry of " << name << "\n";
    os_ << right
        << setw(width_) << "<v, " << symbol << "(x)u>"
        << setw(width_) << "<u, " << symbol << "(x)v>"
        << setw(width_) << "abs error"
        << "\n";
    os_ << scientific << setprecision(precision_) << right
        << setw(width_) << symCheck[0]
        << setw(width_) << symCheck[1]
        << setw(width_) << symCheck[2]
        << "\n";
  }

  os_.copyfmt(oldFormatState);               
  return symCheck;
}

template<typename Real>
vector<Real>
FiniteDifference<Real>::adjoint_consistency_check( f_dderiv_t<Real> A, 
                                                   f_dderiv_t<Real> A_adj,
                                                   f_update_t<Real> A_update,
                                                   const Vector<Real>& u, 
                                                   const Vector<Real>& v, 
                                                   const Vector<Real>& x,
                                                   const string& name,
                                                   const string& symbol ) const {

  auto Au = workspace_.clone(u.dual());
  auto Av = workspace_.clone(v.dual());

  A_update(x);
  A( *Au, u, x );
  A_adj( *Av, v, x );

  Real vAu = v.dot(Au->dual());
  Real uAv = u.dot(Av->dual());
 
  vector<Real> adjCheck(3,0);
  adjCheck[0] = vAu;
  adjCheck[1] = uAv;
  adjCheck[2] = abs(vAu-uAv);

  if (printToStream_) {
      os_ << "\nTest consistency of " << name << " and its adjoint\n";
      os_ << right
          << setw(width_) << "<v, "  << symbol << "(x)u>"
          << setw(width_) << "<adj(" << symbol << "(x))v,u>"
          << setw(width_) << "abs error"
          << "\n";
      os_ << scientific << setprecision(precision_) << right
          << setw(width_) << adjCheck[0]
          << setw(width_) << adjCheck[1]
          << setw(width_) << adjCheck[2]
          << "\n";
    }

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os_);  


  return adjCheck;
}

template<typename Real>
vector<Real>
FiniteDifference<Real>::inverse_check( f_dderiv_t<Real> A, 
                                       f_dderiv_t<Real> A_inv,
                                       f_update_t<Real> A_update,
                                       const Vector<Real>& v, 
                                       const Vector<Real>& x,
                                       const string& name, 
                                       const string& symbol ) const {

  auto Av   = workspace_.clone(v.dual());
  auto AiAv = workspace_.clone(v);

  A_update(x);
  A( *Av, v, x );
  A_inv( *AiAv, *Av, x );

  Real v_norm    = v.norm();
  Real Av_norm   = Av->norm();
  Real AiAv_norm = AiAv->norm();

  AiAv->axpy(-1.0,v);

  Real err = AiAv->norm(); 

  vector<Real> invCheck(4,0);
  invCheck[0] = v_norm;
  invCheck[1] = Av_norm;
  invCheck[2] = AiAv_norm;
  invCheck[3] = err;

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os_);  

if (printToStream_) {
    os_ << "\nTest inverse identity of " << name << "\n";
    os_ << right
        << setw(width_) << "||v||"
        << setw(width_) << "||" << symbol << "v||"
        << setw(width_) << "||inv(" << symbol << ")" << symbol << "v||"
        << setw(width_) << "||v-inv(" << symbol << ")" << symbol << "v||"
        << "\n";
    os_ << scientific << setprecision(precision_) << right
        << setw(width_) << invCheck[0]
        << setw(width_) << invCheck[1]
        << setw(width_) << invCheck[2]
        << setw(width_) << invCheck[3]
        << "\n";
  }

  os_.copyfmt(oldFormatState);               
  return invCheck;
}







} // namespace details

} // namespace ROL

#endif // ROL_FINITEDIFFERENCEDEF_HPP

