// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_VALIDATEFUNCTIONDEF_HPP
#define ROL_VALIDATEFUNCTIONDEF_HPP

namespace ROL {

namespace details { 

using namespace std;

template<typename Real>
ValidateFunction<Real>::ValidateFunction( const int order,
                                          const int numSteps,
                                          const int width,
                                          const int precision,
                                          const bool printToStream,
                                          ostream& os ) : 
  order_(order), numSteps_(numSteps), width_(width), precision_(precision), 
  printToStream_(printToStream), steps_(numSteps_),
  dashline_(string(precision+6,'-')),  os_(os), 
  workspace_(makePtr<VectorWorkspace<Real>>()),
  fd_( order_, workspace_ ) {

  Real fact = 1.0;
  for( auto& e: steps_ ) {
    e = fact;
    fact *= 0.1;
  }
}


template<typename Real>
vector<vector<Real>> 
ValidateFunction<Real>::derivative_check( f_scalar_t<Real> f_value, 
                                          f_vector_t<Real> f_derivative, 
                                          f_update_t<Real> f_update,
                                          const V& g, 
                                          const V& v,
                                          const V& x,
                                          const string& label ) const {
  int numVals = 4;
  vector<Real> tmp(numVals);
  vector<vector<Real>> vCheck(numSteps_,tmp);
  
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(os_);  

  // Evaluate reference scalar at input
  f_update(x);

  // Evaluate vector quantity
  auto r = workspace_->clone(g);
  
  // Compute the derivative in the given direction
  f_derivative( *r, x );
  //Real dr = r->dot(v.dual());
  Real dr = r->apply(v);
  
  for (int i=0; i<numSteps_; i++) {

    vCheck[i][0] = steps_[i];
    vCheck[i][1] = dr;
    vCheck[i][2] = fd_( f_value, f_update, v, x, steps_[i] );
    vCheck[i][3] = abs(vCheck[i][2] - vCheck[i][1]);

    if (printToStream_) {
      if (i==0) {
        os_ << right
             << setw(width_) << "Step size"
             << setw(width_) << label
             << setw(width_) << "FD approx"
             << setw(width_) << "abs error"
             << "\n"
             << setw(width_) << dashline_ 
             << setw(width_) << dashline_
             << setw(width_) << dashline_
             << setw(width_) << dashline_
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

  os_ << endl;

  // Reset format state of outStream.
  os_.copyfmt(oldFormatState);
  return vCheck;
}
                                        
template<typename Real>
vector<vector<Real>> 
ValidateFunction<Real>::derivative_check( f_vector_t<Real> f_value, 
                                          f_dderiv_t<Real> f_derivative, 
                                          f_update_t<Real> f_update,
                                          const V& c, 
                                          const V& v,
                                          const V& x,
                                          const string& label ) const {
  int numVals = 4;
  vector<Real> tmp(numVals);
  vector<vector<Real>> dCheck(numSteps_,tmp);
  
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(os_);  
  
  auto dc = workspace_->clone(c);
  auto jv = workspace_->clone(c);

  // Evaluate reference vector at input
  f_update(x);

  // Compute directional derivative
  f_derivative( *jv, v, x );

  Real norm_jv = jv->norm();  

  for (int i=0; i<numSteps_; i++) {

    fd_( f_value, f_update, *dc, v, x, steps_[i] );

    dCheck[i][0] = steps_[i];
    dCheck[i][1] = norm_jv;
    dCheck[i][2] = dc->norm();

    dc->axpy(-1.0, *jv);
    dCheck[i][3] = dc->norm();

    if (printToStream_) {
      if (i==0) {
        os_ << right
            << setw(width_) << "Step size"
            << setw(width_) << label
            << setw(width_) << "FD approx"
            << setw(width_) << "abs error"
            << "\n"
            << setw(width_) << dashline_ 
            << setw(width_) << dashline_ 
            << setw(width_) << dashline_ 
            << setw(width_) << dashline_ 
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

  os_ << endl;

  // Reset format state of outStream.
  os_.copyfmt(oldFormatState);
  return dCheck;
}


template<typename Real>
vector<Real> 
ValidateFunction<Real>::symmetry_check( f_dderiv_t<Real> A, 
                                        f_update_t<Real> A_update,
                                        const V& u, 
                                        const V& v, 
                                        const V& x,
                                        const string& name,
                                        const string& symbol ) const {
  auto Au = workspace_->clone(u.dual());
  auto Av = workspace_->clone(v.dual());

  A_update(x);
  A( *Au, u, x );
  A( *Av, v, x );

  //Real vAu = v.dot(Au->dual());
  Real vAu = v.apply(*Au);
  //Real uAv = u.dot(Av->dual());
  Real uAv = u.apply(*Av);

  vector<Real> symCheck(3,0);
  symCheck[0] = vAu;
  symCheck[1] = uAv;
  symCheck[2] = abs(vAu-uAv);

  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(os_);  

  string label1 = "<v," + symbol + "(x)u>";
  string label2 = "<u," + symbol + "(x)v>";

  auto width1 = max(width_,static_cast<int>(label1.length())+3);
  auto width2 = max(width_,static_cast<int>(label2.length())+3);

if (printToStream_) {
    os_ << "\nTest symmetry of " << name << "\n";
    os_ << right
        << setw(width1) << label1
        << setw(width2) << label2
        << setw(width_) << "abs error"
        << "\n";
    os_ << scientific << setprecision(precision_) << right
        << setw(width1) << symCheck[0]
        << setw(width2) << symCheck[1]
        << setw(width_) << symCheck[2]
        << "\n";
  }

  os_ << endl;

  os_.copyfmt(oldFormatState);               
  return symCheck;
}

template<typename Real>
vector<Real> 
ValidateFunction<Real>::adjoint_consistency_check( f_dderiv_t<Real> A,
                                                   f_dderiv_t<Real> A_adj,
                                                   f_update_t<Real> A_update,
                                                   const V& u,  
                                                   const V& v,
                                                   const V& x, 
                                                   const string& name,
                                                   const string& symbol ) const {
  auto Au = workspace_->clone(v.dual());
  auto Av = workspace_->clone(u.dual());

  A_update(x);
  A( *Au, u, x );
  A_adj( *Av, v, x );

  //Real vAu = v.dot(Au->dual());
  Real vAu = v.apply(*Au);
  //Real uAv = u.dot(Av->dual());
  Real uAv = u.apply(*Av);
 
  vector<Real> adjCheck(3,0);
  adjCheck[0] = vAu;
  adjCheck[1] = uAv;
  adjCheck[2] = abs(vAu-uAv);

  string label1 = "<v,"   + symbol + "(x)u>";
  string label2 = "<adj(" + symbol + "(x))v,u>";

  auto width1 = max(width_,static_cast<int>(label1.length())+3);
  auto width2 = max(width_,static_cast<int>(label2.length())+3);

  if (printToStream_) {
      os_ << "\nTest consistency of " << name << " and its adjoint\n";
      os_ << right
          << setw(width1) << label1
          << setw(width2) << label2
          << setw(width_) << "abs error"
          << "\n";
      os_ << scientific << setprecision(precision_) << right
          << setw(width1) << adjCheck[0]
          << setw(width2) << adjCheck[1]
          << setw(width_) << adjCheck[2]
          << "\n";
    }

  os_ << endl;

  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(os_);  

  return adjCheck;
}

template<typename Real>
vector<Real> 
ValidateFunction<Real>::inverse_check( f_dderiv_t<Real> A,
                                       f_dderiv_t<Real> A_inv,
                                       f_update_t<Real> A_update,
                                       const V& v,
                                       const V& x, 
                                       const string& name,
                                       const string& symbol ) const {
  auto Av   = workspace_->clone(v.dual());
  auto AiAv = workspace_->clone(v);

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

  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(os_);  

  string label1 = "||[" + symbol + "]v||";
  string label2 = "||[inv(" + symbol + ")" + symbol + "]v||";
  string label3 = "||v-[inv(" + symbol + ")"  + symbol + "]v||";

  auto width1 = max(width_,static_cast<int>(label1.length())+3);
  auto width2 = max(width_,static_cast<int>(label2.length())+3);
  auto width3 = max(width_,static_cast<int>(label3.length())+3);

if (printToStream_) {
    os_ << "\nTest inverse identity of " << name << "\n";
    os_ << right
        << setw(width_) << "||v||"
        << setw(width1) << label1
        << setw(width2) << label2
        << setw(width3) << label3
        << "\n";
    os_ << scientific << setprecision(precision_) << right
        << setw(width_) << invCheck[0]
        << setw(width1) << invCheck[1]
        << setw(width2) << invCheck[2]
        << setw(width3) << invCheck[3]
        << "\n";
  }

  os_ << endl;

  os_.copyfmt(oldFormatState);               
  return invCheck;
}

template<typename Real>
vector<Real>
ValidateFunction<Real>::solve_check( f_solve_t<Real>  solve,
                                     f_vector_t<Real> value,
                                     f_update_t<Real> update,
                                     const V& c,
                                     const V& x,
                                     const string& name ) const {
  auto res1 = workspace_->clone(c);
  auto res2 = workspace_->clone(c);
  auto sol  = workspace_->clone(x);
  sol->set(x);

  solve(*res1,*sol);
  update(*sol);
  value(*res2,*sol);

  Real sol_norm  = sol->norm();
  Real res1_norm = res1->norm();
  Real res2_norm = res2->norm();

  vector<Real> solveCheck(3,0);
  solveCheck[0] = sol_norm;
  solveCheck[1] = res1_norm;
  solveCheck[2] = res2_norm;

  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(os_);  

  string label1 = "Solution Norm";
  string label2 = "Solver Residual";
  string label3 = "True Residual";

  auto width1 = max(width_,static_cast<int>(label1.length())+3);
  auto width2 = max(width_,static_cast<int>(label2.length())+3);
  auto width3 = max(width_,static_cast<int>(label3.length())+3);

if (printToStream_) {
    os_ << endl << "Test solve of " << name << endl;
    os_ << right
        << setw(width1) << label1
        << setw(width2) << label2
        << setw(width3) << label3
        << endl;
    os_ << scientific << setprecision(precision_) << right
        << setw(width1) << solveCheck[0]
        << setw(width2) << solveCheck[1]
        << setw(width3) << solveCheck[2]
        << endl;
  }
  os_ << endl;

  os_.copyfmt(oldFormatState);               
  return solveCheck;
}
       

template<typename Real>
ostream& ValidateFunction<Real>::getStream() const { return os_; }

} // namespace details

} // namespace ROL

#endif // ROL_VALIDATEFUNCTIONDEF_HPP

