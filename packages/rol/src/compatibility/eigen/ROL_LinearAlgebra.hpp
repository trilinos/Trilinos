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



#pragma once
#ifndef ROL_LINEARALGEBRA_HPP
#define ROL_LINEARALGEBRA_HPP

#include <Eigen/Dense>
#include <ostream>

#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"

/** \file  ROL_LinearAlgebra.hpp
    \brief Provides basic capabilities in solving dense
           linear systems and eigenvalue problems using
           Eigen to provide the implementation */


namespace ROL {

namespace LA {

template<typename Real>
using Vector = Eigen::Matrix<Real,1,Eigen::Dynamic>;

template<typename Real>
using Matrix = Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic>;

template<typename> class LinearSolver;

template<typename Real>
class LinearProblem {
private:

  Ptr<LA::Vector<Real>> b_, x_;
  Ptr<LA::Matrix<Real>> A_;
  
public:

  friend class LinearSolver<Real>;

  LinearProblem();
  LinearProblem( const Ptr<LA::Matrix<Real>>& A,
                 const Ptr<LA::Vector<Real>>& x,	
                 const Ptr<LA::Vector<Real>>& b );

  static Ptr<LinearProble> create( const Ptr<LA::Matrix<Real>>& A,
                                   const Ptr<LA::Vector<Real>>& x,	
                                   const Ptr<LA::Vector<Real>>& b );

  void setMatrix( const Ptr<LA::Matrix<Real>>& A );

  void setVectors( const Ptr<LA::Vector<Real>>& x,
                   const Ptr<LA::Vector<Real>>& b );
};




template<typename Real>
class LinearSolver {
private:
  
 Ptr<LA::LinearProblem<Real>> problem_;
 ParameterList                options_;

public:

  void setOptions( ParameterList& opts );
  void setProblem( const Ptr<LA::PLinearroblem<Real>>& problem );
  void solve( std::ostream& outStream );

};

template<typename> class EigenvalueSolver;

template<typename Real>
class EigenvalueProblem {
private:
 
  Ptr<LA::Vector<Real>> d_;         // Vector of eigenvalues
  Ptr<LA::Matrix<Real>> A_,Vl_,Vr_; // Left and right eigenvectors

public:
 
  friend class EigenvalueSolver<Real>;

  EigenvalueProblem();

  EigenvalueProblem( const Ptr<LA::Matrix<Real>>& A,
                     const Ptr<LA::Matrix<Real>>& V,	
                     const Ptr<LA::Vector<Real>>& d );

  EigenvalueProblem( const Ptr<LA::Matrix<Real>>& A,
		     const Ptr<LA::Matrix<Real>>& Vl,	
                     const Ptr<LA::Matrix<Real>>& Vr,	
                     const Ptr<LA::Vector<Real>>& d );

  

  static Ptr<EigenvaluerProble> create( const Ptr<LA::Matrix<Real>>& A,
                                        const Ptr<LA::Matrix<Real>>& V,	
                                        const Ptr<LA::Vector<Real>>& d );

  static Ptr<EigenvaluerProble> create( const Ptr<LA::Matrix<Real>>& A,
                                        const Ptr<LA::Matrix<Real>>& Vl,	
                                        const Ptr<LA::Matrix<Real>>& Vr,	
                                        const Ptr<LA::Vector<Real>>& d );
  
};


template<typename Real>
class EigenvalueSolver {
private:

  Ptr<LA::LinearProblem<Real>> problem_;
  ParameterList                options_;

public:
  void setOptions( ParameterList& opts );
  void setProblem( const Ptr<LA::EigenvalueProblem<Real>>& problem );
  void solve( std::ostream& outStream );
};


} // namespace LA

} // namespace ROL

#endif // ROL_LINEARAPGEBRA_HPP
