// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
ROL2_TypeU_LineSearch_Decl.hpp// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
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
#ifndef ROL2_TYPEU_LINESEARCH_DECL_H
#define ROL2_TYPEU_LINESEARCH_DECL_H

/** \class ROL::LineSearch_U
    \brief Provides interface for and implements line searches.
*/

namespace ROL2 {
namespace TypeU {

template<class Real>
class LineSearch {
public:

  enum class Type : std::int16_t {
     IterationScaling = 0,
     PathBasedTargetLevel,
     Backtracking,
     Bisection,
     GoldenSection,
     CubicInterp,
     Brents,
     UserDefined,
     Last  
  };

  enum class CurvatureCond : std::int16_t {
    Wolfe = 0,
    StrongWolfe,
    GeneralizedWolfe,
    ApproximateWolfe,
    Goldstein,
    Null,
    Last  
  };

  static EnumMap<Type>          type_dict;
  static EnumMap<CurvatureCond> curvature_dict;

  LineSearch() = default;
  virtual ~LineSearch() = default;

  LineSearch( ParameterList& parlist );

  virtual void initialize( const Vector<Real>& x, 
                           const Vector<Real>& g);

  virtual void run(       Real&            alpha, 
                          Real&            fval, 
                          int&             ls_neval, 
                          int&             ls_ngrad,
                    const Real&            gs, 
                    const Vector<Real>&    s, 
                    const Vector<Real>&    x, 
                          Objective<Real>& obj ) = 0;

  // use this function to modify alpha and fval if the maximum number of iterations
  // are reached
  void setMaxitUpdate(       Real& alpha, 
                             Real& fnew, 
                       const Real& fold );


protected:
  virtual bool status( const Type             type, 
                             int&             ls_neval, 
                             int&             ls_ngrad, 
                       const Real             alpha, 
                       const Real             fold, 
                       const Real             sgold, 
                       const Real             fnew, 
                       const Vector<Real>&    x, 
                       const Vector<Real>&    s, 
                             Objective<Real>& obj );

  virtual Real getInitialAlpha(       int&             ls_neval,
                                      int&             ls_ngrad,
                                      Real             fval,
                                      Real             gs,
                                const Vector<Real>&    x,
                                const Vector<Real>&    s,
                                      Objective<Real>& obj);

  void setNextInitialAlpha( Real alpha );
  bool useLocalMinimizer();
  bool takeNoStep();

private:

  Real dirDeriv( const Vector<Real>&    x,      // current iterate
                 const Vector<Real>&    s,      // current trial step
                       Real             alpha,  // current step length
                       Real             fnew,   // f(x+alpha*s)
                       Objective<Real>& obj);

  Ptr<Vector<Real>> xtst_; 

  Real c1_, c2_, c3_, eps_;
  Real fmin_;         // smallest fval encountered
  Real alphaMin_;     // Alpha that yields the smallest fval encountered
  Real alpha0_;
  Real alpha0bnd_;    // Lower bound for initial alpha...if below, set initial alpha to one
 
  CurvatureCond econd_;
  DescentType   edesc_;

  int maxit_;

  bool useralpha_;
  bool usePrevAlpha_; // Use the previous step's accepted alpha as an initial guess
  bool acceptMin_;    // Use smallest fval if sufficient decrease not satisfied
  bool itcond_;       // true if maximum function evaluations reached
  bool FDdirDeriv_;

};

template<class Real>
inline std::string enumToString( LineSearch<Real>::Type e ) { 
  return LineSearch<Real>::type_dict[e]; 
}

template<class Real>
inline std::string enumToString( CurvatureCond e ) { 
  return LineSearch<Real>::curvature_dict[e]; 
}

//template<class Real>
//inline std::string stringToEnum( std::string s, const LineSearch<Real>& ) {
//  return LineSearch<Real>::type_dict[s];
//}
//
//template<class Real>
//inline CurvatureCond stringToEnum( std::string, LineSearch<Real>::CurvatureCond ) {
//  return curvature_dict[s];
//}


template<class Real>
EnumMap<LineSearch<Real>::Type>
LineSearch<Real>::type_dict = { "Iteration Scaling",
                                "Path-Based Target Level",
                                "Backtracking",
                                "Bisection",           
                                "Golden Section",
                                "Cubic Interpolation",
                                "Brent's",    
                                "User Defined" };
template<class Real>
EnumMap<LineSearch<Real>::DescentType>
LineSearch<Real>::curvature_dict = { "Wolfe Conditions",
                                     "Strong Wolfe Conditions",
                                     "Generalized Wolfe Conditions",
                                     "Approximate Wolfe Conditions",
                                     "Goldstein Conditions",
                                     "Null Curvature Condition" };

 




} // namespace ROL2

#endif // ROL2_TYPEU_LINESEARCH_DECL_H
