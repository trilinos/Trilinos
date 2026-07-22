// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test for ROL::TeuchosObjective and
           ROL::TeuchosConstraint
  
     Solves the optimization problem

     \f[ \min_x f(x) = \frac{1}{2} x^\top A x-x^\top b \f]
     
     For simplicity, we take A to be positive definte
 
     Subject to the equality constraint

     \f[ c(x) = Cx-d = 0 \f]
*/
    
#include "ROL_Teuchos_Objective.hpp"
#include "ROL_Teuchos_Constraint.hpp"
#include "ROL_OptimizationSolver.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_SerialDenseSolver.hpp"

#include <random>
#include <iostream>

template<class Ordinal, class Real>
class QuadraticTestObjective : public ROL::TeuchosObjective<Ordinal,Real> {

  template<class T> using RCP = Teuchos::RCP<T>;
  using Vector = Teuchos::SerialDenseVector<Ordinal,Real>;
  using Matrix = Teuchos::SerialDenseMatrix<Ordinal,Real>;
  using Solver = Teuchos::SerialDenseSolver<Ordinal,Real>;

private:

  const RCP<const Matrix> A_;
  const RCP<const Vector> b_;
  Vector Ax_;
  RCP<Vector> scratch_;
  Real   bx_;
  Solver solver_;

  static constexpr Real zero{0.0};
  static constexpr Real one{1.0};
  static constexpr Real half{0.5};

  void applyA( Vector& Av, const Vector& v ) {
    Av.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS, one, *A_, v, zero); 
  }

public:

  QuadraticTestObjective( const RCP<const Matrix>& A, 
                          const RCP<const Vector>& b ) :
    A_{A}, b_{b}, Ax_{*b}, scratch_{Teuchos::rcp( new Vector{*b} )} {
    RCP<Matrix> Af = Teuchos::rcp( new Matrix{*A} );
    solver_.setMatrix(Af);
  }

  void update( const Vector& x, bool flag = true, int iter=-1 ) {
    applyA(Ax_,x);
    bx_ = b_->dot(x);
  }

  Real value( const Vector& x, Real& tol ) {
    return half*x.dot(Ax_)-bx_;
  }
 
  void gradient( Vector& g, const Vector& x, Real& tol ) {
    g.assign(Ax_); 
    g -= *b_;
  }

  void hessVec( Vector& hv, const Vector& v, const Vector& x, Real& tol ) {
    applyA(hv,v);
  }

  void invHessVec( Vector& hv, const Vector& v, const Vector& x, Real& tol ) {
    scratch_->assign(v); 
    auto hvp = Teuchos::rcpFromRef(hv);
    solver_.setVectors(hvp,scratch_);
    solver_.solve();
  }
};


template<class Ordinal, class Real>
class LinearTestConstraint : public ROL::TeuchosConstraint<Ordinal,Real> {

  template<class T> using RCP = Teuchos::RCP<T>;
  using Vector = Teuchos::SerialDenseVector<Ordinal,Real>;
  using Matrix = Teuchos::SerialDenseMatrix<Ordinal,Real>;

private:

  const RCP<const Matrix> C_; 
  const RCP<const Vector> d_;
  Vector c_;

  static constexpr Real one{1.0};
  static constexpr Real zero{0.0};

  void applyC( Vector& Cv, const Vector& v ) {
    Cv.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, *C_, v, zero); 
  }

public:

  LinearTestConstraint( const RCP<const Matrix>& C, 
                        const RCP<const Vector>& d ) :
    C_{C}, d_{d}, c_{ C->numRows() }  {}

  void update( const Vector& x, bool flag = true, int iter = -1 ) {
    applyC(c_,x); c_ -= *d_; 
  }
  
  void value( Vector& c, const Vector& x, Real& tol ) {
    c.assign(c_);
  }

  void applyJacobian( Vector& jv, const Vector& v, const Vector& x, Real& tol ) {
    applyC(jv,v);
  }

  void applyAdjointJacobian( Vector& ajv, const Vector& v, const Vector& x, Real& tol ) {
    ajv.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, one, *C_, v, zero);
  }
  
  void applyAdjointHessian( Vector& ahuv, const Vector& u, 
                            const Vector& v, const Vector& x, Real& tol ) {
    ahuv.putScalar(Real{0});
  }
};




using OrdinalT = int;
using RealT = double;

int main( int argc, char *argv[] ) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  using Teuchos::RCP;
  using Teuchos::rcp;  
  using Vector = Teuchos::SerialDenseVector<OrdinalT,RealT>;
  using Matrix = Teuchos::SerialDenseMatrix<OrdinalT,RealT>;

  // This little trick lets us print to std::cout only if a 
  // (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

//  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  std::default_random_engine gen;
  std::normal_distribution<RealT> dist(0.0,1.0);

  // *** Test body.
  try {

    // Dimension of optimization space
    OrdinalT Nopt = 10;

    // Dimension of constraint space
    OrdinalT Ncon = 3;

    auto A = rcp( new Matrix{Nopt,Nopt} );
    auto b = rcp( new Vector{Nopt} );
    auto C = rcp( new Matrix{Ncon,Nopt} );
    auto d = rcp( new Vector{Ncon} );

    // Create a symmetric random positive definte matrix A
    // random rectangular matrix C, and rectangular vectors b and d

    for( OrdinalT i=0; i<Nopt; ++i ) {
      RealT sum{0};
      for( OrdinalT j=i+1; j<Nopt; ++j ) {
        (*A)(i,j) = dist(gen);
        (*A)(j,i) = (*A)(i,j);
        sum += std::abs((*A)(i,j));
      }
      (*A)(i,i) = 5*sum;
      (*b)(i) = dist(gen);
      for( OrdinalT k=0; k<Ncon; ++k ) {
        (*C)(k,i) = dist(gen);
        if(i==0) { 
          (*d)(k) = dist(gen);
        } 
      }
    }

    *outStream << "\nA = " << printMat(*A);
    *outStream << "\nb = " << printMat(*b);
    *outStream << "\nC = " << printMat(*C);
    *outStream << "\nd = " << printMat(*d);
    

    auto x = rcp( new Vector{Nopt,1} );
    auto l = rcp( new Vector{Ncon,1} );

    // Solution vector
    auto sol = rcp( new ROL::TeuchosVector<OrdinalT,RealT>{x} );
    
    // Equality multiplier
    auto mul = rcp( new ROL::TeuchosVector<OrdinalT,RealT>{l} );

    // Objective
    auto obj = rcp( new QuadraticTestObjective<OrdinalT,RealT>{A,b} );

    // Constraint
    auto con = rcp( new LinearTestConstraint<OrdinalT,RealT>{C,d} );
     
    ROL::OptimizationProblem<RealT> problem(obj,sol,con,mul);

    problem.check(*outStream);

    Teuchos::ParameterList emptyList;

    ROL::OptimizationSolver<RealT> solver(problem,emptyList);
    solver.solve(*outStream);


  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";


  return 0;
}
