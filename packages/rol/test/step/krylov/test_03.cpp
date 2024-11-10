// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_03.cpp
    \brief Test MINRES solver
*/

#include "ROL_StdVector.hpp"
#include "ROL_GMRES.hpp"
#include "ROL_KrylovFactory.hpp"
#include "ROL_RandomVector.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include<iomanip>

// Identity operator for preconditioner
template<class Real> 
class Identity : public ROL::LinearOperator<Real> {
  typedef ROL::Vector<Real> V;
public:
  void apply( V& Hv, const V& v, Real &tol ) const { 
    Hv.set(v); 
  }
}; // class Identity


// Apply a tridiagonal Toeplitz matrix to a ROL::StdVector to test Krylov solvers
template<class Real>
class TridiagonalToeplitzOperator : public ROL::LinearOperator<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;

  typedef typename vector::size_type uint; 

private:

  Real a_; // subdiagonal
  Real b_; // diagonal 
  Real c_; // superdiagonal

  ROL::LAPACK<int,Real> lapack_;

public:

  TridiagonalToeplitzOperator( Real &a, Real &b, Real &c ) : a_(a), b_(b), c_(c) {}

  // Tridiagonal multiplication
  void apply( V &Hv, const V &v, Real &tol ) const {
 
    SV &Hvs = dynamic_cast<SV&>(Hv);
    ROL::Ptr<vector> Hvp = Hvs.getVector();
 
    const SV &vs = dynamic_cast<const SV&>(v);
    ROL::Ptr<const vector> vp = vs.getVector();

    uint n = vp->size();

    (*Hvp)[0] = b_*(*vp)[0] + c_*(*vp)[1];

    for(uint k=1; k<n-1; ++k) {
      (*Hvp)[k]  = a_*(*vp)[k-1] + b_*(*vp)[k] + c_*(*vp)[k+1];  
    } 
  
    (*Hvp)[n-1] = a_*(*vp)[n-2] + b_*(*vp)[n-1];

  }

  // Tridiagonal solve - compare against GMRES
  void applyInverse( V &Hv, const V &v, Real &tol ) const {

      
 
    SV &Hvs = dynamic_cast<SV&>(Hv);
    ROL::Ptr<vector> Hvp = Hvs.getVector();
 
    const SV &vs = dynamic_cast<const SV&>(v);
    ROL::Ptr<const vector> vp = vs.getVector();

    uint n = vp->size();

    const char TRANS = 'N';
    const int NRHS = 1;

    vector dl(n-1,a_);
    vector d(n,b_);
    vector du(n-1,c_);
    vector du2(n-2,0.0);
 
    std::vector<int> ipiv(n);
    int info;

    Hv.set(v); // LAPACK will modify this in place    

    // Do Tridiagonal LU factorization
    lapack_.GTTRF(n,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&info);

    // Solve the system with the LU factors
    lapack_.GTTRS(TRANS,n,NRHS,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&(*Hvp)[0],n,&info);

  }

}; // class TridiagonalToeplitzOperator



typedef double RealT;

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>            vector;
  typedef ROL::StdVector<RealT>         SV; 

  typedef typename vector::size_type    uint;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {

    ROL::ParameterList parlist;
    ROL::ParameterList &gList = parlist.sublist("General");
    ROL::ParameterList &kList = gList.sublist("Krylov");
    
    kList.set("Type","MINRES");
    kList.set("Iteration Limit",20);
    kList.set("Absolute Tolerance",1.e-8);
    kList.set("Relative Tolerance",1.e-6);
    kList.set("Use Initial Guess",false);

    uint dim = 10;

    auto xp = ROL::makePtr<vector>(dim,0.0);
    auto yp = ROL::makePtr<vector>(dim,0.0);
    auto zp = ROL::makePtr<vector>(dim,0.0);
    auto bp = ROL::makePtr<vector>(dim,0.0);

    SV x(xp); // Exact solution
    SV y(yp); // Solution using direct solve
    SV z(zp); // Solution using GMRES

    SV b(bp); // Right-hand-side    
    
    RealT left = -1.0;
    RealT right = 1.0;    

    ROL::RandomizeVector(x,left,right);
   
    RealT sub   =  1.0;
    RealT diag  = -2.0;
    RealT super =  1.0;

    TridiagonalToeplitzOperator<RealT> T(sub,diag,super);
    Identity<RealT> I;

    RealT tol = 0.0;
      
    T.apply(b,x,tol);

    T.applyInverse(y,b,tol);

    auto krylov = ROL::KrylovFactory<RealT>( parlist );

    int iter;
    int flag;
    z.zero();
    krylov->run(z,T,b,I,iter,flag);

    *outStream << std::setw(10) << "Exact"  
               << std::setw(10) << "LAPACK"
               << std::setw(10) << "MINRES" << std::endl;
    *outStream << "---------------------------------" << std::endl;

    for(uint k=0;k<dim;++k) {
      *outStream << std::setw(10) << (*xp)[k] << " " 
                 << std::setw(10) << (*yp)[k] << " "
                 << std::setw(10) << (*zp)[k] << " " << std::endl;
    }  
 
    *outStream << "MINRES performed " << iter << " iterations." << std::endl;

    z.axpy(-1.0,x);

    if( z.norm() > std::sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      ++errorFlag;
    }

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
