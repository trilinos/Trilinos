// @HEADER
// ************************************************************************
// @HEADER


/*! \file  test_02.cpp
    \brief Test Rosenbrock.
*/

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Algorithm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

template<class Real>
class Objective_Rosenbrock : public ROL::Objective<Real> {
private:
  Real alpha_;

public:

  Objective_Rosenbrock(Real alpha = 100.0) : alpha_(alpha) {}

  Real value( const ROL::Vector<Real> &x ) {
    ROL::StdVector<Real> & ex =
      Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();
    int n = xp->size();

    Real val = 0;
    for( int i=0; i<n/2; i++ ) {
      val += alpha_ * pow(pow((*xp)[2*i],2) - (*xp)[2*i+1], 2);
      val += pow((*xp)[2*i] - 1.0, 2);
    }

    return val;
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x ) {

    Teuchos::RCP<const std::vector<Real> > xp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > gp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());

    int n = xp->size();

    for( int i=0; i<n/2; i++ ) {
      (*gp)[2*i]   =  4.0*alpha_*(pow((*xp)[2*i],2) - (*xp)[2*i+1])*(*xp)[2*i] + 2.0*((*xp)[2*i]-1.0); 
      (*gp)[2*i+1] = -2.0*alpha_*(pow((*xp)[2*i],2) - (*xp)[2*i+1]);
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {

    Teuchos::RCP<const std::vector<Real> > xp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<std::vector<Real> > hvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());

    int n = xp->size();

    for( int i=0; i<n/2; i++ ) { 
      Real h11 = 4.0*alpha_*(3.0*pow((*xp)[2*i],2)-(*xp)[2*i+1]) + 2.0;
      Real h12 = -4.0*alpha_*(*xp)[2*i];
      Real h22 = 2.0*alpha_;

      (*hvp)[2*i]   = h11*(*vp)[2*i] + h12*(*vp)[2*i+1]; 
      (*hvp)[2*i+1] = h12*(*vp)[2*i] + h22*(*vp)[2*i+1];
    }
  }
};

template<class Real>
class Objective_SumSquares : public ROL::Objective<Real> {
public:
  Real value( const ROL::Vector<Real> &x ) {
    ROL::StdVector<Real> & ex =
      Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();
    int n = xp->size();

    Real val = 0;
    for (int i=0; i<n; i++) {
      val += pow((*xp)[i], 2);
    }

    return val;
 }

};


int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Test body.

  try {

    int dim = 100;
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> x(x_rcp);
    ROL::StdVector<RealT> y(y_rcp);

    // set x,y
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i] = 1.0;
      (*y_rcp)[i] = i;
    }

    Objective_Rosenbrock<RealT> obj;
    //Objective_SumSquares<RealT> obj;
  
    *outStream << "\nObjective value: " << obj.value(x) << "\n";

    *outStream << "\nDirectional derivative: " << obj.dirDeriv(x,y,1e-2*Teuchos::ScalarTraits<RealT>::eps()) << "\n";

    obj.gradient(y,x);

    *outStream << "\nNorm of finite-difference gradient: " << y.norm() << "\n";

    x.zero();
    y.zero();

    for (int i=0; i<dim/2; i++) {
      (*x_rcp)[2*i]   = -1.2;
      (*x_rcp)[2*i+1] =  1.0;
    }

    ROL::LineSearchStepType LSStype = ROL::LineSearchStep_NewtonKrylov;
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_Secant;
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_Gradient;

    ROL::SecantType Stype = ROL::Secant_lDFP;
    //ROL::SecantType Stype = ROL::Secant_lBFGS;
    //ROL::SecantType Stype = ROL::Secant_BarzilaiBorwein;

    //ROL::LineSearchType LStype = ROL::LineSearchType_Backtracking;
    //ROL::LineSearchType LStype = ROL::LineSearchType_SimpleBacktracking;
    ROL::LineSearchType LStype = ROL::LineSearchType_Brents;
    //ROL::LineSearchType LStype = ROL::LineSearchType_Bisection;
    //ROL::LineSearchType LStype = ROL::LineSearchType_GoldenSection;

    ROL::LineSearchCondition LScond = ROL::LineSearchCondition_Wolfe;
    //ROL::LineSearchCondition LScond = ROL::LineSearchCondition_StrongWolfe;
    //ROL::LineSearchCondition LScond = ROL::LineSearchCondition_Goldstein;
 
    int maxit    = 20;
    RealT rho    = 0.5;
    RealT c1     = 1.e-4;
    RealT c2     = 0.9;
    RealT tol    = 1.e-8;

    int L        = 10;
    int BBtype   = 1;

    RealT CGtol1 = 1.e-4;
    RealT CGtol2 = 1.e-2;
    int maxitCG  = 100;

    ROL::LineSearchStep<RealT> step(LStype,LScond,LSStype,maxit,c1,c2,tol,rho,Stype,L,BBtype,CGtol1,CGtol2,maxitCG);
    ROL::StatusTest<RealT> status(1.e-6,1.e-12,100000);    

    ROL::DefaultAlgorithm<RealT> algo(step,status);

    //Teuchos::RCP<ROL::Algorithm<RealT> > algo;
    //ROL::DefaultAlgorithmFactory<RealT> algoFactory.getAlgo(algo, parlist);
    
    algo.run(x, obj);

    obj.gradient(y,x);

    *outStream << "\nNorm of finite-difference gradient: " << y.norm() << "\n";


  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

