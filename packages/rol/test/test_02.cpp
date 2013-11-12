// @HEADER
// ************************************************************************
// @HEADER


/*! \file  test_02.cpp
    \brief Test Rosenbrock.
*/

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
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

  void invHessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {

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

      (*hvp)[2*i]   = (1.0/(h11*h22-h12*h12))*( h22*(*vp)[2*i] - h12*(*vp)[2*i+1]); 
      (*hvp)[2*i+1] = (1.0/(h11*h22-h12*h12))*(-h12*(*vp)[2*i] + h11*(*vp)[2*i+1]);
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

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x ) {

    Teuchos::RCP<const std::vector<Real> > xp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > gp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());

    int n = xp->size();

    for( int i=0; i<n; i++ ) {
      (*gp)[i] = 2.0*(*xp)[i]; 
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

    for( int i=0; i<n; i++ ) { 
      (*hvp)[i] = 2.0*(*vp)[i]; 
    }
  }

  void invHessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {

    Teuchos::RCP<const std::vector<Real> > xp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<std::vector<Real> > hvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());

    int n = xp->size();

    for( int i=0; i<n; i++ ) { 
      (*hvp)[i] = 0.5*(*vp)[i]; 
    }
  }
};

template<class Real>
class Objective_LeastSquares : public ROL::Objective<Real> {
public:
  Real value( const ROL::Vector<Real> &x ) {
    ROL::StdVector<Real> & ex =
      Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();
    int n = xp->size();

    Real h = 1.0/((Real)n+1.0);

    Real val = 0.0;
    Real res = 0.0;
    for (int i=0; i<n; i++) {
      if ( i == 0 ) {
        res = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i+1]-2.0*(*xp)[i]);
      }  
      else if ( i == n-1 ) {
        res = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]);
      }
      else {
        res = 2.0*h + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]+(*xp)[i+1]);  
      } 
      val += 0.5*res*res; 
    }

    return val;
 }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x ) {

    Teuchos::RCP<const std::vector<Real> > xp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > gp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());

    int n = xp->size();

    std::vector<Real> res(n,0.0);

    Real h = 1.0/((Real)n+1.0);

    for (int i=0; i<n; i++) {
      if ( i == 0 ) {
        res[i] = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i+1]-2.0*(*xp)[i]);
      }  
      else if ( i == n-1 ) {
        res[i] = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]);
      }
      else {
        res[i] = 2.0*h + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]+(*xp)[i+1]);  
      } 
    }

    for (int i=0; i<n; i++) {
      if ( i == 0 ) {
        (*gp)[i] = 1.0/h*(res[i+1]-2.0*res[i]);
      }  
      else if ( i == n-1 ) {
        (*gp)[i] = 1.0/h*(res[i-1]-2.0*res[i]);
      }
      else {
        (*gp)[i] = 1.0/h*(res[i-1]-2.0*res[i]+res[i+1]);  
      } 
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

    std::vector<Real> res(n,0.0);

    Real h = 1.0/((Real)n+1.0);

    for (int i=0; i<n; i++) {
      if ( i == 0 ) {
        res[i] = 1.0/h*((*vp)[i+1]-2.0*(*vp)[i]);
      }  
      else if ( i == n-1 ) {
        res[i] = 1.0/h*((*vp)[i-1]-2.0*(*vp)[i]);
      }
      else {
        res[i] = 1.0/h*((*vp)[i-1]-2.0*(*vp)[i]+(*vp)[i+1]);  
      } 
    }

    for (int i=0; i<n; i++) {
      if ( i == 0 ) {
        (*hvp)[i] = 1.0/h*(res[i+1]-2.0*res[i]);
      }  
      else if ( i == n-1 ) {
        (*hvp)[i] = 1.0/h*(res[i-1]-2.0*res[i]);
      }
      else {
        (*hvp)[i] = 1.0/h*(res[i-1]-2.0*res[i]+res[i+1]);  
      } 
    }
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

    int objFunc = 3;

    int dim = 64;
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > z_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > e_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> x(x_rcp);
    ROL::StdVector<RealT> y(y_rcp);
    ROL::StdVector<RealT> z(z_rcp);
    ROL::StdVector<RealT> e(e_rcp);

    // set x,y
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i] = 1.0;
      (*y_rcp)[i] = i;
      if ( objFunc == 1 ) {
        (*z_rcp)[i] = 1.0;
      }
      else if ( objFunc == 2 ) {
        (*z_rcp)[i] = 0.0;
      }
      else if ( objFunc == 3 ) {
        RealT h  = 1.0/((RealT)dim+1.0);
        RealT pt = (RealT)(i+1)*h;
        (*z_rcp)[i] = pt*(1.0-pt);
      }
    }

    Teuchos::RCP<ROL::Objective<RealT> > obj;
    if ( objFunc == 1 ) {
      *outStream << "\nROSENBROCK'S FUNCTION\n";
      obj = Teuchos::rcp(new Objective_Rosenbrock<RealT>);
    }
    else if ( objFunc == 2 ) {
      *outStream << "\nSUM OF SQUARES OBJECTIVE\n";
      obj = Teuchos::rcp(new Objective_SumSquares<RealT>);
    }
    else if ( objFunc == 3 ) {
      *outStream << "\nLEAST SQUARES OBJECTIVE\n";
      obj = Teuchos::rcp(new Objective_LeastSquares<RealT>);
    }
    obj->checkGradient(x,y,true);

    /* BEGIN SECANT DEFINITION */
    ROL::SecantType Stype = ROL::Secant_lDFP;
    //ROL::SecantType Stype = ROL::Secant_lSR1;
    //ROL::SecantType Stype = ROL::Secant_lBFGS;
    //ROL::SecantType Stype = ROL::Secant_BarzilaiBorwein;

    int L        = 10;
    int BBtype   = 1;
    /* END SECANT DEFINTION */

    /* BEGIN LINE SEARCH STEP DEFINTION */
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_Newton;
    ROL::LineSearchStepType LSStype = ROL::LineSearchStep_NewtonKrylov;
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_NewtonKrylovSecantPreconditioning;
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_Secant;
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_Gradient;

    ROL::LineSearchType LStype = ROL::LineSearchType_Backtracking;
    //ROL::LineSearchType LStype = ROL::LineSearchType_SimpleBacktracking;
    //ROL::LineSearchType LStype = ROL::LineSearchType_Brents;
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

    RealT CGtol1 = 1.e-4;
    RealT CGtol2 = 1.e-2;
    int maxitCG  = 200;

    ROL::LineSearchStep<RealT> LS_step(LStype,LScond,LSStype,maxit,c1,c2,tol,rho,
                                       Stype,L,BBtype,CGtol1,CGtol2,maxitCG);
    /* END LINE SEARCH STEP DEFINITION */

    /* BEGIN TRUST REGION STEP DEFINTION */
    //ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_Newton;
    ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_NewtonKrylov;
    //ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_NewtonKrylovSecantPreconditioning;
    //ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_Secant;
    //ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_Gradient;

    //ROL::TrustRegionType TRtype = ROL::TrustRegionType_CauchyPoint;
    ROL::TrustRegionType TRtype = ROL::TrustRegionType_TruncatedCG;
    //ROL::TrustRegionType TRtype = ROL::TrustRegionType_DoubleDogleg;  
    //ROL::TrustRegionType TRtype = ROL::TrustRegionType_Dogleg;

    maxit        = 200;
    RealT tol1   = 1.e-4;
    RealT tol2   = 1.e-2;
    RealT del    = 100.0;
    RealT delmin = 1.e-8;
    RealT delmax = 5000.0;
    RealT eta0   = 0.05;
    RealT eta1   = 0.05;
    RealT eta2   = 0.9;
    RealT gamma0 = 0.0625;
    RealT gamma1 = 0.25;
    RealT gamma2 = 2.50;
    RealT TRsafe = 1.0;

    ROL::TrustRegionStep<RealT> TR_step(TRtype,TRStype,maxit,tol1,tol2,del,delmin,delmax,
                                        eta0,eta1,eta2,gamma0,gamma1,gamma2,TRsafe,
                                        Stype,L,BBtype);
    /* END TRUST REGION STEP DEFINITION */ 

    ROL::StatusTest<RealT> status(1.e-8,1.e-16,100);    

    for (int i=0; i<dim/2; i++) {
      (*x_rcp)[2*i]   = -1.2;
      (*x_rcp)[2*i+1] =  1.0;
    }
    ROL::DefaultAlgorithm<RealT> LS_algo(LS_step,status);
    LS_algo.run(x, *obj);
    e.set(x);
    e.axpy(-1.0,z);
    *outStream << "\nNorm of Error: " << e.norm() << "\n";

    for (int i=0; i<dim/2; i++) {
      (*x_rcp)[2*i]   = -1.2;
      (*x_rcp)[2*i+1] =  1.0;
    }
    ROL::DefaultAlgorithm<RealT> TR_algo(TR_step,status);
    TR_algo.run(x, *obj);
    e.set(x);
    e.axpy(-1.0,z);
    *outStream << "\nNorm of Error: " << e.norm() << "\n";

    //Teuchos::RCP<ROL::Algorithm<RealT> > algo;
    //ROL::DefaultAlgorithmFactory<RealT> algoFactory.getAlgo(algo, parlist);
    
    obj->gradient(y,x);

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

