// @HEADER
// ************************************************************************
// @HEADER


/*! \file  test_12.cpp
    \brief Check derivative checks.
*/

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
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

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag  = 0;

  // Specify interval on which to generate uniform random numbers.
  RealT left = -1.0, right = 1.0;

  // *** Test body.

  try {

    int dim = 100;
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> x(x_rcp);
    ROL::StdVector<RealT> y(y_rcp);

    // set x,y
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i] = 2.0;
      (*y_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }

    Objective_Rosenbrock<RealT> obj;
    //Objective_SumSquares<RealT> obj;

    std::vector<std::vector<RealT> > gCheck = obj.checkGradient(x, y);

    for (unsigned i=0; i<gCheck.size(); i++) {
      if (i==0) {
       	std::cout << std::right
                  << std::setw(20) << "Step size"
                  << std::setw(20) << "grad'*dir"
                  << std::setw(20) << "FD approx"
                  << std::setw(20) << "abs error"
       	       	  << "\n";
      }
      std::cout << std::scientific << std::setprecision(8) << std::right
                  << std::setw(20) << gCheck[i][0]
                  << std::setw(20) << gCheck[i][1]
                  << std::setw(20) << gCheck[i][2]
                  << std::setw(20) << gCheck[i][3]
                  << "\n";
    }
  
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return 0;

}

