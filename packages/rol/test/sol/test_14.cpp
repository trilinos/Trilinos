// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_Types.hpp"

#include "ROL_ProgressiveHedging.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"

typedef double RealT;

template<class Real> 
class ParametrizedObjectiveEx1 : public ROL::Objective<Real> {
public:
  Real value( const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > ex = 
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    Real quad = 0.0, lin = 0.0;
    std::vector<Real> p = this->getParameter();
    unsigned size = ex->size();
    for ( unsigned i = 0; i < size; i++ ) {
      quad += (*ex)[i]*(*ex)[i]; 
      lin  += (*ex)[i]*p[i+1];
    }
    return std::exp(p[0])*quad + lin + p[size+1];
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > ex = 
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<std::vector<Real> > eg =
      dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    std::vector<Real> p = this->getParameter();
    unsigned size = ex->size();
    for ( unsigned i = 0; i < size; i++ ) {
      (*eg)[i] = 2.0*std::exp(p[0])*(*ex)[i] + p[i+1];
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > ex = 
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<const std::vector<Real> > ev = 
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<std::vector<Real> > ehv =
      dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    std::vector<Real> p = this->getParameter();
    unsigned size = ex->size();
    for ( unsigned i = 0; i < size; i++ ) {
      (*ehv)[i] = 2.0*std::exp(p[0])*(*ev)[i]; 
    } 
  }
};

void setRandomVector(std::vector<RealT> &x) {
  unsigned dim = x.size();
  for ( unsigned i = 0; i < dim; i++ ) {
    x[i] = (RealT)rand()/(RealT)RAND_MAX;
  }
}

void printSolution(const std::vector<RealT> &x,
                   std::ostream & outStream) {
  unsigned dim = x.size();
  outStream << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    outStream << x[i] << ", ";
  }
  outStream << x[dim-1] << ")\n";
}

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input_14.xml";
//    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
//    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
//    Teuchos::ParameterList list = *parlist;
      ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile( filename );
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    unsigned dim = 4;
    ROL::Ptr<std::vector<RealT>> x_ptr  = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT>> x      = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    ROL::Ptr<std::vector<RealT>> xr_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT>> xr     = ROL::makePtr<ROL::StdVector<RealT>>(xr_ptr);
    ROL::Ptr<std::vector<RealT>> xg_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT>> xg     = ROL::makePtr<ROL::StdVector<RealT>>(xg_ptr);
    ROL::Ptr<std::vector<RealT>> d_ptr  = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT>> d      = ROL::makePtr<ROL::StdVector<RealT>>(d_ptr);
    setRandomVector(*d_ptr);
    // Build samplers
    int nSamp = 1000;  
    unsigned sdim = dim + 2;
    std::vector<RealT> tmp(2,0.); tmp[0] = -1.; tmp[1] = 1.;
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT> > bman =
      ROL::makePtr<ROL::BatchManager<RealT>>();
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler =
      ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nSamp,bounds,bman,false,false,100);
    // Build stochastic objective function
    ROL::Ptr<ROL::Objective<RealT> > pObj =
      ROL::makePtr<ParametrizedObjectiveEx1<RealT>>();
    // Build bound constraints
    std::vector<RealT> l(dim,0.0);
    std::vector<RealT> u(dim,1.0);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd = 
      ROL::makePtr<ROL::StdBoundConstraint<RealT>>(l,u);
    bnd->deactivate();
    // Test parametrized objective functions
    *outStream << "Check Derivatives of Parametrized Objective Function\n";
    pObj->setParameter(sampler->getMyPoint(0));
    pObj->checkGradient(*x,*d,true,*outStream);
    pObj->checkHessVec(*x,*d,true,*outStream);
    // Solve using progressive hedging
    ROL::Ptr<ROL::Problem<RealT>> problem
      = ROL::makePtr<ROL::Problem<RealT>>(pObj,x);
    problem->addBoundConstraint(bnd);
    //ROL::Ptr<ROL::OptimizationProblem<RealT>> problem
    //  = ROL::makePtr<ROL::OptimizationProblem<RealT>>(pObj,x,bnd);
    ROL::ProgressiveHedging<RealT> ph_solver(problem,sampler,*parlist);
    ph_solver.check(*outStream);
    ph_solver.run(*outStream);
    // Compute true solution
    RealT denom(0), denom_g(0);
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      denom += sampler->getMyWeight(i)*std::exp(sampler->getMyPoint(i)[0]);
      for (unsigned j = 0; j < dim; ++j) {
        (*xr_ptr)[j] += sampler->getMyWeight(i)*sampler->getMyPoint(i)[j+1];
      }
    }
    sampler->sumAll(&denom,&denom_g,1);
    sampler->sumAll(*xr,*xg);
    for (unsigned j = 0; j < dim; ++j) {
      (*xg_ptr)[j] /= -2.0*denom_g;
    }
    // Print solutions
    *outStream << std::endl << "Progressive Hedging Solution" << std::endl;
    printSolution(*x_ptr,*outStream);
    *outStream << std::endl << "True Solution" << std::endl;
    printSolution(*xg_ptr,*outStream);
    d->set(*x); d->axpy(-1.0,*xg);
    RealT err = d->norm();
    *outStream << std::endl << "Error: " << err << std::endl << std::endl;
    errorFlag += (err > 1e-2 ? 1 : 0);
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
