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
#include "ROL_Stream.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "ROL_StdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_Types.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"

#include "ROL_Objective.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StdEpetraBatchManager.hpp"

#include "ROL_OptimizationProblem.hpp"

typedef double RealT;

template<class Real> 
class ParametrizedObjectiveEx3 : public ROL::Objective<Real> {
public:
  Real value( const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > ex = 
      (dynamic_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &&>(x))).getVector();
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
      (dynamic_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &&>(x))).getVector();
    ROL::Ptr<std::vector<Real> > eg =
      ROL::constPtrCast<std::vector<Real> >((dynamic_cast<ROL::StdVector<Real>&>(g)).getVector());
    std::vector<Real> p = this->getParameter();
    unsigned size = ex->size();
    for ( unsigned i = 0; i < size; i++ ) {
      (*eg)[i] = 2.0*std::exp(p[0])*(*ex)[i] + p[i+1];
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > ex = 
      (dynamic_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &&>(x))).getVector();
    ROL::Ptr<const std::vector<Real> > ev = 
      (dynamic_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &&>(v))).getVector();
    ROL::Ptr<std::vector<Real> > ehv =
      ROL::constPtrCast<std::vector<Real> >((dynamic_cast<ROL::StdVector<Real>&>(hv)).getVector());
    std::vector<Real> p = this->getParameter();
    unsigned size = ex->size();
    for ( unsigned i = 0; i < size; i++ ) {
      (*ehv)[i] = 2.0*std::exp(p[0])*(*ev)[i]; 
    } 
  }
};

void setUpAndSolve(Teuchos::ParameterList &list,
                   ROL::Ptr<ROL::Objective<RealT> > &pObj,
                   ROL::Ptr<ROL::SampleGenerator<RealT> > &sampler,
                   ROL::Ptr<ROL::Vector<RealT> > &x,
                   ROL::Ptr<ROL::Vector<RealT> > &d,
                   ROL::Ptr<ROL::BoundConstraint<RealT> > &bnd,
                   std::ostream & outStream) {
  ROL::OptimizationProblem<RealT> opt(pObj,x,bnd);
  opt.setStochasticObjective(list,sampler);
  outStream << "\nCheck Derivatives of Stochastic Objective Function\n";
  opt.checkObjectiveGradient(*d,true,outStream);
  opt.checkObjectiveHessVec(*d,true,outStream);
  // Run ROL algorithm
  ROL::Ptr<ROL::Step<RealT>>
    step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(list);
  ROL::Ptr<ROL::StatusTest<RealT>>
    status = ROL::makePtr<ROL::StatusTest<RealT>>(list);
  ROL::Algorithm<RealT> algo(step,status,false);
  algo.run(opt,true,outStream);
}

template<class Real>
Real random(const ROL::Ptr<Epetra_Comm> &comm) {
  Real val = 0.0;
  if ( comm->MyPID()==0 ) {
    val = (Real)rand()/(Real)RAND_MAX;
  }
  comm->Broadcast(&val,1,0);
  return val;
}

void setRandomVector(std::vector<RealT> &x,
               const ROL::Ptr<Epetra_Comm> &comm) {
  unsigned dim = x.size();
  for ( unsigned i = 0; i < dim; i++ ) {
    x[i] = random<RealT>(comm);
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
  ROL::Ptr<Epetra_Comm> comm;
#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  comm = ROL::makePtr<Epetra_MpiComm>(MPI_COMM_WORLD);
#else 
  comm = ROL::makePtr<Epetra_SerialComm>();
#endif

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0 && comm->MyPID()==0)
    ROL::makePtrFromRef(std::cout);
  else
    ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input_01.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    Teuchos::ParameterList list = *parlist;
    // Build ROL algorithm
    ROL::Ptr<ROL::Algorithm<RealT> > algo;
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    unsigned dim = 4;
    ROL::Ptr<std::vector<RealT> > x_ptr  = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > x  = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    ROL::Ptr<std::vector<RealT> > d_ptr  = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > d  = ROL::makePtr<ROL::StdVector<RealT>>(d_ptr);
    setRandomVector(*x_ptr,comm);
    setRandomVector(*d_ptr,comm);
    // Build samplers
    int nSamp = 1000;
    unsigned sdim = dim + 2;
    std::vector<RealT> tmp(2,0.); tmp[0] = -1.; tmp[1] = 1.;
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT> > bman =
      ROL::makePtr<ROL::StdEpetraBatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler =
      ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nSamp,bounds,bman,false,false,100);
    // Build risk-averse objective function
    ROL::Ptr<ROL::Objective<RealT> > pObj =
      ROL::makePtr<ParametrizedObjectiveEx3<RealT>>();
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
    /**********************************************************************************************/
    /************************* MEAN VALUE *********************************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN VALUE\n";
    list.sublist("SOL").set("Type","Mean Value"); 
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* RISK NEUTRAL *******************************************************/
    /**********************************************************************************************/
    *outStream << "\nRISK NEUTRAL\n";
    list.sublist("SOL").set("Type","Risk Neutral"); 
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS DEVIATION ************************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS DEVIATION\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Deviation");
    list.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Deviation").set("Deviation Type","Absolute");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS VARIANCE *************************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS VARIANCE\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Variance");
    list.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Variance").set("Deviation Type","Absolute");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS DEVIATION FROM TARGET ************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS DEVIATION FROM TARGET\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Deviation From Target");
    list.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Deviation From Target").set("Deviation Type","Absolute");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS VARIANCE FROM TARGET *************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS VARIANCE FROM TARGET\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Variance From Target");
    list.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Variance From Target").set("Deviation Type","Absolute");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS SEMIDEVIATION ********************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS SEMIDEVIATION\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Deviation");
    list.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Deviation").set("Deviation Type","Upper");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS SEMIVARIANCE *********************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS SEMIVARIANCE\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Variance");
    list.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Variance").set("Deviation Type","Upper");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS SEMIDEVIATION FROM TARGET ********************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS SEMIDEVIATION FROM TARGET\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Deviation From Target");
    list.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Deviation From Target").set("Deviation Type","Upper");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS SEMIVARIANCE FROM TARGET *********************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS SEMIVARIANCE FROM TARGET\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Mean Plus Variance From Target");
    list.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Variance From Target").set("Deviation Type","Upper");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS CVAR *****************************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS CONDITIONAL VALUE AT RISK\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS HMCR *****************************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS HIGHER MOMENT COHERENT RISK MEASURE\n";
    list.sublist("SOL").sublist("Risk Measure").set("Name","HMCR");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    /**********************************************************************************************/
    /************************* EXPONENTIAL UTILITY FUNCTION ***************************************/
    /**********************************************************************************************/
    *outStream << "\nEXPONENTIAL UTILITY FUNCTION\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Entropic Risk");
    setRandomVector(*x_ptr,comm);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
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
