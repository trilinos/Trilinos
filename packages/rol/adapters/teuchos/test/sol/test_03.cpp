// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_ParameterList.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_Types.hpp"

#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"

#include "ROL_OptimizationSolver.hpp"
#include "ROL_RiskMeasureFactory.hpp"

typedef double RealT;

template<class Real>
class ParametrizedObjectiveEx3 : public ROL::Objective<Real> {
public:
  Real value( const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > ex
      = dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    Real quad(0), lin(0);
    std::vector<Real> p = ROL::Objective<Real>::getParameter();
    unsigned size = static_cast<unsigned>(ex->size());
    for ( unsigned i = 0; i < size; i++ ) {
      quad += (*ex)[i]*(*ex)[i];
      lin  += (*ex)[i]*p[i+1];
    }
    return std::exp(p[0])*quad + lin + p[size+1];
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > ex
      = dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<std::vector<Real> > eg
      = dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    std::vector<Real> p = ROL::Objective<Real>::getParameter();
    unsigned size = static_cast<unsigned>(ex->size());
    const Real two(2);
    for ( unsigned i = 0; i < size; i++ ) {
      (*eg)[i] = two*std::exp(p[0])*(*ex)[i] + p[i+1];
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > ev
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<std::vector<Real> > ehv
      = dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    std::vector<Real> p = ROL::Objective<Real>::getParameter();
    unsigned size = static_cast<unsigned>(ev->size());
    const Real two(2);
    for ( unsigned i = 0; i < size; i++ ) {
      (*ehv)[i] = two*std::exp(p[0])*(*ev)[i];
    }
  }
};

void setUpAndSolve(ROL::ParameterList &list,
                   ROL::Ptr<ROL::Objective<RealT> > &pObj,
                   ROL::Ptr<ROL::SampleGenerator<RealT> > &sampler,
                   ROL::Ptr<ROL::Vector<RealT> > &x,
                   ROL::Ptr<ROL::BoundConstraint<RealT> > &bnd,
                   std::ostream & outStream) {
  x->zero();
  ROL::OptimizationProblem<RealT> problem(pObj,x,bnd);
  problem.setStochasticObjective(list,sampler);
  //outStream << "\nCheck Derivatives of Stochastic Objective Function\n";
  //problem.check(outStream);
  ROL::OptimizationSolver<RealT> solver(problem,list);
  solver.solve(outStream);
}

template<class Real>
Real random(const ROL::Ptr<const Teuchos::Comm<int> > &commptr) {
  Real val(0);
  if ( Teuchos::rank<int>(*commptr)==0 ) {
    srand(time(NULL));
    val = (Real)rand()/(Real)RAND_MAX;
  }
  Teuchos::broadcast<int,Real>(*commptr,0,&val);
  return val;
}

void setRandomVector(std::vector<RealT> &x,
               const ROL::Ptr<const Teuchos::Comm<int> > &commptr) {
  unsigned dim = static_cast<unsigned>(x.size());
  for ( unsigned i = 0; i < dim; i++ ) {
    x[i] = random<RealT>(commptr);
  }
}

void printSolution(const std::vector<RealT> &x,
                   std::ostream & outStream) {
  unsigned dim = static_cast<unsigned>(x.size());
  outStream << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    outStream << x[i] << ", ";
  }
  outStream << x[dim-1] << ")\n";
}

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > commptr =
    ROL::toPtr(Teuchos::DefaultComm<int>::getComm());

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0 && commptr->getRank()==0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input_03.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    ROL::ParameterList list = *parlist;
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
    //setRandomVector(*x_ptr,commptr);
    //setRandomVector(*d_ptr,commptr);
    // Build samplers
    int nSamp = 1000;
    unsigned sdim = dim + 2;
    std::vector<RealT> tmp(2,0.); tmp[0] = -1.; tmp[1] = 1.;
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT> > bman =
      ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(commptr);
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
    // Storage for solutions
    std::vector<std::tuple<std::string,std::string,std::vector<RealT>>> solution;
    /**********************************************************************************************/
    /************************* MEAN VALUE *********************************************************/
    /**********************************************************************************************/
    *outStream << "\nMean Value\n";
    list.sublist("SOL").set("Type","Mean Value"); 
    //setRandomVector(*x_ptr,commptr);
    setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("","Mean Value",*x_ptr));
    /**********************************************************************************************/
    /************************* RISK NEUTRAL *******************************************************/
    /**********************************************************************************************/
    *outStream << "\nRisk Neutral\n";
    list.sublist("SOL").set("Type","Risk Neutral"); 
    //setRandomVector(*x_ptr,commptr);
    setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("","Risk Neutral",*x_ptr));
    /**********************************************************************************************/
    /************************* RISK AVERSE ********************************************************/
    /**********************************************************************************************/
    for (ROL::ERiskMeasure er = ROL::RISKMEASURE_CVAR; er != ROL::RISKMEASURE_LAST; er++) {
      std::string name = ROL::ERiskMeasureToString(er);
      *outStream << std::endl << name << std::endl;
      list.sublist("SOL").set("Type","Risk Averse"); 
      list.sublist("SOL").sublist("Risk Measure").set("Name",name);
      if (er == ROL::RISKMEASURE_MEANDEVIATION           ||
          er == ROL::RISKMEASURE_MEANVARIANCE            ||
          er == ROL::RISKMEASURE_MEANDEVIATIONFROMTARGET ||
          er == ROL::RISKMEASURE_MEANVARIANCEFROMTARGET) {
        list.sublist("SOL").sublist("Risk Measure").sublist(name).set("Deviation Type","Absolute");
        //setRandomVector(*x_ptr,commptr);
        setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Risk",name,*x_ptr));

        list.sublist("SOL").sublist("Risk Measure").sublist(name).set("Deviation Type","Upper");
        //setRandomVector(*x_ptr,commptr);
        setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Risk",name,*x_ptr));
      }
      else if (er == ROL::RISKMEASURE_CHEBYSHEVSPECTRAL) {
        list.sublist("SOL").sublist("Risk Measure").sublist(name).set("Weight Type",1);
        //setRandomVector(*x_ptr,commptr);
        setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Risk",name,*x_ptr));

        list.sublist("SOL").sublist("Risk Measure").sublist(name).set("Weight Type",2);
        //setRandomVector(*x_ptr,commptr);
        setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Risk",name,*x_ptr));

        list.sublist("SOL").sublist("Risk Measure").sublist(name).set("Weight Type",3);
        //setRandomVector(*x_ptr,commptr);
        setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Risk",name,*x_ptr));
      }
      else {
        //setRandomVector(*x_ptr,commptr);
        setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Risk",name,*x_ptr));
      }
    }
    /**********************************************************************************************/
    /************************* CONVEX COMBINATION OF RISK MEASURES ********************************/
    /**********************************************************************************************/
    *outStream << "\nConvex Combination if Risk Measures\n";
    list.sublist("SOL").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Convex Combination Risk Measure");
    //setRandomVector(*x_ptr,commptr);
    setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Risk","Convex Combination of Risk Measures",*x_ptr));
    /**********************************************************************************************/
    /************************* DEVIATION **********************************************************/
    /**********************************************************************************************/
    for (ROL::EDeviationMeasure ed = ROL::DEVIATIONMEASURE_MEANVARIANCEQUADRANGLE; ed != ROL::DEVIATIONMEASURE_LAST; ed++) {
      std::string name = ROL::EDeviationMeasureToString(ed);
      *outStream << std::endl << "Deviation: " << name << std::endl;
      list.sublist("SOL").set("Type","Deviation"); 
      list.sublist("SOL").sublist("Deviation Measure").set("Name",name);
      //setRandomVector(*x_ptr);
      setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
      printSolution(*x_ptr,*outStream);
      solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Deviation",name,*x_ptr));
    }
    /**********************************************************************************************/
    /************************* REGRET *************************************************************/
    /**********************************************************************************************/
    for (ROL::ERegretMeasure er = ROL::REGRETMEASURE_MEANABSOLUTELOSS; er != ROL::REGRETMEASURE_LAST; er++) {
      std::string name = ROL::ERegretMeasureToString(er);
      *outStream << std::endl << "Regret: " << name << std::endl;
      list.sublist("SOL").set("Type","Regret"); 
      list.sublist("SOL").sublist("Regret Measure").set("Name",name);
      //setRandomVector(*x_ptr);
      setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
      printSolution(*x_ptr,*outStream);
      solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Regret",name,*x_ptr));
    }
    /**********************************************************************************************/
    /************************* PROBABILITY ********************************************************/
    /**********************************************************************************************/
    for (ROL::EProbability ep = ROL::PROBABILITY_BPOE; ep != ROL::PROBABILITY_LAST; ep++) {
      std::string name = ROL::EProbabilityToString(ep);
      *outStream << std::endl << "Probability: " << name << std::endl;
      list.sublist("SOL").set("Type","Probability"); 
      list.sublist("SOL").sublist("Probability").set("Name",name);
      //setRandomVector(*x_ptr);
      setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
      printSolution(*x_ptr,*outStream);
      solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>>("Probability",name,*x_ptr));
    }

    std::vector<std::tuple<std::string,std::string,std::vector<RealT>>>::iterator it;
    *outStream << std::endl << std::scientific << std::setprecision(6);
    for (it = solution.begin(); it != solution.end(); ++it) {
      *outStream << "  ";
      *outStream << std::setw(20) << std::left << std::get<0>(*it);
      *outStream << std::setw(50) << std::left << std::get<1>(*it);
      for (unsigned i = 0; i < dim; ++i) {
        *outStream << std::setw(18) << std::left << (std::get<2>(*it))[i];
      }
      *outStream << std::endl;
    }
    *outStream << std::endl;
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
