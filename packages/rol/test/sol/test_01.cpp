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

#include "ROL_StdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_Types.hpp"

#include "ROL_OptimizationSolver.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_Solver.hpp"
#include "ROL_RiskMeasureFactory.hpp"
#include "ROL_DeviationMeasureFactory.hpp"
#include "ROL_RegretMeasureFactory.hpp"
#include "ROL_ProbabilityFactory.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"

//#include <fenv.h>

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

bool setUpAndSolve(ROL::ParameterList &list,
                   ROL::Ptr<ROL::Objective<RealT> > &pObj,
                   ROL::Ptr<ROL::SampleGenerator<RealT> > &sampler,
                   ROL::Ptr<ROL::Vector<RealT> > &x,
                   ROL::Ptr<ROL::BoundConstraint<RealT> > &bnd,
                   std::ostream & outStream) {
  x->zero();
  ROL::Ptr<ROL::StochasticProblem<RealT>>
    newprob = ROL::makePtr<ROL::StochasticProblem<RealT>>(pObj,x);
  if (bnd->isActivated()) newprob->addBoundConstraint(bnd);
  newprob->makeObjectiveStochastic(list,sampler);
  outStream << "\nCheck Derivatives of Stochastic Objective Function\n";
  newprob->finalize(false,true,outStream);
  newprob->check(true,outStream);
  ROL::Solver<RealT> solver(newprob,list);
  solver.solve(outStream);
  return solver.getAlgorithmState()->statusFlag == ROL::EXITSTATUS_CONVERGED;
}

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
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  bool print     = iprint > 0;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
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
    std::string filename = "input_01.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    ROL::ParameterList list = *parlist;
    list.sublist("General").set("Output Level",print ? 1 : 0);
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    unsigned dim = 4;
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > x = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    ROL::Ptr<std::vector<RealT> > d_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > d = ROL::makePtr<ROL::StdVector<RealT>>(d_ptr);
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
    //bnd->deactivate();
    // Test parametrized objective functions
    *outStream << "Check Derivatives of Parametrized Objective Function\n";
    pObj->setParameter(sampler->getMyPoint(0));
    pObj->checkGradient(*x,*d,true,*outStream);
    pObj->checkHessVec(*x,*d,true,*outStream);
    // Storage for solutions
    std::vector<std::tuple<std::string,std::string,std::vector<RealT>,bool>> solution;
    bool flag = false;
    /**********************************************************************************************/
    /************************* MEAN VALUE *********************************************************/
    /**********************************************************************************************/
    *outStream << "\nMean Value\n";
    list.sublist("SOL").sublist("Objective").set("Type","Mean Value"); 
    //setRandomVector(*x_ptr);
    flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("","Mean Value",*x_ptr,flag));
    /**********************************************************************************************/
    /************************* RISK NEUTRAL *******************************************************/
    /**********************************************************************************************/
    *outStream << "\nRisk Neutral\n";
    list.sublist("SOL").sublist("Objective").set("Type","Risk Neutral"); 
    //setRandomVector(*x_ptr);
    flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("","Risk Neutral",*x_ptr,flag));
    /**********************************************************************************************/
    /************************* RISK AVERSE ********************************************************/
    /**********************************************************************************************/
    for (ROL::ERiskMeasure er = ROL::RISKMEASURE_CVAR; er != ROL::RISKMEASURE_LAST; er++) {
      std::string name = ROL::ERiskMeasureToString(er);
      *outStream << std::endl << "Risk Averse: " << name << std::endl;
      list.sublist("SOL").sublist("Objective").set("Type","Risk Averse"); 
      list.sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name",name);
      if (er == ROL::RISKMEASURE_MEANDEVIATION           ||
          er == ROL::RISKMEASURE_MEANVARIANCE            ||
          er == ROL::RISKMEASURE_MEANDEVIATIONFROMTARGET ||
          er == ROL::RISKMEASURE_MEANVARIANCEFROMTARGET) {
        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Deviation Type","Absolute");
        //setRandomVector(*x_ptr);
        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));

        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Deviation Type","Upper");
        //setRandomVector(*x_ptr);
        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));
      }
      else if (er == ROL::RISKMEASURE_CHEBYSHEVSPECTRAL) {
        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Weight Type",1);
        //setRandomVector(*x_ptr);
        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));

        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Weight Type",2);
        //setRandomVector(*x_ptr);
        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));

        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Weight Type",3);
        //setRandomVector(*x_ptr);
        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));
      }
      else {
        //setRandomVector(*x_ptr);
        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
        printSolution(*x_ptr,*outStream);
        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));
      }
    }
    /**********************************************************************************************/
    /************************* CONVEX COMBINATION OF RISK MEASURES ********************************/
    /**********************************************************************************************/
    *outStream << "\nRisk Averse: Convex Combination of Risk Measures\n";
    list.sublist("SOL").sublist("Objective").set("Type","Risk Averse"); 
    list.sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name","Convex Combination Risk Measure");
    //setRandomVector(*x_ptr);
    flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
    solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk","Convex Combination of Risk Measures",*x_ptr,flag));
    /**********************************************************************************************/
    /************************* DEVIATION **********************************************************/
    /**********************************************************************************************/
    RealT tol = list.sublist("Status Test").get("Gradient Tolerance",1e-6);
    for (ROL::EDeviationMeasure ed = ROL::DEVIATIONMEASURE_MEANVARIANCEQUADRANGLE; ed != ROL::DEVIATIONMEASURE_LAST; ed++) {
      std::string name = ROL::EDeviationMeasureToString(ed);
      *outStream << std::endl << "Deviation: " << name << std::endl;
      list.sublist("SOL").sublist("Objective").set("Type","Deviation"); 
      list.sublist("SOL").sublist("Objective").sublist("Deviation Measure").set("Name",name);
      list.sublist("Status Test").set("Gradient Tolerance",tol);
      if (ed == ROL::DEVIATIONMEASURE_LOGQUANTILEQUADRANGLE)
        list.sublist("Status Test").set("Gradient Tolerance",std::max(1e-3,tol));
      //setRandomVector(*x_ptr);
      flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
      printSolution(*x_ptr,*outStream);
      solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Deviation",name,*x_ptr,flag));
    }
    /**********************************************************************************************/
    /************************* REGRET *************************************************************/
    /**********************************************************************************************/
    for (ROL::ERegretMeasure er = ROL::REGRETMEASURE_MEANABSOLUTELOSS; er != ROL::REGRETMEASURE_LAST; er++) {
      std::string name = ROL::ERegretMeasureToString(er);
      *outStream << std::endl << "Regret: " << name << std::endl;
      list.sublist("SOL").sublist("Objective").set("Type","Regret"); 
      list.sublist("SOL").sublist("Objective").sublist("Regret Measure").set("Name",name);
      //setRandomVector(*x_ptr);
      flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
      printSolution(*x_ptr,*outStream);
      solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Regret",name,*x_ptr,flag));
    }
    /**********************************************************************************************/
    /************************* PROBABILITY ********************************************************/
    /**********************************************************************************************/
    for (ROL::EProbability ep = ROL::PROBABILITY_BPOE; ep != ROL::PROBABILITY_LAST; ep++) {
      std::string name = ROL::EProbabilityToString(ep);
      *outStream << std::endl << "Probability: " << name << std::endl;
      list.sublist("SOL").sublist("Objective").set("Type","Probability"); 
      list.sublist("SOL").sublist("Objective").sublist("Probability").set("Name",name);
      //setRandomVector(*x_ptr);
      flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
      printSolution(*x_ptr,*outStream);
      solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Probability",name,*x_ptr,flag));
    }

    *outStream << std::endl << std::scientific << std::setprecision(6);
    for (auto it = solution.begin(); it != solution.end(); ++it) {
      *outStream << "  ";
      *outStream << std::setw(20) << std::left << std::get<0>(*it);
      *outStream << std::setw(50) << std::left << std::get<1>(*it);
      for (unsigned i = 0; i < dim; ++i) {
        *outStream << std::setw(18) << std::left << (std::get<2>(*it))[i];
      }
      *outStream << std::setw(10) << std::left << (std::get<3>(*it) ? "CONVERGED" : "NOT CONVERGED");
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
