// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
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

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Stream.hpp"
#include "ROL_NewOptimizationSolver.hpp"

#include "ROL_StdBoundConstraint.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_Bounds.hpp"

#include "ROL_StochasticProblem.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"

#include <numeric>
#include <algorithm>

//#include <fenv.h>

typedef double RealT;

template<typename Real>
class LossEx16 : public ROL::StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    std::vector<Real> p = ROL::StdObjective<Real>::getParameter();
    return std::inner_product(x.begin(),x.end(),p.begin(),static_cast<Real>(0));
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    std::vector<Real> p = ROL::StdObjective<Real>::getParameter();
    g.assign(p.begin(),p.end());
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    std::fill(hv.begin(),hv.end(),static_cast<Real>(0));
  }
};

template<typename Real>
class BudgetEx16 : public ROL::StdConstraint<Real> {
public:
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    c[0] = std::accumulate(x.begin(),x.end(),static_cast<Real>(0)) - static_cast<Real>(1);
  }

  void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    jv[0] = std::accumulate(v.begin(),v.end(),static_cast<Real>(0));
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    std::fill(ajv.begin(),ajv.end(),v[0]);
  }

  void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    std::fill(ahuv.begin(),ahuv.end(),static_cast<Real>(0));
  }
};

void printSolution(const std::vector<RealT> &x,
                   std::ostream & outStream) {
  size_t dim = x.size();
  outStream << std::endl;
  outStream << "x = (";
  for ( size_t i = 0; i < dim-1; i++ ) {
    outStream << x[i] << ", ";
  }
  outStream << x[dim-1] << ")\n";
  outStream << std::endl;
}

int main(int argc, char* argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
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
    std::string filename = "input_16.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    ROL::ParameterList list = *parlist;
    list.sublist("General").set("Output Level",print ? 1 : 0);
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    size_t dim = 10;
    ROL::Ptr<ROL::StdVector<RealT>>
      x = ROL::makePtr<ROL::StdVector<RealT>>(dim,0);
    // Build samplers
    int nsamp = 1000;
    RealT offset(0);
    std::vector<std::vector<RealT>> bounds;
    for (size_t i = 0; i < dim; ++i) {
      RealT li = -static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      RealT ui =  static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      bounds.push_back({li,ui});
      *outStream << "  Asset "        << i << ":"
                 << "  Lower Bound: " << li
                 << "  Upper Bound: " << ui << std::endl;
      offset += static_cast<RealT>(1)/static_cast<RealT>(dim) * ui;
    }
    ROL::Ptr<ROL::BatchManager<RealT>>
      bman = ROL::makePtr<ROL::BatchManager<RealT>>();
    ROL::Ptr<ROL::SampleGenerator<RealT>>
      sampler = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);
    // Build loss function
    ROL::Ptr<ROL::Objective<RealT>>
      loss = ROL::makePtr<LossEx16<RealT>>();
    loss->setParameter(sampler->getMyPoint(0));
    // Build loss constraint
    ROL::Ptr<ROL::Constraint<RealT>>
      losscon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(loss);
    ROL::Ptr<ROL::SingletonVector<RealT>>
      mul_losscon = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
    ROL::Ptr<ROL::SingletonVector<RealT>>
      u_losscon = ROL::makePtr<ROL::SingletonVector<RealT>>(offset);
    ROL::Ptr<ROL::BoundConstraint<RealT>>
      bnd_losscon = ROL::makePtr<ROL::Bounds<RealT>>(*u_losscon,false);
    // Build budget constraint
    ROL::Ptr<ROL::Constraint<RealT>>
      budget = ROL::makePtr<BudgetEx16<RealT>>();
    ROL::Ptr<ROL::StdVector<RealT>>
      mul_budget = ROL::makePtr<ROL::StdVector<RealT>>(1,0);
    // Build optimization problem
    ROL::Ptr<ROL::StochasticProblem<RealT>>
      problem = ROL::makePtr<ROL::StochasticProblem<RealT>>(loss,x);
    ROL::Ptr<ROL::NewOptimizationProblem<RealT>>
      newprob = ROL::staticPtrCast<ROL::NewOptimizationProblem<RealT>>(problem);
    problem->addLinearConstraint("Budget",budget,mul_budget);
    problem->addLinearConstraint("Loss",losscon,mul_losscon,bnd_losscon);

    // Check deterministic problem
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);

    // Markowitz portfolio selection
    problem->edit();
    ROL::ParameterList objlist;
    objlist.sublist("SOL").sublist("Objective").set("Type","Deviation");
    objlist.sublist("SOL").sublist("Objective").sublist("Deviation Measure").set("Name","Variance");
    problem->makeObjectiveStochastic(objlist,sampler);
    ROL::ParameterList conlist;
    conlist.sublist("SOL").sublist("Loss").set("Type","Risk Neutral");
    problem->makeLinearConstraintStochastic("Loss",conlist,sampler,bman);
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);
    ROL::Ptr<ROL::NewOptimizationSolver<RealT>>
      solver = ROL::makePtr<ROL::NewOptimizationSolver<RealT>>(newprob,*parlist);
    solver->solve(*outStream);
    errorFlag += (solver->getAlgorithmState()->statusFlag == ROL::EXITSTATUS_CONVERGED ? 0 : 1);
    printSolution(*x->getVector(),*outStream);
    ROL::Ptr<ROL::Vector<RealT>> xm = x->clone(); xm->set(*x);

    // Markwotiz portfolio selection again
    problem->edit();
    problem->resetStochasticLinearConstraint("Loss");
    conlist.sublist("SOL").sublist("Loss").set("Type","Mean Value");
    problem->makeLinearConstraintStochastic("Loss",conlist,sampler);
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);
    solver = ROL::makePtr<ROL::NewOptimizationSolver<RealT>>(newprob,*parlist);
    solver->solve(*outStream);
    errorFlag += (solver->getAlgorithmState()->statusFlag == ROL::EXITSTATUS_CONVERGED ? 0 : 1);
    printSolution(*x->getVector(),*outStream);
    xm->axpy(-1.0,*x);
    errorFlag += (xm->norm() > std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 1 : 0);

    // bPOE constrained portfolio selection
    problem->edit();
    conlist.sublist("SOL").sublist("Loss").set("Type","Probability");
    conlist.sublist("SOL").sublist("Loss").sublist("Probability").set("Name","bPOE");
    conlist.sublist("SOL").sublist("Loss").sublist("Probability").sublist("bPOE").set("Threshold",0.0);
    conlist.sublist("SOL").sublist("Loss").sublist("Probability").sublist("bPOE").set("Moment Order",2.0);
    problem->removeLinearConstraint("Loss");
    problem->addConstraint("Loss",losscon,mul_losscon,bnd_losscon);
    problem->makeConstraintStochastic("Loss",conlist,sampler);
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);
    solver = ROL::makePtr<ROL::NewOptimizationSolver<RealT>>(newprob,*parlist);
    solver->solve(*outStream);
    errorFlag += (solver->getAlgorithmState()->statusFlag == ROL::EXITSTATUS_CONVERGED ? 0 : 1);
    printSolution(*x->getVector(),*outStream);

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

//
//
//    // Test parametrized objective functions
//    *outStream << "Check Derivatives of Parametrized Objective Function\n";
//    pObj->setParameter(sampler->getMyPoint(0));
//    pObj->checkGradient(*x,*d,true,*outStream);
//    pObj->checkHessVec(*x,*d,true,*outStream);
//    // Storage for solutions
//    std::vector<std::tuple<std::string,std::string,std::vector<RealT>,bool>> solution;
//    bool flag = false;
//    /**********************************************************************************************/
//    /************************* MEAN VALUE *********************************************************/
//    /**********************************************************************************************/
//    *outStream << "\nMean Value\n";
//    list.sublist("SOL").sublist("Objective").set("Type","Mean Value"); 
//    //setRandomVector(*x_ptr);
//    flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//    printSolution(*x_ptr,*outStream);
//    solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("","Mean Value",*x_ptr,flag));
//    /**********************************************************************************************/
//    /************************* RISK NEUTRAL *******************************************************/
//    /**********************************************************************************************/
//    *outStream << "\nRisk Neutral\n";
//    list.sublist("SOL").sublist("Objective").set("Type","Risk Neutral"); 
//    //setRandomVector(*x_ptr);
//    flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//    printSolution(*x_ptr,*outStream);
//    solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("","Risk Neutral",*x_ptr,flag));
//    /**********************************************************************************************/
//    /************************* RISK AVERSE ********************************************************/
//    /**********************************************************************************************/
//    for (ROL::ERiskMeasure er = ROL::RISKMEASURE_CVAR; er != ROL::RISKMEASURE_LAST; er++) {
//      std::string name = ROL::ERiskMeasureToString(er);
//      *outStream << std::endl << "Risk Averse: " << name << std::endl;
//      list.sublist("SOL").sublist("Objective").set("Type","Risk Averse"); 
//      list.sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name",name);
//      if (er == ROL::RISKMEASURE_MEANDEVIATION           ||
//          er == ROL::RISKMEASURE_MEANVARIANCE            ||
//          er == ROL::RISKMEASURE_MEANDEVIATIONFROMTARGET ||
//          er == ROL::RISKMEASURE_MEANVARIANCEFROMTARGET) {
//        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Deviation Type","Absolute");
//        //setRandomVector(*x_ptr);
//        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//        printSolution(*x_ptr,*outStream);
//        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));
//
//        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Deviation Type","Upper");
//        //setRandomVector(*x_ptr);
//        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//        printSolution(*x_ptr,*outStream);
//        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));
//      }
//      else if (er == ROL::RISKMEASURE_CHEBYSHEVSPECTRAL) {
//        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Weight Type",1);
//        //setRandomVector(*x_ptr);
//        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//        printSolution(*x_ptr,*outStream);
//        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));
//
//        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Weight Type",2);
//        //setRandomVector(*x_ptr);
//        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//        printSolution(*x_ptr,*outStream);
//        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));
//
//        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(name).set("Weight Type",3);
//        //setRandomVector(*x_ptr);
//        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//        printSolution(*x_ptr,*outStream);
//        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));
//      }
//      else {
//        //setRandomVector(*x_ptr);
//        flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//        printSolution(*x_ptr,*outStream);
//        solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk",name,*x_ptr,flag));
//      }
//    }
//    /**********************************************************************************************/
//    /************************* CONVEX COMBINATION OF RISK MEASURES ********************************/
//    /**********************************************************************************************/
//    *outStream << "\nRisk Averse: Convex Combination of Risk Measures\n";
//    list.sublist("SOL").sublist("Objective").set("Type","Risk Averse"); 
//    list.sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name","Convex Combination Risk Measure");
//    //setRandomVector(*x_ptr);
//    flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//    printSolution(*x_ptr,*outStream);
//    solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Risk","Convex Combination of Risk Measures",*x_ptr,flag));
//    /**********************************************************************************************/
//    /************************* DEVIATION **********************************************************/
//    /**********************************************************************************************/
//    RealT tol = list.sublist("Status Test").get("Gradient Tolerance",1e-6);
//    for (ROL::EDeviationMeasure ed = ROL::DEVIATIONMEASURE_MEANVARIANCEQUADRANGLE; ed != ROL::DEVIATIONMEASURE_LAST; ed++) {
//      std::string name = ROL::EDeviationMeasureToString(ed);
//      *outStream << std::endl << "Deviation: " << name << std::endl;
//      list.sublist("SOL").sublist("Objective").set("Type","Deviation"); 
//      list.sublist("SOL").sublist("Objective").sublist("Deviation Measure").set("Name",name);
//      list.sublist("Status Test").set("Gradient Tolerance",tol);
//      if (ed == ROL::DEVIATIONMEASURE_LOGQUANTILEQUADRANGLE)
//        list.sublist("Status Test").set("Gradient Tolerance",std::max(1e-3,tol));
//      //setRandomVector(*x_ptr);
//      flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//      printSolution(*x_ptr,*outStream);
//      solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Deviation",name,*x_ptr,flag));
//    }
//    /**********************************************************************************************/
//    /************************* REGRET *************************************************************/
//    /**********************************************************************************************/
//    for (ROL::ERegretMeasure er = ROL::REGRETMEASURE_MEANABSOLUTELOSS; er != ROL::REGRETMEASURE_LAST; er++) {
//      std::string name = ROL::ERegretMeasureToString(er);
//      *outStream << std::endl << "Regret: " << name << std::endl;
//      list.sublist("SOL").sublist("Objective").set("Type","Regret"); 
//      list.sublist("SOL").sublist("Objective").sublist("Regret Measure").set("Name",name);
//      //setRandomVector(*x_ptr);
//      flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//      printSolution(*x_ptr,*outStream);
//      solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Regret",name,*x_ptr,flag));
//    }
//    /**********************************************************************************************/
//    /************************* PROBABILITY ********************************************************/
//    /**********************************************************************************************/
//    for (ROL::EProbability ep = ROL::PROBABILITY_BPOE; ep != ROL::PROBABILITY_LAST; ep++) {
//      std::string name = ROL::EProbabilityToString(ep);
//      *outStream << std::endl << "Probability: " << name << std::endl;
//      list.sublist("SOL").sublist("Objective").set("Type","Probability"); 
//      list.sublist("SOL").sublist("Objective").sublist("Probability").set("Name",name);
//      //setRandomVector(*x_ptr);
//      flag = setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
//      printSolution(*x_ptr,*outStream);
//      solution.push_back(std::tuple<std::string,std::string,std::vector<RealT>,bool>("Probability",name,*x_ptr,flag));
//    }
//
//    *outStream << std::endl << std::scientific << std::setprecision(6);
//    for (auto it = solution.begin(); it != solution.end(); ++it) {
//      *outStream << "  ";
//      *outStream << std::setw(20) << std::left << std::get<0>(*it);
//      *outStream << std::setw(50) << std::left << std::get<1>(*it);
//      for (unsigned i = 0; i < dim; ++i) {
//        *outStream << std::setw(18) << std::left << (std::get<2>(*it))[i];
//      }
//      *outStream << std::setw(10) << std::left << (std::get<3>(*it) ? "CONVERGED" : "NOT CONVERGED");
//      *outStream << std::endl;
//    }
//    *outStream << std::endl;
