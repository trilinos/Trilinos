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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Types.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"

#include "ROL_CVaRVector.hpp"
#include "ROL_ParametrizedObjective.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_AbsoluteValue.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_RiskMeasure.hpp"
#include "ROL_MeanDeviation.hpp"
#include "ROL_MeanVariance.hpp"
#include "ROL_MeanDeviationFromTarget.hpp"
#include "ROL_MeanVarianceFromTarget.hpp"
#include "ROL_CVaR.hpp"
#include "ROL_SmoothCVaRQuad.hpp"
#include "ROL_ExpUtility.hpp"
#include "ROL_RiskAverseObjective.hpp"
#include "ROL_StdEpetraBatchManager.hpp"

template<class Real> 
class ParametrizedObjectiveEx1 : public ROL::ParametrizedObjective<Real> {
public:
  ParametrizedObjectiveEx1( ) {}

  Real value( const ROL::Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > ex = 
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Real ip = 0.0;
    for ( unsigned i = 0; i < ex->size(); i++ ) {
      ip += (*ex)[i]*(this->getParameter())[i];
    }
    return -std::cos(ip);
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > ex = 
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > eg =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());
    Real ip = 0.0;
    for ( unsigned i = 0; i < ex->size(); i++ ) {
      ip += (*ex)[i]*(this->getParameter())[i];
    }
    ip = std::sin(ip);
    for ( unsigned i = 0; i < ex->size(); i++ ) {
      (*eg)[i] = ip*(this->getParameter())[i];
    } 
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > ex = 
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > ev = 
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<std::vector<Real> > ehv =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());
    Real ip1 = 0.0;
    Real ip2 = 0.0;
    for ( unsigned i = 0; i < ex->size(); i++ ) {
      ip1 += (*ex)[i]*(this->getParameter())[i];
      ip2 += (*ev)[i]*(this->getParameter())[i];
    }
    ip1 = std::cos(ip1);
    for ( unsigned i = 0; i < ex->size(); i++ ) {
      (*ehv)[i] = ip1*ip2*(this->getParameter())[i]; 
    } 
  }
};

int main(int argc, char* argv[]) {

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

  /**********************************************************************************************/
  /************************* CONSTRUCT ROL ALGORITHM ********************************************/
  /**********************************************************************************************/
  // Get ROL parameterlist
  std::string filename = "input.xml";
  Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
  Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist) );
  // Build ROL algorithm
  double gtol = parlist->get("Gradient Tolerance",1.e-6);
  double stol = parlist->get("Step Tolerance",1.e-12);
  int maxit   = parlist->get("Maximum Number of Iterations",100);
  ROL::StatusTest<double> status(gtol,stol,maxit);
  //ROL::LineSearchStep<double> step(*parlist);
  Teuchos::RCP<ROL::Step<double> > step;
  Teuchos::RCP<ROL::DefaultAlgorithm<double> > algo;
  /**********************************************************************************************/
  /************************* CONSTRUCT SOL COMPONENTS *******************************************/
  /**********************************************************************************************/
  // Build vectors
  unsigned dim = 2;
  Teuchos::RCP<std::vector<double> > x_rcp = Teuchos::rcp( new std::vector<double>(dim,0.0) );
  ROL::StdVector<double> x(x_rcp);
  Teuchos::RCP<std::vector<double> > d_rcp = Teuchos::rcp( new std::vector<double>(dim,0.0) );
  ROL::StdVector<double> d(d_rcp);
  for ( unsigned i = 0; i < dim; i++ ) {
    (*d_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  // Build samplers
  int nSamp = 1000;  
  std::vector<double> tmp(2,0.0);
  tmp[0] = -1.0; tmp[1] = 1.0;
  std::vector<std::vector<double> > bounds(dim,tmp);
  Teuchos::RCP<ROL::BatchManager<double> > bman =
    Teuchos::rcp(new ROL::BatchManager<double>());
  ROL::MonteCarloGenerator<double> vsampler(nSamp,bounds,bman,false,false,100);
  /*std::vector<double> mean(dim,0.0); mean[1] = 1.0;
  std::vector<double> std(dim,1.0);  std[1]  = 5.0;
  ROL::MonteCarloGenerator<double> gsampler(nSamp,mean,std,bman);
  std::ofstream file;
  std::stringstream name;
  file.open(name.str().c_str());
  for ( unsigned i = 0; i < gsampler.numMySamples(); i++ ) {
    for ( unsigned j = 0; j < dim; j++ ) { 
      file << (gsampler.getMyPoint(i))[j] << "  ";
    }
    file << "\n";
  }
  file.close();*/
  // Build risk-averse objective function
  ParametrizedObjectiveEx1<double> pObj;
  Teuchos::RCP<ROL::RiskMeasure<double> > rm;
  Teuchos::RCP<ROL::RiskAverseObjective<double> > obj;
  // Test parametrized objective functions
  std::cout << "Check Derivatives of Parametrized Objective Function\n";
  pObj.setParameter(vsampler.getMyPoint(0));
  pObj.checkGradient(x,d,true);
  pObj.checkHessVec(x,d,true);
  /**********************************************************************************************/
  /************************* RISK NEUTRAL *******************************************************/
  /**********************************************************************************************/
  std::cout << "\nRISK NEUTRAL\n";
  rm  = Teuchos::rcp( new ROL::RiskMeasure<double>() );
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* MEAN PLUS DEVIATION ************************************************/
  /**********************************************************************************************/
  std::cout << "\nMEAN PLUS DEVIATION\n";
  // Absolute value approximation
  double gamma = 0.0;
  ROL::EAbsoluteValue eav = ROL::ABSOLUTEVALUE_C2;
  Teuchos::RCP<ROL::PositiveFunction<double> > pf = Teuchos::rcp( new ROL::AbsoluteValue<double>(gamma,eav) );
  // Moment vector
  std::vector<double> order(3,0.0); order[0] = 2.0; order[1] = 3.0; order[2] = 4.0;
  // Moment coefficients
  std::vector<double> coeff(3,0.5);
  rm  = Teuchos::rcp( new ROL::MeanDeviation<double>(order,coeff,pf) );
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* MEAN PLUS VARIANCE *************************************************/
  /**********************************************************************************************/
  std::cout << "\nMEAN PLUS VARIANCE\n";
  rm  = Teuchos::rcp( new ROL::MeanVariance<double>(order,coeff,pf) );
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* MEAN PLUS DEVIATION FROM TARGET ************************************/
  /**********************************************************************************************/
  std::cout << "\nMEAN PLUS DEVIATION FROM TARGET\n";
  // Moment targets
  std::vector<double> target(3,-0.5);
  // Risk measure
  rm  = Teuchos::rcp( new ROL::MeanDeviationFromTarget<double>(target,order,coeff,pf) );
  // Risk averse objective
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* MEAN PLUS VARIANCE FROM TARGET *************************************/
  /**********************************************************************************************/
  std::cout << "\nMEAN PLUS VARIANCE FROM TARGET\n";
  // Risk measure
  rm  = Teuchos::rcp( new ROL::MeanVarianceFromTarget<double>(target,order,coeff,pf) );
  // Risk averse objective
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* MEAN PLUS SEMIDEVIATION ********************************************/
  /**********************************************************************************************/
  std::cout << "\nMEAN PLUS SEMIDEVIATION\n";
  // Plus function approximation
  gamma = 1.e2;
  std::vector<double> data2(2,0.0);
  data2[0] = 0.0; data2[1] = 1.0;
  Teuchos::RCP<ROL::Distribution<double> > dist2 =
    Teuchos::rcp(new ROL::Distribution<double>(ROL::DISTRIBUTION_PARABOLIC,data2));
  Teuchos::RCP<ROL::PlusFunction<double> > plusf =
    Teuchos::rcp(new ROL::PlusFunction<double>(dist2,1.0/gamma));
  pf = Teuchos::rcp(new ROL::PlusFunction<double>(dist2,1.0/gamma));
  // Risk measure
  rm  = Teuchos::rcp( new ROL::MeanDeviation<double>(order,coeff,pf) );
  // Risk averse objective
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* MEAN PLUS SEMIVARIANCE *********************************************/
  /**********************************************************************************************/
  std::cout << "\nMEAN PLUS SEMIVARIANCE\n";
  // Risk measure
  rm  = Teuchos::rcp( new ROL::MeanVariance<double>(order,coeff,pf) );
  // Risk averse objective
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* MEAN PLUS SEMIDEVIATION FROM TARGET ********************************/
  /**********************************************************************************************/
  std::cout << "\nMEAN PLUS SEMIDEVIATION FROM TARGET\n";
  // Risk measure
  rm  = Teuchos::rcp( new ROL::MeanDeviationFromTarget<double>(target,order,coeff,pf) );
  // Risk averse objective
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* MEAN PLUS SEMIVARIANCE FROM TARGET *********************************/
  /**********************************************************************************************/
  std::cout << "\nMEAN PLUS SEMIVARIANCE FROM TARGET\n";
  // Risk measure
  rm  = Teuchos::rcp( new ROL::MeanVarianceFromTarget<double>(target,order,coeff,pf) );
  // Risk averse objective
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* MEAN PLUS CVAR *****************************************************/
  /**********************************************************************************************/
  std::cout << "\nMEAN PLUS CONDITIONAL VALUE AT RISK\n";
  double prob = 0.8;
  double c = 0.8;
  rm  = Teuchos::rcp( new ROL::CVaR<double>(prob,c,plusf) );
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  double xv = 10.0*(double)rand()/(double)RAND_MAX-5.0, dv = 10.0*(double)rand()/(double)RAND_MAX-5.0;
  Teuchos::RCP<ROL::Vector<double> > xp = Teuchos::rcp(&x,false);
  Teuchos::RCP<ROL::Vector<double> > dp = Teuchos::rcp(&d,false);
  ROL::CVaRVector<double> xc(xv,xp);
  ROL::CVaRVector<double> dc(dv,dp);
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(xc,dc,true);
  obj->checkHessVec(xc,dc,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(xc,*obj,true);
  // Print Solution
  std::cout << "t = " << xc.getVaR() << "\n";
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* SMOOTHED CVAR QUADRANGLE *******************************************/
  /**********************************************************************************************/
  std::cout << "\nSMOOTHED CONDITIONAL VALUE AT RISK \n";
  prob = 0.9;
  ROL::SmoothCVaRQuad<double> scq(prob,1.0/gamma,plusf);
  scq.checkRegret();
  rm  = Teuchos::rcp(&scq,false);
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  xv = 10.0*(double)rand()/(double)RAND_MAX-5.0;
  dv = 10.0*(double)rand()/(double)RAND_MAX-5.0;
//dv = 0.0;
  xp = Teuchos::rcp(&x,false);
  dp = Teuchos::rcp(&d,false);
//dp->zero();
  ROL::CVaRVector<double> xq(xv,xp);
  ROL::CVaRVector<double> dq(dv,dp);
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(xq,dq,true);
  obj->checkHessVec(xq,dq,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(xq,*obj,true);
  // Print Solution
  std::cout << "t = " << xq.getVaR() << "\n";
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";
  /**********************************************************************************************/
  /************************* EXPONENTIAL UTILITY FUNCTION ***************************************/
  /**********************************************************************************************/
  std::cout << "\nEXPONENTIAL UTILITY FUNCTION\n";
  rm  = Teuchos::rcp( new ROL::ExpUtility<double> );
  obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,*rm,vsampler,vsampler) );
  // Test objective functions
  for ( unsigned i = 0; i < dim; i++ ) {
    (*x_rcp)[i] = 10.0*(double)rand()/(double)RAND_MAX - 5.0;
  }
  std::cout << "\nCheck Derivatives of Risk-Averse Objective Function\n";
  obj->checkGradient(x,d,true);
  obj->checkHessVec(x,d,true);
  // Run ROL algorithm
  step = Teuchos::rcp( new ROL::TrustRegionStep<double>(*parlist) );
  algo = Teuchos::rcp( new ROL::DefaultAlgorithm<double>(*step,status,false) );
  algo->run(x,*obj,true);
  // Print Solution
  std::cout << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    std::cout << (*x_rcp)[i] << ", ";
  }
  std::cout << (*x_rcp)[dim-1] << ")\n";

  return 0;
}
