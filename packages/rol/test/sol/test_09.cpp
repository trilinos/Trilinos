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
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_Types.hpp"
#include "ROL_Algorithm.hpp"

#include "ROL_StochasticProblem.hpp"
#include "ROL_ParametrizedObjective.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"

typedef double RealT;

template<class Real> 
class ParametrizedObjectiveEx8 : public ROL::ParametrizedObjective<Real> {
public:
  Real value( const ROL::Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > ex = 
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
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
    Teuchos::RCP<const std::vector<Real> > ex = 
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > eg =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());
    std::vector<Real> p = this->getParameter();
    unsigned size = ex->size();
    for ( unsigned i = 0; i < size; i++ ) {
      (*eg)[i] = 2.0*std::exp(p[0])*(*ex)[i] + p[i+1];
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > ex = 
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > ev = 
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<std::vector<Real> > ehv =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());
    std::vector<Real> p = this->getParameter();
    unsigned size = ex->size();
    for ( unsigned i = 0; i < size; i++ ) {
      (*ehv)[i] = 2.0*std::exp(p[0])*(*ev)[i]; 
    } 
  }
};

RealT setUpAndSolve(Teuchos::ParameterList &list,
                    Teuchos::RCP<ROL::ParametrizedObjective<RealT> > &pObj,
                    Teuchos::RCP<ROL::SampleGenerator<RealT> > &sampler,
                    Teuchos::RCP<ROL::Vector<RealT> > &x,
                    Teuchos::RCP<ROL::Vector<RealT> > &d,
                    Teuchos::RCP<ROL::BoundConstraint<RealT> > &bnd,
                    std::ostream & outStream) {
  ROL::StochasticProblem<RealT> opt(list,pObj,sampler,x,bnd);
  outStream << "\nCheck Derivatives of Stochastic Objective Function\n";
  opt.checkObjectiveGradient(*d,true,outStream);
  opt.checkObjectiveHessVec(*d,true,outStream);
  // Run ROL algorithm
  ROL::Algorithm<RealT> algo("Trust Region",list,false);
  algo.run(opt,true,outStream);
  Teuchos::RCP<ROL::Objective<RealT> > robj = opt.getObjective();
  RealT tol(1.e-8);
  return robj->value(*(opt.getSolutionVector()),tol);
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

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input_09.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    Teuchos::ParameterList list = *parlist;
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    unsigned dim = 4;
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT>(dim,0.0) );
    Teuchos::RCP<ROL::Vector<RealT> > x = Teuchos::rcp(new ROL::StdVector<RealT>(x_rcp));
    Teuchos::RCP<std::vector<RealT> > xp_rcp = Teuchos::rcp( new std::vector<RealT>(dim,0.0) );
    Teuchos::RCP<ROL::Vector<RealT> > xp = Teuchos::rcp(new ROL::StdVector<RealT>(xp_rcp));
    Teuchos::RCP<std::vector<RealT> > diff_rcp = Teuchos::rcp( new std::vector<RealT>(dim,0.0) );
    Teuchos::RCP<ROL::Vector<RealT> > diff = Teuchos::rcp(new ROL::StdVector<RealT>(diff_rcp));
    Teuchos::RCP<std::vector<RealT> > d_rcp = Teuchos::rcp( new std::vector<RealT>(dim,0.0) );
    Teuchos::RCP<ROL::Vector<RealT> > d = Teuchos::rcp(new ROL::StdVector<RealT>(d_rcp));
    setRandomVector(*d_rcp);
    // Build samplers
    int nSamp = 1000;  
    unsigned sdim = dim + 2;
    std::vector<RealT> tmp(2,0.); tmp[0] = -1.; tmp[1] = 1.;
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman =
      Teuchos::rcp(new ROL::BatchManager<RealT>());
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler =
      Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nSamp,bounds,bman,false,false,100));
    // Build risk-averse objective function
    Teuchos::RCP<ROL::ParametrizedObjective<RealT> > pObj =
      Teuchos::rcp(new ParametrizedObjectiveEx8<RealT>);
    // Build bound constraints
    std::vector<RealT> l(dim,0.0);
    std::vector<RealT> u(dim,1.0);
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd = 
      Teuchos::rcp( new ROL::StdBoundConstraint<RealT>(l,u) );
    bnd->deactivate();
    // Test parametrized objective functions
    *outStream << "Check Derivatives of Parametrized Objective Function\n";
    pObj->setParameter(sampler->getMyPoint(0));
    pObj->checkGradient(*x,*d,true,*outStream);
    pObj->checkHessVec(*x,*d,true,*outStream);
    /**********************************************************************************************/
    /************************* SUPER QUANTILE QUADRANGLE ******************************************/
    /**********************************************************************************************/
    RealT val(0);
    diff->zero(); xp->zero();
    std::vector<RealT> error(20), norm(20), obj(20), objErr(20);
    *outStream << "\nSINGLETON KUSUOKA RISK MEASURE\n";
    list.sublist("SOL").set("Stochastic Optimization Type","Risk Averse"); 
    list.sublist("SOL").sublist("Risk Measure").set("Name","Singleton Kusuoka");
    int nQuadLo = 0;
    int nQuadUp = 21;
    for (int i = nQuadLo; i < nQuadUp; ++i) {
      int order = i+1; //std::pow(2,i+1);
      list.sublist("SOL").sublist("Risk Measure").sublist("Singleton Kusuoka").set("Number of Quadrature Points",order);
      setRandomVector(*x_rcp);
      obj[i] = setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
      printSolution(*x_rcp,*outStream);
      diff->set(*xp); diff->axpy(static_cast<RealT>(-1.0),*x);
      error[i] = diff->norm();
      norm[i] = x->norm();
      objErr[i] = std::abs(val-obj[i]);
      val = obj[i];
      xp->set(*x);
    }
    *outStream << std::right
               << std::setw(20) << "Num quad"
               << std::setw(20) << "norm x"
               << std::setw(20) << "norm diff"
               << std::setw(20) << "obj val"
               << std::setw(20) << "obj diff"
               << std::endl;
    for (int i = nQuadLo; i < nQuadUp; ++i) {
      *outStream << std::fixed << std::setprecision(0) << std::right
                 << std::setw(20) << static_cast<RealT>(i+1)
                 << std::scientific << std::setprecision(11) << std::right
                 << std::setw(20) << norm[i]
                 << std::setw(20) << error[i]
                 << std::setw(20) << obj[i]
                 << std::setw(20) << objErr[i]
                 << std::endl;
    }
    errorFlag += ((objErr[19] > static_cast<RealT>(1.e-3)) ? 1 : 0);
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
