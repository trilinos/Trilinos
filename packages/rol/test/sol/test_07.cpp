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
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_Types.hpp"
#include "ROL_Algorithm.hpp"

#include "ROL_Objective.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"

#include "ROL_StochasticProblem.hpp"

typedef double RealT;

template<class Real> 
class ParametrizedObjectiveEx7 : public ROL::Objective<Real> {
public:
  Real value( const ROL::Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > ex
      = Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();
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
    Teuchos::RCP<const std::vector<Real> > ex
      = Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();
    Teuchos::RCP<std::vector<Real> > eg
      = Teuchos::dyn_cast<ROL::StdVector<Real> >(g).getVector();
    std::vector<Real> p = ROL::Objective<Real>::getParameter();
    unsigned size = static_cast<unsigned>(ex->size());
    const Real two(2);
    for ( unsigned i = 0; i < size; i++ ) {
      (*eg)[i] = two*std::exp(p[0])*(*ex)[i] + p[i+1];
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > ev
      = Teuchos::dyn_cast<const ROL::StdVector<Real> >(v).getVector();
    Teuchos::RCP<std::vector<Real> > ehv
      = Teuchos::dyn_cast<ROL::StdVector<Real> >(hv).getVector();
    std::vector<Real> p = ROL::Objective<Real>::getParameter();
    unsigned size = static_cast<unsigned>(ev->size());
    const Real two(2);
    for ( unsigned i = 0; i < size; i++ ) {
      (*ehv)[i] = two*std::exp(p[0])*(*ev)[i]; 
    } 
  }
};

void setUpAndSolve(Teuchos::ParameterList &list,
                   Teuchos::RCP<ROL::Objective<RealT> > &pObj,
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
}

template<class Real>
Real random(const Teuchos::RCP<const Teuchos::Comm<int> > &commptr) {
  Real val(0);
  if ( Teuchos::rank<int>(*commptr)==0 ) {
    //srand(time(NULL));
    srand(10);
    val = (Real)rand()/(Real)RAND_MAX;
  }
  Teuchos::broadcast<int,Real>(*commptr,0,&val);
  return val;
}

void setRandomVector(std::vector<RealT> &x,
               const Teuchos::RCP<const Teuchos::Comm<int> > &commptr) {
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
  Teuchos::RCP<const Teuchos::Comm<int> > commptr =
    Teuchos::DefaultComm<int>::getComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0 && commptr->getRank()==0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input_07.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    Teuchos::ParameterList list = *parlist;
    // Build ROL algorithm
    Teuchos::RCP<ROL::Algorithm<RealT> > algo;
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    const unsigned dim = 4;
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp(new std::vector<RealT>(dim));
    Teuchos::RCP<std::vector<RealT> > d_rcp = Teuchos::rcp(new std::vector<RealT>(dim));
    Teuchos::RCP<ROL::Vector<RealT> > x = Teuchos::rcp(new ROL::StdVector<RealT>(x_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > d = Teuchos::rcp(new ROL::StdVector<RealT>(d_rcp));
    setRandomVector(*x_rcp,commptr);
    setRandomVector(*d_rcp,commptr);
    // Build samplers
    const RealT zero(0), one(1);
    const int nSamp = 1000;
    const unsigned sdim = dim + 2;
    std::vector<RealT> tmp = {-one, one};
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman =
      Teuchos::rcp(new ROL::StdTeuchosBatchManager<RealT,int>(commptr));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler =
      Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nSamp,bounds,bman,false,false,100));
    // Build risk-averse objective function
    Teuchos::RCP<ROL::Objective<RealT> > pObj =
      Teuchos::rcp(new ParametrizedObjectiveEx7<RealT>);
    // Build bound constraints
    std::vector<RealT> l(dim,zero);
    std::vector<RealT> u(dim,one);
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd = 
      Teuchos::rcp( new ROL::StdBoundConstraint<RealT>(l,u) );
    bnd->deactivate();
    // Test parametrized objective functions
    *outStream << "Check Derivatives of Parametrized Objective Function\n";
    pObj->setParameter(sampler->getMyPoint(0));
    pObj->checkGradient(*x,*d,true,*outStream);
    pObj->checkHessVec(*x,*d,true,*outStream);
    /**********************************************************************************************/
    /************************* MEAN PLUS HMCR *****************************************************/
    /**********************************************************************************************/
    *outStream << "\nMEAN PLUS HIGHER MOMENT COHERENT RISK MEASURE\n";
    setRandomVector(*x_rcp,commptr);
    setUpAndSolve(list,pObj,sampler,x,d,bnd,*outStream);
    printSolution(*x_rcp,*outStream);
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
