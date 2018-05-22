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

#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
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

#include "ROL_OptimizationProblem.hpp"

typedef double RealT;

template<class Real>
class ParametrizedObjectiveEx7 : public ROL::Objective<Real> {
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
  ROL::OptimizationProblem<RealT> opt(pObj,x,bnd);
  opt.setStochasticObjective(list,sampler);
  outStream << "\nCheck Derivatives of Stochastic Objective Function\n";
  opt.check(outStream);
  // Run ROL algorithm
  ROL::Algorithm<RealT> algo("Trust Region",list,false);
  algo.run(opt,true,outStream);
}

template<class Real>
Real random(const ROL::Ptr<const Teuchos::Comm<int> > &commptr) {
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
    std::string filename = "input_07.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    ROL::ParameterList list = *parlist;
    // Build ROL algorithm
    ROL::Ptr<ROL::Algorithm<RealT> > algo;
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    const unsigned dim = 4;
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim);
    ROL::Ptr<std::vector<RealT> > d_ptr = ROL::makePtr<std::vector<RealT>>(dim);
    ROL::Ptr<ROL::Vector<RealT> > x = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    ROL::Ptr<ROL::Vector<RealT> > d = ROL::makePtr<ROL::StdVector<RealT>>(d_ptr);
    setRandomVector(*x_ptr,commptr);
    setRandomVector(*d_ptr,commptr);
    // Build samplers
    const RealT zero(0), one(1);
    const int nSamp = 1000;
    const unsigned sdim = dim + 2;
    std::vector<RealT> tmp = {-one, one};
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT> > bman =
      ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(commptr);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler =
      ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nSamp,bounds,bman,false,false,100);
    // Build risk-averse objective function
    ROL::Ptr<ROL::Objective<RealT> > pObj =
      ROL::makePtr<ParametrizedObjectiveEx7<RealT>>();
    // Build bound constraints
    std::vector<RealT> l(dim,zero);
    std::vector<RealT> u(dim,one);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd = 
      ROL::makePtr<ROL::StdBoundConstraint<RealT>>(l,u);
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
    setRandomVector(*x_ptr,commptr);
    setUpAndSolve(list,pObj,sampler,x,bnd,*outStream);
    printSolution(*x_ptr,*outStream);
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
