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
//         
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_08.cpp
    \brief Interior Point test using Hock & Schittkowski problem 29.
*/

#include "ROL_HS29.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_BoundFletcher.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_OptimizationSolver.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"


typedef double RealT;

int main(int argc, char *argv[]) {


  typedef std::vector<RealT>            vec;
  typedef ROL::StdVector<RealT>         SV;
  typedef ROL::Ptr<ROL::Vector<RealT> > PtrV;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  std::string filename = "input_ex02.xml";
  
  auto parlist = ROL::getParametersFromXmlFile( filename );

  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {

    int xopt_dim  = 3; // Dimension of optimization vectors
    int ci_dim    = 1; // Dimension of inequality constraint

    ROL::Ptr<vec> xopt_ptr = ROL::makePtr<vec>(xopt_dim,1.0); // Feasible initial guess
    ROL::Ptr<vec> li_ptr  = ROL::makePtr<vec>(ci_dim,0.0);

    PtrV xopt = ROL::makePtr<SV>(xopt_ptr);
    PtrV li   = ROL::makePtr<SV>(li_ptr);

    // Original obective
    using ROL::ZOO::Objective_HS29;
    using ROL::ZOO::InequalityConstraint_HS29;
    
    ROL::Ptr<ROL::Objective<RealT> >  obj_hs29 = ROL::makePtr<Objective_HS29<RealT>>();
    ROL::Ptr<ROL::Constraint<RealT> >  incon_hs29 = ROL::makePtr<InequalityConstraint_HS29<RealT>>();

    ROL::Ptr<ROL::Vector<RealT> > bndc = li->clone();
    ROL::Ptr<ROL::BoundConstraint<RealT> > bndcon = ROL::makePtr<ROL::Bounds<RealT> >(*bndc);

    ROL::Ptr<vec> low_ptr = ROL::makePtr<vec>(xopt_dim, ROL::ROL_NINF<RealT>());
    ROL::Ptr<vec> upp_ptr = ROL::makePtr<vec>(xopt_dim, ROL::ROL_INF<RealT>());
    PtrV low = ROL::makePtr<SV>(low_ptr);
    PtrV upp = ROL::makePtr<SV>(upp_ptr);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bndx = ROL::makePtr<ROL::Bounds<RealT> >(low, upp);

    std::string stepname = "Trust Region";

    ROL::OptimizationProblem<RealT> problem( obj_hs29, xopt, bndx, incon_hs29, li, bndcon );  

    ROL::Ptr<ROL::Objective<RealT> > obj = problem.getObjective();    
    ROL::Ptr<ROL::Constraint<RealT> > con = problem.getConstraint();
    ROL::Ptr<ROL::BoundConstraint<RealT> > bndxs = problem.getBoundConstraint();
    ROL::Ptr<ROL::Vector<RealT> > xs = problem.getSolutionVector();
    ROL::Ptr<ROL::Vector<RealT> > lis = problem.getMultiplierVector();

    ROL::Ptr<ROL::BoundFletcher<RealT> > fletcher_penalty = ROL::makePtr<ROL::BoundFletcher<RealT> >(obj, con, bndxs, *xs, *lis, *parlist);

    ROL::Ptr<ROL::Vector<RealT> > v = xs->clone(); v->randomize();
    std::vector<std::vector<RealT> > gCheck = fletcher_penalty->checkGradient(*xs, *v, true );

    ROL::Ptr<ROL::Vector<RealT> > w = xs->clone(); w->randomize();
    std::vector<RealT> hCheck = fletcher_penalty->checkHessSym( *xs, *v, *w, true, *outStream);

    // Define algorithm.
    // ROL::Ptr<ROL::Algorithm<RealT> > algo;
    // algo = ROL::makePtr<ROL::Algorithm<RealT>>(stepname, *parlist);
    // algo->run(*xs, *fletcher_penalty, *bndxs, true, *outStream);   

    ROL::OptimizationSolver<RealT> optSolver(problem, *parlist);
    optSolver.solve(*outStream);    

    *outStream << std::endl << std::setw(20) << "Computed Minimizer" << std::endl;
    for( int i=0;i<xopt_dim;++i ) {   
      *outStream << std::setw(20) << (*xopt_ptr)[i] << std::endl;
    }

    *outStream << "Exact minimizers: x* = (a,b,c), (a,-b,-c), (-a,b,-c), (-a,-b,c)" << std::endl;
    *outStream << "Where a=4, b=" << 2*std::sqrt(2) << ", and c=2" << std::endl;

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
