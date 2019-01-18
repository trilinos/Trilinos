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

#include "test_01.hpp"
#include "ROL_StdTransform_PEBBL.hpp"
#include "ROL_TransformedConstraint_PEBBL.hpp"
#include "ROL_TransformedObjective_PEBBL.hpp"

typedef double RealT;

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
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
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int N = 10;
    ROL::Ptr<std::vector<RealT>> x_ptr    = ROL::makePtr<std::vector<RealT>>(N,0.0);
    ROL::Ptr<ROL::Vector<RealT>> x        = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    ROL::Ptr<std::vector<RealT>> emul_ptr = ROL::makePtr<std::vector<RealT>>(1,0.0);
    ROL::Ptr<ROL::Vector<RealT>> emul     = ROL::makePtr<ROL::StdVector<RealT>>(emul_ptr);
    ROL::Ptr<std::vector<RealT>> xbin_ptr = ROL::makePtr<std::vector<RealT>>(N,0.0);
    ROL::Ptr<ROL::Vector<RealT>> xbin     = ROL::makePtr<ROL::StdVector<RealT>>(xbin_ptr);
    /**********************************************************************************************/
    /************************* CONSTRUCT OBJECTIVE FUNCTION ***************************************/
    /**********************************************************************************************/
    std::vector<RealT> alpha(N,0.0);
    *outStream << std::endl << std::endl;
    *outStream << "alpha =";
    for (int i = 0; i < N; ++i) {
      alpha[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      *outStream << "  " << alpha[i];
    }
    *outStream << std::endl << std::endl;
    ROL::Ptr<ROL::Objective<RealT> > obj
      = ROL::makePtr<Objective_SimpleBinary<RealT>>(alpha);
    /**********************************************************************************************/
    /************************* CONSTRUCT CONSTRAINT ***********************************************/
    /**********************************************************************************************/
    int budget = 3;
    ROL::Ptr<ROL::Constraint<RealT> > econ
      = ROL::makePtr<Constraint_SimpleBinary<RealT>>(budget);
    /**********************************************************************************************/
    /************************* CONSTRUCT BOUND CONSTRAINT *****************************************/
    /**********************************************************************************************/
    ROL::Ptr<std::vector<RealT>> xl_ptr = ROL::makePtr<std::vector<RealT>>(N,0.0);
    ROL::Ptr<ROL::Vector<RealT>> xl     = ROL::makePtr<ROL::StdVector<RealT>>(xl_ptr);
    ROL::Ptr<std::vector<RealT>> xu_ptr = ROL::makePtr<std::vector<RealT>>(N,1.0);
    ROL::Ptr<ROL::Vector<RealT>> xu     = ROL::makePtr<ROL::StdVector<RealT>>(xu_ptr);
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(xl,xu);
    /**********************************************************************************************/
    /************************* SOLVE **************************************************************/
    /**********************************************************************************************/
    ROL::OptimizationProblem<RealT> problem(obj,x,bnd,econ,emul);
    problem.check(*outStream);
    ROL::OptimizationSolver<RealT> solver(problem,*parlist);
    *outStream << "Solve problem with no fixed binary variables"
               << std::endl << std::endl;
    clock_t start = clock();
    solver.solve(*outStream);
    *outStream << "Optimization time: "
               << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC
               << " seconds." << std::endl;
    *outStream << std::endl;
    RealT sum(0);
    *outStream << "x =";
    for (int i = 0; i < N; ++i) {
      *outStream << "  " << (*x_ptr)[i];
      sum += (*x_ptr)[i];
    }
    *outStream << std::endl << std::endl;
    *outStream << "Sum(x) = " << sum << "  Budget = " << budget;
    *outStream << std::endl << std::endl;
    /**********************************************************************************************/
    /************************* ADD BINARY CONSTRAINTS AND SOLVE ***********************************/
    /**********************************************************************************************/
    ROL::Ptr<ROL::Transform_PEBBL<RealT>> trans_bin
      = ROL::makePtr<ROL::StdTransform_PEBBL<RealT>>();
    std::map<int,RealT> fixed;
    fixed.insert(std::pair<int,RealT>(2,0.0));
    fixed.insert(std::pair<int,RealT>(5,0.0));
    fixed.insert(std::pair<int,RealT>(3,1.0));
    fixed.insert(std::pair<int,RealT>(9,1.0));
    trans_bin->add(fixed);
    RealT tol(1e-8);
    trans_bin->value(*xbin,*x,tol);
    ROL::Ptr<ROL::TransformedConstraint_PEBBL<RealT>> econ_trans
      = ROL::makePtr<ROL::TransformedConstraint_PEBBL<RealT>>(econ,trans_bin);
    ROL::Ptr<ROL::TransformedObjective_PEBBL<RealT>> obj_trans
      = ROL::makePtr<ROL::TransformedObjective_PEBBL<RealT>>(obj,trans_bin);
    ROL::OptimizationProblem<RealT> problem_bin(obj_trans,xbin,bnd,econ_trans,emul);
    problem_bin.check(*outStream);
    ROL::OptimizationSolver<RealT> solver_bin(problem_bin,*parlist);
    *outStream << "Solve problem with {2,5} set to 0 and {3,9} set to 1"
               << std::endl << std::endl;
    clock_t start_bin = clock();
    solver_bin.solve(*outStream);
    *outStream << "Optimization time: "
               << (RealT)(clock()-start_bin)/(RealT)CLOCKS_PER_SEC
               << " seconds." << std::endl;
    *outStream << std::endl;
    RealT sum_bin(0);
    *outStream << "x =";
    for (int i = 0; i < N; ++i) {
      *outStream << "  " << (*xbin_ptr)[i];
      sum_bin += (*xbin_ptr)[i];
    }
    *outStream << std::endl << std::endl;
    *outStream << "Sum(x) = " << sum_bin << "  Budget = " << budget;
    *outStream << std::endl << std::endl;

    errorFlag += ((*xbin_ptr)[2]==0.0 ? 0 : 1);
    errorFlag += ((*xbin_ptr)[5]==0.0 ? 0 : 1);
    errorFlag += ((*xbin_ptr)[3]==1.0 ? 0 : 1);
    errorFlag += ((*xbin_ptr)[9]==1.0 ? 0 : 1);
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
