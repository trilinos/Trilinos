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

/*! \file  example_03.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
           full-space methods.
*/

#include "ROL_MonteCarloGenerator.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "test_05.hpp"

typedef double RealT;

typedef H1VectorPrimal<RealT> PrimalStateVector;
typedef H1VectorDual<RealT> DualStateVector;
typedef L2VectorPrimal<RealT> PrimalControlVector;
typedef L2VectorDual<RealT> DualControlVector;
typedef H1VectorDual<RealT> PrimalConstraintVector;
typedef H1VectorPrimal<RealT> DualConstraintVector;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm
    = Teuchos::DefaultComm<int>::getComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  bool print = (iprint>0);
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (print)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.

  try {
    /*************************************************************************/
    /************* INITIALIZE BURGERS FEM CLASS ******************************/
    /*************************************************************************/
    int nx      = 512;   // Set spatial discretization.
    RealT nl    = 1.0;   // Nonlinearity parameter (1 = Burgers, 0 = linear).
    RealT cH1   = 1.0;   // Scale for derivative term in H1 norm.
    RealT cL2   = 0.0;   // Scale for mass term in H1 norm.
    Teuchos::RCP<BurgersFEM<RealT> > fem
      = Teuchos::rcp(new BurgersFEM<RealT>(nx,nl,cH1,cL2));
    fem->test_inverse_mass(*outStream);
    fem->test_inverse_H1(*outStream);
    /*************************************************************************/
    /************* INITIALIZE SIMOPT EQUALITY CONSTRAINT *********************/
    /*************************************************************************/
    bool hess = true;
    Teuchos::RCP<ROL::EqualityConstraint_SimOpt<RealT> > pcon
      = Teuchos::rcp(new EqualityConstraint_BurgersControl<RealT>(fem,hess));
    /*************************************************************************/
    /************* INITIALIZE VECTOR STORAGE *********************************/
    /*************************************************************************/
    // INITIALIZE CONTROL VECTORS
    Teuchos::RCP<std::vector<RealT> > z_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx+2, 0.0) );
    Teuchos::RCP<ROL::Vector<RealT> > zp
      = Teuchos::rcp(new PrimalControlVector(z_rcp,fem));
    // INITIALIZE STATE VECTORS
    Teuchos::RCP<std::vector<RealT> > u_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<ROL::Vector<RealT> > up
      = Teuchos::rcp(new PrimalStateVector(u_rcp,fem));
    // INITIALIZE CONSTRAINT VECTORS
    Teuchos::RCP<std::vector<RealT> > c_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<ROL::Vector<RealT> > cp
      = Teuchos::rcp(new PrimalConstraintVector(c_rcp,fem));
    /*************************************************************************/
    /************* INITIALIZE SAMPLE GENERATOR *******************************/
    /*************************************************************************/
    int dim = 4, nSamp = 10000;
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<RealT> > bounds(dim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new L2VectorBatchManager<RealT,int>(comm));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(
          nSamp,bounds,bman,false,false,100));
    /*************************************************************************/
    /************* CHECK DERIVATIVES AND CONSISTENCY *************************/
    /*************************************************************************/
    RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());
    for (int i = sampler->start(); i < sampler->numMySamples(); i++) {
//      *outStream << "(" << sampler->getMyPoint(i)[0] << ", "
//                        << sampler->getMyPoint(i)[1] << ", "
//                        << sampler->getMyPoint(i)[2] << ", "
//                        << sampler->getMyPoint(i)[3] << ")\n";
      pcon->setParameter(sampler->getMyPoint(i));
      pcon->solve(*cp,*up,*zp,tol);
      RealT rnorm = cp->norm();
      pcon->value(*cp,*up,*zp,tol);
      RealT cnorm = cp->norm();
      errorFlag += ((cnorm > tol) ? 1 : 0) + ((rnorm > tol) ? 1 : 0);
      *outStream << "Sample " << i << "  Rank " << sampler->batchID() << "\n";
      *outStream << std::scientific << std::setprecision(8);
      *outStream << "Test SimOpt solve at feasible (u,z):\n";
      *outStream << "  Solver Residual = " << rnorm << "\n";
      *outStream << "       ||c(u,z)|| = " << cnorm;
      *outStream << "\n\n";
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  comm->barrier();
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
