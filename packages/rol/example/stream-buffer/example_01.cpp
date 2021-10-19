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

/*! \file  example_01.cpp
    \brief Demonstrates how to inherit from std::streambuf to output to screen
           and to a stringstream.
*/

#include "ROL_HS41.hpp"
#include "ROL_HS53.hpp"
#include "ROL_TypeB_TrustRegionSPGAlgorithm.hpp"
#include "ROL_Problem.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
#include <sstream>

typedef double RealT;

template<class charT=char, class Traits=std::char_traits<charT>>
class mybuffer : public std::basic_streambuf<charT,Traits> {
public:
  const std::stringstream& getStringStream() const { return ss; }

protected:
  inline virtual int overflow(int c = Traits::eof()) {
    if (c != Traits::eof())          ss << static_cast<charT>(c);
    if (putchar(c) == Traits::eof()) return Traits::eof();
    return c;
  }

private:
  std::stringstream ss;
};

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a
  // (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  mybuffer<> buf;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtr<std::ostream>(&buf);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {
    //RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());

    ROL::ParameterList list;
    list.sublist("Status Test").set("Gradient Tolerance",1e-8);
    list.sublist("Status Test").set("Constraint Tolerance",1e-8);
    list.sublist("Status Test").set("Step Tolerance",1e-12);
    list.sublist("Status Test").set("Iteration Limit", 250);
    list.sublist("General").set("Output Level",iprint);
    list.sublist("General").sublist("Polyhedral Projection").set("Type","Semismooth Newton");
    list.sublist("General").sublist("Polyhedral Projection").set("Iteration Limit",5000);
    list.sublist("General").sublist("Secant").set("Type","Limited-Memory BFGS");

    ROL::Ptr<ROL::Vector<RealT>>     sol, mul;
    ROL::Ptr<ROL::Objective<RealT>>  obj;
    ROL::Ptr<ROL::Constraint<RealT>> con;
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd;
    ROL::Ptr<ROL::Problem<RealT>> problem;
    ROL::Ptr<ROL::TypeB::TrustRegionSPGAlgorithm<RealT>> algo;

    *outStream << std::endl << "Hock and Schittkowski Problem #41" << std::endl << std::endl;
    ROL::ZOO::getHS41<RealT> HS41;
    obj = HS41.getObjective();
    sol = HS41.getInitialGuess();
    con = HS41.getEqualityConstraint();
    mul = HS41.getEqualityMultiplier();
    bnd = HS41.getBoundConstraint();

    if (mul->dimension() == 1)
      list.sublist("General").sublist("Polyhedral Projection").set("Type","Dai-Fletcher");
    else
      list.sublist("General").sublist("Polyhedral Projection").set("Type","Semismooth Newton");
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,sol);
    problem->addBoundConstraint(bnd);
    problem->addLinearConstraint("LEC",con,mul);
    problem->setProjectionAlgorithm(list);
    algo = ROL::makePtr<ROL::TypeB::TrustRegionSPGAlgorithm<RealT>>(list);
    algo->run(*problem,*outStream);

    if (iprint > 0)
      std::cout << std::endl << "Output from stored stringstream" << std::endl
                << buf.getStringStream().str() << std::endl;
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

