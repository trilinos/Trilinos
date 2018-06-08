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
    \brief Shows how to use Tempus.
*/

#include "Thyra_DefaultSpmdVector.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include <fstream>

#include "example_01.hpp"

typedef double RealT;


int main(int argc, char *argv[]) {

  using Teuchos::RCP;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  try {
    // Play with Thyra vector.
    std::vector<RealT> datavec = {1, 2, 3};
    auto dim = datavec.size();
    Teuchos::ArrayRCP<RealT> datavec_arcp(&datavec[0], 0, dim, false);
    Teuchos::RCP<Thyra::DefaultSpmdVectorSpace<RealT>> vecspace =
      Thyra::defaultSpmdVectorSpace<RealT>(dim);
    Thyra::DefaultSpmdVector<RealT> vec(vecspace, datavec_arcp, 1);

    auto vecnorm = vec.norm_2();
    *outStream << "\nthyra vector norm = " << vecnorm << "\n";

    // Play with Tempus.
    RCP<Teuchos::ParameterList> parList = Teuchos::rcp( new Teuchos::ParameterList );
    Teuchos::updateParametersFromXmlFile("example_01.xml", parList.ptr());
    RCP<Teuchos::ParameterList> tempusParList = sublist(parList, "Tempus", true);
    RCP<SinCosModelEvaluator<RealT>> model = Teuchos::rcp(new SinCosModelEvaluator<RealT>());

    RCP<Tempus::IntegratorBasic<RealT> > integrator = Tempus::integratorBasic<RealT>(tempusParList, model);

    integrator->advanceTime();

    // Test if at 'Final Time'
//    RealT time = integrator->getTime(); // Unused
//    RealT timeFinal = tempusParList->sublist("Demo Integrator").sublist("Time Step Control").get<RealT>("Final Time"); // Unused

    // Output solution.
    std::ofstream ftmp("Tempus_ForwardEuler_SinCos.dat");
    RCP<const Tempus::SolutionHistory<double> > solutionHistory = integrator->getSolutionHistory();
    for (int i=0; i<solutionHistory->getNumStates(); i++) {
      RCP<const Tempus::SolutionState<double> > solutionState = (*solutionHistory)[i];
      double time = solutionState->getTime();
      RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
      ftmp << time << "   " << Thyra::get_ele(*(x_plot), 0) << "   " << Thyra::get_ele(*(x_plot), 1) << std::endl;
    }
    ftmp.close();

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

