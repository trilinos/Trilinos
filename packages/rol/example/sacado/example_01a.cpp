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

/** \file
    \brief An example combining ROL and Sacado to mimize the Zakharov 
           function. The gradient and the action of the Hessian on a given
           vector are computed by Sacado using automatic differentiation.    
                    
           This implementation is far from optimal as vectors of AD type
           are being created repeatedly. A more efficient implementation 
           would encapsulate the functions Zakharov, objgrad, and hessVec
           in an object so that the AD vectors can be instantiated once. 

    \author Created by G. von Winckel
**/

#include <iostream>

#include "ROL_Sacado_Objective.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "example_01a.hpp"

using namespace ROL;

typedef double RealT;

int main(int argc, char **argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
    int iprint     = argc - 1;
    ROL::Ptr<std::ostream> outStream;
    ROL::nullstream bhs; // outputs nothing
    if (iprint > 0)
        outStream = ROL::makePtrFromRef(std::cout);
    else
        outStream = ROL::makePtrFromRef(bhs);

    int errorFlag  = 0;

    // *** Example body.

    try {

        Sacado_Objective<RealT,Zakharov> obj;
    
        int dim = 10; // Set problem dimension. 

        // Load optimizer parameters form XML file
        std::string paramfile = "parameters.xml";
        auto parlist = ROL::getParametersFromXmlFile(paramfile);

        // Define algorithm.
        Algorithm<RealT> algo("Line Search",*parlist);

        // Iteration vector.
        ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
        // Set Initial Guess
        for (int i=0; i<dim; i++) {
            (*x_ptr)[i]   = 2;
        }

        StdVector<RealT> x(x_ptr);

        // Run algorithm.
        algo.run(x, obj, true, *outStream);

        // Get true solution.
        ROL::Ptr<std::vector<RealT> > xtrue_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
        StdVector<RealT> xtrue(xtrue_ptr);
        
        // Compute error.
        x.axpy(-1.0, xtrue);
        RealT abserr = x.norm();
        *outStream << std::scientific << "\n   Absolute Error: " << abserr << std::endl;
        if ( abserr > sqrt(ROL_EPSILON<RealT>()) ) {
            errorFlag += 1;
        }
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
