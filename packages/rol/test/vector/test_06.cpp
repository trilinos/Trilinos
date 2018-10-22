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

/*! \file  test_06.cpp
    \brief Test ProfiledVector interface.
*/

#include "ROL_ProfiledVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Zakharov.hpp"
#include "ROL_Algorithm.hpp"



#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

typedef int    OrdinalT;
typedef double RealT;

template<>
ROL::VectorFunctionCalls<int>
ROL::ProfiledVector<int,RealT>::functionCalls_ = ROL::VectorFunctionCalls<int>();

int main(int argc, char *argv[]) {


  using ROL::ParameterList;

  typedef std::vector<RealT>                  vector;

  typedef ROL::Vector<RealT>                  V;
  typedef ROL::StdVector<RealT>               SV;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);

  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {
    // Dimension of the optimization vector
    int dim = 10;

    std::string paramfile = "parameters.xml";
    auto parlist = ROL::getParametersFromXmlFile(paramfile);

    // Define algorithm.
    ROL::Algorithm<RealT> algo("Trust-Region",*parlist);

    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(dim,1.0);
    ROL::Ptr<vector> k_ptr = ROL::makePtr<vector>(dim);

    for(int i=0;i<dim;++i) {
      (*k_ptr)[i] = 1.0 + i;
    }

    ROL::Ptr<V> xs = ROL::makePtr<SV>(x_ptr);
    ROL::Ptr<V> ks = ROL::makePtr<SV>(k_ptr);

    // Create ProfiledVector objects
    ROL::ProfiledVector<int,RealT> xpf(xs);
    ROL::Ptr<V> kpf = ROL::makePtr<ROL::ProfiledVector<int,RealT>>(ks);

    ROL::ZOO::Objective_Zakharov<RealT> obj(kpf);

    // Run algorithm.
    algo.run(xpf, obj, true, *outStream);

    // Display number of function calls
    ROL::printVectorFunctionCalls(xpf, *outStream);
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
