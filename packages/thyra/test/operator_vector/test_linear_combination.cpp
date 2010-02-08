// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Thyra_VectorImpl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_LinearCombinationTester.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"

using namespace Teuchos;
using namespace Thyra;


#define SKIP_BLOCK_TESTS

#define TEST_FLOATS
#define TEST_COMPLEX

template <class Scalar> inline 
bool runTests(int n, Teuchos::RCP<Teuchos::FancyOStream>& out)
{
  typedef typename Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType Mag;
  typedef Thyra::Ordinal Ordinal;

  Mag epsErr = n * 1.0e2 * ST::prec();
  Mag epsWarn = 1.0e2 * epsErr;

  /* ------- test on a monolithic space ------------ */
  
  *out << "======= Testing on a monolithic vector space ======" << std::endl;
  VectorSpace<Scalar> space(
    Thyra::defaultSpmdVectorSpace<Scalar>(DefaultComm<Ordinal>::getComm(),n,-1));
  
  TestSpecifier<Scalar> spec(true, epsErr, epsWarn);
  
  LinearCombinationTester<Scalar> tester(DefaultComm<Ordinal>::getComm(),
                                         space, out, spec); 
  
  return tester.runAllTests();
}  

int main( int argc, char *argv[] ) 
{
  bool success = false;
  
  GlobalMPISession mpiSession(&argc, &argv);
    
  // Get stream that can print to just root or all streams!
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try {
    
    int  n = 20;
    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    bool verbose = false;
    clp.addOutputSetupOptions(true);
    clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
    clp.setOption( "local-dim", &n, "Local number of elements in each constituent vector." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    if (!verbose) out = rcp(new FancyOStream(rcp(new oblackholestream())));
    
    success = runTests<double>(n, out) ;

#if defined(HAVE_THYRA_FLOAT)
    success = runTests<float>(n, out) && success;
#endif

#if defined(HAVE_THYRA_COMPLEX)
    success = runTests<std::complex<double> >(n, out) && success;
#endif
  }
  
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,success)

    if (success)
      {
        *out << "all tests PASSED!" << std::endl;
        return 0;
      }
    else
      {
        *out << "at least one test FAILED!" << std::endl;
        return 1;
      }
}

