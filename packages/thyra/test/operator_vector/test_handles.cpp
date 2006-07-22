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
#include "Thyra_VectorOpTester.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"

using namespace Teuchos;
using namespace Thyra;

int main( int argc, char *argv[] ) 
{
  
  bool success = false;
  
  GlobalMPISession mpiSession(&argc, &argv);
  typedef Teuchos::ScalarTraits<double> ST;
    
  // Get stream that can print to just root or all streams!
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try {
    
    int  n = 100;
    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "local-dim", &n, "Local number of elements in each constituent vector." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;
    
    RefCountPtr<const VectorSpaceBase<double> >
      sp = rcp(
        new DefaultSpmdVectorSpace<double>(
          Teuchos::DefaultComm<Index>::getComm(),n,-1
          )
        );

    /* ------- test on a monolithic space ------------ */

    *out << "======= Testing on a monolithic vector space ======" << std::endl;
    VectorSpace<double> space = sp;
    
    VectorOpTester<double> tester(space, TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

    tester.runAllTests();

    /* -------- test on a block space ----------------*/
    *out << "======= Testing on a block vector space ======" << std::endl;
    VectorSpace<double> blockSpace = productSpace(space, space);
    
    tester = VectorOpTester<double>(blockSpace, 
                                    TestSpecifier<double>(true, 1.0e-13, 1.0e-10));
    
    tester.runAllTests();           


    

    /* -------- test on a space with recursive block structure -----------*/
    *out << "======= Testing on a recursively blocked vector space ======" << std::endl;
    VectorSpace<double> recSpace = productSpace(space, blockSpace);
    
    tester = VectorOpTester<double>(recSpace, 
                                    TestSpecifier<double>(true, 1.0e-13, 1.0e-10));
    
    tester.runAllTests();           


    


    
    
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,success)

}

