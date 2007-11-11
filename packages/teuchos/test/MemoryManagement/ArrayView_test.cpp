// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#include "Teuchos_ArrayView.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_TestingHelpers.hpp"


//
// Define local macros to make defining tests easier for this particular test
// code.
//
// Note, macros with these types of names should only exist in a *.cpp file
// after all #includes are done!
//


#define TEST_EQUALITY_CONST( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY_CONST( v1, v2, out, success )

#define TEST_EQUALITY( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY( v1, v2, out, success )

#define TEST_ITER_EQUALITY( iter1, iter2 ) \
  TEUCHOS_TEST_ITER_EQUALITY( iter1, iter2, out, success )

#define TEST_ARRAY_ELE_EQUALITY( a, i, val ) \
   TEUCHOS_TEST_ARRAY_ELE_EQUALITY( a, i, val, false, out, local_success )

#define TEST_COMPARE( v1, comp, v2 ) \
  TEUCHOS_TEST_COMPARE( v1, comp, v2, out, success )

#define TEST_COMPARE_ARRAYS( a1, a2 ) \
  { \
    const bool result = compareArrays(a1,#a1,a2,#a2,out); \
    if (!result) success = false; \
  }

#define TEST_THROW( code, ExceptType  ) \
  TEUCHOS_TEST_THROW( code, ExceptType, out, success  )

#define TEST_NOTHROW( code  ) \
  TEUCHOS_TEST_NOTHROW( code, out, success  )


//
// Main templated array test function
//


template<class T>
bool testArrayView( const int n, Teuchos::FancyOStream &out )
{
  
  using Teuchos::ArrayView;
  using Teuchos::setToNull;
  using Teuchos::outArg;
  using Teuchos::NullIteratorTraits;
  using Teuchos::TypeNameTraits;
  using Teuchos::as;
  typedef typename ArrayView<T>::size_type size_type;

  bool success = true;
 
  out
    << "\n***"
    << "\n*** Testing "<<TypeNameTraits<ArrayView<T> >::name()<<" of size = "<<n
    << "\n***\n";
  
  Teuchos::OSTab tab(out);

  TEST_FOR_EXCEPT("ToDo: Implement!");

  return success;

}


//
// Main testing program
//

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;
	
	bool success = true;
  bool result;
 
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  //const int procRank = Teuchos::GlobalMPISession::getRank();
 
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
 
	try {
    
    //
		// Read options from the commandline
    //

    CommandLineProcessor clp(false); // Don't throw exceptions

    int n = 4;
    clp.setOption( "n", &n, "Number of elements in the array" );

		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

		if ( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
			*out << "\nEnd Result: TEST FAILED" << std::endl;
			return parse_return;
		}

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl;
 
    result = testArrayView<int>(n,*out);
    if (!result) success = false;

    /*

    result = testArrayView<float>(n,*out);
    if (!result) success = false;

    result = testArrayView<double>(n,*out);
    if (!result) success = false;

    result = testArrayView<std::complex<double> >(n,*out);
    if (!result) success = false;

    */
 
    // ToDo: Fill in the rest of the types!
 
	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
 
  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;
 
  return ( success ? 0 : 1 );
 
}
