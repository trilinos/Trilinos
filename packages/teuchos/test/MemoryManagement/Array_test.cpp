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

#include "Teuchos_Array.new.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"


namespace Teuchos {

/** \brief Base traits class for getting a properly initialized null pointer.
 *
 * This default traits class simply
 */
template<typename Ptr>
class NullPtr {
public:
  static Ptr get() { return Ptr(0); }
};


template<typename T>
class NullPtr<RCP<T> > {
public:
  static RCP<T> get() { return null; }
};


template<typename T>
class NullPtr<ArrayRCP<T> > {
public:
  static ArrayRCP<T> get() { return null; }
};


template<typename Ptr>
void setToNull( Ptr *ptr ) {
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==ptr);
#endif
  *ptr = NullPtr<Ptr>::get();
}


} // namespace Teuchos


template<class T>
bool test_Array( const int n, Teuchos::FancyOStream &out )
{
  
  using Teuchos::Array;

  bool success = false, result;
 
  out << "\nTesting Array<"<<Teuchos::TypeNameTraits<T>::name()<<"> of size = "<<n<<"\n";
  
  Teuchos::OSTab tab(out);

  Array<T> a(n);
 
  {
    out << "\nInitializing data ...\n";
    for( int i = 0; i < n; ++i )
      a[i] = i;
  }

  {
    out << "\nTest that a[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      const T t_i = i;
      result = ( t_i == a[i] );
      if (!result) {
        out << "\nError, a["<<i<<"] = " << a[i] << " == i = " << i << ": failed!";
      }
      if (!result) local_success = false;
    }
    if (local_success) {
      out << "passed\n";
    }
    else {
      success = false;
    }
  }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest that a[-1] throws ... ";
    try {
      T val = a[-1];
      val=0;
      success = false;
      out << "failed\n"; 
    }
    catch (const Teuchos::RangeError& except) {
      out << "passed\n"; 
    }
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest that a[n] throws ... ";
    try {
      T val = a[n];
      val=0;
      success = false;
      out << "failed\n"; 
    }
    catch (const Teuchos::RangeError& except) {
      out << "passed\n"; 
    }
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  typedef typename Array<T>::iterator iter_t;

  iter_t itr = Teuchos::NullPtr<iter_t>::get();

  TEST_FOR_EXCEPT("ToDo: Implement more tests!");

  return success;

}

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;
  using Teuchos::Array;
	
	bool success = true;
 
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  //const int procRank = Teuchos::GlobalMPISession::getRank();
 
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
 
	try {
    
    //
		// Read options from the commandline
    //

    CommandLineProcessor clp(false); // Don't throw exceptions

    int n = 10;
    clp.setOption( "n", &n, "Number of elements in the array" );

		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

		if ( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
			*out << "\nEnd Result: TEST FAILED" << std::endl;
			return parse_return;
		}

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl;
 
    test_Array<int>(n,*out);
    test_Array<float>(n,*out);
    test_Array<double>(n,*out);
 
    // ToDo: Fill in the rest of the types!
 
    *out << "\nAll tests for Array seem to check out!\n";
 
	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
 
  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;	
 
  return ( success ? 0 : 1 );
 
}
