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
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


template<class Container1, class Container2>
bool compareContainers(
  const Container1 &c1, const std::string &c1_name,
  const Container2 &c2, const std::string &c2_name,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::as;
  bool success = true;

  out << "\nComparing " << c1_name << " == " << c2_name << " ... ";

  const int n = c1.size();

  // Compare sizes
  if (as<int>(c2.size()) != n) {
    out << "\nError, "<<c1_name<<".size() = "<<c1.size()<<" == " 
        << c2_name<<".size() = "<<c2.size()<<" : failed!\n";
    return false;
  }
  
  // Compare elements
  for( int i = 0; i < n; ++i ) {
    const bool result = ( c1[i] == c2[i] );
    if (!result) {
      out << "\nError, "<<c1_name<<"["<<i<<"] = "<<c1[i]<<" == "
          << c2_name<<"["<<i<<"] = "<<c2[i]<<": failed!\n";
      success = false;
    }
  }
  if (success) {
    out << "passed\n";
  }

  return success;

}


template<class T>
bool testArray( const int n, Teuchos::FancyOStream &out )
{
  
  using Teuchos::Array;
  using Teuchos::setToNull;
  using Teuchos::outArg;
  using Teuchos::getConst;
  using Teuchos::NullIteratorTraits;
  using Teuchos::TypeNameTraits;

  bool success = true, result;
 
  out << "\nTesting "<<TypeNameTraits<Array<T> >::name()<<" of size = "<<n<<"\n";
  
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

  {
    
    out << "\nTest nonconst forward iterator access ... ";

    bool local_success = true;

    typedef typename Array<T>::iterator iter_t;

    iter_t iter = a.begin();

    for ( int i = 0; i < n; ++i, ++iter ) {
      const T t_i = a[i];
      result = ( t_i == *iter );
      if (!result) {
        out << "\nError, a["<<i<<"] = " << a[i] << " == *iter = " << (*iter) << ": failed!";
      }
      if (!result) local_success = false;
    }

    iter = NullIteratorTraits<iter_t>::getNull();

    if (local_success) {
      out << "passed\n";
    }
    else {
      success = false;
    }
    
  }

  {
    
    out << "\nTest const forward iterator access ... ";

    bool local_success = true;

    typedef typename Array<T>::const_iterator iter_t;

    iter_t iter = getConst(a).begin();

    for ( int i = 0; i < n; ++i, ++iter ) {
      const T t_i = a[i];
      result = ( t_i == *iter );
      if (!result) {
        out << "\nError, a["<<i<<"] = " << a[i] << " == *iter = " << (*iter) << ": failed!";
      }
      if (!result) local_success = false;
    }

    iter = NullIteratorTraits<iter_t>::getNull();

    if (local_success) {
      out << "passed\n";
    }
    else {
      success = false;
    }
    
  }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest that *a2.begin() throws for empty a2 ... ";
    try {
      Array<T> a2;
      T val = *a2.begin();
      val=0;
      success = false;
      out << "failed\n"; 
    }
    catch (const Teuchos::NullReferenceError& except) {
      out << "passed\n"; 
    }
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest that a2.begin() == a2.end() for empty a2 ... ";
    Array<T> a2;
    result = ( a2.begin() == a2.end() );
    if (result) {
      out << "passed\n";
    }
    else {
      out << "failed\n";
      success = false;
    }
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest that *(a.begin()-1) throws ... ";
    try {
      T val = *(a.begin()-1);
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
    out << "\nTest that *a.end() throws ... ";
    try {
      T val = *a.end();
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
    out << "\nTest that *(getConst(a).begin()-1) throws ... ";
    try {
      T val = *(getConst(a).begin()-1);
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
    out << "\nTest that *getConst(a).end() throws ... ";
    try {
      T val = *getConst(a).end();
      val=0;
      success = false;
      out << "failed\n"; 
    }
    catch (const Teuchos::RangeError& except) {
      out << "passed\n"; 
    }
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  {
    
    out << "\nTest nonconst reverse iterator access ... ";

    bool local_success = true;

    typedef typename Array<T>::reverse_iterator iter_t;

    iter_t iter = a.rbegin();

    for ( int i = 0; i < n; ++i, ++iter ) {
      const T t_i = a[n-i-1];
      result = ( t_i == *iter );
      if (!result) {
        out << "\nError, a["<<n-i-1<<"] = " << t_i << " == *iter = " << (*iter) << ": failed!";
      }
      if (!result) local_success = false;
    }

    iter = NullIteratorTraits<iter_t>::getNull();

    if (local_success) {
      out << "passed\n";
    }
    else {
      success = false;
    }
    
  }

  {
    
    out << "\nTest const reverse iterator access ... ";

    bool local_success = true;

    typedef typename Array<T>::const_reverse_iterator iter_t;

    iter_t iter = getConst(a).rbegin();

    for ( int i = 0; i < n; ++i, ++iter ) {
      const T t_i = a[n-i-1];
      result = ( t_i == *iter );
      if (!result) {
        out << "\nError, a["<<n-i-1<<"] = " << t_i << " == *iter = " << (*iter) << ": failed!";
      }
      if (!result) local_success = false;
    }

    iter = NullIteratorTraits<iter_t>::getNull();

    if (local_success) {
      out << "passed\n";
    }
    else {
      success = false;
    }
    
  }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest that *a2.rbegin() throws for empty a2 ... ";
    try {
      Array<T> a2;
      T val = *a2.rbegin();
      val=0;
      success = false;
      out << "failed\n"; 
    }
    catch (const Teuchos::NullReferenceError& except) {
      out << "passed\n"; 
    }
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest that a2.rbegin() == a2.rend() for empty a2 ... ";
    Array<T> a2;
    result = ( a2.rbegin() == a2.rend() );
    if (result) {
      out << "passed\n";
    }
    else {
      out << "failed\n";
      success = false;
    }
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest that *(a.rbegin()-1) throws ... ";
    try {
      T val = *(a.rbegin()-1);
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
    out << "\nTest that *a.rend() throws ... ";
    try {
      T val = *a.rend();
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
    out << "\nTest that *(getConst(a).rbegin()-1) throws ... ";
    try {
      T val = *(getConst(a).rbegin()-1);
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
    out << "\nTest that *getConst(a).rend() throws ... ";
    try {
      T val = *getConst(a).rend();
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
    out << "\nTest that an iterator reference set to null does not throw ... ";
    try {
      typedef typename Array<T>::iterator iter_t;
      iter_t iter = NullIteratorTraits<iter_t>::getNull();
      {
        Array<T> a2(n);
        iter = a2.begin();
        setToNull(outArg(iter));
      }
      out << "passed\n"; 
    }
    catch (const Teuchos::DanglingReferenceError& except) {
      success = false;
      out << "failed\n"; 
    }
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest that a dangling iterator reference throws exception ... ";
    try {
      typedef typename Array<T>::iterator iter_t;
      iter_t iter = NullIteratorTraits<iter_t>::getNull();
      {
        Array<T> a2(n);
        iter = a2.begin();
        // We did not set to null before a2 was deleted!
      }
      success = false;
      out << "failed\n"; 
    }
    catch (const Teuchos::DanglingReferenceError& except) {
      out << "passed\n"; 
    }
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  {
    out << "\nTest copy conversion to and from Teuchos::Array and std::vector ...\n";
    std::vector<T> v2 = a.toVector();
    Array<T> a2(v2);
    result = compareContainers(a2,"a2",a,"a",out);
    if (!result) success = false;
  }

  {
    out << "\nTest assignment operator taking an std::vector ...\n";
    std::vector<T> v2 = a.toVector();
    Array<T> a2;
    a2 = v2;
    result = compareContainers(a2,"a2",a,"a",out);
    if (!result) success = false;
  }

  out << "\nToDo: Test insert and delete operations!\n";
  
  out << "\nToDo: Test dangling iterators when doing insert and delete operations!\n";


  out << "\nToDo: Implement more tests!\n";

  return success;

}


int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;
  using Teuchos::Array;
	
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

    int n = 10;
    clp.setOption( "n", &n, "Number of elements in the array" );

		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

		if ( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
			*out << "\nEnd Result: TEST FAILED" << std::endl;
			return parse_return;
		}

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl;
 
    result = testArray<int>(n,*out);
    if (!result) success = false;

    result = testArray<float>(n,*out);
    if (!result) success = false;

    result = testArray<double>(n,*out);
    if (!result) success = false;
 
    // ToDo: Fill in the rest of the types!
 
	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
 
  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;
 
  return ( success ? 0 : 1 );
 
}
