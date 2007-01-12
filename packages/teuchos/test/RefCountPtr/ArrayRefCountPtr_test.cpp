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

#include "Teuchos_ArrayRefCountPtr.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"

// Temporarily uncomment any or all of these macros to see compilation
// failures for code that is rightfully not supposed to compile (which is a
// wonderful thing)!  The fact that this code does not compile show that the
// design of the Teuchos::ArrayRefCountPtr class supports full support of
// const projection in all of its forms when dealing with arrays of objects.
//#define SHOW_COMPILE_FAILURE_1
//#define SHOW_COMPILE_FAILURE_2
//#define SHOW_COMPILE_FAILURE_3

template<class T>
void test_ArrayRefCountPtr_iterators(
  const Teuchos::ArrayRefCountPtr<T>   &ptr
  ,const bool                          verbose
  ,Teuchos::FancyOStream               &out_arg
  )
{

  using Teuchos::ArrayRefCountPtr;
  using Teuchos::null;
  using Teuchos::arcp;
  using Teuchos::arcp_reinterpret_cast;

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::rcp(&out_arg,false);
  
  if(verbose)
    *out << "\nTesting iterators and accessors for ptr = " << ptr <<"\n";

  Teuchos::OSTab tab(out);

  const int size = ptr.size();

  // Pointer ++
  
  {
    if(verbose)
      *out << "\nChecking ++itr and < ...\n";
    ArrayRefCountPtr<T> itr = ptr;
    for( int i = 0; itr < ptr+size; ++i, ++itr )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking itr++ and <= ...\n";
    ArrayRefCountPtr<T> itr = ptr;
    for( int i = 0;  itr <= ptr+size-1; ++i, itr++ )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking itr+=1 and != ...\n";
    ArrayRefCountPtr<T> itr = ptr;
    for( int i = 0; itr != ptr+size; ++i, itr+=1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking itr=itr+1 and == ...\n";
    ArrayRefCountPtr<T> itr = ptr;
    for( int i = 0; !( itr == ptr+size ); ++i, itr=itr+1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }

  // Pointer --
  
  {
    if(verbose)
      *out << "\nChecking --itr and >= ...\n";
    ArrayRefCountPtr<T> itr = ptr+size-1;
    for( int i = size-1; itr >= ptr; --i, --itr )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking itr-- and > ...\n";
    ArrayRefCountPtr<T> itr = ptr+size-1;
    for( int i = size-1; itr+1 > ptr; i--, itr-- )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking itr-=1 and != ...\n";
    ArrayRefCountPtr<T> itr = ptr+size-1;
    for( int i = size-1; itr+1 != ptr; i--, itr-=1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking itr=itr-1 and == ...\n";
    ArrayRefCountPtr<T> itr = ptr+size-1;
    for( int i = size-1; !( itr+1 == ptr ); i--, itr=itr-1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }

  // Iterator ++
  
  {
    if(verbose)
      *out << "\nChecking iterator ++itr and < ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin();
    for( int i = 0; itr < ptr.end(); ++i, ++itr )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking iterator itr++ and <= ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin();
    for( int i = 0;  itr <= ptr.end()-1; ++i, itr++ )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking iterator itr+=1 and != ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin();
    for( int i = 0; itr != ptr.end(); ++i, itr+=1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking iterator itr=itr+1 and == ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin();
    for( int i = 0; !( itr == ptr.end() ); ++i, itr=itr+1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }

  // Iterator --
  
  {
    if(verbose)
      *out << "\nChecking iterator --itr and >= ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin()+size-1;
    for( int i = size-1; itr >= ptr.begin(); --i, --itr )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking iterator itr-- and > ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin()+size-1;
    for( int i = size-1; itr+1 > ptr.begin(); i--, itr-- )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking iterator itr-=1 and != ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin()+size-1;
    for( int i = size-1; itr+1 != ptr.begin(); i--, itr-=1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  {
    if(verbose)
      *out << "\nChecking iterator itr=itr-1 and == ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin()+size-1;
    for( int i = size-1; !( itr+1 == ptr.begin() ); i--, itr=itr-1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }

}

template<class T>
void test_ArrayRefCountPtr(
  const Teuchos::ArrayRefCountPtr<T>   &ptr
  ,const bool                          verbose
  ,Teuchos::FancyOStream               &out_arg
  )
{

  using Teuchos::ArrayRefCountPtr;
  using Teuchos::null;
  using Teuchos::arcp;
  using Teuchos::arcp_reinterpret_cast;

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::rcp(&out_arg,false);
  
  if(verbose)
    *out << "\nTesting ptr = " << ptr <<"\n";

  Teuchos::OSTab tab(out);

  const int size = ptr.size();
  
  {
    if(verbose)
      *out << "\nInitializing data ...\n";
    for( int i = 0; i < size; ++i )
      ptr[i] = i;
  }

  TEST_FOR_EXCEPT( !(&*ptr == ptr.get()) )

  test_ArrayRefCountPtr_iterators(ptr,verbose,out_arg);

}

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;
  using Teuchos::null;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::ArrayRefCountPtr;
  using Teuchos::arcp;
  using Teuchos::arcp_reinterpret_cast;
	
	bool success = true, verbose = true;
  
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  //const int procRank = Teuchos::GlobalMPISession::getRank();
  
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
	try {

		// Read options from the commandline
    int   num_ints = 10;
    int   num_doubles = 10;
    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "num-ints", &num_ints, "Number of ints to allocate space for" );
    clp.setOption( "num-doubles", &num_doubles, "Number of doubles to allocate space for" );
		clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
			*out << "\nEnd Result: TEST FAILED" << endl;
			return parse_return;
		}

    const int sizeOfDouble = sizeof(double);
    const int sizeOfInt = sizeof(int);

    const int total_bytes = num_doubles*sizeOfDouble + num_ints*sizeOfInt;

    if(verbose)
      *out << std::endl << Teuchos::Teuchos_Version() << std::endl;
    
		if(verbose)
			*out << "\nTesting basic ArrayRefCountPtr functionality ...\n";

    ArrayRefCountPtr<char>
      char_ptr1 = arcp<char>(total_bytes);

    if(verbose)
      *out << "\nchar_ptr1 = " << char_ptr1 << "\n";

    TEST_FOR_EXCEPT( !(char_ptr1.size() == total_bytes) );
    TEST_FOR_EXCEPT( !(char_ptr1.lowerOffset() == 0) );
    TEST_FOR_EXCEPT( !(char_ptr1.upperOffset() == total_bytes-1) );
    TEST_FOR_EXCEPT( !(char_ptr1.count() == 1) );
    test_ArrayRefCountPtr(char_ptr1,verbose,*out);

    ArrayRefCountPtr<char>
      char_ptr2 = null;

    if(verbose)
      *out << "\nchar_ptr2 = " << char_ptr2 << "\n";

    TEST_FOR_EXCEPT( !(char_ptr2.size() == 0) );
    TEST_FOR_EXCEPT( !(char_ptr2.get() == NULL) );
    TEST_FOR_EXCEPT( !(char_ptr2.count() == 0) );

    ArrayRefCountPtr<char>
      char_ptr2b(char_ptr1); // excplicitly test copy constructor

    if(verbose)
      *out << "\nchar_ptr2b = " << char_ptr2b << "\n";

    TEST_FOR_EXCEPT( !(char_ptr2b.size() == total_bytes) );
    TEST_FOR_EXCEPT( !(char_ptr2b.lowerOffset() == 0) );
    TEST_FOR_EXCEPT( !(char_ptr2b.upperOffset() == total_bytes-1) );
    TEST_FOR_EXCEPT( !(char_ptr2b.count() == 2) );
    test_ArrayRefCountPtr(char_ptr2b,verbose,*out);

    char_ptr2b = null;

    TEST_FOR_EXCEPT( !(char_ptr2b.size() == 0) );
    TEST_FOR_EXCEPT( !(char_ptr2.get() == NULL) );
    TEST_FOR_EXCEPT( !(char_ptr2b.count() == 0) );
    TEST_FOR_EXCEPT( !(char_ptr1.count() == 1) );

    ArrayRefCountPtr<char>
      char_ptr3 = char_ptr1.subview(total_bytes/2,total_bytes/2);

    TEST_FOR_EXCEPT( !(char_ptr1.count() == 2) );
    TEST_FOR_EXCEPT( !(char_ptr3.count() == 2) );
    TEST_FOR_EXCEPT( !(char_ptr3.lowerOffset() == 0) );
    TEST_FOR_EXCEPT( !(char_ptr3.upperOffset() == total_bytes/2-1) );
    test_ArrayRefCountPtr(char_ptr3,verbose,*out);

    if(verbose)
      *out << "\nchar_ptr3 = " << char_ptr3 << "\n";

    if(verbose)
      *out << "\nBreak up char_ptr1 into views of double and int data\n";

    int offset = 0;
    
    ArrayRefCountPtr<double>
      double_ptr1 = arcp_reinterpret_cast<double>(char_ptr1.subview(offset,sizeOfDouble*num_doubles));

    if(verbose)
      *out
        << "\ndouble_ptr1 = " << double_ptr1 << "\n";

    TEST_FOR_EXCEPT( !(char_ptr1.count() == 3) );
    TEST_FOR_EXCEPT( !(double_ptr1.count() == 3) );
    TEST_FOR_EXCEPT( !(double_ptr1.size() == num_doubles) );

    test_ArrayRefCountPtr(double_ptr1,verbose,*out);
    
    offset += sizeOfDouble*num_doubles;

    ArrayRefCountPtr<int>
      int_ptr1 = arcp_reinterpret_cast<int>(char_ptr1.subview(offset,sizeOfInt*num_ints));

    if(verbose)
      *out
        << "\nint_ptr1 = " << int_ptr1 << "\n";

    TEST_FOR_EXCEPT( !(char_ptr1.count() == 4) );
    TEST_FOR_EXCEPT( !(int_ptr1.count() == 4) );
    TEST_FOR_EXCEPT( !(int_ptr1.size() == num_ints) );

    test_ArrayRefCountPtr(int_ptr1,verbose,*out);

    if(verbose)
      *out << "\nCreating a constant view of double_ptr1\n";
    
    ArrayRefCountPtr<const double>
      double_ptr2 = double_ptr1.getConst();

    test_ArrayRefCountPtr_iterators(double_ptr2,verbose,*out);

#ifdef SHOW_COMPILE_FAILURE_1
    // This will not compile since this function tries to use operator[] to
    // change data but it can't since it returns a reference to a const
    // double!
    for( int i = 0; i < double_ptr2.size(); ++i ) {
      double_ptr2[i] = 1.0; // Error, you can change the value!
    }
#endif

    if(verbose)
      *out << "\nCreating an array of RefCountPtr objects!\n";

    ArrayRefCountPtr<RefCountPtr<double> >
      rcp_ptr1 = arcp<RefCountPtr<double> >(num_doubles);

    for( int i = 0; i < num_doubles; ++i )
      rcp_ptr1[i] = rcp(new double(i));

    test_ArrayRefCountPtr_iterators(rcp_ptr1,verbose,*out);

    if(verbose)
      *out << "\nCreating a const view of rcp_ptr1\n";

    ArrayRefCountPtr<const RefCountPtr<double> >
      rcp_ptr2 = rcp_ptr1.getConst();

    test_ArrayRefCountPtr_iterators(rcp_ptr2,verbose,*out);

    if(verbose)
      *out << "\nCreating an ARCP<double*> object doubleptr_ptr1 and dynamically allocation each element\n";

    ArrayRefCountPtr<double*>
      doubleptr_ptr1 = arcp<double*>(total_bytes);

    for( int i = 0; i < doubleptr_ptr1.size(); ++i )
      doubleptr_ptr1[i] = new double(i);

    test_ArrayRefCountPtr_iterators(doubleptr_ptr1,verbose,*out);

    if(verbose)
      *out << "\nCreating an ARCP<double*const> view of a doubleptr_ptr1\n";
    
    ArrayRefCountPtr<double*const>
      doubleptr_ptr2 = doubleptr_ptr1.getConst();

    test_ArrayRefCountPtr_iterators(doubleptr_ptr2,verbose,*out);

#ifdef SHOW_COMPILE_FAILURE_2
    // This will not compile since this function tries to use operator[] to
    // change data but it can't since it returns a reference to a double*const
    // object!
    for( int i = 0; i < doubleptr_ptr2.size(); ++i ) {
      *doubleptr_ptr2[i] = 1.0; // Fine, you can change the value that is being pointed to for this entry!
      doubleptr_ptr2[i] = NULL; // Error, you can't change the pointer entry!
    }
#endif

    if(verbose)
      *out << "\nCreating an ARCP<const double * const> view of a doubleptr_ptr1\n";
    
    ArrayRefCountPtr<const double*const>
      doubleptr_ptr3 = Teuchos::arcp_implicit_cast<const double*const>(doubleptr_ptr1);

    test_ArrayRefCountPtr_iterators(doubleptr_ptr3,verbose,*out);

#ifdef SHOW_COMPILE_FAILURE_3
    // This will not compile since this function tries to use operator[] to
    // change data but it can't since it returns a reference to a double*const
    // object!
    for( int i = 0; i < doubleptr_ptr3.size(); ++i ) {
      *doubleptr_ptr3[i] = 1.0; // Error, you can't change the value that is being pointed to!
      doubleptr_ptr3[i] = NULL; // Error, you can't change the pointer either!
    }
#endif

    for( int i = 0; i < doubleptr_ptr1.size(); ++i )
      delete doubleptr_ptr1[i];

    if(verbose)
      *out << "\nWrapping std::vector<T> objects as ArrayRefCount objects ...\n";

    ArrayRefCountPtr<char>
      vchar_ptr1 = arcp(rcp(new std::vector<char>(total_bytes)));

    if(verbose)
      *out << "\nvchar_ptr1 = " << vchar_ptr1 << "\n";
    
    test_ArrayRefCountPtr(vchar_ptr1,verbose,*out);
    
    ArrayRefCountPtr<const char>
      vchar_ptr2 = arcp(
        Teuchos::rcp_implicit_cast<const std::vector<char> >(
          Teuchos::get_std_vector(vchar_ptr1)
          )
        );
    
    if(verbose)
      *out << "\nvchar_ptr2 = " << vchar_ptr2 << "\n";

    test_ArrayRefCountPtr_iterators(vchar_ptr2,verbose,*out);

#ifndef __sun
    // RAB: 2006/07/12: The sun compiler declares this call to
    // get_std_vector(...) to be ambiguous (which is nonsense based on
    // everything I know about C++)!
    TEST_FOR_EXCEPT( Teuchos::get_std_vector(vchar_ptr2)->size() != static_cast<size_t>(total_bytes) );
#endif
    
    // ToDo: Fill in the rest of the tests!
    
		if(verbose)
			*out << "\nAll tests for ArrayRefCountPtr seem to check out!\n";
    
	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);
  
  if(success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;	
  
  return ( success ? 0 : 1 );
  
}
