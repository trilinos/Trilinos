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

//#define SHOW_COMPILE_FAILURE_1;

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

  const int dim = ptr.dim();

  // Pointer ++
  
  if(1) {
    if(verbose)
      *out << "\nChecking ++itr and < ...\n";
    ArrayRefCountPtr<T> itr = ptr;
    for( int i = 0; itr < ptr+dim; ++i, ++itr )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking itr++ and <= ...\n";
    ArrayRefCountPtr<T> itr = ptr;
    for( int i = 0;  itr <= ptr+dim-1; ++i, itr++ )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking itr+=1 and != ...\n";
    ArrayRefCountPtr<T> itr = ptr;
    for( int i = 0; itr != ptr+dim; ++i, itr+=1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking itr=itr+1 and == ...\n";
    ArrayRefCountPtr<T> itr = ptr;
    for( int i = 0; !( itr == ptr+dim ); ++i, itr=itr+1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }

  // Pointer --
  
  if(1) {
    if(verbose)
      *out << "\nChecking --itr and >= ...\n";
    ArrayRefCountPtr<T> itr = ptr+dim-1;
    for( int i = dim-1; itr >= ptr; --i, --itr )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking itr-- and > ...\n";
    ArrayRefCountPtr<T> itr = ptr+dim-1;
    for( int i = dim-1; itr+1 > ptr; i--, itr-- )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking itr-=1 and != ...\n";
    ArrayRefCountPtr<T> itr = ptr+dim-1;
    for( int i = dim-1; itr+1 != ptr; i--, itr-=1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking itr=itr-1 and == ...\n";
    ArrayRefCountPtr<T> itr = ptr+dim-1;
    for( int i = dim-1; !( itr+1 == ptr ); i--, itr=itr-1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }

  // Iterator ++
  
  if(1) {
    if(verbose)
      *out << "\nChecking iterator ++itr and < ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin();
    for( int i = 0; itr < ptr+dim; ++i, ++itr )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking iterator itr++ and <= ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin();
    for( int i = 0;  itr <= ptr+dim-1; ++i, itr++ )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking iterator itr+=1 and != ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin();
    for( int i = 0; itr != ptr+dim; ++i, itr+=1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking iterator itr=itr+1 and == ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin();
    for( int i = 0; !( itr == ptr+dim ); ++i, itr=itr+1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }

  // Iterator --
  
  if(1) {
    if(verbose)
      *out << "\nChecking iterator --itr and >= ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin()+dim-1;
    for( int i = dim-1; itr >= ptr; --i, --itr )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking iterator itr-- and > ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin()+dim-1;
    for( int i = dim-1; itr+1 > ptr; i--, itr-- )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking iterator itr-=1 and != ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin()+dim-1;
    for( int i = dim-1; itr+1 != ptr; i--, itr-=1 )
      TEST_FOR_EXCEPT( !(*itr == ptr[i]) );
  }
  
  if(1) {
    if(verbose)
      *out << "\nChecking iterator itr=itr-1 and == ...\n";
    typename ArrayRefCountPtr<T>::const_iterator itr = ptr.begin()+dim-1;
    for( int i = dim-1; !( itr+1 == ptr ); i--, itr=itr-1 )
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

  const int dim = ptr.dim();
  
  if(1) {
    if(verbose)
      *out << "\nInitializing data ...\n";
    for( int i = 0; i < dim; ++i )
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

    TEST_FOR_EXCEPT( !(char_ptr1.dim() == total_bytes) );
    TEST_FOR_EXCEPT( !(char_ptr1.lowerOffset() == 0) );
    TEST_FOR_EXCEPT( !(char_ptr1.upperOffset() == total_bytes-1) );
    TEST_FOR_EXCEPT( !(char_ptr1.count() == 1) );
    test_ArrayRefCountPtr(char_ptr1,verbose,*out);

    ArrayRefCountPtr<char>
      char_ptr2 = null;

    if(verbose)
      *out << "\nchar_ptr2 = " << char_ptr2 << "\n";

    TEST_FOR_EXCEPT( !(char_ptr2.dim() == 0) );
    TEST_FOR_EXCEPT( !(char_ptr2.get() == NULL) );
    TEST_FOR_EXCEPT( !(char_ptr2.count() == 0) );

    ArrayRefCountPtr<char>
      char_ptr2b(char_ptr1); // excplicitly test copy constructor

    if(verbose)
      *out << "\nchar_ptr2b = " << char_ptr2b << "\n";

    TEST_FOR_EXCEPT( !(char_ptr2b.dim() == total_bytes) );
    TEST_FOR_EXCEPT( !(char_ptr2b.lowerOffset() == 0) );
    TEST_FOR_EXCEPT( !(char_ptr2b.upperOffset() == total_bytes-1) );
    TEST_FOR_EXCEPT( !(char_ptr2b.count() == 2) );
    test_ArrayRefCountPtr(char_ptr2b,verbose,*out);

    char_ptr2b = null;

    TEST_FOR_EXCEPT( !(char_ptr2b.dim() == 0) );
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
    TEST_FOR_EXCEPT( !(double_ptr1.dim() == num_doubles) );

    test_ArrayRefCountPtr(double_ptr1,verbose,*out);
    
    offset += sizeOfDouble*num_doubles;

    ArrayRefCountPtr<int>
      int_ptr1 = arcp_reinterpret_cast<int>(char_ptr1.subview(offset,sizeOfInt*num_ints));

    if(verbose)
      *out
        << "\nint_ptr1 = " << int_ptr1 << "\n";

    TEST_FOR_EXCEPT( !(char_ptr1.count() == 4) );
    TEST_FOR_EXCEPT( !(int_ptr1.count() == 4) );
    TEST_FOR_EXCEPT( !(int_ptr1.dim() == num_ints) );

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
    test_ArrayRefCountPtr(double_ptr2,verbose,*out);
#endif

    if(verbose)
      *out << "\nCreating an array of RefCountPtr objects!\n";

    ArrayRefCountPtr<RefCountPtr<double> >
      rcp_ptr1 = arcp<RefCountPtr<double> >(num_doubles);

    for( int i = 0; i < num_doubles; ++i )
      rcp_ptr1[i] = rcp(new double(i));

    test_ArrayRefCountPtr_iterators(rcp_ptr1,verbose,*out);

    if(verbose)
      *out << "\nCreating a constant view of double_ptr1\n";

    ArrayRefCountPtr<const RefCountPtr<double> >
      rcp_ptr2 = rcp_ptr1.getConst();

    test_ArrayRefCountPtr_iterators(rcp_ptr2,verbose,*out);

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

    TEST_FOR_EXCEPT( Teuchos::get_std_vector(vchar_ptr2)->size() != static_cast<size_t>(total_bytes) );
    
    // ToDo: Fill in the rest of the tests!
    
		if(verbose)
			*out << "\nAll tests for ArrayRefCountPtr seem to check out!\n";
    
	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);
  
  if(success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;	
  
  return ( success ? 0 : 1 );
  
}
