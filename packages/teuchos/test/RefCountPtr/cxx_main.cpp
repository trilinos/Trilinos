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

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Version.hpp"

// Return constants from class functions
const int
A_g_return  = 1,
	A_f_return  = 2,
	B1_g_return = 3,
	B1_f_return = 4,
	B2_g_return = 5,
	B2_f_return = 6,
	C_g_return  = 7,
	C_f_return  = 8,
	D_g_return  = 9,
	D_f_return  = 10,
	E_g_return  = 11,
	E_f_return  = 12;

/*

 Polymorphic multiple inheritance example

            -----
           |  A  |
            -----
             /|\
              | 
         ------------
        |            |
      -----        ------
     |  B1 |      |  B2  |
      -----        ------
       /|\          /|\
        |            |
         ------------
              |
            -----
           |  C  |
            -----

*/

class A {
	int A_g_, A_f_;
public:
	A() : A_g_(A_g_return), A_f_(A_f_return) {}
	virtual ~A() {}
	virtual int A_g() { return A_g_; }
	virtual int A_f() const { return A_f_; }
};

class B1 : virtual public A {
	int B1_g_, B1_f_;
public:
	B1() : B1_g_(B1_g_return), B1_f_(B1_f_return) {}
	virtual int B1_g() { return B1_g_; }
	virtual int B1_f() const { return B1_f_; }
};

class B2 : virtual public A {
	int B2_g_, B2_f_;
public:
	B2() : B2_g_(B2_g_return), B2_f_(B2_f_return) {}
	virtual int B2_g() { return B2_g_; }
	virtual int B2_f() const { return B2_f_; }
};

class C : virtual public B1, virtual public B2
{
	int C_g_, C_f_;
public:
	C() : C_g_(C_g_return), C_f_(C_f_return) {}
	virtual int C_g() { return C_g_; }
	virtual int C_f() const { return C_f_; }
	
};

/*

 Non-polymophic classes hiearchy examlpe

            -----
           |  D  |
            -----
             /|\
              | 
            -----
           |  E  |
            -----

*/

class D 
{
	int D_g_, D_f_;
public:
	D() : D_g_(D_g_return), D_f_(D_f_return) {}
	int D_g() { return D_g_; }
	int D_f() const { return D_f_; }
};

class E : public D
{
	int E_g_, E_f_;
public:
	E() : E_g_(E_g_return), E_f_(E_f_return) {}
	int E_g() { return E_g_; }
	int E_f() const { return E_f_; }
};

// //////////////////////////

// 
// Uncomment these macros to see example errors
//

//#define SHOW_COMPILE_TIME_ERRORS
//#define SHOW_RUN_TIME_ERROR_1
//#define SHOW_RUN_TIME_ERROR_2
//#define SHOW_RUN_TIME_ERROR_3
//#define SHOW_RUN_TIME_ERROR_4
#define SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS
#define SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS_PRINT
//#define SHOW_MEMORY_LEAK_1

//
// This program prints minimal output to standard error
//

int main( int argc, char* argv[] ) {

	using Teuchos::RefCountPtr;
	using Teuchos::DeallocDelete;
	using Teuchos::rcp;
	using Teuchos::rcp_implicit_cast;
	using Teuchos::rcp_const_cast;
	using Teuchos::rcp_static_cast;
	using Teuchos::rcp_dynamic_cast;
	using Teuchos::set_extra_data;
	using Teuchos::get_extra_data;
	using Teuchos::get_dealloc;
	using Teuchos::CommandLineProcessor;
	
	bool verbose = true;

	try {

		// Read options from the commandline
		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

	if (verbose)
		cout << Teuchos::Teuchos_Version() << endl << endl;

		// Create some smart pointers

		RefCountPtr<A>       a_ptr1  = rcp(new C);
#ifndef __sun
		// RAB: 2003/11/24: The Sun compiler ("Forte Developer 7 C++
		// 5.4 2002/03/09" returned from CC -V) does not seem to be
		// following the standard when it comes to the handling of
		// temporary objects and therefore the count() is not currect.
		// In the above statement, at least one and perhaps two
		// temporary RefCountPtr objects are created before the object
		// a_ptr1 is initialized.  However, the standard says that the
		// lifetime of temprary objects must not extend past the
		// statement in which it was created (see section 10.4.10 in
		// Stroustroup, 3ed edition).  This compiler stinks!!!!!
		assert( a_ptr1.count()  == 1 );
#endif
		assert( a_ptr1.get() != NULL );
		RefCountPtr<D>       d_ptr1  = rcp(new E);
#ifndef __sun
		assert( d_ptr1.count()  == 1 );
#endif
		assert( d_ptr1.get() != NULL);

		if(1) {

			// Create some more smart points (no new memory!)

			const RefCountPtr<const A> ca_ptr1 = rcp_const_cast<const A>(a_ptr1); 
#ifndef __sun
			assert( a_ptr1.count()  == 2 );
#endif
			assert( ca_ptr1.get() != NULL );
#ifndef __sun
			assert( ca_ptr1.count() == 2 );
#endif
			const RefCountPtr<const D> cd_ptr1 = rcp_const_cast<const D>(d_ptr1);
#ifndef __sun
			assert( d_ptr1.count()  == 2 );
#endif
			assert( cd_ptr1.get() != NULL );
#ifndef __sun
			assert( cd_ptr1.count() == 2 );
#endif

#ifdef SHOW_RUN_TIME_ERROR_1
			// Conversion using get() is a no no!  When a_ptr2 is deleted so will the allocated
			// object and then a_ptr1 will be corrupted after this block ends!
			const RefCountPtr<A> a_ptr2 = a_ptr1.get();
#endif

			// Test assignment functions

			a_ptr1 = rcp_const_cast<A>(ca_ptr1); // Should be okay, assignment to self

#ifdef SHOW_COMPILE_TIME_ERRORS
			ca_ptr1 = ca_ptr1; // Should not compile since ca_ptr1 is declared constant
			ca_ptr1->A_g();    // Should not compile since A_g() is a non-const member function
#endif

			// Test function calls through operaor->(...)

			assert( a_ptr1->A_g()  == A_g_return );
			assert( a_ptr1->A_f()  == A_f_return );
			assert( ca_ptr1->A_f() == A_f_return );
			assert( d_ptr1->D_g()  == D_g_return );
			assert( d_ptr1->D_f()  == D_f_return );
			assert( cd_ptr1->D_f() == D_f_return );
		
			// Test funciton calls through operator*(...)

			assert( (*a_ptr1).A_g()  == A_g_return );
			assert( (*a_ptr1).A_f()  == A_f_return );
			assert( (*ca_ptr1).A_f() == A_f_return );
			assert( (*d_ptr1).D_g()  == D_g_return );
			assert( (*d_ptr1).D_f()  == D_f_return );
			assert( (*cd_ptr1).D_f() == D_f_return );

			// Test dynamic and static conversions

			// Cast down the inheritance hiearchy (const A -> const B1)
			const RefCountPtr<const B1> cb1_ptr1 = rcp_dynamic_cast<const B1>(ca_ptr1);
			assert( cb1_ptr1.get() != NULL );
#ifndef __sun
			assert( cb1_ptr1.count() == 3 );
			assert( ca_ptr1.count()  == 3 );
			assert( a_ptr1.count()   == 3 );
#endif

			// Cast up the inheritance hiearchy (const B1 -> const A)
			assert( rcp_implicit_cast<const A>(cb1_ptr1)->A_f()  == A_f_return );
			assert( RefCountPtr<const A>(cb1_ptr1)->A_f()        == A_f_return );
			// Implicit cast from const to non-const (A -> const A)
			assert( rcp_implicit_cast<const A>(a_ptr1)->A_f()    == A_f_return );
			assert( RefCountPtr<const A>(a_ptr1)->A_f()          == A_f_return );
			// Cast away constantness (const B1 -> B1)
			assert( rcp_const_cast<B1>(cb1_ptr1)->B1_g()         == B1_g_return );
			// Cast across the inheritance hiearchy (const B1 -> const B2)
			assert( rcp_dynamic_cast<const B2>(cb1_ptr1)->B2_f() == B2_f_return );
			// Cast down the inheritance hiearchy (const B1 -> const C)
			assert( rcp_dynamic_cast<const C>(cb1_ptr1)->C_f()   == C_f_return );

			// Cast away constantness (const C -> C)
			const RefCountPtr<C>
				c_ptr1 = rcp_const_cast<C>(rcp_dynamic_cast<const C>(ca_ptr1));
			assert( c_ptr1.get() != NULL );
#ifndef __sun
			assert( c_ptr1.count()   == 4 );
			assert( ca_ptr1.count()  == 4 );
			assert( a_ptr1.count()   == 4 );
#endif

			// Cast down the inheritance hiearchy using static_cast<...> (const D -> const E)
			const RefCountPtr<const E>
				ce_ptr1 = rcp_static_cast<const E>(cd_ptr1); // This is not checked at runtime!
			assert( ce_ptr1.get() != NULL);
#ifndef __sun
			assert( ce_ptr1.count()  == 3 );
			assert( cd_ptr1.count()  == 3 );
			assert( d_ptr1.count()   == 3 );
#endif

			// Cast up the inheritance hiearchy (const E -> const D)
			assert( rcp_implicit_cast<const D>(ce_ptr1)->D_f()   == D_f_return ); 
			// Cast away constantness (const E -> E)
			assert( rcp_const_cast<E>(ce_ptr1)->E_g()            == E_g_return );
			assert( ce_ptr1->D_f()                                    == D_f_return );

#ifdef SHOW_COMPILE_TIME_ERRORS
			// Try to cast down inheritance hiearchy using dynamic_cast<...> (const D -> const E)
			rcp_dynamic_cast<const E>( cd_ptr1 )->E_f();  // This should not compile since D and E are not polymophic
#endif

#ifndef _INTEL
			try {
				// Try to cast form one interface to another that is not supported (B2 -> B1).
				// The RefCountPtr<B1> returned from rcp_dynamic_cast<...> should be null!
				// Note that RefCountPtr<...>::optertor->() should throw an exception but even
				// so no memory leak occurs.  If you don't believe me then step through with a
				// debugger and see for yourself.
				assert( rcp_dynamic_cast<B1>( rcp(new B2) )->B1_g() == B1_g_return );
				return -1; // Should not be executed!
			}
			catch( const std::logic_error )
			{}
#endif

			// Manually clean up some memory

			delete d_ptr1.release();  // Now d_ptr1.get() no longer points to a valid object but okay
			// as long as no other access to this object is attempted! (see below)
#ifdef SHOW_RUN_TIME_ERROR_2
			assert( d_ptr1->D_g() == D_g_return ); // Should cause a segmentation fault since d_ptr.get() was deleted!
#endif

#ifdef SHOW_MEMORY_LEAK_1
			a_ptr1.release(); // If we release but do not delete manually then this is a memory leak!
#endif
		
			// Here at the end of the block, all of the other smart pointers are deleted!
		}
		// Check that all of the other references where removed but these
#ifndef __sun
		assert( a_ptr1.count() == 1 );
		assert( d_ptr1.count() == 1 );
#endif

		// Assign some other dynamically created objects.
	
		a_ptr1 = rcp(new A);  // In each case the current dynamically allocated object is deleted ...
		a_ptr1 = rcp(new B1); // before the new reference is set.
		a_ptr1 = rcp(new B2); // ""
		a_ptr1 = rcp(new C);  // ""
		d_ptr1 = rcp(new D);  // ""
		d_ptr1 = rcp(new E);  // ""

		// Assign pointers to some automatic objects that do not need deleted.
		// We can do this but we need to remove ownership of the pointer
		// from the smart pointer objects so that they do not try to
		// delete them.  If we forget then delete will be called on these
		// pointers and will cause a runtime error.

		C c; // Automatic object what will be deleted by compiler at end of block
		a_ptr1 = rcp(&c);
#ifndef SHOW_RUN_TIME_ERROR_3
		// Release ownership so that a_ptr1 will not try to delete &c when a_ptr1 goes out of scope
		a_ptr1.release();
#endif	

		E e; // Automatic object what will be deleted by compiler at end of block
		d_ptr1 = rcp(&e);
#ifndef SHOW_RUN_TIME_ERROR_4
		// Release ownership so that d_ptr1 will not try to delete &e when a_ptr1 goes out of scope
		d_ptr1.release();
#endif

#ifdef SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS
		// Allocate an using new and then store the non-base address in in
		// a RefCountPtr and then try to delete (this is a no-no usually).
		C *c_ptr5 = new C;      // Okay, no type info lost and address should be same as returned from malloc(...)
#ifdef SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS_PRINT
		const void *c_ptr5_base = dynamic_cast<void*>(c_ptr5);
		if(verbose) {
			std::cout << "\nSize of C = " << sizeof(C) << std::endl;
			std::cout << "Base address of object of type C        = " << dynamic_cast<void*>(c_ptr5) << std::endl;
			std::cout << "Offset to address of object of type C   = " << ((long int)c_ptr5                   - (long int)c_ptr5_base) << std::endl;
			std::cout << "Offset of B1 object in object of type C = " << ((long int)static_cast<B1*>(c_ptr5) - (long int)c_ptr5_base) << std::endl;
			std::cout << "Offset of B2 object in object of type C = " << ((long int)static_cast<B2*>(c_ptr5) - (long int)c_ptr5_base) << std::endl;
			std::cout << "Offset of A object in object of type C  = " << ((long int)static_cast<A*>(c_ptr5)  - (long int)c_ptr5_base) << std::endl;
		}
#endif // SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS_PRINT
		A *a_rptr5 = c_ptr5;    // Here the address has changed and is no longer the same as the base address
		a_ptr1 = rcp(a_rptr5);  // This is a no-no and could cause trouble!
		a_ptr1 = Teuchos::null; // This will cause a segmentation fault in free(...) on many platforms
#endif // SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS
	
		// Test out getting the deallocator object
		a_ptr1 = rcp( new C, DeallocDelete<C>(), true );
		get_dealloc<DeallocDelete<C> >(a_ptr1);
    
		// Test storing extra data and then getting it out again
		set_extra_data( int(-5), "int", &a_ptr1 );
		assert( get_extra_data<int>(a_ptr1,"int") == -5 );
		set_extra_data( rcp(new B1), "B1", &a_ptr1 );
		assert( get_extra_data<RefCountPtr<B1> >(a_ptr1,"B1")->B1_f() == B1_f_return );
		assert( get_extra_data<int>(const_cast<const RefCountPtr<A>&>(a_ptr1),"int") == -5 ); // test const version

		// Set pointers to null to force releasing any owned memory
		a_ptr1 = Teuchos::null;
		d_ptr1 = Teuchos::null;

		if(verbose)
			std::cout << "RefCountPtr<...> seems to check out!\n";

	} // end try
	catch( const std::exception &excpt ) {
		if(verbose)
			std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
		return 1;
	}
	catch( ... ) {
		if(verbose)
			std::cerr << "*** Caught an unknown exception\n";
		return 1;
	}

	return 0;

}
