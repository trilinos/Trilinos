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
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"

#ifdef HAVE_TEUCHOS_BOOST
#  include "Teuchos_RefCountPtrBoostSharedPtrConversions.hpp"
#endif

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

class C;

class A {
	int A_g_, A_f_;
public:
	A() : A_g_(A_g_return), A_f_(A_f_return) {}
	virtual ~A(); // See below
	virtual int A_g() { return A_g_; }
	virtual int A_f() const { return A_f_; }
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
private:
  Teuchos::RefCountPtr<C> c_;
public:
  void set_C(const Teuchos::RefCountPtr<C> &c ) { c_ = c; }
#endif // TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
};

class B1 : virtual public A {
	int B1_g_, B1_f_;
public:
	B1() : B1_g_(B1_g_return), B1_f_(B1_f_return) {}
	~B1() { B1_g_ = -1; B1_f_ = -1; }
	virtual int B1_g() { return B1_g_; }
	virtual int B1_f() const { return B1_f_; }
};

class B2 : virtual public A {
	int B2_g_, B2_f_;
public:
	B2() : B2_g_(B2_g_return), B2_f_(B2_f_return) {}
	~B2() { B2_g_ = -1; B2_f_ = -1; }
	virtual int B2_g() { return B2_g_; }
	virtual int B2_f() const { return B2_f_; }
};

class C : virtual public B1, virtual public B2
{
	int C_g_, C_f_;
public:
	C() : C_g_(C_g_return),C_f_(C_f_return) {}
	~C() { C_g_ = -1; C_f_ = -1; }
	virtual int C_g() { return C_g_; }
	virtual int C_f() const { return C_f_; }
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
private:
  Teuchos::RefCountPtr<A> a_;
public:
  void set_A(const Teuchos::RefCountPtr<A> &a ) { a_ = a; }
#endif // TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
};

// Need to put this here if we have circular references
A::~A() { A_g_ = -1; A_f_ = -1; }

class Get_A_f_return {
  const A *a_;
  int *a_f_return_;
  Get_A_f_return();
public:
  Get_A_f_return( const A *a, int *a_f_return ) : a_(a), a_f_return_(a_f_return) {}
  ~Get_A_f_return() { *a_f_return_ = a_->A_f(); }
};

void deallocA(A* ptr)
{
  std::cout << "\nCalled deallocA(...)!\n";
  delete ptr;
}

void deallocHandleA(A** handle)
{
  std::cout << "\nCalled deallocHandleA(...)!\n";
  A *ptr = *handle;
  delete ptr;
  *handle = 0;
}

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
	using Teuchos::deallocFunctorDelete;
	using Teuchos::deallocFunctorHandleDelete;
	using Teuchos::null;
	using Teuchos::rcp;
	using Teuchos::is_null;
	using Teuchos::rcp_implicit_cast;
	using Teuchos::rcp_const_cast;
	using Teuchos::rcp_static_cast;
	using Teuchos::rcp_dynamic_cast;
	using Teuchos::set_extra_data;
	using Teuchos::get_extra_data;
	using Teuchos::get_optional_extra_data;
	using Teuchos::get_dealloc;
	using Teuchos::get_optional_dealloc;
	using Teuchos::CommandLineProcessor;
	
	bool success = true, verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const int procRank = Teuchos::GlobalMPISession::getRank();

  Teuchos::oblackholestream blackhole;
  std::ostream &out = ( procRank == 0 ? std::cout : blackhole );

	try {

		// Read options from the commandline
		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
			out << "\nEnd Result: TEST FAILED" << endl;
			return parse_return;
		}

    blackhole << "\nThis should not print anywhere.\n";

    if(verbose)
      out << std::endl << Teuchos::Teuchos_Version() << std::endl;

		if(verbose)
			out << "\nTesting basic RefCountPtr functionality ...\n";

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
		TEST_FOR_EXCEPT( a_ptr1.count()  != 1 );
#endif
		TEST_FOR_EXCEPT( a_ptr1.get() == NULL );
		TEST_FOR_EXCEPT( a_ptr1 == null );
		TEST_FOR_EXCEPT( !(a_ptr1 != null) );
		TEST_FOR_EXCEPT( is_null(a_ptr1) );
		RefCountPtr<D>       d_ptr1  = rcp(new E);
#ifndef __sun
		TEST_FOR_EXCEPT( d_ptr1.count()  != 1 );
#endif
		TEST_FOR_EXCEPT( d_ptr1.get() == NULL);

		if(1) {

			// Create some more smart points (no new memory!)

			const RefCountPtr<const A> ca_ptr1 = rcp_const_cast<const A>(a_ptr1); 
      TEST_FOR_EXCEPT( !(ca_ptr1 == a_ptr1) );
      TEST_FOR_EXCEPT( ca_ptr1 != a_ptr1 );
#ifndef __sun
			TEST_FOR_EXCEPT( a_ptr1.count()  != 2 );
#endif
			TEST_FOR_EXCEPT( ca_ptr1.get() == NULL );
#ifndef __sun
			TEST_FOR_EXCEPT( ca_ptr1.count() != 2 );
#endif
			const RefCountPtr<const D> cd_ptr1 = rcp_const_cast<const D>(d_ptr1);
#ifndef __sun
			TEST_FOR_EXCEPT( d_ptr1.count()  != 2 );
#endif
			TEST_FOR_EXCEPT( cd_ptr1.get() == NULL );
#ifndef __sun
			TEST_FOR_EXCEPT( cd_ptr1.count() != 2 );
#endif

#ifdef SHOW_RUN_TIME_ERROR_1
			// Conversion using get() is a no no!  When a_ptr2 is deleted so will the allocated
			// object and then a_ptr1 will be corrupted after this block ends!
			const RefCountPtr<A> a_ptr2 = a_ptr1.get();
#endif

			// Test assignment functions

			a_ptr1 = rcp_const_cast<A>(ca_ptr1.assert_not_null()); // Should be okay, assignment to self

#ifdef SHOW_COMPILE_TIME_ERRORS
			ca_ptr1 = ca_ptr1; // Should not compile since ca_ptr1 is declared constant
			ca_ptr1->A_g();    // Should not compile since A_g() is a non-const member function
#endif

			// Test function calls through operaor->(...)

			TEST_FOR_EXCEPT( a_ptr1->A_g()  != A_g_return );
			TEST_FOR_EXCEPT( a_ptr1->A_f()  != A_f_return );
			TEST_FOR_EXCEPT( ca_ptr1->A_f() != A_f_return );
			TEST_FOR_EXCEPT( d_ptr1->D_g()  != D_g_return );
			TEST_FOR_EXCEPT( d_ptr1->D_f()  != D_f_return );
			TEST_FOR_EXCEPT( cd_ptr1->D_f() != D_f_return );
		
			// Test funciton calls through operator*(...)

			TEST_FOR_EXCEPT( (*a_ptr1).A_g()  != A_g_return );
			TEST_FOR_EXCEPT( (*a_ptr1).A_f()  != A_f_return );
			TEST_FOR_EXCEPT( (*ca_ptr1).A_f() != A_f_return );
			TEST_FOR_EXCEPT( (*d_ptr1).D_g()  != D_g_return );
			TEST_FOR_EXCEPT( (*d_ptr1).D_f()  != D_f_return );
			TEST_FOR_EXCEPT( (*cd_ptr1).D_f() != D_f_return );

			// Test dynamic and static conversions

			// Cast down the inheritance hiearchy (const A -> const B1)
			const RefCountPtr<const B1> cb1_ptr1 = rcp_dynamic_cast<const B1>(ca_ptr1);
			TEST_FOR_EXCEPT( cb1_ptr1.get() == NULL );
#ifndef __sun
			TEST_FOR_EXCEPT( cb1_ptr1.count() != 3 );
			TEST_FOR_EXCEPT( ca_ptr1.count()  != 3 );
			TEST_FOR_EXCEPT( a_ptr1.count()   != 3 );
#endif

			// Cast up the inheritance hiearchy (const B1 -> const A)
			TEST_FOR_EXCEPT( rcp_implicit_cast<const A>(cb1_ptr1)->A_f()  != A_f_return );
			TEST_FOR_EXCEPT( RefCountPtr<const A>(cb1_ptr1)->A_f()        != A_f_return );
			// Implicit cast from const to non-const (A -> const A)
			TEST_FOR_EXCEPT( rcp_implicit_cast<const A>(a_ptr1)->A_f()    != A_f_return );
			TEST_FOR_EXCEPT( RefCountPtr<const A>(a_ptr1)->A_f()          != A_f_return );
			// Cast away constantness (const B1 -> B1)
			TEST_FOR_EXCEPT( rcp_const_cast<B1>(cb1_ptr1)->B1_g()         != B1_g_return );
			// Cast across the inheritance hiearchy (const B1 -> const B2)
			TEST_FOR_EXCEPT( rcp_dynamic_cast<const B2>(cb1_ptr1)->B2_f() != B2_f_return );
			// Cast down the inheritance hiearchy (const B1 -> const C)
			TEST_FOR_EXCEPT( rcp_dynamic_cast<const C>(cb1_ptr1)->C_f()   != C_f_return );

			// Cast away constantness (const C -> C)
			const RefCountPtr<C>
				c_ptr1 = rcp_const_cast<C>(rcp_dynamic_cast<const C>(ca_ptr1));
			TEST_FOR_EXCEPT( c_ptr1.get() == NULL );
#ifndef __sun
			TEST_FOR_EXCEPT( c_ptr1.count()   != 4 );
			TEST_FOR_EXCEPT( ca_ptr1.count()  != 4 );
			TEST_FOR_EXCEPT( a_ptr1.count()   != 4 );
#endif

			// Cast down the inheritance hiearchy using static_cast<...> (const D -> const E)
			const RefCountPtr<const E>
				ce_ptr1 = rcp_static_cast<const E>(cd_ptr1); // This is not checked at runtime!
			TEST_FOR_EXCEPT( ce_ptr1.get() == NULL);
#ifndef __sun
			TEST_FOR_EXCEPT( ce_ptr1.count()  != 3 );
			TEST_FOR_EXCEPT( cd_ptr1.count()  != 3 );
			TEST_FOR_EXCEPT( d_ptr1.count()   != 3 );
#endif

			// Cast up the inheritance hiearchy (const E -> const D)
			TEST_FOR_EXCEPT( rcp_implicit_cast<const D>(ce_ptr1)->D_f()   != D_f_return ); 
			// Cast away constantness (const E -> E)
			TEST_FOR_EXCEPT( rcp_const_cast<E>(ce_ptr1)->E_g()            != E_g_return );
			TEST_FOR_EXCEPT( ce_ptr1->D_f()                               != D_f_return );

#ifdef SHOW_COMPILE_TIME_ERRORS
			// Try to cast down inheritance hiearchy using dynamic_cast<...> (const D -> const E)
			rcp_dynamic_cast<const E>( cd_ptr1 )->E_f();  // This should not compile since D and E are not polymophic
#endif

#ifndef _INTEL // Intel compiler does not seem to be doing dynamic cast correctly?
#ifdef TEUCHOS_DEBUG  // operator->() only throws exception when TEUCHOS_DEBUG is defined
			try {
				// Try to cast form one interface to another that is not supported (B2 -> B1).
				// The RefCountPtr<B1> returned from rcp_dynamic_cast<...> should be null!
				// Note that RefCountPtr<...>::optertor->() should throw an exception in debug
				// mode (i.e. TEUCHOS_DEBUG is defined) but even so no memory leak occurs.  If you
				// don't believe me then step through with a debugger and see for yourself.
				TEST_FOR_EXCEPT( rcp_dynamic_cast<B1>( rcp(new B2) )->B1_g() != B1_g_return );
				return -1; // Should not be executed!
			}
			catch( const std::logic_error &excpt )
			{}
#endif
			try {
				// Try to cast form one interface to another that is not supported (B2 -> B1).
				// Note that rcp_dynamic_cast<B1>(...,true) should throw an exception but even
				// so no memory leak occurs.  If you don't believe me then step through with a
				// debugger and see for yourself.
				rcp_dynamic_cast<B1>( rcp(new B2), true );
				return -1; // Should not be executed!
			}
			catch( const std::bad_cast &excpt )
			{}
#endif

			// Manually clean up some memory

			delete d_ptr1.release();  // Now d_ptr1.get() no longer points to a valid object but okay
			// as long as no other access to this object is attempted! (see below)
#ifdef SHOW_RUN_TIME_ERROR_2
			TEST_FOR_EXCEPT( d_ptr1->D_g() == D_g_return ); // Should cause a segmentation fault since d_ptr.get() was deleted!
#endif

#ifdef SHOW_MEMORY_LEAK_1
			a_ptr1.release(); // If we release but do not delete manually then this is a memory leak!
#endif
		
			// Here at the end of the block, all of the other smart pointers are deleted!
		}
		// Check that all of the other references where removed but these
#ifndef __sun
		TEST_FOR_EXCEPT( a_ptr1.count() != 1 );
		TEST_FOR_EXCEPT( d_ptr1.count() != 1 );
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
			out << "\nSize of C = " << sizeof(C) << std::endl;
			out << "Base address of object of type C        = " << dynamic_cast<void*>(c_ptr5) << std::endl;
			out << "Offset to address of object of type C   = " << ((long int)c_ptr5                   - (long int)c_ptr5_base) << std::endl;
			out << "Offset of B1 object in object of type C = " << ((long int)static_cast<B1*>(c_ptr5) - (long int)c_ptr5_base) << std::endl;
			out << "Offset of B2 object in object of type C = " << ((long int)static_cast<B2*>(c_ptr5) - (long int)c_ptr5_base) << std::endl;
			out << "Offset of A object in object of type C  = " << ((long int)static_cast<A*>(c_ptr5)  - (long int)c_ptr5_base) << std::endl;
		}
#endif // SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS_PRINT
		A *a_rptr5 = c_ptr5;    // Here the address has changed and is no longer the same as the base address
		a_ptr1 = rcp(a_rptr5);  // This is a no-no and could cause trouble!
		a_ptr1 = null; // This will cause a segmentation fault in free(...) on many platforms
#endif // SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS
	
		// Test out getting the deallocator object
		a_ptr1 = rcp( new C, DeallocDelete<C>(), true );
		get_dealloc<DeallocDelete<C> >(a_ptr1);
    TEST_FOR_EXCEPT( get_optional_dealloc<DeallocDelete<C> >(a_ptr1)==NULL );
    TEST_FOR_EXCEPT( get_optional_dealloc<DeallocDelete<A> >(a_ptr1)!=NULL );
    TEST_FOR_EXCEPT( get_optional_dealloc<DeallocDelete<C> >(const_cast<const RefCountPtr<A>&>(a_ptr1))==NULL );
    TEST_FOR_EXCEPT( get_optional_dealloc<DeallocDelete<A> >(const_cast<const RefCountPtr<A>&>(a_ptr1))!=NULL );
    
		// Test storing extra data and then getting it out again
		TEST_FOR_EXCEPT( get_optional_extra_data<RefCountPtr<B1> >(a_ptr1,"blahblah") != NULL );
		TEST_FOR_EXCEPT( get_optional_extra_data<int>(const_cast<const RefCountPtr<A>&>(a_ptr1),"blahblah") != NULL ); // test const version
		set_extra_data( int(-5), "int", &a_ptr1 );
		TEST_FOR_EXCEPT( get_extra_data<int>(a_ptr1,"int") != -5 );
		set_extra_data( rcp(new B1), "B1", &a_ptr1 );
		TEST_FOR_EXCEPT( get_extra_data<RefCountPtr<B1> >(a_ptr1,"B1")->B1_f() != B1_f_return );
		TEST_FOR_EXCEPT( get_extra_data<int>(const_cast<const RefCountPtr<A>&>(a_ptr1),"int") != -5 ); // test const version
		TEST_FOR_EXCEPT( (*get_optional_extra_data<RefCountPtr<B1> >(a_ptr1,"B1"))->B1_f() != B1_f_return );
		TEST_FOR_EXCEPT( *get_optional_extra_data<int>(const_cast<const RefCountPtr<A>&>(a_ptr1),"int") != -5 ); // test const version
		TEST_FOR_EXCEPT( get_optional_extra_data<RefCountPtr<B1> >(a_ptr1,"blahblah") != NULL );
		TEST_FOR_EXCEPT( get_optional_extra_data<int>(const_cast<const RefCountPtr<A>&>(a_ptr1),"blahblah") != NULL ); // test const version

    // Test pre-destruction of extra data
    int a_f_return = -2;
    set_extra_data( Teuchos::rcp(new Get_A_f_return(&*a_ptr1,&a_f_return)), "a_f_return", &a_ptr1, Teuchos::PRE_DESTROY );

		// Set pointers to null to force releasing any owned memory
		a_ptr1 = null;
		d_ptr1 = null;

#ifndef __sun
		// RAB: 2004/08/12: It appears that SUN compiler is not deleting the piece of extra
		// data properly and therefore the destructor of the above Get_A_f_return object
		// is not being called (which sets the value of af_return).  This compiler stinks!
    TEST_FOR_EXCEPT( a_f_return != A_f_return ); // Should be been called in destructor of a_ptr1 but before the A object is destroyed!
#endif

    // Testing the deallocFunctorDelete function and DeallocFunctorDelete class
    a_ptr1 = rcp( new C, deallocFunctorDelete<A>(deallocA), true );
    a_ptr1 = null;

    // Testing the deallocFunctorHandleDelete function and DeallocFunctorHandleDelete class
    a_ptr1 = rcp( new C, deallocFunctorHandleDelete<A>(deallocHandleA), true );
    a_ptr1 = null;

#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
    
		if(verbose)
			out << "\nCreate a circular reference that will case a memory leak! ...\n";
    if(1) {
      RefCountPtr<A> a = rcp(new A());
      RefCountPtr<C> c = rcp(new C());
      a->set_C(c);
      c->set_A(a);
    }

#endif // TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES

#ifdef HAVE_TEUCHOS_BOOST

		if(verbose)
			out << "\nTesting basic RefCountPtr compatibility with boost::shared_ptr ...\n";

    boost::shared_ptr<A>  a_sptr1(new C());
    RefCountPtr<A>        a_rsptr1 = rcp(a_sptr1);
		TEST_FOR_EXCEPT( a_rsptr1.get() != a_sptr1.get() );
    boost::shared_ptr<A>  a_sptr2 = shared_pointer(a_rsptr1);
		TEST_FOR_EXCEPT( a_sptr2.get() != a_sptr1.get() );
    RefCountPtr<A>        a_rsptr2 = rcp(a_sptr2);
		TEST_FOR_EXCEPT( a_rsptr2.get() != a_rsptr1.get() );
		//TEST_FOR_EXCEPT( a_rsptr2 != a_rsptr1 );  // This should work if boost::get_deleter() works correctly!
    boost::shared_ptr<A>  a_sptr3 = shared_pointer(a_rsptr2);
		TEST_FOR_EXCEPT( a_sptr3.get() != a_rsptr2.get() );

#endif // HAVE_TEUCHOS_BOOST

		if(verbose)
			out << "\nAll tests for RefCountPtr seem to check out!\n";

	} // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);
  
  try {
    // This should show that the A and C RCP objects are still around!
    Teuchos::PrivateUtilityPack::print_active_RefCountPtr_nodes(out);
	} // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);
  
  if(success)
    out << "\nEnd Result: TEST PASSED" << std::endl;	

  return ( success ? 0 : 1 );

}
