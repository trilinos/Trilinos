// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_Version.hpp"

#ifdef HAVE_TEUCHOS_BOOST
#  include "Teuchos_RCPBoostSharedPtrConversions.hpp"
#endif

#include "TestClasses.hpp"


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

  using Teuchos::RCP;
  using Teuchos::DeallocDelete;
  using Teuchos::deallocFunctorDelete;
  using Teuchos::deallocFunctorHandleDelete;
  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::inOutArg;
  using Teuchos::is_null;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_static_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;
  using Teuchos::get_extra_data;
  using Teuchos::get_nonconst_extra_data;
  using Teuchos::get_optional_extra_data;
  using Teuchos::get_optional_nonconst_extra_data;
  using Teuchos::get_dealloc;
  using Teuchos::get_nonconst_dealloc;
  using Teuchos::get_optional_dealloc;
  using Teuchos::get_optional_nonconst_dealloc;
  using Teuchos::rcpWithEmbeddedObj;
  using Teuchos::rcpWithEmbeddedObjPreDestroy;
  using Teuchos::rcpWithEmbeddedObjPostDestroy;
  using Teuchos::getEmbeddedObj;
  using Teuchos::getNonconstEmbeddedObj;
  using Teuchos::getConst;
  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool createCircRefs = false; // Don't create memory leak by default!

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const int procRank = Teuchos::GlobalMPISession::getRank();

  Teuchos::oblackholestream blackhole;
  std::ostream &out = ( procRank == 0 ? std::cout : blackhole );

  try {

    // Read options from the commandline
    CommandLineProcessor clp(false); // Don't throw exceptions
    clp.setOption( "create-circ-refs", "no-create-circ-refs", &createCircRefs,
      "Set if output is printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
      out << "\nEnd Result: TEST FAILED" << std::endl;
      return parse_return;
    }

    blackhole << "\nThis should not print anywhere.\n";

    out << std::endl << Teuchos::Teuchos_Version() << std::endl;

    out << "\nTesting basic RCP functionality ...\n";

    // Create some smart pointers

    RCP<A> a_ptr1 = rcp(new C);
    out << "\na_ptr1 = " << a_ptr1 << "\n";
    // RAB: 2003/11/24: The Sun compiler ("Forte Developer 7 C++
    // 5.4 2002/03/09" returned from CC -V) does not seem to be
    // following the standard when it comes to the handling of
    // temporary objects and therefore the count() is not currect.
    // In the above statement, at least one and perhaps two
    // temporary RCP objects are created before the object
    // a_ptr1 is initialized. However, the standard says that the
    // lifetime of temprary objects must not extend past the
    // statement in which it was created (see section 10.4.10 in
    // Stroustroup, 3ed edition). This compiler stinks!!!!!
    TEUCHOS_TEST_FOR_EXCEPT( a_ptr1.strong_count() != 1 );
    TEUCHOS_TEST_FOR_EXCEPT( !a_ptr1.shares_resource(a_ptr1) );
    TEUCHOS_TEST_FOR_EXCEPT( a_ptr1.ptr() == null );
    TEUCHOS_TEST_FOR_EXCEPT( a_ptr1 == null );
    TEUCHOS_TEST_FOR_EXCEPT( !(a_ptr1 != null) );
    TEUCHOS_TEST_FOR_EXCEPT( is_null(a_ptr1) );
    RCP<D> d_ptr1 = rcp(new E);
    TEUCHOS_TEST_FOR_EXCEPT( d_ptr1.shares_resource(a_ptr1) );
    TEUCHOS_TEST_FOR_EXCEPT( d_ptr1.strong_count() != 1 );
    TEUCHOS_TEST_FOR_EXCEPT( d_ptr1.get() == NULL);
    TEUCHOS_TEST_FOR_EXCEPT( d_ptr1.getRawPtr() == NULL);

    {

      // Create some more smart points (no new memory!)

      const RCP<const A> ca_ptr1 = rcp_const_cast<const A>(a_ptr1);
      TEUCHOS_TEST_FOR_EXCEPT( !(ca_ptr1 == a_ptr1) );
      TEUCHOS_TEST_FOR_EXCEPT( ca_ptr1 != a_ptr1 );
      TEUCHOS_TEST_FOR_EXCEPT( !ca_ptr1.shares_resource(a_ptr1) );
      TEUCHOS_TEST_FOR_EXCEPT( a_ptr1.strong_count() != 2 );
      TEUCHOS_TEST_FOR_EXCEPT( ca_ptr1.ptr() == null );
      TEUCHOS_TEST_FOR_EXCEPT( ca_ptr1.strong_count() != 2 );
      const RCP<const D> cd_ptr1 = rcp_const_cast<const D>(d_ptr1);
      TEUCHOS_TEST_FOR_EXCEPT( d_ptr1.strong_count() != 2 );
      TEUCHOS_TEST_FOR_EXCEPT( cd_ptr1.ptr() == null );
      TEUCHOS_TEST_FOR_EXCEPT( cd_ptr1.strong_count() != 2 );

#ifdef SHOW_RUN_TIME_ERROR_1
      // Conversion using get() is a no no! When a_ptr2 is deleted so will the allocated
      // object and then a_ptr1 will be corrupted after this block ends!
      const RCP<A> a_ptr2 = a_ptr1.get();
#endif

      // Test assignment functions

      a_ptr1 = rcp_const_cast<A>(ca_ptr1.assert_not_null()); // Should be okay, assignment to self

#ifdef SHOW_COMPILE_TIME_ERRORS
      ca_ptr1 = ca_ptr1; // Should not compile since ca_ptr1 is declared constant
      ca_ptr1->A_g(); // Should not compile since A_g() is a non-const member function
#endif

      // Test function calls through operaor->(...)

      TEUCHOS_TEST_FOR_EXCEPT( a_ptr1->A_g() != A_g_return );
      TEUCHOS_TEST_FOR_EXCEPT( a_ptr1->A_f() != A_f_return );
      TEUCHOS_TEST_FOR_EXCEPT( ca_ptr1->A_f() != A_f_return );
      TEUCHOS_TEST_FOR_EXCEPT( d_ptr1->D_g() != D_g_return );
      TEUCHOS_TEST_FOR_EXCEPT( d_ptr1->D_f() != D_f_return );
      TEUCHOS_TEST_FOR_EXCEPT( cd_ptr1->D_f() != D_f_return );

      // Test funciton calls through operator*(...)

      TEUCHOS_TEST_FOR_EXCEPT( (*a_ptr1).A_g() != A_g_return );
      TEUCHOS_TEST_FOR_EXCEPT( (*a_ptr1).A_f() != A_f_return );
      TEUCHOS_TEST_FOR_EXCEPT( (*ca_ptr1).A_f() != A_f_return );
      TEUCHOS_TEST_FOR_EXCEPT( (*d_ptr1).D_g() != D_g_return );
      TEUCHOS_TEST_FOR_EXCEPT( (*d_ptr1).D_f() != D_f_return );
      TEUCHOS_TEST_FOR_EXCEPT( (*cd_ptr1).D_f() != D_f_return );

      // Test dynamic and static conversions

      // Cast down the inheritance hiearchy (const A -> const B1)
      const RCP<const B1> cb1_ptr1 = rcp_dynamic_cast<const B1>(ca_ptr1);
      TEUCHOS_TEST_FOR_EXCEPT( cb1_ptr1.ptr() == null );
      TEUCHOS_TEST_FOR_EXCEPT( cb1_ptr1.strong_count() != 3 );
      TEUCHOS_TEST_FOR_EXCEPT( ca_ptr1.strong_count() != 3 );
      TEUCHOS_TEST_FOR_EXCEPT( a_ptr1.strong_count() != 3 );

      // Cast up the inheritance hiearchy (const B1 -> const A)
      TEUCHOS_TEST_FOR_EXCEPT( rcp_implicit_cast<const A>(cb1_ptr1)->A_f() != A_f_return );
      TEUCHOS_TEST_FOR_EXCEPT( RCP<const A>(cb1_ptr1)->A_f() != A_f_return );
      // Implicit cast from const to non-const (A -> const A)
      TEUCHOS_TEST_FOR_EXCEPT( rcp_implicit_cast<const A>(a_ptr1)->A_f() != A_f_return );
      TEUCHOS_TEST_FOR_EXCEPT( RCP<const A>(a_ptr1)->A_f() != A_f_return );
      // Cast away constantness (const B1 -> B1)
      TEUCHOS_TEST_FOR_EXCEPT( rcp_const_cast<B1>(cb1_ptr1)->B1_g() != B1_g_return );
      // Cast across the inheritance hiearchy (const B1 -> const B2)
      TEUCHOS_TEST_FOR_EXCEPT( rcp_dynamic_cast<const B2>(cb1_ptr1)->B2_f() != B2_f_return );
      // Cast down the inheritance hiearchy (const B1 -> const C)
      TEUCHOS_TEST_FOR_EXCEPT( rcp_dynamic_cast<const C>(cb1_ptr1)->C_f() != C_f_return );

      // Cast away constantness (const C -> C)
      const RCP<C>
        c_ptr1 = rcp_const_cast<C>(rcp_dynamic_cast<const C>(ca_ptr1));
      TEUCHOS_TEST_FOR_EXCEPT( c_ptr1.ptr() == null );
      TEUCHOS_TEST_FOR_EXCEPT( c_ptr1.strong_count() != 4 );
      TEUCHOS_TEST_FOR_EXCEPT( ca_ptr1.strong_count() != 4 );
      TEUCHOS_TEST_FOR_EXCEPT( a_ptr1.strong_count() != 4 );

      // Cast down the inheritance hiearchy using static_cast<...> (const D -> const E)
      const RCP<const E>
        ce_ptr1 = rcp_static_cast<const E>(cd_ptr1); // This is not checked at runtime!
      TEUCHOS_TEST_FOR_EXCEPT( ce_ptr1.ptr() == null);
      TEUCHOS_TEST_FOR_EXCEPT( ce_ptr1.strong_count() != 3 );
      TEUCHOS_TEST_FOR_EXCEPT( cd_ptr1.strong_count() != 3 );
      TEUCHOS_TEST_FOR_EXCEPT( d_ptr1.strong_count() != 3 );

      // Cast up the inheritance hiearchy (const E -> const D)
      TEUCHOS_TEST_FOR_EXCEPT( rcp_implicit_cast<const D>(ce_ptr1)->D_f() != D_f_return );
      // Cast away constantness (const E -> E)
      TEUCHOS_TEST_FOR_EXCEPT( rcp_const_cast<E>(ce_ptr1)->E_g() != E_g_return );
      TEUCHOS_TEST_FOR_EXCEPT( ce_ptr1->D_f() != D_f_return );

#ifdef SHOW_COMPILE_TIME_ERRORS
      // Try to cast down inheritance hiearchy using dynamic_cast<...> (const D -> const E)
      rcp_dynamic_cast<const E>( cd_ptr1 )->E_f(); // This should not compile since D and E are not polymophic
#endif

#ifndef _INTEL // Intel compiler does not seem to be doing dynamic cast correctly?
#ifdef TEUCHOS_DEBUG // operator->() only throws std::exception when TEUCHOS_DEBUG is defined
      try {
        // Try to cast form one interface to another that is not supported (B2 -> B1).
        // The RCP<B1> returned from rcp_dynamic_cast<...> should be null!
        // Note that RCP<...>::optertor->() should throw an std::exception in debug
        // mode (i.e. TEUCHOS_DEBUG is defined) but even so no memory leak occurs. If you
        // don't believe me then step through with a debugger and see for yourself.
        TEUCHOS_TEST_FOR_EXCEPT( rcp_dynamic_cast<B1>( rcp(new B2) )->B1_g() != B1_g_return );
        return -1; // Should not be executed!
      }
      catch( const std::logic_error &excpt )
      {}
#endif
      try {
        // Try to cast form one interface to another that is not supported (B2 -> B1).
        // Note that rcp_dynamic_cast<B1>(...,true) should throw an std::exception but even
        // so no memory leak occurs. If you don't believe me then step through with a
        // debugger and see for yourself.
        rcp_dynamic_cast<B1>( rcp(new B2), true );
        return -1; // Should not be executed!
      }
      catch( const std::bad_cast &excpt )
      {}
#endif

      // Manually clean up some memory

      delete d_ptr1.release().get(); // Now d_ptr1.get() no longer points to a valid object but okay
      // as long as no other access to this object is attempted! (see below)
#ifdef SHOW_RUN_TIME_ERROR_2
      TEUCHOS_TEST_FOR_EXCEPT( d_ptr1->D_g() == D_g_return ); // Should cause a segmentation fault since d_ptr.get() was deleted!
#endif

#ifdef SHOW_MEMORY_LEAK_1
      a_ptr1.release(); // If we release but do not delete manually then this is a memory leak!
#endif

      // Here at the end of the block, all of the other smart pointers are deleted!
    }
    // Check that all of the other references where removed but these
    TEUCHOS_TEST_FOR_EXCEPT( a_ptr1.strong_count() != 1 );
    TEUCHOS_TEST_FOR_EXCEPT( d_ptr1.strong_count() != 1 );

    // Assign some other dynamically created objects.

    a_ptr1 = rcp(new A); // In each case the current dynamically allocated object is deleted ...
    a_ptr1 = rcp(new B1); // before the new reference is set.
    a_ptr1 = rcp(new B2); // ""
    a_ptr1 = rcp(new C); // ""
    d_ptr1 = rcp(new D); // ""
    d_ptr1 = rcp(new E); // ""

    // Assign pointers to some automatic objects that do not need deleted.
    // We can do this but we need to remove ownership of the pointer
    // from the smart pointer objects so that they do not try to
    // delete them. If we forget then delete will be called on these
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
    // a RCP and then try to delete (this is a no-no usually).
    C *c_ptr5 = new C; // Okay, no type info lost and address should be same as returned from malloc(...)
#ifdef SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS_PRINT
    const void *c_ptr5_base = dynamic_cast<void*>(c_ptr5);
    out << "\nSize of C = " << sizeof(C) << std::endl;
    out << "Base address of object of type C = " << dynamic_cast<void*>(c_ptr5) << std::endl;
    out << "Offset to address of object of type C = " << ((long int)c_ptr5 - (long int)c_ptr5_base) << std::endl;
    out << "Offset of B1 object in object of type C = " << ((long int)static_cast<B1*>(c_ptr5) - (long int)c_ptr5_base) << std::endl;
    out << "Offset of B2 object in object of type C = " << ((long int)static_cast<B2*>(c_ptr5) - (long int)c_ptr5_base) << std::endl;
    out << "Offset of A object in object of type C = " << ((long int)static_cast<A*>(c_ptr5) - (long int)c_ptr5_base) << std::endl;
#endif // SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS_PRINT
    A *a_rptr5 = c_ptr5; // Here the address has changed and is no longer the same as the base address
    a_ptr1 = rcp(a_rptr5); // This is a no-no and could cause trouble!
    a_ptr1 = null; // This will cause a segmentation fault in free(...) on many platforms
#endif // SHOW_RUN_TIME_ERROR_VIRTUAL_BASE_CLASS

    // Test out getting the deallocator object
    a_ptr1 = rcpWithDealloc( new C, DeallocDelete<C>() );
    get_dealloc<DeallocDelete<C> >(a_ptr1);
    get_nonconst_dealloc<DeallocDelete<C> >(a_ptr1);
    TEUCHOS_TEST_FOR_EXCEPT( get_optional_nonconst_dealloc<DeallocDelete<C> >(a_ptr1)==null );
    TEUCHOS_TEST_FOR_EXCEPT( get_optional_nonconst_dealloc<DeallocDelete<A> >(a_ptr1)!=null );
    TEUCHOS_TEST_FOR_EXCEPT( get_optional_dealloc<DeallocDelete<C> >(const_cast<const RCP<A>&>(a_ptr1))==null );
    TEUCHOS_TEST_FOR_EXCEPT( get_optional_dealloc<DeallocDelete<A> >(const_cast<const RCP<A>&>(a_ptr1))!=null );

    // Test storing extra data and then getting it out again
    TEUCHOS_TEST_FOR_EXCEPT( get_optional_nonconst_extra_data<RCP<B1> >(a_ptr1,"blahblah") != null );
    TEUCHOS_TEST_FOR_EXCEPT( get_optional_extra_data<int>(const_cast<const RCP<A>&>(a_ptr1),"blahblah") != null ); // test const version
    set_extra_data( int(-5), "int", inOutArg(a_ptr1) );
    TEUCHOS_TEST_FOR_EXCEPT( get_extra_data<int>(a_ptr1,"int") != -5 );
    TEUCHOS_TEST_FOR_EXCEPT( get_nonconst_extra_data<int>(a_ptr1,"int") != -5 );
    set_extra_data( rcp(new B1), "B1", inOutArg(a_ptr1) );
    TEUCHOS_TEST_FOR_EXCEPT( get_extra_data<RCP<B1> >(a_ptr1,"B1")->B1_f() != B1_f_return );
    TEUCHOS_TEST_FOR_EXCEPT( get_extra_data<int>(const_cast<const RCP<A>&>(a_ptr1),"int") != -5 ); // test const version
    TEUCHOS_TEST_FOR_EXCEPT( (*get_optional_extra_data<RCP<B1> >(a_ptr1,"B1"))->B1_f() != B1_f_return );
    TEUCHOS_TEST_FOR_EXCEPT( *get_optional_extra_data<int>(const_cast<const RCP<A>&>(a_ptr1),"int") != -5 ); // test const version
    TEUCHOS_TEST_FOR_EXCEPT( get_optional_extra_data<RCP<B1> >(a_ptr1,"blahblah") != null );
    TEUCHOS_TEST_FOR_EXCEPT( get_optional_extra_data<int>(const_cast<const RCP<A>&>(a_ptr1),"blahblah") != null ); // test const version

    // Test storage of extra data as embedded objects and then getting it out
    // again

    {
      RCP<A> a_ptr = rcpWithEmbeddedObj(new C,int(-5));
      const int intRtn1 = getEmbeddedObj<C,int>(a_ptr);
      TEUCHOS_TEST_FOR_EXCEPT( intRtn1 != -5 );
      getNonconstEmbeddedObj<C,int>(a_ptr) = -4;
      const int intRtn2 = getEmbeddedObj<C,int>(a_ptr);
      TEUCHOS_TEST_FOR_EXCEPT( intRtn2 != -4 );
    }

    {
      RCP<A> a_ptr = rcpWithEmbeddedObjPreDestroy(new C,int(-5));
      const int intRtn1 = getEmbeddedObj<C,int>(a_ptr);
      TEUCHOS_TEST_FOR_EXCEPT( intRtn1 != -5 );
      getNonconstEmbeddedObj<C,int>(a_ptr) = -4;
      const int intRtn2 = getEmbeddedObj<C,int>(a_ptr);
      TEUCHOS_TEST_FOR_EXCEPT( intRtn2 != -4 );
    }

    {
      RCP<A> a_ptr = rcpWithEmbeddedObjPostDestroy(new C,int(-5));
      const int intRtn1 = getEmbeddedObj<C,int>(a_ptr);
      TEUCHOS_TEST_FOR_EXCEPT( intRtn1 != -5 );
      getNonconstEmbeddedObj<C,int>(a_ptr) = -4;
      const int intRtn2 = getEmbeddedObj<C,int>(a_ptr);
      TEUCHOS_TEST_FOR_EXCEPT( intRtn2 != -4 );
    }

    // Test pre-destruction of extra data
    int a_f_return = -2;
    set_extra_data( rcp(new Get_A_f_return(&*a_ptr1,&a_f_return)),
      "a_f_return", inOutArg(a_ptr1), Teuchos::PRE_DESTROY );

    // Set pointers to null to force releasing any owned memory
    a_ptr1 = null;
    d_ptr1 = null;

    // RAB: 2004/08/12: It appears that SUN compiler is not deleting the piece of extra
    // data properly and therefore the destructor of the above Get_A_f_return object
    // is not being called (which sets the value of af_return). This compiler stinks!
    TEUCHOS_TEST_FOR_EXCEPT( a_f_return != A_f_return ); // Should be been called in destructor of a_ptr1 but before the A object is destroyed!

    // Testing the deallocFunctorDelete function and DeallocFunctorDelete class
    a_ptr1 = rcpWithDealloc( new C, deallocFunctorDelete<A>(deallocA) );
    a_ptr1 = null;

    // Testing the deallocFunctorHandleDelete function and DeallocFunctorHandleDelete class
    a_ptr1 = rcpWithDealloc( new C, deallocFunctorHandleDelete<A>(deallocHandleA) );
    a_ptr1 = null;

#ifdef TEUCHOS_DEBUG

    if (createCircRefs) {
      out << "\nCreate a circular reference that will cause a memory leak! ...\n";
#  if !defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)
      // Only trun on tracing if you have to
      Teuchos::RCPNodeTracer::setTracingActiveRCPNodes(true);
#  endif
      RCP<A> a = rcp(new A());
      RCP<C> c2 = rcp(new C());
      a->set_C(c2);
      c2->set_A(a);
    }

#endif // TEUCHOS_DEBUG

#ifndef TEUCHOS_DEBUG

    out << "\nTesting using RCP to wrap an undefined opaque object (no TNT) ...\n";
    {
      RCP<UndefinedType> op_ptr =
        rcpWithDeallocUndef (createOpaque (),
                             deallocFunctorHandleDelete<UndefinedType> (destroyOpaque),
                             true);
      TEUCHOS_ASSERT_EQUALITY( getOpaqueValue(&*op_ptr), getOpaqueValue_return );
    }
    // 2008/08/01: rabartl: Above, we can only wrap an undefined type in
    // nondebug mode since there is no TypeNameTraits class defined for it and
    // the default uses typeid(...) which you can't call on an undefined type.
    // If you define a specialization of TypeNameTraits for this class, then
    // it will compile just fine!  This is related to bug 4016.

#endif // not TEUCHOS_DEBUG

    out << "\nTesting using RCP to wrap an undefined opaque object (with TNT) ...\n";
    {
      RCP<UndefinedType2> op_ptr = rcpWithDeallocUndef( createOpaque2(),
        deallocFunctorHandleDelete<UndefinedType2>(destroyOpaque2) );
      TEUCHOS_ASSERT_EQUALITY( getOpaque2Value(&*op_ptr), getOpaque2Value_return );
    }
    // 2008/08/01: rabartl: Above, we can wrap an undefined type in debug mode
    // as long as we have a TypeNameTraits specialization of it to avoid
    // calling typeid(...).

#ifdef HAVE_TEUCHOS_BOOST

    out << "\nTesting basic RCP compatibility with boost::shared_ptr ...\n";

    boost::shared_ptr<A> a_sptr1(new C());
    RCP<A> a_rsptr1 = rcp(a_sptr1);
    TEUCHOS_TEST_FOR_EXCEPT( a_rsptr1.get() != a_sptr1.get() );
    TEUCHOS_TEST_FOR_EXCEPT( a_rsptr1.getRawPtr() != a_sptr1.get() );
    TEUCHOS_TEST_FOR_EXCEPT( a_rsptr1.get() != a_rsptr1.getRawPtr() );
    boost::shared_ptr<A> a_sptr2 = shared_pointer(a_rsptr1);
    // There seems no standard way to test that a shared_ptr shares the same node
    //TEUCHOS_TEST_FOR_EXCEPT( a_sptr2._internal_equiv(a_sptr1) != true );
    RCP<A> a_rsptr2 = rcp(a_sptr2);
    TEUCHOS_TEST_FOR_EXCEPT( a_rsptr2.ptr() != a_rsptr1.ptr() );
    //TEUCHOS_TEST_FOR_EXCEPT( a_rsptr2 != a_rsptr1 ); // This should work if boost::get_deleter() works correctly!
    boost::shared_ptr<A> a_sptr3 = shared_pointer(a_rsptr2);
    TEUCHOS_TEST_FOR_EXCEPT( a_sptr3.get() != a_rsptr2.get() );

    out << "\nCompatibility with boost::shared_ptr passed ...\n";

#endif // HAVE_TEUCHOS_BOOST

    out << "\nAll tests for RCP seem to check out!\n";

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  try {
    // In debug mode, this should show that the A and C RCP objects are still
    // around!
    if (createCircRefs) {
      out << "\nPrinting the active nodes just to see them!\n";
      Teuchos::RCPNodeTracer::printActiveRCPNodes(out);
#if defined(TEUCHOS_DEBUG) && !defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)
      TEUCHOS_ASSERT_EQUALITY( 2, Teuchos::RCPNodeTracer::numActiveRCPNodes() );
#endif
    }
  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  if(success)
    out << "\nEnd Result: TEST PASSED" << std::endl;

  return ( success ? 0 : 1 );

}
