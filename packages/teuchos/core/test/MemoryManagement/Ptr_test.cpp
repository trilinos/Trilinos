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

#include "TestClasses.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


int main( int argc, char* argv[] ) {

  using Teuchos::Ptr;
  using Teuchos::ptr;
  using Teuchos::ptrFromRef;
  using Teuchos::constPtr;
  using Teuchos::outArg;
  using Teuchos::inOutArg;
  using Teuchos::inoutArg;
  using Teuchos::optInArg;
  using Teuchos::constOptInArg;
  using Teuchos::CommandLineProcessor;

        bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  //const int procRank = Teuchos::GlobalMPISession::getRank();

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

        try {

    //
                // Read options from the commandline
    //

    CommandLineProcessor clp(false); // Don't throw exceptions

                CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

                if ( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
                        *out << "\nEnd Result: TEST FAILED" << std::endl;
                        return parse_return;
                }

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl;

    *out << "\nTesting Teuchos::Ptr class ...\n";

    {
      // Test null construction
      Ptr<A> a_ptr;
      *out << "\nNull a_ptr = " << a_ptr << "\n";
      TEUCHOS_ASSERT_EQUALITY( 0, a_ptr.get() );
      TEUCHOS_ASSERT_EQUALITY( 0, a_ptr.getRawPtr() );
#ifdef TEUCHOS_DEBUG
      try {
        A &a = *a_ptr; // Should throw!
        a.A_g();
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
          "Error, Ptr::operator*() on null Ptr should have thrown exception!" );
      }
      catch( const Teuchos::NullReferenceError &except ) {
        // Caught expected exception!
      }
#endif
#ifdef TEUCHOS_DEBUG
      try {
        a_ptr->A_g(); // Should throw!
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
          "Error, Ptr::operator->() on null Ptr should have thrown exception!" );
      }
      catch( const Teuchos::NullReferenceError &except ) {
        // Caught expected exception!
      }
#endif
    }

    {
      // Test basic construction of Ptr
      A a;
      Ptr<A> a_ptr(&a);
      *out << "\nNon-null a_ptr = " << a_ptr << "\n";
      TEUCHOS_ASSERT_EQUALITY( &a, &*a_ptr );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.get() );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.getRawPtr() );
    }

    {
      // Test copy constructor for Ptr
      A a;
      Ptr<A> a_ptr1(&a);
      Ptr<A> a_ptr2(a_ptr1);
      TEUCHOS_ASSERT_EQUALITY( &*a_ptr1, &*a_ptr2 );
    }

    {
      // Test implicit copy conversion
      C c;
      Ptr<C> c_ptr(&c);
      Ptr<A> a_ptr(c_ptr);
      TEUCHOS_ASSERT_EQUALITY( &*a_ptr, &*c_ptr );
    }

    {
      // Test assignment operator
      C c;
      Ptr<C> c_ptr(&c);
      Ptr<A> a_ptr;
      a_ptr = c_ptr;
      TEUCHOS_ASSERT_EQUALITY( &*a_ptr, &*c_ptr );
    }

    {
      // Test construction of Ptr from ptr()
      A a;
      Ptr<A> a_ptr = ptr(&a);
      TEUCHOS_ASSERT_EQUALITY( &a, &*a_ptr );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.get() );
    }

    {
      // Test construction of Ptr from ptrFromRef()
      A a;
      Ptr<A> a_ptr = ptrFromRef(a);
      TEUCHOS_ASSERT_EQUALITY( &a, &*a_ptr );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.get() );
    }

    {
      // Test construction of Ptr from constPtr()
      A a;
      Ptr<const A> a_ptr = constPtr(a);
      TEUCHOS_ASSERT_EQUALITY( &a, &*a_ptr );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.get() );
    }

    {
      // Test construction of Ptr from outArg()
      A a;
      Ptr<A> a_ptr = outArg(a);
      TEUCHOS_ASSERT_EQUALITY( &a, &*a_ptr );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.get() );
    }

    {
      // Test construction of Ptr from inOutArg()
      A a;
      Ptr<A> a_ptr = inOutArg(a);
      TEUCHOS_ASSERT_EQUALITY( &a, &*a_ptr );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.get() );
    }

    {
      // Test construction of Ptr from inOutArg()
      A a;
      Ptr<A> a_ptr = inoutArg(a);
      TEUCHOS_ASSERT_EQUALITY( &a, &*a_ptr );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.get() );
    }

    {
      // Test construction of Ptr from optInArg()
      A a;
      Ptr<const A> a_ptr = optInArg(a);
      TEUCHOS_ASSERT_EQUALITY( &a, &*a_ptr );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.get() );
    }

    {
      // Test construction of Ptr from optInArg()
      A a;
      Ptr<const A> a_ptr = constOptInArg(a);
      TEUCHOS_ASSERT_EQUALITY( &a, &*a_ptr );
      TEUCHOS_ASSERT_EQUALITY( &a, a_ptr.get() );
    }

    {
      // Test ptr_implicit_cast()
      C c;
      Ptr<C> c_ptr(&c);
      Ptr<A> a_ptr1 = c_ptr;
      Ptr<A> a_ptr2 = Teuchos::ptr_implicit_cast<A>(c_ptr);
      TEUCHOS_ASSERT_EQUALITY( &*a_ptr1, &*a_ptr2 );
    }

    {
      // Test ptr_static_cast()
      E e;
      Ptr<D> d_ptr(&e);
      Ptr<E> e_ptr = Teuchos::ptr_static_cast<E>(d_ptr);
      TEUCHOS_ASSERT_EQUALITY( &*e_ptr, &e );
    }

    {
      // Test ptr_const_cast()
      C c;
      Ptr<const C> c_ptr1(&c);
      Ptr<C> c_ptr2 = Teuchos::ptr_const_cast<C>(c_ptr1);
      TEUCHOS_ASSERT_EQUALITY( &*c_ptr2, &*c_ptr1 );
    }

    {
      // Test null ptr_dynamic_cast()
      Ptr<A> a_ptr;
      Ptr<C> c_ptr = Teuchos::ptr_dynamic_cast<C>(a_ptr);
      TEUCHOS_ASSERT_EQUALITY( c_ptr.get(), 0 );
    }

    {
      // Test non-throw non-null ptr_dynamic_cast()
      C c;
      Ptr<A> a_ptr(&c);
      Ptr<C> c_ptr = Teuchos::ptr_dynamic_cast<C>(a_ptr);
      TEUCHOS_ASSERT_EQUALITY( &*c_ptr, &c );
    }

    {
      // Test good throwing non-null ptr_dynamic_cast()
      C c;
      Ptr<A> a_ptr(&c);
      Ptr<C> c_ptr = Teuchos::ptr_dynamic_cast<C>(a_ptr,true);
      TEUCHOS_ASSERT_EQUALITY( &*c_ptr, &c );
    }

    {
      // Test bad throwing non-null ptr_dynamic_cast()
      B1 b1;
      Ptr<A> a_ptr(&b1);
      try {
        Ptr<C> b2_ptr = Teuchos::ptr_dynamic_cast<C>(a_ptr,true);
        (void) b2_ptr; // Silence "set but not used" compiler warning.
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
          "If you get here then the test failed!" );
      }
      catch ( const Teuchos::m_bad_cast &except ) {
        // Test passed!
      }
    }

        }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;

  return ( success ? 0 : 1 );

}
