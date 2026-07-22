// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Array_UnitTest_helpers.hpp"
#include "TestClasses.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_getRawPtr.hpp"

namespace {

using ArrayUnitTestHelpers::n;
using ArrayUnitTestHelpers::generateArray;

typedef Teuchos_Ordinal Ordinal;
using Teuchos::getRawPtr;
using Teuchos::as;
using Teuchos::null;
using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::arcp;
using Teuchos::arcpCloneNode;
using Teuchos::arcp_reinterpret_cast;
using Teuchos::arcp_reinterpret_cast_nonpod;
using Teuchos::ArrayView;
using Teuchos::getConst;
using Teuchos::DuplicateOwningRCPError;
using Teuchos::NullReferenceError;
using Teuchos::DanglingReferenceError;
using Teuchos::RangeError;
using Teuchos::RCP_STRONG;
using Teuchos::RCP_WEAK;
using Teuchos::implicit_ptr_cast;
using Teuchos::getRawPtr;


//
// Non templated unit tests
//


TEUCHOS_UNIT_TEST( ArrayRCP, memberPointer )
{
  ArrayRCP<A> a_arcp = arcp<A>(1);
  TEST_EQUALITY_CONST( a_arcp->A_f(), A_f_return );
}


TEUCHOS_UNIT_TEST( ArrayRCP, getConst_null )
{
  const ArrayRCP<A> a1_arcp;
  const ArrayRCP<const A> a2_arcp = a1_arcp.getConst();
  TEST_ASSERT(is_null(a2_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, operator_parenth_ArrayView_null )
{
  const ArrayRCP<A> a_arcp;
  const ArrayView<A> av = a_arcp();
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST( ArrayRCP, operator_parenth_ArrayView_const_null )
{
  const ArrayRCP<const A> a_arcp;
  const ArrayView<const A> av = a_arcp();
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST( ArrayRCP, null_zero_ArrayView_operator )
{
  const ArrayRCP<const A> a_arcp;
  const ArrayView<const A> av = a_arcp(0, 0);
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST( ArrayRCP, null_zero_ArrayView_view_func )
{
  const ArrayRCP<const A> a_arcp;
  const ArrayView<const A> av = a_arcp.view(0, 0);
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST( ArrayRCP, null_zero_ArrayView_persistingView )
{
  const ArrayRCP<const A> a_arcp;
  const ArrayRCP<const A> a_arcp2 = a_arcp.persistingView(0, 0);
  TEST_ASSERT(is_null(a_arcp2));
}


TEUCHOS_UNIT_TEST( ArrayRCP, raw_ptr_nonowning_self_view )
{
  A *data = new A[10];
  ArrayRCP<A> arcp_view(data, 0, 10, false);
  ArrayView<A> view = arcp_view(0, 5);
  arcp_view = null;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW(view.size(), DanglingReferenceError);
#endif
  delete [] data;
}


TEUCHOS_UNIT_TEST( ArrayRCP, implicit_ArrayRCP_const )
{
  const ArrayRCP<A> a_arcp;
  const ArrayRCP<const A> ac_arcp = a_arcp;
  TEST_ASSERT(is_null(ac_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, ArrayRCP_void_throws )
{
  TEST_THROW( const ArrayRCP<      void>  v_arcp, std::logic_error );
  TEST_THROW( const ArrayRCP<const void> cv_arcp, std::logic_error );
}


TEUCHOS_UNIT_TEST( ArrayRCP, release )
{
  ArrayRCP<A> a_arcp = arcp<A>(1);
  delete [] a_arcp.release();
}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_null )
{
  ArrayRCP<A> a_arcp = arcp<A>(0, 0, -1, false);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_dealloc_null )
{
  ArrayRCP<A> a_arcp = arcp<A, Teuchos::DeallocNull<A> >(0, 0, -1,
    Teuchos::DeallocNull<A>(), false);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, convert_from_vector_null )
{
  const RCP<std::vector<int> > v_rcp;
  const ArrayRCP<int> a_arcp = arcp(v_rcp);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, convert_from_const_vector_null )
{
  const RCP<const std::vector<int> > v_rcp;
  const ArrayRCP<const int> a_arcp = arcp(v_rcp);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, convert_from_vector_unsized )
{
  const RCP<std::vector<int> > v_rcp = rcp(new std::vector<int>);
  const ArrayRCP<int> a_arcp = arcp(v_rcp);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, convert_from_const_vector_unsized )
{
  const RCP<const std::vector<int> > v_rcp = rcp(new std::vector<int>);
  const ArrayRCP<const int> a_arcp = arcp(v_rcp);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, arcpWithEmbeddedObj )
{
  const ArrayRCP<const int> a_arcp =
    Teuchos::arcpWithEmbeddedObj<int>(new int[1], 0, 1, as<int>(1), true);
  const int embeddedObj = Teuchos::getEmbeddedObj<int,int>(a_arcp);
  TEST_EQUALITY_CONST( embeddedObj, as<int>(1) );
}


TEUCHOS_UNIT_TEST( ArrayRCP, nonnull )
{
  ECHO(ArrayRCP<int> a_arcp = arcp<int>(10));
  TEST_EQUALITY_CONST(is_null(a_arcp), false);
  TEST_EQUALITY_CONST(nonnull(a_arcp), true);
  ECHO(a_arcp = null);
  TEST_EQUALITY_CONST(is_null(a_arcp), true);
  TEST_EQUALITY_CONST(nonnull(a_arcp), false);
}


TEUCHOS_UNIT_TEST( ArrayRCP, weak_strong )
{

  ECHO(ArrayRCP<int> arcp1 = arcp<int>(10));
  TEST_EQUALITY_CONST( arcp1.strength(), RCP_STRONG );

  ECHO(ArrayRCP<int> arcp2 = arcp1.create_weak());

  TEST_EQUALITY_CONST( arcp2.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( arcp1.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp1.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp2.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp2.weak_count(), 1 );

  ECHO(ArrayRCP<int> arcp3 = arcp2.create_strong());

  TEST_EQUALITY_CONST( arcp3.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( arcp1.strong_count(), 2 );
  TEST_EQUALITY_CONST( arcp1.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp2.strong_count(), 2 );
  TEST_EQUALITY_CONST( arcp2.weak_count(), 1 );

  // This will make the underlying object A gets deleted!
  ECHO(arcp1 = null);
  ECHO(arcp3 = null);

  ECHO(arcp2 = null); // Should make the underlying node go away

}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_reinterpret_cast_null )
{
  ECHO(ArrayRCP<char> arcp_char = null);
  ECHO(ArrayRCP<int> arcp_int = arcp_reinterpret_cast<int>(arcp_char));
  TEST_EQUALITY_CONST(arcp_int, null);
}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_reinterpret_cast_char_to_int )
{

  const int sizeOfInt = sizeof(int);
  const int sizeOfChar = sizeof(char);
  const int num_ints = n;
  const int num_chars = (num_ints*sizeOfInt)/sizeOfChar;
  out << "num_ints = " << num_ints << "\n";
  out << "num_chars = " << num_chars << "\n";

  ECHO(ArrayRCP<char> arcp_char = arcp<char>(num_chars));
  ECHO(ArrayRCP<int> arcp_int = arcp_reinterpret_cast<int>(arcp_char));
  TEST_EQUALITY(arcp_int.size(), num_ints);
  TEST_EQUALITY(implicit_ptr_cast<void>(&arcp_int[0]),
    implicit_ptr_cast<void>(&arcp_char[0]));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-1])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1])+1));

  ECHO(arcp_char+=sizeOfInt);
  ECHO(arcp_int = arcp_reinterpret_cast<int>(arcp_char));
  TEST_EQUALITY(arcp_int.size(), num_ints);
  TEST_EQUALITY_CONST( arcp_int.lowerOffset(), -1);
  TEST_EQUALITY( arcp_int.upperOffset(), num_ints-2);
  TEST_EQUALITY( implicit_ptr_cast<void>(&arcp_int[-1]),
    implicit_ptr_cast<void>(&arcp_char[-sizeOfInt])
    );
  TEST_EQUALITY( implicit_ptr_cast<void>((&arcp_int[num_ints-2])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1-sizeOfInt])+1));

}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_reinterpret_cast_int_to_char )
{

  const int sizeOfInt = sizeof(int);
  const int sizeOfChar = sizeof(char);
  const int num_ints = n;
  const int num_chars = (num_ints*sizeOfInt)/sizeOfChar;
  out << "num_ints = " << num_ints << "\n";
  out << "num_chars = " << num_chars << "\n";

  ECHO(ArrayRCP<int> arcp_int = arcp<int>(num_ints));
  ECHO(ArrayRCP<char> arcp_char = arcp_reinterpret_cast<char>(arcp_int));
  TEST_EQUALITY(arcp_char.size(), num_chars);
  TEST_EQUALITY(implicit_ptr_cast<void>(&arcp_int[0]),
    implicit_ptr_cast<void>(&arcp_char[0]));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-1])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1])+1));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-1])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1])+1));

  ECHO(++arcp_int);
  ECHO(arcp_char = arcp_reinterpret_cast<char>(arcp_int));
  TEST_EQUALITY(as<int>(arcp_char.lowerOffset()), as<int>(-sizeOfInt));
  TEST_EQUALITY(as<int>(arcp_char.upperOffset()), as<int>(num_chars-1-sizeOfInt));
  TEST_EQUALITY(implicit_ptr_cast<void>(&arcp_int[-1]),
    implicit_ptr_cast<void>(&arcp_char[-sizeOfInt]));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-2])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1-sizeOfInt])+1));

}


TEUCHOS_UNIT_TEST( ArrayRCP, evil_reinterpret_cast )
{
  ECHO(ArrayRCP<ArrayRCP<int> > arcp1 = arcp<ArrayRCP<int> >(n));
  ECHO(ArrayRCP<ArrayRCP<const int> > arcp2 =
    arcp_reinterpret_cast<ArrayRCP<const int> >(arcp1));
  TEST_EQUALITY(arcp2.size(), arcp1.size());
  TEST_EQUALITY(implicit_ptr_cast<const void>(&arcp1[0]),
    implicit_ptr_cast<const void>(&arcp2[0]));
  ECHO(ArrayRCP<const ArrayRCP<const int> > arcp3 = arcp2);
  TEST_EQUALITY(arcp3.size(), arcp1.size());
  TEST_EQUALITY(implicit_ptr_cast<const void>(&arcp1[0]),
    implicit_ptr_cast<const void>(&arcp3[0]));
  out << "arcp3 = " << arcp3 << "\n";
}


//
// Test arcpCloneNode(...)
//


TEUCHOS_UNIT_TEST( ArrayRCP, arcpCloneNode_null )
{
  ECHO(ArrayRCP<ArrayRCP<int> > arcp1 = null);
  ECHO(ArrayRCP<ArrayRCP<int> > arcp2 = arcpCloneNode(arcp1));
  TEST_EQUALITY(arcp2, null);
}


TEUCHOS_UNIT_TEST( ArrayRCP, arcpCloneNode_basic )
{

  ECHO(ArrayRCP<int> arcp1 = arcp<int>(n));

  ECHO(ArrayRCP<int> arcp2 = arcpCloneNode(arcp1));
  TEST_ASSERT(nonnull(arcp2));
  TEST_EQUALITY(arcp1.strong_count(), 2);
  TEST_EQUALITY(arcp2.strong_count(), 1);

  ECHO(ArrayRCP<int> arcp3 = arcp2);
  TEST_EQUALITY(arcp1.strong_count(), 2);
  TEST_EQUALITY(arcp2.strong_count(), 2);
  TEST_EQUALITY(arcp3.strong_count(), 2);

  ECHO(ArrayRCP<int> arcp4 = arcp1);
  TEST_EQUALITY(arcp1.strong_count(), 3);
  TEST_EQUALITY(arcp2.strong_count(), 2);
  TEST_EQUALITY(arcp3.strong_count(), 2);

  ECHO(arcp4 = null);
  TEST_EQUALITY(arcp1.strong_count(), 2);
  TEST_EQUALITY(arcp2.strong_count(), 2);
  TEST_EQUALITY(arcp3.strong_count(), 2);
  TEST_EQUALITY(arcp4.strong_count(), 0);

  ECHO(arcp1 = null);
  TEST_EQUALITY(arcp1.strong_count(), 0);
  TEST_EQUALITY(arcp2.strong_count(), 2);
  TEST_EQUALITY(arcp3.strong_count(), 2);
  TEST_EQUALITY(arcp4.strong_count(), 0);

  ECHO(arcp2 = null);
  TEST_EQUALITY(arcp2.strong_count(), 0);
  TEST_EQUALITY(arcp3.strong_count(), 1);

  ECHO(arcp3 = null);
  TEST_EQUALITY(arcp3.strong_count(), 0);

}


//
// Test arcp_reinterpret_cast_nonpod(...)
//


class MockObject {
  int member_;
public:

  MockObject(int member_in = -1) : member_(member_in) { ++(numConstructorsCalled()); }
  MockObject(const MockObject &mo) : member_(mo.member_) { ++(numCopyConstructorsCalled()); }
  ~MockObject() { ++(numDestructorsCalled()); }
  int member() const { return member_; }

  static int & numConstructorsCalled()
    { static int s_numConstructorsCalled = 0; return s_numConstructorsCalled; }
  static int & numCopyConstructorsCalled()
    { static int s_numCopyConstructorsCalled = 0; return s_numCopyConstructorsCalled; }
  static int & numDestructorsCalled()
    { static int s_numDestructorsCalled = 0; return s_numDestructorsCalled; }
  static void reset() { numConstructorsCalled() = numCopyConstructorsCalled() = numDestructorsCalled() = 0; }

};


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_reinterpret_cast_nonpod_default_construct )
{

  const int sizeOfMockObject = sizeof(MockObject);
  const int sizeOfChar = sizeof(char);
  const int num_objs = n;
  const int num_chars = (num_objs*sizeOfMockObject)/sizeOfChar;
  out << "num_objs = " << num_objs << "\n";
  out << "num_chars = " << num_chars << "\n";

  ECHO(ArrayRCP<char> arcp_chars = arcp<char>(num_chars));

  ECHO(MockObject::reset());
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 0);
  ECHO(ArrayRCP<MockObject> arcp_objs =
    arcp_reinterpret_cast_nonpod<MockObject>(arcp_chars));
  TEST_EQUALITY(arcp_objs.size(), num_objs);
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 1);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), num_objs);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 1);
  {
    int sum = 0; for (int i=0; i < num_objs; ++i) sum += arcp_objs[i].member();
    TEST_EQUALITY(sum, -num_objs);
  }

  ECHO(ArrayRCP<MockObject> arcp_objs2 = arcp_objs);
  TEST_EQUALITY(arcp_objs.size(), num_objs);
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 1);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), num_objs);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 1);
  {
    int sum = 0; for (int i=0; i < num_objs; ++i) sum += arcp_objs[i].member();
    TEST_EQUALITY(sum, -num_objs);
  }

  ECHO(arcp_objs = null);
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 1);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), num_objs);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 1);

  ECHO(arcp_objs2 = null);
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 1);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), num_objs);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), num_objs + 1);

}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_reinterpret_cast_nonpod_copy_construct )
{

  const int sizeOfMockObject = sizeof(MockObject);
  const int sizeOfChar = sizeof(char);
  const int num_objs = n;
  const int num_chars = (num_objs*sizeOfMockObject)/sizeOfChar;
  out << "num_objs = " << num_objs << "\n";
  out << "num_chars = " << num_chars << "\n";

  ECHO(ArrayRCP<char> arcp_chars = arcp<char>(num_chars));

  ECHO(MockObject::reset());
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 0);
  ECHO(const MockObject mockObj(1));
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 1);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 0);

  ECHO(MockObject::reset());
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 0);
  ECHO(ArrayRCP<MockObject> arcp_objs =
    arcp_reinterpret_cast_nonpod(arcp_chars, mockObj));
  TEST_EQUALITY(arcp_objs.size(), num_objs);
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), num_objs);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 0);
  {
    int sum = 0; for (int i=0; i < num_objs; ++i) sum += arcp_objs[i].member();
    TEST_EQUALITY(sum, num_objs);
  }

  ECHO(ArrayRCP<MockObject> arcp_objs2 = arcp_objs);
  TEST_EQUALITY(arcp_objs.size(), num_objs);
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), num_objs);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 0);
  {
    int sum = 0; for (int i=0; i < num_objs; ++i) sum += arcp_objs[i].member();
    TEST_EQUALITY(sum, num_objs);
  }

  ECHO(arcp_objs = null);
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), num_objs);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), 0);

  ECHO(arcp_objs2 = null);
  TEST_EQUALITY(MockObject::numConstructorsCalled(), 0);
  TEST_EQUALITY(MockObject::numCopyConstructorsCalled(), num_objs);
  TEST_EQUALITY(MockObject::numDestructorsCalled(), num_objs);

}


//
// Test catching of duplicate owning ArrayRCP objects
//


TEUCHOS_UNIT_TEST( ArrayRCP, duplicate_arcp_owning )
{
  SET_RCPNODE_TRACING();
  ECHO(A *a_ptr = new A[n]);
  ECHO(ArrayRCP<A> a_arcp1 = arcp(a_ptr, 0, n)); // Okay
#if defined(TEUCHOS_DEBUG)
  // With node tracing turned on, the implementation knows that an RCPNode
  // already exists pointing to this same underlying array and will therefore
  // throw.
  TEST_THROW(ArrayRCP<A> a_arcp2 = arcp(a_ptr, 0, n), DuplicateOwningRCPError);
#else
  // Will not determine they are point to the same object!
  ECHO(ArrayRCP<A> a_arcp2 = arcp(a_ptr, 0, n));
  TEST_EQUALITY(a_arcp2.getRawPtr(), a_ptr);
  ECHO(a_arcp2.release()); // Better or we will get a segfault!
#endif
}


TEUCHOS_UNIT_TEST( ArrayRCP, dangling_nonowning )
{
  SET_RCPNODE_TRACING();
  ECHO(A *a_ptr = new A[n]);
  ECHO(ArrayRCP<A> a_arcp1 = arcp(a_ptr, 0, n)); // Okay
  ECHO(ArrayRCP<A> a_arcp2 = arcp(a_ptr, 0, n, false)); // Okay
  a_arcp1 = null;
#if defined(TEUCHOS_DEBUG)
  // With node tracing turned on, the implementation knows that the original
  // array is deleted and this is a dangling reference!
  TEST_THROW(a_arcp2.getRawPtr(), DanglingReferenceError);
#else
  // With node tracing turned off, the implemetation does not know the
  // original array is deleted and therefore it will return a now invalid
  // pointer.
  TEST_NOTHROW(a_arcp2.getRawPtr());
#endif
}


class WeirdDealloc {
  int size_;
  RCP<std::ostream > out_;
public:
  WeirdDealloc(int size, const RCP<std::ostream> &out) : size_(size), out_(out) {}
  void free(void *ptr) const
    {
      int * const int_ptr = reinterpret_cast<int*>(ptr);
      {
        // Create an ArrayView that points to the same memory that is being
        // deallocated by the owning ArrayRCP.  Here, if RCPNode tracing is
        // enabled, this will thrown and there is really no way around it.
        ArrayView<const int> tmpav(int_ptr, size_, Teuchos::RCP_DISABLE_NODE_LOOKUP);
        assert(tmpav[0] == int_ptr[0]);
        *out_ << tmpav << std::endl;
        // Create a copy of the ArrayView and make sure that it does not do
        // node tracing either.
        ArrayView<const int> tmpav2(tmpav);
        assert(tmpav2[0] == int_ptr[0]);
        *out_ << tmpav2 << std::endl;
        // Assign the ArrayView and make sure that it does not do node tracing
        // either.
        ArrayView<const int> tmpav3;
        tmpav3 = tmpav;
        assert(tmpav3[0] == int_ptr[0]);
        *out_ << tmpav2 << std::endl;
      }
      delete [] int_ptr;
    }
};


TEUCHOS_UNIT_TEST( ArrayRCP, weirdDealloc )
{
  using Teuchos::rcpFromRef;
  const int size = 4;
  const bool ownsMem = true;
  int *int_ptr = new int[size];
  std::fill_n(int_ptr, size, 0);
  ArrayRCP<int> a = arcp<int>( int_ptr , 0, size,
    WeirdDealloc(size, rcpFromRef(out)), ownsMem );
  a = Teuchos::null;
}


//
// Templated unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, construct_n, T )
{
  std::vector<T> a(n, as<T>(1));
  ArrayRCP<T> a_arcp(n, as<T>(1));
  TEST_COMPARE_ARRAYS(a, a_arcp);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, assignSelf, T )
{
  ArrayRCP<T> a_arcp;
  a_arcp = a_arcp;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, assign_n_val, T )
{
  const T val = as<T>(1);
  std::vector<T> a;
  a.assign(n, val);
  ArrayRCP<T> a_arcp;
  a_arcp.assign(as<Ordinal>(n), val);
  TEST_COMPARE_ARRAYS(a, a_arcp);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, assign_begin_end, T )
{
  const T val = as<T>(1);
  std::vector<T> a;
  a.assign(n, val);
  ArrayRCP<T> a_arcp;
  a_arcp.assign(a.begin(), a.end());
  TEST_COMPARE_ARRAYS(a, a_arcp);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, print_iterators, T )
{
  typedef typename ArrayRCP<T>::const_iterator const_iterator;
  ECHO(ArrayRCP<T> a_arcp = arcp<T>(n));
  ECHO(const_iterator itr = a_arcp.begin());
  out << "itr = " << itr << "\n";
  TEST_EQUALITY(itr, a_arcp.begin());
  ECHO(itr += n);
  out << "itr = " << itr << "\n";
  TEST_EQUALITY(itr, a_arcp.end());
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, deepCopy, T )
{
  const T val = as<T>(1);
  std::vector<T> a;
  a.assign(n, val);
  ArrayRCP<T> a_arcp = arcp<T>(n);
  ArrayRCP<T> a_arcp_cpy = a_arcp;
  a_arcp.deepCopy(Teuchos::arrayViewFromVector(a));
  TEST_COMPARE_ARRAYS(a, a_arcp);
  TEST_EQUALITY(a_arcp.getRawPtr(), a_arcp_cpy.getRawPtr());
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, resize, T )
{
  const T val1 = as<T>(1);
  const T val2 = as<T>(2);

  std::vector<T> a;
  ArrayRCP<T> a_arcp;

  out << "\nChecking resize(n, val1) ...\n";
  a.resize(n, val1);
  a_arcp.resize(n, val1);
  TEST_COMPARE_ARRAYS(a, a_arcp);

  out << "\nChecking resize(2*n, val2) ...\n";
  a.resize(2*n, val2);
  a_arcp.resize(2*n, val2);
  TEST_COMPARE_ARRAYS(a, a_arcp);

  out << "\nChecking resize(n/2) ...\n";
  a.resize(n/2);
  a_arcp.resize(n/2);
  TEST_COMPARE_ARRAYS(a, a_arcp);

  out << "\nChecking resize(0) ...\n";
  a.resize(0);
  a_arcp.resize(0);
  TEST_COMPARE_ARRAYS(a, a_arcp);

#ifdef TEUCHOS_DEBUG
  a_arcp = arcp<T>(n);
  ++a_arcp;
  TEST_THROW(a_arcp.resize(1), std::out_of_range);
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, clear, T )
{
  ArrayRCP<T> a_arcp = arcp<T>(n);
  TEST_EQUALITY( a_arcp.size(), n );
  a_arcp.clear();
  TEST_EQUALITY( a_arcp.size(), 0 );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, nullIterator, T )
{
  typedef ArrayRCP<T> iter_t;
  ArrayRCP<T> arcp1 = Teuchos::NullIteratorTraits<iter_t>::getNull();
  TEST_EQUALITY_CONST(arcp1, Teuchos::null);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, implicitConversions, T )
{

  ECHO(ArrayRCP<T> arcp1 = arcp<T>(n));
  ECHO(ArrayRCP<const T> arcp2 = arcp1);
  TEST_ASSERT(arcp1.shares_resource(arcp2));

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, weakDelete, T )
{

  ECHO(ArrayRCP<T> arcp_strong = arcp<T>(n));

  TEST_EQUALITY_CONST( arcp_strong.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( arcp_strong.is_null(), false );
  TEST_EQUALITY_CONST( arcp_strong.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_strong.weak_count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.total_count(), 1 );

  ECHO(ArrayRCP<T> arcp_weak1 = arcp_strong.create_weak());

  TEST_EQUALITY_CONST( arcp_weak1.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( arcp_weak1.is_null(), false );
  TEST_EQUALITY_CONST( arcp_weak1.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak1.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak1.total_count(), 2 );

  TEST_EQUALITY_CONST( arcp_strong.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_strong.is_null(), false );
  TEST_EQUALITY_CONST( arcp_strong.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp_strong.total_count(), 2 );

  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_strong), true );

  TEST_EQUALITY( arcp_weak1.get(), arcp_weak1.getRawPtr() );
  TEST_EQUALITY( arcp_weak1.get(), arcp_strong.get() );
  TEST_EQUALITY( arcp_weak1.getRawPtr(), arcp_strong.getRawPtr() );

  ECHO(ArrayRCP<T> arcp_weak2 = arcp_weak1);

  TEST_EQUALITY_CONST( arcp_weak2.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( arcp_weak2.is_null(), false );
  TEST_EQUALITY_CONST( arcp_weak2.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak2.weak_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak2.total_count(), 3 );

  TEST_EQUALITY_CONST( arcp_strong.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_strong.is_null(), false );
  TEST_EQUALITY_CONST( arcp_strong.weak_count(), 2 );
  TEST_EQUALITY_CONST( arcp_strong.total_count(), 3 );

  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_strong), true );
  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_weak2), true );
  TEST_EQUALITY_CONST( arcp_weak2.shares_resource(arcp_strong), true );

  TEST_EQUALITY( arcp_weak2.get(), arcp_strong.get() );
  TEST_EQUALITY( arcp_weak2.getRawPtr(), arcp_strong.getRawPtr() );

  ECHO(arcp_strong = null); // This deletes the underlying object of type T!

  TEST_EQUALITY_CONST( arcp_strong.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( arcp_strong.is_null(), true );
  TEST_EQUALITY_CONST( arcp_strong.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.weak_count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.total_count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( arcp_strong.shares_resource(arcp_weak1), false );
  TEST_EQUALITY_CONST( arcp_strong.shares_resource(arcp_weak2), false );

  TEST_EQUALITY_CONST( arcp_weak1.has_ownership(), true );
  TEST_EQUALITY_CONST( arcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.weak_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak1.total_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak1.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( arcp_weak2.has_ownership(), true );
  TEST_EQUALITY_CONST( arcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak2.weak_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak2.total_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak2.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_weak2), true );

  ECHO(arcp_weak1.assert_not_null()); // Does not throw!
  ECHO(arcp_weak2.assert_not_null()); // Does not throw!

  TEST_THROW( arcp_weak1.assert_valid_ptr(), DanglingReferenceError );
#ifdef TEUCHOS_DEBUG
  TEST_THROW( arcp_weak1.operator->(), DanglingReferenceError );
  TEST_THROW( *arcp_weak1, DanglingReferenceError );
  TEST_THROW( arcp_weak1.create_weak(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.get(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.getRawPtr(), DanglingReferenceError );
  TEST_THROW( arcp_weak1[0], DanglingReferenceError );
  TEST_THROW( ++arcp_weak1, DanglingReferenceError );
  TEST_THROW( arcp_weak1++, DanglingReferenceError );
  TEST_THROW( --arcp_weak1, DanglingReferenceError );
  TEST_THROW( arcp_weak1--, DanglingReferenceError );
  TEST_THROW( arcp_weak1+=1, DanglingReferenceError );
  TEST_THROW( arcp_weak1-=1, DanglingReferenceError );
  TEST_THROW( arcp_weak1+1, DanglingReferenceError );
  TEST_THROW( arcp_weak1-1, DanglingReferenceError );
  TEST_THROW( arcp_weak1.getConst(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.persistingView(0,n), DanglingReferenceError );
  TEST_THROW( arcp_weak1.lowerOffset(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.upperOffset(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.size(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.begin(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.end(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.view(0,n), DanglingReferenceError );
  TEST_THROW( arcp_weak1(0,n), DanglingReferenceError );
  TEST_THROW( arcp_weak1(), DanglingReferenceError );
  TEST_THROW( {ArrayView<T> av = arcp_weak1();}, DanglingReferenceError );
  TEST_THROW( {ArrayRCP<const T> ap = getConst(arcp_weak1);},
    DanglingReferenceError );
  TEST_THROW( arcp_weak1.release(), DanglingReferenceError );
#endif // TEUCHOS_DEBUG

  ECHO(arcp_weak1 = null); // Just deicrements weak count!

  TEST_EQUALITY_CONST( arcp_weak1.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( arcp_weak1.is_null(), true );
  TEST_EQUALITY_CONST( arcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.weak_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.total_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( arcp_weak2.has_ownership(), true );
  TEST_EQUALITY_CONST( arcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak2.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak2.total_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak2.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_weak2), false );

  TEST_THROW( arcp_weak2.assert_valid_ptr(), DanglingReferenceError );
#ifdef TEUCHOS_DEBUG
  // ToDo: Fill in
#endif // TEUCHOS_DEBUG

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, danglingArrayView, T )
{
  ArrayView<T> av;
  {
    ArrayRCP<T> arcp1 = arcp<T>(n);
    av = arcp1();
  }
#ifdef TEUCHOS_DEBUG
  TEST_THROW( av.size(), DanglingReferenceError );
  TEST_THROW( av.toString(), DanglingReferenceError );
  TEST_THROW( av.getRawPtr(), DanglingReferenceError );
  TEST_THROW( av[0], DanglingReferenceError );
  TEST_THROW( av.front(), DanglingReferenceError );
  TEST_THROW( av.back(), DanglingReferenceError );
  TEST_THROW( av.view(0, n), DanglingReferenceError );
  TEST_THROW( av(0, n), DanglingReferenceError );
  TEST_THROW( av(), DanglingReferenceError );
  TEST_THROW( av.getConst(), DanglingReferenceError );
  TEST_THROW( av.begin(), DanglingReferenceError );
  TEST_THROW( av.end(), DanglingReferenceError );
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, getRawPtr, T )
{
  ArrayRCP<const T> cptr;
  ArrayRCP<T> ptr;
  TEST_EQUALITY_CONST( getRawPtr(cptr), (const T*)NULL );
  TEST_EQUALITY_CONST( getRawPtr(ptr), (T*)NULL );
  cptr = arcp<T>(n);
  ptr  = arcp<T>(n);
  TEST_EQUALITY( getRawPtr(cptr), &cptr[0]);
  TEST_EQUALITY( getRawPtr(ptr),  &ptr[0] );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CPtr, getRawPtr, T )
{
  const T *cptr = NULL;
  T *ptr = NULL;
  TEST_EQUALITY_CONST( getRawPtr(cptr), (const T*)NULL );
  TEST_EQUALITY_CONST( getRawPtr(ptr),  (T*)NULL );
  cptr = new T[n];
  ptr  = new T[n];
  TEST_EQUALITY( getRawPtr(cptr), &cptr[0]);
  TEST_EQUALITY( getRawPtr(ptr),  &ptr[0] );
  delete [] cptr;
  delete [] ptr;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, arcp_zero, T )
{
  ArrayRCP<T> arcp_strong = arcp<T>(0);
  TEST_EQUALITY(arcp_strong.size(), as<Ordinal>(0));
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, arcpFromArrayView, T )
{
  Array<T> a = generateArray<T>(n);
  ArrayView<T> av = a;
  ArrayRCP<T> arcp1 = Teuchos::arcpFromArrayView(av);
  TEST_COMPARE_ARRAYS(arcp1, av);
}


#ifdef TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, arcp_neg, T )
{
  TEST_THROW(ArrayRCP<T> arcp_strong = arcp<T>(-1),
    std::out_of_range);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, outOfBounds, T )
{
  ECHO(ArrayRCP<T> arcp1 = arcp<T>(n));
  TEST_THROW(arcp1(-1,n), RangeError);
  TEST_THROW(arcp1(0,n+1), RangeError);
  TEST_THROW(arcp1(0,-1), RangeError);
  TEST_THROW(arcp1.view(-1,n), RangeError);
  TEST_THROW(arcp1.view(0,n+1), RangeError);
  TEST_THROW(arcp1.view(0,-1), RangeError);
  TEST_THROW(arcp1.persistingView(-1,n), RangeError);
  TEST_THROW(arcp1.persistingView(0,n+1), RangeError);
  TEST_THROW(arcp1.persistingView(0,-1), RangeError);
}


#endif // TEUCHOS_DEBUG


//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, arcp_neg, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, outOfBounds, T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, construct_n, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, assignSelf, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, assign_n_val, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, assign_begin_end, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, print_iterators, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, deepCopy, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, resize, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, clear, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, nullIterator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, implicitConversions, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, weakDelete, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, danglingArrayView, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, getRawPtr, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CPtr, getRawPtr, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, arcp_zero, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, arcpFromArrayView, T ) \
  DEBUG_UNIT_TEST_GROUP(T)


UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(double)
UNIT_TEST_GROUP(float)


} // namespace
