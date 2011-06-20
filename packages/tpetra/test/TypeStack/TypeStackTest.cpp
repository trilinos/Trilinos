#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TypeTraits.hpp>

#include "TpetraExt_TypeStack.hpp"

template <class TS>
void recurseTypes(Teuchos::FancyOStream &os) 
{
  using Teuchos::TypeNameTraits;
  //
  typedef typename TS::type type;
  os << TypeNameTraits<type>::name();
  if (TS::bottom) {
    os << std::endl;
  }
  else {
    os << " -> ";
    recurseTypes<typename TS::next>(os);
  }
}

TEUCHOS_UNIT_TEST( TypeStack, FullStack4 ) {
  // test the type stack
  TPETRAEXT_TYPESTACK4(FullStack,double,
                       float,
                       int,
                       short)
  TEST_EQUALITY_CONST( FullStack::height, 4 );
  recurseTypes<FullStack>(out);
  // test that the macro is equivalent to the manual instantiation
  using Tpetra::Ext::TypeStack;
  typedef TypeStack< double, 
          TypeStack< float, 
          TypeStack< int , 
                     short  > > > 
                     FullStackManual;
  const bool same = Teuchos::TypeTraits::is_same<FullStack,FullStackManual>::value;
  TEST_EQUALITY_CONST( same, true );
}

TEUCHOS_UNIT_TEST( TypeStack, FullStack3 ) {
  // test the type stack
  TPETRAEXT_TYPESTACK3(FullStack,double,
                       float,
                       int)
  TEST_EQUALITY_CONST( FullStack::height, 3 );
  recurseTypes<FullStack>(out);
  // test that the macro is equivalent to the manual instantiation
  using Tpetra::Ext::TypeStack;
  typedef TypeStack< double, 
          TypeStack< float, 
                     int > >
                     FullStackManual;
  const bool same = Teuchos::TypeTraits::is_same<FullStack,FullStackManual>::value;
  TEST_EQUALITY_CONST( same, true );
}

TEUCHOS_UNIT_TEST( TypeStack, FullStack2 ) {
  // test the type stack
  TPETRAEXT_TYPESTACK2(FullStack,double,
                       float)
  TEST_EQUALITY_CONST( FullStack::height, 2 );
  recurseTypes<FullStack>(out);
  // test that the macro is equivalent to the manual instantiation
  using Tpetra::Ext::TypeStack;
  typedef TypeStack< double, 
                     float >
                     FullStackManual;
  const bool same = Teuchos::TypeTraits::is_same<FullStack,FullStackManual>::value;
  TEST_EQUALITY_CONST( same, true );
}

TEUCHOS_UNIT_TEST( TypeStack, FullStack1 ) {
  // test a boring type stack of one type
  // unlike above, this isn't actually a TypeStack object
  TPETRAEXT_TYPESTACK1(FullStack, int)
  TEST_EQUALITY_CONST( FullStack::height, 1 );
  recurseTypes<FullStack>(out);
  // test that the macro is equivalent to the manual instantiation
  using Tpetra::Ext::TypeStackBottom;
  typedef TypeStackBottom<int> 
                          FullStackManual;
  const bool same = Teuchos::TypeTraits::is_same<FullStack,FullStackManual>::value;
  TEST_EQUALITY_CONST( same, true );
}
