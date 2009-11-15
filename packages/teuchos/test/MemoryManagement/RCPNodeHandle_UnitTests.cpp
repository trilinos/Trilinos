#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCPNode.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "TestClasses.hpp"


namespace Teuchos {


template<class T>
class MockRCP {
public:
  T* access_private_ptr() const
    {
      return (T*)(0x777777); // Just some bogus address printed to out
    }
};


template<class T>
RCPNode* basicRCPNodeNoAlloc(T* p, const bool has_ownership)
{
  RCPNodeTmpl<T,DeallocDelete<T> > *rcpNode = 
    new RCPNodeTmpl<T,DeallocDelete<T> >(p, DeallocDelete<T>(), has_ownership);
  return rcpNode;
}


template<class T>
RCPNode* basicRCPNode(const bool has_ownership, T **p_out = 0)
{
  T *p = new T;
  if (p_out)
    *p_out = p;
  RCPNode *rcpNode = basicRCPNodeNoAlloc<T>(p, has_ownership);
  return rcpNode;
}


void deleteRCPNode( RCPNode **node )
{
  TEUCHOS_ASSERT(node);
  TEUCHOS_ASSERT(*node);
  (*node)->delete_obj();
  delete (*node);
  *node = 0;
}


template<class T>
RCPNodeHandle basicRCPNodeHandle(const bool has_ownership, T **p_out = 0)
{
  using Teuchos::typeName;
  using Teuchos::concreteTypeName;
  T *p = 0;
  RCPNode *rcpNode = basicRCPNode(has_ownership, &p);
  if (p_out)
    *p_out = p;
#ifdef TEUCHOS_DEBUG
  return RCPNodeHandle(rcpNode, p, typeName(*p), concreteTypeName(*p), has_ownership);
#else
  return RCPNodeHandle(rcpNode);
#endif
}


TEUCHOS_STATIC_SETUP()
{
  // Set the ordering of the tests to make sure the run in the specified order.
  // This is important because we want to build up tests from basic to more
  // advanced in order to catch errors.
  UnitTestRepository::setTestOrdering(true);
}


//
// Non-templated tests
//


TEUCHOS_UNIT_TEST( RCPNodeHandle, assignSelf )
{
  RCPNodeHandle nodeRef;
  nodeRef = nodeRef;
}


TEUCHOS_UNIT_TEST( RCPNodeHandle, defaultConstruct)
{
  RCPNodeHandle nodeRef;
  TEST_EQUALITY_CONST( nodeRef.count(), 0 );
  TEST_EQUALITY_CONST( nodeRef.has_ownership(), false );
  nodeRef.has_ownership(true);
  TEST_EQUALITY_CONST( nodeRef.has_ownership(), false );
#ifdef TEUCHOS_DEBUG
  TEST_EQUALITY_CONST( nodeRef.get_base_obj_map_key_void_ptr(), static_cast<void*>(0) );
  TEST_THROW({nodeRef.set_extra_data(any(),"", Teuchos::PRE_DESTROY, true);},
    NullReferenceError);
  TEST_THROW({any &a = nodeRef.get_extra_data("int","blob"); (void)a;},
    NullReferenceError);
  TEST_THROW({const any &a = getConst(nodeRef).get_extra_data("int","blob"); (void)a;},
    NullReferenceError);
  TEST_THROW({any *a = nodeRef.get_optional_extra_data("int","blob"); (void)a;},
    NullReferenceError);
  TEST_THROW({const any *a = getConst(nodeRef).get_optional_extra_data("int","blob"); (void)a;},
    NullReferenceError);
#endif // TEUCHOS_DEBUG
}


#ifdef TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST( RCPNodeHandle, add_New_RCPNode_basic )
{
  SET_RCPNODE_TRACING();
  RCPNode *node = basicRCPNode<A>(true);
  const int numActiveNodesBase = numActiveRCPNodes();
  ECHO(add_new_RCPNode(node, "dummy"));
  TEST_EQUALITY(numActiveRCPNodes(), numActiveNodesBase+1);
  ECHO(remove_RCPNode(node));
  TEST_EQUALITY(numActiveRCPNodes(), numActiveNodesBase);
  deleteRCPNode(&node);
}


TEUCHOS_UNIT_TEST( RCPNodeHandle, add_New_RCPNode_add_twice_error )
{
  SET_RCPNODE_TRACING();
  RCPNode *node = basicRCPNode<A>(true);
  const int numActiveNodesBase = numActiveRCPNodes();
  ECHO(add_new_RCPNode(node, "dummy"));
  TEST_EQUALITY(numActiveRCPNodes(), numActiveNodesBase+1);
  TEST_THROW(add_new_RCPNode(node, "dummy"), DuplicateOwningRCPError);
  ECHO(remove_RCPNode(node));
  TEST_EQUALITY(numActiveRCPNodes(), numActiveNodesBase);
  deleteRCPNode(&node);
}


TEUCHOS_UNIT_TEST( RCPNodeHandle, add_New_RCPNode_add_two_nodes_same_obj )
{
  SET_RCPNODE_TRACING();
  ECHO(C *c_ptr = new C);
  ECHO(RCPNode *node_c = basicRCPNodeNoAlloc<C>(c_ptr, true));
  ECHO(RCPNode *node_b1 = basicRCPNodeNoAlloc<B1>(c_ptr, true));
  ECHO(const int numActiveNodesBase = numActiveRCPNodes());
  ECHO(add_new_RCPNode(node_c, "dummy"));
  TEST_EQUALITY(numActiveRCPNodes(), numActiveNodesBase+1);
#ifdef HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR
  // We can detect that these are the same object!
  TEST_THROW(add_new_RCPNode(node_b1, "dummy"), DuplicateOwningRCPError);
#else
  // We can not detect if these are the same object!
  ECHO(add_new_RCPNode(node_b1, "dummy"));
  TEST_EQUALITY(numActiveRCPNodes(), numActiveNodesBase+2);
  ECHO(remove_RCPNode(node_b1));
#endif
  TEST_EQUALITY(numActiveRCPNodes(), numActiveNodesBase+1);
  ECHO(remove_RCPNode(node_c));
  TEST_EQUALITY(numActiveRCPNodes(), numActiveNodesBase);
  ECHO(node_b1->has_ownership(false));
  ECHO(deleteRCPNode(&node_b1));
  ECHO(deleteRCPNode(&node_c));
}


#ifdef HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING
TEUCHOS_UNIT_TEST( RCPNodeHandle, remove_RCPNode_missing_node )
{
  SET_RCPNODE_TRACING();
  RCPNode *node = basicRCPNode<A>(true);
  TEST_THROW(remove_RCPNode(node), std::logic_error);
  deleteRCPNode(&node);
}
#endif // HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING


#endif // TEUCHOS_DEBUG


//
// Templated tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, basicConstruct_owns_mem, T )
{
  T *p = 0;
  RCPNodeHandle nodeRef(basicRCPNodeHandle<T>(true, &p));
  TEST_EQUALITY_CONST( nodeRef.count(), 1 );
  TEST_EQUALITY_CONST( nodeRef.has_ownership(), true );
  nodeRef.has_ownership(false);
  TEST_EQUALITY_CONST( nodeRef.has_ownership(), false );
#ifdef TEUCHOS_DEBUG
  TEST_INEQUALITY_CONST( nodeRef.get_base_obj_map_key_void_ptr(), static_cast<void*>(0) );
  TEST_EQUALITY( nodeRef.get_base_obj_map_key_void_ptr(), static_cast<void*>(p) );
  TEST_THROW({any &a = nodeRef.get_extra_data("int","blob"); (void)a;},
    std::invalid_argument);
  TEST_THROW({const any &a = getConst(nodeRef).get_extra_data("int","blob"); (void)a;},
    std::invalid_argument);
#endif // TEUCHOS_DEBUG
  TEST_EQUALITY_CONST(nodeRef.get_optional_extra_data("int","blob"), 0);
  TEST_EQUALITY_CONST(getConst(nodeRef).get_optional_extra_data("int","blob"), 0);
  nodeRef.has_ownership(true);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, basicConstruct_no_owns_mem, T )
{
  RCPNodeHandle nodeRef(basicRCPNodeHandle<T>(false));
  TEST_EQUALITY_CONST( nodeRef.count(), 1 );
  TEST_EQUALITY_CONST( nodeRef.has_ownership(), false );
  nodeRef.has_ownership(true);
  TEST_EQUALITY_CONST( nodeRef.has_ownership(), true );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, weakPtr_basic_1, T )
{

  ECHO(RCPNodeHandle nodeRef1(basicRCPNodeHandle<T>(true)));
  TEST_EQUALITY_CONST( nodeRef1.strength(), RCP_STRONG );

  ECHO(RCPNodeHandle nodeRef2 = nodeRef1.create_weak());

  TEST_EQUALITY_CONST( nodeRef2.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( nodeRef1.strong_count(), 1 );
  TEST_EQUALITY_CONST( nodeRef1.weak_count(), 1 );
  TEST_EQUALITY_CONST( nodeRef2.strong_count(), 1 );
  TEST_EQUALITY_CONST( nodeRef2.weak_count(), 1 );

  ECHO(RCPNodeHandle nodeRef3 = nodeRef2.create_strong());

  TEST_EQUALITY_CONST( nodeRef3.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( nodeRef1.strong_count(), 2 );
  TEST_EQUALITY_CONST( nodeRef1.weak_count(), 1 );
  TEST_EQUALITY_CONST( nodeRef2.strong_count(), 2 );
  TEST_EQUALITY_CONST( nodeRef2.weak_count(), 1 );
  
  MockRCP<T> mockRCP;
  ECHO(nodeRef2.debug_assert_valid_ptr(mockRCP)); // Should not throw!

  // This will make the underlying object T get deleted!
  ECHO(nodeRef1 = null);
  ECHO(nodeRef3 = null);

  TEST_EQUALITY_CONST( nodeRef1.node_ptr()==0, true );
  TEST_EQUALITY_CONST( nodeRef1.is_node_null(), true );
  TEST_EQUALITY_CONST( nodeRef1.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( nodeRef2.node_ptr()!=0, true );
  TEST_EQUALITY_CONST( nodeRef2.is_node_null(), false );
  TEST_EQUALITY_CONST( nodeRef2.is_valid_ptr(), false );

#ifdef TEUCHOS_DEBUG
  TEST_THROW( nodeRef2.debug_assert_valid_ptr(mockRCP),
    DanglingReferenceError );
#endif

  ECHO(nodeRef2 = null); // Should make the underlying node go away

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, weakPtr_basic_2, T )
{

  ECHO(RCPNodeHandle nodeRef1(basicRCPNodeHandle<T>(true)));
  TEST_EQUALITY_CONST( nodeRef1.strength(), RCP_STRONG );

  ECHO(RCPNodeHandle nodeRef2 = nodeRef1.create_weak());
  TEST_EQUALITY_CONST( nodeRef2.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( nodeRef1.strong_count(), 1 );
  TEST_EQUALITY_CONST( nodeRef1.weak_count(), 1 );
  TEST_EQUALITY_CONST( nodeRef2.strong_count(), 1 );
  TEST_EQUALITY_CONST( nodeRef2.weak_count(), 1 );
  
  MockRCP<T> mockRCP;

  ECHO(nodeRef2.debug_assert_valid_ptr(mockRCP)); // Should not throw!

  ECHO(nodeRef2 = null); // The underlying object stays alive!

  TEST_EQUALITY_CONST( nodeRef2.node_ptr()==0, true );
  TEST_EQUALITY_CONST( nodeRef2.is_node_null(), true );
  TEST_EQUALITY_CONST( nodeRef2.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( nodeRef1.node_ptr()!=0, true );
  TEST_EQUALITY_CONST( nodeRef1.is_node_null(), false );
  TEST_EQUALITY_CONST( nodeRef1.is_valid_ptr(), true );

  nodeRef1.debug_assert_valid_ptr(mockRCP); // Should not throw!

}

//
// Test behavior of RCP node tracing but only if it is off by default
//
// NOTE: If node tracing is on by default then we can't control how many nodes
// get created in other code not in the unit test.
//

#if defined(TEUCHOS_DEBUG) && !defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)
#  define DO_RCPNODE_TRACING_TESTS 1
#endif


#ifdef DO_RCPNODE_TRACING_TESTS


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, debugWithNodeTracing, T )
{

  TEST_EQUALITY_CONST(Teuchos::isTracingActiveRCPNodes(), false);
  Teuchos::setTracingActiveRCPNodes(true);
  TEST_EQUALITY_CONST(Teuchos::isTracingActiveRCPNodes(), true);

  {

    T *p = new T; // Never do this in production code!
    const std::string T_name = "T_name";
    const std::string concreateT_name = "concreateT_name";
    const bool has_ownership = true;
    RCPNode *node = new RCPNodeTmpl<T,DeallocDelete<T> >(
      p, DeallocDelete<T>(), has_ownership);

    RCPNodeHandle nodeRef(node, p, T_name, concreateT_name, has_ownership);

    TEST_EQUALITY_CONST(Teuchos::numActiveRCPNodes(), 1);

    out << "\nMake sure output is printed when there is an active node with tracing ...\n";

    const void* rcpNodeKey = getRCPNodeBaseObjMapKeyVoidPtr(p);

    std::ostringstream expendedOutput_oss;
    expendedOutput_oss
        << "\n***"
        << "\n*** Warning! The following Teuchos::RCPNode objects were created but have"
        << "\n*** not been destroyed yet.  This may be an indication that these objects may"
        << "\n*** be involved in a circular dependency!  A memory checking tool may complain"
        << "\n*** that these objects are not destroyed correctly."
        << "\n***\n"
        << "\n"
        << "  0: RCPNode (map_key_void_ptr=" << rcpNodeKey << ")\n"
        << "       Information = {T="<<T_name<<", ConcreteT="<<concreateT_name<<", p="<<p<<", has_ownership="<<has_ownership<<"}\n"
        << "       RCPNode address = " << node << "\n"
        << "       Call number = " << get_add_new_RCPNode_call_number()
        << "\n";

    std::ostringstream printActiveRCPNodes_out;
    Teuchos::printActiveRCPNodes(printActiveRCPNodes_out);
    TEST_EQUALITY( printActiveRCPNodes_out.str(), expendedOutput_oss.str() );

    // NOTE: The above test basically copied and pasted the ouptut stream code
    // from Teuchos::printActiveRCPNodes(...) and will need to be maintained
    // with this code.  However, this is a good test because it ensures that
    // the arguments node, p, T_name, and concreateT_name are passed, stored,
    // and retrieved correctly.  It is also a good test because it ensures
    // that output is printed when node tracing is turned on.
    //
    // This is the very definition of a "white box" test but that is just fine
    // for a unit test.

  }

  TEST_EQUALITY_CONST(Teuchos::numActiveRCPNodes(), 0);

  out << "\nMake sure no output is printed when there are no active nodes ...\n";
  const std::string expendedOutput = "";
  std::ostringstream printActiveRCPNodes_out;
  Teuchos::printActiveRCPNodes(printActiveRCPNodes_out);
  TEST_EQUALITY( printActiveRCPNodes_out.str(), expendedOutput );

  Teuchos::setTracingActiveRCPNodes(false);;
  TEST_EQUALITY_CONST(Teuchos::isTracingActiveRCPNodes(), false);

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, debugWithoutNodeTracing, T )
{

  TEST_EQUALITY_CONST(Teuchos::isTracingActiveRCPNodes(), false);
  Teuchos::setTracingActiveRCPNodes(false);
  TEST_EQUALITY_CONST(Teuchos::isTracingActiveRCPNodes(), false);

  T *p = new T; // Never do this in production code!
  const std::string T_name = "T_name";
  const std::string concreateT_name = "concreateT_name";
  const bool has_ownership = true;
  RCPNode *node = new RCPNodeTmpl<T,DeallocDelete<T> >(
    p, DeallocDelete<T>(), has_ownership);
  
  RCPNodeHandle nodeRef(node, p, T_name, concreateT_name, has_ownership);

  TEST_EQUALITY_CONST(Teuchos::numActiveRCPNodes(), 0);
  
  out << "\nMake sure not output is printed when there is an active node without tracing ...\n";
  const std::string expendedOutput = "";
  std::ostringstream printActiveRCPNodes_out;
  Teuchos::printActiveRCPNodes(printActiveRCPNodes_out);
  TEST_EQUALITY( printActiveRCPNodes_out.str(), expendedOutput );

}


#endif // DO_RCPNODE_TRACING_TESTS


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, copyConstruct, T )
{
  RCPNodeHandle nodeRef1(basicRCPNodeHandle<T>(true));
  RCPNodeHandle nodeRef2(nodeRef1);
  TEST_EQUALITY_CONST( nodeRef1.count(), 2 );
  TEST_EQUALITY_CONST( nodeRef2.count(), 2 );
  TEST_EQUALITY_CONST( nodeRef1.has_ownership(), true );
  TEST_EQUALITY_CONST( nodeRef2.has_ownership(), true );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, assignmentOperator, T )
{
  RCPNodeHandle nodeRef1(basicRCPNodeHandle<T>(true));
  RCPNodeHandle nodeRef2;
  nodeRef2 = nodeRef1;
  TEST_EQUALITY_CONST( nodeRef1.count(), 2 );
  TEST_EQUALITY_CONST( nodeRef2.count(), 2 );
  TEST_EQUALITY_CONST( nodeRef1.has_ownership(), true );
  TEST_EQUALITY_CONST( nodeRef2.has_ownership(), true );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, extraData_basic, T )
{

  RCPNodeHandle nodeRef(basicRCPNodeHandle<T>(true));

  const int v1 = 2;
  const any a1(v1);
  nodeRef.set_extra_data(a1, "a1", Teuchos::PRE_DESTROY, true); 

  any &a2 = nodeRef.get_extra_data(a1.typeName(), "a1");
  TEST_EQUALITY_CONST( a1.same(a2), true );
  TEST_EQUALITY( any_cast<int>(a2), v1 );

  any *a3 = nodeRef.get_optional_extra_data(a1.typeName(), "a1");
  TEST_EQUALITY_CONST( a3!=0, true );
  TEST_EQUALITY( &a2, a3 );
  TEST_EQUALITY_CONST( a3->same(a1), true );

  RCPNodeHandle nodeRef2 = nodeRef;

  const int v2 = 3;
  a2 = v2;
  TEST_EQUALITY( any_cast<int>(a1), v1 );
  TEST_EQUALITY( any_cast<int>(*a3), v2 );

  any &a4 = nodeRef2.get_extra_data(a1.typeName(), "a1");
  TEST_EQUALITY( &a4, &a2 );
  TEST_EQUALITY( &a4, a3 );
  TEST_EQUALITY( any_cast<int>(a4), v2 );
  
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, extraData_basic_const, T )
{

  RCPNodeHandle nodeRef(basicRCPNodeHandle<T>(true));

  const int v1 = 2;
  const any a1(v1);
  nodeRef.set_extra_data(a1, "a1", Teuchos::PRE_DESTROY, true); 
  
  const RCPNodeHandle nodeRef2 = nodeRef;

  const any &a2 = nodeRef2.get_extra_data(a1.typeName(), "a1");
  TEST_EQUALITY_CONST( a1.same(a2), true );
  TEST_EQUALITY( any_cast<int>(a2), v1 );

  const any *a3 = nodeRef2.get_optional_extra_data(a1.typeName(), "a1");
  TEST_EQUALITY_CONST( a3!=0, true );
  TEST_EQUALITY( &a2, a3 );
  TEST_EQUALITY_CONST( a3->same(a1), true );
  
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, extraData_failed, T )
{

  RCPNodeHandle nodeRef(basicRCPNodeHandle<T>(true));

  const int v1 = 2;
  const any a1(v1);
  nodeRef.set_extra_data(a1, "a1", Teuchos::PRE_DESTROY, true); 

#ifdef TEUCHOS_DEBUG

  TEST_THROW({nodeRef.get_extra_data("wrong type", "a1");},
    std::invalid_argument);

  TEST_THROW({nodeRef.get_extra_data(a1.typeName(), "wrong name");},
    std::invalid_argument);

#endif // TEUCHOS_DEBUG

  any *a2 = nodeRef.get_optional_extra_data("wrong type", "a1");
  TEST_EQUALITY_CONST( a2, 0 );

  any *a3 = nodeRef.get_optional_extra_data(a1.typeName(), "wrong name");
  TEST_EQUALITY_CONST( a3, 0 );
  
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCPNodeHandle, extraData_failed_const, T )
{

  RCPNodeHandle nodeRef(basicRCPNodeHandle<T>(true));

  const int v1 = 2;
  const any a1(v1);
  nodeRef.set_extra_data(a1, "a1", Teuchos::PRE_DESTROY, true); 
  
  const RCPNodeHandle nodeRef2 = nodeRef;

#ifdef TEUCHOS_DEBUG

  TEST_THROW({nodeRef2.get_extra_data("wrong type", "a1");},
    std::invalid_argument);

  TEST_THROW({nodeRef2.get_extra_data(a1.typeName(), "wrong name");},
    std::invalid_argument);

#endif // TEUCHOS_DEBUG

  const any *a2 = nodeRef2.get_optional_extra_data("wrong type", "a1");
  TEST_EQUALITY_CONST( a2, 0 );

  const any *a3 = nodeRef2.get_optional_extra_data(a1.typeName(), "wrong name");
  TEST_EQUALITY_CONST( a3, 0 );
  
}


//
// Instantiations
//


#ifdef DO_RCPNODE_TRACING_TESTS

#  define DEBUG_UNIT_TEST_GROUP( T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, debugWithNodeTracing, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, debugWithoutNodeTracing, T )

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, basicConstruct_owns_mem, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, basicConstruct_no_owns_mem, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, weakPtr_basic_1, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, weakPtr_basic_2, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, copyConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, assignmentOperator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, extraData_basic, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, extraData_basic_const, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, extraData_failed, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCPNodeHandle, extraData_failed_const, T ) \
  DEBUG_UNIT_TEST_GROUP(T)


UNIT_TEST_GROUP(A)
//UNIT_TEST_GROUP(B1)
//UNIT_TEST_GROUP(B2)
UNIT_TEST_GROUP(C)
//UNIT_TEST_GROUP(D)
UNIT_TEST_GROUP(E)

// 2008/09/22: rabartl: Above: We don't need to test with all of these classes
// in order to test this functionality.


} // namespace Teuchos 
