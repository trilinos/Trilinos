#ifndef stk_adapt_IAdapter_hpp
#define stk_adapt_IAdapter_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk_classic {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * An IAdapter is an abstract base class that provides all functionality of a Refiner but requires sub-classes to 
     *   overload the refineMethodApply() and buildUnrefineList() methods.  The refineMethodApply() method is used by the driving code during
     *   refinement/parallel communication, while the buildUnrefineList() supplies the driving code with a list of
     *   candidate elements that are requested to be unrefined.  Generally, a user will not overload IAdapter directly,
     *   but instead will use one of its pre-existing subclasses (see below).
     *
     * There are two flavors of derived classes that can be used to develop your own adapter: 
     *  a) overload "mark" method on abstract base class IEdgeAdapter or IElementAdapter
     *  b) or use the predicate-based approach where you supply simple structs with an operator() that supplies
     *       information on whether an edge or element is to be refined/unrefined (classes named *Predicate*.*pp )
     *
     * Details on the two flavors of adapters:
     *
     * a) 
     * IEdgeAdapter is a base class for your class that provide a mark(edge_info) method.
     * IElementAdapter allows overloading mark(element_info) method.
     * In either case, the mark methods tell the driving code to refine any (or all edges) of the mesh, unrefine an
     * element (note that all "siblings" of the element must be marked before the collection is deleted and the parent
     * is reinstated to replace the children), or refine an element by subsequently marking all its edges.
     *
     * @see IEdgeAdapter, IElementAdapter, 
     *   unit_tests/TestLocalRefinerTri_N_3_IEdgeAdapter, 
     *   unit_tests/TestLocalRefinerTri_N_3_IElementAdapter, 
     *   and these tests in  unit_tests/UnitTestLocalRefiner:
     *     STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_3_1_IEdgeAdapter)
     *     STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_3_1_IElementAdapter)
     *   
     * b)
     * Predicate based adapters are similar to Edge/IElementAdapter but are constructed with user-supplied structs that must
     * provide an operator() that takes edge- or element-based info, respectively, depending on the sub-flavor.  
     *
     *   PredicateBasedIElementAdapter takes two structs defining a refine and unrefine operation on an element.
     *   PredicateBasedEdgeAdapter takes structs defining a refine operation on an edge, and unrefine on element.
     *
     *   The structs are required to satisfy signatures of ElementRefinePredicate (for element-based ops) and
     *   IEdgeBasedAdapterPredicate respectively.
     *
     * @see ElementRefinePredicate, IEdgeBasedAdapterPredicate, PredicateBasedEdgeAdapter, PredicateBasedIElementAdapter
     *   and this tests in unit_tests/UnitTestLocalRefiner.cpp:
     *      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_5_ElementBased)
     *
     *   and the various tests in regression_tests/RegressionTestUniformRefiner.cpp.
     *
     *
     * 
     */

    /// Note: these can be bitwise or'ed to say this element can be refined or unrefined which is useful during testing
    ///   (for example, during refinement we may want to refine a region, then immediately unrefine, so elements
    ///    could have both refine and unrefine marks).
    /// It can also be used as an error check that enforces that both options can't be set at the same time.

    enum AdaptInstruction {
      DO_NOTHING  = 0,
      DO_REFINE   = 1 << 0,
      DO_UNREFINE = 1 << 1
    };

    class IAdapter : public Refiner
    {
    public:


      virtual ElementUnrefineCollection  buildUnrefineList() = 0;

    protected:
      IAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk_classic::mesh::FieldBase *proc_rank_field=0);

      virtual void 
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, 
                                              vector<NeededEntityType>& needed_entity_ranks) = 0;


    };

    IAdapter::IAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk_classic::mesh::FieldBase *proc_rank_field) : 
      Refiner(eMesh, bp, proc_rank_field)
    {
    }


  }
}
#endif
