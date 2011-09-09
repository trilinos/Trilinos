#ifndef stk_adapt_Marker_hpp
#define stk_adapt_Marker_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A Marker is an abstract base class that provides all functionality of a Refiner but requires sub-classes to 
     *   overload the apply() and buildUnrefineList() methods.  The apply() method is used by the driving code during
     *   refinement/parallel communication, while the buildUnrefineList() supplies the driving code with a list of
     *   candidate elements that are requested to be unrefined.  Generally, a user will not overload Marker directly,
     *   but instead will use one of its pre-existing subclasses (see below).
     *
     * There are two flavors of derived classes that can be used to develop your own marker: 
     *  a) overload "mark" method on abstract base class EdgeMarker or ElementMarker
     *  b) or use the predicate-based approach where you supply simple structs with an operator() that supplies
     *       information on whether an edge or element is to be refined/unrefined (classes named *Predicate*.*pp )
     *
     * Details on the two flavors of markers:
     *
     * a) 
     * EdgeMarker is a base class for your class that provide a mark(edge_info) and markUnrefine(element_info) methods.
     * ElementMarker allows overloading mark(element_info) and markUnrefine(element_info) methods.
     * In either case, the mark methods tell the driving code to refine any (or all edges) of the mesh, unrefine an
     * element (note that all "siblings" of the element must be marked before the collection is deleted and the parent
     * is reinstated to replace the children), or refine an element by subsequently marking all its edges.
     *
     * @see EdgeMarker, ElementMarker, 
     *   unit_tests/TestLocalRefinerTri_N_3_EdgeMarker, 
     *   unit_tests/TestLocalRefinerTri_N_3_ElementMarker, 
     *   and these tests in  unit_tests/UnitTestLocalRefiner:
     *     STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_3_1_EM)
     *     STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_3_1_ElementMarker)
     *   
     * b)
     * Predicate based markers are similar to Edge/ElementMarker but are constructed with user-supplied struct that must
     * provide an operator() that takes an edge or element info, depending on the sub-flavor.  
     *
     *   PredicateBasedMarker takes two structs defining a refine and unrefine operation on an element, and 
     *   Similarly, PredicateBasedEdgeMarker takes structs defining a refine operation on an edge, and unrefine on element.
     *
     *   The structs are required to satisfy signatures of ElementFieldBasedRefinePredicate (for element-based ops) and
     *   EdgeBasedMarkerPredicate respectively.
     *
     * @see ElementFieldBasedRefinePredicate, EdgeBasedMarkerPredicate, PredicateBasedEdgeMarker, PredicateBasedMarker
     *   and this tests in unit_tests/UnitTestLocalRefiner.cpp:
     *      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_5_FieldBased)
     *
     *   and the various tests in regression_tests/RegressionTestUniformRefiner.cpp.
     *
     *
     * 
     */
    class Marker : public Refiner
    {
    public:
      Marker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      virtual ElementUnrefineCollection  buildUnrefineList() = 0;

    protected:

      virtual void 
      apply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, 
                                              vector<NeededEntityType>& needed_entity_ranks) = 0;


    };

    Marker::Marker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) : 
      Refiner(eMesh, bp, proc_rank_field)
    {
    }


  }
}
#endif
