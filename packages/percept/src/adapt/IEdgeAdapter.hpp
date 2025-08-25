// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_IEdgeAdapter_hpp
#define adapt_IEdgeAdapter_hpp

#include <adapt/IAdapter.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * An IEdgeAdapter is an abstract base class for derived classes that are required to overload the mark method,
     *    which provides info such as the element the edge belongs to, which edge ordinal it is, the nodes of the edge
     *    and the edge coordinates.
     */
    class IEdgeAdapter : public IAdapter
    {
    public:
      IEdgeAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0)
      : IAdapter(eMesh, bp, proc_rank_field) {}

      /// can be overriden
      virtual void  buildUnrefineList(ElementUnrefineCollection&) override;

    public:

      /// Client supplies these methods - given an element, which edge, and the nodes on the edge, return instruction on what to do to the edge,
      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE

      virtual int markEdge(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
                           double *coord0, double *coord1, std::vector<int>* existing_edge_marks) = 0;

      /// This convenience method calls mark and if all edges are marked for unrefine, it returns DO_UNREFINE to unrefine the element.
      /// This method can be overriden to allow for an "element-based" determination that doesn't need to visit edges.
      virtual int markUnrefine(const stk::mesh::Entity element);

      /// This convenience method calls mark and and returns how many edges are marked for refine
      int markCountRefinedEdges(const stk::mesh::Entity element);

      virtual void
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
                                              vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data) override;


    };
  }

#endif
