// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_QualityElementAdapter_hpp
#define adapt_QualityElementAdapter_hpp

#include <functional>
#include <cmath>

#include <adapt/PredicateBasedElementAdapter.hpp>

  namespace percept {


    template<class RefinePredicate>
    class QualityElementAdapter : public PredicateBasedElementAdapter<RefinePredicate>
    {

      // parameters:
      double m_criterion_bad_element_quality;
      double m_criterion_bad_edge_quality;

      double m_criterion_bad_transition_element_quality;
      double m_criterion_bad_transition_edge_quality;

      typedef PredicateBasedElementAdapter<RefinePredicate> Base;
    public:

      QualityElementAdapter(RefinePredicate& predicate_refine,
                           percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0) :
        Base(predicate_refine, eMesh, bp, proc_rank_field)
      {
        // note, angles only valid for triangle quality; however, they give a reasonable way to compute ratios for 3D as well
        double bad_angle = 30.0;  // degrees
        m_criterion_bad_element_quality = 1./(2.*std::sin(0.5*bad_angle*M_PI/180.0));  // ratio of max-edge to min-edge is > than this value, it's a bad element
        m_criterion_bad_edge_quality = m_criterion_bad_element_quality;     // ratio of local edge to min-edge > than this value, it's a bad edge

        // we can have different transition element requirements, here relaxed
        double bad_angle_transition_elements = 20.0;  // degrees
        m_criterion_bad_transition_element_quality = 1./(2.*std::sin(0.5*bad_angle_transition_elements*M_PI/180.0));
        m_criterion_bad_transition_edge_quality = m_criterion_bad_transition_element_quality;

        //std::cout << "m_criterion_bad_transition_element_quality= " << m_criterion_bad_transition_element_quality
        //          << " m_criterion_bad_element_quality = " << m_criterion_bad_element_quality << std::endl;
      }

      /// set quality criteria using angles (in degrees) - if primary (marked) elements
      /// have angle less than @param bad_angle, or if secondary (transition) elements
      /// have angle less than @param bad_angle_transition_elements,try to improve
      /// by different marking strategy (mark longest edges)
      void setQualityCriteriaAngle(double bad_angle, double bad_angle_transition_elements)
      {
        m_criterion_bad_element_quality = 1./(2.*std::sin(0.5*bad_angle*M_PI/180.0));  // ratio of max-edge to min-edge is > than this value, it's a bad element
        m_criterion_bad_edge_quality = m_criterion_bad_element_quality;     // ratio of local edge to min-edge > than this value, it's a bad edge

        m_criterion_bad_transition_element_quality = 1./(2.*std::sin(0.5*bad_angle_transition_elements*M_PI/180.0));
        m_criterion_bad_transition_edge_quality = m_criterion_bad_transition_element_quality;
      }
      /// set quality criteria directly - criteria are maximum allowed ratio of max to min edge lengths
      /// if primary/marked elements (@param criterion_bad_element_quality) or secondary/transition
      /// (@param criterion_bad_transition_element_quality) have max/min edge length greater than
      /// that specified, try to improve quality with different marking
      void setQualityCriteria(double criterion_bad_element_quality, double criterion_bad_transition_element_quality)
      {
        m_criterion_bad_element_quality = criterion_bad_element_quality;
        m_criterion_bad_edge_quality = m_criterion_bad_element_quality;

        m_criterion_bad_transition_element_quality = criterion_bad_transition_element_quality;
        m_criterion_bad_transition_edge_quality = m_criterion_bad_transition_element_quality;
      }

      /// fine-grained control - advanced use - see code
      void setQualityCriteriaAll(double criterion_bad_element_quality, double criterion_bad_edge_quality,
                                 double criterion_bad_transition_element_quality, double criterion_bad_transition_edge_quality)
      {
        m_criterion_bad_element_quality = criterion_bad_element_quality;
        m_criterion_bad_edge_quality    = criterion_bad_edge_quality;

        m_criterion_bad_transition_element_quality = criterion_bad_transition_element_quality;
        m_criterion_bad_transition_edge_quality    = criterion_bad_transition_edge_quality;
      }

    protected:

      virtual void
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
            vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data)
      {
        const CellTopologyData * const cell_topo_data = Base::m_eMesh.get_cell_topology(element);
        CellTopology cell_topo(cell_topo_data);

        bool markInfo = (Base::m_predicate_refine(element) & DO_REFINE);

        for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
          {
            unsigned numSubDimNeededEntities = 0;
            stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

            if (needed_entity_rank == Base::m_eMesh.edge_rank())
              {
                numSubDimNeededEntities = cell_topo_data->edge_count;
              }
            else if (needed_entity_rank == Base::m_eMesh.face_rank())
              {
                numSubDimNeededEntities = cell_topo_data->side_count;
              }
            else if (needed_entity_rank == Base::m_eMesh.element_rank())
              {
                numSubDimNeededEntities = 1;
              }

            // setup to mark all edges
            std::vector<bool> markInfoVec(numSubDimNeededEntities, markInfo);

            setMarkFromQuality(element, needed_entity_rank, markInfo, markInfoVec);

            for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
              {
                (Base::m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, markInfoVec[iSubDimOrd],bucket_topo_data);
              }
          }
      }


      void setMarkAllEdgesFromQuality(stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, bool markInfo, std::vector<bool>& markInfoVec,
                                      double criterion_bad_element_quality, double criterion_bad_edge_quality)
      {

        // Alg A: find longest edge, if not marked, mark it
        // Alg B: mark all edges longer than X-factor * min edge

        if (needed_entity_rank != Base::m_eMesh.edge_rank())
          return;

        const CellTopologyData * const cell_topo_data = Base::m_eMesh.get_cell_topology(element);
        int spatialDimension = Base::m_eMesh.get_spatial_dim();
        CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (Base::m_eMesh, element, Base::m_eMesh.node_rank());
        //VERIFY_OP_ON(elem_nodes.size(), ==, 3, "only for tris");

        stk::mesh::FieldBase* coordField = Base::m_eMesh.get_coordinates_field();

        unsigned numSubDimNeededEntities = 0;
        numSubDimNeededEntities = cell_topo_data->edge_count;

        //int iedge_max = -1;
        double edge_len_max = 0.0;
        double edge_len_min = std::numeric_limits<double>::max();
        double edge_lengths[numSubDimNeededEntities];
        //bool marked[numSubDimNeededEntities];
        for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
          {
            stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
            stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
            double * const coord0 = Base::m_eMesh.field_data( *coordField , node0 );
            double * const coord1 = Base::m_eMesh.field_data( *coordField , node1 );
            double edge_len = 0.0;
            for (int i=0; i < spatialDimension; ++i)
              {
                edge_len += (coord0[i] - coord1[i])*(coord0[i] - coord1[i]);
              }
            edge_len = std::sqrt(edge_len);

            // if the edge is marked, the new element will have half edge length
            bool is_marked = markInfo & DO_REFINE;
            bool is_not_marked = false;
            if (!is_marked)
              {
                Base::getEdgeMarkInfo(element, needed_entity_rank, iSubDimOrd, is_marked, is_not_marked);
              }
            //marked[iSubDimOrd] = is_marked;
            if (is_marked)
              {
                edge_len *= 0.5;
              }

            if (edge_len > edge_len_max)
              {
                edge_len_max = edge_len;
                //iedge_max = iSubDimOrd;
              }
            edge_len_min = std::min(edge_len_min, edge_len);
            edge_lengths[iSubDimOrd] = edge_len;
          }

        if (edge_len_max/edge_len_min < criterion_bad_element_quality)
          {
            if (markInfo & DO_REFINE)
              {
                for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                  {
                    markInfoVec[iSubDimOrd] = true;
                  }
              }
            // else, it's a transition element, predicted quality is ok, no additional marks needed
          }
        else
          {
            for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
              {
                if (edge_lengths[iSubDimOrd] / edge_len_min >= criterion_bad_edge_quality)
                  {
                    markInfoVec[iSubDimOrd] = true;
                  }
                else
                  {
                    markInfoVec[iSubDimOrd] = false;
                  }
              }
          }
      }

      void setMarkFromQuality(stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, bool markInfo, std::vector<bool>& markInfoVec)
      {
        // Alg A: find longest edge, if not marked, mark it
        // Alg B: mark all edges longer than X-factor * min edge
        //
        // now using Alg B
        //
        bool useAlgB = true;

        if (needed_entity_rank != Base::m_eMesh.edge_rank())
          return;

        // if it is marked for refine, and thus not a transition element, mark edges based on quality
        if (markInfo & DO_REFINE)
          {
            setMarkAllEdgesFromQuality(element, needed_entity_rank, markInfo, markInfoVec, m_criterion_bad_element_quality, m_criterion_bad_edge_quality);
            return;
          }

        // it's a possible transition element...
        const CellTopologyData * const cell_topo_data = Base::m_eMesh.get_cell_topology(element);
        int spatialDimension = Base::m_eMesh.get_spatial_dim();
        CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (Base::m_eMesh, element,  Base::m_eMesh.node_rank());
        //VERIFY_OP_ON(elem_nodes.size(), ==, 3, "only for tris");

        stk::mesh::FieldBase* coordField = Base::m_eMesh.get_coordinates_field();

        unsigned numSubDimNeededEntities = 0;
        numSubDimNeededEntities = cell_topo_data->edge_count;

        bool at_least_one_edge_marked = false;
        for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
          {
            bool is_marked = false;
            bool is_not_marked = false;
            Base::getEdgeMarkInfo(element, needed_entity_rank, iSubDimOrd, is_marked, is_not_marked);
            if (is_marked)
              {
                at_least_one_edge_marked = true;
                break;
              }
          }
        // not a transition element
        if (!at_least_one_edge_marked)
          return;

        if (useAlgB)
          {
            setMarkAllEdgesFromQuality(element, needed_entity_rank, markInfo, markInfoVec,
                                       m_criterion_bad_transition_element_quality,
                                       m_criterion_bad_transition_edge_quality);
            return;
          }

        int iedge_max = -1;
        double edge_len_max = 0.0;
        for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
          {
            stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
            stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
            double * const coord0 = Base::m_eMesh.field_data( *coordField , node0 );
            double * const coord1 = Base::m_eMesh.field_data( *coordField , node1 );
            double edge_len = 0.0;
            for (int i=0; i < spatialDimension; ++i)
              {
                edge_len += (coord0[i] - coord1[i])*(coord0[i] - coord1[i]);
              }
            edge_len = std::sqrt(edge_len);
            if (edge_len > edge_len_max)
              {
                edge_len_max = edge_len;
                iedge_max = iSubDimOrd;
              }
          }

        bool is_marked = false;
        bool is_not_marked = false;
        // FIXME - is it necessary to get existing marks? can't we just mark it regardless?
        Base::getEdgeMarkInfo(element, needed_entity_rank, (unsigned)iedge_max, is_marked, is_not_marked);
        //if (!is_marked && is_not_marked)
        if (!is_marked)
          {
            //std::cout << "P[" << m_eMesh.get_rank() << "] setMarkFromQuality:: marking an unmarked edge" << std::endl;
            markInfoVec[iedge_max] = true;
          }

      }

    };



  }

#endif
