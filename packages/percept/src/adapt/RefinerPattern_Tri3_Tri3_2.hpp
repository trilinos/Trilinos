// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Tri3_Tri3_2_sierra_hpp
#define adapt_RefinerPattern_Tri3_Tri3_2_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"

  namespace percept {

    /// this is for testing only - or could be unsed in future for a bisection-based refinement scheme

    template <>
    class RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 2 > : public URP<shards::Triangle<3>,shards::Triangle<3>  >
    {

      UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > * m_edge_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Triangle<3>, shards::Triangle<3>  >(eMesh),
                                                                                                    m_edge_breaker(0)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        if (m_eMesh.get_spatial_dim() == 2)
          {
            m_edge_breaker =  new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) ;
          }

      }

      ~RefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
      }


      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
      {
        EXCEPTWATCH;
        bp.resize(2);

        if (eMesh.get_spatial_dim() == 2)
          {
            bp[0] = this;
            bp[1] = m_edge_breaker;
          }
        else
          {
            bp.resize(0);
          }
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[0].second = 1u;
      }

      virtual unsigned getNumNewElemPerElem() override { return 2; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& /*nodeRegistry*/,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        typedef std::array<stk::mesh::EntityId, 3> tri_tuple_type;
        static vector<tri_tuple_type> elems(2);

        shards::CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;
        add_parts = m_toParts;


        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 3; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }
        if (num_edges_marked > 1)
          {
            throw std::runtime_error("RefinerPattern_Tri3_Tri3_2 can only refine element with one marked edge");
          }
        if (num_edges_marked == 0)
          return;

        double tmp_x[3];
        for (int iedge = 0; iedge < 3; iedge++)
          {
            double * mp = midPoint(EDGE_COORD(iedge,0), EDGE_COORD(iedge,1), eMesh.get_spatial_dim(), tmp_x);

            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                eMesh.createOrGetNode(EDGE_N(iedge), mp);
                elems[0] = {VERT_N(iedge), EDGE_N(iedge), VERT_N((iedge+2)%3) };
                elems[1] = {EDGE_N(iedge), VERT_N((iedge+1)%3), VERT_N((iedge+2)%3) };
              }
          }

        //nodeRegistry.prolongateCoords(*const_cast<stk::mesh::Entity>(&element), stk::topology::ELEMENT_RANK, 0u);


        bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity newElement = stk::mesh::Entity();
            if (!use_declare_element_side)
              newElement = *element_pool;
            //std::cout << "tmp newElement id = " << m_eMesh.identifier(newElement) << std::endl;

            stk::mesh::Entity nodes[3] = {eMesh.createOrGetNode(elems[ielem][0]),
                                          eMesh.createOrGetNode(elems[ielem][1]),
                                          eMesh.createOrGetNode(elems[ielem][2])};
            create_side_element(eMesh, use_declare_element_side, nodes, 3, newElement);

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(m_eMesh.get_rank());
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            eMesh.get_bulk_data()->change_entity_parts( newElement, add_parts, remove_parts );

            // 3 nodes of the new tris
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem][0]), 0);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem][1]), 1);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem][2]), 2);

            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            ft_element_pool++;
            if (!use_declare_element_side)
              element_pool++;

          }


      }

    };

  }

#endif
