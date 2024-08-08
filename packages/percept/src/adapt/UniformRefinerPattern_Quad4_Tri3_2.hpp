// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Quad4_Tri3_2_hpp
#define adapt_UniformRefinerPattern_Quad4_Tri3_2_hpp


//#include "UniformRefinerPattern.hpp"

  namespace percept {

    //struct Specialization {};

    template <>
    class UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 2  > : public URP<shards::Quadrilateral<4> , shards::Triangle<3> >
    {
    public:

      virtual bool edgeMarkIsEnough() { return false; }
      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Quadrilateral<4> , shards::Triangle<3> >(eMesh)
       {
         EXCEPTWATCH;
         m_primaryEntityRank = eMesh.face_rank();
         if (m_eMesh.get_spatial_dim() == 2)
           m_primaryEntityRank = stk::topology::ELEMENT_RANK;

         setNeededParts(eMesh, block_names, false);

       }

      virtual void doBreak() {}

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(0);
        //needed_entities[0] = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK :  m_eMesh.face_rank());
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 2; }

      virtual size_t estimateNumberOfNewElements(percept::PerceptMesh& eMesh, stk::mesh::EntityRank rank, NodeRegistry& nodeRegistry, size_t num_elem_not_ghost)
      {

        unsigned num_elem_needed = num_elem_not_ghost * this->getNumNewElemPerElem();
        return num_elem_needed;
      }

      virtual StringStringMap fixSurfaceAndEdgeSetNamesMap()
      {
        StringStringMap str_map;
        str_map["wedge6"] = "tet4";
        str_map["pyramid5"] = "tet4";
        str_map["hex8"] = "tet4";
        str_map["quad4"] = "tri3";
        return str_map;
      }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        typedef std::array<stk::mesh::EntityId, 3> tri_tuple_type;
        static vector<tri_tuple_type> elems(2);

        CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

        add_parts = m_toParts;
        {
          unsigned globalIqf  = VERT_N(0);
          unsigned minVal     = globalIqf;
          unsigned indxMinVal = 0;
          for (unsigned iFaceNodeOrd=1; iFaceNodeOrd < 4; iFaceNodeOrd++)
            {
              globalIqf = VERT_N(iFaceNodeOrd);
              if (globalIqf < minVal)
                {
                  minVal = globalIqf;
                  indxMinVal = iFaceNodeOrd;
                }
            }

          unsigned istart = indxMinVal;
          elems[0] = {VERT_N((0 + istart) % 4), VERT_N((1 + istart) % 4), VERT_N((2 + istart) % 4)};
          elems[1] = {VERT_N((0 + istart) % 4), VERT_N((2 + istart) % 4), VERT_N((3 + istart) % 4)};
        }

        bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity newElement = stk::mesh::Entity();
            if (!use_declare_element_side)
              newElement = *element_pool;

            // 3 nodes of the new quads
            stk::mesh::Entity nodes[3] = {
              eMesh.createOrGetNode(elems[ielem][0]),
              eMesh.createOrGetNode(elems[ielem][1]),
              eMesh.createOrGetNode(elems[ielem][2]) };

            create_side_element(eMesh, use_declare_element_side, nodes, 3, newElement);

            if (proc_rank_field && m_eMesh.entity_rank(element) == stk::topology::ELEMENT_RANK)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            change_entity_parts(eMesh, element, newElement);

            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            std::vector<stk::mesh::Entity> elements(1,element);
            eMesh.prolongateElementFields( elements, newElement);

            ft_element_pool++;
            if (!use_declare_element_side)
              element_pool++;

          }

      }

    };

  }

#endif
