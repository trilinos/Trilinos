// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Line2_Line3_1_sierra_hpp
#define adapt_UniformRefinerPattern_Line2_Line3_1_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <percept/PerceptBoostArray.hpp>


  namespace percept {

    template <>
    class UniformRefinerPattern< shards::Line<2>, shards::Line<3>, 1, SierraPort > : public URP<shards::Line<2> , shards::Line<3> >
    {

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Line<2> , shards::Line<3> >(eMesh)
      {
        m_primaryEntityRank = m_eMesh.edge_rank();
        if (m_eMesh.get_spatial_dim() == 1)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, false);
        Elem::StdMeshObjTopologies::bootstrap();

      }


      virtual void doBreak() override {}

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
      {
        EXCEPTWATCH;
        bp.resize(1);

        if (eMesh.get_spatial_dim() == 1)
          {
            bp[0] = this;
          }
        else
          {
            bp.resize(0);
          }

      }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(1);
        //needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[0].first = (m_eMesh.get_spatial_dim() == 1 ? stk::topology::ELEMENT_RANK : m_eMesh.edge_rank());
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 1; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        typedef std::array<stk::mesh::EntityId,3> quadratic_type;
        static vector<quadratic_type> elems(1);

        CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

        add_parts = m_toParts;

#if 0
        double coord_x[3];
        for (int iedge = 0; iedge < 1; iedge++)
          {
            double * mp = midPoint(VERT_COORD(0), VERT_COORD(1), eMesh.get_spatial_dim(), coord_x);
            //double * mp = midPoint(EDGE_COORD(iedge,0), EDGE_COORD(iedge,1), eMesh.get_spatial_dim(), coord_x);

            if (!EDGE_N(iedge))
              {
                std::cout << "P[" << eMesh.get_rank() << " nid ## = 0  " << std::endl;
              }
            eMesh.createOrGetNode(EDGE_N(iedge), mp);
          }
        // FIXME
        nodeRegistry.prolongateCoords(element, m_primaryEntityRank, 0u);
        nodeRegistry.prolongateFields(element, m_primaryEntityRank, 0u);
#endif
        nodeRegistry.addToExistingParts(element, m_primaryEntityRank, 0u);

#define CENTROID_N NN(m_primaryEntityRank,0)


        {
          quadratic_type& EN = elems[0];

          for (unsigned ind = 0; ind < 2; ind++)
            {
              stk::mesh::EntityId inode = VERT_N(ind);
              EN[ind] = inode;
            }

          EN[2] = CENTROID_N;
        }

#undef CENTROID_N

        // write a diagram of the refinement pattern as a vtk file, or a latex/tikz/pgf file
#define WRITE_DIAGRAM 0
#if WRITE_DIAGRAM

#endif

        bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity newElement = stk::mesh::Entity();
            if (!use_declare_element_side)
              newElement = *element_pool;

            stk::mesh::Entity nodes[3];
            for (int inode=0; inode < 3; inode++)
              {
                stk::mesh::EntityId eid = elems[ielem][inode];
                stk::mesh::Entity node = eMesh.createOrGetNode(eid);
                nodes[inode] = node;
              }
            create_side_element(eMesh, use_declare_element_side, nodes, 3, newElement);

            // FIXME
            if (0 && proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            change_entity_parts(eMesh, element, newElement);


            {
              if (!elems[ielem][0])
                {
                  std::cout << "P[" << eMesh.get_rank() << " nid = 0  " << std::endl;
                  exit(1);
                }

            }

            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            ft_element_pool++;
            if (!use_declare_element_side)
              element_pool++;

          }

      }

    };

  }

#endif
