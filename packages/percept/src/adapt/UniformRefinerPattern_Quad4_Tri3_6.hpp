// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Quad4_Tri3_6_hpp
#define adapt_UniformRefinerPattern_Quad4_Tri3_6_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > : public URP<shards::Quadrilateral<4> , shards::Triangle<3> >
    {
#define EDGE_BREAKER_Q4_T3_6_S 1
#if EDGE_BREAKER_Q4_T3_6_S
      UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > * m_edge_breaker;
#endif

    public:

      virtual bool edgeMarkIsEnough() { return false; }
      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Quadrilateral<4> , shards::Triangle<3> >(eMesh)
       {
         m_primaryEntityRank = m_eMesh.face_rank();
         if (m_eMesh.get_spatial_dim() == 2)
           m_primaryEntityRank = stk::topology::ELEMENT_RANK;

         setNeededParts(eMesh, block_names, false);

        Elem::StdMeshObjTopologies::bootstrap();

#if EDGE_BREAKER_Q4_T3_6_S
        if (m_eMesh.get_spatial_dim() == 2)
          m_edge_breaker =  new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) ;
        else
          m_edge_breaker = 0;
#endif

       }

      ~UniformRefinerPattern()
      {
#if EDGE_BREAKER_Q4_T3_6_S
        if (m_edge_breaker) delete m_edge_breaker;
#endif
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp.resize(2);

        if (eMesh.get_spatial_dim() == 2)
          {
            bp[0] = this;
#if EDGE_BREAKER_Q4_T3_6_S
            bp[1] = m_edge_breaker;
#endif
          }
        else if (eMesh.get_spatial_dim() == 3)
          {
            bp.resize(0);
            // FIXME
            //             std::cout << "ERROR" ;
            //             exit(1);
          }

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank(); // edges have 2 nodes
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 6; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        typedef boost::tuple<stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId> tri_tuple_type;
        static vector<tri_tuple_type> elems(6);

        CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        //stk::mesh::Part & active = mesh->ActivePart();
        //stk::mesh::Part & quad4  = mesh->QuadPart();

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

        //add_parts.push_back( &active );
        //FIXME
        //add_parts.push_back( const_cast<mesh::Part*>( eMesh.getPart(m_toTopoPartName) ));
        add_parts = m_toParts;

        /**
           \node[above] at (p4.side 1){2};
           \node[left] at (p4.side 2){3};
           \node[below] at (p4.side 3){0};
           \node[right] at (p4.side 4){1};
        */



        double tmp_x[3];
        for (int iedge = 0; iedge < 4; iedge++)
          {
            double * mp = midPoint(EDGE_COORD(iedge,0), EDGE_COORD(iedge,1), eMesh.get_spatial_dim(), tmp_x);

            if (!EDGE_N(iedge))
              {
                std::cout << "P[" << eMesh.get_rank() << " nid ## = 0 << " << std::endl;
              }
            eMesh.createOrGetNode(EDGE_N(iedge), mp);

          }


        elems[0] = tri_tuple_type(VERT_N(0), EDGE_N(0), EDGE_N(3));
        elems[1] = tri_tuple_type(VERT_N(1), EDGE_N(1), EDGE_N(0));
        elems[2] = tri_tuple_type(EDGE_N(0), EDGE_N(1), EDGE_N(3));

        elems[3] = tri_tuple_type(VERT_N(2), EDGE_N(2), EDGE_N(1));
        elems[4] = tri_tuple_type(VERT_N(3), EDGE_N(3), EDGE_N(2));
        elems[5] = tri_tuple_type(EDGE_N(2), EDGE_N(3), EDGE_N(1));

        // write a diagram of the refinement pattern as a vtk file, or a latex/tikz/pgf file
#define WRITE_DIAGRAM 0
#if WRITE_DIAGRAM

#endif

        bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            //stk::mesh::Entity newElement = eMesh.get_bulk_data()->declare_entity(Element, *element_id_pool, eMesh.getPart(interface_table::shards_Triangle_3) );
            //stk::mesh::Entity newElement = eMesh.get_bulk_data()->declare_entity(Element, *element_id_pool, eMesh.getPart(interface_table::shards_Triangle_3) );

            stk::mesh::Entity newElement = stk::mesh::Entity();
            if (!use_declare_element_side)
              newElement = *element_pool;

            stk::mesh::Entity nodes[3] = {
              eMesh.createOrGetNode(elems[ielem].get<0>()),
              eMesh.createOrGetNode(elems[ielem].get<1>()),
              eMesh.createOrGetNode(elems[ielem].get<2>())};

            create_side_element(eMesh, use_declare_element_side, nodes, 3, newElement);

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(m_eMesh.get_rank());
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            //eMesh.get_bulk_data()->change_entity_parts( newElement, add_parts, remove_parts );
            change_entity_parts(eMesh, element, newElement);

            {
              if (!elems[ielem].get<0>())
                {
                  std::cout << "P[" << eMesh.get_rank() << " nid = 0 << " << std::endl;
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
