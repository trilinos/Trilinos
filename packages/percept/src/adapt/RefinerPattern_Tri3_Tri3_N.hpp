// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Tri3_Tri3_N_sierra_hpp
#define adapt_RefinerPattern_Tri3_Tri3_N_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <adapt/TriangulateTri.hpp>

//#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"
#include "RefinerPattern_Line2_Line2_N.hpp"

  namespace percept {

    typedef TriangulateTri::tri_tuple_type_local tri_tuple_type_local;
    typedef TriangulateTri::tri_tuple_type tri_tuple_type;

    /// general refinement pattern

    // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
    template <>
    class RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > : public URP<shards::Triangle<3>,shards::Triangle<3>  >
    {

      RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;
      TransitionElementType *m_transition_element_field;
      bool m_transition_element_field_set;
    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Triangle<3>, shards::Triangle<3>  >(eMesh),
                                                                                                    m_edge_breaker(0),
                                                                                                    m_transition_element_field(0), m_transition_element_field_set(false)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        if (m_eMesh.get_spatial_dim() == 2)
          {
            m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;
          }

      }

      ~RefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp.resize(2);
        bp.assign(2, (UniformRefinerPatternBase *)0);

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

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[0].second = 1u;
      }

      // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
      virtual unsigned getNumNewElemPerElem() { return 4; }


      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        if (0 && eMesh.check_entity_duplicate(element))
          {
            throw std::logic_error("RefinerPattern_Tri3_Tri3_N::createNewElements bad duplicate element of PARENT!");
          }

        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        static vector<tri_tuple_type> elems(4);
        static vector<tri_tuple_type_local> elems_local(4);
        unsigned num_new_elems=0;

        shards::CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
        //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

        if (!m_transition_element_field_set)
          {
            m_transition_element_field_set = true;
            m_transition_element_field = eMesh.get_transition_element_field();
            if(m_transition_element_field && m_transition_element_field->name() != "transition_element") m_transition_element_field = nullptr;
          }


        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;
        add_parts = m_toParts;

        unsigned edge_marks[3] = {0,0,0};
        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 3; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                edge_marks[iedge] = 1;
                ++num_edges_marked;
              }
          }
        if (num_edges_marked == 0 && !allow_single_refine)
          return;

        stk::mesh::Entity elem_nodes_local[3] = {stk::mesh::Entity()};
        for (int inode=0; inode < 3; inode++)
          {
            elem_nodes_local[inode] = elem_nodes[inode].entity();
          }
        TriangulateTri tq;
        tq.triangulate_face(eMesh, elem_nodes_local, edge_marks, elems_local);

#define CV_EV(i) ( i < 3 ? VERT_N(i) : EDGE_N(i-3) )

        num_new_elems = elems_local.size();
        elems.resize(num_new_elems);
        for (unsigned ielem=0; ielem < num_new_elems; ielem++)
          {
            elems[ielem] = {CV_EV(elems_local[ielem][0] ), CV_EV(elems_local[ielem][1] ), CV_EV(elems_local[ielem][2] )};
          }

        //std::cout << "tmp RefinerPattern_Tri3_Tri3_N::num_edges_marked= " << num_edges_marked << std::endl;

        //nodeRegistry.prolongateCoords(*const_cast<stk::mesh::Entity>(&element), stk::topology::ELEMENT_RANK, 0u);

        bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity newElement = stk::mesh::Entity();
            if (!use_declare_element_side)
              newElement = *element_pool;

            stk::mesh::Entity nodes[3] = {eMesh.createOrGetNode(elems[ielem][0]),
                                          eMesh.createOrGetNode(elems[ielem][1]),
                                          eMesh.createOrGetNode(elems[ielem][2])};

            create_side_element(eMesh, use_declare_element_side, nodes, 3, newElement);

            if (m_transition_element_field)
              {
                int *transition_element = 0;
                transition_element = stk::mesh::field_data( *m_transition_element_field , newElement );
                VERIFY_OP_ON(transition_element, !=, 0, "transition_element");
                if (num_edges_marked == 3)
                  {
                    transition_element[0] = 0;
                  }
                else
                  {
                    transition_element[0] = 1;
                  }
                if (0)
                  std::cout << "tmp srk found transition_element = " << m_eMesh.identifier(newElement) << " num_edges_marked= " << num_edges_marked
                            << " val= " << transition_element[0]
                            << " element= " << m_eMesh.identifier(element)
                            << std::endl;
              }

            if (proc_rank_field && eMesh.get_spatial_dim() == 2)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(m_eMesh.get_rank());
                if (fdata)
                  fdata[0] = double(eMesh.owner_rank(newElement));
              }

            //eMesh.get_bulk_data()->change_entity_parts( newElement, add_parts, remove_parts );
            change_entity_parts(eMesh, element, newElement);


            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            std::vector<stk::mesh::Entity> elements(1,element);
            eMesh.prolongateElementFields( elements, newElement);

            // FIXME tmp - could be slow

            if (0 && eMesh.check_entity_duplicate(newElement))
              {
                /*
                if (eMesh.check_entity_duplicate(element))
                  {
                    std::cout << "RefinerPattern_Tri3_Tri3_N::createNewElements bad duplicate element of PARENT 2!" << std::endl;
                  }
                */
                std::cout << "RefinerPattern_Tri3_Tri3_N bad duplicate element= " << element << " newElement= " << newElement << " elems.size() = " << elems.size() << std::endl;
                std::cout << "===> newElement.ischild, is parent = " << eMesh.isChildElement(newElement) << " " << eMesh.isParentElement(newElement) << std::endl;
                throw std::logic_error("RefinerPattern_Tri3_Tri3_N::createNewElements bad duplicate element");
              }

            ft_element_pool++;
            if (!use_declare_element_side)
              element_pool++;

          }


      }

    };

  }

#endif
