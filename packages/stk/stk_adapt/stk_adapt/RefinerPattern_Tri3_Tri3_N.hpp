#ifndef stk_adapt_RefinerPattern_Tri3_Tri3_N_sierra_hpp
#define stk_adapt_RefinerPattern_Tri3_Tri3_N_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"

namespace stk {
  namespace adapt {

    typedef boost::tuple<unsigned, unsigned, unsigned> tri_tuple_type_local;
    typedef boost::tuple<stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId> tri_tuple_type;

    /// general refinement pattern
    
    // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
    template <>
    class RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > : public URP<shards::Triangle<3>,shards::Triangle<3>  >
    {

      UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > * m_edge_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Triangle<3>, shards::Triangle<3>  >(eMesh),
                                                                                                    m_edge_breaker(0)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.getSpatialDim() == 2)
          m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        if (m_eMesh.getSpatialDim() == 2)
          {
            m_edge_breaker =  new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) ;
          }

      }

      ~RefinerPattern() 
      {
        if (m_edge_breaker) delete m_edge_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);

        if (eMesh.getSpatialDim() == 2)
          {
            bp[0] = this;
            if (m_eMesh.getSpatialDim() == 2)
              {
                bp[1] = m_edge_breaker;
              }
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

      /**
       *
       *   Convention: input is the element's nodes and the marks on the edges.  Output is an array
       *     of "elements" defined as local id's of nodes forming those elements, where {0,1,2} represent
       *     the original vertices and {3,4,5} are the edges:
       *
       *                2  
       *                o
       *               / \
       *              /   \
       *             /     \
       *          5 *       * 4
       *           /         \
       *          /           \
       *         /             \ 
       *        o-------*-------o
       *       0        3        1
       */

      // Note: this will form the basis of triangulating faces in 3D, so it is generalized to a 
      //   generic method.
      // Note: code below doesn't orient the face except for a rotation - we need a polarity flip check as
      //   well for the general, 3D face case
      //

#define T_VERT_N(i) (i)
#define T_EDGE_N(i) ((i)+3)

      static void triangulate_face(PerceptMesh& eMesh, stk::mesh::Entity *elem_nodes[3], unsigned edge_marks[3], 
                                   vector<tri_tuple_type_local>& elems)
      {
        elems.resize(0);

        const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< shards::Triangle<3> >();

        shards::CellTopology cell_topo(cell_topo_data);
        //const stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK); /NLM
        VectorFieldType* coordField = eMesh.getCoordinatesField();

        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 3; iedge++)
          {
            unsigned num_nodes_on_edge = edge_marks[iedge];
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }

        //std::cout << "tmp RefinerPattern_Tri3_Tri3_N::num_edges_marked= " << num_edges_marked << std::endl;

        if (num_edges_marked == 3)
          {
            elems.resize(4);

            elems[0] = tri_tuple_type( T_VERT_N(0),    T_EDGE_N(0), T_EDGE_N(2) );
            elems[1] = tri_tuple_type( T_VERT_N(1),    T_EDGE_N(1), T_EDGE_N(0) );
            elems[2] = tri_tuple_type( T_VERT_N(2),    T_EDGE_N(2), T_EDGE_N(1) );
            elems[3] = tri_tuple_type( T_EDGE_N(0),    T_EDGE_N(1), T_EDGE_N(2) );
          }
        else if (num_edges_marked == 2)
          {
            /**
             *
             *   case 1: jedge == max length edge
             *
             *                i2  
             *                o
             *               /|\
             *              / | \
             *             /  |  \
             *            /   |   * jedgep
             *           /    |  / \
             *          /     | /   \
             *         /      |/     \ 
             *        o-------*-------o
             *       i0      jedge     i1
             *
             *
             *   case 2: jedge+1 == max length edge
             *
             *                i2 
             *                o
             *               / \ 
             *              /   \
             *             /     \
             *            /     _.* jedgep
             *           /   _.* / \
             *          / _.*   /   \
             *         /.*     /     \
             *        o-------*-------o
             *       i0      jedge     i1
             *
             */

            elems.resize(3);

            // find first of two marked edges in sequence (get in "standard" orientation), and longest marked edge
            int jedge = -1;
            int jedge_max_edge = -1;
            double max_edge_length = -1.0;
            unsigned id_diff_0 = 0u;
            unsigned id_diff_1 = 0u;
            for (int iedge = 0; iedge < 3; iedge++)
              {

                unsigned num_nodes_on_edge = edge_marks[iedge];
                unsigned num_nodes_on_edge_p = edge_marks[(iedge+1)%3];
                if (num_nodes_on_edge && num_nodes_on_edge_p)
                  {
                    jedge = iedge;
                  }

                if (num_nodes_on_edge)
                  {
                    stk::mesh::Entity * node_0 = elem_nodes[cell_topo_data->edge[iedge].node[0]];
                    stk::mesh::Entity * node_1 = elem_nodes[cell_topo_data->edge[iedge].node[1]];

                    bool reverse = false;
                    // ensure edge_len is computed identically, independent of edge orientation
                    if (node_0->identifier() > node_1->identifier())
                      {
                        reverse = true;
                        stk::mesh::Entity *node_temp = node_0;
                        node_0 = node_1;
                        node_1 = node_temp;
                      }

                    double * const coord_0 = stk::mesh::field_data( *coordField , *node_0 );
                    double * const coord_1 = stk::mesh::field_data( *coordField , *node_1 );
                    double edge_len_squared = 0.0;

                    edge_len_squared = 
                      (coord_0[0] - coord_1[0])*(coord_0[0] - coord_1[0])+
                      (coord_0[1] - coord_1[1])*(coord_0[1] - coord_1[1])+
                      (eMesh.getSpatialDim() == 2 ? 0 : 
                       (coord_0[2] - coord_1[2])*(coord_0[2] - coord_1[2]) );

                    if (edge_len_squared > max_edge_length)
                      {
                        id_diff_0 = node_0->identifier();
                        id_diff_1 = node_1->identifier();
                        max_edge_length = edge_len_squared;
                        jedge_max_edge = iedge;
                      }
                    // intentional floating-point comparison (tie-break)
                    else if (edge_len_squared == max_edge_length)
                      {
                        unsigned loc_id_diff_0 = node_0->identifier();
                        unsigned loc_id_diff_1 = node_1->identifier();
                        bool lexical_less = false;
                        if (loc_id_diff_0 < id_diff_0)
                          {
                            lexical_less = true;
                          }
                        else if (loc_id_diff_0 == id_diff_0 && loc_id_diff_1 < id_diff_1)
                          {
                            lexical_less = true;
                          }
                        if (!lexical_less)
                          {
                            max_edge_length = edge_len_squared;
                            jedge_max_edge = iedge;
                          }
                      }
                  }
              }

            if (jedge < 0 || jedge_max_edge < 0) 
              {
                std::cout << "jedge = " << jedge << " jedge_max_edge = " << jedge_max_edge << std::endl;
                throw std::runtime_error("RefinerPattern_Tri3_Tri3_N jedge < 0");
              }

            //stk::mesh::Entity & node0 = *elem_nodes[iii].entity();
            int i0 = cell_topo_data->edge[jedge].node[0];
            if (i0 != jedge)
              {
                std::cout << "i0 = " << i0 << " jedge= " << jedge << std::endl;
                throw std::runtime_error("RefinerPattern_Tri3_Tri3_N i0 != jedge");
              }

            int i1 = (i0+1)%3;
            int i2 = (i0+2)%3;
            int jedgep = (jedge+1)%3;
            if (jedge_max_edge == jedge)
              {
                elems[0] = tri_tuple_type( T_VERT_N(i0),    T_EDGE_N(jedge),         T_VERT_N(i2)       );
                elems[1] = tri_tuple_type( T_EDGE_N(jedge), T_EDGE_N( jedgep ),      T_VERT_N(i2)       );
                elems[2] = tri_tuple_type( T_EDGE_N(jedge), T_VERT_N(i1),            T_EDGE_N( jedgep ) );
              }
            else
              {
                elems[0] = tri_tuple_type( T_VERT_N(i0),    T_EDGE_N(jedge),         T_EDGE_N( jedgep ) );
                elems[1] = tri_tuple_type( T_VERT_N(i0),    T_EDGE_N( jedgep ),      T_VERT_N(i2)       );
                elems[2] = tri_tuple_type( T_EDGE_N(jedge), T_VERT_N(i1),            T_EDGE_N( jedgep ) );
              }
          }
        else if (num_edges_marked == 1)
          {
            elems.resize(2);
            for (int iedge = 0; iedge < 3; iedge++)
              {
                unsigned num_nodes_on_edge = edge_marks[iedge];
                if (num_nodes_on_edge)
                  {
                    elems[0] = tri_tuple_type(T_VERT_N(iedge), T_EDGE_N(iedge), T_VERT_N((iedge+2)%3) );
                    elems[1] = tri_tuple_type(T_EDGE_N(iedge), T_VERT_N((iedge+1)%3), T_VERT_N((iedge+2)%3) );
                    break;
                  }
              }
          }
        else if (num_edges_marked == 0)
          {
#if 0
            // this allows each level to be at the same hierarchical level by having a single parent to single child
            elems.resize(1);
            elems[0] = tri_tuple_type(T_VERT_N(0), T_VERT_N(1), T_VERT_N(2) );
#else
            if (elems.size() != 0)
              {
                std::cout << "tmp num_edges_marked= 0 " << elems.size() << std::endl;
                throw std::logic_error("hmmmmmmmmmmm");
              }

            return;
#endif
          }

      }



      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity *>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        if (0 && eMesh.check_entity_duplicate(element))
          {
            throw std::logic_error("RefinerPattern_Tri3_Tri3_N::createNewElements bad duplicate element of PARENT!");
          }

        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::tuple<stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId> tri_tuple_type;
        typedef boost::tuple<int, int, int> tri_tuple_type_int;
        static vector<tri_tuple_type> elems(4);
        static vector<tri_tuple_type_local> elems_local(4);
        unsigned num_new_elems=0;

        shards::CellTopology cell_topo(cell_topo_data);
        const stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
        //VectorFieldType* coordField = eMesh.getCoordinatesField();

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
        if (num_edges_marked == 0)
          return;

        stk::mesh::Entity *elem_nodes_local[3] = {0,0,0};
        for (int inode=0; inode < 3; inode++)
          {
            elem_nodes_local[inode] = elem_nodes[inode].entity();
          }
        triangulate_face(eMesh, elem_nodes_local, edge_marks, elems_local);
        
#define CV_EV(i) ( i < 3 ? VERT_N(i) : EDGE_N(i-3) )

        num_new_elems = elems_local.size();
        elems.resize(num_new_elems);
        for (unsigned ielem=0; ielem < num_new_elems; ielem++)
          {
            elems[ielem] = tri_tuple_type( CV_EV(elems_local[ielem].get<0>() ), CV_EV(elems_local[ielem].get<1>() ), CV_EV(elems_local[ielem].get<2>() ) );
          }

        //std::cout << "tmp RefinerPattern_Tri3_Tri3_N::num_edges_marked= " << num_edges_marked << std::endl;

        //nodeRegistry.makeCentroidCoords(*const_cast<stk::mesh::Entity *>(&element), m_eMesh.element_rank(), 0u);
        
        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity& newElement = *(*element_pool);

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(m_eMesh.getRank());
                if (fdata)
                  fdata[0] = double(newElement.owner_rank());
              }

            //eMesh.getBulkData()->change_entity_parts( newElement, add_parts, remove_parts );
            change_entity_parts(eMesh, element, newElement);

            set_parent_child_relations(eMesh, element, newElement, ielem);

            interpolateElementFields(eMesh, element, newElement);

            {
              if (!elems[ielem].get<0>())
                {
                  std::cout << "P[" << eMesh.getRank() << "] nid = 0 << " << std::endl;
                  //exit(1);
                }
            }

            // 3 nodes of the new tris
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<0>()), 0);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<1>()), 1);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<2>()), 2);

            // FIXME tmp - could be slow

            if (0 && eMesh.check_entity_duplicate(newElement))
              {
                if (eMesh.check_entity_duplicate(element))
                  {
                    std::cout << "RefinerPattern_Tri3_Tri3_N::createNewElements bad duplicate element of PARENT 2!" << std::endl;
                  }
                std::cout << "RefinerPattern_Tri3_Tri3_N bad duplicate element= " << element << " newElement= " << newElement << " elems.size() = " << elems.size() << std::endl;
                std::cout << "===> newElement.ischild, is parent = " << eMesh.isChildElement(newElement) << " " << eMesh.isParentElement(newElement) << std::endl;
                throw std::logic_error("RefinerPattern_Tri3_Tri3_N::createNewElements bad duplicate element");
              }
            
            element_pool++;

          }

      
      }
      
    };

  }
}
#endif
