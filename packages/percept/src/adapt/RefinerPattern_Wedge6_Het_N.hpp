// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Wedge6_Het_N_hpp
#define adapt_RefinerPattern_Wedge6_Het_N_hpp

#include <percept/mesh/geometry/volume/VolumeUtil.hpp>

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Quad4_Quad4_Transition.hpp"
#include "RefinerPattern_Tri3_Tri3_HangingNode.hpp"
#include "RefinerPattern_Quad4_Het_N.hpp"

#include "RefinerPattern_NoBreak.hpp"

#include <percept/PerceptBoostArray.hpp>
#include <adapt/DiscretizeWedge.hpp>

namespace percept {

  extern bool s_do_transition_break;
  extern bool s_allow_special_wedge_refine;  // for wedge boundary-layer refine

  // Some explanation: Wedge local (non-hanging-node) refinement pattern creates a heterogeneous mesh, so we create three
  //   sub-patterns to deal with the three different resulting topologies (wedge, tet and pyramid).  A third
  //   (parent) pattern is created to refer to the two sub-patterns, akin to URP_Heterogeneous_3D.

  // this pattern is only intended to be used for "transition" element breaking, that is, all hanging
  // nodes are removed by filling the interior of the wedge with sub-wedges, tets or pyramids

  //================================================================================================================================================================



  //================================================================================================================================================================
  // Heterogeneous child patterns
  //================================================================================================================================================================

  // Wedge to Tet
  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  struct WedgeTetPartial {};
  template <>
  class RefinerPattern<shards::Wedge<6>, shards::Tetrahedron<4>, -1, WedgeTetPartial > : public URP<shards::Wedge<6>, shards::Tetrahedron<4>  >
  {
    typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition > FaceBreaker;
    typedef RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1, TriHangingNode >  TriFaceBreakerType;
    FaceBreaker * m_face_breaker;
    TriFaceBreakerType * m_face_breaker_tri;

    TransitionElementType *m_transition_element_field;
    bool m_transition_element_field_set;

  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Wedge<6>, shards::Tetrahedron<4>  >(eMesh),
                                                                                                  m_face_breaker(0),m_face_breaker_tri(0),
                                                                                                  m_transition_element_field(0), m_transition_element_field_set(false)
    {
      m_primaryEntityRank = m_eMesh.element_rank();
      bool sameTopology = false;
      setNeededParts(eMesh, block_names, sameTopology);
      Elem::StdMeshObjTopologies::bootstrap();

      m_face_breaker =  new FaceBreaker(eMesh, block_names) ;
      m_face_breaker_tri =  new TriFaceBreakerType(eMesh, block_names) ;
    }

    ~RefinerPattern()
    {
      if (m_face_breaker) delete m_face_breaker;
      if (m_face_breaker_tri) delete m_face_breaker_tri;
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
    {
      EXCEPTWATCH;
      bp.resize(3);

      bp[0] = this;
      bp[1] = m_face_breaker;
      bp[2] = m_face_breaker_tri;
    }

    virtual void doBreak() {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
    {
      needed_entities.resize(3);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      needed_entities[2].first = m_eMesh.element_rank();
      setToOne(needed_entities);

      int faces[5] = {1,1,1,0,0};
      needed_entities[1].third.resize(5);
      needed_entities[1].third.assign(faces,faces+5);
    }

    // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
    virtual unsigned getNumNewElemPerElem() { return 48; }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                      vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
      static vector<wedge_to_tet_tuple_type> elems_tet(48);
      static vector<wedge_to_tet_tuple_type_local> elems_tet_local(48);
      static vector<wedge_to_pyr_tuple_type> elems_pyr(24);
      static vector<wedge_to_pyr_tuple_type_local> elems_pyr_local(24);
      static vector<wedge_to_wedge_tuple_type> elems_wedge(24);
      static vector<wedge_to_wedge_tuple_type_local> elems_wedge_local(24);
      unsigned num_new_elems=0;

      shards::CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
      //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

      if (!m_transition_element_field_set)
        {
          m_transition_element_field_set = true;
          m_transition_element_field = eMesh.get_transition_element_field();
          if(m_transition_element_field && m_transition_element_field->name() != "transition_element_3") m_transition_element_field = nullptr;
        }

      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Part*> remove_parts;
      add_parts = m_toParts;

      unsigned edge_marks[DiscretizeWedge::nedges] = { 0,0,0, 0,0,0, 0,0,0 };
      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < DiscretizeWedge::nedges; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (num_nodes_on_edge)
            {
              edge_marks[iedge] = 1;
              ++num_edges_marked;
            }
        }

      //unsigned num_faces_marked = 0;
      unsigned face_marks[DiscretizeWedge::nfaces] = {0,0,0,0,0};
      stk::mesh::EntityRank rank = m_eMesh.face_rank();

      for (int iface = 0; iface < DiscretizeWedge::nfaces; iface++)
        {
          if ( new_sub_entity_nodes[rank].size() )
            {
              if (new_sub_entity_nodes[rank][iface].size())
                {
                  face_marks[iface] = 1;
                  //++num_faces_marked;
                }
            }
        }

      //if (num_edges_marked == 0 && num_faces_marked == 0)
      if (num_edges_marked == 0)
        return;

      DiscretizeWedge::discretize(eMesh, element,
                                  edge_marks, face_marks,
                                  elems_tet_local, elems_pyr_local, elems_wedge_local, s_allow_special_wedge_refine);

#define Q_CENTROID_N_EV (useAltNode ? m_eMesh.identifier(new_nodes[0]) : NN(m_primaryEntityRank, 0) )

#define Q_CV_EV(i) \
      ( DiscretizeWedge::is_mapped_v(i) ? VERT_N(i) :    \
        ( DiscretizeWedge::is_mapped_e(i) ? EDGE_N(i - DiscretizeWedge::edge_offset) : \
          ( DiscretizeWedge::is_mapped_f(i) ? FACE_N(i - DiscretizeWedge::face_offset) : Q_CENTROID_N_EV ) ) )

      std::vector<stk::mesh::Entity> new_nodes;
      bool useAltNode = false;
      if (!new_sub_entity_nodes[m_primaryEntityRank][0].size())
        {
          //throw std::runtime_error("no centroid node available");
           // std::cout << " bad element: " << std::endl;
           // m_eMesh.print(element);
           // m_eMesh.dump_vtk("bad.vtk");
          if (1)
            {
              if (!m_eMesh.getEntitiesUsingIdServerNewNodes( 1, new_nodes))
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] node_rank entity pool deplenished in RefinerPattern_Wedge6_Het_N" << std::endl;
                  throw std::logic_error("node rank entity pool deplenished");
                }
              useAltNode = true;
            }
        }

      if (useAltNode)
        {
          nodeRegistry.forceInMap(new_nodes, NodeRegistry::NR_MARK, element, stk::topology::ELEMENT_RANK, 0u);
          nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);
          new_sub_entity_nodes[m_primaryEntityRank][0].resize(1);
          new_sub_entity_nodes[m_primaryEntityRank][0][0] = m_eMesh.identifier(new_nodes[0]);
        }

      num_new_elems = elems_tet_local.size();
      elems_tet.resize(num_new_elems);
      for (unsigned ielem=0; ielem < num_new_elems; ielem++)
        {
          wedge_to_tet_tuple_type ht;
          for (unsigned ii=0; ii < 4; ++ii)
            {
              unsigned ee = elems_tet_local[ielem][ii];
              bool isf = DiscretizeWedge::is_mapped_f(ee);
              if (isf)
                {
                  int ifaceOrd = ee - DiscretizeWedge::face_offset;
                  if (new_sub_entity_nodes[m_eMesh.face_rank()].size() &&
                      new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd].size() &&
                      new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd][0])
                    {
                      stk::mesh::EntityId id = new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd][0];
                      if (0) std::cout << "no error, id= " << id << std::endl;
                    }
                  else
                    {
                      std::cout << "error WedgeTetPartial" << std::endl;

                      std::set<stk::mesh::Entity> list;
                      list.insert(element);
                      m_eMesh.dump_vtk("err-face.vtk", false, &list);

                      SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
                      nodeRegistry.getSubDimEntity(subDimEntity, element, m_eMesh.face_rank(), ifaceOrd);
                      std::vector<stk::mesh::Entity> nv;
                      for (unsigned jj=0; jj < subDimEntity.size(); ++jj)
                        {
                          nv.push_back(subDimEntity[jj]);
                        }
                      m_eMesh.dump_vtk(nv, "err-nodes.vtk");
                    }
                }
              ht[ii] = Q_CV_EV(ee);
            }
          elems_tet[ielem] = ht;
        }

      for (unsigned ielem=0; ielem < elems_tet.size(); ielem++)
        {
          stk::mesh::Entity newElement = *element_pool;

          if (m_transition_element_field)
            {
              int *transition_element = 0;
              transition_element = stk::mesh::field_data( *m_transition_element_field , newElement );
              //if (num_edges_marked == 12 && num_faces_marked == 6)
              if (num_edges_marked == 9)
                {
                  transition_element[0] = 0;
                  // FIXME error???
                  //throw std::runtime_error("logic error in RefinerPattern_Wedge6_Het_N");
                }
              else
                {
                  transition_element[0] = 1;
                }
            }

          if (proc_rank_field)
            {
              double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
              if (fdata)
                fdata[0] = double(eMesh.owner_rank(newElement));
            }

          change_entity_parts(eMesh, element, newElement);

          {
            if (!elems_tet[ielem][0])
              {
                std::cout << "P[" << eMesh.get_rank() << "] nid = 0 << " << std::endl;
                //exit(1);
              }
          }

          // 4 nodes of the new tets
          for (unsigned ii=0; ii < 4; ++ii)
            {
              stk::mesh::Entity n0 = eMesh.createOrGetNode(elems_tet[ielem][ii]);
              eMesh.get_bulk_data()->declare_relation(newElement, n0, ii);
            }

          unsigned nchild = eMesh.numChildren(element);
          //set_parent_child_relations(eMesh, element, newElement, ielem);
          set_parent_child_relations(eMesh, element, newElement,  *ft_element_pool, nchild);

          if (0) std::cout << "tmp createNewElements: isParentElement= " << m_eMesh.isParentElement(element, false) << std::endl;
          std::vector<stk::mesh::Entity> elements(1,element);
          eMesh.prolongateElementFields( elements, newElement);

          ++ft_element_pool;
          ++element_pool;

        }
    }

  };

  // Wedge to Pyramid
  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  struct WedgePyrPartial {};
  template <>
  class RefinerPattern<shards::Wedge<6>, shards::Pyramid<5>, -1, WedgePyrPartial > : public URP< shards::Wedge<6>, shards::Pyramid<5>  >
  {
    typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition > QuadFaceBreaker;
    typedef RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > TriFaceBreaker;
    QuadFaceBreaker * m_quad_face_breaker;
    TriFaceBreaker * m_tri_face_breaker;

    TransitionElementType *m_transition_element_field;
    bool m_transition_element_field_set;

  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Wedge<6>, shards::Pyramid<5>  >(eMesh),
                                                                                                  m_quad_face_breaker(0), m_tri_face_breaker(0),
                                                                                                  m_transition_element_field(0), m_transition_element_field_set(false)
    {
      m_primaryEntityRank = m_eMesh.element_rank();
      bool sameTopology = false;
      setNeededParts(eMesh, block_names, sameTopology);
      Elem::StdMeshObjTopologies::bootstrap();

      m_quad_face_breaker =  new QuadFaceBreaker(eMesh, block_names) ;
      m_tri_face_breaker =  new TriFaceBreaker(eMesh, block_names) ;
    }

    ~RefinerPattern()
    {
      if (m_quad_face_breaker) delete m_quad_face_breaker;
      if (m_tri_face_breaker) delete m_tri_face_breaker;
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
    {
      EXCEPTWATCH;
      bp.resize(3);

      bp[0] = this;
      bp[1] = m_quad_face_breaker;
      bp[2] = m_tri_face_breaker;
    }

    virtual void doBreak() {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
    {
      needed_entities.resize(3);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      needed_entities[2].first = m_eMesh.element_rank();
      setToOne(needed_entities);

      int faces[5] = {1,1,1,0,0};
      needed_entities[1].third.resize(5);
      needed_entities[1].third.assign(faces,faces+5);
    }

    // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
    virtual unsigned getNumNewElemPerElem() { return 24; }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                      vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
      static vector<wedge_to_tet_tuple_type> elems_tet(48);
      static vector<wedge_to_tet_tuple_type_local> elems_tet_local(48);
      static vector<wedge_to_pyr_tuple_type> elems_pyr(24);
      static vector<wedge_to_pyr_tuple_type_local> elems_pyr_local(24);
      static vector<wedge_to_wedge_tuple_type> elems_wedge(24);
      static vector<wedge_to_wedge_tuple_type_local> elems_wedge_local(24);
      unsigned num_new_elems=0;

      shards::CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
      //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

      if (!m_transition_element_field_set)
        {
          m_transition_element_field_set = true;
          m_transition_element_field = eMesh.get_transition_element_field();
          if(m_transition_element_field->name() != "transition_element_3") m_transition_element_field = nullptr;
        }

      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Part*> remove_parts;
      add_parts = m_toParts;

      unsigned edge_marks[9] = { 0,0,0, 0,0,0, 0,0,0 };
      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < 9; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (num_nodes_on_edge)
            {
              edge_marks[iedge] = 1;
              ++num_edges_marked;
            }
        }

      unsigned face_marks[5] = {0,0,0,0,0};
      unsigned num_faces_marked = 0;
      stk::mesh::EntityRank rank = m_eMesh.face_rank();

      for (int iface = 0; iface < 5; iface++)
        {
          if ( new_sub_entity_nodes[rank].size() )
            {
              if (new_sub_entity_nodes[rank][iface].size())
                {
                  face_marks[iface] = 1;
                  ++num_faces_marked;
                }
            }
        }

      if (num_edges_marked == 0 && num_faces_marked == 0)
        return;

      DiscretizeWedge::discretize(eMesh, element,
                                  edge_marks, face_marks,
                                  elems_tet_local, elems_pyr_local, elems_wedge_local, s_allow_special_wedge_refine);

      std::vector<stk::mesh::Entity> new_nodes;
      bool useAltNode = false;
      if (!new_sub_entity_nodes[m_primaryEntityRank][0].size())
        {
          if (1)
            {
              if (!m_eMesh.getEntitiesUsingIdServerNewNodes( 1, new_nodes))
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] node_rank entity pool deplenished in RefinerPattern_Wedge6_Het_N" << std::endl;
                  throw std::logic_error("node rank entity pool deplenished");
                }
              useAltNode = true;
            }
        }

      if (useAltNode)
        {
          nodeRegistry.forceInMap(new_nodes, NodeRegistry::NR_MARK, element, stk::topology::ELEMENT_RANK, 0u);
          nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);
          new_sub_entity_nodes[m_primaryEntityRank][0].resize(1);
          new_sub_entity_nodes[m_primaryEntityRank][0][0] = m_eMesh.identifier(new_nodes[0]);
        }

      num_new_elems = elems_pyr_local.size();
      elems_pyr.resize(num_new_elems);
      for (unsigned ielem=0; ielem < num_new_elems; ielem++)
        {
          wedge_to_pyr_tuple_type hp;
          for (unsigned ii=0; ii < 5; ++ii)
            {
              unsigned ee = elems_pyr_local[ielem][ii];

              bool isf = DiscretizeWedge::is_mapped_f(ee);
              if (isf)
                {
                  int ifaceOrd = ee - DiscretizeWedge::face_offset;
                  if (new_sub_entity_nodes[m_eMesh.face_rank()].size() &&
                      new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd].size() &&
                      new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd][0])
                    {
                      stk::mesh::EntityId id = new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd][0];
                      if (0) std::cout << "no error, id= " << id << std::endl;
                    }
                  else
                    {
                      std::cout << "error WedgePyrPartial" << std::endl;

                      if (1)
                        {
                          std::vector<unsigned> em(&edge_marks[0], &edge_marks[0]+9);
                          std::vector<unsigned> fm(&face_marks[0], &face_marks[0]+5);
                          std::cout << "WedgePyrPartial:createNewElements start:: element= " << m_eMesh.identifier(element)
                                    << " num_edges_marked= " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                                    << "\n em= " << em
                                    << "\n fm= " << fm
                                    << "\n ielem= " << ielem
                                    << " num_new_elems= " << num_new_elems << " ee= " << ee << " ifaceOrd= " << ifaceOrd
                                    << " new_sub_entity_nodes[m_eMesh.face_rank()].size= " << new_sub_entity_nodes[m_eMesh.face_rank()].size()
                                    << std::endl;
                          eMesh.print_entity(element);

                          stk::mesh::Entity faceNeigh = eMesh.get_face_neighbor(element, ifaceOrd);
                          std::cout << "faceNeigh= " << eMesh.id(faceNeigh) << std::endl;
                          eMesh.print_entity(faceNeigh);
                        }

                      DiscretizeWedge::discretize(eMesh, element,
                                                  edge_marks, face_marks,
                                                  elems_tet_local, elems_pyr_local, elems_wedge_local, s_allow_special_wedge_refine, true);

                      std::set<stk::mesh::Entity> list;
                      list.insert(element);
                      m_eMesh.dump_vtk("err-face.vtk", false, &list);

                      SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
                      nodeRegistry.getSubDimEntity(subDimEntity, element, m_eMesh.face_rank(), ifaceOrd);
                      std::vector<stk::mesh::Entity> nv;
                      for (unsigned jj=0; jj < subDimEntity.size(); ++jj)
                        {
                          nv.push_back(subDimEntity[jj]);
                        }
                      m_eMesh.dump_vtk(nv, "err-nodes.vtk");
                    }
                }
              hp[ii] = Q_CV_EV(ee);
            }
          elems_pyr[ielem] = hp;
        }

      for (unsigned ielem=0; ielem < elems_pyr.size(); ielem++)
        {
          stk::mesh::Entity newElement = *element_pool;

          if (m_transition_element_field)
            {
              int *transition_element = 0;
              transition_element = stk::mesh::field_data( *m_transition_element_field , newElement );
              transition_element[0] = 1;
            }

          if (proc_rank_field)
            {
              double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
              if (fdata)
                fdata[0] = double(eMesh.owner_rank(newElement));
            }

          change_entity_parts(eMesh, element, newElement);

          {
            if (!elems_pyr[ielem][0])
              {
                std::cout << "P[" << eMesh.get_rank() << "] nid = 0 << " << std::endl;
                //exit(1);
              }
          }

          // 5 nodes of the new pyramids
          for (unsigned ii=0; ii < 5; ++ii)
            {
              stk::mesh::Entity n0 = eMesh.createOrGetNode(elems_pyr[ielem][ii]);
              eMesh.get_bulk_data()->declare_relation(newElement, n0, ii);
            }

          unsigned nchild = eMesh.numChildren(element);
          set_parent_child_relations(eMesh, element, newElement,  *ft_element_pool, nchild);

          std::vector<stk::mesh::Entity> elements(1,element);
          eMesh.prolongateElementFields( elements, newElement);

          ++ft_element_pool;
          ++element_pool;
        }
    }

  };

  // Wedge to Wedge
  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  struct WedgeWedgePartial {};
  template <>
  class RefinerPattern<shards::Wedge<6>, shards::Wedge<6>, -1, WedgeWedgePartial > : public URP<shards::Wedge<6>, shards::Wedge<6>  >
  {
    typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition > FaceBreaker;
    typedef RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1, TriHangingNode >  TriFaceBreakerType;
    FaceBreaker * m_face_breaker;
    TriFaceBreakerType * m_face_breaker_tri;

    TransitionElementType *m_transition_element_field;
    bool m_transition_element_field_set;

  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Wedge<6>, shards::Wedge<6>  >(eMesh),
                                                                                                  m_face_breaker(0), m_face_breaker_tri(0),
                                                                                                  m_transition_element_field(0), m_transition_element_field_set(false)
    {
      m_primaryEntityRank = m_eMesh.element_rank();
      bool sameTopology = true;
      setNeededParts(eMesh, block_names, sameTopology);
      Elem::StdMeshObjTopologies::bootstrap();

      m_face_breaker =  new FaceBreaker(eMesh, block_names) ;
      m_face_breaker_tri =  new TriFaceBreakerType(eMesh, block_names) ;
    }

    ~RefinerPattern()
    {
      if (m_face_breaker) delete m_face_breaker;
      if (m_face_breaker_tri) delete m_face_breaker_tri;
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
    {
      EXCEPTWATCH;
      bp.resize(3);

      bp[0] = this;
      bp[1] = m_face_breaker;
      bp[2] = m_face_breaker_tri;
    }

    virtual void doBreak() {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
    {
      needed_entities.resize(3);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      needed_entities[2].first = m_eMesh.element_rank();
      setToOne(needed_entities);

      int faces[5] = {1,1,1,0,0};
      needed_entities[1].third.resize(5);
      needed_entities[1].third.assign(faces,faces+5);
    }

    // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
    virtual unsigned getNumNewElemPerElem() { return 48; }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                      vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
      static vector<wedge_to_tet_tuple_type> elems_tet(48);
      static vector<wedge_to_tet_tuple_type_local> elems_tet_local(48);
      static vector<wedge_to_pyr_tuple_type> elems_pyr(24);
      static vector<wedge_to_pyr_tuple_type_local> elems_pyr_local(24);
      static vector<wedge_to_wedge_tuple_type> elems_wedge(24);
      static vector<wedge_to_wedge_tuple_type_local> elems_wedge_local(24);
      unsigned num_new_elems=0;

      //std::cout << "s_allow_special_wedge_refine= " << s_allow_special_wedge_refine << std::endl;
      if (!s_allow_special_wedge_refine)
        return;

      shards::CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
      //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

      if (!m_transition_element_field_set)
        {
          m_transition_element_field_set = true;
          m_transition_element_field = eMesh.get_transition_element_field();
          if(m_transition_element_field->name() != "transition_element_3") m_transition_element_field = nullptr;
        }

      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Part*> remove_parts;
      add_parts = m_toParts;

      unsigned edge_marks[DiscretizeWedge::nedges] = { 0,0,0, 0,0,0, 0,0,0 };
      unsigned num_edges_marked=0;
      bool special_iso = true;
      int special_te_case_2 = -1;
      int special_te = -1;
      for (int iedge = 0; iedge < DiscretizeWedge::nedges; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (iedge < 6 && !num_nodes_on_edge) special_iso=false;
          if (iedge >= 6 && num_nodes_on_edge) special_iso=false;
          if (num_nodes_on_edge)
            {
              edge_marks[iedge] = 1;
              ++num_edges_marked;
            }
        }
      if (special_iso)
        {
          VERIFY_OP_ON(num_edges_marked, ==, 6, "special_iso bad");
        }

      unsigned num_faces_marked = 0;
      unsigned face_marks[DiscretizeWedge::nfaces] = {0,0,0,0,0};
      stk::mesh::EntityRank rank = m_eMesh.face_rank();

      for (int iface = 0; iface < DiscretizeWedge::nfaces; iface++)
        {
          if ( new_sub_entity_nodes[rank].size() )
            {
              if (new_sub_entity_nodes[rank][iface].size())
                {
                  if (iface < 3) special_iso = false;
                  face_marks[iface] = 1;
                  ++num_faces_marked;
                }
            }
        }

      if (num_edges_marked == 2) // && num_faces_marked == 0)
        {
          special_te = -1;
          for (unsigned ii=0; ii < 3; ii++)
            {
              if (edge_marks[ii] && edge_marks[ii+3])
                {
                  special_te = ii;
                  break;
                }
            }
        }

      // special case: 2 edges of each tri face are marked
      if (num_edges_marked == 4) // && num_faces_marked == 0)
        {
          special_te_case_2 = -1;
          int nfound = 0;
          for (unsigned ii=0; ii < 3; ii++)
            {
              if (edge_marks[ii] && edge_marks[ii+3])
                {
                  ++nfound;
                  if (special_te_case_2 < 0)
                    special_te_case_2 = 0;
                  special_te_case_2 += 1 << (ii);
                }
            }
          if (nfound != 2)
            special_te_case_2 = -1;
        }

      if (num_edges_marked == 0)
        return;

      DiscretizeWedge::discretize(eMesh, element,
                                  edge_marks, face_marks,
                                  elems_tet_local, elems_pyr_local, elems_wedge_local, s_allow_special_wedge_refine);

      std::vector<stk::mesh::Entity> new_nodes;
      bool useAltNode = false;
      if (!new_sub_entity_nodes[m_primaryEntityRank][0].size())
        {
          if (1)
            {
              if (!m_eMesh.getEntitiesUsingIdServerNewNodes( 1, new_nodes))
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] node_rank entity pool deplenished in RefinerPattern_Wedge6_Het_N" << std::endl;
                  throw std::logic_error("node rank entity pool deplenished");
                }
              useAltNode = true;
            }
        }

      if (useAltNode)
        {
          nodeRegistry.forceInMap(new_nodes, NodeRegistry::NR_MARK, element, stk::topology::ELEMENT_RANK, 0u);
          nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);
          new_sub_entity_nodes[m_primaryEntityRank][0].resize(1);
          new_sub_entity_nodes[m_primaryEntityRank][0][0] = m_eMesh.identifier(new_nodes[0]);
        }

      num_new_elems = elems_wedge_local.size();
      if (num_new_elems && s_allow_special_wedge_refine)
        {
          if (!(special_iso || (special_te >= 0) || (special_te_case_2 >= 0)))
            {
              std::vector<unsigned> em(&edge_marks[0], &edge_marks[0]+DiscretizeWedge::nedges);
              std::vector<unsigned> fm(&face_marks[0], &face_marks[0]+DiscretizeWedge::nfaces);
              std::cout << "WedgeTetPartial:createNewElements:: element= " << m_eMesh.identifier(element)
                        << " num_edges_marked= " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                        << "\n em= " << em
                        << "\n fm= " << fm
                        << std::endl;
            }

          VERIFY_OP_ON((special_iso || (special_te >= 0) || (special_te_case_2 >= 0)), ==, true, "bad wedge marks");
        }
      elems_wedge.resize(num_new_elems);
      for (unsigned ielem=0; ielem < num_new_elems; ielem++)
        {
          wedge_to_wedge_tuple_type hw;
          for (unsigned ii=0; ii < 6; ++ii)
            {
              unsigned ee = elems_wedge_local[ielem][ii];
              bool isf = DiscretizeWedge::is_mapped_f(ee);
              VERIFY_OP_ON(isf, == , false, "can't happen");
              hw[ii] = Q_CV_EV(ee);
            }
          elems_wedge[ielem] = hw;
        }

      for (unsigned ielem=0; ielem < elems_wedge.size(); ielem++)
        {
          stk::mesh::Entity newElement = *element_pool;

          if (m_transition_element_field)
            {
              int *transition_element = 0;
              transition_element = stk::mesh::field_data( *m_transition_element_field , newElement );
              //if (num_edges_marked == 12 && num_faces_marked == 6)
              // if (s_allow_special_wedge_refine && (special_iso || special_te >= 0 || special_te_case_2 >= 0))
              //   {
              //     transition_element[0] = 1;
              //   }
              // else if (num_edges_marked == 9) // && num_faces_marked == 3)
              if (num_edges_marked == 9 || (s_allow_special_wedge_refine && special_iso))
                {
                  transition_element[0] = 0;
                  // FIXME error???
                  //throw std::runtime_error("logic error in RefinerPattern_Wedge6_Het_N");
                }
              else
                {
                  transition_element[0] = 1;
                }
            }

          if (proc_rank_field)
            {
              double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
              if (fdata)
                fdata[0] = double(eMesh.owner_rank(newElement));
            }

          change_entity_parts(eMesh, element, newElement);

          {
            if (!elems_wedge[ielem][0])
              {
                std::cout << "P[" << eMesh.get_rank() << "] nid = 0 << " << std::endl;
                //exit(1);
              }
          }

          // 6 nodes of the new wedges
          for (unsigned ii=0; ii < 6; ++ii)
            {
              stk::mesh::Entity n0 = eMesh.createOrGetNode(elems_wedge[ielem][ii]);
              eMesh.get_bulk_data()->declare_relation(newElement, n0, ii);
            }

          unsigned nchild = eMesh.numChildren(element);
          //set_parent_child_relations(eMesh, element, newElement, ielem);
          set_parent_child_relations(eMesh, element, newElement,  *ft_element_pool, nchild);

          if (0) std::cout << "tmp createNewElements: isParentElement= " << m_eMesh.isParentElement(element, false) << std::endl;
          std::vector<stk::mesh::Entity> elements(1,element);
          eMesh.prolongateElementFields( elements, newElement);

          ++ft_element_pool;
          ++element_pool;

        }
    }

  };

  //================================================================================================================================================================
  // Parent pattern
  //================================================================================================================================================================

  struct WedgeHet {};

  template<>
  class RefinerPattern<shards::Wedge<6>, shards::Wedge<6>, 20, WedgeHet > : public UniformRefinerPatternBase
  {
    std::vector<UniformRefinerPatternBase *> m_bp;
  public:
    std::vector<UniformRefinerPatternBase *> m_bp_exported;
    typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition > FaceBreaker;
    FaceBreaker * m_face_breaker;

  protected:

    percept::PerceptMesh& m_eMesh;

  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  m_eMesh(eMesh)
    {
      m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      Elem::StdMeshObjTopologies::bootstrap();

      m_bp.resize(6);

      m_bp[0] = new  RefinerPattern<shards::Wedge<6>, shards::Wedge<6>, -1, WedgeWedgePartial > (eMesh, block_names) ;
      m_bp[1] = new  RefinerPattern<shards::Wedge<6>,  shards::Tetrahedron<4>, -1, WedgeTetPartial >         (eMesh, block_names) ;
      m_bp[2] = new  RefinerPattern<shards::Wedge<6>,  shards::Pyramid<5>,     -1, WedgePyrPartial >         (eMesh, block_names) ;
      m_bp[3] = new  RefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, -1, TetTempPartialNoBreak > (eMesh, block_names) ;
      m_bp[4] = new  RefinerPattern<shards::Pyramid<5>,     shards::Pyramid<5>,     -1, PyrTempPartialNoBreak > (eMesh, block_names) ;
      m_bp[5] = new  RefinerPattern<shards::Wedge<6>,     shards::Wedge<6>,         -1, WedgeTempPartialNoBreak > (eMesh, block_names) ;

      m_bp[0]->setNeededParts(eMesh, block_names, true);
      m_bp[1]->setNeededParts(eMesh, block_names, false);
      m_bp[2]->setNeededParts(eMesh, block_names, false);
      m_bp[3]->setNeededParts(eMesh, block_names, true);
      m_bp[4]->setNeededParts(eMesh, block_names, true);
      m_bp[5]->setNeededParts(eMesh, block_names, true);

      // repeat setNeededParts to catch newly created parts
      bool skipConvertedParts = false;
      m_bp[0]->setNeededParts(eMesh, block_names, true, skipConvertedParts);
      m_bp[1]->setNeededParts(eMesh, block_names, false, skipConvertedParts);
      m_bp[2]->setNeededParts(eMesh, block_names, false, skipConvertedParts);
      m_bp[3]->setNeededParts(eMesh, block_names, true, skipConvertedParts);
      m_bp[4]->setNeededParts(eMesh, block_names, true, skipConvertedParts);
      m_bp[5]->setNeededParts(eMesh, block_names, true, skipConvertedParts);

      m_bp_exported.resize(0);
      m_bp_exported.push_back(m_bp[3]);
      m_bp_exported.push_back(m_bp[4]);
      m_bp_exported.push_back(m_bp[5]);

      m_face_breaker =  new FaceBreaker(eMesh, block_names) ;
    }

    ~RefinerPattern()
    {
      if (m_face_breaker) delete m_face_breaker;
      for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
        {
          if (m_bp[ibp]) delete m_bp[ibp];
        }
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
    {
      EXCEPTWATCH;
      bp.resize(0);
      bp = m_bp;
      //bp.push_back(m_face_breaker);
    }

    virtual void doBreak() {
      throw std::runtime_error("shouldn't call RefinerPattern_Wedge6_Het_N::doBreak()");
    }

    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
    {
      needed_entities.resize(3);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      needed_entities[2].first = m_eMesh.element_rank();
      setToOne(needed_entities);

      int faces[5] = {1,1,1,0,0};
      needed_entities[1].third.resize(5);
      needed_entities[1].third.assign(faces,faces+5);
    }

    virtual unsigned getNumNewElemPerElem() { return 20; }

    virtual unsigned getFromTypeKey()
    {
      return shards::Wedge<6>::key;
    }

    // this is a bit bogus, but need to return something
    virtual unsigned getToTypeKey()
    {
      return shards::Wedge<6>::key;
    }

    virtual std::string getFromTopoPartName() {
      shards::CellTopology cell_topo(getFromTopology());
      return cell_topo.getName();
    }
    virtual std::string getToTopoPartName() {
      shards::CellTopology cell_topo(getToTopology());
      return cell_topo.getName();
    }

    virtual const CellTopologyData * getFromTopology() { return shards::getCellTopologyData< shards::Wedge<6> >(); }
    virtual const CellTopologyData * getToTopology() { return shards::getCellTopologyData< shards::Wedge<6> >(); }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0)
    {
      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < 9; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (num_nodes_on_edge)
            {
              ++num_edges_marked;
            }
        }
      //unsigned num_faces_marked = 0;
      //stk::mesh::EntityRank rank = m_eMesh.face_rank();

      //for (int iface = 0; iface < 6; iface++)
      //  {
      //    if ( new_sub_entity_nodes[rank].size() )
      //      {
      //        if (new_sub_entity_nodes[rank][iface].size())
      //          ++num_faces_marked;
      //      }
      //  }

      if ( num_edges_marked == 9 )
        //if ( num_edges_marked == 9 && num_faces_marked == 5)
        {
          // vector<stk::mesh::Entity>::iterator start = element_pool;
          // m_bp[0]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, proc_rank_field);
          // ptrdiff_t n_gen = element_pool - start;
          // if (0)
          //   std::cout << "n_gen = " << n_gen << std::endl;
          throw std::runtime_error("not ready");
        }
      else
        {
          for (unsigned ii=0; ii < m_bp.size(); ++ii)
            {
              //m_bp[ii]->m_ep_begin = m_ep_begin;
              //m_bp[ii]->m_ep_end = m_ep_end;
              m_bp[ii]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
            }
        }

    }
  };

}

#undef Q_VERT_N
#undef Q_EDGE_N
#undef Q_CENTROID_N
#undef Q_CENTROID_N_EV
#undef Q_CV_EV

#endif
