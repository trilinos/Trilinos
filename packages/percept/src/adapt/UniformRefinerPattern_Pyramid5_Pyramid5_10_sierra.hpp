// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Pyramid5_Pyramid5_10_sierra_hpp
#define adapt_UniformRefinerPattern_Pyramid5_Pyramid5_10_sierra_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Quad4_Quad4_4_sierra.hpp"
#include "UniformRefinerPattern_Tri3_Tri3_4_sierra.hpp"

#include <percept/PerceptBoostArray.hpp>

  namespace percept {

#define DEBUG_Pyramid5_Pyramid5_10 0

    // Some explanation: Pyramid refinement pattern creates a heterogeneous mesh, so we create two
    //   sub-patterns to deal with the two different resulting topologies (pyramid and tet).  A third
    //   (parent) pattern is created to refer to the two sub-patterns, akin to URP_Heterogeneous_3D.

    //================================================================================================================================================================
    template <>
    class UniformRefinerPattern<shards::Pyramid<5>, shards::Pyramid<5>, 6, SierraPort > : public URP<shards::Pyramid<5>, shards::Pyramid<5>  >
    {


    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType /*block_names*/ = BlockNamesType()) :  URP<shards::Pyramid<5>, shards::Pyramid<5>  >(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;
        Elem::StdMeshObjTopologies::bootstrap();
        m_do_strip_hashes=false;
      }

      ~UniformRefinerPattern()
      {
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;
        bp.resize(1);
        bp[0] = this;
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(3);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = m_eMesh.face_rank();
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 6; }


      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool, 
                                        proc_rank_field);
      }

    };

    //================================================================================================================================================================
    template <>
    class UniformRefinerPattern<shards::Pyramid<5>, shards::Tetrahedron<4>, 4, SierraPort > : public URP<shards::Pyramid<5>, shards::Tetrahedron<4>  >
    {

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType /*block_names*/ = BlockNamesType()) :  URP<shards::Pyramid<5>, shards::Tetrahedron<4>  >(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        Elem::StdMeshObjTopologies::bootstrap();
        m_do_strip_hashes=false;

      }

      ~UniformRefinerPattern()
      {
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;
        bp.resize(1);

        bp[0] = this;
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(3);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = m_eMesh.face_rank();
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 4; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool, 
                                        proc_rank_field);
      }

    };

    //================================================================================================================================================================
    template <>
    //class UniformRefinerPattern<shards::Pyramid<5>, shards::Pyramid<5>, 6, SierraPort > : public URP<shards::Pyramid<5>, shards::Pyramid<5>  >
    class UniformRefinerPattern<shards::Pyramid<5>, shards::Pyramid<5>, 10, SierraPort > : public UniformRefinerPatternBase
    {
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > * m_face_breaker;
      UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > * m_face_breaker_tri;

    protected:

      percept::PerceptMesh& m_eMesh;

    public:
      std::vector<UniformRefinerPatternBase *> m_bp;  // FIXME - protect with suitable design

      //UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Pyramid<5>, shards::Pyramid<5>  >(eMesh), m_eMesh(eMesh)
      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : m_eMesh(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;
        m_do_strip_hashes=false;

        Elem::StdMeshObjTopologies::bootstrap();

        // list all types of known break patterns to be used here
        m_bp.resize(0);

        m_bp.push_back(new  UniformRefinerPattern<shards::Pyramid<5>,       shards::Pyramid<5>,      6, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Pyramid<5>,       shards::Tetrahedron<4>,          4, SierraPort > (eMesh, block_names));

        bool sameTopology = false;
        //setNeededParts(eMesh, block_names, sameTopology);
        m_bp[0]->setNeededParts(eMesh, block_names, true); // don't force a new part for pyramids
        if (DEBUG_Pyramid5_Pyramid5_10)
          {
            std::cout << "tmp Pyramid5_Pyramid5_10 printParts m_bp[0]= " ; printParts(m_bp[0]);
          }
        m_bp[1]->setNeededParts(eMesh, block_names, sameTopology);
        if (DEBUG_Pyramid5_Pyramid5_10)
          {
            std::cout << "tmp Pyramid5_Pyramid5_10 printParts m_bp[1]= " ; printParts(m_bp[1]);
          }

        // repeat to catch new parts
        bool skipConvertedParts = false;
        m_bp[0]->setNeededParts(eMesh, block_names, true, skipConvertedParts); // don't force a new part for pyramids
        m_bp[1]->setNeededParts(eMesh, block_names, sameTopology, skipConvertedParts);

        for (int ibp=0; ibp < 2; ibp++)
          {
            bool merge=true;
            mergeOrAddParts(m_bp[ibp], this, merge);
          }

        if (DEBUG_Pyramid5_Pyramid5_10)
          {
            std::cout << "tmp Pyramid5_Pyramid5_10 printParts this= " ;
            printParts(this);
          }

        m_face_breaker =  new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > (eMesh, block_names) ;
        m_face_breaker_tri = new UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > (eMesh, block_names);
      }

      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
        if (m_face_breaker_tri) delete m_face_breaker_tri;
        for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
          {
            if (m_bp[ibp]) delete m_bp[ibp];
          }
      }

      void setSubPatternsForSetNeededParts( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        //m_transition_breaker->setSubPatterns(bp, eMesh);
        bp = m_bp;
        bp.push_back(this);
        bp.push_back( m_face_breaker);
        bp.push_back( m_face_breaker_tri );
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;

        bp.resize(3u);
        bp[0] = this;
        bp[1] = m_face_breaker;
        bp[2] = m_face_breaker_tri;
      }

      virtual void doBreak() override
      {
        throw std::runtime_error("shouldn't call URP_Pyramid5_Pyramid5::doBreak()");
      }
      virtual unsigned getFromTypeKey() override
      {
        //throw std::runtime_error("shouldn't call URP_Pyramid5_Pyramid5::getFromTypeKey()");
        return shards::Pyramid<5>::key;
      }
      virtual unsigned getToTypeKey() override
      {
        //throw std::runtime_error("shouldn't call URP_Pyramid5_Pyramid5::getToTypeKey()");
        return shards::Pyramid<5>::key;
      }

      virtual std::string getFromTopoPartName() override {
        shards::CellTopology cell_topo(getFromTopology());
        return cell_topo.getName();
      }
      virtual std::string getToTopoPartName() override {
        shards::CellTopology cell_topo(getToTopology());
        return cell_topo.getName();
      }

      virtual const CellTopologyData * getFromTopology() override { return shards::getCellTopologyData< shards::Pyramid<5> >(); }
      virtual const CellTopologyData * getToTopology() override { return shards::getCellTopologyData< shards::Pyramid<5> >(); }

//       virtual const CellTopologyData *  getFromTopology()
//       {
//         throw std::runtime_error("shouldn't call URP_Pyramid5_Pyramid5::getFromTopology()");
//       }

//       virtual const CellTopologyData *  getToTopology() {
//         throw std::runtime_error("shouldn't call URP_Pyramid5_Pyramid5::getToTopology()");
//       }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        //throw std::runtime_error("shouldn't call URP_Pyramid5_Pyramid5::fillNeededEntities()");
        needed_entities.resize(3);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = m_eMesh.face_rank();
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override
      {
        //throw std::runtime_error("shouldn't call URP_Pyramid5_Pyramid5::getNumNewElemPerElem()");
        return 10;
      }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        //throw std::runtime_error("shouldn't call URP_Pyramid5_Pyramid5::createNewElements()");

        // pyramids
        m_bp[0]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
        // tets
        m_bp[1]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
      }

    };

  }

#endif
