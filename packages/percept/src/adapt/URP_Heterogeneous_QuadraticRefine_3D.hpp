// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_URP_Heterogeneous_QuadraticRefine_3D_hpp
#define adapt_URP_Heterogeneous_QuadraticRefine_3D_hpp


#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <Shards_CellTopologyData.h>

#define FACE_BREAKER_HETERO_QUADRATIC_REFINE_3D 0
#if FACE_BREAKER_HETERO_QUADRATIC_REFINE_3D
#include "UniformRefinerPattern_Quad4_Quad4_4_sierra.hpp"
#include "UniformRefinerPattern_Tri3_Tri3_4_sierra.hpp"
#endif

  namespace percept {

    class URP_Heterogeneous_QuadraticRefine_3D : public UniformRefinerPatternBase
    {

      std::vector<UniformRefinerPatternBase *> m_bp;

#if FACE_BREAKER_HETERO_QUADRATIC_REFINE_3D
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > * m_face_breaker;
#endif

    protected:

      percept::PerceptMesh& m_eMesh;

    public:

      URP_Heterogeneous_QuadraticRefine_3D(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  m_eMesh(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        Elem::StdMeshObjTopologies::bootstrap();

        // list all types of known break patterns to be used here
        m_bp.resize(0);

        // refine
        m_bp.push_back(  new UniformRefinerPattern<shards::Hexahedron<27>,    shards::Hexahedron<27>,    8, SierraPort > (eMesh, block_names) );
        m_bp.push_back(  new UniformRefinerPattern<shards::Wedge<15>,         shards::Wedge<15>,         8, SierraPort > (eMesh, block_names) );
        m_bp.push_back(  new UniformRefinerPattern<shards::Tetrahedron<10>,   shards::Tetrahedron<10>,   8, SierraPort > (eMesh, block_names) );

        // FIXME
        //m_bp.push_back(  new UniformRefinerPattern<shards::ShellQuadrilateral<9>,  shards::ShellQuadrilateral<9>,   4, SierraPort > (eMesh, block_names) );
        //m_bp.push_back(  new UniformRefinerPattern<shards::ShellTriangle<6>,       shards::ShellTriangle<6>,   4, SierraPort > (eMesh, block_names) );

        // enrich
//         m_bp.push_back ( new UniformRefinerPattern< shards::Wedge<6>, shards::Wedge<15>, 1, SierraPort >                 (eMesh, block_names) );
//         m_bp.push_back ( new UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<10>,  1, SierraPort >   (eMesh, block_names) );
//         m_bp.push_back ( new UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<27>,   1, SierraPort >   (eMesh, block_names) );

#if FACE_BREAKER_HETERO_QUADRATIC_REFINE_3D

        m_bp.push_back(  new UniformRefinerPattern<shards::Quadrilateral<9>, shards::Quadrilateral<9>, 4, SierraPort > (eMesh, block_names) );
        m_bp.push_back(  new UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > (eMesh, block_names) );

#endif

      }

      ~URP_Heterogeneous_QuadraticRefine_3D()
      {
        for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
          {
            if (m_bp[ibp]) delete m_bp[ibp];
          }
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;

        bp = m_bp;

      }

      virtual std::string getFromTopoPartName() override {
        shards::CellTopology cell_topo(getFromTopology());
        return cell_topo.getName();
      }
      virtual std::string getToTopoPartName() override {
        shards::CellTopology cell_topo(getToTopology());
        return cell_topo.getName();
      }

      virtual void doBreak() override
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_QuadraticRefine_3D::doBreak()");

      }
      virtual unsigned getFromTypeKey() override
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_QuadraticRefine_3D::getFromTypeKey()");
      }
      virtual unsigned getToTypeKey() override
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_QuadraticRefine_3D::getToTypeKey()");
      }

      virtual const CellTopologyData *  getFromTopology() override
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_QuadraticRefine_3D::getFromTopology()");
      }
      virtual const CellTopologyData *  getToTopology() override
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_3D::getToTopology()");
      }


      void fillNeededEntities(std::vector<NeededEntityType>& /*needed_entities*/) override
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_QuadraticRefine_3D::fillNeededEntities()");
      }

      virtual unsigned getNumNewElemPerElem() override
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_QuadraticRefine_3D::getNumNewElemPerElem()");
        return 8;
      }

      void
      createNewElements(percept::PerceptMesh& /*eMesh*/, NodeRegistry& /*nodeRegistry*/,
                        stk::mesh::Entity /*element*/,  NewSubEntityNodesType& /*new_sub_entity_nodes*/, vector<stk::mesh::Entity>::iterator& /*element_pool*/,
                                        vector<stk::mesh::Entity>::iterator& /*ft_element_pool*/,
                        stk::mesh::FieldBase */*proc_rank_field*/=0) override
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_QuadraticRefine_3D::createNewElements()");
      }
    };

  }

#endif
