// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_URP_Heterogeneous_Enrich_3D_hpp
#define adapt_URP_Heterogeneous_Enrich_3D_hpp


#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <Shards_CellTopologyData.h>


#include "UniformRefinerPattern_Quad4_Quad9_1_sierra.hpp"
#include "UniformRefinerPattern_Tri3_Tri6_1_sierra.hpp"
#include "UniformRefinerPattern_ShellLine2_ShellLine3_1_sierra.hpp"
#include "UniformRefinerPattern_Beam2_Beam3_1_sierra.hpp"
#include "UniformRefinerPattern_ShellQuad4_ShellQuad9_1_sierra.hpp"
#include "UniformRefinerPattern_ShellTri3_ShellTri6_1_sierra.hpp"

  namespace percept {

    class URP_Heterogeneous_Enrich_3D : public UniformRefinerPatternBase
    {

      std::vector<UniformRefinerPatternBase *> m_bp;


    protected:

      percept::PerceptMesh& m_eMesh;

    public:

      URP_Heterogeneous_Enrich_3D(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  m_eMesh(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        Elem::StdMeshObjTopologies::bootstrap();

        // list all types of known break patterns to be used here
        m_bp.resize(0);

        int spatialDim = eMesh.get_spatial_dim();

        // refine

        // put them in reverse topological rank order

        if (spatialDim != 3)
          {
            throw std::runtime_error("URP_Heterogeneous_Enrich_3D is only for 3D meshes");
          }


        // enrich
        m_bp.push_back ( new UniformRefinerPattern<shards::Wedge<6>,        shards::Wedge<15>, 1, SierraPort >          (eMesh, block_names) );
        m_bp.push_back ( new UniformRefinerPattern<shards::Tetrahedron<4>,  shards::Tetrahedron<10>,  1, SierraPort >   (eMesh, block_names) );
        m_bp.push_back ( new UniformRefinerPattern<shards::Hexahedron<8>,   shards::Hexahedron<27>,   1, SierraPort >   (eMesh, block_names) );
        m_bp.push_back ( new UniformRefinerPattern<shards::Pyramid<5>,      shards::Pyramid<13>, 1, SierraPort >        (eMesh, block_names) );

        m_bp.push_back(  new UniformRefinerPattern<shards::ShellQuadrilateral<4>,  shards::ShellQuadrilateral<9>,   1, SierraPort > (eMesh, block_names) );
        m_bp.push_back(  new UniformRefinerPattern<shards::ShellTriangle<3>,       shards::ShellTriangle<6>,   1, SierraPort >      (eMesh, block_names) );
        m_bp.push_back(  new UniformRefinerPattern<shards::Beam<2>,       shards::Beam<3>,   1, SierraPort >      (eMesh, block_names) );

        m_bp.push_back(  new UniformRefinerPattern<shards::Quadrilateral<4>,       shards::Quadrilateral<9>, 1, SierraPort >        (eMesh, block_names) );
        m_bp.push_back(  new UniformRefinerPattern<shards::Triangle<3>,            shards::Triangle<6>, 1, SierraPort >             (eMesh, block_names) );

        m_bp.push_back(  new UniformRefinerPattern<shards::ShellLine<2>,       shards::ShellLine<3>,   1, SierraPort >      (eMesh, block_names) );

      }

      ~URP_Heterogeneous_Enrich_3D()
      {
        for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
          {
            if (m_bp[ibp]) delete m_bp[ibp];
          }
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;

        bp = m_bp;

      }

      virtual std::string getFromTopoPartName() {
        shards::CellTopology cell_topo(getFromTopology());
        return cell_topo.getName();
      }
      virtual std::string getToTopoPartName() {
        shards::CellTopology cell_topo(getToTopology());
        return cell_topo.getName();
      }

      virtual void doBreak()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::doBreak()");

      }
      virtual unsigned getFromTypeKey()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::getFromTypeKey()");
      }
      virtual unsigned getToTypeKey()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::getFromTypeKey()");
      }

      virtual const CellTopologyData *  getFromTopology()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::getFromTopology()");
      }
      virtual const CellTopologyData *  getToTopology()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_3D::getToTopology()");
      }


      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::fillNeededEntities()");
      }

      virtual unsigned getNumNewElemPerElem()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::getNumNewElemPerElem()");
        return 8;
      }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::createNewElements()");
      }
    };

  }

#endif
