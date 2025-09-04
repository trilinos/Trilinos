// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_URP_Tet4_Wedge6_Hex8_hpp
#define adapt_URP_Tet4_Wedge6_Hex8_hpp


#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <Shards_CellTopologyData.h>

  namespace percept {

    class URP_Tet4_Wedge6_Hex8 : public UniformRefinerPatternBase
    {

      std::vector<UniformRefinerPatternBase *> m_bp;

    protected:

      percept::PerceptMesh& m_eMesh;

    public:

      URP_Tet4_Wedge6_Hex8(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  m_eMesh(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        //!setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        // list all types of known break patterns to be used here
        m_bp.resize(0);

        int spatialDim = eMesh.get_spatial_dim();

        // refine

        // put them in reverse topological rank order

        if (spatialDim == 3)
          {
            m_bp.push_back(new  UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Hexahedron<8>,   4 > (eMesh, block_names));
            m_bp.push_back(new  UniformRefinerPattern<shards::Wedge<6>,         shards::Hexahedron<8>,   6 > (eMesh, block_names));
          }

        m_bp.push_back(new  UniformRefinerPattern<shards::Triangle<3>,      shards::Quadrilateral<4>, 3, Specialization > (eMesh, block_names) );
        m_bp.push_back(new  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > (eMesh, block_names));

        m_bp.push_back(new  UniformRefinerPattern<shards::Line<2>,          shards::Line<2>,          2, SierraPort > (eMesh, block_names));

      }

      ~URP_Tet4_Wedge6_Hex8()
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

      virtual void doBreak() override
      {
        throw std::runtime_error("shouldn't call URP_Tet4_Wedge6_Hex8::doBreak()");

      }
      virtual unsigned getFromTypeKey() override
      {
        throw std::runtime_error("shouldn't call URP_Tet4_Wedge6_Hex8::getFromTypeKey()");

      }
      virtual unsigned getToTypeKey() override
      {
        throw std::runtime_error("shouldn't call URP_Tet4_Wedge6_Hex8::getToTypeKey()");

      }

      virtual std::string getFromTopoPartName() override {
        shards::CellTopology cell_topo(getFromTopology());
        return cell_topo.getName();
      }
      virtual std::string getToTopoPartName() override {
        shards::CellTopology cell_topo(getToTopology());
        return cell_topo.getName();
      }

      virtual const CellTopologyData *  getFromTopology() override
      {
        throw std::runtime_error("shouldn't call URP_Tet4_Wedge6_Hex8::getFromTopology()");
      }

      virtual const CellTopologyData *  getToTopology() override {
        throw std::runtime_error("shouldn't call URP_Tet4_Wedge6_Hex8::getToTopology()");
      }

      void fillNeededEntities(std::vector<NeededEntityType>& /*needed_entities*/) override
      {
        throw std::runtime_error("shouldn't call URP_Tet4_Wedge6_Hex8::fillNeededEntities()");
      }

      virtual unsigned getNumNewElemPerElem() override
      {
        throw std::runtime_error("shouldn't call URP_Tet4_Wedge6_Hex8::getNumNewElemPerElem()");
        return 0;
      }

      void
      createNewElements(percept::PerceptMesh& /*eMesh*/, NodeRegistry& /*nodeRegistry*/,
                        stk::mesh::Entity /*element*/,  NewSubEntityNodesType& /*new_sub_entity_nodes*/, vector<stk::mesh::Entity>::iterator& /*element_pool*/,
                        vector<stk::mesh::Entity>::iterator& /*ft_element_pool*/,

                        stk::mesh::FieldBase */*proc_rank_field*/=0) override
      {
        throw std::runtime_error("shouldn't call URP_Tet4_Wedge6_Hex8::createNewElements()");
      }
    };

  }

#endif
