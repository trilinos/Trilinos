// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"

#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"
#include "Intrepid2_HDIV_QUAD_In_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"
#include "Intrepid2_HDIV_HEX_In_FEM.hpp"

#include <vector>
#include <tuple>
#include <utility>

using Teuchos::rcp;
using Teuchos::RCP;
using namespace panzer;

// This test covers DG support for the GeometricAggFieldPattern and
// the FieldAggPattern. The base field patterns don't know anything
// about CG/DG. Only the internal aggregates understand the
// differences.

namespace panzer {

  // Tests a CG only basis to make sure DG hasn't broken things. This
  // is purposely a Q1 basis so that the cell interior dimension has
  // no DOFs. This uncovered a bug related to assuming the cell dim
  // always had a DOF.
  TEUCHOS_UNIT_TEST(tFieldPattern_DG, HGRAD_CG)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD = rcp(new Intrepid2::Basis_HGRAD_LINE_Cn_FEM<PHX::Device,double,double>(1));
    RCP<const FieldPattern> HGRAD = rcp(new Intrepid2FieldPattern(bHGRAD));

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;
    FieldVec f;
    f.push_back(std::make_tuple(0,FieldType::CG,HGRAD));
    FieldAggPattern afp(f);
    out << "HGRAD QUAD CG" << std::endl;
    out << "dimension = " << afp.getDimension() << std::endl;
    out << "GeomAggFP = \n" << *afp.getGeometricAggFieldPattern() << std::endl;
    out << "AggFP = \n" << afp << std::endl;
    TEST_EQUALITY(afp.getGeometricAggFieldPattern()->getSubcellIndices(1,0).size(),0);
    TEST_EQUALITY(afp.getSubcellIndices(1,0).size(),0);
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, HGRAD_1D_LINE)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD2 = rcp(new Intrepid2::Basis_HGRAD_LINE_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD2 = rcp(new Intrepid2FieldPattern(bHGRAD2));

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;
    FieldVec f;
    f.push_back(std::make_tuple(0,FieldType::DG,HGRAD2));
    FieldAggPattern afp(f);
    out << "HGRAD LINE DG" << std::endl;
    out << "dimension = " << afp.getDimension() << std::endl;
    out << "GeomAggFP = \n" << *afp.getGeometricAggFieldPattern() << std::endl;
    out << "AggFP = \n" << afp << std::endl;
    TEST_EQUALITY(afp.getGeometricAggFieldPattern()->getSubcellIndices(1,0).size(),1);
    // Nodes
    for (int i=0; i < 2; ++i) {
      TEST_EQUALITY(afp.getSubcellIndices(0,i).size(),0);
    }
    // Cell (Edge)
    TEST_EQUALITY(afp.getSubcellIndices(1,0).size(),3);
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, HGRAD_2D_QUAD)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD2 = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD2 = rcp(new Intrepid2FieldPattern(bHGRAD2));

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;
    FieldVec f;
    f.push_back(std::make_tuple(0,FieldType::DG,HGRAD2));
    FieldAggPattern afp(f);
    out << "HGRAD QUAD DG" << std::endl;
    out << "dimension = " << afp.getDimension() << std::endl;
    out << "GeomAggFP = \n" << *afp.getGeometricAggFieldPattern() << std::endl;
    out << "AggFP = \n" << afp << std::endl;
    TEST_EQUALITY(afp.getGeometricAggFieldPattern()->getSubcellIndices(2,0).size(),1);
    // Nodes
    for (int i=0; i < 4; ++i) {
      TEST_EQUALITY(afp.getSubcellIndices(0,i).size(),0);
    }
    // Edges
    for (int i=0; i < 4; ++i) {
      TEST_EQUALITY(afp.getSubcellIndices(1,i).size(),0);
    }
    // Cell (Face)
    TEST_EQUALITY(afp.getSubcellIndices(2,0).size(),9);
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, HGRAD_3D_HEX)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD2 = rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD2 = rcp(new Intrepid2FieldPattern(bHGRAD2));

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;
    FieldVec f;
    f.push_back(std::make_tuple(0,FieldType::DG,HGRAD2));
    FieldAggPattern afp(f);
    out << "HGRAD HEX DG" << std::endl;
    out << "dimension = " << afp.getDimension() << std::endl;
    out << "GeomAggFP = \n" << *afp.getGeometricAggFieldPattern() << std::endl;
    out << "AggFP = \n" << afp << std::endl;
    TEST_EQUALITY(afp.getGeometricAggFieldPattern()->getSubcellIndices(3,0).size(),1);
    // Nodes
    for (int i=0; i < 8; ++i) {
      TEST_EQUALITY(afp.getSubcellIndices(0,i).size(),0);
    }
    // Edges
    for (int i=0; i < 12; ++i) {
      TEST_EQUALITY(afp.getSubcellIndices(1,i).size(),0);
    }
    // Faces
    for (int i=0; i < 6; ++i) {
      TEST_EQUALITY(afp.getSubcellIndices(2,i).size(),0);
    }
    // Cell
    TEST_EQUALITY(afp.getSubcellIndices(3,0).size(),27);
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, Mixed_CG_DG_2D)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHCURL2 = rcp(new Intrepid2::Basis_HCURL_QUAD_In_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HCURL2 = rcp(new Intrepid2FieldPattern(bHCURL2));
    out << "HCURL2\n" << *HCURL2 << std::endl;

    RCP<Basis> bHDIV2 = rcp(new Intrepid2::Basis_HDIV_QUAD_In_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HDIV2 = rcp(new Intrepid2FieldPattern(bHDIV2));
    out << "HDIV2\n" << *HDIV2 << std::endl;

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    FieldVec f;
    f.push_back(std::make_tuple(1,FieldType::DG,HCURL2));
    f.push_back(std::make_tuple(2,FieldType::CG,HDIV2));
    FieldAggPattern afp(f);
    out << "dimension = " << afp.getDimension() << std::endl;
    out << "GeomAggFP = \n" << *afp.getGeometricAggFieldPattern() << std::endl;
    out << "AggFP = \n" << afp << std::endl;
    
    // Sizing
    {
      TEST_EQUALITY(afp.getGeometricAggFieldPattern()->getSubcellIndices(2,0).size(),1);
      // Nodes
      for (int i=0; i < 4; ++i) {
        TEST_EQUALITY(afp.getSubcellIndices(0,i).size(),0);
      }
      // Edges
      for (int i=0; i < 4; ++i) {
        TEST_EQUALITY(afp.getSubcellIndices(1,i).size(),2);
      }
      // Cell (Face)
      TEST_EQUALITY(afp.getSubcellIndices(2,0).size(),16);
    }
    
    // Check that all DG come after all CG
    {
      // CG
      const auto& offsets_2 = afp.localOffsets(2);
      TEST_EQUALITY(offsets_2.size(),12);
      for (auto& i : offsets_2) {
        TEST_ASSERT(i >= 0);
        TEST_ASSERT(i < 12);
      }
      // DG
      const auto& offsets_1 = afp.localOffsets(1);
      TEST_EQUALITY(offsets_2.size(),12);
      for (auto& i : offsets_1) {
        TEST_ASSERT(i >= 12);
        TEST_ASSERT(i < 24);
      }
    }
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, Mixed_CG_DG_3D)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHCURL2 = rcp(new Intrepid2::Basis_HCURL_HEX_In_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HCURL2 = rcp(new Intrepid2FieldPattern(bHCURL2));
    out << "HCURL2\n" << *HCURL2 << std::endl;

    RCP<Basis> bHDIV2 = rcp(new Intrepid2::Basis_HDIV_HEX_In_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HDIV2 = rcp(new Intrepid2FieldPattern(bHDIV2));
    out << "HDIV2\n" << *HDIV2 << std::endl;

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    FieldVec f;
    f.push_back(std::make_tuple(1,FieldType::DG,HCURL2));
    f.push_back(std::make_tuple(2,FieldType::CG,HDIV2));
    FieldAggPattern afp(f);
    out << "dimension = " << afp.getDimension() << std::endl;
    out << "GeomAggFP = \n" << *afp.getGeometricAggFieldPattern() << std::endl;
    out << "AggFP = \n" << afp << std::endl;
    
    // Sizing
    {
      TEST_EQUALITY(afp.getGeometricAggFieldPattern()->getSubcellIndices(3,0).size(),1);
      // Nodes
      for (int i=0; i < 8; ++i) {
        TEST_EQUALITY(afp.getSubcellIndices(0,i).size(),0);
      }
      // Edges
      for (int i=0; i < 12; ++i) {
        TEST_EQUALITY(afp.getSubcellIndices(1,i).size(),0);
      }
      // Faces
      for (int i=0; i < 6; ++i) {
        TEST_EQUALITY(afp.getSubcellIndices(2,i).size(),4);
      }
      // Cell
      TEST_EQUALITY(afp.getSubcellIndices(3,0).size(),66); // 12 CG + 54 DG
    }
    
    // Check that all DG come after all CG
    {
      // CG
      const auto& offsets_2 = afp.localOffsets(2);
      TEST_EQUALITY(offsets_2.size(),36);
      for (auto& i : offsets_2) {
        TEST_ASSERT(i >= 0);
        TEST_ASSERT(i < 36);
      }
      // DG
      const auto& offsets_1 = afp.localOffsets(1);
      TEST_EQUALITY(offsets_1.size(),54);
      for (auto& i : offsets_1) {
        TEST_ASSERT(i >= 36);
        TEST_ASSERT(i < 90);
      }
    }
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, Mixed_CG_DG_INTERLEAVED)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD = rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD = rcp(new Intrepid2FieldPattern(bHGRAD));
    out << "HGRAD\n" << *HGRAD << std::endl;

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    FieldVec f;
    f.push_back(std::make_tuple(1,FieldType::DG,HGRAD));
    f.push_back(std::make_tuple(2,FieldType::CG,HGRAD));
    f.push_back(std::make_tuple(4,FieldType::DG,HGRAD));
    f.push_back(std::make_tuple(5,FieldType::CG,HGRAD));
    FieldAggPattern afp(f);
    out << "dimension = " << afp.getDimension() << std::endl;
    out << "GeomAggFP = \n" << *afp.getGeometricAggFieldPattern() << std::endl;
    out << "AggFP = \n" << afp << std::endl;
    // Sizes
    {
      TEST_EQUALITY(afp.getGeometricAggFieldPattern()->getSubcellIndices(3,0).size(),1);
      // Nodes
      for (int i=0; i < 8; ++i) {
        TEST_EQUALITY(afp.getSubcellIndices(0,i).size(),2);
      }
      // Edges
      for (int i=0; i < 12; ++i) {
        TEST_EQUALITY(afp.getSubcellIndices(1,i).size(),2);
      }
      // Faces
      for (int i=0; i < 6; ++i) {
        TEST_EQUALITY(afp.getSubcellIndices(2,i).size(),2);
      }
      // Cell
      TEST_EQUALITY(afp.getSubcellIndices(3,0).size(),56); // 2 CG + 54 DG
    }
    
    // Check that all DG come after all CG
    {
      // CG
      const auto& offsets_2 = afp.localOffsets(2);
      TEST_EQUALITY(offsets_2.size(),27);
      for (auto& i : offsets_2) {
        TEST_ASSERT(i >= 0);
        TEST_ASSERT(i < 54);
      }
      const auto& offsets_5 = afp.localOffsets(5);
      TEST_EQUALITY(offsets_5.size(),27);
      for (auto& i : offsets_5) {
        TEST_ASSERT(i >= 0);
        TEST_ASSERT(i < 54);
      }
      // DG
      const auto& offsets_1 = afp.localOffsets(1);
      TEST_EQUALITY(offsets_1.size(),27);
      for (auto& i : offsets_1) {
        TEST_ASSERT(i >= 54);
        TEST_ASSERT(i < 108);
      }
      const auto& offsets_4 = afp.localOffsets(4);
      TEST_EQUALITY(offsets_4.size(),27);
      for (auto& i : offsets_4) {
        TEST_ASSERT(i >= 54);
        TEST_ASSERT(i < 108);
      }
    }
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, Mixed_CG_DG_OFFSET_CLOSURE)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD = rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD = rcp(new Intrepid2FieldPattern(bHGRAD));
    out << "HGRAD\n" << *HGRAD << std::endl;

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    FieldVec f;
    f.push_back(std::make_tuple(1,FieldType::DG,HGRAD));
    f.push_back(std::make_tuple(2,FieldType::CG,HGRAD));
    FieldAggPattern afp(f);
    out << "dimension = " << afp.getDimension() << std::endl;
    out << "GeomAggFP = \n" << *afp.getGeometricAggFieldPattern() << std::endl;
    out << "AggFP = \n" << afp << std::endl;
    
    // Check the DG faces
    for (int i=0; i < afp.getFieldPattern(1)->getSubcellCount(2); ++i) {
      const auto& offsets = afp.localOffsets_closure(1,2,0);
      const auto& o1 = offsets.first;
      const auto& o2 = offsets.second;      
      TEST_EQUALITY(o1.size(), 9);
      TEST_EQUALITY(o1.size(), o2.size());
    }

    // Check the DG edges
    for (int i=0; i < afp.getFieldPattern(1)->getSubcellCount(1); ++i) {
      const auto& offsets = afp.localOffsets_closure(1,1,0);
      const auto& o1 = offsets.first;
      const auto& o2 = offsets.second;
      TEST_EQUALITY(o1.size(), 3);
      TEST_EQUALITY(o1.size(), o2.size());
    }

    // Check that the offset numbering is correct and unique. NOTE:
    // this check uses the knowledge of the underlying numbering in
    // shards topology and intrepid2 bases/field pattern. If either
    // library changes the underlying ordering, this check will break.

    // Offsets for a DOF from HGRAD pattern. The internal
    // implementation orders the dofs based on increasing subcell
    // indices.
    std::vector<int> v2_gold(27);
    v2_gold[0] = 0;
    v2_gold[1] = 8;
    v2_gold[2] = 1;
    v2_gold[3] = 11;
    v2_gold[4] = 24;
    v2_gold[5] = 9;
    v2_gold[6] = 3;
    v2_gold[7] = 10;    
    v2_gold[8]  = 2;
    v2_gold[9]  = 16;
    v2_gold[10] = 20;
    v2_gold[11] = 17;
    v2_gold[12] = 23;
    v2_gold[13] = 26;
    v2_gold[14] = 21;
    v2_gold[15] = 19;
    v2_gold[16] = 22;
    v2_gold[17] = 18;
    v2_gold[18] = 4;
    v2_gold[19] = 12;
    v2_gold[20] = 5;
    v2_gold[21] = 15;
    v2_gold[22] = 25; 
    v2_gold[23] = 13;
    v2_gold[24] = 7;
    v2_gold[25] = 14;
    v2_gold[26] = 6;

    // v2 is CG field
    const auto& v2 = afp.localOffsets(2);    
    for(std::size_t i=0;i<v2.size();i++)
      TEST_EQUALITY(v2[i],v2_gold[i]);

    // v1 is DG field. Same pattern as v2 but offset by number of v2
    // DOFs.
    const auto& v1 = afp.localOffsets(1);    
    for(std::size_t i=0;i<v1.size();i++)
      TEST_EQUALITY(v1[i],v2_gold[i]+27);
  }

}
