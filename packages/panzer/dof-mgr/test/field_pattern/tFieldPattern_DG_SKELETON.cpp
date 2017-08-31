// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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

// This test covers DG and SKELETON support in the
// GeometricAggFieldPattern and the AggField Pattern. The base field
// patterns don't know anything about CG/DG/SKELETON. Only the
// internal aggregates understand the differences.

namespace panzer {

  // Tests a CG only basis to make sure DG/SKELETON hasn't broken
  // things. This is purposely a Q1 basis so that the cell interior
  // dimension has no DOFs. This uncovered a bug related to assuming
  // the cell dim always had a DOF.
  TEUCHOS_UNIT_TEST(tFieldPattern_DG, HGRAD_CG)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD = rcp(new Intrepid2::Basis_HGRAD_LINE_Cn_FEM<PHX::Device,double,double>(1));
    RCP<const FieldPattern> HGRAD = rcp(new Intrepid2FieldPattern(bHGRAD));

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    // DG only
    {
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
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, HGRAD_1D_LINE)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD2 = rcp(new Intrepid2::Basis_HGRAD_LINE_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD2 = rcp(new Intrepid2FieldPattern(bHGRAD2));

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    // DG only
    {
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

    // SKELETON only
    {
      FieldVec f;
      f.push_back(std::make_tuple(0,FieldType::SKELETON,HGRAD2));
      FieldAggPattern afp(f);
      out << "HGRAD LINE SKELETON" << std::endl;
      out << "dimension = " << afp.getDimension() << std::endl;
      out << "GeomAggFP = \n" << *afp.getGeometricAggFieldPattern() << std::endl;
      out << "AggFP = \n" << afp << std::endl;
      TEST_EQUALITY(afp.getGeometricAggFieldPattern()->getSubcellIndices(1,0).size(),1);
      // Nodes
      for (int i=0; i < 2; ++i) {
        TEST_EQUALITY(afp.getSubcellIndices(0,i).size(),0);
      }
      // Cell (Edge)
      TEST_EQUALITY(afp.getSubcellIndices(1,0).size(),2);
    }

  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, HGRAD_2D_QUAD)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD2 = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD2 = rcp(new Intrepid2FieldPattern(bHGRAD2));

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    // DG only
    {
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

    // SKELETON only
    {
      FieldVec f;
      f.push_back(std::make_tuple(0,FieldType::SKELETON,HGRAD2));
      FieldAggPattern afp(f);
      out << "HGRAD QUAD SKELETON" << std::endl;
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
      TEST_EQUALITY(afp.getSubcellIndices(2,0).size(),8);
    }

  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, HGRAD_3D_HEX)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD2 = rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD2 = rcp(new Intrepid2FieldPattern(bHGRAD2));

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    // DG only
    {
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

    // SKELETON only
    {
      FieldVec f;
      f.push_back(std::make_tuple(0,FieldType::SKELETON,HGRAD2));
      FieldAggPattern afp(f);
      out << "HGRAD HEX SKELETON" << std::endl;
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
      TEST_EQUALITY(afp.getSubcellIndices(3,0).size(),26);
    }

  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, Mixed_CG_DG_SKELETON_2D)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD2 = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD2 = rcp(new Intrepid2FieldPattern(bHGRAD2));
    out << "HGRAD2\n" << *HGRAD2 << std::endl;

    RCP<Basis> bHCURL2 = rcp(new Intrepid2::Basis_HCURL_QUAD_In_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HCURL2 = rcp(new Intrepid2FieldPattern(bHCURL2));
    out << "HCURL2\n" << *HCURL2 << std::endl;

    RCP<Basis> bHDIV2 = rcp(new Intrepid2::Basis_HDIV_QUAD_In_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HDIV2 = rcp(new Intrepid2FieldPattern(bHDIV2));
    out << "HDIV2\n" << *HDIV2 << std::endl;

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    {
      FieldVec f;
      // Ordering is important. Check that all DG come after all CG
      // and that all SKELETON come after all CG+DG.
      f.push_back(std::make_tuple(0,FieldType::SKELETON,HGRAD2));
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
        TEST_EQUALITY(afp.getSubcellIndices(2,0).size(),24);
      }

      // Blocking
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
        // SKELETON
        const auto& offsets_0 = afp.localOffsets(0);
        TEST_EQUALITY(offsets_0.size(),8);
        for (auto& i : offsets_0) {
          TEST_ASSERT(i >= 24);
          TEST_ASSERT(i < 32);
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, Mixed_CG_DG_SKELETON_3D)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD2 = rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD2 = rcp(new Intrepid2FieldPattern(bHGRAD2));
    out << "HGRAD2\n" << *HGRAD2 << std::endl;

    RCP<Basis> bHCURL2 = rcp(new Intrepid2::Basis_HCURL_HEX_In_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HCURL2 = rcp(new Intrepid2FieldPattern(bHCURL2));
    out << "HCURL2\n" << *HCURL2 << std::endl;

    RCP<Basis> bHDIV2 = rcp(new Intrepid2::Basis_HDIV_HEX_In_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HDIV2 = rcp(new Intrepid2FieldPattern(bHDIV2));
    out << "HDIV2\n" << *HDIV2 << std::endl;

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    {
      FieldVec f;
      // Ordering is important. Check that all DG come after all CG
      // and that all SKELETON come after all CG+DG.
      f.push_back(std::make_tuple(0,FieldType::SKELETON,HGRAD2));
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
        TEST_EQUALITY(afp.getSubcellIndices(3,0).size(),92); // 26 CG + 54 DG + 12 SKELETON
      }

      // Blocking
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
        // SKELETON
        const auto& offsets_0 = afp.localOffsets(0);
        TEST_EQUALITY(offsets_0.size(),26);
        for (auto& i : offsets_0) {
          TEST_ASSERT(i >= 90);
          TEST_ASSERT(i < 116);
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, Mixed_CG_DG_SKELETON_INTERLEAVED)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD = rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD = rcp(new Intrepid2FieldPattern(bHGRAD));
    out << "HGRAD\n" << *HGRAD << std::endl;

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    {
      FieldVec f;
      // Ordering is important. Check that all DG come after all CG
      // and that all SKELETON come after all CG+DG.
      f.push_back(std::make_tuple(0,FieldType::SKELETON,HGRAD));
      f.push_back(std::make_tuple(1,FieldType::DG,HGRAD));
      f.push_back(std::make_tuple(2,FieldType::CG,HGRAD));
      f.push_back(std::make_tuple(3,FieldType::SKELETON,HGRAD));
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
        TEST_EQUALITY(afp.getSubcellIndices(3,0).size(),108); // 2 CG + 54 DG + 52 SKELETON
      }

      // Blocking
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
        // SKELETON
        const auto& offsets_0 = afp.localOffsets(0);
        TEST_EQUALITY(offsets_0.size(),26);
        for (auto& i : offsets_0) {
          TEST_ASSERT(i >= 108);
          TEST_ASSERT(i < 160);
        }
        const auto& offsets_3 = afp.localOffsets(3);
        TEST_EQUALITY(offsets_3.size(),26);
        for (auto& i : offsets_3) {
          TEST_ASSERT(i >= 108);
          TEST_ASSERT(i < 160);
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST(tFieldPattern_DG, Mixed_CG_DG_SKELETON_OFFSET_CLOSURE)
  {
    using Basis = Intrepid2::Basis<PHX::Device,double,double>;

    RCP<Basis> bHGRAD = rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> HGRAD = rcp(new Intrepid2FieldPattern(bHGRAD));
    out << "HGRAD\n" << *HGRAD << std::endl;

    using FieldVec = std::vector<std::tuple<int,FieldType,RCP<const FieldPattern>>>;

    {
      FieldVec f;
      // Ordering is important. Check that all DG come after all CG
      // and that all SKELETON come after all CG+DG.
      f.push_back(std::make_tuple(0,FieldType::SKELETON,HGRAD));
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
        // for (const auto& o : o1)
        //   out << "dg o1=" << o << std::endl;
        // for (const auto& o : o2)
        //   out << "dg o2=" << o << std::endl;
        
        TEST_EQUALITY(o1.size(), 9);
        TEST_EQUALITY(o1.size(), o2.size());
      }

      // Check the DG edges
      for (int i=0; i < afp.getFieldPattern(1)->getSubcellCount(1); ++i) {
        const auto& offsets = afp.localOffsets_closure(1,1,0);
        const auto& o1 = offsets.first;
        const auto& o2 = offsets.second;
        // for (const auto& o : o1)
        //   out << "dg o1=" << o << std::endl;
        // for (const auto& o : o2)
        //   out << "dg o2=" << o << std::endl;
        
        TEST_EQUALITY(o1.size(), 3);
        TEST_EQUALITY(o1.size(), o2.size());
      }

      // Check the SKELETON faces
      for (int i=0; i < afp.getFieldPattern(0)->getSubcellCount(2); ++i) {
        const auto& offsets = afp.localOffsets_closure(0,2,i);
        const auto& o1 = offsets.first;
        const auto& o2 = offsets.second;
        // for (const auto& o : o1)
        //   out << "sk o1=" << o << std::endl;
        // for (const auto& o : o2)
        //   out << "sk o2=" << o << std::endl;
        TEST_EQUALITY(o1.size(), 9);
        TEST_EQUALITY(o1.size(), o2.size());
      }

      // Check the SKELETON edges
      for (int i=0; i < afp.getFieldPattern(0)->getSubcellCount(1); ++i) {
        const auto& offsets = afp.localOffsets_closure(0,1,i);
        const auto& o1 = offsets.first;
        const auto& o2 = offsets.second;
        // for (const auto& o : o1)
        //   out << "sk o1=" << o << std::endl;
        // for (const auto& o : o2)
        //   out << "sk o2=" << o << std::endl;
        TEST_EQUALITY(o1.size(), 3);
        TEST_EQUALITY(o1.size(), o2.size());
      }
      
      

    }
  }

}
