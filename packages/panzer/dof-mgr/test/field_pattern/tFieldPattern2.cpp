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

#include "Panzer_FieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_EdgeFieldPattern.hpp"
#include "Panzer_FaceFieldPattern.hpp"
#include "Panzer_ElemFieldPattern.hpp"

#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

  std::string note = "***   NOTE: UNIT TEST BASED ON SEPT 2010   ***\n"
                     "***   INTREPID AND SHARDS Trilinos-dev     ***\n"
                     "***   DOXYGEN WEBSITE                      ***\n";

  /////////////////////////////////////////////
  // 2D tests
  /////////////////////////////////////////////

  bool intrepid_equals(const RCP<Intrepid2::Basis<PHX::Device,double,double> > & basisA,
                       const RCP<Intrepid2::Basis<PHX::Device,double,double> > & basisB,
                       const std::string & file,int lineNo);
  bool intrepid_same_geom(const RCP<Intrepid2::Basis<PHX::Device,double,double> > & basisA,
                          const RCP<Intrepid2::Basis<PHX::Device,double,double> > & basisB,
                          const std::string & file,int lineNo);

  // triangle tests
  TEUCHOS_UNIT_TEST(tFieldPattern, test_equals)
  {
    out << note << std::endl;

    RCP<Intrepid2::Basis<PHX::Device,double,double> > basisA;
    RCP<Intrepid2::Basis<PHX::Device,double,double> > basisB;

    // test for same order
    basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(4));
    basisB = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(4));

    TEST_ASSERT(intrepid_same_geom(basisA,basisB,__FILE__,__LINE__));
    TEST_ASSERT(intrepid_equals(basisA,basisB,__FILE__,__LINE__));

    // test for different  order
    basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(4));
    basisB = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(3));

    TEST_ASSERT(intrepid_same_geom(basisA,basisB,__FILE__,__LINE__));
    TEST_ASSERT(not intrepid_equals(basisA,basisB,__FILE__,__LINE__));

    // test basis accessor
    panzer::Intrepid2FieldPattern fp(basisA);
    TEST_ASSERT(basisA.get() == fp.getIntrepidBasis().get());
  }

  // Checks that the clone function and compare operator are
  // functioning as expected.
  TEUCHOS_UNIT_TEST(tFieldPattern, CloneAndEqualsOP)
  {
    // Build nodal Q1 and Q2
    auto basis1 = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(1));
    auto basis2 = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(2));
    RCP<const FieldPattern> fp1 = rcp(new Intrepid2FieldPattern(basis1));
    RCP<const FieldPattern> fp2 = rcp(new Intrepid2FieldPattern(basis2));

    TEST_ASSERT(*fp1 != *fp2);
    auto clone1 = fp1->clone();
    TEST_ASSERT(fp1 != clone1); // pointers not equal
    TEST_ASSERT(*fp1 == *clone1); // values equal

    // Check differences between Q1 and Q2
    TEST_EQUALITY(fp1->getSubcellCount(0),4);
    TEST_EQUALITY(fp1->getSubcellCount(1),4);
    TEST_EQUALITY(fp1->getSubcellCount(2),1);
    TEST_EQUALITY(fp1->getSubcellCount(3),0);
    TEST_EQUALITY(fp2->getSubcellCount(0),4);
    TEST_EQUALITY(fp2->getSubcellCount(1),4);
    TEST_EQUALITY(fp2->getSubcellCount(2),1);
    TEST_EQUALITY(fp2->getSubcellCount(3),0);

    TEST_EQUALITY(fp1->getSubcellIndices(1,0).size(),0); // no edge dof
    TEST_EQUALITY(fp2->getSubcellIndices(1,0).size(),1); // one edge dof

    std::vector<std::pair<panzer::FieldType,RCP<const panzer::FieldPattern>>> geometric_patterns{{panzer::FieldType::CG,fp1},
                                                                                                 {panzer::FieldType::CG,fp2}};
    RCP<const FieldPattern> fp_geom_agg = Teuchos::make_rcp<panzer::GeometricAggFieldPattern>(geometric_patterns);
    TEST_EQUALITY(fp_geom_agg->getSubcellIndices(1,0).size(),1); // one edge dof
    auto fpgeom_agg_clone = fp_geom_agg->clone();
    TEST_ASSERT(fp_geom_agg != fpgeom_agg_clone); // pointers not equal
    TEST_ASSERT(*fp_geom_agg == *fpgeom_agg_clone); // values equal

    std::vector<std::tuple<int,panzer::FieldType,RCP<const panzer::FieldPattern>>> field_patterns{{0,panzer::FieldType::CG,fp1},
                                                                                                  {1,panzer::FieldType::CG,fp2},
                                                                                                  {2,panzer::FieldType::CG,fp2}};
    RCP<const FieldPattern> fp_agg = Teuchos::make_rcp<panzer::FieldAggPattern>(field_patterns);
    TEST_EQUALITY(fp_agg->getSubcellIndices(0,0).size(),3); // three vertex dof
    TEST_EQUALITY(fp_agg->getSubcellIndices(1,0).size(),2); // two edge dof

    auto fp_geom_agg2 = Teuchos::rcp_dynamic_cast<const panzer::FieldAggPattern>(fp_agg)->getGeometricAggFieldPattern();
    TEST_ASSERT(fp_geom_agg2 != fpgeom_agg_clone); // pointers not equal
    TEST_ASSERT(*fp_geom_agg2 == *fpgeom_agg_clone); // values equal

    // Let's make some aggregates that are not the same as above
    std::vector<std::pair<panzer::FieldType,RCP<const panzer::FieldPattern>>> geometric_patterns_alt{{panzer::FieldType::CG,fp1}};
    RCP<const FieldPattern> fp_geom_agg_alt = Teuchos::make_rcp<panzer::GeometricAggFieldPattern>(geometric_patterns_alt);
    TEST_ASSERT(*fp_geom_agg_alt != *fpgeom_agg_clone); // values not equal

    std::vector<std::tuple<int,panzer::FieldType,RCP<const panzer::FieldPattern>>> field_patterns_alt{{0,panzer::FieldType::CG,fp2}};
    RCP<const FieldPattern> fp_agg_alt = Teuchos::make_rcp<panzer::FieldAggPattern>(field_patterns_alt);
    TEST_INEQUALITY(*fp_agg,*fp_agg_alt);

    // Geometric agg for both fp_agg are the same
    TEST_EQUALITY(*Teuchos::rcp_dynamic_cast<const panzer::FieldAggPattern>(fp_agg)->getGeometricAggFieldPattern(),
                  *Teuchos::rcp_dynamic_cast<const panzer::FieldAggPattern>(fp_agg_alt)->getGeometricAggFieldPattern());

    RCP<const FieldPattern> nodal_fp = Teuchos::make_rcp<panzer::NodalFieldPattern>(basis1->getBaseCellTopology());
    auto nodal_fp_clone = nodal_fp->clone();
    TEST_ASSERT(nodal_fp != nodal_fp_clone);
    TEST_ASSERT(*nodal_fp == *nodal_fp_clone);

    auto edge_fp = Teuchos::make_rcp<panzer::EdgeFieldPattern>(basis1->getBaseCellTopology());
    auto edge_fp_clone = edge_fp->clone();
    TEST_ASSERT(edge_fp != edge_fp_clone);
    TEST_ASSERT(*edge_fp == *edge_fp_clone);

    auto face_fp = Teuchos::make_rcp<panzer::FaceFieldPattern>(basis1->getBaseCellTopology());
    auto face_fp_clone = face_fp->clone();
    TEST_ASSERT(face_fp != face_fp_clone);
    TEST_ASSERT(*face_fp == *face_fp_clone);

    auto elem_fp = Teuchos::make_rcp<panzer::ElemFieldPattern>(basis1->getBaseCellTopology());
    auto elem_fp_clone = elem_fp->clone();
    TEST_ASSERT(elem_fp != elem_fp_clone);
    TEST_ASSERT(*elem_fp == *elem_fp_clone);
  }

  bool intrepid_equals(const RCP<Intrepid2::Basis<PHX::Device,double,double> > & basisA,
                       const RCP<Intrepid2::Basis<PHX::Device,double,double> > & basisB,
                       const std::string & file,int lineNo)
  {
    // notice if Intrepid2FieldPattern implements "equals" then this functionality
    // is not quite correct!
    Teuchos::RCP<const FieldPattern> fpA = rcp(new Intrepid2FieldPattern(basisA));
    Teuchos::RCP<const FieldPattern> fpB= rcp(new Intrepid2FieldPattern(basisB));

    bool forward = fpA->equals(*fpB);
    bool bckward = fpB->equals(*fpA);

    std::cout << " intrepid_equals \n";
    fpA->print(std::cout);
    fpB->print(std::cout);

    if(bckward!=forward) { // if true, this is a big problem
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Test \"intrepid_equals\" leads to an"
                                 "inconsistency with the FieldPattern::equals "
                                 "functionality: " << file << ": " << lineNo );
    }
    TEUCHOS_ASSERT(bckward==forward);

    return forward; // they are equal so that it shouldn't matter at this point
  }

  bool intrepid_same_geom(const RCP<Intrepid2::Basis<PHX::Device,double,double> > & basisA,
                          const RCP<Intrepid2::Basis<PHX::Device,double,double> > & basisB,
                          const std::string & file,int lineNo)
  {
    // notice if Intrepid2FieldPattern implements "sameGeometry" then this functionality
    // is not quite correct!
    Teuchos::RCP<const FieldPattern> fpA = rcp(new Intrepid2FieldPattern(basisA));
    Teuchos::RCP<const FieldPattern> fpB= rcp(new Intrepid2FieldPattern(basisB));

    bool forward = fpA->sameGeometry(*fpB);
    bool bckward = fpB->sameGeometry(*fpA);

    std::cout << " intrepid_geom_equals \n";
    fpA->print(std::cout);
    fpB->print(std::cout);

    if(bckward!=forward) { // if true, this is a big problem
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Test \"intrepid_samegeom\" leads to an"
                                 "inconsistency with the FieldPattern::sameGeometry "
                                 "functionality: " << file << ": " << lineNo );
    }
    TEUCHOS_ASSERT(bckward==forward);

    return forward; // they are equal so that it shouldn't matter at this point
  }

}
