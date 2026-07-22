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

  // tests:
  //    - getGeometricAggPattern
  //    - getFieldPattern
  //    - getDimension 
  TEUCHOS_UNIT_TEST(tFieldAggPattern, testA)
  {

    out << note << std::endl;

    // basis to build patterns from
    const int maxOrder = Intrepid2::Parameters::MaxOrder;
    const int basisA_degree = std::min(3, maxOrder);
    const int basisB_degree = std::min(5, maxOrder);
    const int basisC_degree = std::min(1, maxOrder);
    RCP<Intrepid2::Basis<PHX::Device,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(basisA_degree));
    RCP<Intrepid2::Basis<PHX::Device,double,double> > basisB = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(basisB_degree));
    RCP<Intrepid2::Basis<PHX::Device,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(basisC_degree));

    RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
    RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
    RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisC));

    std::vector<int> closureIndices;
    std::vector<std::pair<FieldType,RCP<const FieldPattern>>> patternV;
    std::vector<std::tuple<int,panzer::FieldType,RCP<const FieldPattern> > > patternM;

    patternV.push_back(std::make_pair(FieldType::CG,patternA));
    patternV.push_back(std::make_pair(FieldType::CG,patternB));
    patternV.push_back(std::make_pair(FieldType::CG,patternC));

    GeometricAggFieldPattern geom(patternV);

    patternM.push_back(std::make_tuple(7,FieldType::CG,patternA));
    patternM.push_back(std::make_tuple(3,FieldType::CG,patternB));
    patternM.push_back(std::make_tuple(4,FieldType::CG,patternC));

    // test build of geometric field pattern
    {
      FieldAggPattern agg; 
   
      TEST_THROW(agg.getSubcellClosureIndices(2,0,closureIndices),std::logic_error);
      // TEST_THROW(*agg.getGeometricAggFieldPattern(),Teuchos::NullReferenceError);

      TEST_EQUALITY(agg.getGeometricAggFieldPattern(),Teuchos::null);

      agg.buildPattern(patternM);

      bool equality = false;
      TEST_NOTHROW(equality = geom.equals(*agg.getGeometricAggFieldPattern()));
      TEST_ASSERT(equality);

      TEST_EQUALITY(geom.getDimension(),agg.getGeometricAggFieldPattern()->getDimension());
    }

    // test build of geometric field pattern
    {
      FieldAggPattern agg(patternM); 

      TEST_THROW(agg.getSubcellClosureIndices(2,0,closureIndices),std::logic_error);
   
      bool equality = false;
      TEST_NOTHROW(equality = geom.equals(*agg.getGeometricAggFieldPattern()));
      TEST_ASSERT(equality);

      TEST_EQUALITY(geom.getDimension(),agg.getGeometricAggFieldPattern()->getDimension());
    }

    {
      FieldAggPattern agg(patternM); 
   
      TEST_THROW(agg.getFieldPattern(5),std::logic_error); // no field 5

      TEST_NOTHROW(agg.getFieldPattern(3));
      TEST_NOTHROW(agg.getFieldPattern(7));
      TEST_NOTHROW(agg.getFieldPattern(4));

      TEST_EQUALITY(agg.getFieldPattern(3),patternB);
      TEST_EQUALITY(agg.getFieldPattern(7),patternA);
      TEST_EQUALITY(agg.getFieldPattern(4),patternC);
    }
  }

  // tests:
  //    - buildPattern --> different and same geometries
  TEUCHOS_UNIT_TEST(tFieldAggPattern, testB)
  {
    out << note << std::endl;
   
    const int maxOrder = Intrepid2::Parameters::MaxOrder; 
    const int poly_A = std::min(3, maxOrder),
              poly_B = std::min(5, maxOrder);
    
    // basis to build patterns from
    RCP<Intrepid2::Basis<PHX::Device,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(poly_A));
    RCP<Intrepid2::Basis<PHX::Device,double,double> > basisB = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(poly_B));

    const int ndof_A[] = { basisA->getDofTag(basisA->getDofOrdinal(0,0,0))(3),
                           basisA->getDofTag(basisA->getDofOrdinal(1,0,0))(3),
                           basisA->getDofTag(basisA->getDofOrdinal(2,0,0))(3) };
    const int ndof_B[] = { basisB->getDofTag(basisB->getDofOrdinal(0,0,0))(3),
                           basisB->getDofTag(basisB->getDofOrdinal(1,0,0))(3),
                           basisB->getDofTag(basisB->getDofOrdinal(2,0,0))(3) };

    RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
    RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));

    // test that single construction gives the same pattern
    {
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,patternB));

      FieldAggPattern agg(patternM); 
      agg.print(out); // for debugging purposes

      TEST_EQUALITY(agg.getDimension(),patternB->getDimension()); 
      TEST_EQUALITY((int) agg.fieldIds().size(),patternB->numberIds());

      // this cannot match in high order case
      // TEST_EQUALITY((int) agg.numFieldsPerId().size(),patternB->numberIds());

      const std::vector<int> & fieldsPerId = agg.numFieldsPerId();
      const std::vector<int> & fieldIds = agg.fieldIds();

      int numberFields = 0;
      for(std::size_t i=0;i<fieldsPerId.size();i++)
        numberFields += fieldsPerId[i];

      TEST_EQUALITY((int) fieldIds.size(),numberFields);


      // check fields per ID
      {
        bool allMatch = true;
        const int fields_B[] = { ndof_B[0], ndof_B[0], ndof_B[0], ndof_B[0], 
                                 ndof_B[1], ndof_B[1], ndof_B[1], ndof_B[1], 
                                 ndof_B[2] };
        int size_fields_B = 0;
        for (int i=0;i<9;++i) 
          size_fields_B += (fields_B[i] > 0);

        TEST_ASSERT(static_cast<int>(fieldsPerId.size()) == size_fields_B);
        for(std::size_t i=0;i<fieldsPerId.size();i++) {
          allMatch &= (fieldsPerId[i]==fields_B[i]);
        }
        TEST_ASSERT(allMatch);
      }

      // check fields ids
      {
        bool allIds = true;
        for(std::size_t i=0;i<fieldIds.size();i++)
          allIds &= (fieldIds[i]==3);
        TEST_ASSERT(allIds);
      }

      for(int i=0;i<agg.getDimension();i++)
        TEST_EQUALITY(agg.getSubcellCount(i),patternB->getSubcellCount(i)); 
      // this cannot be compared anymore
      // agg create a pattern based on the geometrical order while patternB is arbitrary given from intrepid
      // TEST_ASSERT(agg.equals(*patternB)); 
    }

    // test that single construction gives the same pattern
    {
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,patternB));
      patternM.push_back(std::make_tuple(7,FieldType::CG,patternA));

      FieldAggPattern agg(patternM); 
      agg.print(out); // for debugging purposes

      TEST_ASSERT(patternB->sameGeometry(agg));
      TEST_EQUALITY(agg.getDimension(),patternB->getDimension()); 
      TEST_EQUALITY((int) agg.numFieldsPerId().size(),agg.getGeometricAggFieldPattern()->numberIds());

      const std::vector<int> & fieldsPerId = agg.numFieldsPerId();
      const std::vector<int> & fieldIds = agg.fieldIds();

      int numberFields = 0;
      for(std::size_t i=0;i<fieldsPerId.size();i++)
        numberFields += fieldsPerId[i];

      TEST_EQUALITY((int) fieldIds.size(),numberFields);
      TEST_EQUALITY((int) fieldsPerId.size(),agg.getGeometricAggFieldPattern()->numberIds());

      RCP<const FieldPattern> geomPattern = agg.getGeometricAggFieldPattern();

      // check out fieldsPerId
      for (int offset=0, dimension=0;dimension<3;++dimension) {
        const int nsubcell = geomPattern->getSubcellCount(dimension);
        for(int entities=0;entities<nsubcell;entities++) {
          TEST_EQUALITY(fieldsPerId[offset],ndof_A[dimension]+ndof_B[dimension]); 
        }
        offset += nsubcell;
      }

      // this tests length of fieldsPerId is not longer than expected
      TEST_EQUALITY(geomPattern->getSubcellCount(0)+ 
                    geomPattern->getSubcellCount(1)+ 
                    geomPattern->getSubcellCount(2),
                    (int) fieldsPerId.size());

      // check out fieldIds
      for (int offset=0, dimension=0;dimension<3;++dimension) {
        const int nsubcell = geomPattern->getSubcellCount(dimension);
        for(int entities=0;entities<nsubcell;entities++) {
          for (int ndof=0;ndof<ndof_B[dimension];++ndof,++offset) {
            TEST_EQUALITY(fieldIds[offset],3);
          }
          for (int ndof=0;ndof<ndof_A[dimension];++ndof,++offset) {
            TEST_EQUALITY(fieldIds[offset],7); 
          }
        }
      }

      // this tests length of fieldIds is not longer than expected
      TEST_EQUALITY((ndof_A[0]+ndof_B[0])*geomPattern->getSubcellCount(0)+ 
                    (ndof_A[1]+ndof_B[1])*geomPattern->getSubcellCount(1)+ 
                    (ndof_A[2]+ndof_B[2])*geomPattern->getSubcellCount(2),
                    (int) fieldIds.size());

      // check out subcell information
      for (int index=0, dimension=0;dimension<3;++dimension) {
        for(int nodes=0;nodes<geomPattern->getSubcellCount(dimension);++nodes) {
          const std::vector<int> & indices = agg.getSubcellIndices(dimension,nodes);
          TEST_EQUALITY((int) indices.size(),ndof_A[dimension]+ndof_B[dimension]);
          for(std::size_t i=0;i<indices.size();i++,index++) 
            TEST_EQUALITY(indices[i],index);
        }
      }

    }
  }

  // tests:
  //    - localOffsets --> different and same geometries
  TEUCHOS_UNIT_TEST(tFieldAggPattern, testC)
  {
    out << note << std::endl;

    // basis to build patterns from

    {

      const int maxOrder = Intrepid2::Parameters::MaxOrder;
      const int poly = std::min(4,maxOrder); 
      RCP<Intrepid2::Basis<PHX::Device,double,double> > basis = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(poly));

      const int ndof[] = { basis->getDofTag(basis->getDofOrdinal(0,0,0))(3),
                           basis->getDofTag(basis->getDofOrdinal(1,0,0))(3),
                           basis->getDofTag(basis->getDofOrdinal(2,0,0))(3) };
      
      RCP<const FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,pattern));

      FieldAggPattern agg(patternM);
      agg.print(out);

      TEST_THROW(agg.localOffsets(7),std::logic_error); // check for failure if ID does not exist

      RCP<const FieldPattern> geomPattern = agg.getGeometricAggFieldPattern();

      std::vector<int> v_true, v_tmp(basis->getCardinality());
      for (int offset=0, dimension=0;dimension<3;++dimension) {
        const int nsubcell = geomPattern->getSubcellCount(dimension);
        for (int entities=0;entities<nsubcell;entities++) 
          for (int k=0;k<ndof[dimension];++k,++offset)
            v_tmp[basis->getDofOrdinal(dimension,entities,k)] = offset;
      }
      v_true = v_tmp;
      
      std::vector<int> v = agg.localOffsets(3);
      TEST_EQUALITY(v.size(),v_true.size());
      for(std::size_t i=0;i<v.size();i++)
        TEST_EQUALITY(v[i],v_true[i]);
    }

    {
      const int maxOrder = Intrepid2::Parameters::MaxOrder; 
      const int polyC1 = std::min(3, maxOrder),
                polyC2 = std::min(4, maxOrder);
      RCP<Intrepid2::Basis<PHX::Device,double,double> > basisC1 = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(polyC1));
      RCP<Intrepid2::Basis<PHX::Device,double,double> > basisC2 = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(polyC2));
      
      const int ndofC1[] = { basisC1->getDofTag(basisC1->getDofOrdinal(0,0,0))(3),
                             basisC1->getDofTag(basisC1->getDofOrdinal(1,0,0))(3),
                             basisC1->getDofTag(basisC1->getDofOrdinal(2,0,0))(3) };
      const int ndofC2[] = { basisC2->getDofTag(basisC2->getDofOrdinal(0,0,0))(3),
                             basisC2->getDofTag(basisC2->getDofOrdinal(1,0,0))(3),
                             basisC2->getDofTag(basisC2->getDofOrdinal(2,0,0))(3) };

      RCP<const FieldPattern> patternC1 = rcp(new Intrepid2FieldPattern(basisC1));
      RCP<const FieldPattern> patternC2 = rcp(new Intrepid2FieldPattern(basisC2));

      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(0,FieldType::CG,patternC2));
      patternM.push_back(std::make_tuple(1,FieldType::CG,patternC2));
      patternM.push_back(std::make_tuple(2,FieldType::CG,patternC2));
      patternM.push_back(std::make_tuple(3,FieldType::CG,patternC1));

      FieldAggPattern agg(patternM);
      agg.print(out);

      TEST_THROW(agg.localOffsets(7),std::logic_error); // check for failure if ID does not exist

      RCP<const FieldPattern> geomPattern = agg.getGeometricAggFieldPattern();

      {
        out << "\nTesting Q4-Q2 NS: Pressure" << std::endl;
        std::vector<int> v_true, v_tmp(basisC1->getCardinality());
        for (int offset=0, dimension=0;dimension<3;++dimension) {
          const int nsubcell = geomPattern->getSubcellCount(dimension);
          for (int entities=0;entities<nsubcell;entities++) {
            offset += (3*ndofC2[dimension]);
            for (int k=0;k<ndofC1[dimension];++k,++offset)
              v_tmp[basisC1->getDofOrdinal(dimension,entities,k)] = offset;
          }
        }
        v_true = v_tmp;

        std::vector<int> v = agg.localOffsets(3);
        TEST_EQUALITY(v.size(),v_true.size());
        for(std::size_t i=0;i<v.size();i++)
          TEST_EQUALITY(v[i],v_true[i]);
      }

      // for(int d=0;d<3;d++) {
      //   out << "\nTesting Q4-Q2 NS: Velocity d=" << d << std::endl;
      //   std::vector<int> v_true, v_tmp(basisC2->getCardinality());
      //   for (int offset=0, dimension=0;dimension<3;++dimension) {
      //     const int nsubcell = geomPattern->getSubcellCount(dimension);
      //     for (int entities=0;entities<nsubcell;entities++) {
      //       offset += (d*ndofC2[dimension]);
      //       for (int k=0;k<ndofC1[dimension];++k,++offset)
      //         v_tmp[basisC1->getDofOrdinal(dimension,entities,k)] = offset;
      //       offset += (2-d)*ndofC2[dimension];
      //       offset += ndofC1[dimension];
      //     }
      //   }
      //   v_true = v_tmp;
        
      //   std::vector<int> v = agg.localOffsets(d);
      //   TEST_EQUALITY(v.size(),v_true.size());
      //   for(std::size_t i=0;i<v.size();i++)
      //     TEST_EQUALITY(v[i],v_true[i]);
      // }
    }
    
    //    {
    //       out << "Crazy HEX basis test" << std::endl;
 
    //       RCP<Intrepid2::Basis<double,FieldContainer> > basisC1 = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer>);
    //       RCP<Intrepid2::Basis<double,FieldContainer> > basisDivI1 = rcp(new Intrepid2::Basis_HDIV_HEX_I1_FEM<double,FieldContainer>);
    //       RCP<Intrepid2::Basis<double,FieldContainer> > basisCurlI1 = rcp(new Intrepid2::Basis_HCURL_HEX_I1_FEM<double,FieldContainer>);
    //       RCP<const FieldPattern> patternC1 = rcp(new Intrepid2FieldPattern(basisC1));
    //       RCP<const FieldPattern> patternDivI1 = rcp(new Intrepid2FieldPattern(basisDivI1));
    //       RCP<const FieldPattern> patternCurlI1 = rcp(new Intrepid2FieldPattern(basisCurlI1));
    //       std::vector<std::pair<int,RCP<const FieldPattern> > > patternM;
    //       patternM.push_back(std::make_pair(3,patternC1));
    //       patternM.push_back(std::make_pair(9,patternDivI1));
    //       patternM.push_back(std::make_pair(1,patternCurlI1));

    //       FieldAggPattern agg(patternM);
    //       agg.print(out);

    //       TEST_THROW(agg.localOffsets(7),std::logic_error); // check for failure if ID does not exist

    //       {
    //          out << "\nTesting C1" << std::endl;
    //          std::vector<int> v = agg.localOffsets(3);
    //          for(std::size_t i=0;i<v.size();i++)
    //             TEST_EQUALITY(v[i],(int) i);
    //       }

    //       {
    //          out << "\nTesting Div-I1" << std::endl;
    //          std::vector<int> v = agg.localOffsets(9);
    //          for(std::size_t i=0;i<v.size();i++)
    //             TEST_EQUALITY(v[i],(int) i+20);
    //       }

    //       {
    //          out << "\nTesting Curl-I1" << std::endl;
    //          std::vector<int> v = agg.localOffsets(1);
    //          for(std::size_t i=0;i<v.size();i++)
    //             TEST_EQUALITY(v[i],(int) i+8);
    //       }
    //    }

    //    {
    //       out << "Highorder HEX basis test: FAILS!!!! DISABLED FOR NOW!!!!!" << std::endl;
 
    // /*
    //       RCP<Intrepid2::Basis<double,FieldContainer> > basis = rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<double,FieldContainer>(4,Intrepid2::POINTTYPE_EQUISPACED));
    //       RCP<const FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));
    //       std::vector<std::pair<int,RCP<const FieldPattern> > > patternM;
    //       patternM.push_back(std::make_pair(3,pattern));

    //       FieldAggPattern agg(patternM);
    //       out << "4th order hex" << std::endl;
    //       pattern->print(out);

    //       out << "Aggregated hex" << std::endl;
    //       agg.print(out);

    //       out << "Geometric hex" << std::endl;
    //       agg.getGeometricAggFieldPattern()->print(out);

    //       TEST_THROW(agg.localOffsets(7),std::logic_error); // check for failure if ID does not exist
    // */
    //    }
  }

  // // tests:
  // //    - localOffsets_closure --> different and same geometries
  // TEUCHOS_UNIT_TEST(tFieldAggPattern, testD)
  // {
  //    {
  //       RCP<Intrepid2::Basis<double,FieldContainer> > basis = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<double,FieldContainer>);
  //       RCP<const FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));
  //       std::vector<std::pair<int,RCP<const FieldPattern> > > patternM;
  //       patternM.push_back(std::make_pair(3,pattern));

  //       FieldAggPattern agg(patternM);
  //       agg.print(out);

  //       out << "Testing throws of localOffsets_closure" << std::endl;
  //       TEST_THROW(agg.localOffsets_closure(7,0,1),std::logic_error); // check for failure if ID does not exist
  //       TEST_THROW(agg.localOffsets_closure(3,2,6),std::logic_error); // check for failure if sub cell doesn't exist

  //       std::vector<int> v_true;
  //       std::vector<int> v;

  //       // test (2,0) closure
  //       out << "Testing localOffsets_closure(3,2,0) on a C2 HEX" << std::endl;
  //       v_true.resize(9);
  //       v = agg.localOffsets_closure(3,2,0).first;

  //       // nodes
  //       v_true[0] = 0; v_true[1] = 1; v_true[2] = 4; v_true[3] = 5;
 
  //       // edges     
  //       v_true[4]  = 8; v_true[5]  = 17; v_true[6] = 12; v_true[7] = 16;

  //       // areas
  //       v_true[8] = 20;

  //       TEST_EQUALITY(v.size(),v_true.size());

  //       std::sort(v.begin(),v.end());
  //       std::sort(v_true.begin(),v_true.end());
  //       for(std::size_t i=0;i<v.size();i++)
  //          TEST_EQUALITY(v[i],v_true[i]);

  //       // test (1,7) closure
  //       out << "Testing localOffsets_closure(3,1,7) on a C2 HEX" << std::endl;
  //       v_true.resize(3);
  //       v = agg.localOffsets_closure(3,1,7).first;

  //       // nodes
  //       v_true[0] = 4; v_true[1] = 7;
 
  //       // edges     
  //       v_true[2]  = 15; 

  //       TEST_EQUALITY(v.size(),v_true.size());

  //       std::sort(v.begin(),v.end());
  //       std::sort(v_true.begin(),v_true.end());
  //       for(std::size_t i=0;i<v.size();i++)
  //          TEST_EQUALITY(v[i],v_true[i]);

  //       // test all dimension zero closures
  //       out << "Testing localOffsets_closure(3,0,*) on a C2 HEX" << std::endl;
  //       v_true.resize(1);
  //       for(int i=0;i<8;i++) {
  //          v = agg.localOffsets_closure(3,0,i).first;

  //          // nodes
  //          v_true[0] = i;
 
  //          TEST_EQUALITY(v.size(),v_true.size());
  //          TEST_EQUALITY(v[0],v_true[0]);
  //       }
  //    }

  //    {
  //       RCP<Intrepid2::Basis<double,FieldContainer> > basisC1 = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer>);
  //       RCP<Intrepid2::Basis<double,FieldContainer> > basisC2 = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<double,FieldContainer>);
  //       RCP<const FieldPattern> patternP = rcp(new Intrepid2FieldPattern(basisC1));
  //       RCP<const FieldPattern> patternU = rcp(new Intrepid2FieldPattern(basisC2));
  //       std::vector<std::pair<int,RCP<const FieldPattern> > > patternM;
  //       int numP = 8;
  //       int numU = 4;
  //       patternM.push_back(std::make_pair(numP,patternP));
  //       patternM.push_back(std::make_pair(numU,patternU));

  //       FieldAggPattern agg(patternM);
  //       agg.print(out);

  //       out << "Testing throws of localOffsets_closure" << std::endl;
  //       TEST_THROW(agg.localOffsets_closure(7,0,1),std::logic_error); // check for failure if ID does not exist
  //       TEST_THROW(agg.localOffsets_closure(3,2,6),std::logic_error); // check for failure if sub cell doesn't exist

  //       TEST_EQUALITY(agg.localOffsets_closure(numP,1,6).first.size(),2); // pressure basis only has nodes
  //       TEST_EQUALITY(agg.localOffsets_closure(numP,2,1).first.size(),4); // pressure basis only has nodes

  //       std::vector<int> v_true;
  //       std::vector<int> v;

  //       // test numP first
  //       ///////////////////////////////////////

  //       // test (2,2) closure
  //       out << "Testing localOffsets_closure(numP,2,2) on a C2 HEX" << std::endl;
  //       v_true.resize(4);
  //       v = agg.localOffsets_closure(numP,2,2).first;

  //       // nodes
  //       v_true[0] = 2*2; v_true[1] = 2*3; v_true[2] = 2*7; v_true[3] = 2*6;
 
  //       TEST_EQUALITY(v.size(),v_true.size());

  //       std::sort(v.begin(),v.end());
  //       std::sort(v_true.begin(),v_true.end());
  //       for(std::size_t i=0;i<v.size();i++)
  //          TEST_EQUALITY(v[i],v_true[i]);

  //       // test (1,11) closure
  //       out << "Testing localOffsets_closure(numP,1,11) on a C2 HEX" << std::endl;
  //       v_true.resize(2);
  //       v = agg.localOffsets_closure(numP,1,11).first;

  //       // nodes
  //       v_true[0] = 2*3; v_true[1] = 2*7;
 
  //       TEST_EQUALITY(v.size(),v_true.size());

  //       std::sort(v.begin(),v.end());
  //       std::sort(v_true.begin(),v_true.end());
  //       for(std::size_t i=0;i<v.size();i++)
  //          TEST_EQUALITY(v[i],v_true[i]);

  //       // test all dimension zero closures
  //       out << "Testing localOffsets_closure(numP,0,*) on a C2 HEX" << std::endl;
  //       v_true.resize(1);
  //       for(int i=0;i<8;i++) {
  //          v = agg.localOffsets_closure(numP,0,i).first;

  //          // nodes
  //          v_true[0] = 2*i;
 
  //          TEST_EQUALITY(v.size(),v_true.size());
  //          TEST_EQUALITY(v[0],v_true[0]);
  //       }

  //       // test numU second
  //       ///////////////////////////////////////

  //       // test (2,4) closure
  //       out << "Testing localOffsets_closure(numU,2,4) on a C2 HEX" << std::endl;
  //       v_true.resize(9);
  //       v = agg.localOffsets_closure(numU,2,4).first;

  //       // nodes
  //       v_true[0] = 2*0+1; v_true[1] = 2*1+1; v_true[2] = 2*2+1; v_true[3] = 2*3+1;
 
  //       // edges     
  //       v_true[4]  = 16; v_true[5]  = 17; v_true[6] = 18; v_true[7] = 19;

  //       // areas
  //       v_true[8] = 32;

  //       TEST_EQUALITY(v.size(),v_true.size());

  //       std::sort(v.begin(),v.end());
  //       std::sort(v_true.begin(),v_true.end());
  //       for(std::size_t i=0;i<v.size();i++)
  //          TEST_EQUALITY(v[i],v_true[i]);

  //       // test (1,10) closure
  //       out << "Testing localOffsets_closure(numU,1,10) on a C2 HEX" << std::endl;
  //       v_true.resize(3);
  //       v = agg.localOffsets_closure(numU,1,10).first;

  //       // nodes
  //       v_true[0] = 2*2+1; v_true[1] = 2*6+1;
 
  //       // edges     
  //       v_true[2]  = 26; 

  //       TEST_EQUALITY(v.size(),v_true.size());

  //       std::sort(v.begin(),v.end());
  //       std::sort(v_true.begin(),v_true.end());
  //       for(std::size_t i=0;i<v.size();i++)
  //          TEST_EQUALITY(v[i],v_true[i]);

  //       // test all dimension zero closures
  //       out << "Testing localOffsets_closure(numU,0,*) on a C2 HEX" << std::endl;
  //       v_true.resize(1);
  //       for(int i=0;i<8;i++) {
  //          v = agg.localOffsets_closure(numU,0,i).first;

  //          // nodes
  //          v_true[0] = 2*i+1;
 
  //          TEST_EQUALITY(v.size(),v_true.size());
  //          TEST_EQUALITY(v[0],v_true[0]);
  //       }
  //    }
  // }

  // // Tests the HCURL case where the geometry includes the nodes, but the pattern does not
  // TEUCHOS_UNIT_TEST(tFieldAggPattern, testE)
  // {
  //    out << note << std::endl;

  //    // basis to build patterns from
  //    RCP<Intrepid2::Basis<double,FieldContainer> > basisA = rcp(new Intrepid2::Basis_HCURL_QUAD_I1_FEM<double,FieldContainer>);

  //    RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
  //    RCP<const FieldPattern> patternNode = rcp(new NodalFieldPattern(basisA->getBaseCellTopology()));

  //    std::vector<int> closureIndices;
  //    std::vector<RCP<const FieldPattern> > patternV;
  //    std::vector<std::pair<int,RCP<const FieldPattern> > > patternM;

  //    patternV.push_back(patternA);
  //    patternV.push_back(patternNode);

  //    GeometricAggFieldPattern geom(patternV);

  //    patternM.push_back(std::make_pair(7,patternA));

  //    // test build of geometric field pattern
  //    {
  //       FieldAggPattern agg; 
   
  //       agg.buildPattern(patternM,Teuchos::rcpFromRef(geom));

  //       bool equality = false;
  //       TEST_NOTHROW(equality = geom.equals(*agg.getGeometricAggFieldPattern()));
  //       TEST_ASSERT(equality);

  //       TEST_EQUALITY(geom.getDimension(),agg.getGeometricAggFieldPattern()->getDimension());

  //       out << agg << std::endl;
  //    }
  // }

}
