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

#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"

// 3D basis 
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"

#include "Shards_BasicTopologies.hpp"

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

// triangle tests
TEUCHOS_UNIT_TEST(tIntrepid2FieldPattern, test2d_tri_c1)
{

   out << note << std::endl;

   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis;

   basis = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);

   Teuchos::RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));
  
   out << "output stream test: " << std::endl;
   out << *pattern << std::endl;

   TEST_ASSERT(pattern->consistentSubcells());
   TEST_EQUALITY(pattern->getDimension(),2);
   TEST_EQUALITY(pattern->numberIds(),3);

   TEST_EQUALITY(pattern->getSubcellCount(0),3); // nodes
   TEST_EQUALITY(pattern->getSubcellCount(1),3); // edges
   TEST_EQUALITY(pattern->getSubcellCount(2),1); // cells
  
   std::vector<int> v;
   v = pattern->getSubcellIndices(0,0);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],0);

   v = pattern->getSubcellIndices(0,1);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],1);

   v = pattern->getSubcellIndices(0,2);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],2);

   TEST_EQUALITY(pattern->getSubcellIndices(1,0).size(),0);
   TEST_EQUALITY(pattern->getSubcellIndices(1,1).size(),0);
   TEST_EQUALITY(pattern->getSubcellIndices(1,2).size(),0);
   TEST_EQUALITY(pattern->getSubcellIndices(2,0).size(),0);

   // test subcell closure functionality: useful for dirchlet boundary conditions
   std::vector<int> bndryIndices;
   std::vector<int> bndryIndices_true(2);

   //
   pattern->getSubcellClosureIndices(1,0,bndryIndices); 
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 1;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,1,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 1;
   bndryIndices_true[1] = 2;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,2,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 2;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));
}


   // test quadratic triangle
TEUCHOS_UNIT_TEST(tIntrepid2FieldPattern, test2d_tri_c2)
{

   out << note << std::endl;

   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis;

   basis = rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM<PHX::exec_space,double,double>);

   Teuchos::RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

   TEST_ASSERT(pattern->consistentSubcells());
   TEST_EQUALITY(pattern->getDimension(),2);
   TEST_EQUALITY(pattern->numberIds(),6);

   TEST_EQUALITY(pattern->getSubcellCount(0),3); // nodes
   TEST_EQUALITY(pattern->getSubcellCount(1),3); // edges
   TEST_EQUALITY(pattern->getSubcellCount(2),1); // cells
  
   std::vector<int> v;
   v = pattern->getSubcellIndices(0,0);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],0);

   v = pattern->getSubcellIndices(0,1);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],1);

   v = pattern->getSubcellIndices(0,2);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],2);

   v = pattern->getSubcellIndices(1,0);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],3);

   v = pattern->getSubcellIndices(1,1);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],4);

   v = pattern->getSubcellIndices(1,2);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],5);

   TEST_EQUALITY(pattern->getSubcellIndices(2,0).size(),0);

   // test subcell closure functionality: useful for dirchlet boundary conditions
   std::vector<int> bndryIndices;
   std::vector<int> bndryIndices_true(3);

   //
   pattern->getSubcellClosureIndices(1,0,bndryIndices); 
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 1;
   bndryIndices_true[2] = 3;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,1,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 1;
   bndryIndices_true[1] = 2;
   bndryIndices_true[2] = 4;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,2,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 2;
   bndryIndices_true[2] = 5;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));
}

// quad tests
TEUCHOS_UNIT_TEST(tIntrepid2FieldPattern, test2d_quad_c1)
{

   out << note << std::endl;

   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis;
   basis = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double>);

   Teuchos::RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

   TEST_ASSERT(pattern->consistentSubcells());
   TEST_EQUALITY(pattern->getDimension(),2);
   TEST_EQUALITY(pattern->numberIds(),4);

   TEST_EQUALITY(pattern->getSubcellCount(0),4); // nodes
   TEST_EQUALITY(pattern->getSubcellCount(1),4); // edges
   TEST_EQUALITY(pattern->getSubcellCount(2),1); // cells
  
   std::vector<int> v;
   v = pattern->getSubcellIndices(0,0);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],0);

   v = pattern->getSubcellIndices(0,1);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],1);

   v = pattern->getSubcellIndices(0,2);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],2);

   v = pattern->getSubcellIndices(0,3);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],3);

   TEST_EQUALITY(pattern->getSubcellIndices(1,0).size(),0);
   TEST_EQUALITY(pattern->getSubcellIndices(1,1).size(),0);
   TEST_EQUALITY(pattern->getSubcellIndices(1,2).size(),0);
   TEST_EQUALITY(pattern->getSubcellIndices(1,3).size(),0);
   TEST_EQUALITY(pattern->getSubcellIndices(2,0).size(),0);

   // test subcell closure functionality: useful for dirchlet boundary conditions
   std::vector<int> bndryIndices;
   std::vector<int> bndryIndices_true(2);

   //
   pattern->getSubcellClosureIndices(1,0,bndryIndices); 
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 1;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,1,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 1;
   bndryIndices_true[1] = 2;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,2,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 2;
   bndryIndices_true[1] = 3;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,3,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 3;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));
}

TEUCHOS_UNIT_TEST(tIntrepid2FieldPattern, test2d_quad_c2)
{

   out << note << std::endl;

   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis;
   basis = rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double>);

   Teuchos::RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

   TEST_ASSERT(pattern->consistentSubcells());
   TEST_EQUALITY(pattern->getDimension(),2);
   TEST_EQUALITY(pattern->numberIds(),9);

   TEST_EQUALITY(pattern->getSubcellCount(0),4); // nodes
   TEST_EQUALITY(pattern->getSubcellCount(1),4); // edges
   TEST_EQUALITY(pattern->getSubcellCount(2),1); // cells

   std::vector<int> v;
   v = pattern->getSubcellIndices(0,0);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],0);

   v = pattern->getSubcellIndices(0,1);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],1);

   v = pattern->getSubcellIndices(0,2);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],2);

   v = pattern->getSubcellIndices(0,3);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],3);

   v = pattern->getSubcellIndices(1,0);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],4);

   v = pattern->getSubcellIndices(1,1);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],5);

   v = pattern->getSubcellIndices(1,2);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],6);

   v = pattern->getSubcellIndices(1,3);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],7);

   v = pattern->getSubcellIndices(2,0);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],8);

   // test subcell closure functionality: useful for dirchlet boundary conditions
   std::vector<int> bndryIndices;
   std::vector<int> bndryIndices_true(3);

   //
   pattern->getSubcellClosureIndices(1,0,bndryIndices); 
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 1;
   bndryIndices_true[2] = 4;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,1,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 1;
   bndryIndices_true[1] = 2;
   bndryIndices_true[2] = 5;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,2,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 2;
   bndryIndices_true[1] = 3;
   bndryIndices_true[2] = 6;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));
  
   //
   pattern->getSubcellClosureIndices(1,3,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 3;
   bndryIndices_true[2] = 7;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   // test Cell
   pattern->getSubcellClosureIndices(2,0,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   for(std::size_t i=0;i<bndryIndices_true.size();++i)
      bndryIndices_true[i] = i;
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));
}

/////////////////////////////////////////////
// 3D tests
/////////////////////////////////////////////

// hex tests
TEUCHOS_UNIT_TEST(tIntrepid2FieldPattern, test3d_hex_c1)
{

   out << note << std::endl;

   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis;
   basis = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);

   Teuchos::RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

   TEST_ASSERT(pattern->consistentSubcells());
   TEST_EQUALITY(pattern->getDimension(),3);
   TEST_EQUALITY(pattern->numberIds(),8);

   TEST_EQUALITY(pattern->getSubcellCount(0),8); // nodes
   TEST_EQUALITY(pattern->getSubcellCount(1),12); // edges
   TEST_EQUALITY(pattern->getSubcellCount(2),6); // faces
   TEST_EQUALITY(pattern->getSubcellCount(3),1); // cells

   std::vector<int> v;
   for(int i=0;i<8;i++) {
      v = pattern->getSubcellIndices(0,i);
      TEST_EQUALITY(v.size(),1);
      TEST_EQUALITY(v[0],i);
   }

   // test subcell closure functionality: useful for dirchlet boundary conditions
   std::vector<int> bndryIndices;
   std::vector<int> bndryIndices_true(2);

   //
   pattern->getSubcellClosureIndices(1,1,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 1;
   bndryIndices_true[1] = 2;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,2,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 2;
   bndryIndices_true[1] = 3;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));
  
   //
   pattern->getSubcellClosureIndices(1,6,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 6;
   bndryIndices_true[1] = 7;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   // test face
   bndryIndices_true.resize(4);

   pattern->getSubcellClosureIndices(2,3,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 3;
   bndryIndices_true[2] = 4;
   bndryIndices_true[3] = 7;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   pattern->getSubcellClosureIndices(2,1,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 1;
   bndryIndices_true[1] = 2;
   bndryIndices_true[2] = 5;
   bndryIndices_true[3] = 6;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

//    // test cell
//    bndryIndices_true.resize(8);
// 
//    pattern->getSubcellClosureIndices(3,0,bndryIndices);
//    std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision
// 
//    for(std::size_t i=0;i<bndryIndices_true.size();++i)
//       bndryIndices_true[i] = i;
//    TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));
}

TEUCHOS_UNIT_TEST(tIntrepid2FieldPattern, test3d_hex_c2)
{

   out << note << std::endl;

   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis;
   basis = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);

   Teuchos::RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

   TEST_ASSERT(pattern->consistentSubcells());
   TEST_EQUALITY(pattern->getDimension(),3);
   TEST_EQUALITY(pattern->numberIds(),27);

   TEST_EQUALITY(pattern->getSubcellCount(0),8); // nodes
   TEST_EQUALITY(pattern->getSubcellCount(1),12); // edges
   TEST_EQUALITY(pattern->getSubcellCount(2),6); // faces
   TEST_EQUALITY(pattern->getSubcellCount(3),1); // cells

   std::vector<int> v;
   for(int i=0;i<8;i++) {
      v = pattern->getSubcellIndices(0,i);
      TEST_EQUALITY(v.size(),1);
      TEST_EQUALITY(v[0],i);
   }

   for(int i=0;i<4;i++) {
      v = pattern->getSubcellIndices(1,i);
      TEST_EQUALITY(v.size(),1);
      TEST_EQUALITY(v[0],i+8);
   }

   for(int i=0;i<4;i++) {
      v = pattern->getSubcellIndices(1,4+i);
      TEST_EQUALITY(v.size(),1);
      TEST_EQUALITY(v[0],i+16);
   }

   for(int i=0;i<4;i++) {
      v = pattern->getSubcellIndices(1,8+i);
      TEST_EQUALITY(v.size(),1);
      TEST_EQUALITY(v[0],i+12);
   }

   std::vector<int> v_true(6);
   v_true[0] = 25;
   v_true[1] = 24;
   v_true[2] = 26;
   v_true[3] = 23;
   v_true[4] = 21;
   v_true[5] = 22;
   for(std::size_t i=0;i<6;i++) {
      v = pattern->getSubcellIndices(2,i);
      TEST_EQUALITY(v.size(),1);
      TEST_EQUALITY(v[0],v_true[i]);
   }

   v = pattern->getSubcellIndices(3,0);
   TEST_EQUALITY(v.size(),1);
   TEST_EQUALITY(v[0],20);


   // test subcell closure functionality: useful for dirchlet boundary conditions
   out << "Begin Dim 1 subcell closures" << std::endl;
   std::vector<int> bndryIndices;
   std::vector<int> bndryIndices_true(3);

   //
   pattern->getSubcellClosureIndices(1,1,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 1;
   bndryIndices_true[1] = 2;
   bndryIndices_true[2] = 9;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   //
   pattern->getSubcellClosureIndices(1,2,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 2;
   bndryIndices_true[1] = 3;
   bndryIndices_true[2] = 10;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));
  
   //
   pattern->getSubcellClosureIndices(1,6,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 6;
   bndryIndices_true[1] = 7;
   bndryIndices_true[2] = 18;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));

   // test face
   out << "Begin Dim 2 subcell closures" << std::endl;
   bndryIndices_true.resize(9);

   pattern->getSubcellClosureIndices(2,3,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 3;
   bndryIndices_true[2] = 4;
   bndryIndices_true[3] = 7;
   bndryIndices_true[4] = 11;
   bndryIndices_true[5] = 12;
   bndryIndices_true[6] = 15;
   bndryIndices_true[7] = 19;
   bndryIndices_true[8] = 23;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices.begin(),bndryIndices.end(),bndryIndices_true.begin()));

   pattern->getSubcellClosureIndices(2,1,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 1;
   bndryIndices_true[1] = 2;
   bndryIndices_true[2] = 5;
   bndryIndices_true[3] = 6;
   bndryIndices_true[4] = 9;
   bndryIndices_true[5] = 13;
   bndryIndices_true[6] = 14;
   bndryIndices_true[7] = 17;
   bndryIndices_true[8] = 24;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices.begin(),bndryIndices.end(),bndryIndices_true.begin()));

   out << "EXPECTING NEW SUB CELL CLOSURE!" << std::endl;
   bndryIndices.clear();
   pattern->getSubcellClosureIndices(2,0,bndryIndices);
   std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision

   bndryIndices_true[0] = 0;
   bndryIndices_true[1] = 1;
   bndryIndices_true[2] = 4;
   bndryIndices_true[3] = 5;
   bndryIndices_true[4] = 8;
   bndryIndices_true[5] = 12;
   bndryIndices_true[6] = 13;
   bndryIndices_true[7] = 16;
   bndryIndices_true[8] = 25;
   TEST_ASSERT(bndryIndices.size()==bndryIndices_true.size());
   TEST_ASSERT(std::equal(bndryIndices.begin(),bndryIndices.end(),bndryIndices_true.begin()));

   out << "BEGIN COMPARISON" << std::endl;
   for(std::size_t i=0;i<bndryIndices_true.size();++i)
      out << bndryIndices[i] << " ?==" << bndryIndices_true[i] << std::endl;
//    // test cell
//    bndryIndices_true.resize(8);
// 
//    pattern->getSubcellClosureIndices(3,0,bndryIndices);
//    std::sort(bndryIndices.begin(),bndryIndices.end()); // sort for comparision
// 
//    for(std::size_t i=0;i<bndryIndices_true.size();++i)
//       bndryIndices_true[i] = i;
//    TEST_ASSERT(std::equal(bndryIndices_true.begin(),bndryIndices_true.end(),bndryIndices.begin()));
}

TEUCHOS_UNIT_TEST(tIntrepid2FieldPattern, test2d_tri_hdivi1)
{

   out << note << std::endl;

   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HDIV_TRI_I1_FEM<PHX::exec_space,double,double>);
   RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

   pattern->print(out);

   TEST_EQUALITY(pattern->getDimension(),2);
   TEST_EQUALITY(pattern->numberIds(),3);
}

///////////////////////////////////////////
// constant tests
///////////////////////////////////////////

TEUCHOS_UNIT_TEST(tIntrepid2FieldPattern, constant_basis)
{

   out << note << std::endl;

   shards::CellTopology tet(shards::getCellTopologyData<shards::Tetrahedron<4> >());
   shards::CellTopology hex(shards::getCellTopologyData<shards::Hexahedron<8> >());
   shards::CellTopology tri(shards::getCellTopologyData<shards::Triangle<3> >());
   shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >()); // i'm irritated there are 4 letters in quad!

   {
     RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HVOL_C0_FEM<PHX::exec_space,double,double>(tet));
     RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

     pattern->print(out);

     TEST_EQUALITY(pattern->getDimension(),3);
     TEST_EQUALITY(pattern->numberIds(),1);
     TEST_EQUALITY(pattern->getSubcellCount(0),4);
     TEST_EQUALITY(pattern->getSubcellCount(1),6);
     TEST_EQUALITY(pattern->getSubcellCount(2),4);
     TEST_EQUALITY(pattern->getSubcellCount(3),1);

     // make sure unknowns are distributed correctly
     TEST_EQUALITY(pattern->getSubcellIndices(0,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,2).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,3).size(),0);

     TEST_EQUALITY(pattern->getSubcellIndices(1,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,2).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,3).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,4).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,5).size(),0);

     TEST_EQUALITY(pattern->getSubcellIndices(2,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(2,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(2,2).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(2,3).size(),0);

     const std::vector<int> & vec = pattern->getSubcellIndices(3,0);
     TEST_EQUALITY(vec.size(),1);
     TEST_EQUALITY(vec[0],0);
   }

   {
     RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HVOL_C0_FEM<PHX::exec_space,double,double>(hex));
     RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

     pattern->print(out);

     TEST_EQUALITY(pattern->getDimension(),3);
     TEST_EQUALITY(pattern->numberIds(),1);
     TEST_EQUALITY(pattern->getSubcellCount(0),8);
     TEST_EQUALITY(pattern->getSubcellCount(1),12);
     TEST_EQUALITY(pattern->getSubcellCount(2),6);
     TEST_EQUALITY(pattern->getSubcellCount(3),1);

     // make sure unknowns are distributed correctly
     TEST_EQUALITY(pattern->getSubcellIndices(0,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,2).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,3).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,4).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,5).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,6).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,7).size(),0);

     TEST_EQUALITY(pattern->getSubcellIndices(1,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,2).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,3).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,4).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,5).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,6).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,7).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,8).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,9).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,10).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,11).size(),0);

     TEST_EQUALITY(pattern->getSubcellIndices(2,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(2,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(2,2).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(2,3).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(2,4).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(2,5).size(),0);

     const std::vector<int> & vec = pattern->getSubcellIndices(3,0);
     TEST_EQUALITY(vec.size(),1);
     TEST_EQUALITY(vec[0],0);
   }

   {
     RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HVOL_C0_FEM<PHX::exec_space,double,double>(tri));
     RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

     pattern->print(out);

     TEST_EQUALITY(pattern->getDimension(),2);
     TEST_EQUALITY(pattern->numberIds(),1);
     TEST_EQUALITY(pattern->getSubcellCount(0),3);
     TEST_EQUALITY(pattern->getSubcellCount(1),3);
     TEST_EQUALITY(pattern->getSubcellCount(2),1);

     // make sure unknowns are distributed correctly
     TEST_EQUALITY(pattern->getSubcellIndices(0,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,2).size(),0);

     TEST_EQUALITY(pattern->getSubcellIndices(1,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,2).size(),0);

     const std::vector<int> & vec = pattern->getSubcellIndices(2,0);
     TEST_EQUALITY(vec.size(),1);
     TEST_EQUALITY(vec[0],0);
   }

   {
     RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HVOL_C0_FEM<PHX::exec_space,double,double>(quad));
     RCP<FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

     pattern->print(out);

     TEST_EQUALITY(pattern->getDimension(),2);
     TEST_EQUALITY(pattern->numberIds(),1);
     TEST_EQUALITY(pattern->getSubcellCount(0),4);
     TEST_EQUALITY(pattern->getSubcellCount(1),4);
     TEST_EQUALITY(pattern->getSubcellCount(2),1);

     // make sure unknowns are distributed correctly
     TEST_EQUALITY(pattern->getSubcellIndices(0,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,2).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(0,3).size(),0);

     TEST_EQUALITY(pattern->getSubcellIndices(1,0).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,1).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,2).size(),0);
     TEST_EQUALITY(pattern->getSubcellIndices(1,3).size(),0);

     const std::vector<int> & vec = pattern->getSubcellIndices(2,0);
     TEST_EQUALITY(vec.size(),1);
     TEST_EQUALITY(vec[0],0);
   }
}

}
