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

// include some intrepid basis functions
// 2D basis 
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"

// 3D basis 
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

std::string note = "***   NOTE: UNIT TEST BASED ON SEPT 2010   ***\n"
                   "***   INTREPID AND SHARDS Trilinos-dev     ***\n"
                   "***   DOXYGEN WEBSITE                      ***\n";

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

/////////////////////////////////////////////
// 2D tests
/////////////////////////////////////////////

bool intrepid_equals(const RCP<Intrepid2::Basis<PHX::exec_space,double,double> > & basisA,
                     const RCP<Intrepid2::Basis<PHX::exec_space,double,double> > & basisB,
                     const std::string & file,int lineNo);
bool intrepid_same_geom(const RCP<Intrepid2::Basis<PHX::exec_space,double,double> > & basisA,
                       const RCP<Intrepid2::Basis<PHX::exec_space,double,double> > & basisB,
                       const std::string & file,int lineNo);

// triangle tests
TEUCHOS_UNIT_TEST(tFieldPattern, test_equals)
{

   out << note << std::endl;

   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA;
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB;

   basisA = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
   basisB = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
   TEST_ASSERT(intrepid_same_geom(basisA,basisB,__FILE__,__LINE__));
   TEST_ASSERT(intrepid_equals(basisA,basisB,__FILE__,__LINE__));

   basisA = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
   basisB = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);
   TEST_ASSERT(intrepid_same_geom(basisA,basisB,__FILE__,__LINE__));
   TEST_ASSERT(not intrepid_equals(basisA,basisB,__FILE__,__LINE__));

   basisA = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);
   basisB = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);
   TEST_ASSERT(intrepid_same_geom(basisA,basisB,__FILE__,__LINE__));
   TEST_ASSERT(intrepid_equals(basisA,basisB,__FILE__,__LINE__));

   basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double>);
   basisB = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
   TEST_ASSERT(not intrepid_same_geom(basisA,basisB,__FILE__,__LINE__));
   TEST_ASSERT(not intrepid_equals(basisA,basisB,__FILE__,__LINE__));

   basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double>);
   basisB = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double>);
   TEST_ASSERT(intrepid_same_geom(basisA,basisB,__FILE__,__LINE__));
   TEST_ASSERT(intrepid_equals(basisA,basisB,__FILE__,__LINE__));

   basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double>);
   basisB = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double>);
   TEST_ASSERT(intrepid_same_geom(basisA,basisB,__FILE__,__LINE__));
   TEST_ASSERT(not intrepid_equals(basisA,basisB,__FILE__,__LINE__));

   basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double>);
   basisB = rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double>);
   TEST_ASSERT(intrepid_same_geom(basisA,basisB,__FILE__,__LINE__));
   TEST_ASSERT(intrepid_equals(basisA,basisB,__FILE__,__LINE__));
}

bool intrepid_equals(const RCP<Intrepid2::Basis<PHX::exec_space,double,double> > & basisA,
                     const RCP<Intrepid2::Basis<PHX::exec_space,double,double> > & basisB,
                     const std::string & file,int lineNo)
{
   // notice if Intrepid2FieldPattern implements "equals" then this functionality
   // is not quite correct!
   Teuchos::RCP<const FieldPattern> fpA = rcp(new Intrepid2FieldPattern(basisA));
   Teuchos::RCP<const FieldPattern> fpB= rcp(new Intrepid2FieldPattern(basisB));

   bool forward = fpA->equals(*fpB); 
   bool bckward = fpB->equals(*fpA); 

   if(bckward!=forward) { // if true, this is a big problem
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                         "Test \"intrepid_equals\" leads to an"
                         "inconsistency with the FieldPattern::equals "
                         "functionality: " << file << ": " << lineNo );
   }
   TEUCHOS_ASSERT(bckward==forward);

   return forward; // they are equal so that it shouldn't matter at this point 
}

bool intrepid_same_geom(const RCP<Intrepid2::Basis<PHX::exec_space,double,double> > & basisA,
                       const RCP<Intrepid2::Basis<PHX::exec_space,double,double> > & basisB,
                       const std::string & file,int lineNo)
{
   // notice if Intrepid2FieldPattern implements "sameGeometry" then this functionality
   // is not quite correct!
   Teuchos::RCP<const FieldPattern> fpA = rcp(new Intrepid2FieldPattern(basisA));
   Teuchos::RCP<const FieldPattern> fpB= rcp(new Intrepid2FieldPattern(basisB));

   bool forward = fpA->sameGeometry(*fpB); 
   bool bckward = fpB->sameGeometry(*fpA); 

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
