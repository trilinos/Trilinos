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

#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"

// 3D basis 
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"

#include "Intrepid_HDIV_TRI_I1_FEM.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

/////////////////////////////////////////////
// 2D tests
/////////////////////////////////////////////

// triangle tests
TEUCHOS_UNIT_TEST(tNodalFieldPattern, test2d_tri_c1)
{
   typedef Intrepid::FieldContainer<double> FieldContainer;
   RCP<Intrepid::Basis<double,FieldContainer> > basis1, basis2;

   basis1 = rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double,FieldContainer>);
   basis2 = rcp(new Intrepid::Basis_HGRAD_TRI_C2_FEM<double,FieldContainer>);

   Teuchos::RCP<FieldPattern> pattern1 = rcp(new IntrepidFieldPattern(basis1));
   Teuchos::RCP<FieldPattern> pattern2 = rcp(new IntrepidFieldPattern(basis2));
   Teuchos::RCP<FieldPattern> nodalPattern = rcp(new NodalFieldPattern(pattern1->getCellTopology()));
  
   TEST_ASSERT(nodalPattern->equals(*pattern1))
   TEST_ASSERT(!nodalPattern->equals(*pattern2))
}

TEUCHOS_UNIT_TEST(tNodalFieldPattern, test3d_HEX_c1)
{
   typedef Intrepid::FieldContainer<double> FieldContainer;
   RCP<Intrepid::Basis<double,FieldContainer> > basis1, basis2;

   basis1 = rcp(new Intrepid::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer>);
   basis2 = rcp(new Intrepid::Basis_HGRAD_HEX_C2_FEM<double,FieldContainer>);

   Teuchos::RCP<FieldPattern> pattern1 = rcp(new IntrepidFieldPattern(basis1));
   Teuchos::RCP<FieldPattern> pattern2 = rcp(new IntrepidFieldPattern(basis2));
   Teuchos::RCP<FieldPattern> nodalPattern = rcp(new NodalFieldPattern(pattern1->getCellTopology()));
  
   TEST_ASSERT(nodalPattern->equals(*pattern1))
   TEST_ASSERT(!nodalPattern->equals(*pattern2))
}

}
