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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Panzer_DOFManagerFEI.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_Traits.hpp"

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY);

typedef Intrepid::FieldContainer<double> FieldContainer;

const double tolerance = 1.0e-14;

inline double dog_func(double x, double y) { return 4.0*x*x+y; }
inline double cat_func(double x, double y) { return y*y+3.0*x+5.0; }

TEUCHOS_UNIT_TEST(tSolutionReader, test)
{
   RCP<Epetra_Comm> comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

   RCP<panzer::ConnManager<int,int> > connManager 
         = rcp(new panzer_stk::STKConnManager(mesh));
   RCP<panzer::DOFManagerFEI<int,int> > dofManager 
         = rcp(new panzer::DOFManagerFEI<int,int>(connManager,MPI_COMM_WORLD));

   RCP<const panzer::FieldPattern> linPattern = rcp(new panzer::IntrepidFieldPattern(
        rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer>)));
   dofManager->addField("eblock-0_0","cat",linPattern);
   dofManager->addField("dog",linPattern);
   dofManager->buildGlobalUnknowns();

   panzer::EpetraLinearObjFactory<panzer::Traits,int> elof(comm.getConst(),dofManager);
   RCP<Epetra_Vector> vec = rcp(new Epetra_Vector(*elof.getGhostedMap()));

   read_solution_data(*dofManager,*mesh,*vec);

   // grab coordinates and local indices
   std::vector<std::size_t> e0indices, e1indices;
   FieldContainer e0_coords, e1_coords; 
   panzer_stk::workset_utils::getIdsAndVertices(*mesh,"eblock-0_0",e0indices,e0_coords);
   panzer_stk::workset_utils::getIdsAndVertices(*mesh,"eblock-1_0",e1indices,e1_coords);

   const std::vector<int> & dog0_offsets =
         dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("dog"));
   const std::vector<int> & dog1_offsets =
         dofManager->getGIDFieldOffsets("eblock-1_0",dofManager->getFieldNum("dog"));
   const std::vector<int> & cat0_offsets =
         dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("cat"));

   // check the output in the 0,0 block
   for(std::size_t cell=0;cell<e0indices.size();cell++) {
      std::vector<int> GIDs;
      std::size_t localId = e0indices[cell]; 
      dofManager->getElementGIDs(localId,GIDs);

      for(std::size_t b=0;b<4;b++) { // for p/w bilinear quads
         double x = e0_coords(cell,b,0);
         double y = e0_coords(cell,b,1);

         int dog_lid = vec->Map().LID(GIDs[dog0_offsets[b]]);
         int cat_lid = vec->Map().LID(GIDs[cat0_offsets[b]]);

         TEST_FLOATING_EQUALITY((*vec)[dog_lid],dog_func(x,y),tolerance); 
         TEST_FLOATING_EQUALITY((*vec)[cat_lid],cat_func(x,y),tolerance); 
      }
   }

   // check the output in the 1,0 block
   for(std::size_t cell=0;cell<e1indices.size();cell++) {
      std::vector<int> GIDs;
      std::size_t localId = e1indices[cell]; 
      dofManager->getElementGIDs(localId,GIDs);

      for(std::size_t b=0;b<4;b++) { // for p/w bilinear quads
         double x = e1_coords(cell,b,0);
         double y = e1_coords(cell,b,1);

         int dog_lid = vec->Map().LID(GIDs[dog1_offsets[b]]);

         TEST_FLOATING_EQUALITY((*vec)[dog_lid],dog_func(x,y),tolerance); 
      }
   }
}

Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY)
{
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;
  typedef panzer_stk::STK_Interface::VectorFieldType CoordinateField;

  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  pl->set("X Blocks",2);
  pl->set("Y Blocks",1);
  pl->set("X Elements",elemX);
  pl->set("Y Elements",elemY);
  
  panzer_stk::SquareQuadMeshFactory factory;
  factory.setParameterList(pl);
  RCP<panzer_stk::STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
 
  // add in some fields
  mesh->addSolutionField("dog","eblock-0_0");
  mesh->addSolutionField("dog","eblock-1_0");

  mesh->addSolutionField("cat","eblock-0_0");

  factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD); 

  VariableField * dog_field = mesh->getMetaData()->get_field<VariableField>("dog");
  VariableField * cat_field = mesh->getMetaData()->get_field<VariableField>("cat");
  CoordinateField * cField = mesh->getMetaData()->get_field<CoordinateField>("coordinates");
  TEUCHOS_ASSERT(dog_field!=0);
  TEUCHOS_ASSERT(cat_field!=0);
  TEUCHOS_ASSERT(cField!=0);

  // fill the fields with data
  stk::mesh::Part * block0_part = mesh->getElementBlockPart("eblock-0_0");
  TEUCHOS_ASSERT(block0_part!=0);
  const std::vector<stk::mesh::Bucket*> nodeData 
      = mesh->getBulkData()->buckets(mesh->getNodeRank());
  for(std::size_t b=0;b<nodeData.size();++b) {
     stk::mesh::Bucket * bucket = nodeData[b];


     // build all buckets
     for(stk::mesh::Bucket::iterator itr=bucket->begin();
         itr!=bucket->end();++itr) {

        stk::mesh::EntityArray<CoordinateField> coordinates(*cField,*itr);
        stk::mesh::EntityArray<VariableField> dog_array(*dog_field,*itr);

        double x = coordinates(0);
        double y = coordinates(1);

        dog_array() = dog_func(x,y);
        if(bucket->member(*block0_part)) {
           stk::mesh::EntityArray<VariableField> cat_array(*cat_field,*itr);
           cat_array() = cat_func(x,y);
        }
        
     }
  }
  
  if(mesh->isWritable())
     mesh->writeToExodus("SolutionReader.exo");

  return mesh;
}
