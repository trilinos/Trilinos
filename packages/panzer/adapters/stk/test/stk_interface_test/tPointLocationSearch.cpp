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

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

// STK_search objects
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search_util/stk_mesh/CreateBoundingBox.hpp>

typedef stk::mesh::Field<double, stk::mesh::Cartesian>  VectorField;

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tPointLocationSearch, basic)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",1);
   pl->set("Z Elements",1);
   pl->set("X0",0.0);
   pl->set("Xf",1.0);
   pl->set("Y0",0.0);
   pl->set("Yf",1.0);
   pl->set("Z0",0.0);
   pl->set("Zf",1.0);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   // Build (domain) bounding boxes for all cells in mesh
   RCP<stk::mesh::fem::FEMMetaData> meta_data = mesh->getMetaData();
   RCP<stk::mesh::BulkData> bulk_data = mesh->getBulkData(); 
   const stk::mesh::Field<double, stk::mesh::Cartesian>* domain_coord_field = &(mesh->getCoordinatesField());
   stk::ParallelMachine comm = bulk_data->parallel();
   // NOTE: the create bounding boxes call has specific typedefs on data.  We need to rewrite for general case.
   std::vector<AxisAlignedBoundingBox3D> domain_vector;
   

   Teuchos::FancyOStream os(Teuchos::rcpFromRef(std::cout));
   os.setShowProcRank(true);
   os.setProcRankAndSize(stk::parallel_machine_rank(comm),stk::parallel_machine_size(comm));



     /*

   std::vector<stk::mesh::Entity*> my_elements;
   mesh->getMyElements(my_elements);
   for (std::vector<stk::mesh::Entity*>::const_iterator e=my_elements.begin(); e!=my_elements.end();++e) {
     os << "element id = " << (*e)->identifier() << std::endl;
   

     AxisAlignedBoundingBox3D bbox;
     std::vector<double> box(6);
     // loop over nodes to get min/max coordinates and center

     // set on bbox 

     // set the key for bbox

     // add bbox to vector
     
   }
   
     */


   stk::search_util::build_axis_aligned_bbox(*bulk_data,
					     mesh->getElementRank(),
					     const_cast<stk::mesh::Field<double, stk::mesh::Cartesian>* >(domain_coord_field),
					     domain_vector);
   






   os << "\nsize of domain = " << domain_vector.size() << std::endl;
   for (std::vector<AxisAlignedBoundingBox3D>::const_iterator i=domain_vector.begin(); i != domain_vector.end(); ++i)
     os << i->key.ident.id() << "," << i->key.proc << " = " 
	<< "X(" << i->box[0] << "," << i->box[0+AxisAlignedBoundingBox3D::DIMENSION] << ") "
	<< "Y(" << i->box[1] << "," << i->box[1+AxisAlignedBoundingBox3D::DIMENSION] << ") "
	<< "Z(" << i->box[2] << "," << i->box[2+AxisAlignedBoundingBox3D::DIMENSION] << ") "
	<< std::endl;
   

   std::vector<PointBoundingBox3D> pts;
   {   
     PointBoundingBox3D p;
     double center[3] = {0.25,0.25,0.25};
     
     stk::mesh::EntityKey pt_key(0,451);
     IdentProc id(pt_key,stk::parallel_machine_rank(comm));
     p.key = id;
     p.set_center(center);
     
     pts.push_back(p);
   }

   stk::search::FactoryOrder order;
   order.m_communicator = comm;
   order.m_algorithm = stk::search::FactoryOrder::BIHTREE;
   
   IdentProcRelation relation;
   
   stk::search::coarse_search(relation, pts, domain_vector, order);
   
   for (IdentProcRelation::const_iterator i=relation.begin(); i != relation.end(); ++i)
     os << "Relation Domain(" << i->first.ident.id() << "," << i->first.proc << ") " 
	<< "Range(" << i->second.ident.id() << "," << i->second.proc << ")"
	<< std::endl;
   
   

}

}
