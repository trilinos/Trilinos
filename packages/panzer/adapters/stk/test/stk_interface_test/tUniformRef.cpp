#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

/*
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_percept/PerceptMesh.hpp"
#include "stk_adapt/UniformRefiner.hpp"
#include "stk_adapt/UniformRefinerPattern.hpp"
*/
#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/UniformRefiner.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tUniformRef, stk_fixture)
{
   typedef stk::mesh::Field<int> ProcIdFieldType;

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,3,3,3);
   ProcIdFieldType & processorIdField = hf.m_fem_meta.declare_field<ProcIdFieldType>("PROC_ID");
   stk::mesh::put_field( processorIdField , 3, hf.m_fem_meta.universal_part());

   stk::percept::PerceptMesh eMesh(&hf.m_fem_meta,&hf.m_bulk_data,false);
   stk::adapt::Hex8_Hex8_8 break_quad(eMesh);

   hf.m_fem_meta.commit();
   hf.generate_mesh();

   const std::vector<stk::mesh::Bucket*> & buckets = hf.m_bulk_data.buckets(3);
   for(std::size_t i=0;i<buckets.size();++i) {
      stk::mesh::Bucket & b = *buckets[i]; 
      for(std::size_t j=0;j<b.size();++j) {
         stk::mesh::Entity & element = b[j];
         // set processor rank
         int * procId = stk::mesh::field_data(processorIdField,element);
         procId[0] = hf.m_bulk_data.parallel_rank();
      }
   }

   stk::adapt::UniformRefiner breaker(eMesh, break_quad, &processorIdField);
   breaker.setRemoveOldElements(true);

   breaker.doBreak();
}

#if 0
TEUCHOS_UNIT_TEST(tUniformRef, refTest)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",3);
   pl->set("Z Elements",3);

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   // SquareQuadMeshFactory factory; 
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
   mesh->instantiateBulkData(MPI_COMM_WORLD);

   RCP<stk::mesh::fem::FEMMetaData> meta = mesh->getMetaData();
   // meta->FEM_initialize(2);
   RCP<stk::mesh::BulkData> bulk = mesh->getBulkData();

   stk::percept::PerceptMesh eMesh(&*meta,&*bulk,false);
   // stk::adapt::Quad4_Quad4_4 break_quad(eMesh);
   stk::adapt::Hex8_Hex8_8 break_quad(eMesh);

   // mesh->initialize(MPI_COMM_WORLD);
   factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);


   stk::mesh::FieldBase* proc_rank_field = eMesh.getField("PROC_ID");

   stk::adapt::UniformRefiner breaker(eMesh, break_quad, proc_rank_field);
   breaker.setRemoveOldElements(true);

   breaker.doBreak();

/*
   const stk::mesh::FieldBase::Restriction & r = mesh->getCoordinatesField().restriction(stk::mesh::fem::FEMMetaData::NODE_RANK, meta->universal_part());
   unsigned dataStride = r.dimension();
   std::cout << "DS = " << dataStride << std::endl;
 
   factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
   stk::percept::PerceptMesh eMesh(&*meta, &*bulk, true);
   stk::adapt::Quad4_Quad4_4 quad(eMesh);

   int scalarDimension = 0;
   stk::adapt::UniformRefiner breaker(eMesh, quad, proc_rank_field);
   breaker.doBreak();
*/
}
#endif

}
