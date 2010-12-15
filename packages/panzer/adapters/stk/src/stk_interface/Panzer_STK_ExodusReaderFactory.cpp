
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_Interface.hpp"

#ifdef HAVE_IOSS 

#include <stk_io/util/UseCase_mesh.hpp>
#include <Ionit_Initializer.h>

#include <stk_io/IossBridge.hpp>

namespace panzer_stk {

STK_ExodusReaderFactory::STK_ExodusReaderFactory()
   : fileName_("")
{ }

STK_ExodusReaderFactory::STK_ExodusReaderFactory(const std::string & fileName)
   : fileName_(fileName)
{ }

Teuchos::RCP<STK_Interface> STK_ExodusReaderFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;

   RCP<STK_Interface> mesh = rcp(new STK_Interface(2));
   RCP<stk::mesh::MetaData> metaData = mesh->getMetaData();

   // read in meta data
   Ioss::Init::Initializer io;
   stk::io::util::MeshData meshData;
   stk::io::util::create_input_mesh("exodusii", fileName_, "", parallelMach,
                                    *metaData, meshData, false); 
   // build spactial dimension and some topological meta data
   unsigned spaceDim = metaData->get_field<VectorFieldType>("coordinates")->max_size(stk::mesh::Node);
   // mesh->setDimension(spaceDim);

   // build element blocks
   registerElementBlocks(*mesh);
   registerSidesets(*mesh);

   // this calls commit on meta data
   mesh->initialize(parallelMach,false);

   mesh->beginModification();
   RCP<stk::mesh::BulkData> bulkData = mesh->getBulkData();
   stk::io::util::populate_bulk_data(*bulkData, meshData, "exodusii");
   // bulkData->modification_end();
   mesh->endModification();
   
   return mesh;
}

//! From ParameterListAcceptor
void STK_ExodusReaderFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> STK_ExodusReaderFactory::getValidParameters() const
{
   return Teuchos::null;
}

void STK_ExodusReaderFactory::registerElementBlocks(STK_Interface & mesh) const 
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::mesh::PartVector & parts = metaData->get_parts();

   stk::mesh::PartVector::const_iterator partItr;
   for(partItr=parts.begin();partItr!=parts.end();++partItr) {
      const stk::mesh::Part * part = *partItr;
      const CellTopologyData * ct = stk::mesh::fem::get_cell_topology(*part).getCellTopologyData();

      // if an element part ==> this is element block
      if(part->primary_entity_rank()==stk::mesh::Element) {
         TEUCHOS_ASSERT(ct!=0);
         mesh.addElementBlock(part->name(),ct);
      }
   }
}

void STK_ExodusReaderFactory::registerSidesets(STK_Interface & mesh) const
{
   using Teuchos::RCP;
   unsigned spatialDim = mesh.getDimension();

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::mesh::PartVector & parts = metaData->get_parts();

   stk::mesh::PartVector::const_iterator partItr;
   for(partItr=parts.begin();partItr!=parts.end();++partItr) {
      const stk::mesh::Part * part = *partItr;
      const CellTopologyData * ct = stk::mesh::fem::get_cell_topology(*part).getCellTopologyData();

      // if an side part ==> this is a sideset: only take null cell topologies
      if(part->primary_entity_rank()==spatialDim-1 && ct==0)
         mesh.addSideset(part->name());
   }
}

}

#endif
