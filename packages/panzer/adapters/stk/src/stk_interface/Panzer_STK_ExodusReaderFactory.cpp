
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

   RCP<STK_Interface> mesh = rcp(new STK_Interface());
   RCP<stk::mesh::MetaData> metaData = mesh->getMetaData();

   // read in meta data
   Ioss::Init::Initializer io;
   stk::io::util::MeshData meshData;
   stk::io::util::create_input_mesh("exodusii", fileName_, "", parallelMach,
                                    *metaData, meshData, false); 

   // build spactial dimension and some topological meta data
   unsigned spaceDim = metaData->get_field<VectorFieldType>("coordinates")->max_size(stk::mesh::Node);
   stk::mesh::TopologicalMetaData md(*metaData,spaceDim);

   mesh->setDimension(md.spatial_dimension);

   // build element blocks
   registerElementBlocks(*mesh,md);
   registerSidesets(*mesh,md);

   // this calls commit on meta data
   mesh->initialize(parallelMach,false);
   
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

void STK_ExodusReaderFactory::registerElementBlocks(STK_Interface & mesh,const stk::mesh::TopologicalMetaData & md) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::mesh::PartVector & parts = metaData->get_parts();

   stk::mesh::PartVector::const_iterator partItr;
   for(partItr=parts.begin();partItr!=parts.end();++partItr) {
      const stk::mesh::Part * part = *partItr;

      // if an element part ==> this is element block
      if(part->primary_entity_rank()==stk::mesh::Element) {
         const CellTopologyData * ct = stk::mesh::get_cell_topology(*part);
         mesh.addElementBlock(part->name(),ct);
      }
   }
}

void STK_ExodusReaderFactory::registerSidesets(STK_Interface & mesh,const stk::mesh::TopologicalMetaData & md) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::mesh::PartVector & parts = metaData->get_parts();

   stk::mesh::PartVector::const_iterator partItr;
   for(partItr=parts.begin();partItr!=parts.end();++partItr) {
      const stk::mesh::Part * part = *partItr;

      // if an element part ==> this is element block
      if(part->primary_entity_rank()==md.spatial_dimension-1)
         mesh.addSideset(part->name());
   }
}

}

#endif
