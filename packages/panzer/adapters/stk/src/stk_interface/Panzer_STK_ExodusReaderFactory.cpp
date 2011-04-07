
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_Interface.hpp"

#ifdef HAVE_IOSS 

#include <stk_io/util/UseCase_mesh.hpp>
#include <Ionit_Initializer.h>
#include <Ioss_ElementBlock.h>

#include <stk_io/IossBridge.hpp>

#include "Teuchos_StandardParameterEntryValidators.hpp"

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
   
   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // in here you would add your fields...but it is better to use
   // the two step construction

   // this calls commit on meta data
   mesh->initialize(parallelMach,false);

   completeMeshConstruction(*mesh,parallelMach); 

   return mesh;
}

/** This builds all the meta data of the mesh. Does not call metaData->commit.
  * Allows user to add solution fields and other pieces. The mesh can be "completed"
  * by calling <code>completeMeshConstruction</code>.
  */
Teuchos::RCP<STK_Interface> STK_ExodusReaderFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const 
{ 
   using Teuchos::RCP;
   using Teuchos::rcp;

   RCP<STK_Interface> mesh = rcp(new STK_Interface());
   RCP<stk::mesh::fem::FEMMetaData> femMetaData = mesh->getMetaData();
   stk::mesh::MetaData & metaData = stk::mesh::fem::FEMMetaData::get_meta_data(*femMetaData);

   // read in meta data
   Ioss::Init::Initializer io;
   stk::io::util::MeshData * meshData = new stk::io::util::MeshData;
   stk::io::util::create_input_mesh("exodusii", fileName_, "", parallelMach,
                                    *femMetaData, *meshData, false); 

   // store mesh data pointer for later use in initializing 
   // bulk data
   metaData.declare_attribute_with_delete(meshData);

   mesh->initializeFromMetaData();

   // build element blocks
   registerElementBlocks(*mesh,*meshData);
   registerSidesets(*mesh,*meshData);

   mesh->addPeriodicBCs(periodicBCVec_);

   return mesh; 
}

void STK_ExodusReaderFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   // grab mesh data pointer to build the bulk data
   stk::mesh::MetaData & metaData = stk::mesh::fem::FEMMetaData::get_meta_data(*mesh.getMetaData());
   stk::io::util::MeshData * meshData = 
         const_cast<stk::io::util::MeshData *>(metaData.get_attribute<stk::io::util::MeshData>());
         // if const_cast is wrong ... why does it feel so right?
         // I believe this is safe since we are basically hiding this object under the covers
         // until the mesh construction can be completed...below I cleanup the object myself.
   TEUCHOS_ASSERT(metaData.remove_attribute(meshData)); 
      // remove the MeshData attribute

   RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   mesh.beginModification();
   stk::io::util::populate_bulk_data(*bulkData, *meshData, "exodusii");
   mesh.endModification();

   mesh.buildSubcells();
   mesh.buildLocalElementIDs();

   // clean up mesh data object
   delete meshData;
}

//! From ParameterListAcceptor
void STK_ExodusReaderFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   TEST_FOR_EXCEPTION_PURE_MSG(!paramList->isParameter("File Name"),
        Teuchos::Exceptions::InvalidParameterName,
        "Error, the parameter {name=\"File Name\","
        "type=\"string\""
        "\nis required in parameter (sub)list \""<< paramList->name() <<"\"."
        "\n\nThe parsed parameter parameter list is: \n" << paramList->currentParametersString()
   );
      
   paramList->validateParametersAndSetDefaults(*getValidParameters(),0); 

   setMyParamList(paramList);

   fileName_ = paramList->get<std::string>("File Name");

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> STK_ExodusReaderFactory::getValidParameters() const
{
   static Teuchos::RCP<Teuchos::ParameterList> validParams;

   if(validParams==Teuchos::null) {
      validParams = Teuchos::rcp(new Teuchos::ParameterList);
      validParams->set<std::string>("File Name","<file name not set>","Name of exodus file to be read", 
                                    Teuchos::rcp(new Teuchos::FileNameValidator));

      Teuchos::ParameterList & bcs = validParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
   }

   return validParams.getConst();
}

void STK_ExodusReaderFactory::registerElementBlocks(STK_Interface & mesh,stk::io::util::MeshData & meshData) const 
{
   using Teuchos::RCP;

   RCP<stk::mesh::fem::FEMMetaData> femMetaData = mesh.getMetaData();
   stk::mesh::MetaData & metaData = stk::mesh::fem::FEMMetaData::get_meta_data(*femMetaData);

   // here we use the Ioss interface because they don't add
   // "bonus" element blocks and its easier to determine
   // "real" element blocks versus STK-only blocks
   const Ioss::ElementBlockContainer & elem_blocks = meshData.m_region->get_element_blocks();
   for(Ioss::ElementBlockContainer::const_iterator itr=elem_blocks.begin();itr!=elem_blocks.end();++itr) {
      Ioss::GroupingEntity * entity = *itr;
      const std::string & name = entity->name(); 

      const stk::mesh::Part * part = femMetaData->get_part(name);
      // const CellTopologyData * ct = stk::mesh::fem::get_cell_topology(*part).getCellTopologyData();
      const CellTopologyData * ct = femMetaData->get_cell_topology(*part).getCellTopologyData();

      TEUCHOS_ASSERT(ct!=0);
      mesh.addElementBlock(part->name(),ct);
   }
}

template <typename SetType>
void buildSetNames(const SetType & setData,std::vector<std::string> & names)
{
   // pull out all names for this set
   for(typename SetType::const_iterator itr=setData.begin();itr!=setData.end();++itr) {
      Ioss::GroupingEntity * entity = *itr;
      names.push_back(entity->name());
   }
}

void STK_ExodusReaderFactory::registerSidesets(STK_Interface & mesh,stk::io::util::MeshData & meshData) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::fem::FEMMetaData> metaData = mesh.getMetaData();
   const stk::mesh::PartVector & parts = metaData->get_parts();

   std::cout << std::endl;
   stk::mesh::PartVector::const_iterator partItr;
   for(partItr=parts.begin();partItr!=parts.end();++partItr) {
      const stk::mesh::Part * part = *partItr;
      const stk::mesh::PartVector & subsets = part->subsets();
      // const CellTopologyData * ct = stk::mesh::fem::get_cell_topology(*part).getCellTopologyData();
      const CellTopologyData * ct = metaData->get_cell_topology(*part).getCellTopologyData();

      /*
      if(part->primary_entity_rank()==mesh.getSideRank() && ct==0) {
         std::cout << "found side: \"" << part->name() << "\": ct = " << ct << " subsets = ";
         for(std::size_t i = 0;i<subsets.size();i++)
            std::cout << "\"" << subsets[i]->name() << "\" ";
         std::cout << std::endl;
      }
      */

      // if a side part ==> this is a sideset: now storage is recursive
      // on part contains all sub parts with consistent topology
      if(part->primary_entity_rank()==mesh.getSideRank() && ct==0 && subsets.size()>0) {
         TEST_FOR_EXCEPTION(subsets.size()!=1,std::runtime_error,
                            "STK_ExodusReaderFactory::registerSidesets error - part \"" << part->name() << 
                            "\" has more than one subset"); 

         // grab cell topology and name of subset part
         const stk::mesh::Part * ss_part = subsets[0];
         // const CellTopologyData * ss_ct = stk::mesh::fem::get_cell_topology(*ss_part).getCellTopologyData();
         const CellTopologyData * ss_ct = metaData->get_cell_topology(*ss_part).getCellTopologyData();
 
         // only add subset parts that have no topology
         if(ss_ct!=0) 
            mesh.addSideset(part->name(),ss_ct);
      }
   }
}

}

#endif
