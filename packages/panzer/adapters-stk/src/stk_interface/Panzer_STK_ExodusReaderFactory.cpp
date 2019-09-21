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


#include <Teuchos_TimeMonitor.hpp>
#include <PanzerAdaptersSTK_config.hpp>

#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_Interface.hpp"

#ifdef PANZER_HAVE_IOSS

#include <Ionit_Initializer.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_Region.h>
#include <GetBuckets.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace panzer_stk {

int getMeshDimension(const std::string & meshStr,
                     stk::ParallelMachine parallelMach,
                     const bool isExodus)
{
  stk::io::StkMeshIoBroker meshData(parallelMach);
  meshData.property_add(Ioss::Property("LOWER_CASE_VARIABLE_NAMES", false));
  if (isExodus)
    meshData.add_mesh_database(meshStr, "exodusII", stk::io::READ_MESH);
  else
    meshData.add_mesh_database(meshStr, "pamgen", stk::io::READ_MESH);
  meshData.create_input_mesh();
  return Teuchos::as<int>(meshData.meta_data_rcp()->spatial_dimension());
}

STK_ExodusReaderFactory::STK_ExodusReaderFactory()
  : fileName_(""), restartIndex_(0), isExodus_(true), userMeshScaling_(false), meshScaleFactor_(0.0),
    levelsOfRefinement_(0)
{ }

STK_ExodusReaderFactory::STK_ExodusReaderFactory(const std::string & fileName,
                                                 const int restartIndex,
                                                 const bool isExodus)
  : fileName_(fileName), restartIndex_(restartIndex), isExodus_(isExodus), userMeshScaling_(false), meshScaleFactor_(0.0),
    levelsOfRefinement_(0)
{ }

Teuchos::RCP<STK_Interface> STK_ExodusReaderFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::STK_ExodusReaderFactory::buildMesh()");

   using Teuchos::RCP;
   using Teuchos::rcp;

   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // in here you would add your fields...but it is better to use
   // the two step construction

   // this calls commit on meta data
   const bool buildRefinementSupport = levelsOfRefinement_ > 0 ? true : false;
   mesh->initialize(parallelMach,false,buildRefinementSupport);

   completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

/** This builds all the meta data of the mesh. Does not call metaData->commit.
  * Allows user to add solution fields and other pieces. The mesh can be "completed"
  * by calling <code>completeMeshConstruction</code>.
  */
Teuchos::RCP<STK_Interface> STK_ExodusReaderFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::STK_ExodusReaderFactory::buildUncomittedMesh()");

   using Teuchos::RCP;
   using Teuchos::rcp;

   // read in meta data
   stk::io::StkMeshIoBroker* meshData = new stk::io::StkMeshIoBroker(parallelMach);
   meshData->property_add(Ioss::Property("LOWER_CASE_VARIABLE_NAMES", false));

   // add in "FAMILY_TREE" entity for doing refinement
   std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
   entity_rank_names.push_back("FAMILY_TREE");
   meshData->set_rank_name_vector(entity_rank_names);

   if (isExodus_)
     meshData->add_mesh_database(fileName_, "exodusII", stk::io::READ_MESH);
   else
     meshData->add_mesh_database(fileName_, "pamgen", stk::io::READ_MESH);

   meshData->create_input_mesh();
   RCP<stk::mesh::MetaData> metaData = meshData->meta_data_rcp();

   RCP<STK_Interface> mesh = rcp(new STK_Interface(metaData));
   mesh->initializeFromMetaData();
   mesh->instantiateBulkData(parallelMach);
   meshData->set_bulk_data(mesh->getBulkData());

   // read in other transient fields, these will be useful later when
   // trying to read other fields for use in solve
   meshData->add_all_mesh_fields_as_input_fields();

   // store mesh data pointer for later use in initializing
   // bulk data
   mesh->getMetaData()->declare_attribute_with_delete(meshData);

   // build element blocks
   registerElementBlocks(*mesh,*meshData);
   registerSidesets(*mesh);
   registerNodesets(*mesh);

   mesh->addPeriodicBCs(periodicBCVec_);

   return mesh;
}

void STK_ExodusReaderFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::STK_ExodusReaderFactory::completeMeshConstruction()");

   using Teuchos::RCP;
   using Teuchos::rcp;


   if(not mesh.isInitialized()) {
     const bool buildRefinementSupport = levelsOfRefinement_ > 0 ? true : false;
     mesh.initialize(parallelMach,true,buildRefinementSupport);
   }

   // grab mesh data pointer to build the bulk data
   stk::mesh::MetaData & metaData = *mesh.getMetaData();
   stk::mesh::BulkData & bulkData = *mesh.getBulkData();
   stk::io::StkMeshIoBroker * meshData =
     const_cast<stk::io::StkMeshIoBroker *>(metaData.get_attribute<stk::io::StkMeshIoBroker>());
         // if const_cast is wrong ... why does it feel so right?
         // I believe this is safe since we are basically hiding this object under the covers
         // until the mesh construction can be completed...below I cleanup the object myself.
   TEUCHOS_ASSERT(metaData.remove_attribute(meshData));
      // remove the MeshData attribute

   // build mesh bulk data
   meshData->populate_bulk_data();

   const bool deleteParentElements = true;
   if (levelsOfRefinement_ > 0)
     mesh.refineMesh(levelsOfRefinement_,deleteParentElements);

   // The following section of code is applicable if mesh scaling is
   // turned on from the input file.
   if (userMeshScaling_)
   {
     stk::mesh::Field<double,stk::mesh::Cartesian>* coord_field =
       metaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");

     stk::mesh::Selector select_all_local = metaData.locally_owned_part() | metaData.globally_shared_part();
     stk::mesh::BucketVector const& my_node_buckets = bulkData.get_buckets(stk::topology::NODE_RANK, select_all_local);

     int mesh_dim = mesh.getDimension();

     // Scale the mesh
     const double inv_msf = 1.0/meshScaleFactor_;
     for (size_t i=0; i < my_node_buckets.size(); ++i)
     {
       stk::mesh::Bucket& b = *(my_node_buckets[i]);
       double* coordinate_data = field_data( *coord_field, b );

       for (size_t j=0; j < b.size(); ++j) {
         for (int k=0; k < mesh_dim; ++k) {
           coordinate_data[mesh_dim*j + k] *= inv_msf;
         }
       }
     }
   }

   // put in a negative index and (like python) the restart will be from the back
   // (-1 is the last time step)
   int restartIndex = restartIndex_;
   if(restartIndex<0) {
     std::pair<int,double> lastTimeStep = meshData->get_input_io_region()->get_max_time();
     restartIndex = 1+restartIndex+lastTimeStep.first;
   }

   // populate mesh fields with specific index
   meshData->read_defined_input_fields(restartIndex);

   mesh.buildSubcells();
   mesh.buildLocalElementIDs();

   if (userMeshScaling_) {
     stk::mesh::Field<double,stk::mesh::Cartesian>* coord_field =
       metaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");
     std::vector< const stk::mesh::FieldBase *> fields;
     fields.push_back(coord_field);

     stk::mesh::communicate_field_data(bulkData, fields);
   }

   if(restartIndex>0) // process_input_request is a no-op if restartIndex<=0 ... thus there would be no inital time
      mesh.setInitialStateTime(meshData->get_input_io_region()->get_state_time(restartIndex));
   else
      mesh.setInitialStateTime(0.0); // no initial time to speak, might as well use 0.0

   // clean up mesh data object
   delete meshData;

   // calls Stk_MeshFactory::rebalance
   this->rebalance(mesh);
}

//! From ParameterListAcceptor
void STK_ExodusReaderFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(!paramList->isParameter("File Name"),
        Teuchos::Exceptions::InvalidParameterName,
        "Error, the parameter {name=\"File Name\","
        "type=\"string\""
        "\nis required in parameter (sub)list \""<< paramList->name() <<"\"."
        "\n\nThe parsed parameter parameter list is: \n" << paramList->currentParametersString()
   );

   // Set default values here. Not all the params should be set so this
   // has to be done manually as opposed to using
   // validateParametersAndSetDefaults().
   if(!paramList->isParameter("Restart Index"))
     paramList->set<int>("Restart Index", -1);

   if(!paramList->isParameter("File Type"))
     paramList->set("File Type", "Exodus");

   if(!paramList->isSublist("Periodic BCs"))
     paramList->sublist("Periodic BCs");

   Teuchos::ParameterList& p_bcs = paramList->sublist("Periodic BCs");
   if (!p_bcs.isParameter("Count"))
     p_bcs.set<int>("Count", 0);

   if(!paramList->isParameter("Levels of Uniform Refinement"))
     paramList->set<int>("Levels of Uniform Refinement", 0);

   paramList->validateParameters(*getValidParameters(),0);

   setMyParamList(paramList);

   fileName_ = paramList->get<std::string>("File Name");

   restartIndex_ = paramList->get<int>("Restart Index");

   {
     const auto fileType = paramList->get<std::string>("File Type");
     isExodus_ = fileType == "Exodus";
   }

   // get any mesh scale factor
   if (paramList->isParameter("Scale Factor"))
   {
     meshScaleFactor_ = paramList->get<double>("Scale Factor");
     userMeshScaling_ = true;
   }

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_);

   levelsOfRefinement_ = paramList->get<int>("Levels of Uniform Refinement");
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> STK_ExodusReaderFactory::getValidParameters() const
{
   static Teuchos::RCP<Teuchos::ParameterList> validParams;

   if(validParams==Teuchos::null) {
      validParams = Teuchos::rcp(new Teuchos::ParameterList);
      validParams->set<std::string>("File Name","<file name not set>","Name of exodus file to be read",
                                    Teuchos::rcp(new Teuchos::FileNameValidator));

      validParams->set<int>("Restart Index",-1,"Index of solution to read in",
			    Teuchos::rcp(new Teuchos::AnyNumberParameterEntryValidator(Teuchos::AnyNumberParameterEntryValidator::PREFER_INT,Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes(true))));

      Teuchos::setStringToIntegralParameter<int>("File Type",
                                                 "Exodus",
                                                 "Choose input file type - either \"Exodus\" or \"Pamgen\"",
                                                 Teuchos::tuple<std::string>("Exodus","Pamgen"),
                                                 validParams.get()
                                                 );

      validParams->set<double>("Scale Factor", 1.0, "Scale factor to apply to mesh after read",
                               Teuchos::rcp(new Teuchos::AnyNumberParameterEntryValidator(Teuchos::AnyNumberParameterEntryValidator::PREFER_DOUBLE,Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes(true))));

      Teuchos::ParameterList & bcs = validParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions

      validParams->set("Levels of Uniform Refinement",0,"Number of levels of inline uniform mesh refinement");
   }

   return validParams.getConst();
}

void STK_ExodusReaderFactory::registerElementBlocks(STK_Interface & mesh,stk::io::StkMeshIoBroker & meshData) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> femMetaData = mesh.getMetaData();

   // here we use the Ioss interface because they don't add
   // "bonus" element blocks and its easier to determine
   // "real" element blocks versus STK-only blocks
   const Ioss::ElementBlockContainer & elem_blocks = meshData.get_input_io_region()->get_element_blocks();
   for(Ioss::ElementBlockContainer::const_iterator itr=elem_blocks.begin();itr!=elem_blocks.end();++itr) {
      Ioss::GroupingEntity * entity = *itr;
      const std::string & name = entity->name();

      const stk::mesh::Part * part = femMetaData->get_part(name);
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

void STK_ExodusReaderFactory::registerSidesets(STK_Interface & mesh) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::mesh::PartVector & parts = metaData->get_parts();

   stk::mesh::PartVector::const_iterator partItr;
   for(partItr=parts.begin();partItr!=parts.end();++partItr) {
      const stk::mesh::Part * part = *partItr;
      const stk::mesh::PartVector & subsets = part->subsets();
      // const CellTopologyData * ct = stk::mesh::get_cell_topology(*part).getCellTopologyData();
      const CellTopologyData * ct = metaData->get_cell_topology(*part).getCellTopologyData();

      // if a side part ==> this is a sideset: now storage is recursive
      // on part contains all sub parts with consistent topology
      if(part->primary_entity_rank()==mesh.getSideRank() && ct==0 && subsets.size()>0) {
         TEUCHOS_TEST_FOR_EXCEPTION(subsets.size()!=1,std::runtime_error,
                            "STK_ExodusReaderFactory::registerSidesets error - part \"" << part->name() <<
                            "\" has more than one subset");

         // grab cell topology and name of subset part
         const stk::mesh::Part * ss_part = subsets[0];
         // const CellTopologyData * ss_ct = stk::mesh::get_cell_topology(*ss_part).getCellTopologyData();
         const CellTopologyData * ss_ct = metaData->get_cell_topology(*ss_part).getCellTopologyData();

         // only add subset parts that have no topology
         if(ss_ct!=0)
            mesh.addSideset(part->name(),ss_ct);
      }
   }
}

void STK_ExodusReaderFactory::registerNodesets(STK_Interface & mesh) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::mesh::PartVector & parts = metaData->get_parts();

   stk::mesh::PartVector::const_iterator partItr;
   for(partItr=parts.begin();partItr!=parts.end();++partItr) {
      const stk::mesh::Part * part = *partItr;
      const CellTopologyData * ct = metaData->get_cell_topology(*part).getCellTopologyData();

      // if a side part ==> this is a sideset: now storage is recursive
      // on part contains all sub parts with consistent topology
      if(part->primary_entity_rank()==mesh.getNodeRank() && ct==0) {

         // only add subset parts that have no topology
         if(part->name()!=STK_Interface::nodesString)
            mesh.addNodeset(part->name());
      }
   }
}

}

#endif
