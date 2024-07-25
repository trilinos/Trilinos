// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Panzer_STK_SculptMeshFactory.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <PanzerAdaptersSTK_config.hpp>

#include "elsa.h"
using Teuchos::RCP;
using Teuchos::rcp;



namespace panzer_stk {

SculptMeshFactory::SculptMeshFactory()
{
   initializeWithDefaults();
}

//! Destructor
SculptMeshFactory::~SculptMeshFactory()
{
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> SculptMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::SculptMeshFactory::buildMesh()");

   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // commit meta data
   mesh->initialize(parallelMach);

   // build bulk data
   completeMeshConstruction(*mesh,parallelMach);

   // wrtie exodus file
   //mesh->writeToExodus("STKSculptMesh.exo");
 
   return mesh;
}

Teuchos::RCP<STK_Interface> SculptMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::SculptMeshFactory::buildUncomittedMesh()");

   RCP<STK_Interface> mesh = rcp(new STK_Interface(3));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   procTuple_ = procRankToProcTuple(machRank_);

   //if( machRank_ == 0 )
   {
       // call Sculptor
       char diatom_file[1000];
       writeDiatomFile( stlFileDir_, stlFileName_, diatom_file  );
  
       callSculptor( parallelMach, diatom_file );
   
        // build meta information: blocks and side set setups
        buildMetaData(parallelMach,*mesh);

        mesh->addPeriodicBCs(periodicBCVec_);
        mesh->setBoundingBoxSearchFlag(useBBoxSearch_);

   }

//   if( machRank_ == 0 )
//                if(mesh->isWritable())
//                               mesh->writeToExodus("STKSculptMesh.exo");
            

   return mesh;
}

int SculptMeshFactory::writeDiatomFile( std::string stl_path, std::string stl_filename, char *diatom_file ) const
{
 
  strcpy( diatom_file, stl_path.c_str() );
  strcat( diatom_file, "stl.diatom" );
  FILE *fp = fopen( diatom_file, "w" );
  if ( !fp )
  {
    printf( "ERROR: Unable to open %s for writing\n", diatom_file );
    return 0;
  }

  char stl_fullfile[1000];
  strcpy( stl_fullfile, stl_path.c_str() );
  strcat( stl_fullfile, stl_filename.c_str() );

  fprintf( fp, "  diatom\n" );
  fprintf( fp, "    package \'box\'\n" );
  fprintf( fp, "      material 1\n" );
  fprintf( fp, "      insert stl\n" );
  fprintf( fp, "        FILE = \'%s\'\n", stl_fullfile );
  fprintf( fp, "      endinsert\n" );
  fprintf( fp, "    endpackage\n" );
  fprintf( fp, "  enddiatom\n" );
  fclose( fp );
  
  return 1;

}
int SculptMeshFactory::callSculptor(stk::ParallelMachine parallelMach, char *diatom_file_name ) const 
{

  char * base_exodus_file_name = NULL;
  char * base_vfrac_file_name = NULL;
  int nelx, nely, nelz; 

  nelx = xInterval_;
  nely = yInterval_;
  nelz = zInterval_;

  int mesh_void = 0;
 
  double gmin[3];
  double gmax[3];

  gmin[0] = xMin_;
  gmin[1] = yMin_;
  gmin[2] = zMin_;
  gmax[0] = xMax_;
  gmax[1] = yMax_;
  gmax[2] = zMax_;

  int stair = 0; 
  int smooth = 1;
  int smooth_iterations = 7;
 
  int gen_sidesets = 4; //for stl based sidesets
  int adaptive_grid = 0;  
  int adapt_level = 2;
  int adapt_type = 0;
   
  printf("\n Sculpt BBox Min ( %lf, %lf, %lf )\n", xMin_, yMin_, zMin_ );
  printf("\n Sculpt BBox Max ( %lf, %lf, %lf )\n", xMax_, yMax_, zMax_ );
  
  int cr_result = Create_Sculptor_Mesh(diatom_file_name,
				       base_exodus_file_name,
                                       base_vfrac_file_name,
				       0, //vfac_input 
				       machSize_, //comm.size(), 
				       machRank_, //comm.rank(),
                                       1,
				       nelx,
				       nely,
				       nelz,
				       gmin,
				       gmax,
				       stair, 
				       smooth, 
				       10,/*num_laplac_iters*/
                                       0, // max opt iters
				       .4,/*opt_threshold*/
                                       0, // max pcol iters
                                       .4, // pcol threshold
				       mesh_void, 
				       gen_sidesets,
				       adapt_type, /* adatptive type*/
                                       adaptive_grid,/*adaptive_grid*/
                                       adapt_level, /* adapt level */
				       0, // max deg iter
                                       0.0,/*double htet_threshold*/
				       0,/*int pillow*/
                                       0, // capture
                                       0, //micro_expand
                                       0, //align
                                       0, //cell_size
				       NULL,/*char * quality_filename*/
				       NULL,/*char * comm_maps_file_name*/
                                       0, // write geom
				       0/*int quiet 1 is quiet*/
				       );



   if (cr_result == 1){
      if(machRank_ == 0)
	printf("Error Generating Sculptor Mesh\n");
        return 1;
   }

  return 0;
}

void SculptMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::SculptMeshFactory::completeMeshConstruction()");

   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   buildElements(parallelMach,mesh);

   mesh.buildSubcells();
   mesh.buildLocalElementIDs();
   mesh.buildLocalEdgeIDs();
   mesh.buildLocalFaceIDs();

   addSideSets(mesh);
   addNodeSets(mesh);
   addEdgeBlocks(mesh);
   addFaceBlocks(mesh);

   this->rebalance(mesh);
}

//! From ParameterListAcceptor
void SculptMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   paramList->validateParametersAndSetDefaults(*getValidParameters(),0);

   setMyParamList(paramList);

   xInterval_ = paramList->get<int>("xInterval");
   yInterval_ = paramList->get<int>("yInterval");
   zInterval_ = paramList->get<int>("zInterval");


   xMin_ = paramList->get<double>("xMin");
   yMin_ = paramList->get<double>("yMin");
   zMin_ = paramList->get<double>("zMin");

   xMax_ = paramList->get<double>("xMax");
   yMax_ = paramList->get<double>("yMax");
   zMax_ = paramList->get<double>("zMax");
 
   stlFileDir_ = paramList->get<std::string>("stlFileDir");
   stlFileName_ = paramList->get<std::string>("stlFileName"); 

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_,useBBoxSearch_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> SculptMeshFactory::getValidParameters() const
{
   static RCP<Teuchos::ParameterList> defaultParams;

   // fill with default values
   if(defaultParams == Teuchos::null) {
      defaultParams = rcp(new Teuchos::ParameterList);

      defaultParams->set<int>("xInterval",10);
      defaultParams->set<int>("yInterval",10);
      defaultParams->set<int>("zInterval",10);

      defaultParams->set<double>("xMin",0.0);
      defaultParams->set<double>("yMin",0.0);
      defaultParams->set<double>("zMin",0.0);

      defaultParams->set<double>("xMax",1.0);
      defaultParams->set<double>("yMax",1.0);
      defaultParams->set<double>("zMax",1.0);
     
      defaultParams->set<std::string>("stlFileDir", "NULL");
      defaultParams->set<std::string>("stlFileName", "NULL");

      Teuchos::ParameterList & bcs = defaultParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
   }


   return defaultParams;
}

void SculptMeshFactory::initializeWithDefaults()
{
   // get valid parameters
   RCP<Teuchos::ParameterList> validParams = rcp(new Teuchos::ParameterList(*getValidParameters()));

   // set that parameter list
   setParameterList(validParams);
   
}

void SculptMeshFactory::buildMetaData(stk::ParallelMachine parallelMach, STK_Interface & mesh) const
{
   struct MeshStorageStruct *mss = get_sculpt_mesh();

   int nBlocks_ = mss->num_elem_blk;
   int nSidesets_ = mss->num_side_sets;
   int nNodesets_ = mss->num_node_sets;


   typedef shards::Hexahedron<8> HexTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<HexTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(2,0);

   const CellTopologyData * edge_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);
   const CellTopologyData * face_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(2,0);


   // build meta data
   //mesh.setDimension(3);
   for( int b = 0; b < nBlocks_; b++){
      std::stringstream ebPostfix;
      ebPostfix << "-" << mss->block_id[b];
      mesh.addElementBlock("eblock"+ebPostfix.str(),ctd);
   }


   // add sidesets 
     int side_set_id;
     machRank_ = stk::parallel_machine_rank(parallelMach); 
     for(int ict = 0;ict < nSidesets_;ict ++){
        std::stringstream sPostfix;
        sPostfix << "-" << mss->side_set_id[ict];
        mesh.addSideset("Sideset"+sPostfix.str(),side_ctd);
    }

   // add nodesets
   for(int nx=0;nx<nNodesets_;nx++) {
     std::stringstream nPostfix;
     nPostfix << "-" << nx;
     mesh.addNodeset("Nodeset"+nPostfix.str());
   }

   mesh.addEdgeBlock(panzer_stk::STK_Interface::edgeBlockString, edge_ctd);
   mesh.addFaceBlock(panzer_stk::STK_Interface::faceBlockString, face_ctd);
}

void SculptMeshFactory::buildNodes( stk::ParallelMachine paralleMach, STK_Interface &mesh ) const
{
   struct MeshStorageStruct *mss = get_sculpt_mesh();
   int num_nodes = mss->num_nodes;


  int dimensionality = 3;

  if (num_nodes){
    int global_node_numbers;
    for(int ict = 0; ict < num_nodes; ict ++){
      global_node_numbers = mss->global_node_numbers[ict];
      std::vector<double> coord(3, 0.0);
      coord[0] = mss->coord[0*num_nodes+ict];
      coord[1] = mss->coord[1*num_nodes+ict];
      coord[2] = mss->coord[2*num_nodes+ict];
      mesh.addNode(global_node_numbers, coord );

      //std::cout<<"Node "<<global_node_numbers<<": ( "<<coord[0]<<", "<<coord[1]<<", "<<coord[2]<<" )"<<std::endl;      

    }
  }
 

} 

void SculptMeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
{
   struct MeshStorageStruct *mss = get_sculpt_mesh();
   int num_blocks  = mss->num_elem_blk;
  

   int *block_id = new int[num_blocks];
   //char ** element_types = new std::string[num_blocks];
   int *elements = new int[num_blocks];
   int *nodes_per_element = new int[num_blocks];
   int *element_attributes = new int[num_blocks];
   int **elmt_node_linkage = new int*[num_blocks];

   for(int b = 0; b < num_blocks; b++){
      block_id[b] = mss->block_id[b];
    //  element_types[b] = mss->element_types[b];
      elements[b] = mss->elements[b];
      nodes_per_element[b] = mss->nodes_per_element[b];
      element_attributes[b] = mss->element_attributes[b];
   }


   int elm_start = 1;
   mesh.beginModification();
   // build each block
   for(int ib=0;ib<num_blocks;ib++) {
     buildBlock(parallelMach,mesh, ib, block_id, elm_start, elements, nodes_per_element, element_attributes, elmt_node_linkage );
     elm_start += elements[ib];
   }
   mesh.endModification();
}

void SculptMeshFactory::buildBlock(stk::ParallelMachine parallelMach,STK_Interface & mesh, int block_index, int *block_id, int elm_start, int *elements, int *nodes_per_element, int *elem_attributes, int **elmt_node_linkage ) const
{

  struct MeshStorageStruct *mss = get_sculpt_mesh();

   // add blocks     
   std::stringstream blockName;
   blockName << "eblock-" << block_id[block_index];
   stk::mesh::Part * block = mesh.getElementBlockPart(blockName.str());


   buildNodes( parallelMach, mesh );

 
    // read element block properties
    //read element connectivity information into a temporary array
      if(elements[block_index]) {
       int maximum_nodes = elements[block_index] * nodes_per_element[block_index];
       elmt_node_linkage[block_index]        = new int[maximum_nodes];
       for(int ict = 0;ict < elements[block_index]; ict ++){
         std::vector<stk::mesh::EntityId> nodes(nodes_per_element[block_index]);
         //std::cout<<"Element id = "<<elm_start+ ict<<std::endl;
         //std::cout<<"Element global id = "<<mss->global_element_numbers[elm_start+ ict-1]<<std::endl;
         for(int nct = 0; nct < nodes_per_element[block_index]; nct++){
            elmt_node_linkage[block_index][(ict*nodes_per_element[block_index])+nct] = mss->elmt_node_linkage[block_index][(ict*nodes_per_element[block_index])+nct];
            nodes[nct] =  mss->global_node_numbers[elmt_node_linkage[block_index][(ict*nodes_per_element[block_index])+nct]-1];
            //std::cout<<" Node linkage id = "<<elmt_node_linkage[block_index][(ict*nodes_per_element[block_index])+nct]<<std::endl;
            //std::cout<<" Node global  id = "<<nodes[nct]<<std::endl;
         }

         stk::mesh::EntityId gid = mss->global_element_numbers[elm_start+ ict-1];
         RCP<ElementDescriptor> ed = rcp(new ElementDescriptor(gid,nodes));
         mesh.addElement(ed,block);
       }
      }
      else {
        elmt_node_linkage[block_index] = NULL;
     }
}

const stk::mesh::Relation * SculptMeshFactory::getRelationByID(unsigned ID,stk::mesh::PairIterRelation relations) const
{
   for(std::size_t i=0;i<relations.size();i++) 
      if(relations[i].identifier()==ID)
         return &relations[i];

   return 0;
}



void SculptMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

    struct MeshStorageStruct *mss = get_sculpt_mesh();
    int num_side_sets  = mss->num_side_sets;

    int *side_set_id = new int[num_side_sets];
    int *num_elements_in_side_set = new int[num_side_sets];
    int *num_nodes_in_side_set = new int[num_side_sets];
    int *num_df_in_side_set = new int[num_side_sets];
    int **side_set_elements = new int*[num_side_sets];
    int **side_set_faces = new int*[num_side_sets];
    //Element_Type **side_set_element_type = new Element_Type*[num_side_sets];
    int **side_set_node_counter = new int*[num_side_sets];
    int **side_set_nodes = new int*[num_side_sets];
    double **side_set_df = new double*[num_side_sets];
    
    for(int ict = 0;ict < num_side_sets;ict ++){
        side_set_id[ict] = mss->side_set_id[ict];
    }

   for(int i = 0; i < num_side_sets; i++) {
     
      std::stringstream sidesetName;
      sidesetName << "Sideset-" << mss->side_set_id[i];
      stk::mesh::Part * sideset = mesh.getSideset(sidesetName.str());


      num_elements_in_side_set[i] = mss->num_elements_in_side_set[i];
      num_df_in_side_set[i] = mss->num_df_in_side_set[i];
      
      int ne = num_elements_in_side_set[i];
      side_set_elements[i] = new int[ne];
      side_set_faces[i] = new int[ne];
      //side_set_element_type[i] = new Element_Type[ne];
      side_set_node_counter[i] = new int[ne];
      side_set_df[i] = new double[num_df_in_side_set[i]];

     
      if(ne) {

        for(int nct = 0; nct < ne; nct ++){

          std::vector<stk::mesh::EntityId> nodes(4);

          int sculpt_elem_id =  mss->global_element_numbers[ mss->side_set_elements[i][nct]-1 ]; 
          int sculpt_face_id = -1 ;

          std::vector<stk::mesh::Entity> localElmts;
          mesh.getMyElements(localElmts);

          std::vector<stk::mesh::Entity>::const_iterator itr;
          for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
            stk::mesh::Entity element = (*itr);

            if( element->identifier() == sculpt_elem_id )
            { 
              sculpt_face_id =  mss->side_set_faces[i][nct];

              stk::mesh::EntityId gid = element->identifier();
 
              stk::mesh::PairIterRelation relations = element->relations(mesh.getSideRank());

              stk::mesh::Entity side = getRelationByID(sculpt_face_id-1,relations)->entity();
 
              if( side != NULL )
              {
                if(side->owner_rank()==machRank_)
                   mesh.addEntityToSideset(*side,sideset );
              }
            }
         }

         if( sculpt_face_id == -1 )
            printf(" ERROR:  Could not find the element id for a sideset face \n");
      }
    }
  }
  mesh.endModification();
}

void SculptMeshFactory::addNodeSets(STK_Interface & mesh) const
{
    mesh.beginModification();

    struct MeshStorageStruct *mss = get_sculpt_mesh();
    int num_node_sets  = mss->num_node_sets;
 

    if (num_node_sets) {
      int *node_set_id = new int[num_node_sets];
      int *num_nodes_in_node_set = new int[num_node_sets];
      int *num_df_in_node_set = new int[num_node_sets];
      int **node_set_nodes = new int*[num_node_sets];
      double **node_set_df = new double*[num_node_sets];

      for(int ict = 0; ict < num_node_sets; ict ++){
        node_set_id[ict] = mss->node_set_id[ict];
        num_nodes_in_node_set[ict] = mss->num_nodes_in_node_set[ict];
        num_df_in_node_set[ict] = mss->num_df_in_node_set[ict];
      }

      for(int i = 0; i < num_node_sets; i++) {
        node_set_nodes[i] = new int[num_nodes_in_node_set[i]];
        node_set_df[i] = NULL;
        if(num_nodes_in_node_set[i]) {
          for(int nct = 0; nct < num_nodes_in_node_set[i];nct ++){
            node_set_nodes[i][nct] = mss->node_set_nodes[i][nct];
          }
        }
      }
    

      for(int i = 0; i < num_node_sets; i++) {
    
        std::stringstream nodesetName;
        nodesetName << "Nodeset-" << mss->node_set_id[i];
        stk::mesh::Part * nodeset = mesh.getNodeset(nodesetName.str());

        for( int j = 0; j < num_nodes_in_node_set[i]; j++ )
        {
           int node_id = node_set_nodes[i][j];
           Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
           if(machRank_==0)
           {  
              stk::mesh::Entity node = bulkData->get_entity(mesh.getNodeRank(),node_id);
              mesh.addEntityToNodeset(*node, nodeset);
           }
        }
      }

    }
    mesh.endModification();
}

void ScupltMeshFactory::addEdgeBlocks(STK_Interface & mesh) const
{
   mesh.beginModification();

   stk::mesh::Part * edge_block = mesh.getEdgeBlock(panzer_stk::STK_Interface::edgeBlockString);

   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();

   std::vector<stk::mesh::Entity> edges;
   bulkData->get_entities(mesh.getEdgeRank(),metaData->locally_owned_part(),edges);
   for(auto edge : edges) {
      mesh.addEntityToEdgeBlock(edge, edge_block);
   }

   mesh.endModification();
}

void ScupltMeshFactory::addFaceBlocks(STK_Interface & mesh) const
{
   mesh.beginModification();

   stk::mesh::Part * face_block = mesh.getFaceBlock(panzer_stk::STK_Interface::faceBlockString);

   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();

   std::vector<stk::mesh::Entity> faces;
   bulkData->get_entities(mesh.getFaceRank(),metaData->locally_owned_part(),faces);
   for(auto face : faces) {
      mesh.addEntityToFaceBlock(face, face_block);
   }

   mesh.endModification();
}

//! Convert processor rank to a tuple
Teuchos::Tuple<std::size_t,2> SculptMeshFactory::procRankToProcTuple(std::size_t procRank) const
{
   std::size_t i=0,j=0;

   j = procRank/machSize_; 
   procRank = procRank % machSize_;
   i = procRank;

   return Teuchos::tuple(i,j);
}

} // end panzer_stk
