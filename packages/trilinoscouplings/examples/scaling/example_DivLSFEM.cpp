// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   example_DivLSFEM.cpp
    \brief  Example of a div-curl system on a hexadedral mesh using div-conforming elements.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.

    \remark Code requires an xml input file with Pamgen mesh description and material parameters
            named DivLSFEMin.xml.
*/

#undef DEBUG_PRINTING
#define DEBUG_PRINTING

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "../mesh_spec_lt/im_exodusII_l.h"
#include "../mesh_spec_lt/im_ne_nemesisI_l.h"

// ML Includes
#include "ml_epetra_utils.h"

using namespace std;
using namespace Intrepid;

class topo_entity{
public:
  topo_entity(){
    local_id = -1;
    global_id = -1;
    owned = true;
  };
  void add_node(long long the_val,long long * global_nids){
    local_node_ids.push_back(the_val);
    sorted_local_node_ids.push_back(the_val);
    sorted_global_node_ids.push_back(global_nids[the_val-1]);
  }
  void sort(){
    sorted_local_node_ids.sort();
    sorted_global_node_ids.sort();
  }
  ~topo_entity(){};
  std::list <long long > local_node_ids;
  std::list <long long > sorted_local_node_ids;
  std::list <long long > sorted_global_node_ids;
  long long local_id;
  long long global_id;
  bool owned;
};

/*******************************************************************************/
inline bool compare_sorted_global_node_ids ( topo_entity* const x,  topo_entity* const y )
/*******************************************************************************/
{
  assert(x->sorted_global_node_ids.size() == y->sorted_global_node_ids.size());
  if(x->sorted_global_node_ids < y->sorted_global_node_ids)return true;    
  return false;
}

struct fecomp{
  bool operator () ( topo_entity* x,  topo_entity*  y )const
  {
    if(x->sorted_local_node_ids < y->sorted_local_node_ids)return true;    
    return false;
  }
};

void  Conform_Boundary_IDS_topo_entity(std::vector < std:: vector < topo_entity * > > & topo_entities,
					long long * proc_ids, 
					long long  rank);

void  Conform_Boundary_IDS(long long ** comm_entities,
			   long long * entity_counts,
			   long long * proc_ids, 
			   long long * data_array,
			   long long num_comm_pairs,
			   long long  rank);

/*******************************************************************************/
void calc_global_node_ids(long long * globalNodeIds,
			  bool * nodeIsOwned,
			  long long numNodes,
			  long long num_node_comm_maps,
			  long long * node_cmap_node_cnts,
			  long long * node_comm_proc_ids,
			  long long * * comm_node_ids,
			  int rank)
/*******************************************************************************/
{
  for(long long i = 0; i < numNodes; i ++){
    globalNodeIds[i] = 1l;
    nodeIsOwned[i] = true;
  }
  for(long long j = 0; j < num_node_comm_maps; j++) {
    for(long long k = 0; k < node_cmap_node_cnts[j] ; k ++){
      if(node_comm_proc_ids[j] < rank){
	globalNodeIds[comm_node_ids[j][k]-1] = -1;	
	nodeIsOwned[comm_node_ids[j][k]-1] = false;
      }
    }
  }
  long long num_unique_nodes = 0;
  for(long long  i = 0 ; i < numNodes; i ++)if(globalNodeIds[i] == 1l)num_unique_nodes ++;
  long long start_id = 0;
  MPI_Scan(&num_unique_nodes,&start_id,1,
	   MPI_LONG_LONG_INT,
	   MPI_SUM,
	   MPI_COMM_WORLD);
  start_id -= num_unique_nodes;

  int num_assigned = 0;
  for(long long  i = 0 ; i < numNodes; i ++)if(globalNodeIds[i] == 1l){
    globalNodeIds[i] = num_assigned + start_id;
    num_assigned ++;
  }

  //Conforms global nodal ids
  Conform_Boundary_IDS(comm_node_ids,
		       node_cmap_node_cnts,
		       node_comm_proc_ids, 
		       globalNodeIds,
		       num_node_comm_maps,
		       rank);

}


/*******************************************************************************/
void calc_global_ids(std::vector < topo_entity * > eof_vec,
		long long **comm_node_ids,
		long long * node_comm_proc_ids,
		long long * node_cmap_node_cnts,
		int num_node_comm_maps,
		int rank,
		std::string fname_string)
/*******************************************************************************/
{
  std::vector < std:: vector < topo_entity *> > topo_entities;
  int nncm = num_node_comm_maps;
  // make a vector of sets of comm nodes
  std::vector < std::set < long long> > comm_vector_set;
  for(int i = 0; i < nncm; i ++){
    std::vector < topo_entity * >v;
    topo_entities.push_back(v);
    std::set < long long > as;
    for(int j = 0; j < node_cmap_node_cnts[i];j ++){
      as.insert(comm_node_ids[i][j]);
    }
    comm_vector_set.push_back(as);
  }
/*run over all edges, faces*/

  for(unsigned tec = 0;tec != eof_vec.size();tec ++){
    topo_entity * teof = eof_vec[tec];

    for(int i = 0; i < nncm; i ++){
      bool found = true;
      std::list <long long > :: iterator lit;
      for( lit = teof->local_node_ids.begin();
	   lit != teof->local_node_ids.end() && found == true;
	   lit ++){
	if(comm_vector_set[i].find(*lit) == comm_vector_set[i].end())found = false;
      }
      //if all component nodes found then add face,edge to comm lists
      if(found){
	topo_entities[i].push_back(teof);
	if(node_comm_proc_ids[i] < rank)teof->owned = false;//not owned if found on lower processor
      }
      else{
      }
    }
  }

  //need to sort the edges_or_face vectors by their sorted global node ids
  for(unsigned i = 0; i < topo_entities.size(); i ++){
    if(!topo_entities[i].empty()){
      std::sort(topo_entities[i].begin(),topo_entities[i].end(),compare_sorted_global_node_ids);
    }
  }
#ifdef DEBUG_PRINTING
 std::stringstream aname;
  aname << fname_string;
  aname << rank;
  aname << ".txt";
  ofstream fout(aname.str().c_str());
#endif

  //need to sort the edges_or_face vectors by their sorted global node ids
  for(unsigned i = 0; i < topo_entities.size(); i ++){
    if(!topo_entities[i].empty()){
      std::sort(topo_entities[i].begin(),topo_entities[i].end(),compare_sorted_global_node_ids);
    }
#ifdef DEBUG_PRINTING
    fout << " from proc rank " << rank 
	 << " to proc rank " << node_comm_proc_ids[i] 
	 << " has " << topo_entities[i].size() 
	 << " entries " << std::endl;

    if(!topo_entities[i].empty()){
      for(unsigned j = 0; j < topo_entities[i].size();j ++){
	topo_entity * eof = topo_entities[i][j];
	for(std::list < long long > :: iterator klit = eof->sorted_global_node_ids.begin();
	    klit != eof->sorted_global_node_ids.end(); klit ++){
	  fout << (*klit) << " " ;
	}
	fout << endl;
      }
    }
#endif
  }

  //count the number of entities owned;
  long long owned_entities = 0;
  for(unsigned ict = 0; ict < eof_vec.size(); ict ++){
    if(eof_vec[ict]->owned)owned_entities ++;
  }
#ifdef DEBUG_PRINTING
  fout << " proc " << rank << " owns " << owned_entities << " edges " << std::endl;
#endif

  long long start_id = 0;
  MPI_Scan(&owned_entities,&start_id,1,
	   MPI_LONG_LONG_INT,
	   MPI_SUM,
	   MPI_COMM_WORLD);
  start_id -= owned_entities;
#ifdef DEBUG_PRINTING
  fout << " proc " << rank << " start_id " << start_id << std::endl;
#endif
  //DMH
  long long num_assigned = 0;
  for(unsigned ict = 0; ict < eof_vec.size(); ict ++){
    if(eof_vec[ict]->owned){
      eof_vec[ict]->global_id = num_assigned + start_id;
      start_id ++;
    }
  }

  Conform_Boundary_IDS_topo_entity(topo_entities,
				    node_comm_proc_ids, 
				    rank);

#ifdef DEBUG_PRINTING
  for(unsigned ict = 0; ict < eof_vec.size(); ict ++){
    fout << "on proc " << rank << " entity " << ict << " has id " << eof_vec[ict]->global_id << std::endl;
  }

  fout.close();
#endif

}

// Functions to evaluate exact solution and its derivatives
int evalu(double & uExact0, 
          double & uExact1, 
          double & uExact2, 
          double & x, 
          double & y, 
          double & z);

double evalDivu(double & x, 
                double & y, 
                double & z);

int evalCurlu(double & curlu0, 
              double & curlu1, 
              double & curlu2, 
              double & x, 
              double & y, 
              double & z, 
              double & mu);

int evalCurlCurlu(double & curlCurlu0, 
                  double & curlCurlu1, 
                  double & curlCurlu2, 
                  double & x, 
                  double & y, 
                  double & z, 
                  double & mu);

int main(int argc, char *argv[]) {

  int error = 0;
#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int rank=mpiSession.getRank();
  int numProcs=mpiSession.getNProc();
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  int MyPID = Comm.MyPID();
#else
  int rank=0;
  int numProcs=1;
  int MyPID = 0;
  Epetra_SerialComm Comm;
#endif
  Epetra_Time Time(Comm);

   //Check number of arguments
    TEST_FOR_EXCEPTION( ( argc < 1 ),
                      std::invalid_argument,
                      ">>> ERROR (example_01): Invalid number of arguments. See code listing for requirements.");
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 1)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|          Example: Div-Curl System on Hexahedral Mesh                        |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";

#ifdef HAVE_MPI
  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;

  std::set < topo_entity * , fecomp > edge_set;
  std::set < topo_entity * , fecomp > face_set;

  std::vector < topo_entity * > edge_vector;
  std::vector < topo_entity * > face_vector;

  std::vector < int > edge_comm_procs;
#endif

  int dim = 3;

// ************************************ GET INPUTS **************************************

  // Input file
    std::string   xmlInFileName = "DivLSFEMin.xml";

  // Read xml file into parameter list
    Teuchos::ParameterList inputList;

   if(xmlInFileName.length()) {
      std::cout << "\nReading parameter list from the XML file \""<<xmlInFileName<<"\" ...\n\n";
      Teuchos::updateParametersFromXmlFile(xmlInFileName,&inputList);
      inputList.print(std::cout,2,true,true);
      std::cout << "\n";
    }
    else
    {
      std::cout << "Cannot read input file: " << xmlInFileName << "\n";
      return 0;
    }

  // Get pamgen mesh definition
    std::string meshInput = Teuchos::getParameter<std::string>(inputList,"meshInput");
 

// *********************************** CELL TOPOLOGY **********************************

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >() );

   // Get dimensions 
    int numNodesPerElem = hex_8.getNodeCount();
    int numEdgesPerElem = hex_8.getEdgeCount();
    int numFacesPerElem = hex_8.getSideCount();
    int numNodesPerEdge = 2;
    int numNodesPerFace = 4;
    int numEdgesPerFace = 4;
    int spaceDim = hex_8.getDimension();

   // Build reference element edge to node map
    FieldContainer<int> refEdgeToNode(numEdgesPerElem,numNodesPerEdge);
    for (int i=0; i<numEdgesPerElem; i++){
        refEdgeToNode(i,0)=hex_8.getNodeMap(1, i, 0);
        refEdgeToNode(i,1)=hex_8.getNodeMap(1, i, 1);
    }

   // Build reference element face to node map
    FieldContainer<int> refFaceToNode(numFacesPerElem,numNodesPerFace);
    for (int i=0; i<numFacesPerElem; i++){
        refFaceToNode(i,0)=hex_8.getNodeMap(2, i, 0);
        refFaceToNode(i,1)=hex_8.getNodeMap(2, i, 1);
        refFaceToNode(i,2)=hex_8.getNodeMap(2, i, 2);
        refFaceToNode(i,3)=hex_8.getNodeMap(2, i, 3);
    }

   // Build reference element face to edge map (Hardcoded for now)
    FieldContainer<int> refFaceToEdge(numFacesPerElem,numEdgesPerFace);
        refFaceToEdge(0,0)=0;
        refFaceToEdge(0,1)=9;
        refFaceToEdge(0,2)=4;
        refFaceToEdge(0,3)=8;
        refFaceToEdge(1,0)=1;
        refFaceToEdge(1,1)=10;
        refFaceToEdge(1,2)=5;
        refFaceToEdge(1,3)=9;
        refFaceToEdge(2,0)=2;
        refFaceToEdge(2,1)=11;
        refFaceToEdge(2,2)=6;
        refFaceToEdge(2,3)=10;
        refFaceToEdge(3,0)=3;
        refFaceToEdge(3,1)=8;
        refFaceToEdge(3,2)=7;
        refFaceToEdge(3,3)=11;
        refFaceToEdge(4,0)=3;
        refFaceToEdge(4,1)=2;
        refFaceToEdge(4,2)=1;
        refFaceToEdge(4,3)=0;
        refFaceToEdge(5,0)=4;
        refFaceToEdge(5,1)=5;
        refFaceToEdge(5,2)=6;
        refFaceToEdge(5,3)=7;

// *********************************** GENERATE MESH ************************************

  if (MyPID == 0) {
    std::cout << "Generating mesh ... \n\n";
  }

   // Generate mesh with Pamgen

    long long maxInt = 9223372036854775807LL;
    Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);
    
   // Get mesh size info
    char title[100];
    long long numDim;
    long long numNodes;
    long long numElems;
    long long numElemBlk;
    long long numNodeSets;
    long long numSideSets;
    int id = 0;

    im_ex_get_init_l(id, title, &numDim, &numNodes, 
                                &numElems, &numElemBlk, &numNodeSets,
                                &numSideSets);
    long long * block_ids = new long long [numElemBlk];
    error += im_ex_get_elem_blk_ids_l(id, block_ids);




    long long  *nodes_per_element   = new long long [numElemBlk];
    long long  *element_attributes  = new long long [numElemBlk];
    long long  *elements            = new long long [numElemBlk];
    char      **element_types       = new char * [numElemBlk];
    long long **elmt_node_linkage   = new long long * [numElemBlk];


    for(long long i = 0; i < numElemBlk; i ++){
      element_types[i] = new char [MAX_STR_LENGTH + 1];
      error += im_ex_get_elem_block_l(id, 
				      block_ids[i], 
				      element_types[i],
				      (long long*)&(elements[i]),
				      (long long*)&(nodes_per_element[i]), 
				      (long long*)&(element_attributes[i]));
    }

    /*connectivity*/
    for(long long b = 0; b < numElemBlk; b++){
      elmt_node_linkage[b] =  new long long [nodes_per_element[b]*elements[b]];
      error += im_ex_get_elem_conn_l(id,block_ids[b],elmt_node_linkage[b]);
    }

   // Get mu value for each block of elements from parameter list
    double  *mu = new double [numElemBlk];
    for(int b = 0; b < numElemBlk; b++){
       stringstream muBlock;
       muBlock.clear();
       muBlock << "mu" << b;
       mu[b] = inputList.get(muBlock.str(),1.0);
    }

  // Get node-element connectivity and set element mu value
    int telct = 0;
    FieldContainer<int> elemToNode(numElems,numNodesPerElem);
    FieldContainer<double> muVal(numElems);
    for(long long b = 0; b < numElemBlk; b++){
      for(long long el = 0; el < elements[b]; el++){
	for (int j=0; j<numNodesPerElem; j++) {
	  elemToNode(telct,j) = elmt_node_linkage[b][el*numNodesPerElem + j]-1;
	}
        muVal(telct) = mu[b];     
	telct ++;
      }
    }
 

   // Read node coordinates and place in field container
    FieldContainer<double> nodeCoord(numNodes,dim);
    double * nodeCoordx = new double [numNodes];
    double * nodeCoordy = new double [numNodes];
    double * nodeCoordz = new double [numNodes];
    im_ex_get_coord_l(id,nodeCoordx,nodeCoordy,nodeCoordz);
    for (int i=0; i<numNodes; i++) {          
      nodeCoord(i,0)=nodeCoordx[i];
      nodeCoord(i,1)=nodeCoordy[i];
      nodeCoord(i,2)=nodeCoordz[i];
    }

#ifdef HAVE_MPI
    /*parallel info*/
    long long num_internal_nodes;
    long long num_border_nodes;
    long long num_external_nodes;
    long long num_internal_elems;
    long long num_border_elems;
    long long num_node_comm_maps;
    long long num_elem_comm_maps;
    im_ne_get_loadbal_param_l( id, 
			       &num_internal_nodes,
			       &num_border_nodes, 
			       &num_external_nodes,
			       &num_internal_elems, 
			       &num_border_elems,
			       &num_node_comm_maps,
			       &num_elem_comm_maps,
			       0/*unused*/ );

    if(num_node_comm_maps > 0){
      node_comm_proc_ids   = new long long  [num_node_comm_maps];
      node_cmap_node_cnts  = new long long  [num_node_comm_maps];
      node_cmap_ids        = new long long  [num_node_comm_maps];
      comm_node_ids        = new long long* [num_node_comm_maps];
      comm_node_proc_ids   = new long long* [num_node_comm_maps];
  
      long long *  elem_cmap_ids        = new long long [num_elem_comm_maps];
      long long *  elem_cmap_elem_cnts  = new long long [num_elem_comm_maps];


      if ( im_ne_get_cmap_params_l( id, 
				  node_cmap_ids,
				  (long long*)node_cmap_node_cnts, 
				  elem_cmap_ids,
				  (long long*)elem_cmap_elem_cnts, 
				  0/*not used proc_id*/ ) < 0 )++error;
      
      for(long long j = 0; j < num_node_comm_maps; j++) {
	comm_node_ids[j]       = new long long [node_cmap_node_cnts[j]];
	comm_node_proc_ids[j]  = new long long [node_cmap_node_cnts[j]];
	if ( im_ne_get_node_cmap_l( id, 
				  node_cmap_ids[j], 
				  comm_node_ids[j], 
				  comm_node_proc_ids[j],
				  0/*not used proc_id*/ ) < 0 )++error;
	node_comm_proc_ids[j] = comm_node_proc_ids[j][0];
      }
    }


    //Calculate global node ids
    long long * globalNodeIds = new long long[numNodes];
    bool * nodeIsOwned = new bool[numNodes];

    calc_global_node_ids(globalNodeIds,
			 nodeIsOwned,
			 numNodes,
			 num_node_comm_maps,
			 node_cmap_node_cnts,
			 node_comm_proc_ids,
			 comm_node_ids,
			 rank);    

#endif 

    FieldContainer<int> elemToEdge(numElems,numEdgesPerElem);
    FieldContainer<int> elemToFace(numElems,numFacesPerElem);

   // calculate edge and face ids
    int elct = 0;
    for(long long b = 0; b < numElemBlk; b++){
      if(nodes_per_element[b] == 4){
      }
      else if (nodes_per_element[b] == 8){
	//loop over all elements and push their edges onto a set if they are not there already
	for(long long el = 0; el < elements[b]; el++){
	  std::set< topo_entity * > ::iterator fit;
	  for (int i=0; i < numEdgesPerElem; i++){
	    topo_entity * teof = new topo_entity;
	    for(int j = 0; j < numNodesPerEdge;j++){
	      teof->add_node(elmt_node_linkage[b][el*numNodesPerElem + refEdgeToNode(i,j)],globalNodeIds);
	    }
	    teof->sort();
	    fit = edge_set.find(teof);
	    if(fit == edge_set.end()){
	      teof->local_id = edge_vector.size();
	      edge_set.insert(teof);
	      elemToEdge(elct,i)= edge_vector.size();
	      edge_vector.push_back(teof);
	    }
	    else{
	      elemToEdge(elct,i) = (*fit)->local_id;
	      delete teof;
	    }
	  }
	  for (int i=0; i < numFacesPerElem; i++){
	    topo_entity * teof = new topo_entity;
	    for(int j = 0; j < numNodesPerFace;j++){
	      teof->add_node(elmt_node_linkage[b][el*numNodesPerElem + refFaceToNode(i,j)],globalNodeIds);
	    }
	    teof->sort();
	    fit = face_set.find(teof);
	    if(fit == face_set.end()){
	      teof->local_id = face_vector.size();
	      face_set.insert(teof);
	      elemToFace(elct,i)= face_vector.size();
	      face_vector.push_back(teof);
	    }
	    else{
	      elemToFace(elct,i) = (*fit)->local_id;
	      delete teof;
	    }
	  }
	  elct ++;
	}
      }
    }
    
   // Edge to Node connectivity
    FieldContainer<int> edgeToNode(edge_vector.size(), numNodesPerEdge);
    for(unsigned ect = 0; ect != edge_vector.size(); ect++){
      std::list<long long>::iterator elit;
      int nct = 0;
      for(elit  = edge_vector[ect]->local_node_ids.begin();
	  elit != edge_vector[ect]->local_node_ids.end();
	  elit ++){
	edgeToNode(ect,nct) = *elit-1;
	nct++;
      }
    }

   // Face to Node connectivity
    FieldContainer<int> faceToNode(face_vector.size(), numNodesPerFace);
    for(unsigned fct = 0; fct != face_vector.size(); fct++){
      std::list<long long>::iterator flit;
      int nct = 0;
      for(flit  = face_vector[fct]->local_node_ids.begin();
	  flit != face_vector[fct]->local_node_ids.end();
	  flit ++){
	faceToNode(fct,nct) = *flit-1;
	nct++;
      }
    }

   // Face to Edge connectivity
    FieldContainer<int> faceToEdge(face_vector.size(), numEdgesPerFace);
    FieldContainer<bool> faceDone(face_vector.size());
    for (int ielem = 0; ielem < numElems; ielem++){
       for (int iface = 0; iface < numFacesPerElem; iface++){
         if (!faceDone(elemToFace(ielem,iface))){
           for (int iedge = 0; iedge < numEdgesPerFace; iedge++){
              faceToEdge(elemToFace(ielem,iface),iedge) = 
                           elemToEdge(ielem,refFaceToEdge(iface,iedge));
              faceDone(elemToFace(ielem,iface))=1;
           }
         }
       }
    }   

    int numEdges = edge_vector.size();
    int numFaces = face_vector.size();

  if (MyPID == 0) {
    std::cout << " Number of Elements: " << numElems << " \n";
    std::cout << "    Number of Nodes: " << numNodes << " \n";
    std::cout << "    Number of Edges: " << numEdges << " \n";
    std::cout << "    Number of Faces: " << numFaces << " \n\n";
  }
   
#ifdef HAVE_MPI
    std::string doing_type;
    doing_type = "EDGES";
    calc_global_ids(edge_vector,
	       comm_node_ids,
	       node_comm_proc_ids, 
	       node_cmap_node_cnts,
	       num_node_comm_maps,
	       rank,
	       doing_type);


    doing_type = "FACES";
    calc_global_ids(face_vector,
	       comm_node_ids,
	       node_comm_proc_ids, 
	       node_cmap_node_cnts,
	       num_node_comm_maps,
	       rank,
	       doing_type);
#endif

 
   // Output element to face connectivity
    ofstream el2fout("elem2face.dat");
    ofstream el2nout("elem2node.dat");
    for (int i=0; i<numElems; i++) {
      for (int l=0; l<numFacesPerElem; l++) {
         el2fout << elemToFace(i,l) << "  ";
      } 
      el2fout << "\n";
      for (int m=0; m<numNodesPerElem; m++) {
        el2nout << elemToNode(i,m) << "  ";
      } 
      el2nout << "\n";
    }
    el2fout.close();
    el2nout.close();

   // Output face to edge and face to node connectivity
    ofstream f2edout("face2edge.dat");
    ofstream f2nout("face2node.dat");
    for (int k=0; k<numFaces; k++) {
       for (int i=0; i<numEdgesPerFace; i++) {
           f2edout << faceToEdge(k,i) << "  ";
       } 
       for (int j=0; j<numNodesPerFace; j++) {
           f2nout << faceToNode(k,j) << "  ";
       } 
       f2edout << "\n";
       f2nout << "\n";
    }
    f2edout.close();
    f2nout.close();

   // Container indicating whether a face is on the boundary (1-yes 0-no)
    FieldContainer<int> edgeOnBoundary(numEdges);
    FieldContainer<int> faceOnBoundary(numFaces);

 // Get boundary (side set) information
    long long * sideSetIds = new long long [numSideSets];
    long long numSidesInSet;
    long long numDFinSet;
    im_ex_get_side_set_ids_l(id,sideSetIds);
    for (int i=0; i<numSideSets; i++) {
        im_ex_get_side_set_param_l(id,sideSetIds[i],&numSidesInSet,&numDFinSet);
        if (numSidesInSet > 0){
          long long * sideSetElemList = new long long [numSidesInSet];
          long long * sideSetSideList = new long long [numSidesInSet];
          im_ex_get_side_set_l(id,sideSetIds[i],sideSetElemList,sideSetSideList);
          for (int j=0; j<numSidesInSet; j++) {
             int iface = sideSetSideList[j]-1;
             faceOnBoundary(elemToFace(sideSetElemList[j]-1,iface))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),0))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),1))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),2))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),3))=1;
          }
          delete [] sideSetElemList;
          delete [] sideSetSideList;
       }
    }

    delete [] sideSetIds;

   //TEMP
    ofstream fFaceout("faceOnBndy.dat");
    for (int i=0; i<numFaces; i++){
       fFaceout << faceOnBoundary(i) <<"\n";
    }
    fFaceout.close();

    ofstream fEdgeout("edgeOnBndy.dat");
    for (int i=0; i<numEdges; i++){
       fEdgeout << edgeOnBoundary(i) <<"\n";
    }
    fEdgeout.close();

    // Print coords
    FILE *f=fopen("coords.dat","w");
    for (int i=0; i < numNodes; i++) {
       fprintf(f,"%22.16e %22.16e %22.16e\n",nodeCoord(i,0),nodeCoord(i,1),nodeCoord(i,2));
    }
    fclose(f);


// **************************** INCIDENCE MATRIX **************************************

   // Edge to face incidence matrix
  if (MyPID == 0) {
    std::cout << "Building incidence matrix ... \n\n";
  }

    Epetra_Map globalMapD(numFaces, 0, Comm);
    Epetra_Map globalMapC(numEdges, 0, Comm);
    Epetra_Map globalMapG(numNodes, 0, Comm);
    Epetra_FECrsMatrix DCurl(Copy, globalMapD, globalMapC, 4);
    Epetra_FECrsMatrix DGrad(Copy, globalMapC, globalMapG, 2);

    double edgevals[2];
    edgevals[0]=-0.5; edgevals[1]=0.5;
    for (int j=0; j<numEdges; j++){
        int rowNum = j;
        int colNum[2];
        colNum[0] = edgeToNode(j,0);
        colNum[1] = edgeToNode(j,1);
        DGrad.InsertGlobalValues(1, &rowNum, 2, colNum, edgevals);
    }

    double vals[4];
    for (int j=0; j<numFaces; j++){
        int rowNum = j;
        int colNum[4];
        colNum[0] = faceToEdge(j,0);
        colNum[1] = faceToEdge(j,1);
        colNum[2] = faceToEdge(j,2);
        colNum[3] = faceToEdge(j,3);

        if (edgeToNode(faceToEdge(j,0),1) == edgeToNode(faceToEdge(j,1),0) ||
               edgeToNode(faceToEdge(j,0),1) == edgeToNode(faceToEdge(j,1),1)){
           vals[0]=0.5;}
        else vals[0]=-0.5;

        if (edgeToNode(faceToEdge(j,1),1) == edgeToNode(faceToEdge(j,2),0) ||
               edgeToNode(faceToEdge(j,1),1) == edgeToNode(faceToEdge(j,2),1)){
           vals[1]=0.5;
        }
        else vals[1]=-0.5;

        if (edgeToNode(faceToEdge(j,2),1) == edgeToNode(faceToEdge(j,3),0) ||
               edgeToNode(faceToEdge(j,2),1) == edgeToNode(faceToEdge(j,3),1)){
           vals[2]=0.5;
        }
        else vals[2]=-0.5;

        if (edgeToNode(faceToEdge(j,3),1) == edgeToNode(faceToEdge(j,0),0) ||
               edgeToNode(faceToEdge(j,3),1) == edgeToNode(faceToEdge(j,0),1)){
           vals[3]=0.5;
        }
        else vals[3]=-0.5;

        DCurl.InsertGlobalValues(1, &rowNum, 4, colNum, vals);
    }


// ************************************ CUBATURE **************************************

   // Get numerical integration points and weights
  if (MyPID == 0) {
    std::cout << "Getting cubature ... \n\n";
  }

    DefaultCubatureFactory<double>  cubFactory;                                   
    int cubDegree = 2;
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(hex_8, cubDegree); 

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);

   // Get numerical integration points and weights for hexahedron face
    //             (needed for rhs boundary term)

    // Define topology of the face parametrization domain as [-1,1]x[-1,1]
    CellTopology paramQuadFace(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

    // Define cubature
    DefaultCubatureFactory<double>  cubFactoryFace;
    Teuchos::RCP<Cubature<double> > hexFaceCubature = cubFactoryFace.create(paramQuadFace, 3);
    int cubFaceDim    = hexFaceCubature -> getDimension();
    int numFacePoints = hexFaceCubature -> getNumPoints();

    // Define storage for cubature points and weights on [-1,1]x[-1,1]
    FieldContainer<double> paramGaussWeights(numFacePoints);
    FieldContainer<double> paramGaussPoints(numFacePoints,cubFaceDim);

    // Define storage for cubature points on workset faces
    hexFaceCubature -> getCubature(paramGaussPoints, paramGaussWeights);


// ************************************** BASIS ***************************************

   // Define basis 
  if (MyPID == 0) {
    std::cout << "Getting basis ... \n\n";
  }
    Basis_HCURL_HEX_I1_FEM<double, FieldContainer<double> > hexHCurlBasis;
    Basis_HDIV_HEX_I1_FEM<double, FieldContainer<double> > hexHDivBasis;
    Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;

    int numFieldsC = hexHCurlBasis.getCardinality();
    int numFieldsD = hexHDivBasis.getCardinality();
    int numFieldsG = hexHGradBasis.getCardinality();

  // Evaluate basis at cubature points
     FieldContainer<double> hexCVals(numFieldsC, numCubPoints, spaceDim); 
     FieldContainer<double> hexDVals(numFieldsD, numCubPoints, spaceDim); 
     FieldContainer<double> hexDivs(numFieldsD, numCubPoints); 
     FieldContainer<double> hexGVals(numFieldsG, numCubPoints); 
     FieldContainer<double> worksetDVals(numFieldsD, numFacePoints, spaceDim); 

     hexHCurlBasis.getValues(hexCVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(hexDVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(hexDivs, cubPoints, OPERATOR_DIV);
     hexHGradBasis.getValues(hexGVals, cubPoints, OPERATOR_VALUE);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL MASS and STIFFNESS MATRICES *************


  if (MyPID == 0) {
    std::cout << "Building mass and stiffness matrices ... \n\n";
  }

 // Settings and data structures for mass and stiffness matrices
    typedef CellTools<double>  CellTools;
    typedef FunctionSpaceTools fst;
    int numCells = 1; 

   // Containers for nodes, edge and face signs 
    FieldContainer<double> hexNodes(numCells, numNodesPerElem, spaceDim);
    FieldContainer<double> hexEdgeSigns(numCells, numFieldsC);
    FieldContainer<double> hexFaceSigns(numCells, numFieldsD);
   // Containers for Jacobian
    FieldContainer<double> hexJacobian(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobDet(numCells, numCubPoints);
   // Containers for element HCURL mass matrix
    FieldContainer<double> massMatrixC(numCells, numFieldsC, numFieldsC);
    FieldContainer<double> weightedMeasure(numCells, numCubPoints);
    FieldContainer<double> weightedMeasureMu(numCells, numCubPoints);
    FieldContainer<double> hexCValsTransformed(numCells, numFieldsC, numCubPoints, spaceDim);
    FieldContainer<double> hexCValsTransformedWeighted(numCells, numFieldsC, numCubPoints, spaceDim);
   // Containers for element HDIV mass matrix
    FieldContainer<double> massMatrixD(numCells, numFieldsD, numFieldsD);
    FieldContainer<double> hexDValsTransformed(numCells, numFieldsD, numCubPoints, spaceDim);
    FieldContainer<double> hexDValsTransformedWeighted(numCells, numFieldsD, numCubPoints, spaceDim);
   // Containers for element HDIV stiffness matrix
    FieldContainer<double> stiffMatrixD(numCells, numFieldsD, numFieldsD);
    FieldContainer<double> hexDivsTransformed(numCells, numFieldsD, numCubPoints);
    FieldContainer<double> hexDivsTransformedWeighted(numCells, numFieldsD, numCubPoints);
   // Containers for element HGRAD mass matrix
    FieldContainer<double> massMatrixG(numCells, numFieldsG, numFieldsG);
    FieldContainer<double> hexGValsTransformed(numCells, numFieldsG, numCubPoints);
    FieldContainer<double> hexGValsTransformedWeighted(numCells, numFieldsG, numCubPoints);
   // Containers for right hand side vectors
    FieldContainer<double> rhsDatag(numCells, numCubPoints, cubDim);
    FieldContainer<double> rhsDatah(numCells, numCubPoints);
    FieldContainer<double> gD(numCells, numFieldsD);
    FieldContainer<double> hD(numCells, numFieldsD);
    FieldContainer<double> gDBoundary(numCells, numFieldsD);
    FieldContainer<double> refGaussPoints(numFacePoints,spaceDim);
    FieldContainer<double> worksetGaussPoints(numCells,numFacePoints,spaceDim);
    FieldContainer<double> worksetJacobians(numCells, numFacePoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobDet(numCells, numFacePoints);
    FieldContainer<double> worksetFaceTu(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetFaceTv(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetFaceN(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetVFieldVals(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetDValsTransformed(numCells, numFieldsD, numFacePoints, spaceDim);
    FieldContainer<double> curluFace(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetDataCrossField(numCells, numFieldsD, numFacePoints, spaceDim);
   // Container for cubature points in physical space
    FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);

    
   // Global arrays in Epetra format
    Epetra_FECrsMatrix MassC(Copy, globalMapC, numFieldsC);
    Epetra_FECrsMatrix MassD(Copy, globalMapD, numFieldsD);
    Epetra_FECrsMatrix MassG(Copy, globalMapG, numFieldsG);
    Epetra_FECrsMatrix StiffD(Copy, globalMapD, numFieldsD);
    Epetra_FEVector rhsD(globalMapD);

    ofstream fSignsout("faceSigns.dat");

 // *** Element loop ***
    for (int k=0; k<numElems; k++) {

     // Physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

     // Face signs
      for (int j=0; j<numFacesPerElem; j++) {
         hexFaceSigns(0,j) = -1.0;
         for (int i=0; i<numNodesPerFace; i++) {
           int indf=i+1;
           if (indf > numNodesPerFace) indf=0;
           if (elemToNode(k,refFaceToNode(j,0))==faceToNode(elemToFace(k,j),i) &&
               elemToNode(k,refFaceToNode(j,1))==faceToNode(elemToFace(k,j),indf))
                hexFaceSigns(0,j) = 1.0;
          }
         fSignsout << hexFaceSigns(0,j) << "  ";
       }
       fSignsout << "\n";

     // Edge signs
      for (int j=0; j<numEdgesPerElem; j++) {
          if (elemToNode(k,refEdgeToNode(j,0))==edgeToNode(elemToEdge(k,j),0) &&
              elemToNode(k,refEdgeToNode(j,1))==edgeToNode(elemToEdge(k,j),1))
              hexEdgeSigns(0,j) = 1.0;
          else 
              hexEdgeSigns(0,j) = -1.0;
       }

    // Compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInv, hexJacobian );
       CellTools::setJacobianDet(hexJacobDet, hexJacobian );


// ************************** Compute element HCurl mass matrices *******************************

     // transform to physical coordinates 
      fst::HCURLtransformVALUE<double>(hexCValsTransformed, hexJacobInv, 
                                   hexCVals);

     // compute weighted measure
      fst::computeMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

     // combine mu value with weighted measure
      for (int nC = 0; nC < numCells; nC++){
        for (int nPt = 0; nPt < numCubPoints; nPt++){
          weightedMeasureMu(nC,nPt) = weightedMeasure(nC,nPt) * muVal(k);
        }
      }

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexCValsTransformedWeighted,
                                   weightedMeasureMu, hexCValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixC,
                             hexCValsTransformed, hexCValsTransformedWeighted,
                             COMP_BLAS);
     // apply edge signs
      fst::applyLeftFieldSigns<double>(massMatrixC, hexEdgeSigns);
      fst::applyRightFieldSigns<double>(massMatrixC, hexEdgeSigns);

     // assemble into global matrix
      for (int row = 0; row < numFieldsC; row++){
        for (int col = 0; col < numFieldsC; col++){
            int rowIndex = elemToEdge(k,row);
            int colIndex = elemToEdge(k,col);
            double val = massMatrixC(0,row,col);
            MassC.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ************************** Compute element HDiv mass matrices *******************************

     // transform to physical coordinates 
      fst::HDIVtransformVALUE<double>(hexDValsTransformed, hexJacobian, hexJacobDet,
                                   hexDVals);

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexDValsTransformedWeighted,
                                   weightedMeasure, hexDValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixD,
                             hexDValsTransformed, hexDValsTransformedWeighted,
                             COMP_BLAS);
     // apply face signs
      fst::applyLeftFieldSigns<double>(massMatrixD, hexFaceSigns);
      fst::applyRightFieldSigns<double>(massMatrixD, hexFaceSigns);

     // assemble into global matrix
      for (int row = 0; row < numFieldsD; row++){
        for (int col = 0; col < numFieldsD; col++){
            int rowIndex = elemToFace(k,row);
            int colIndex = elemToFace(k,col);
            double val = massMatrixD(0,row,col);
            MassD.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ************************ Compute element HDiv stiffness matrices *****************************

      // transform to physical coordinates 
      fst::HDIVtransformDIV<double>(hexDivsTransformed, hexJacobDet,
                                    hexDivs);

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexDivsTransformedWeighted,
                                   weightedMeasure, hexDivsTransformed);

     // integrate to compute element stiffness matrix
      fst::integrate<double>(stiffMatrixD,
                             hexDivsTransformed, hexDivsTransformedWeighted,
                             COMP_BLAS);

     // apply face signs
      fst::applyLeftFieldSigns<double>(stiffMatrixD, hexFaceSigns);
      fst::applyRightFieldSigns<double>(stiffMatrixD, hexFaceSigns);

     // assemble into global matrix
      for (int row = 0; row < numFieldsD; row++){
        for (int col = 0; col < numFieldsD; col++){
            int rowIndex = elemToFace(k,row);
            int colIndex = elemToFace(k,col);
            double val = stiffMatrixD(0,row,col);
            StiffD.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }
// ************************** Compute element HGrad mass matrices *******************************

     // transform to physical coordinates
      fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);

     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGValsTransformedWeighted,
                                   weightedMeasure, hexGValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixG,
                             hexGValsTransformed, hexGValsTransformedWeighted, COMP_BLAS);

      // assemble into global matrix
      for (int row = 0; row < numFieldsG; row++){
        for (int col = 0; col < numFieldsG; col++){
            int rowIndex = elemToNode(k,row);
            int colIndex = elemToNode(k,col);
            double val = massMatrixG(0,row,col);
            MassG.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }


// ******************************* Build right hand side ************************************

      // transform integration points to physical points
       FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);
       CellTools::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);

      // evaluate right hand side functions at physical points
       FieldContainer<double> rhsDatag(numCells, numCubPoints, cubDim);
       FieldContainer<double> rhsDatah(numCells, numCubPoints);
       for (int nPt = 0; nPt < numCubPoints; nPt++){

          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);
          double du1, du2, du3;

          evalCurlCurlu(du1, du2, du3, x, y, z, muVal(k));
          rhsDatag(0,nPt,0) = du1;
          rhsDatag(0,nPt,1) = du2;
          rhsDatag(0,nPt,2) = du3;
         
          rhsDatah(0,nPt) = evalDivu(x, y, z);
       }

     // integrate (g,curl w) term
      fst::integrate<double>(gD, rhsDatag, hexDValsTransformedWeighted,
                             COMP_BLAS);

     // integrate (h,div w) term
      fst::integrate<double>(hD, rhsDatah, hexDivsTransformedWeighted,
                             COMP_BLAS);

     // apply signs
      fst::applyFieldSigns<double>(gD, hexFaceSigns);
      fst::applyFieldSigns<double>(hD, hexFaceSigns);

     // calculate boundary term
      for (int i = 0; i < numFacesPerElem; i++){
        if (faceOnBoundary(elemToFace(k,i))){

         // map Gauss points on quad to reference face: paramGaussPoints -> refGaussPoints
            CellTools::mapToReferenceSubcell(refGaussPoints,
                                   paramGaussPoints,
                                   2, i, hex_8);

         // get basis values at points on reference cell
            hexHDivBasis.getValues(worksetDVals, refGaussPoints, OPERATOR_VALUE);

         // compute Jacobians at Gauss pts. on reference face for all parent cells
            CellTools::setJacobian(worksetJacobians, refGaussPoints,
                         hexNodes, hex_8);
            CellTools::setJacobianDet(worksetJacobDet, worksetJacobians);

         // transform to physical coordinates
            fst::HDIVtransformVALUE<double>(worksetDValsTransformed, worksetJacobians,
                                   worksetJacobDet, worksetDVals);

         // map Gauss points on quad from ref. face to face workset: refGaussPoints -> worksetGaussPoints
            CellTools::mapToPhysicalFrame(worksetGaussPoints,
                                refGaussPoints,
                                hexNodes, hex_8);

         // compute face tangents
            CellTools::getPhysicalFaceTangents(worksetFaceTu,
                                     worksetFaceTv,
                                     paramGaussPoints,
                                     worksetJacobians,
                                     i, hex_8);

         // face outer normals (relative to parent cell) are uTan x vTan
            RealSpaceTools<double>::vecprod(worksetFaceN, worksetFaceTu, worksetFaceTv);


         // evaluate curl u at face points
           for(int nPt = 0; nPt < numFacePoints; nPt++){

             double x = worksetGaussPoints(0, nPt, 0);
             double y = worksetGaussPoints(0, nPt, 1);
             double z = worksetGaussPoints(0, nPt, 2);

             evalCurlu(curluFace(0,nPt,0), curluFace(0,nPt,1), curluFace(0,nPt,2), x, y, z, muVal(k));
           }

         // compute the cross product of curluFace with basis and multiply by weights
           for (int nF = 0; nF < numFieldsD; nF++){
              for(int nPt = 0; nPt < numFacePoints; nPt++){
                  worksetDataCrossField(0,nF,nPt,0) = (curluFace(0,nPt,1)*worksetDValsTransformed(0,nF,nPt,2)
                                 - curluFace(0,nPt,2)*worksetDValsTransformed(0,nF,nPt,1))
                                  * paramGaussWeights(nPt);
                  worksetDataCrossField(0,nF,nPt,1) = (curluFace(0,nPt,2)*worksetDValsTransformed(0,nF,nPt,0)
                                 - curluFace(0,nPt,0)*worksetDValsTransformed(0,nF,nPt,2))
                                  * paramGaussWeights(nPt);
                  worksetDataCrossField(0,nF,nPt,2) = (curluFace(0,nPt,0)*worksetDValsTransformed(0,nF,nPt,1)
                                 - curluFace(0,nPt,1)*worksetDValsTransformed(0,nF,nPt,0))
                                  *paramGaussWeights(nPt);
              } //nPt
           } //nF

          // integrate
           fst::integrate<double>(gDBoundary, worksetFaceN, worksetDataCrossField,
                             COMP_CPP);

          // apply signs
           fst::applyFieldSigns<double>(gDBoundary, hexFaceSigns);

          // add into hC term
            for (int nF = 0; nF < numFieldsD; nF++){
                gD(0,nF) = gD(0,nF) - gDBoundary(0,nF);
            }

        } // if faceOnBoundary
      } // numFaces

    // assemble into global vector
     for (int row = 0; row < numFieldsD; row++){
           int rowIndex = elemToFace(k,row);
           double val = hD(0,row)+gD(0,row);
           rhsD.SumIntoGlobalValues(1, &rowIndex, &val);
     }
 
     
 } // *** end element loop ***

  // Assemble over multiple processors, if necessary
   DCurl.GlobalAssemble(); DCurl.FillComplete(MassD.RowMap(),MassC.RowMap()); 
   DGrad.GlobalAssemble(); DGrad.FillComplete(MassC.RowMap(),MassG.RowMap()); 
   MassC.GlobalAssemble();  MassC.FillComplete();
   MassD.GlobalAssemble();  MassD.FillComplete();
   MassG.GlobalAssemble();  MassG.FillComplete();
   StiffD.GlobalAssemble(); StiffD.FillComplete();
   rhsD.GlobalAssemble();

 // Build the inverse diagonal for MassC
   Epetra_CrsMatrix MassCinv(Copy,MassC.RowMap(),MassC.RowMap(),1);
   Epetra_Vector DiagC(MassC.RowMap());

   DiagC.PutScalar(1.0);
   MassC.Multiply(false,DiagC,DiagC);
   for(int i=0; i<DiagC.MyLength(); i++) {
     DiagC[i]=1.0/DiagC[i];
   }
   for(int i=0; i<DiagC.MyLength(); i++) {
     int CID=MassC.GCID(i);
     MassCinv.InsertGlobalValues(MassC.GRID(i),1,&(DiagC[i]),&CID);
   }
   MassCinv.FillComplete();

  // Set value to zero on diagonal that corresponds to boundary edge
   for(int i=0;i<numEdges;i++) {
     if (edgeOnBoundary(i)){
      double val=0.0;
      MassCinv.ReplaceGlobalValues(i,1,&val,&i);
     }
   }

  // Adjust matrices and rhs due to boundary conditions
     int numBCFaces=0;
     for (int i=0; i<numFaces; i++){
         if (faceOnBoundary(i)){
             numBCFaces++;
         }
      }
      int * BCFaces = new int [numBCFaces];
      int indbc=0;
      for (int i=0; i<numFaces; i++){
         if (faceOnBoundary(i)){
            BCFaces[indbc]=i;
            indbc++;
            rhsD[0][i]=0;
         }
      }
      ML_Epetra::Apply_OAZToMatrix(BCFaces, numBCFaces, StiffD);
      ML_Epetra::Apply_OAZToMatrix(BCFaces, numBCFaces, MassD);
      delete [] BCFaces;

  // Dump matrices to disk
   EpetraExt::RowMatrixToMatlabFile("mag_m1_matrix.dat",MassC);
   EpetraExt::RowMatrixToMatlabFile("mag_m1inv_matrix.dat",MassCinv);
   EpetraExt::RowMatrixToMatlabFile("mag_m2_matrix.dat",MassD);
   EpetraExt::RowMatrixToMatlabFile("mag_k2_matrix.dat",StiffD);
   EpetraExt::RowMatrixToMatlabFile("mag_t0_matrix.dat",DGrad);
   EpetraExt::RowMatrixToMatlabFile("mag_t1_matrix.dat",DCurl);
   EpetraExt::MultiVectorToMatrixMarketFile("rhs2_vector.dat",rhsD,0,0,false);

   fSignsout.close();


 // delete mesh
 Delete_Pamgen_Mesh();

 //clean up
  for(long long b = 0; b < numElemBlk; b++){
     delete [] elmt_node_linkage[b];
     delete [] element_types[b];
   }
   delete [] block_ids;
   delete [] nodes_per_element;
   delete [] element_attributes;
   delete [] element_types;
   delete [] elmt_node_linkage;
   delete [] elements;
   delete [] nodeCoordx;
   delete [] nodeCoordy;
   delete [] nodeCoordz;

#ifdef HAVE_MPI
   delete [] globalNodeIds;
   delete [] nodeIsOwned;
   if(num_node_comm_maps > 0){
      delete [] node_comm_proc_ids;
      delete [] node_cmap_node_cnts;
      delete [] node_cmap_ids;
      for(long long i=0;i<num_node_comm_maps;i++){
        delete [] comm_node_ids[i];
        delete [] comm_node_proc_ids[i];
      }

      delete [] comm_node_ids;
      delete [] comm_node_proc_ids;
   }
#endif


 // reset format state of std::cout
 std::cout.copyfmt(oldFormatState);

   MPI_Finalize();

   exit(0);
}


// Calculates value of exact solution u
 int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z)
 {

   // function 1
    uExact0 = exp(y+z)*(x+1.0)*(x-1.0);
    uExact1 = exp(x+z)*(y+1.0)*(y-1.0);
    uExact2 = exp(x+y)*(z+1.0)*(z-1.0);
  
/*
   // function 2
    uExact0 = cos(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
    uExact1 = cos(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
    uExact2 = cos(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
 
*/
 /*
   // function 3
    uExact0 = x*x-1.0;
    uExact1 = y*y-1.0;
    uExact2 = z*z-1.0;

   // function 4
    uExact0 = sin(M_PI*x);
    uExact1 = sin(M_PI*y);
    uExact2 = sin(M_PI*z);
 */  

   return 0;
 }

// Calculates divergence of exact solution u
 double evalDivu(double & x, double & y, double & z)
 {
   
   // function 1
    double divu = 2.0*x*exp(y+z)+2.0*y*exp(x+z)+2.0*z*exp(x+y);

   // function 2
   //double divu = 2.0*x*cos(M_PI*y)*cos(M_PI*z) + 2.0*y*cos(M_PI*x)*cos(M_PI*z)
   //               + 2.0*z*cos(M_PI*x)*cos(M_PI*y);
   
   // function 3
   // double divu = 2.0*(x + y + z);
   
   // function 4
   // double divu = M_PI*(cos(M_PI*x)+cos(M_PI*y)+cos(M_PI*z));

   return divu;
 }


// Calculates curl of exact solution u
 int evalCurlu(double & curlu0, double & curlu1, double & curlu2, 
                double & x, double & y, double & z, double & mu)
 {
  
   // function 1
    double duxdy = exp(y+z)*(x+1.0)*(x-1.0);
    double duxdz = exp(y+z)*(x+1.0)*(x-1.0);
    double duydx = exp(x+z)*(y+1.0)*(y-1.0);
    double duydz = exp(x+z)*(y+1.0)*(y-1.0);
    double duzdx = exp(x+y)*(z+1.0)*(z-1.0);
    double duzdy = exp(x+y)*(z+1.0)*(z-1.0);
 

  /*
   // function 2
    double duxdy = -M_PI*sin(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
    double duxdz = -M_PI*sin(M_PI*z)*cos(M_PI*y)*(x+1.0)*(x-1.0);
    double duydx = -M_PI*sin(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
    double duydz = -M_PI*sin(M_PI*z)*cos(M_PI*x)*(y+1.0)*(y-1.0);
    double duzdx = -M_PI*sin(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
    double duzdy = -M_PI*sin(M_PI*y)*cos(M_PI*x)*(z+1.0)*(z-1.0);
  */

    curlu0 = (duzdy - duydz)/mu;
    curlu1 = (duxdz - duzdx)/mu;
    curlu2 = (duydx - duxdy)/mu;
 
  /*
   // function 3 and 4
    curlu0 = 0;
    curlu1 = 0;
    curlu2 = 0;
  */
  
   return 0;
 }

// Calculates curl of the curl of exact solution u
 int evalCurlCurlu(double & curlCurlu0, double & curlCurlu1, double & curlCurlu2, 
                    double & x, double & y, double & z, double & mu)
{
   
   // function 1
    double dcurlu0dy = exp(x+y)*(z+1.0)*(z-1.0) - 2.0*y*exp(x+z);
    double dcurlu0dz = 2.0*z*exp(x+y) - exp(x+z)*(y+1.0)*(y-1.0); 
    double dcurlu1dx = 2.0*x*exp(y+z) - exp(x+y)*(z+1.0)*(z-1.0); 
    double dcurlu1dz = exp(y+z)*(x+1.0)*(x-1.0) - 2.0*z*exp(x+y);
    double dcurlu2dx = exp(x+z)*(y+1.0)*(y-1.0) - 2.0*x*exp(y+z);
    double dcurlu2dy = 2.0*y*exp(x+z) - exp(y+z)*(x+1.0)*(x-1.0);
                       

 /*
   // function 2
    double dcurlu0dy = -M_PI*M_PI*cos(M_PI*y)*cos(M_PI*x)*(z+1.0)*(z-1.0)
                           + 2.0*y*M_PI*sin(M_PI*z)*cos(M_PI*x);
    double dcurlu0dz = -2.0*z*M_PI*sin(M_PI*y)*cos(M_PI*x)
                          + M_PI*M_PI*cos(M_PI*z)*cos(M_PI*x)*(y+1.0)*(y-1.0);
    double dcurlu1dx = -2.0*x*M_PI*sin(M_PI*z)*cos(M_PI*y)
                          + M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
    double dcurlu1dz = -M_PI*M_PI*cos(M_PI*z)*cos(M_PI*y)*(x+1.0)*(x-1.0)
                           + 2.0*z*M_PI*sin(M_PI*x)*cos(M_PI*y);
    double dcurlu2dx = -M_PI*M_PI*cos(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0)
                           + 2.0*x*M_PI*sin(M_PI*y)*cos(M_PI*z);
    double dcurlu2dy = -2.0*y*M_PI*sin(M_PI*x)*cos(M_PI*z)
                          + M_PI*M_PI*cos(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
 */
                       
    curlCurlu0 = (dcurlu2dy - dcurlu1dz)/mu;
    curlCurlu1 = (dcurlu0dz - dcurlu2dx)/mu;
    curlCurlu2 = (dcurlu1dx - dcurlu0dy)/mu;
 
 /*
   // function 3 and 4
    curlCurlu0 = 0.0;
    curlCurlu1 = 0.0;
    curlCurlu2 = 0.0;
 */

    return 0;
}

/*****************************************************************************/
void  Conform_Boundary_IDS(long long ** comm_entities,
			   long long * entity_counts,
			   long long * proc_ids, 
			   long long * data_array,
			   long long num_comm_pairs,
			   long long  rank)
/*****************************************************************************/
{
  //relies on data array having a -1 for unassigned values
  unsigned nncm = num_comm_pairs;

  //Will load an ownership flag along with the data to allow the conform
  MPI_Request * req = new MPI_Request[nncm];
  
  long long ** send_buffer = new long long * [nncm];
  long long ** receive_buffer = new long long * [nncm];
  for(unsigned i = 0; i < nncm; i ++)send_buffer[i]    = new long long [entity_counts[i]];
  for(unsigned i = 0; i < nncm; i ++)receive_buffer[i] = new long long [entity_counts[i]];

  // load up send buffer
  for(unsigned i = 0; i < nncm; i ++){
    for(unsigned j = 0; j < entity_counts[i];j++){
      send_buffer[i][j] = data_array[comm_entities[i][j]-1];
    }
  }

  //communicate

  for(unsigned i = 0; i < nncm ;i ++){
    int size = entity_counts[i];
    int proc = proc_ids[i];
    MPI_Irecv(receive_buffer[i],size, MPI_LONG_LONG_INT, proc, 1, MPI_COMM_WORLD, req + i);
  }

  for(unsigned i = 0; i < nncm ;i ++){
    int size = entity_counts[i];
    int proc = proc_ids[i];
    MPI_Send(send_buffer[i], size, MPI_LONG_LONG_INT, proc, 1,MPI_COMM_WORLD);
  }

  for(unsigned i = 0; i < nncm ;i ++){
    MPI_Status stat;
    MPI_Wait(req + i, &stat);
  }

  for(unsigned i = 0; i < nncm; i ++){
    for(unsigned j = 0; j < entity_counts[i];j++){
      if(receive_buffer[i][j] >= 0)data_array[comm_entities[i][j]-1] = receive_buffer[i][j];
    }
  }
  
  for(unsigned i = 0; i < nncm; i ++){
    if(send_buffer[i])   delete [] send_buffer[i];
    if(receive_buffer[i])delete [] receive_buffer[i];
  }
  delete [] send_buffer;
  delete [] receive_buffer;
}


/*****************************************************************************/
void  Conform_Boundary_IDS_topo_entity(std::vector < std:: vector < topo_entity * > > & topo_entities,
				       long long * proc_ids, 
				       long long  rank)
/*****************************************************************************/
{
  //relies on data array having a -1 for unassigned values
  unsigned nncm = topo_entities.size();

  //Will load an ownership flag along with the data to allow the conform
  MPI_Request * req = new MPI_Request[nncm];
  
  long long ** send_buffer = new long long * [nncm];
  long long ** receive_buffer = new long long * [nncm];
  for(unsigned i = 0; i < nncm; i ++){
    if(topo_entities[i].size() > 0){
      send_buffer[i]    = new long long [topo_entities[i].size()];
      receive_buffer[i] = new long long [topo_entities[i].size()];
    }
    else{
      send_buffer[i] = NULL;
      receive_buffer[i] = NULL;
    }
  }

  // load up send buffer
  for(unsigned i = 0; i < nncm; i ++){
    for(unsigned j = 0; j < topo_entities[i].size();j++){
      send_buffer[i][j] = topo_entities[i][j]->global_id;
    }
  }

  //communicate

  for(unsigned i = 0; i < nncm ;i ++){
    int size = topo_entities[i].size();
    if(size > 0){
      int proc = proc_ids[i];
      MPI_Irecv(receive_buffer[i],size, MPI_LONG_LONG_INT, proc, 1, MPI_COMM_WORLD, req + i);
    }
  }

  for(unsigned i = 0; i < nncm ;i ++){
    int size = topo_entities[i].size();
    if(size > 0){
      int proc = proc_ids[i];
      MPI_Send(send_buffer[i], size, MPI_LONG_LONG_INT, proc, 1,MPI_COMM_WORLD);
    }
  }

  for(unsigned i = 0; i < nncm ;i ++){
    MPI_Status stat;
    int size = topo_entities[i].size();
    if(size > 0){
      MPI_Wait(req + i, &stat);
    }
  }

  for(unsigned i = 0; i < nncm; i ++){
    for(unsigned j = 0; j < topo_entities[i].size();j++){
      if(receive_buffer[i][j] >= 0)topo_entities[i][j]->global_id = receive_buffer[i][j];
    }
  }

  for(unsigned i = 0; i < nncm; i ++){
    if(send_buffer[i])   delete [] send_buffer[i];
    if(receive_buffer[i])delete [] receive_buffer[i];
  }
  delete [] send_buffer;
  delete [] receive_buffer;
}
