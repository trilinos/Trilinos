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

/** \file   example_Poisson.cpp
    \brief  Example solution of a Poisson equation (div grad u = f) using Trilinos packages.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
 
    \remark Sample command line
    \code   ./example_Poisson.exe 10 10 10 \endcode
*/

#undef DEBUG_PRINTING
// #define DEBUG_PRINTING

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Import.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "../mesh_spec_lt/im_exodusII_l.h"
#include "../mesh_spec_lt/im_ne_nemesisI_l.h"

// AztecOO includes
#include "AztecOO.h"

// ML Includes
#include "ml_MultiLevelPreconditioner.h"


int TestMultiLevelPreconditionerLaplace(char ProblemType[],
				 Teuchos::ParameterList   & MLList,
                                 Epetra_CrsMatrix   & A,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol);


using namespace std;
using namespace Intrepid;
using namespace shards;

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

// Functions to evaluate exact solution and derivatives
double evalu(double & x, double & y, double & z);
int evalGradu(double & x, double & y, double & z, double & gradu1, double & gradu2, double & gradu3);
double evalDivGradu(double & x, double & y, double & z);

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int error = 0;
  int rank=mpiSession.getRank();
  int numProcs=mpiSession.getNProc();
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  int MyPID = Comm.MyPID();

   //Check number of arguments
    TEST_FOR_EXCEPTION( ( argc < 4 ),
                      std::invalid_argument,
                      ">>> ERROR (example_Poisson): Invalid number of arguments. See code listing for requirements.");
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 4)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|          Example: Poisson Equation on Hexahedral Mesh                       |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";

  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;
  
  std::set < topo_entity * , fecomp > edge_set;
  std::set < topo_entity * , fecomp > face_set;

  std::vector < topo_entity * > node_vector;
  std::vector < topo_entity * > edge_vector;
  std::vector < topo_entity * > face_vector;

  std::vector < int > edge_comm_procs;


// ************************************ GET INPUTS **************************************

    int NX            = atoi(argv[1]);  // num intervals in x direction (assumed box domain, 0,1)
    int NY            = atoi(argv[2]);  // num intervals in y direction (assumed box domain, 0,1)
    int NZ            = atoi(argv[3]);  // num intervals in z direction (assumed box domain, 0,1)

// *********************************** CELL TOPOLOGY **********************************

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hex_8(shards::getCellTopologyData<Hexahedron<8> >() );

   // Get dimensions 
    int numNodesPerElem = hex_8.getNodeCount();
    int spaceDim = hex_8.getDimension();
    int dim = 3;

// *********************************** GENERATE MESH ************************************

  if (MyPID == 0) {
    std::cout << "Generating mesh ... \n\n";

    std::cout << "    NX" << "   NY" << "   NZ\n";
    std::cout << std::setw(5) << NX <<
                 std::setw(5) << NY <<
                 std::setw(5) << NZ << "\n\n";
  }

   // Cube
    double leftX = 0.0, rightX = 1.0;
    double leftY = 0.0, rightY = 1.0;
    double leftZ = 0.0, rightZ = 1.0;

   // Create Pamgen input file
    stringstream ss;
    ss.clear();
    ss << "mesh \n";
    ss << "  rectilinear \n"; 
    ss << "     nx = " << NX << "\n";
    ss << "     ny = " << NY << "\n"; 
    ss << "     nz = " << NZ << "\n"; 
    ss << "     bx = 1\n";
    ss << "     by = 1\n"; 
    ss << "     bz = 1\n"; 
    ss << "     gmin = " << leftX << " " << leftY << " " << leftZ << "\n";
    ss << "     gmax = " << rightX << " " << rightY << " " << rightZ << "\n";
    ss << "  end \n";
    ss << "  set assign \n";
    ss << "     sideset, ilo, 1\n"; 
    ss << "     sideset, jlo, 2\n"; 
    ss << "     sideset, klo, 3\n"; 
    ss << "     sideset, ihi, 4\n"; 
    ss << "     sideset, jhi, 5\n"; 
    ss << "     sideset, khi, 6\n"; 
    ss << "  end \n";
    ss << "end \n";

    string meshInput = ss.str();
  if (MyPID == 0) {
    std::cout << meshInput <<"\n";
  }

   // Generate mesh with Pamgen

    long int maxInt = 9223372036854775807LL;
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
    long long numNodesGlobal;
    long long numElemsGlobal;
    long long numElemBlkGlobal;
    long long numNodeSetsGlobal;
    long long numSideSetsGlobal;

    im_ne_get_init_global_l(id, &numNodesGlobal, &numElemsGlobal, 
                         &numElemBlkGlobal, &numNodeSetsGlobal,
                         &numSideSetsGlobal);

    long long * block_ids = new long long [numElemBlk];
    error += im_ex_get_elem_blk_ids_l(id, block_ids);

/*    int edgePerElem = 4;
    if(dim == 3)edgePerElem = 12;
    int facePerElem = 0;
    if(dim == 3)facePerElem = 6;
    FieldContainer<int> elemToEdge(numElems,edgePerElem);
    FieldContainer<int> elemToFace(numElems,facePerElem); */

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


  // Get node-element connectivity
    int telct = 0;
    FieldContainer<int> elemToNode(numElems,numNodesPerElem);
    for(long long b = 0; b < numElemBlk; b++){
      for(long long el = 0; el < elements[b]; el++){
	for (int j=0; j<numNodesPerElem; j++) {
	  elemToNode(telct,j) = elmt_node_linkage[b][el*numNodesPerElem + j]-1;
	}
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
    delete [] nodeCoordx;
    delete [] nodeCoordy;
    delete [] nodeCoordz;

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

      delete [] elem_cmap_ids;
      delete [] elem_cmap_elem_cnts;      
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

  
    //create edges and calculate edge ids
    /*connectivity*/
/*
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
    
    FieldContainer<int> edgeToNode(edge_vector.size(), numNodesPerEdge);
    for(unsigned ect = 0; ect != edge_vector.size(); ect++){
      std::list<long long>::iterator elit;
      int nct = 0;
      for(elit  = edge_vector[ect]->local_node_ids.begin();
	  elit != edge_vector[ect]->local_node_ids.end();
	  elit ++){
	edgeToNode(ect,nct) = *elit;
	nct++;
      }
    }
    int numEdges = edge_vector.size();
    int numFaces = face_vector.size();
   */
 //   std::cout << "    Number of Edges: " << numEdges << " \n";
 //  std::cout << "    Number of Faces: " << numFaces << " \n\n";
   
/*    std::string doing_type;
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
 */ 

//    MPI_Finalize();
//    exit(0);


   // Get boundary (side set) information
   // Side set 1 - left,  Side set 2 - front, Side set 3 - bottom,
   // Side set 4 - right, Side set 5 - back,  Side set 6 - top
    long long * sideSetIds = new long long [numSideSets];
    FieldContainer<int> numElemsOnBoundary(numSideSets);
    long long numSidesinSet;
    long long numDFinSet;
    im_ex_get_side_set_ids_l(id,sideSetIds);
    for (int i=0; i<numSideSets; i++) {
        im_ex_get_side_set_param_l(id,sideSetIds[i],&numSidesinSet,&numDFinSet);
        numElemsOnBoundary(i)=numSidesinSet;
     }

   // Container indicating whether a node is on the boundary (1-yes 0-no)
    FieldContainer<int> nodeOnBoundary(numNodes);

   // Side set 1: left
    if (numElemsOnBoundary(0) > 0){
     long long * sideSetElemList1 = new long long [numElemsOnBoundary(0)];
     long long * sideSetSideList1 = new long long [numElemsOnBoundary(0)];
     im_ex_get_side_set_l(id,sideSetIds[0],sideSetElemList1,sideSetSideList1);
     for (int i=0; i<numElemsOnBoundary(0); i++) {
          nodeOnBoundary(elemToNode(sideSetElemList1[i]-1,0))=1;
          nodeOnBoundary(elemToNode(sideSetElemList1[i]-1,3))=1;
          nodeOnBoundary(elemToNode(sideSetElemList1[i]-1,4))=1;
          nodeOnBoundary(elemToNode(sideSetElemList1[i]-1,7))=1;
     }
     delete [] sideSetElemList1;
     delete [] sideSetSideList1;
   }

   // Side set 2: front
    if (numElemsOnBoundary(1) > 0){
     long long * sideSetElemList2 = new long long [numElemsOnBoundary(1)];
     long long * sideSetSideList2 = new long long [numElemsOnBoundary(1)];
     im_ex_get_side_set_l(id,sideSetIds[1],sideSetElemList2,sideSetSideList2);
     for (int i=0; i<numElemsOnBoundary(1); i++) {
          nodeOnBoundary(elemToNode(sideSetElemList2[i]-1,0))=1;
          nodeOnBoundary(elemToNode(sideSetElemList2[i]-1,1))=1;
          nodeOnBoundary(elemToNode(sideSetElemList2[i]-1,5))=1;
          nodeOnBoundary(elemToNode(sideSetElemList2[i]-1,4))=1;
     }
     delete [] sideSetElemList2;
     delete [] sideSetSideList2;
    }

   // Side set 3: bottom
    if (numElemsOnBoundary(2) > 0){
     long long * sideSetElemList3 = new long long [numElemsOnBoundary(2)];
     long long * sideSetSideList3 = new long long [numElemsOnBoundary(2)];
     im_ex_get_side_set_l(id,sideSetIds[2],sideSetElemList3,sideSetSideList3);
     for (int i=0; i<numElemsOnBoundary(2); i++) {
          nodeOnBoundary(elemToNode(sideSetElemList3[i]-1,0))=1;
          nodeOnBoundary(elemToNode(sideSetElemList3[i]-1,1))=1;
          nodeOnBoundary(elemToNode(sideSetElemList3[i]-1,2))=1;
          nodeOnBoundary(elemToNode(sideSetElemList3[i]-1,3))=1;
     }
     delete [] sideSetElemList3;
     delete [] sideSetSideList3;
    }

   // Side set 4: right
    if (numElemsOnBoundary(3) > 0){
     long long * sideSetElemList4 = new long long [numElemsOnBoundary(3)];
     long long * sideSetSideList4 = new long long [numElemsOnBoundary(3)];
     im_ex_get_side_set_l(id,sideSetIds[3],sideSetElemList4,sideSetSideList4);
     for (int i=0; i<numElemsOnBoundary(3); i++) {
          nodeOnBoundary(elemToNode(sideSetElemList4[i]-1,1))=1;
          nodeOnBoundary(elemToNode(sideSetElemList4[i]-1,2))=1;
          nodeOnBoundary(elemToNode(sideSetElemList4[i]-1,5))=1;
          nodeOnBoundary(elemToNode(sideSetElemList4[i]-1,6))=1;
     }
     delete [] sideSetElemList4;
     delete [] sideSetSideList4;
    }

   // Side set 5: back
    if (numElemsOnBoundary(4) > 0){
     long long * sideSetElemList5 = new long long [numElemsOnBoundary(4)];
     long long * sideSetSideList5 = new long long [numElemsOnBoundary(4)];
     im_ex_get_side_set_l(id,sideSetIds[4],sideSetElemList5,sideSetSideList5);
     for (int i=0; i<numElemsOnBoundary(4); i++) {
          nodeOnBoundary(elemToNode(sideSetElemList5[i]-1,2))=1;
          nodeOnBoundary(elemToNode(sideSetElemList5[i]-1,3))=1;
          nodeOnBoundary(elemToNode(sideSetElemList5[i]-1,6))=1;
          nodeOnBoundary(elemToNode(sideSetElemList5[i]-1,7))=1;
     }
     delete [] sideSetElemList5;
     delete [] sideSetSideList5;
    }

   // Side set 6: top
    if (numElemsOnBoundary(5) > 0){
     long long * sideSetElemList6 = new long long [numElemsOnBoundary(5)];
     long long * sideSetSideList6 = new long long [numElemsOnBoundary(5)];
     im_ex_get_side_set_l(id,sideSetIds[5],sideSetElemList6,sideSetSideList6);
     for (int i=0; i<numElemsOnBoundary(5); i++) {
          nodeOnBoundary(elemToNode(sideSetElemList6[i]-1,4))=1;
          nodeOnBoundary(elemToNode(sideSetElemList6[i]-1,5))=1;
          nodeOnBoundary(elemToNode(sideSetElemList6[i]-1,6))=1;
          nodeOnBoundary(elemToNode(sideSetElemList6[i]-1,7))=1;
     }
     delete [] sideSetElemList6;
     delete [] sideSetSideList6;
    }

    delete [] sideSetIds;


// ************************************ CUBATURE **************************************

  if (MyPID == 0) {
    std::cout << "Getting cubature ... \n\n";
  }

   // Get numerical integration points and weights
    DefaultCubatureFactory<double>  cubFactory;                                   
    int cubDegree = 2;
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(hex_8, cubDegree); 

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);


// ************************************** BASIS ***************************************

  if (MyPID == 0) {
     std::cout << "Getting basis ... \n\n";
  }

   // Define basis 
     Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;
     int numFieldsG = hexHGradBasis.getCardinality();
     FieldContainer<double> hexGVals(numFieldsG, numCubPoints); 
     FieldContainer<double> hexGrads(numFieldsG, numCubPoints, spaceDim); 

  // Evaluate basis values and gradients at cubature points
     hexHGradBasis.getValues(hexGVals, cubPoints, OPERATOR_VALUE);
     hexHGradBasis.getValues(hexGrads, cubPoints, OPERATOR_GRAD);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL STIFFNESS MATRIX *************

  if (MyPID == 0) {
    std::cout << "Building stiffness matrix and right hand side ... \n\n";
  }

 // Settings and data structures for mass and stiffness matrices
    typedef CellTools<double>  CellTools;
    typedef FunctionSpaceTools fst;
    int numCells = 1; 

   // Container for nodes
    FieldContainer<double> hexNodes(numCells, numNodesPerElem, spaceDim);
   // Containers for Jacobian
    FieldContainer<double> hexJacobian(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobDet(numCells, numCubPoints);
   // Containers for element HGRAD stiffness matrix
    FieldContainer<double> localStiffMatrix(numCells, numFieldsG, numFieldsG);
    FieldContainer<double> weightedMeasure(numCells, numCubPoints);
    FieldContainer<double> hexGradsTransformed(numCells, numFieldsG, numCubPoints, spaceDim);
    FieldContainer<double> hexGradsTransformedWeighted(numCells, numFieldsG, numCubPoints, spaceDim);
   // Containers for right hand side vectors
    FieldContainer<double> rhsData(numCells, numCubPoints, cubDim);
    FieldContainer<double> localRHS(numCells, numFieldsG);
    FieldContainer<double> hexGValsTransformed(numCells, numFieldsG, numCubPoints);
    FieldContainer<double> hexGValsTransformedWeighted(numCells, numFieldsG, numCubPoints);
   // Container for cubature points in physical space
    FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);

    // Count owned nodes
    int ownedNodes=0;
    for(int i=0;i<numNodes;i++)
      if(nodeIsOwned[i]) ownedNodes++;

    // Build a list of the OWNED global ids...
    // NTS: will need to switch back to long long
    int *ownedGIDs=new int[ownedNodes];    
    int oidx=0;
    for(int i=0;i<numNodes;i++)
      if(nodeIsOwned[i]){
        ownedGIDs[oidx]=(int)globalNodeIds[i];
        oidx++;
      }
    
    // Global arrays in Epetra format
    Epetra_Map globalMapG(-1,ownedNodes,ownedGIDs,0,Comm);
    Epetra_FECrsMatrix StiffMatrix(Copy, globalMapG, numFieldsG);
    Epetra_FEVector rhs(globalMapG);

    int int_ne=(int)numElems;
    int globalElems;
    Comm.SumAll(&int_ne,&globalElems,1);
    if(!Comm.MyPID()){
      std::cout << " Number of Global Elements: " << globalElems << " \n";
      std::cout << "    Number of Global Nodes: " << globalMapG.NumGlobalElements() << " \n\n";
    }
    
 // *** Element loop ***
    for (int k=0; k<numElems; k++) {

     // Physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

    // Compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInv, hexJacobian );
       CellTools::setJacobianDet(hexJacobDet, hexJacobian );

// ************************** Compute element HGrad stiffness matrices *******************************
  
     // transform to physical coordinates 
      fst::HGRADtransformGRAD<double>(hexGradsTransformed, hexJacobInv, hexGrads);
      
     // compute weighted measure
      fst::computeMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGradsTransformedWeighted,
                                   weightedMeasure, hexGradsTransformed);

     // integrate to compute element stiffness matrix
      fst::integrate<double>(localStiffMatrix,
                             hexGradsTransformed, hexGradsTransformedWeighted, COMP_CPP);

      // assemble into global matrix
      int err = 0;
      for (int row = 0; row < numFieldsG; row++){
        for (int col = 0; col < numFieldsG; col++){
            int rowIndex = globalNodeIds[elemToNode(k,row)];
            int colIndex = globalNodeIds[elemToNode(k,col)];
            double val = localStiffMatrix(0,row,col);
            StiffMatrix.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ******************************* Build right hand side ************************************

      // transform integration points to physical points
       FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);
       CellTools::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);

      // evaluate right hand side functions at physical points
       FieldContainer<double> rhsData(numCells, numCubPoints);
       for (int nPt = 0; nPt < numCubPoints; nPt++){

          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);

          rhsData(0,nPt) = evalDivGradu(x, y, z);
       }

     // transform basis values to physical coordinates 
      fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);
      
     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGValsTransformedWeighted,
                                   weightedMeasure, hexGValsTransformed);

     // integrate rhs term
      fst::integrate<double>(localRHS, rhsData, hexGValsTransformedWeighted,
                             COMP_CPP);

    // assemble into global vector
     for (int row = 0; row < numFieldsG; row++){
           int rowIndex = globalNodeIds[elemToNode(k,row)];
           double val = -localRHS(0,row);
           err = rhs.SumIntoGlobalValues(1, &rowIndex, &val);
      }
 
     
 } // *** end element loop ***


  // Assemble over multiple processors
   StiffMatrix.GlobalAssemble(); StiffMatrix.FillComplete();
   rhs.GlobalAssemble();

 
  // Adjust stiffness matrix and rhs based on boundary conditions
   for (int row = 0; row<numNodes; row++){
       if (nodeOnBoundary(row)) {
          int rowindex = globalNodeIds[row];
          for (int col=0; col<numNodesGlobal; col++){
              double val = 0.0;
              int colindex = col;
              StiffMatrix.ReplaceGlobalValues(1, &rowindex, 1, &colindex, &val);
          }
          double val = 1.0;
          StiffMatrix.ReplaceGlobalValues(1, &rowindex, 1, &rowindex, &val);
          val = 0.0;
          rhs.ReplaceGlobalValues(1, &rowindex, &val);
       }
    }

   
  // Dump matrices to disk
     //   EpetraExt::RowMatrixToMatlabFile("stiff_matrix.dat",StiffMatrix);
     //   EpetraExt::MultiVectorToMatlabFile("rhs_vector.dat",rhs);


   // Run the solver
   Teuchos::ParameterList MLList;
   ML_Epetra::SetDefaults("SA",MLList);
   Epetra_FEVector xexact(globalMapG);
   Epetra_FEVector uh(globalMapG);
   double TotalErrorResidual=0, TotalErrorExactSol=0;

   // Get exact solution at nodes
    for (int i = 0; i<numNodes; i++) {
       if (nodeIsOwned[i]){
          double x = nodeCoord(i,0);
          double y = nodeCoord(i,1);
          double z = nodeCoord(i,2);
          double exactu = evalu(x, y, z);
          int rowindex=globalNodeIds[i];
          xexact.SumIntoGlobalValues(1, &rowindex, &exactu);
       }
    }
    xexact.GlobalAssemble();
       
  //  EpetraExt::MultiVectorToMatlabFile("uexact.dat",xexact);
   
    TestMultiLevelPreconditionerLaplace("laplace",MLList,StiffMatrix,xexact,rhs,uh,
                                       TotalErrorResidual, TotalErrorExactSol);

   // ********  Calculate Error in Solution *************** 
     double L2err = 0.0;
     double L2errTot = 0.0;
     double H1err = 0.0;
     double H1errTot = 0.0;
     double Linferr = 0.0;
     double LinferrTot = 0.0;

   // Import solution onto current processor
     Epetra_Map  solnMap(numNodesGlobal, numNodesGlobal, 0, Comm);
     Epetra_Import  solnImporter(solnMap, globalMapG);
     Epetra_Vector  uCoeff(solnMap);
     uCoeff.Import(uh, solnImporter, Insert);

   // Get cubature points and weights for error calc (may be different from previous)
     DefaultCubatureFactory<double>  cubFactoryErr;                                   
     int cubDegErr = 3;
     Teuchos::RCP<Cubature<double> > hexCubErr = cubFactoryErr.create(hex_8, cubDegErr); 
     int cubDimErr       = hexCubErr->getDimension();
     int numCubPointsErr = hexCubErr->getNumPoints();
     FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
     FieldContainer<double> cubWeightsErr(numCubPointsErr);
     hexCubErr->getCubature(cubPointsErr, cubWeightsErr);

   // Containers for Jacobian
     FieldContainer<double> hexJacobianE(numCells, numCubPointsErr, spaceDim, spaceDim);
     FieldContainer<double> hexJacobInvE(numCells, numCubPointsErr, spaceDim, spaceDim);
     FieldContainer<double> hexJacobDetE(numCells, numCubPointsErr);
     FieldContainer<double> weightedMeasureE(numCells, numCubPointsErr);

  // Evaluate basis values and gradients at cubature points
     FieldContainer<double> uhGVals(numFieldsG, numCubPointsErr); 
     FieldContainer<double> uhGValsTrans(numCells,numFieldsG, numCubPointsErr); 
     FieldContainer<double> uhGrads(numFieldsG, numCubPointsErr, spaceDim); 
     FieldContainer<double> uhGradsTrans(numCells, numFieldsG, numCubPointsErr, spaceDim); 
     hexHGradBasis.getValues(uhGVals, cubPointsErr, OPERATOR_VALUE);
     hexHGradBasis.getValues(uhGrads, cubPointsErr, OPERATOR_GRAD);


   // Loop over elements
    for (int k=0; k<numElems; k++){

      double L2errElem = 0.0;
      double H1errElem = 0.0;
      double uExact; 
      double graduExact1, graduExact2, graduExact3;

     // physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

    // compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobianE, cubPointsErr, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInvE, hexJacobianE );
       CellTools::setJacobianDet(hexJacobDetE, hexJacobianE );

      // transform integration points to physical points
       FieldContainer<double> physCubPoints(numCells,numCubPointsErr, cubDimErr);
       CellTools::mapToPhysicalFrame(physCubPoints, cubPointsErr, hexNodes, hex_8);

      // transform basis values to physical coordinates 
       fst::HGRADtransformVALUE<double>(uhGValsTrans, uhGVals);
       fst::HGRADtransformGRAD<double>(uhGradsTrans, hexJacobInvE, uhGrads);

      // compute weighted measure
       fst::computeMeasure<double>(weightedMeasureE, hexJacobDetE, cubWeightsErr);

      // loop over cubature points
       for (int nPt = 0; nPt < numCubPoints; nPt++){

         // get exact solution and gradients
          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);
          uExact = evalu(x, y, z);
          evalGradu(x, y, z, graduExact1, graduExact2, graduExact3);

         // calculate approximate solution and gradients
          double uApprox = 0.0;
          double graduApprox1 = 0.0;
          double graduApprox2= 0.0;
          double graduApprox3 = 0.0;
          for (int i = 0; i < numFieldsG; i++){
             int rowIndex = globalNodeIds[elemToNode(k,i)];
             double uh1 = uCoeff.Values()[rowIndex];
             uApprox += uh1*uhGValsTrans(0,i,nPt); 
             graduApprox1 += uh1*uhGradsTrans(0,i,nPt,0); 
             graduApprox2 += uh1*uhGradsTrans(0,i,nPt,1); 
             graduApprox3 += uh1*uhGradsTrans(0,i,nPt,2); 
          }

         // evaluate the error at cubature points
          Linferr = max(Linferr, abs(uExact - uApprox));

          L2errElem+=(uExact - uApprox)*(uExact - uApprox)*weightedMeasureE(0,nPt);
          H1errElem+=((graduExact1 - graduApprox1)*(graduExact1 - graduApprox1))
                     *weightedMeasureE(0,nPt);
          H1errElem+=((graduExact2 - graduApprox2)*(graduExact2 - graduApprox2))
                     *weightedMeasureE(0,nPt);
          H1errElem+=((graduExact3 - graduApprox3)*(graduExact3 - graduApprox3))
                     *weightedMeasureE(0,nPt);
        }

       L2err+=L2errElem;
       H1err+=H1errElem;
     }

   // sum over all processors
    Comm.SumAll(&L2err,&L2errTot,1);
    Comm.SumAll(&H1err,&H1errTot,1);
    Comm.MaxAll(&Linferr,&LinferrTot,1);

   if (MyPID == 0) {
    std::cout << "\n" << "L2 Error:  " << sqrt(L2errTot) <<"\n";
    std::cout << "H1 Error:  " << sqrt(H1errTot+L2errTot) <<"\n";
    std::cout << "LInf Error:  " << LinferrTot <<"\n";
   }

   // Cleanup
   for(long long b = 0; b < numElemBlk; b++){     
     delete [] elmt_node_linkage[b];
     delete [] element_types[b];
   }
   delete [] block_ids;
   delete [] nodes_per_element;
   delete [] element_attributes;
   delete [] element_types;
   delete [] elmt_node_linkage;
   delete [] ownedGIDs;
   delete [] elements;
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
   

   // delete mesh
   Delete_Pamgen_Mesh();
   
   // reset format state of std::cout
   std::cout.copyfmt(oldFormatState);
   
   //   MPI_Finalize();
 
   exit(0);

}


// Calculates value of exact solution u
 double evalu(double & x, double & y, double & z)
 {
   // u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)
  // double exactu = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);

  // or

   // u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x+y+z)
   double exactu = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);

   return exactu;
 }

// Calculates gradient of exact solution u
 int evalGradu(double & x, double & y, double & z, double & gradu1, double & gradu2, double & gradu3)
 {
   //  for u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)
  /*      gradu1 = M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
        gradu2 = M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
        gradu3 = M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
  */

  // or

   // for u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x+y+z)
     gradu1 = (M_PI*cos(M_PI*x)+sin(M_PI*x))
                  *sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);
     gradu2 = (M_PI*cos(M_PI*y)+sin(M_PI*y))
                  *sin(M_PI*x)*sin(M_PI*z)*exp(x+y+z);
     gradu3 = (M_PI*cos(M_PI*z)+sin(M_PI*z))
                  *sin(M_PI*x)*sin(M_PI*y)*exp(x+y+z);
  
   return 0;
 }

// Calculates Laplacian of exact solution u
 double evalDivGradu(double & x, double & y, double & z)
 {
   //  for u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)
   //double divGradu = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);

  // or

   // for u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x+y+z)
   double divGradu = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*cos(M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*exp(x+y+z)
                    + 3.0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);
  
   
   return divGradu;
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
  delete [] req;
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


// Test ML
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
				 Teuchos::ParameterList   & MLList,
                                 Epetra_CrsMatrix   & A,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
				 double & TotalErrorExactSol)
{
  Epetra_MultiVector x(xexact);
  x.PutScalar(0.0);
  
  Epetra_LinearProblem Problem(&A,&x,&b); 
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  
  Epetra_Time Time(A.Comm());
  
  // =================== //
  // call ML and AztecOO //
  // =================== //
  
  AztecOO solver(Problem);  
  ML_Epetra::MultiLevelPreconditioner *MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 32);

  solver.Iterate(200, 1e-10);
  
  delete MLPrec;

  uh = *lhs;
  
  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //  
  double d = 0.0, d_tot = 0.0;  
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - xexact[0][i]) * ((*lhs)[0][i] - xexact[0][i]);
  
  A.Comm().SumAll(&d,&d_tot,1);
  
  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A.Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);
  
  string msg = ProblemType;
  
  if (A.Comm().MyPID() == 0) {
    cout << msg << endl << "......Using " << A.Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
  
  return( solver.NumIters() );
  
}


